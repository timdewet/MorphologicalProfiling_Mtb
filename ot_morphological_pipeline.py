"""
Optimal Transport Morphological Profiling Pipeline
===================================================
Computes pairwise Wasserstein (Sinkhorn) distances between M. tuberculosis
conditions (CRISPRi knockdowns and drug treatments) using cell-level
morphological features from MicrobeJ.  Ranks drug–gene matches by phenotypic
similarity and optionally tests significance via permutation.

Usage:
    python ot_morphological_pipeline.py
"""

import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ot
import pandas as pd
import seaborn as sns
import umap
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import cosine, euclidean
from scipy.stats import false_discovery_control
from sklearn.decomposition import PCA

# ==============================================================================
# CONFIG
# ==============================================================================

CONFIG = {
    # --- Input ---
    "input_files": [
        "input_data/data_extraction_18_03_25.csv",
        "input_data/all_morphology_combined.csv",
    ],
    "name_separator": "__",
    "control_labels": ["NT", "No_drug"],
    "exclude_reporters": [],

    # --- Name corrections ---
    # Applied after parsing.  Each dict must have "experiment" (exact match on
    # the raw EXPERIMENT string) plus any fields to overwrite.
    "name_corrections": [
        {"experiment": "WT_Reporters_+_drug__imiB__Inn",     "reporter": "iniB", "knockdown": "INH"},
        {"experiment": "WT_Reporters_+_drug__imiB__EMB",     "reporter": "iniB"},
        {"experiment": "WT_Reporters_+_drug__imiB__No_drug", "reporter": "iniB"},
        {"experiment": "WT_Reporters_+_drug__imiB__RIF",     "reporter": "iniB"},
        {"experiment": "ATC_Strains__recA__dnaW2",           "knockdown": "dnaN1_rep2"},
    ],

    # --- Features (29 SHAPE columns, matching R pipeline) ---
    "shape_features": [
        "SHAPE.angularity",
        "SHAPE.angularity.amplitude",
        "SHAPE.angularity.max",
        "SHAPE.angularity.median",
        "SHAPE.angularity.mid",
        "SHAPE.angularity.min",
        "SHAPE.angularity.stdev",
        "SHAPE.angularity.variation",
        "SHAPE.area",
        "SHAPE.aspectRatio",
        "SHAPE.circularity",
        "SHAPE.curvature",
        "SHAPE.feret",
        "SHAPE.feret.max",
        "SHAPE.feret.min",
        "SHAPE.length",
        "SHAPE.perimeter",
        "SHAPE.pole",
        "SHAPE.roundness",
        "SHAPE.sinuosity",
        "SHAPE.solidity",
        "SHAPE.width",
        "SHAPE.width.amplitude",
        "SHAPE.width.max",
        "SHAPE.width.median",
        "SHAPE.width.mid",
        "SHAPE.width.min",
        "SHAPE.width.stdev",
        "SHAPE.width.variation",
    ],
    "include_fluorescence": False,
    "fluorescence_col": "INTENSITY.ch1.mean",

    # --- Normalisation ---
    "normalize_within_reporter": True,
    # Subtract each condition's centroid so OT compares the *shape* of the
    # perturbation (direction + spread) rather than absolute position.
    # Without this, mild perturbations all match each other simply because
    # they sit near the control origin.
    "center_conditions": False,

    # --- Distance mode ---
    # "cell_ot": Sinkhorn on raw cell distributions (original approach)
    # "sscore_cosine": S-scores (mean + CV per feature per condition,
    #   standardised vs controls), then cosine distance between condition
    #   profiles.  Mirrors the R pipeline's approach and avoids the
    #   within-condition-variance domination problem.
    # "sscore_euclidean": Same S-scores but Euclidean distance.
    "distance_mode": "sscore_cosine",

    # --- Dimensionality reduction (only used in cell_ot mode) ---
    # PCA on all non-control cells before OT.  Decorrelates features and
    # reduces curse-of-dimensionality effects on the cost matrix.
    "pca_before_ot": True,
    "pca_variance_threshold": 0.85,  # retain PCs explaining this fraction

    # --- OT parameters (only used in cell_ot mode) ---
    "sinkhorn_reg": 0.05,
    "sinkhorn_max_iter": 1000,
    "subsample_n": 500,
    "random_seed": 42,

    # --- Drug identification ---
    "drug_experiment_pattern": "drug",  # matched case-insensitively on experiment_type

    # --- Permutation test (set n_permutations > 0 to enable) ---
    "n_permutations": 1000,
    "permutation_top_k": 5,       # only test top-K gene matches per drug
    "permutation_subsample_n": 200,
    "alpha": 0.05,

    # --- Output ---
    "fig_dir": "figures",
    "output_csv": "ot_ranked_matches.csv",
    "distance_csv": "ot_distance_matrix.csv",
    "fig_dpi": 300,
}


# ==============================================================================
# DATA LOADING & PARSING
# ==============================================================================

def load_and_parse_data(config):
    """Load CSVs, parse EXPERIMENT column, apply corrections, filter."""
    frames = []
    for fpath in config["input_files"]:
        df = pd.read_csv(fpath)
        df["_source_file"] = Path(fpath).name
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    print(f"Loaded {len(df):,} cells from {len(config['input_files'])} file(s)")

    # Parse EXPERIMENT into components
    sep = config["name_separator"]
    parts = df["EXPERIMENT"].str.split(sep, n=2, expand=True)
    df["experiment_type"] = parts[0]
    df["reporter"] = parts[1]
    df["knockdown"] = parts[2]

    # Apply name corrections
    for corr in config["name_corrections"]:
        mask = df["EXPERIMENT"] == corr["experiment"]
        for field in ("experiment_type", "reporter", "knockdown"):
            if field in corr:
                df.loc[mask, field] = corr[field]

    # Flag controls
    ctrl_pat = "^(" + "|".join(re.escape(c) for c in config["control_labels"]) + ")(_.+)?$"
    df["is_control"] = df["knockdown"].str.match(ctrl_pat)

    # Flag drugs
    df["is_drug"] = df["experiment_type"].str.contains(
        config["drug_experiment_pattern"], case=False, na=False
    ) & ~df["is_control"]

    # Exclude reporters
    if config["exclude_reporters"]:
        before = len(df)
        df = df[~df["reporter"].isin(config["exclude_reporters"])].copy()
        print(f"Excluded reporters {config['exclude_reporters']}: "
              f"{before - len(df):,} cells removed, {len(df):,} remaining")

    # Condition label = knockdown name (pools across reporters)
    df["condition_label"] = df["knockdown"]

    # Summary
    n_ctrl = df["is_control"].sum()
    n_drug = df["is_drug"].sum()
    n_gene = (~df["is_control"] & ~df["is_drug"]).sum()
    conditions = df.loc[~df["is_control"], "condition_label"].nunique()
    print(f"Conditions: {conditions} non-control "
          f"({df.loc[df['is_drug'], 'condition_label'].nunique()} drugs, "
          f"{df.loc[~df['is_control'] & ~df['is_drug'], 'condition_label'].nunique()} genes)")
    print(f"Cells: {n_gene:,} gene-KD, {n_drug:,} drug, {n_ctrl:,} control")

    return df


# ==============================================================================
# NORMALISATION
# ==============================================================================

def normalize_features(df, features, config):
    """Z-score normalise features using control cells as baseline."""
    df = df.copy()

    if config["normalize_within_reporter"]:
        for reporter, grp in df.groupby("reporter"):
            ctrl_mask = grp["is_control"]
            if ctrl_mask.sum() == 0:
                warnings.warn(f"No control cells for reporter '{reporter}'; "
                              "using global controls for this group.")
                ctrl_data = df.loc[df["is_control"], features]
            else:
                ctrl_data = grp.loc[ctrl_mask, features]
            mu = ctrl_data.mean()
            sigma = ctrl_data.std().replace(0, 1)
            df.loc[grp.index, features] = (grp[features] - mu) / sigma
        print("Normalised features within each reporter background")
    else:
        ctrl_data = df.loc[df["is_control"], features]
        mu = ctrl_data.mean()
        sigma = ctrl_data.std().replace(0, 1)
        df[features] = (df[features] - mu) / sigma
        print("Normalised features using pooled controls")

    return df


# ==============================================================================
# S-SCORE COMPUTATION
# ==============================================================================

def _cv(x):
    """Coefficient of variation, handling edge cases."""
    x = x.dropna()
    if len(x) == 0:
        return np.nan
    mu = x.mean()
    if mu == 0 or not np.isfinite(mu):
        return np.nan
    return x.std() / mu


def compute_s_scores(df, features, config):
    """Compute S-scores (mean + CV) per condition, standardised vs controls.

    For each condition, computes the mean and CV of each raw feature across
    its cells.  These are then z-scored using the distribution of means/CVs
    across control samples, yielding a profile vector per condition.

    Returns a DataFrame with one row per non-control condition_label and
    columns for each S-score feature (feature_mean, feature_CV).
    """
    # Group by EXPERIMENT (original sample) for mean/CV computation,
    # since condition_label pools replicates
    sample_means = df.groupby("EXPERIMENT")[features].mean()
    sample_cvs = df.groupby("EXPERIMENT")[features].apply(
        lambda g: g.apply(_cv)
    )

    # Get metadata per EXPERIMENT
    meta = df.drop_duplicates("EXPERIMENT").set_index("EXPERIMENT")[
        ["condition_label", "is_control", "is_drug", "reporter"]
    ]

    # Control baselines
    ctrl_exps = meta[meta["is_control"]].index

    if config.get("normalize_within_reporter", False):
        # Compute baselines per reporter
        scored_parts = []
        for reporter, rep_meta in meta.groupby("reporter"):
            rep_ctrl = rep_meta[rep_meta["is_control"]].index
            if len(rep_ctrl) == 0:
                rep_ctrl = ctrl_exps  # fallback to global

            for stat_name, stat_df in [("mean", sample_means), ("CV", sample_cvs)]:
                ctrl_vals = stat_df.loc[stat_df.index.isin(rep_ctrl)]
                mu = ctrl_vals.mean()
                sigma = ctrl_vals.std().replace(0, 1)

                rep_exps = rep_meta.index
                scored = (stat_df.loc[stat_df.index.isin(rep_exps)] - mu) / sigma
                scored.columns = [f"{f}_{stat_name}" for f in features]
                scored_parts.append(scored)

        all_scored = pd.concat(scored_parts, axis=1)
        # Average duplicate columns (same feature scored from different reporters)
        all_scored = all_scored.T.groupby(level=0).mean().T
    else:
        scored_parts = []
        for stat_name, stat_df in [("mean", sample_means), ("CV", sample_cvs)]:
            ctrl_vals = stat_df.loc[stat_df.index.isin(ctrl_exps)]
            mu = ctrl_vals.mean()
            sigma = ctrl_vals.std().replace(0, 1)
            scored = (stat_df - mu) / sigma
            scored.columns = [f"{f}_{stat_name}" for f in features]
            scored_parts.append(scored)
        all_scored = pd.concat(scored_parts, axis=1)

    # Attach condition_label and aggregate by it (pools across reporters)
    all_scored["condition_label"] = meta.loc[all_scored.index, "condition_label"].values
    all_scored["is_control"] = meta.loc[all_scored.index, "is_control"].values
    all_scored["is_drug"] = meta.loc[all_scored.index, "is_drug"].values

    # Remove controls, average across experiments within same condition
    non_ctrl = all_scored[~all_scored["is_control"]].copy()
    score_cols = [c for c in non_ctrl.columns
                  if c not in ("condition_label", "is_control", "is_drug")]
    profiles = non_ctrl.groupby("condition_label")[score_cols].mean()

    # Preserve is_drug flag
    drug_flag = non_ctrl.drop_duplicates("condition_label").set_index("condition_label")["is_drug"]
    profiles["is_drug"] = drug_flag

    n_features = len(score_cols)
    print(f"S-scores: {n_features} features "
          f"({len(features)} means + {len(features)} CVs) "
          f"for {len(profiles)} conditions")

    return profiles, score_cols


def compute_sscore_distances(profiles, score_cols, metric="cosine"):
    """Pairwise distances between S-score condition profiles."""
    conditions = profiles.index.tolist()
    n = len(conditions)
    dist_func = cosine if metric == "cosine" else euclidean

    mat = profiles[score_cols].values
    # Replace NaN with 0 for distance computation
    mat = np.nan_to_num(mat, nan=0.0)

    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = dist_func(mat[i], mat[j])
            dist_matrix[i, j] = dist_matrix[j, i] = d

    dist_df = pd.DataFrame(dist_matrix, index=conditions, columns=conditions)
    print(f"Computed {n*(n-1)//2} pairwise {metric} distances")
    return dist_df


# ==============================================================================
# OPTIMAL TRANSPORT DISTANCES
# ==============================================================================

def compute_ot_distances(df, features, conditions, config):
    """Pairwise Sinkhorn distances between all non-control conditions."""
    rng = np.random.default_rng(config["random_seed"])
    n_sub = config["subsample_n"]
    n = len(conditions)

    # Pre-extract and subsample
    arrays = {}
    for cond in conditions:
        X = df.loc[df["condition_label"] == cond, features].values
        if len(X) == 0:
            raise ValueError(f"No cells for condition '{cond}'")
        if len(X) >= n_sub:
            idx = rng.choice(len(X), n_sub, replace=False)
        else:
            idx = rng.choice(len(X), n_sub, replace=True)
        arrays[cond] = X[idx].astype(np.float64)

    # Centre each condition's distribution at the origin so OT compares
    # the *shape* of the perturbation (spread, skew, multimodality) rather
    # than absolute position in feature space.
    if config.get("center_conditions", False):
        for cond in conditions:
            arrays[cond] = arrays[cond] - arrays[cond].mean(axis=0)
        print("Centred each condition's distribution at the origin")

    dist_matrix = np.zeros((n, n))
    n_pairs = n * (n - 1) // 2
    print(f"Computing {n_pairs} pairwise Sinkhorn distances "
          f"({n} conditions, {n_sub} cells each) ...")

    pair_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            X_a, X_b = arrays[conditions[i]], arrays[conditions[j]]
            M = ot.dist(X_a, X_b, metric="sqeuclidean")
            M_max = M.max()
            if M_max > 0:
                M /= M_max
            a = np.ones(len(X_a)) / len(X_a)
            b = np.ones(len(X_b)) / len(X_b)
            w = ot.sinkhorn2(
                a, b, M,
                reg=config["sinkhorn_reg"],
                numItermax=config["sinkhorn_max_iter"],
                warn=False,
            )
            dist_matrix[i, j] = dist_matrix[j, i] = float(w)
            pair_count += 1
            if pair_count % 50 == 0 or pair_count == n_pairs:
                print(f"  {pair_count}/{n_pairs} pairs done")

    return pd.DataFrame(dist_matrix, index=conditions, columns=conditions)


# ==============================================================================
# DRUG–GENE MATCH RANKING
# ==============================================================================

def rank_matches(dist_df, drug_conditions, gene_conditions):
    """For each drug, rank gene knockdowns by ascending Wasserstein distance."""
    rows = []
    for drug in sorted(drug_conditions):
        distances = dist_df.loc[drug, gene_conditions].sort_values()
        for rank, (gene, dist) in enumerate(distances.items(), 1):
            rows.append({
                "drug": drug,
                "rank": rank,
                "gene": gene,
                "distance": dist,
            })
    return pd.DataFrame(rows)


# ==============================================================================
# PERMUTATION TEST
# ==============================================================================

def permutation_test(df, features, drug_conditions, gene_conditions,
                     observed_dist_df, config):
    """Permutation test for top-K drug–gene matches."""
    rng = np.random.default_rng(config["random_seed"])
    top_k = config["permutation_top_k"]
    n_perm = config["n_permutations"]
    n_sub = config["permutation_subsample_n"]

    # Identify pairs to test
    pairs_to_test = []
    for drug in drug_conditions:
        top_genes = observed_dist_df.loc[drug, gene_conditions].nsmallest(top_k).index
        for gene in top_genes:
            pairs_to_test.append((drug, gene))

    observed = {(d, g): observed_dist_df.loc[d, g] for d, g in pairs_to_test}
    counts = {pair: 0 for pair in pairs_to_test}

    # Pool non-control cells
    non_ctrl = df[~df["is_control"]].copy()
    all_conditions = list(set(drug_conditions) | set(gene_conditions))
    condition_sizes = {c: (non_ctrl["condition_label"] == c).sum()
                       for c in all_conditions}

    print(f"Running {n_perm} permutations for {len(pairs_to_test)} drug–gene pairs ...")

    for perm_i in range(n_perm):
        # Shuffle condition labels
        shuffled = non_ctrl.copy()
        new_labels = np.concatenate([
            np.repeat(c, condition_sizes[c]) for c in all_conditions
        ])
        rng.shuffle(new_labels)
        shuffled["condition_label"] = new_labels

        # Subsample and compute distances for tested pairs only
        arrays = {}
        needed = set()
        for d, g in pairs_to_test:
            needed.add(d)
            needed.add(g)

        for cond in needed:
            X = shuffled.loc[shuffled["condition_label"] == cond, features].values
            if len(X) >= n_sub:
                idx = rng.choice(len(X), n_sub, replace=False)
            else:
                idx = rng.choice(len(X), n_sub, replace=True)
            arrays[cond] = X[idx].astype(np.float64)

        for drug, gene in pairs_to_test:
            X_a, X_b = arrays[drug], arrays[gene]
            M = ot.dist(X_a, X_b, metric="sqeuclidean")
            M_max = M.max()
            if M_max > 0:
                M /= M_max
            a = np.ones(n_sub) / n_sub
            b = np.ones(n_sub) / n_sub
            w = float(ot.sinkhorn2(a, b, M, reg=config["sinkhorn_reg"],
                                   numItermax=config["sinkhorn_max_iter"],
                                   warn=False))
            if w <= observed[(drug, gene)]:
                counts[(drug, gene)] += 1

        if (perm_i + 1) % 50 == 0 or perm_i == n_perm - 1:
            print(f"  permutation {perm_i + 1}/{n_perm}")

    # Compute p-values and FDR
    p_values = {pair: (counts[pair] + 1) / (n_perm + 1) for pair in pairs_to_test}
    pv_df = pd.DataFrame([
        {"drug": d, "gene": g, "p_value": p_values[(d, g)]}
        for d, g in pairs_to_test
    ])
    pv_df["p_adj"] = false_discovery_control(pv_df["p_value"].values, method="bh")
    pv_df["significant"] = pv_df["p_adj"] < config["alpha"]

    n_sig = pv_df["significant"].sum()
    print(f"Permutation test: {n_sig}/{len(pv_df)} pairs significant "
          f"(FDR < {config['alpha']})")

    return pv_df


def permutation_test_sscore(df, features, drug_conditions, gene_conditions,
                            observed_dist_df, config, metric="cosine"):
    """Permutation test for S-score distances.

    Shuffles condition labels at the cell level, recomputes S-scores from
    the shuffled data, then computes pairwise distances.  Tests whether
    observed drug–gene distances are smaller than expected by chance.
    """
    rng = np.random.default_rng(config["random_seed"])
    top_k = config["permutation_top_k"]
    n_perm = config["n_permutations"]
    dist_func = cosine if metric == "cosine" else euclidean

    # Identify pairs to test (top-K gene matches per drug)
    pairs_to_test = []
    for drug in drug_conditions:
        top_genes = observed_dist_df.loc[drug, gene_conditions].nsmallest(top_k).index
        for gene in top_genes:
            pairs_to_test.append((drug, gene))

    observed = {(d, g): observed_dist_df.loc[d, g] for d, g in pairs_to_test}
    counts = {pair: 0 for pair in pairs_to_test}

    # Work with non-control cells only
    non_ctrl = df[~df["is_control"]].copy()
    all_conditions = list(set(drug_conditions) | set(gene_conditions))
    condition_sizes = {c: (non_ctrl["condition_label"] == c).sum()
                       for c in all_conditions}

    # Pre-compute the EXPERIMENT-level structure: each EXPERIMENT maps to
    # a condition label.  We'll shuffle at the EXPERIMENT level to preserve
    # within-sample correlations.
    exp_to_cond = non_ctrl.drop_duplicates("EXPERIMENT").set_index("EXPERIMENT")["condition_label"]
    experiment_list = exp_to_cond.index.tolist()
    exp_condition_list = np.array(exp_to_cond.values.tolist())

    print(f"Running {n_perm} S-score permutations for "
          f"{len(pairs_to_test)} drug–gene pairs ...")

    needed = set()
    for d, g in pairs_to_test:
        needed.add(d)
        needed.add(g)

    for perm_i in range(n_perm):
        # Shuffle condition labels across experiments
        shuffled_labels = exp_condition_list.copy()
        rng.shuffle(shuffled_labels)
        perm_exp_to_cond = dict(zip(experiment_list, shuffled_labels))

        # Recompute S-scores with shuffled labels
        # Mean per experiment (unchanged), but group by shuffled condition
        exp_means = non_ctrl.groupby("EXPERIMENT")[features].mean()
        exp_cvs = non_ctrl.groupby("EXPERIMENT")[features].apply(
            lambda g: g.apply(_cv)
        )

        # Map experiments to shuffled conditions
        exp_cond_series = pd.Series(perm_exp_to_cond)

        # For each needed condition, compute its S-score profile
        profiles = {}
        for cond in needed:
            cond_exps = exp_cond_series[exp_cond_series == cond].index
            if len(cond_exps) == 0:
                continue
            mean_vec = exp_means.loc[cond_exps].mean().values
            cv_vec = exp_cvs.loc[cond_exps].mean().values
            profiles[cond] = np.concatenate([mean_vec, cv_vec])

        # Replace NaN
        for cond in profiles:
            profiles[cond] = np.nan_to_num(profiles[cond], nan=0.0)

        # Compute distances for tested pairs
        for drug, gene in pairs_to_test:
            if drug not in profiles or gene not in profiles:
                continue
            d = dist_func(profiles[drug], profiles[gene])
            if d <= observed[(drug, gene)]:
                counts[(drug, gene)] += 1

        if (perm_i + 1) % 200 == 0 or perm_i == n_perm - 1:
            print(f"  permutation {perm_i + 1}/{n_perm}")

    # Compute p-values and FDR
    p_values = {pair: (counts[pair] + 1) / (n_perm + 1) for pair in pairs_to_test}
    pv_df = pd.DataFrame([
        {"drug": d, "gene": g, "p_value": p_values[(d, g)]}
        for d, g in pairs_to_test
    ])
    pv_df["p_adj"] = false_discovery_control(pv_df["p_value"].values, method="bh")
    pv_df["significant"] = pv_df["p_adj"] < config["alpha"]

    n_sig = pv_df["significant"].sum()
    print(f"Permutation test: {n_sig}/{len(pv_df)} pairs significant "
          f"(FDR < {config['alpha']})")

    return pv_df


# ==============================================================================
# VISUALISATION
# ==============================================================================

def generate_plots(dist_df, match_table, drug_conditions, gene_conditions, config,
                    dist_label=None):
    """Generate heatmap, UMAP, and ranked-match bar chart."""
    if dist_label is None:
        dist_label = config["distance_mode"]
    fig_dir = Path(config["fig_dir"])
    fig_dir.mkdir(exist_ok=True)
    dpi = config["fig_dpi"]

    # --- 1. Distance heatmap with clustering ---
    print("Generating distance heatmap ...")
    cond_type = pd.Series("gene", index=dist_df.index)
    cond_type[cond_type.index.isin(drug_conditions)] = "drug"
    palette = {"gene": "#4C72B0", "drug": "#DD8452"}
    row_colors = cond_type.map(palette).rename("type")

    g = sns.clustermap(
        dist_df,
        method="ward",
        cmap="viridis_r",
        figsize=(12, 10),
        row_colors=row_colors,
        col_colors=row_colors,
        linewidths=0,
        xticklabels=True,
        yticklabels=True,
        dendrogram_ratio=(0.12, 0.12),
        cbar_pos=(0.02, 0.82, 0.03, 0.15),
    )
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.tick_params(labelsize=8)
    # Legend for condition type
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=palette["gene"], label="Gene KD"),
                       Patch(facecolor=palette["drug"], label="Drug")]
    g.ax_heatmap.legend(handles=legend_elements, loc="upper left",
                        bbox_to_anchor=(1.02, 1), frameon=False, fontsize=9)
    g.savefig(fig_dir / "ot_distance_heatmap.png", dpi=dpi, bbox_inches="tight")
    plt.close()

    # --- 2. UMAP of conditions ---
    print("Generating condition UMAP ...")
    n_conditions = len(dist_df)
    n_neighbors = min(5, n_conditions - 1)

    reducer = umap.UMAP(
        metric="precomputed",
        n_neighbors=n_neighbors,
        min_dist=0.3,
        random_state=config["random_seed"],
    )
    embedding = reducer.fit_transform(dist_df.values)

    fig, ax = plt.subplots(figsize=(10, 8))
    for ctype, color in palette.items():
        mask = cond_type == ctype
        ax.scatter(embedding[mask, 0], embedding[mask, 1],
                   c=color, s=80, label="Gene KD" if ctype == "gene" else "Drug",
                   edgecolors="white", linewidth=0.5, zorder=3)
    for i, label in enumerate(dist_df.index):
        ax.annotate(label, (embedding[i, 0], embedding[i, 1]),
                    fontsize=7, ha="center", va="bottom",
                    xytext=(0, 5), textcoords="offset points")
    ax.legend(frameon=False)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"Condition-level UMAP ({dist_label} distances)")
    sns.despine(ax=ax)
    fig.savefig(fig_dir / "ot_umap_conditions.png", dpi=dpi, bbox_inches="tight")
    plt.close()

    # --- 3. Top matches per drug ---
    print("Generating ranked match chart ...")
    top_n = min(5, len(gene_conditions))
    top_matches = match_table[match_table["rank"] <= top_n].copy()
    n_drugs = top_matches["drug"].nunique()

    fig, axes = plt.subplots(1, n_drugs, figsize=(4 * n_drugs, 5), sharey=False)
    if n_drugs == 1:
        axes = [axes]

    for ax, (drug, grp) in zip(axes, top_matches.groupby("drug")):
        grp = grp.sort_values("distance")
        colors = ["#DD8452"] * len(grp)
        if "significant" in grp.columns:
            colors = ["#55A868" if sig else "#CCCCCC"
                      for sig in grp["significant"]]
        ax.barh(grp["gene"], grp["distance"], color=colors, edgecolor="white")
        ax.set_xlabel(f"Distance ({dist_label})")
        ax.set_title(drug, fontweight="bold")
        ax.invert_yaxis()
        sns.despine(ax=ax)

    fig.suptitle("Top gene matches per drug", fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(fig_dir / "ot_ranked_matches.png", dpi=dpi, bbox_inches="tight")
    plt.close()

    print(f"Figures saved to {fig_dir}/")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 60)
    print("OT Morphological Profiling Pipeline")
    print("=" * 60)

    # 1. Load and parse
    df = load_and_parse_data(CONFIG)

    # 2. Select features
    features = list(CONFIG["shape_features"])
    if CONFIG["include_fluorescence"]:
        features.append(CONFIG["fluorescence_col"])

    # Drop rows with NaN in feature columns
    n_before = len(df)
    df = df.dropna(subset=features).copy()
    if len(df) < n_before:
        print(f"Dropped {n_before - len(df):,} rows with missing feature values")

    # 3. Normalise
    df = normalize_features(df, features, CONFIG)

    # 4. Identify non-control conditions
    non_ctrl = df[~df["is_control"]]
    drug_conditions = sorted(non_ctrl.loc[non_ctrl["is_drug"], "condition_label"].unique())
    gene_conditions = sorted(non_ctrl.loc[~non_ctrl["is_drug"], "condition_label"].unique())
    all_conditions = gene_conditions + drug_conditions

    print(f"\nDrug conditions ({len(drug_conditions)}): {drug_conditions}")
    print(f"Gene conditions ({len(gene_conditions)}): {gene_conditions}")

    # Cell counts per condition
    print("\nCells per condition:")
    for cond in all_conditions:
        n = (non_ctrl["condition_label"] == cond).sum()
        ctype = "drug" if cond in drug_conditions else "gene"
        print(f"  {cond:20s}  {n:>6,}  ({ctype})")

    # 5. Compute distances (mode-dependent)
    mode = CONFIG["distance_mode"]
    print(f"\nDistance mode: {mode}")

    if mode.startswith("sscore"):
        # S-score approach: condition-level profiles, pairwise distance
        profiles, score_cols = compute_s_scores(df, features, CONFIG)
        metric = "cosine" if mode == "sscore_cosine" else "euclidean"
        dist_df = compute_sscore_distances(profiles, score_cols, metric=metric)

        # Derive drug/gene lists from profiles
        drug_conditions = sorted(
            profiles[profiles["is_drug"]].index.tolist()
        )
        gene_conditions = sorted(
            profiles[~profiles["is_drug"]].index.tolist()
        )
        all_conditions = gene_conditions + drug_conditions

    else:
        # Cell-level OT approach
        # Optional PCA
        if CONFIG.get("pca_before_ot", False):
            threshold = CONFIG["pca_variance_threshold"]
            pca = PCA().fit(non_ctrl[features].values)
            cumvar = np.cumsum(pca.explained_variance_ratio_)
            n_components = int(np.searchsorted(cumvar, threshold) + 1)
            pca = PCA(n_components=n_components)
            pc_cols = [f"PC{i+1}" for i in range(n_components)]
            df[pc_cols] = pca.fit_transform(df[features].values)
            non_ctrl = df[~df["is_control"]]
            features = pc_cols
            print(f"PCA: {n_components} components explain "
                  f"{cumvar[n_components-1]:.1%} of variance")

        dist_df = compute_ot_distances(non_ctrl, features, all_conditions, CONFIG)

    dist_df.to_csv(CONFIG["distance_csv"])
    print(f"Distance matrix saved to {CONFIG['distance_csv']}")

    # 6. Rank drug–gene matches
    match_table = rank_matches(dist_df, drug_conditions, gene_conditions)

    # 7. Permutation test (if enabled)
    if CONFIG["n_permutations"] > 0:
        if mode == "cell_ot":
            pv_df = permutation_test(
                df, features, drug_conditions, gene_conditions, dist_df, CONFIG
            )
        else:
            metric = "cosine" if mode == "sscore_cosine" else "euclidean"
            pv_df = permutation_test_sscore(
                df, features, drug_conditions, gene_conditions, dist_df,
                CONFIG, metric=metric
            )
        match_table = match_table.merge(pv_df, on=["drug", "gene"], how="left")

    # 8. Save ranked matches
    match_table.to_csv(CONFIG["output_csv"], index=False)
    print(f"\nRanked matches saved to {CONFIG['output_csv']}")
    print("\nTop match per drug:")
    for _, row in match_table[match_table["rank"] == 1].iterrows():
        sig = ""
        if "p_adj" in row and pd.notna(row["p_adj"]):
            sig = f"  (p_adj={row['p_adj']:.4f})"
        print(f"  {row['drug']:15s} → {row['gene']:15s}  "
              f"dist={row['distance']:.4f}{sig}")

    # 9. Generate plots
    generate_plots(dist_df, match_table, drug_conditions, gene_conditions, CONFIG)

    print("\n" + "=" * 60)
    print("Pipeline complete.")
    print("=" * 60)


if __name__ == "__main__":
    main()
