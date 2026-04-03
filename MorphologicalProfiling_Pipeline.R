# ==============================================================================
# UMAP + consensus clustering pipeline (generalised)
# - Loads a single dataset and parses sample names
# - Computes S-scores (mean + CV) vs NT controls for selected features
# - Runs multiple UMAP+HDBSCAN fits, builds a co-association matrix
# - Hierarchical clustering on consensus distances
# - Two output modes:
#   A) Predictive: UMAP trained on knockdown samples only, drug samples projected
#   B) Combined:   UMAP trained on all samples, drugs identified separately
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, readr, stringr, ggplot2, tibble, tidyr
  library(dbscan)      # hdbscan
  library(uwot)        # umap
  library(ggrepel)     # geom_text_repel
  library(cluster)     # silhouette
  library(plotly)      # ggplotly (optional)
})

# Optional: custom theme (only if present)
if (file.exists("Scripts/Theme.R")) source("Scripts/Theme.R")

set.seed(10)

# ----------------------------- configuration ----------------------------------

input_files <- c(
  "input_data/data_extraction_18_03_25.csv",
  "input_data/all_morphology_combined.csv"
)

# Sample name format: ExperimentType__Reporter__Knockdown
# e.g. ATC_Strains__cydA__menH, WT_Reporters_+_drug__cydA__BDQ
name_separator <- "__"

# Which values in the knockdown field identify untreated/control samples?
# These are used to compute S-score baselines
control_labels <- c("NT", "No_drug")

# Helper: match control labels including numbered replicates (NT_1, NT_2, ...)
is_control_label <- function(x) {
  patt <- paste0("^(", paste(control_labels, collapse = "|"), ")(_.+)?$")
  str_detect(x, patt)
}

# Reporter backgrounds to exclude from analysis (empty vector = include all)
# e.g. c("cydA") to drop all samples with the cydA reporter background
exclude_reporters <- c("cydA")

# --- Sample name corrections ---
# Applied to parsed metadata AFTER parsing but BEFORE analysis.
# Each entry: list(experiment = "original_name", field = "corrected_value", ...)
# Supported fields: experiment_type, reporter, knockdown
name_corrections <- list(
  list(experiment = "WT_Reporters_+_drug__imiB__Inn", reporter = "iniB", knockdown = "INH"),
  list(experiment = "ATC_Strains__recA__dnaW2",       knockdown = "dnaN1_rep2")
)

# --- Replicate groups ---
# Define which experiments are replicates of the same biological condition.
# Each entry: canonical_name = c("experiment1", "experiment2", ...)
# When merge_replicates = TRUE, cells are pooled into one sample per group.
# When FALSE, replicates stay separate but are labelled _rep1, _rep2, etc.
replicate_groups <- list(
  dnaN = c("ATC_Strains__recA__dnaN1", "ATC_Strains__recA__dnaW2")
)
merge_replicates <- FALSE

# Include control samples in UMAP? If FALSE, controls are used for S-score
# baseline only and excluded from the embedding
include_controls_in_umap <- TRUE

# Include reporter fluorescence (INTENSITY.ch1.mean) as an additional feature?
# When TRUE, a per-reporter normalised fluorescence S-score is added to the
# feature matrix. Normalisation uses the control for the SAME reporter background.
include_fluorescence <- FALSE
fluorescence_col     <- "INTENSITY.ch1.mean"

# Which experiment_type values identify drug-treated samples?
# Samples matching this pattern are separated from knockdown samples
drug_experiment_pattern <- "drug"   # matched with str_detect (case-insensitive)

variables_of_interest <- c(
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
  "SHAPE.width.variation"
)

# Set to TRUE to run parameter selection diagnostics before the main analysis.
# Once you've chosen parameters, set to FALSE and fill in the values below.
run_parameter_selection <- FALSE

# Parameter grid for selection (only used when run_parameter_selection = TRUE)
# Actual values are capped to dataset size at runtime (see parameter selection section)
param_grid <- list(
  n_neighbors_fracs = c(0.15, 0.25, 0.4, 0.55, 0.75),  # as fraction of n samples
  minPts_fracs      = c(0.1, 0.2, 0.3, 0.45, 0.6),     # as fraction of n samples
  k_range           = 2:8                                # consensus k values to evaluate
)
n_umap_runs_param <- 20   # fewer runs during parameter search (faster)

# Final parameters (used for main analysis)
umap_params    <- list(n_neighbors = 3, min_dist = 0)
hdbscan_params <- list(minPts = 3)
n_umap_runs    <- 50
consensus_k    <- 4  # cutree() k

# --- Figure output ---
fig_dir    <- "figures"
fig_dpi    <- 300
fig_width  <- 10    # inches
fig_height <- 7     # inches

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ----------------------------- helpers ----------------------------------------

cv <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mu <- mean(x)
  if (!is.finite(mu) || mu == 0) return(NA_real_)
  sd(x) / mu
}

require_cols <- function(df, cols, df_name = deparse(substitute(df))) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "%s is missing required columns: %s",
      df_name, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
}

# Parse sample names into components
# Expected format: ExperimentType__Reporter__Knockdown
parse_sample_name <- function(name, sep = "__") {
  parts <- str_split_fixed(name, fixed(sep), n = 3)
  tibble(
    EXPERIMENT      = name,
    experiment_type = parts[, 1],
    reporter        = parts[, 2],
    knockdown       = parts[, 3]
  )
}

# Compute S-scores for mean and CV vs control samples:
# S = (value - mean(control_values)) / sd(control_values)
compute_s_scores <- function(df, variables, control_col = "is_control") {
  require_cols(df, c("EXPERIMENT", control_col), "df")
  require_cols(df, variables, "df")

  sample_means <- df %>%
    group_by(EXPERIMENT) %>%
    summarise(across(all_of(variables), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  sample_cvs <- df %>%
    group_by(EXPERIMENT) %>%
    summarise(across(all_of(variables), ~ cv(.x)), .groups = "drop")

  long_means <- sample_means %>%
    pivot_longer(-EXPERIMENT, names_to = "feature", values_to = "value") %>%
    mutate(stat = "mean")

  long_cvs <- sample_cvs %>%
    pivot_longer(-EXPERIMENT, names_to = "feature", values_to = "value") %>%
    mutate(stat = "cv")

  long_all <- bind_rows(long_means, long_cvs)

  # Identify control sample names
  control_names <- df %>%
    filter(.data[[control_col]]) %>%
    distinct(EXPERIMENT) %>%
    pull(EXPERIMENT)

  baseline <- long_all %>%
    filter(EXPERIMENT %in% control_names) %>%
    group_by(stat, feature) %>%
    summarise(mu = mean(value, na.rm = TRUE),
              sigma = sd(value, na.rm = TRUE),
              .groups = "drop")

  scored <- long_all %>%
    left_join(baseline, by = c("stat", "feature")) %>%
    mutate(s = (value - mu) / sigma,
           out_feature = case_when(
             stat == "mean" ~ paste0(feature, "_mean"),
             stat == "cv"   ~ paste0(feature, "_CV"),
             TRUE           ~ paste0(feature, "_", stat)
           )) %>%
    dplyr::select(EXPERIMENT, out_feature, s) %>%
    pivot_wider(names_from = out_feature, values_from = s)

  # Rejoin metadata
  df %>%
    distinct(EXPERIMENT, experiment_type, reporter, knockdown, is_control, is_drug) %>%
    left_join(scored, by = "EXPERIMENT")
}

make_umap_hdbscan <- function(mat, seed, n_neighbors, min_dist, minPts,
                              ret_model = FALSE) {
  set.seed(seed)
  emb <- uwot::umap(mat, n_neighbors = n_neighbors, min_dist = min_dist,
                     ret_model = ret_model)
  if (ret_model) {
    coords <- as.data.frame(emb$embedding)
  } else {
    coords <- as.data.frame(emb)
  }
  fit <- dbscan::hdbscan(coords[, 1:2, drop = FALSE], minPts = minPts)
  result <- list(embedding = coords, cluster = as.character(fit$cluster))
  if (ret_model) result$model <- emb
  result
}

consensus_from_clusters <- function(cluster_matrix) {
  n_runs <- ncol(cluster_matrix)
  n <- nrow(cluster_matrix)

  assoc_counts <- matrix(0, n, n)
  for (r in seq_len(n_runs)) {
    assoc_counts <- assoc_counts + outer(cluster_matrix[, r], cluster_matrix[, r], FUN = "==")
  }
  assoc_prop <- assoc_counts / n_runs
  as.dist(1 - assoc_prop)
}

# ----------------------------- load + parse data ------------------------------

cell_data <- bind_rows(lapply(input_files, readr::read_csv, show_col_types = FALSE))

# Parse sample names from EXPERIMENT column
sample_info <- cell_data %>%
  distinct(EXPERIMENT) %>%
  pull(EXPERIMENT) %>%
  map_dfr(~ parse_sample_name(.x, sep = name_separator))

# Apply name corrections
if (length(name_corrections) > 0) {
  message("--- Applying name corrections ---")
  for (correction in name_corrections) {
    idx <- which(sample_info$EXPERIMENT == correction$experiment)
    if (length(idx) == 1) {
      fields <- setdiff(names(correction), "experiment")
      for (f in fields) {
        old_val <- sample_info[[f]][idx]
        sample_info[[f]][idx] <- correction[[f]]
        message(sprintf("  Corrected %s: %s -> %s (sample: %s)",
                        f, old_val, correction[[f]], correction$experiment))
      }
    } else if (length(idx) == 0) {
      warning(sprintf("Name correction: experiment '%s' not found in data",
                      correction$experiment))
    }
  }

  # Rebuild EXPERIMENT names from corrected components and propagate to cell_data
  sample_info <- sample_info %>%
    mutate(EXPERIMENT_new = paste(experiment_type, reporter, knockdown,
                                  sep = name_separator))
  correction_map <- sample_info %>%
    filter(EXPERIMENT != EXPERIMENT_new) %>%
    dplyr::select(EXPERIMENT, EXPERIMENT_new)
  if (nrow(correction_map) > 0) {
    cell_data <- cell_data %>%
      mutate(EXPERIMENT = as.character(EXPERIMENT)) %>%
      left_join(correction_map, by = "EXPERIMENT") %>%
      mutate(EXPERIMENT = if_else(!is.na(EXPERIMENT_new),
                                  EXPERIMENT_new, EXPERIMENT)) %>%
      dplyr::select(-EXPERIMENT_new)
    message(sprintf("  Updated EXPERIMENT names for %d corrected sample(s)",
                    nrow(correction_map)))
  }
  sample_info <- sample_info %>%
    mutate(EXPERIMENT = EXPERIMENT_new) %>%
    dplyr::select(-EXPERIMENT_new)
}

# Apply replicate grouping
if (length(replicate_groups) > 0) {
  rep_lookup <- tibble(EXPERIMENT = character(), rep_group = character(),
                       rep_label = character(), merged_name = character())
  for (group_name in names(replicate_groups)) {
    exps <- replicate_groups[[group_name]]
    # Build a canonical merged EXPERIMENT name from the group members
    group_info <- sample_info %>% filter(EXPERIMENT %in% exps)
    exp_type <- group_info$experiment_type[1]
    reporters <- sort(unique(group_info$reporter))
    reporter_str <- paste(reporters, collapse = "+")
    canonical_name <- paste(exp_type, reporter_str, group_name, sep = name_separator)

    for (i in seq_along(exps)) {
      rep_lookup <- bind_rows(rep_lookup, tibble(
        EXPERIMENT  = exps[i],
        rep_group   = group_name,
        rep_label   = if (merge_replicates) group_name else paste0(group_name, "_rep", i),
        merged_name = canonical_name
      ))
    }
  }

  if (merge_replicates) {
    # Relabel EXPERIMENT in cell_data so grouped experiments share a canonical name
    cell_data <- cell_data %>%
      mutate(EXPERIMENT = as.character(EXPERIMENT)) %>%
      left_join(rep_lookup %>% dplyr::select(EXPERIMENT, merged_name), by = "EXPERIMENT") %>%
      mutate(EXPERIMENT = if_else(!is.na(merged_name), merged_name, EXPERIMENT)) %>%
      dplyr::select(-merged_name)
    # Re-parse sample_info for merged names
    sample_info <- cell_data %>%
      distinct(EXPERIMENT) %>%
      pull(EXPERIMENT) %>%
      map_dfr(~ parse_sample_name(.x, sep = name_separator))
    message(sprintf("Merged %d replicate group(s) into single samples",
                    length(replicate_groups)))
  } else {
    # Just relabel knockdown for clarity in plots
    sample_info <- sample_info %>%
      left_join(rep_lookup %>% dplyr::select(EXPERIMENT, rep_label), by = "EXPERIMENT") %>%
      mutate(knockdown = if_else(!is.na(rep_label), rep_label, knockdown)) %>%
      dplyr::select(-rep_label)
    message(sprintf("Replicate groups defined (%d) — kept as separate samples with _repN labels",
                    length(replicate_groups)))
  }
}

# Join parsed info back and flag controls / drug samples
cell_data <- cell_data %>%
  mutate(EXPERIMENT = as.character(EXPERIMENT)) %>%
  left_join(sample_info, by = "EXPERIMENT") %>%
  mutate(
    is_control = is_control_label(knockdown),
    is_drug    = str_detect(experiment_type, regex(drug_experiment_pattern, ignore_case = TRUE))
  )

# Print summary of parsed samples
message("--- Sample summary ---")
message(sprintf("Total samples: %d", n_distinct(cell_data$EXPERIMENT)))
message(sprintf("Control samples (knockdown in {%s}): %d",
                paste(control_labels, collapse = ", "),
                n_distinct(cell_data$EXPERIMENT[cell_data$is_control])))
message(sprintf("Drug-treated samples: %d",
                n_distinct(cell_data$EXPERIMENT[cell_data$is_drug])))
message(sprintf("Knockdown samples: %d",
                n_distinct(cell_data$EXPERIMENT[!cell_data$is_drug & !cell_data$is_control])))
message(sprintf("Experiment types: %s",
                paste(unique(cell_data$experiment_type), collapse = ", ")))
message(sprintf("Reporters: %s",
                paste(unique(cell_data$reporter), collapse = ", ")))

# ----------------------------- compute S-scores --------------------------------

s_values <- compute_s_scores(cell_data, variables_of_interest, control_col = "is_control")

# ---- Reporter-specific fluorescence S-scores --------------------------------
# INTENSITY.ch1.mean reflects the transcriptional reporter (recA/cydA/iniB).
# Each sample's fluorescence is normalised to the control for the SAME reporter,
# because baseline fluorescence differs between reporters.

if (include_fluorescence && fluorescence_col %in% names(cell_data)) {

  message("--- Computing reporter-specific fluorescence S-scores ---")
  message("  Each reporter is a separate feature (different biological pathways)")

  reporters <- sort(unique(cell_data$reporter))
  message(sprintf("  Reporters detected: %s", paste(reporters, collapse = ", ")))

  # Per-experiment mean and CV of fluorescence
  fluor_per_exp <- cell_data %>%
    group_by(EXPERIMENT, reporter, is_control) %>%
    summarise(fluor_mean = mean(.data[[fluorescence_col]], na.rm = TRUE),
              fluor_cv   = cv(.data[[fluorescence_col]]),
              .groups = "drop")

  # Baseline: mean + sd of control experiments within each reporter
  fluor_baseline <- fluor_per_exp %>%
    filter(is_control) %>%
    group_by(reporter) %>%
    summarise(
      ctrl_mu_mean    = mean(fluor_mean, na.rm = TRUE),
      ctrl_sigma_mean = sd(fluor_mean,   na.rm = TRUE),
      ctrl_mu_cv      = mean(fluor_cv,   na.rm = TRUE),
      ctrl_sigma_cv   = sd(fluor_cv,     na.rm = TRUE),
      n_ctrl          = n(),
      .groups = "drop"
    )

  # If only one control per reporter, sd is NA — fall back to pooled sd
  pooled_sigma_mean <- sd(fluor_per_exp$fluor_mean[fluor_per_exp$is_control], na.rm = TRUE)
  pooled_sigma_cv   <- sd(fluor_per_exp$fluor_cv[fluor_per_exp$is_control],   na.rm = TRUE)

  fluor_baseline <- fluor_baseline %>%
    mutate(
      ctrl_sigma_mean = if_else(is.na(ctrl_sigma_mean) | ctrl_sigma_mean == 0,
                                pooled_sigma_mean, ctrl_sigma_mean),
      ctrl_sigma_cv   = if_else(is.na(ctrl_sigma_cv) | ctrl_sigma_cv == 0,
                                pooled_sigma_cv, ctrl_sigma_cv)
    )
  message("  Reporter baselines (n controls):")
  for (r in seq_len(nrow(fluor_baseline))) {
    message(sprintf("    %s: n=%d, mu=%.2f, sigma=%.2f",
                    fluor_baseline$reporter[r], fluor_baseline$n_ctrl[r],
                    fluor_baseline$ctrl_mu_mean[r], fluor_baseline$ctrl_sigma_mean[r]))
  }

  # Compute S-scores and pivot to wide: one column per reporter
  fluor_scores_long <- fluor_per_exp %>%
    left_join(fluor_baseline %>% dplyr::select(-n_ctrl), by = "reporter") %>%
    mutate(
      s_mean = (fluor_mean - ctrl_mu_mean) / ctrl_sigma_mean,
      s_cv   = (fluor_cv   - ctrl_mu_cv)   / ctrl_sigma_cv
    ) %>%
    dplyr::select(EXPERIMENT, reporter, s_mean, s_cv)

  # Pivot: each reporter becomes its own pair of columns
  # e.g. cydA_fluor_mean, cydA_fluor_CV, iniB_fluor_mean, ...
  fluor_wide <- fluor_scores_long %>%
    pivot_wider(
      id_cols     = EXPERIMENT,
      names_from  = reporter,
      values_from = c(s_mean, s_cv),
      names_glue  = "{reporter}_fluor_{ifelse(.value == 's_mean', 'mean', 'CV')}"
    )

  # Each sample only has its own reporter measured — fill NAs with 0
  # (S-score of 0 = no change from control baseline, biologically appropriate)
  fluor_cols <- setdiff(names(fluor_wide), "EXPERIMENT")
  fluor_wide <- fluor_wide %>%
    mutate(across(all_of(fluor_cols), ~ replace_na(.x, 0)))

  message(sprintf("  Fluorescence features added: %s", paste(fluor_cols, collapse = ", ")))
  message("  NAs imputed with 0 (= baseline, no reporter response)")

  # Merge into s_values
  s_values <- s_values %>%
    left_join(fluor_wide, by = "EXPERIMENT")

} else if (include_fluorescence) {
  warning(sprintf("Fluorescence column '%s' not found in data. Skipping.", fluorescence_col))
}

# Build numeric matrix for UMAP — optionally exclude control samples
meta_cols <- c("EXPERIMENT", "experiment_type", "reporter", "knockdown", "is_control", "is_drug")

if (include_controls_in_umap) {
  s_values_umap <- s_values
} else {
  s_values_umap <- s_values %>% filter(!is_control)
  message("Excluding control samples from UMAP (used for S-score baseline only)")
}

if (length(exclude_reporters) > 0) {
  n_excluded <- sum(s_values_umap$reporter %in% exclude_reporters)
  s_values_umap <- s_values_umap %>% filter(!reporter %in% exclude_reporters)
  message(sprintf("Excluding reporter(s) from analysis: %s (%d samples removed)",
                  paste(exclude_reporters, collapse = ", "), n_excluded))
}

feature_mat <- s_values_umap %>%
  dplyr::select(-all_of(meta_cols)) %>%
  as.matrix()
rownames(feature_mat) <- as.character(s_values_umap$EXPERIMENT)

# Drop columns that are all NA or contain any NA
feature_mat <- feature_mat[, colSums(is.na(feature_mat)) != nrow(feature_mat), drop = FALSE]
feature_mat <- feature_mat[, colSums(is.na(feature_mat)) == 0, drop = FALSE]

# Split into knockdown (+ controls if included) and drug matrices
is_drug_sample <- s_values_umap$is_drug
feature_mat_kd   <- feature_mat[!is_drug_sample, , drop = FALSE]
feature_mat_drug <- feature_mat[is_drug_sample, , drop = FALSE]

samples_kd   <- rownames(feature_mat_kd)
samples_drug <- rownames(feature_mat_drug)

# ==============================================================================
# PARAMETER SELECTION (optional — set run_parameter_selection = TRUE)
# ==============================================================================
#
# Strategy:
#   1. Grid search over n_neighbors x minPts using the knockdown matrix
#   2. For each combination, run a reduced consensus clustering
#   3. Evaluate across a range of k values using:
#      - Mean silhouette width on the consensus distance matrix
#      - Consensus CDF area-under-curve (PAC: proportion of ambiguous clustering)
#        Lower PAC = cleaner consensus = better separation
#   4. Also report average HDBSCAN noise proportion and stability
#   5. Produce diagnostic plots to guide parameter choice

if (run_parameter_selection) {

  message("=== Running parameter selection ===")

  # Helper: proportion of ambiguous clustering (PAC)
  # Fraction of consensus values in (lower, upper) — lower is better
  compute_pac <- function(consensus_dist, lower = 0.1, upper = 0.9) {
    vals <- as.numeric(1 - as.matrix(consensus_dist))  # back to similarity
    vals <- vals[upper.tri(matrix(0, attr(consensus_dist, "Size"),
                                     attr(consensus_dist, "Size")))]
    mean(vals > lower & vals < upper)
  }

  # Compute actual integer values from fractions, capped to dataset size
  n_kd <- nrow(feature_mat_kd)
  nn_values  <- sort(unique(pmax(2L, pmin(n_kd - 1L,
                    as.integer(round(param_grid$n_neighbors_fracs * n_kd))))))
  mpt_values <- sort(unique(pmax(2L, pmin(n_kd - 1L,
                    as.integer(round(param_grid$minPts_fracs * n_kd))))))
  message(sprintf("  Dataset has %d knockdown samples", n_kd))
  message(sprintf("  n_neighbors grid: %s", paste(nn_values, collapse = ", ")))
  message(sprintf("  minPts grid:      %s", paste(mpt_values, collapse = ", ")))

  # Run grid search
  param_results <- list()
  grid <- expand.grid(
    n_neighbors = nn_values,
    minPts      = mpt_values,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(grid))) {
    nn  <- grid$n_neighbors[i]
    mpt <- grid$minPts[i]
    message(sprintf("  n_neighbors = %d, minPts = %d ...", nn, mpt))

    # Consensus runs
    runs_param <- map(1:n_umap_runs_param, ~ make_umap_hdbscan(
      mat = feature_mat_kd, seed = .x,
      n_neighbors = nn, min_dist = 0, minPts = mpt
    ))

    # HDBSCAN diagnostics (from individual runs)
    noise_props <- map_dbl(runs_param, ~ mean(.x$cluster == "0"))
    mean_noise  <- mean(noise_props)

    # Build consensus
    clust_mat_param <- do.call(cbind, lapply(runs_param, `[[`, "cluster"))
    dist_cons_param <- consensus_from_clusters(clust_mat_param)
    hc_param        <- hclust(dist_cons_param)

    pac <- compute_pac(dist_cons_param)

    # Evaluate silhouette across k range
    for (k in param_grid$k_range) {
      clust_k <- cutree(hc_param, k = k)

      # Silhouette requires >= 2 clusters actually present
      n_actual <- length(unique(clust_k))
      if (n_actual >= 2 && n_actual < nrow(feature_mat_kd)) {
        sil <- silhouette(clust_k, dist_cons_param)
        mean_sil <- mean(sil[, "sil_width"])
      } else {
        mean_sil <- NA_real_
      }

      param_results[[length(param_results) + 1]] <- tibble(
        n_neighbors = nn,
        minPts      = mpt,
        k           = k,
        mean_silhouette = mean_sil,
        pac             = pac,
        mean_noise_prop = mean_noise
      )
    }
  }

  param_results_df <- bind_rows(param_results)

  # --- Diagnostic plots ---

  # 1. Silhouette heatmap for each k (faceted), across n_neighbors x minPts
  p_sil <- ggplot(param_results_df,
                  aes(x = factor(n_neighbors), y = factor(minPts),
                      fill = mean_silhouette)) +
    geom_tile() +
    geom_text(aes(label = round(mean_silhouette, 2)), size = 2.5) +
    scale_fill_viridis_c(na.value = "grey80") +
    facet_wrap(~ paste("k =", k)) +
    labs(x = "n_neighbors", y = "minPts", fill = "Mean\nSilhouette",
         title = "Silhouette width across parameter grid") +
    theme_minimal()
  print(p_sil)
  ggsave(file.path(fig_dir, "param_silhouette_heatmap.png"), p_sil,
         width = fig_width, height = fig_height, dpi = fig_dpi)

  # 2. PAC across n_neighbors x minPts (independent of k)
  pac_df <- param_results_df %>%
    distinct(n_neighbors, minPts, pac)

  p_pac <- ggplot(pac_df,
                  aes(x = factor(n_neighbors), y = factor(minPts), fill = pac)) +
    geom_tile() +
    geom_text(aes(label = round(pac, 2)), size = 3) +
    scale_fill_viridis_c(direction = -1, na.value = "grey80") +
    labs(x = "n_neighbors", y = "minPts", fill = "PAC\n(lower = better)",
         title = "Proportion of Ambiguous Clustering (PAC)") +
    theme_minimal()
  print(p_pac)
  ggsave(file.path(fig_dir, "param_pac_heatmap.png"), p_pac,
         width = fig_width, height = fig_height, dpi = fig_dpi)

  # 3. Noise proportion across n_neighbors x minPts
  noise_df <- param_results_df %>%
    distinct(n_neighbors, minPts, mean_noise_prop)

  p_noise <- ggplot(noise_df,
                    aes(x = factor(n_neighbors), y = factor(minPts),
                        fill = mean_noise_prop)) +
    geom_tile() +
    geom_text(aes(label = round(mean_noise_prop, 2)), size = 3) +
    scale_fill_viridis_c(na.value = "grey80") +
    labs(x = "n_neighbors", y = "minPts", fill = "Mean\nNoise %",
         title = "HDBSCAN noise proportion") +
    theme_minimal()
  print(p_noise)
  ggsave(file.path(fig_dir, "param_noise_heatmap.png"), p_noise,
         width = fig_width, height = fig_height, dpi = fig_dpi)

  # 4. Silhouette vs k for each n_neighbors/minPts combo (line plot)
  p_sil_k <- ggplot(param_results_df,
                    aes(x = k, y = mean_silhouette,
                        color = factor(minPts))) +
    geom_line() +
    geom_point() +
    facet_wrap(~ paste("n_neighbors =", n_neighbors)) +
    labs(x = "k (number of clusters)", y = "Mean silhouette width",
         color = "minPts",
         title = "Silhouette vs k across parameter grid") +
    theme_minimal()
  print(p_sil_k)
  ggsave(file.path(fig_dir, "param_silhouette_vs_k.png"), p_sil_k,
         width = fig_width, height = fig_height, dpi = fig_dpi)

  # Print top parameter combinations
  best_params <- param_results_df %>%
    filter(!is.na(mean_silhouette)) %>%
    arrange(desc(mean_silhouette)) %>%
    head(10)
  message("\n--- Top 10 parameter combinations by silhouette ---")
  print(as.data.frame(best_params))

  message("\nUpdate umap_params, hdbscan_params, and consensus_k based on the above,")
  message("then set run_parameter_selection <- FALSE to run the main analysis.")
  stop("Parameter selection complete. Review diagnostics and update config.",
       call. = FALSE)
}

# ==============================================================================
# MODE A: Predictive — UMAP trained on knockdowns, drugs projected in
# ==============================================================================

message("=== Mode A: Predictive UMAP (knockdowns train, drugs projected) ===")
message(sprintf("Training on %d knockdown samples, projecting %d drug samples (%d features)",
                nrow(feature_mat_kd), nrow(feature_mat_drug), ncol(feature_mat_kd)))

# --- consensus clustering on knockdown samples only ---
runs_kd <- map(
  1:n_umap_runs,
  ~ make_umap_hdbscan(
    mat         = feature_mat_kd,
    seed        = .x,
    n_neighbors = umap_params$n_neighbors,
    min_dist    = umap_params$min_dist,
    minPts      = hdbscan_params$minPts
  )
)

cluster_mat_kd <- do.call(cbind, lapply(runs_kd, `[[`, "cluster"))
rownames(cluster_mat_kd) <- samples_kd
colnames(cluster_mat_kd) <- paste0("run_", seq_len(n_umap_runs))

dist_consensus_kd <- consensus_from_clusters(cluster_mat_kd)
hc_kd <- hclust(dist_consensus_kd)
clusters_consensus_kd <- cutree(hc_kd, k = consensus_k)

# --- final UMAP with model retained for projection ---
final_kd <- make_umap_hdbscan(
  mat         = feature_mat_kd,
  seed        = 10,
  n_neighbors = umap_params$n_neighbors,
  min_dist    = umap_params$min_dist,
  minPts      = hdbscan_params$minPts,
  ret_model   = TRUE
)

# Project drug samples into the knockdown-trained UMAP space
drug_emb <- uwot::umap_transform(feature_mat_drug, final_kd$model)
drug_emb <- as.data.frame(drug_emb)

# Build knockdown UMAP data
d_umap_kd <- final_kd$embedding %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(
    EXPERIMENT  = samples_kd,
    cluster     = as.character(clusters_consensus_kd),
    sample_type = "Knockdown"
  ) %>%
  left_join(sample_info, by = "EXPERIMENT") %>%
  mutate(
    sample_type = if_else(is_control_label(knockdown), "Control", sample_type),
    cluster     = if_else(is_control_label(knockdown), "Control", cluster)
  )

# Build projected drug UMAP data
d_umap_drug <- drug_emb %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(
    EXPERIMENT  = samples_drug,
    cluster     = "Drug",
    sample_type = "Drug"
  ) %>%
  left_join(sample_info, by = "EXPERIMENT")

d_umap_predictive <- bind_rows(d_umap_kd, d_umap_drug) %>%
  mutate(cluster = factor(cluster))

# --- Mode A plots ---
p_predictive <- ggplot(d_umap_predictive, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(data = filter(d_umap_predictive, sample_type == "Knockdown"), alpha = 0.6) +
  geom_point(data = filter(d_umap_predictive, sample_type == "Drug"),
             shape = 17, size = 3) +
  geom_text_repel(aes(label = knockdown), size = 3, max.overlaps = 20) +
  labs(x = "UMAP1", y = "UMAP2", color = "Cluster",
       title = "Mode A: Drugs projected onto knockdown-trained UMAP") +
  theme_minimal()

if (exists("theme_Publication")) {
  p_predictive <- p_predictive +
    theme_Publication() +
    theme(legend.position = "right")
}

p_predictive
ggsave(file.path(fig_dir, "modeA_umap_predictive.png"), p_predictive,
       width = fig_width, height = fig_height, dpi = fig_dpi)

# ==============================================================================
# MODE B: Combined — all samples train UMAP, drugs identified separately
# ==============================================================================

message("=== Mode B: Combined UMAP (all samples, drugs highlighted) ===")
message(sprintf("Training on all %d samples (%d features)",
                nrow(feature_mat), ncol(feature_mat)))

# --- consensus clustering on all samples ---
runs_all <- map(
  1:n_umap_runs,
  ~ make_umap_hdbscan(
    mat         = feature_mat,
    seed        = .x,
    n_neighbors = umap_params$n_neighbors,
    min_dist    = umap_params$min_dist,
    minPts      = hdbscan_params$minPts
  )
)

cluster_mat_all <- do.call(cbind, lapply(runs_all, `[[`, "cluster"))
rownames(cluster_mat_all) <- rownames(feature_mat)
colnames(cluster_mat_all) <- paste0("run_", seq_len(n_umap_runs))

dist_consensus_all <- consensus_from_clusters(cluster_mat_all)
hc_all <- hclust(dist_consensus_all)
clusters_consensus_all <- cutree(hc_all, k = consensus_k)

# --- final embedding ---
final_all <- make_umap_hdbscan(
  mat         = feature_mat,
  seed        = 10,
  n_neighbors = umap_params$n_neighbors,
  min_dist    = umap_params$min_dist,
  minPts      = hdbscan_params$minPts
)

d_umap_combined <- final_all$embedding %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(
    EXPERIMENT = rownames(feature_mat),
    cluster    = as.character(clusters_consensus_all)
  ) %>%
  left_join(sample_info, by = "EXPERIMENT") %>%
  left_join(s_values %>% dplyr::select(EXPERIMENT, is_drug, is_control), by = "EXPERIMENT") %>%
  mutate(
    sample_type = case_when(
      is_control ~ "Control",
      is_drug    ~ "Drug",
      TRUE       ~ "Knockdown"
    ),
    cluster = if_else(is_control, "Control", cluster),
    cluster = factor(cluster)
  )

# --- Mode B plots ---

# Coloured by cluster, drugs shown as triangles
p_combined <- ggplot(d_umap_combined, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(data = filter(d_umap_combined, sample_type != "Drug"), alpha = 0.6) +
  geom_point(data = filter(d_umap_combined, sample_type == "Drug"),
             shape = 17, size = 3) +
  geom_text_repel(aes(label = knockdown), size = 3, max.overlaps = 20) +
  labs(x = "UMAP1", y = "UMAP2", color = "Cluster",
       title = "Mode B: All samples (drugs = triangles)") +
  theme_minimal()

if (exists("theme_Publication")) {
  p_combined <- p_combined +
    theme_Publication() +
    theme(legend.position = "right")
}

p_combined
ggsave(file.path(fig_dir, "modeB_umap_combined.png"), p_combined,
       width = fig_width, height = fig_height, dpi = fig_dpi)

# Coloured by sample type
p_combined_type <- ggplot(d_umap_combined, aes(UMAP1, UMAP2, color = sample_type)) +
  geom_point() +
  geom_text_repel(aes(label = knockdown), size = 3, max.overlaps = 20) +
  labs(x = "UMAP1", y = "UMAP2", color = "Sample Type",
       title = "Mode B: All samples coloured by type") +
  theme_minimal()

if (exists("theme_Publication")) {
  p_combined_type <- p_combined_type +
    theme_Publication() +
    theme(legend.position = "right")
}

p_combined_type
ggsave(file.path(fig_dir, "modeB_umap_by_type.png"), p_combined_type,
       width = fig_width, height = fig_height, dpi = fig_dpi)

# plotly::ggplotly(p_combined)
