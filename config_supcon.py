"""
Shared Configuration for SupCon + OT Morphological Profiling Pipeline
=====================================================================
Imported by both extract_crops.py (Phase 1) and supcon_ot_pipeline.py (Phase 2).
"""

from pathlib import Path

# Base paths
_BASE = Path("/Users/timdewet/Library/CloudStorage/Dropbox-MMRU/Timothy de Wet"
             "/Science/Admin/2026")
_IMAGING_PIPELINE = _BASE / "Code" / "ImagingPipeline"
_MORPH_PROFILING = _BASE / "Code" / "MorphologicalProfiling_Mtb"
_MTB_LIBRARY = _BASE / "MtbLibrary"

CONFIG = {
    # ── Data sources ─────────────────────────────────────────────────────────
    "tiff_dirs": [
        str(_MTB_LIBRARY / "segmented_and_classified"),  # 47 TIFFs (QC-filtered)
    ],

    # ── Condition parsing ────────────────────────────────────────────────────
    # TIFF filenames follow: ExperimentType__Reporter__Knockdown.tif
    "name_separator": "__",
    "control_labels": ["NT", "No_drug"],
    "drug_experiment_pattern": "drug",  # case-insensitive match on experiment_type

    # Name corrections — applied after parsing filename stem.
    # Each dict matches on the full stem and overwrites specified fields.
    "name_corrections": [
        {"stem": "WT_Reporters_+_drug__imiB__Inn",     "reporter": "iniB", "knockdown": "INH"},
        {"stem": "WT_Reporters_+_drug__imiB__EMB",     "reporter": "iniB"},
        {"stem": "WT_Reporters_+_drug__imiB__No_drug", "reporter": "iniB"},
        {"stem": "WT_Reporters_+_drug__imiB__RIF",     "reporter": "iniB"},
        {"stem": "ATC_Strains__recA__dnaW2",           "knockdown": "dnaN1_rep2"},
    ],

    # ── Crop extraction (Phase 1) ───────────────────────────────────────────
    "crop_size": 96,
    "crop_pad": 10,  # pixels of padding around bounding box
    "imaging_pipeline_dir": str(_IMAGING_PIPELINE),
    "crop_output_h5": str(_MORPH_PROFILING / "crops" / "all_crops.h5"),
    "min_cells_per_condition": 50,

    # ── SupCon training (Phase 2) ────────────────────────────────────────────
    "supcon_input_channels": 2,   # phase + mask
    "embedding_dim": 512,         # ResNet-18 output before projection
    "projection_dim": 128,        # SupCon projection head output
    "backbone": "resnet18",
    "pretrained": True,
    "batch_size": 128,
    "epochs": 100,
    "lr": 1e-4,
    "weight_decay": 1e-4,
    "warmup_epochs": 5,           # linear LR warmup
    "temperature": 0.05,          # SupCon loss temperature
    "val_fraction": 0.15,
    "patience": 20,               # early stopping patience (epochs)
    "freeze_layers": 6,           # freeze conv1, bn1, relu, maxpool, layer1, layer2
    "exclude_controls_from_training": False,
    "max_samples_per_class": 5000,  # cap per class per epoch for balance

    # ── OT matching ──────────────────────────────────────────────────────────
    "sinkhorn_reg": 0.05,
    "sinkhorn_max_iter": 1000,
    "subsample_n": 500,
    "n_permutations": 10000,
    "permutation_top_k": 1,
    "permutation_subsample_n": 200,
    "pca_dims": 50,               # PCA reduction before OT (None to skip)
    "alpha": 0.05,
    "random_seed": 42,

    # ── Pathway mapping ───────────────────────────────────────────────────
    "gene_pathway_csv": str(_MORPH_PROFILING / "gene_pathway_map.csv"),

    # ── Output paths ─────────────────────────────────────────────────────────
    "checkpoint_dir": str(_MORPH_PROFILING / "checkpoints"),
    "embedding_dir": str(_MORPH_PROFILING / "embeddings"),
    "fig_dir": str(_MORPH_PROFILING / "figures"),
    "output_csv": str(_MORPH_PROFILING / "ot_supcon_ranked_matches.csv"),
    "pathway_output_csv": str(_MORPH_PROFILING / "ot_supcon_pathway_matches.csv"),
    "distance_csv": str(_MORPH_PROFILING / "ot_supcon_distance_matrix.csv"),
    "fig_dpi": 300,
}
