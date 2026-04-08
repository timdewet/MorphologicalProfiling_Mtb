"""
Shared Configuration for SupCon + OT Morphological Profiling Pipeline
=====================================================================
M. smegmatis branch.

Imported by both extract_crops.py (Phase 1) and supcon_ot_pipeline.py (Phase 2).
"""

from pathlib import Path

# Base paths
_BASE = Path("/Users/timdewet/Library/CloudStorage/Dropbox-MMRU/Timothy de Wet"
             "/Science/Admin/2026")
_IMAGING_PIPELINE = _BASE / "Code" / "ImagingPipeline"
_MORPH_PROFILING = _BASE / "Code" / "MorphologicalProfiling_Mtb"
_SMEG_DATA = Path("/Users/timdewet/Library/CloudStorage/Dropbox-MMRU"
                   "/Network Data/Microscopy/Tim/2026/MSM_Reanalysis")

CONFIG = {
    # ── Data sources ─────────────────────────────────────────────────────────
    "tiff_dirs": [
        str(_SMEG_DATA),  # 533 TIFFs (267 mutants, 20 drugs, 26 plasmid controls)
    ],

    # ── Condition parsing ────────────────────────────────────────────────────
    # All filenames start with "Labelled__" prefix, then:
    #   Mutants:  Labelled__MSMEG_XXXX_RY
    #   Drugs:    Labelled__Drugs__ZZZ_AX_RY  (concentration as fraction of MIC)
    #   Controls: Labelled__Plasmid_RY        (empty vector, 26 replicates)
    "control_labels": ["Plasmid"],  # empty vector controls (26 replicates)

    # Name corrections — applied after parsing filename stem.
    # Each dict matches on the full stem and overwrites specified fields.
    "name_corrections": [],  # TODO: add corrections as needed

    # ── Crop extraction (Phase 1) ───────────────────────────────────────────
    "crop_size": 128,  # larger than Mtb (96) — M. smegmatis cells are bigger
    "crop_pad": 10,    # pixels of padding around bounding box
    "imaging_pipeline_dir": str(_IMAGING_PIPELINE),
    "crop_output_h5": str(_MORPH_PROFILING / "output" / "crops" / "all_crops.h5"),
    "min_cells_per_condition": 50,

    # ── SupCon training (Phase 2) ────────────────────────────────────────────
    "supcon_input_channels": 3,   # phase + ParB fluorescence + mask
    "embedding_dim": 512,         # ResNet-18 output before projection
    "projection_dim": 128,        # SupCon projection head output
    "backbone": "resnet18",
    "pretrained": True,
    "batch_size": 256,
    "epochs": 100,
    "lr": 1e-4,
    "weight_decay": 1e-4,
    "warmup_epochs": 5,           # linear LR warmup
    "temperature": 0.05,          # SupCon loss temperature
    "val_fraction": 0.15,
    "patience": 15,               # early stopping patience (epochs)
    "freeze_layers": 0,           # 0 = no freezing
    "exclude_controls_from_training": False,
    "max_samples_per_class": 5000,  # cap per class per epoch for balance

    # ── OT matching ──────────────────────────────────────────────────────────
    "sinkhorn_reg": 0.05,
    "sinkhorn_max_iter": 1000,
    "subsample_n": 500,
    "n_permutations": 10000,
    "permutation_top_k": 1,
    "permutation_subsample_n": 200,
    "pca_dims": 50,              # PCA reduction before OT (None to skip)
    "alpha": 0.05,
    "random_seed": 42,

    # ── Gene annotation ─────────────────────────────────────────────────
    # Accession → gene name mapping from DetailedAll_UpdatedAnnotations.csv
    # Used for display labels in figures (accession is kept as internal ID)
    "annotation_csv": str(_MORPH_PROFILING / "input_data" / "DetailedAll_UpdatedAnnotations.csv"),
    "gene_pathway_csv": str(_MORPH_PROFILING / "gene_pathway_map.csv"),

    # ── Output paths ─────────────────────────────────────────────────────────
    "checkpoint_dir": str(_MORPH_PROFILING / "output" / "checkpoints"),
    "embedding_dir": str(_MORPH_PROFILING / "output" / "embeddings"),
    "fig_dir": str(_MORPH_PROFILING / "figures"),
    "output_csv": str(_MORPH_PROFILING / "output" / "ot_supcon_ranked_matches.csv"),
    "pathway_output_csv": str(_MORPH_PROFILING / "output" / "ot_supcon_pathway_matches.csv"),
    "distance_csv": str(_MORPH_PROFILING / "output" / "ot_supcon_distance_matrix.csv"),
    "fig_dpi": 300,
}
