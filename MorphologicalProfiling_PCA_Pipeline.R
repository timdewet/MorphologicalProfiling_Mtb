# ==============================================================================
# PCA + hierarchical clustering pipeline (generalised)
# - Loads a single dataset and parses sample names
# - Computes S-scores (mean + CV) vs NT controls for selected features
# - Reporter-specific fluorescence S-scores (one feature per reporter)
# - PCA for dimensionality reduction, Ward's D2 hierarchical clustering
# - k selected by silhouette + gap statistic
# - Two output modes:
#   A) Predictive: PCA trained on knockdown samples only, drug samples projected
#   B) Combined:   PCA trained on all samples, drugs identified separately
# - Robustness: pvclust + cell-level bootstrap
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, readr, stringr, ggplot2, tibble, tidyr
  library(ggrepel)     # geom_text_repel
  library(cluster)     # silhouette, clusGap, pam
  library(uwot)        # umap (display only)
  library(pvclust)     # bootstrap dendrogram p-values
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

# Include control samples in the analysis? If FALSE, controls are used for
# S-score baseline only and excluded from PCA/clustering
include_controls_in_analysis <- TRUE

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

# PCA variance threshold — retain PCs explaining this fraction of total variance
pca_variance_threshold <- 0.85

# Range of k values to evaluate for clustering
k_range <- 2:8

# Manual override for k (set to NULL to auto-select via silhouette + gap)
manual_k <- NULL

# Number of cell-level bootstrap iterations
n_bootstrap <- 200

# --- Figure output ---
fig_dir    <- "figures"
fig_dpi    <- 300
fig_width  <- 10    # inches
fig_height <- 7     # inches

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ----------------------------- helpers ----------------------------------------

# Helper to save base R plots (dendrograms, heatmaps, etc.) to PNG
save_base_plot <- function(filename, expr, width = fig_width, height = fig_height,
                           dpi = fig_dpi) {
  png(file.path(fig_dir, filename),
      width = width, height = height, units = "in", res = dpi)
  eval(expr)
  dev.off()
}

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

# Compute reporter-specific fluorescence S-scores and merge into s_values
add_fluorescence_features <- function(cell_data, s_values, fluorescence_col,
                                      cv_fn = cv) {
  if (!(fluorescence_col %in% names(cell_data))) {
    warning(sprintf("Fluorescence column '%s' not found in data. Skipping.",
                    fluorescence_col))
    return(s_values)
  }

  message("--- Computing reporter-specific fluorescence S-scores ---")
  message("  Each reporter is a separate feature (different biological pathways)")

  reporters <- sort(unique(cell_data$reporter))
  message(sprintf("  Reporters detected: %s", paste(reporters, collapse = ", ")))

  # Per-experiment mean and CV of fluorescence
  fluor_per_exp <- cell_data %>%
    group_by(EXPERIMENT, reporter, is_control) %>%
    summarise(fluor_mean = mean(.data[[fluorescence_col]], na.rm = TRUE),
              fluor_cv   = cv_fn(.data[[fluorescence_col]]),
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
  pooled_sigma_mean <- sd(fluor_per_exp$fluor_mean[fluor_per_exp$is_control],
                          na.rm = TRUE)
  pooled_sigma_cv   <- sd(fluor_per_exp$fluor_cv[fluor_per_exp$is_control],
                          na.rm = TRUE)

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
                    fluor_baseline$ctrl_mu_mean[r],
                    fluor_baseline$ctrl_sigma_mean[r]))
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

  message(sprintf("  Fluorescence features added: %s",
                  paste(fluor_cols, collapse = ", ")))
  message("  NAs imputed with 0 (= baseline, no reporter response)")

  s_values %>%
    left_join(fluor_wide, by = "EXPERIMENT")
}

# Build feature matrix from s_values, dropping metadata and NA columns
build_feature_matrix <- function(s_values, meta_cols) {
  mat <- s_values %>%
    dplyr::select(-all_of(meta_cols)) %>%
    as.matrix()
  rownames(mat) <- as.character(s_values$EXPERIMENT)

  # Drop columns that are all NA or contain any NA

  mat <- mat[, colSums(is.na(mat)) != nrow(mat), drop = FALSE]
  mat <- mat[, colSums(is.na(mat)) == 0, drop = FALSE]
  mat
}

# Run PCA, select components, return list with model + scores + diagnostics
run_pca <- function(mat, variance_threshold = 0.85, label = "") {
  pca <- prcomp(mat, center = TRUE, scale. = TRUE)
  cum_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  n_pcs <- which(cum_var >= variance_threshold)[1]

  # Fallback: at least 2 PCs for visualisation
  n_pcs <- max(2L, n_pcs)
  # Cap at n-1 (max meaningful PCs)
  n_pcs <- min(n_pcs, nrow(mat) - 1L, ncol(mat))

  pca_scores <- pca$x[, 1:n_pcs, drop = FALSE]

  message(sprintf("--- PCA %s---", if (label != "") paste0("(", label, ") ") else ""))
  message(sprintf("  Features: %d, Samples: %d", ncol(mat), nrow(mat)))
  message(sprintf("  PCs retained: %d (%.1f%% variance at threshold %.0f%%)",
                  n_pcs, 100 * cum_var[n_pcs], 100 * variance_threshold))
  for (i in 1:min(5, n_pcs)) {
    message(sprintf("    PC%d: %.1f%% (cumulative: %.1f%%)",
                    i, 100 * pca$sdev[i]^2 / sum(pca$sdev^2), 100 * cum_var[i]))
  }

  list(
    model     = pca,
    scores    = pca_scores,
    cum_var   = cum_var,
    n_pcs     = n_pcs,
    var_explained = pca$sdev^2 / sum(pca$sdev^2)
  )
}

# Select k using silhouette + gap statistic
select_k <- function(scores, k_range, dist_mat = NULL, manual_k = NULL) {
  if (!is.null(manual_k)) {
    message(sprintf("  Using manual k = %d", manual_k))
    return(manual_k)
  }

  if (is.null(dist_mat)) dist_mat <- dist(scores, method = "euclidean")

  hc <- hclust(dist_mat, method = "ward.D2")

  # Cap k_range at n-1
  k_range <- k_range[k_range < nrow(scores)]

  # Silhouette scan
  sil_widths <- sapply(k_range, function(k) {
    cl <- cutree(hc, k = k)
    n_actual <- length(unique(cl))
    if (n_actual >= 2 && n_actual < nrow(scores)) {
      mean(silhouette(cl, dist_mat)[, "sil_width"])
    } else {
      NA_real_
    }
  })
  names(sil_widths) <- k_range

  # Gap statistic
  gap <- clusGap(scores, FUNcluster = function(x, k) {
    list(cluster = cutree(hclust(dist(x), method = "ward.D2"), k = k))
  }, K.max = max(k_range), B = 500, verbose = FALSE)

  best_k_gap <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"])
  best_k_sil <- k_range[which.max(sil_widths)]

  message(sprintf("  Best k by silhouette: %d (width = %.3f)",
                  best_k_sil, max(sil_widths, na.rm = TRUE)))
  message(sprintf("  Best k by gap statistic: %d", best_k_gap))

  # Prefer silhouette; report if they disagree
  best_k <- best_k_sil
  if (best_k_gap != best_k_sil) {
    message(sprintf("  NOTE: silhouette and gap disagree (sil=%d, gap=%d). Using silhouette.",
                    best_k_sil, best_k_gap))
  }

  list(
    best_k     = best_k,
    sil_widths = sil_widths,
    gap        = gap,
    best_k_sil = best_k_sil,
    best_k_gap = best_k_gap
  )
}

# Cell-level bootstrap: resample cells, recompute S-scores -> PCA -> cluster
run_cell_bootstrap <- function(cell_data, variables, fluorescence_col,
                               include_fluorescence, feature_meta_cols,
                               include_controls, var_threshold,
                               best_k, n_boot = 200, exclude_drug = TRUE) {

  message(sprintf("--- Cell-level bootstrap (%d iterations) ---", n_boot))

  # Get sample names for the target matrix
  s0 <- compute_s_scores(cell_data, variables, control_col = "is_control")
  if (include_fluorescence) {
    s0 <- add_fluorescence_features(cell_data, s0, fluorescence_col, cv_fn = cv)
  }
  if (!include_controls) s0 <- s0 %>% filter(!is_control)
  if (exclude_drug) s0 <- s0 %>% filter(!is_drug)
  sample_names <- as.character(s0$EXPERIMENT)
  n_samples <- length(sample_names)

  boot_assignments <- matrix(NA_integer_, nrow = n_samples, ncol = n_boot)
  rownames(boot_assignments) <- sample_names

  for (b in seq_len(n_boot)) {
    if (b %% 50 == 0) message(sprintf("  Bootstrap iteration %d / %d", b, n_boot))

    # Resample cells within each experiment (with replacement)
    boot_data <- cell_data %>%
      group_by(EXPERIMENT) %>%
      slice_sample(prop = 1, replace = TRUE) %>%
      ungroup()

    # Recompute S-scores
    s_boot <- tryCatch(
      compute_s_scores(boot_data, variables, control_col = "is_control"),
      error = function(e) NULL
    )
    if (is.null(s_boot)) next

    if (include_fluorescence) {
      s_boot <- tryCatch(
        add_fluorescence_features(boot_data, s_boot, fluorescence_col, cv_fn = cv),
        error = function(e) s_boot
      )
    }

    if (!include_controls) s_boot <- s_boot %>% filter(!is_control)
    if (exclude_drug) s_boot <- s_boot %>% filter(!is_drug)

    mat_boot <- build_feature_matrix(s_boot, feature_meta_cols)

    # Ensure same samples in same order
    if (!all(sample_names %in% rownames(mat_boot))) next
    mat_boot <- mat_boot[sample_names, , drop = FALSE]

    # PCA + cluster
    pca_boot <- tryCatch(
      prcomp(mat_boot, center = TRUE, scale. = TRUE),
      error = function(e) NULL
    )
    if (is.null(pca_boot)) next

    cum_var_boot <- cumsum(pca_boot$sdev^2) / sum(pca_boot$sdev^2)
    n_pcs_boot <- max(2L, which(cum_var_boot >= var_threshold)[1])
    n_pcs_boot <- min(n_pcs_boot, nrow(mat_boot) - 1L, ncol(mat_boot))

    scores_boot <- pca_boot$x[, 1:n_pcs_boot, drop = FALSE]
    hc_boot <- hclust(dist(scores_boot), method = "ward.D2")
    boot_assignments[, b] <- cutree(hc_boot, k = best_k)
  }

  # Remove failed iterations
  valid <- colSums(!is.na(boot_assignments)) == n_samples
  boot_assignments <- boot_assignments[, valid, drop = FALSE]
  message(sprintf("  Successful iterations: %d / %d", sum(valid), n_boot))

  # Build co-assignment matrix
  n_valid <- ncol(boot_assignments)
  coassoc <- matrix(0, n_samples, n_samples)
  rownames(coassoc) <- colnames(coassoc) <- sample_names
  for (b in seq_len(n_valid)) {
    coassoc <- coassoc + outer(boot_assignments[, b], boot_assignments[, b], "==")
  }
  coassoc <- coassoc / n_valid

  list(
    coassignment = coassoc,
    assignments  = boot_assignments,
    n_valid      = n_valid
  )
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
    is_drug    = str_detect(experiment_type,
                            regex(drug_experiment_pattern, ignore_case = TRUE))
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
                n_distinct(cell_data$EXPERIMENT[!cell_data$is_drug &
                                                  !cell_data$is_control])))
message(sprintf("Experiment types: %s",
                paste(unique(cell_data$experiment_type), collapse = ", ")))
message(sprintf("Reporters: %s",
                paste(unique(cell_data$reporter), collapse = ", ")))

# ----------------------------- compute S-scores --------------------------------

s_values <- compute_s_scores(cell_data, variables_of_interest,
                             control_col = "is_control")

# ---- Reporter-specific fluorescence S-scores --------------------------------

if (include_fluorescence) {
  s_values <- add_fluorescence_features(cell_data, s_values, fluorescence_col,
                                        cv_fn = cv)
}

# ---- Build feature matrices -------------------------------------------------

meta_cols <- c("EXPERIMENT", "experiment_type", "reporter", "knockdown",
               "is_control", "is_drug")

if (include_controls_in_analysis) {
  s_values_analysis <- s_values
} else {
  s_values_analysis <- s_values %>% filter(!is_control)
  message("Excluding control samples from analysis (used for S-score baseline only)")
}

if (length(exclude_reporters) > 0) {
  n_excluded <- sum(s_values_analysis$reporter %in% exclude_reporters)
  s_values_analysis <- s_values_analysis %>% filter(!reporter %in% exclude_reporters)
  message(sprintf("Excluding reporter(s) from analysis: %s (%d samples removed)",
                  paste(exclude_reporters, collapse = ", "), n_excluded))
}

feature_mat <- build_feature_matrix(s_values_analysis, meta_cols)

# Split into knockdown (+ controls if included) and drug matrices
is_drug_sample   <- s_values_analysis$is_drug
feature_mat_kd   <- feature_mat[!is_drug_sample, , drop = FALSE]
feature_mat_drug <- feature_mat[is_drug_sample, , drop = FALSE]

samples_kd   <- rownames(feature_mat_kd)
samples_drug <- rownames(feature_mat_drug)

message(sprintf("\nFeature matrix: %d samples x %d features",
                nrow(feature_mat), ncol(feature_mat)))
message(sprintf("  Knockdown/control: %d samples", nrow(feature_mat_kd)))
message(sprintf("  Drug-treated: %d samples", nrow(feature_mat_drug)))

# ==============================================================================
# MODE A: Predictive — PCA trained on knockdowns, drugs projected in
# ==============================================================================

message("\n=== Mode A: Predictive (knockdowns train, drugs projected) ===")

# --- PCA on knockdown samples ---
pca_kd <- run_pca(feature_mat_kd, pca_variance_threshold, label = "Mode A")

# --- k selection ---
dist_kd <- dist(pca_kd$scores, method = "euclidean")
hc_kd   <- hclust(dist_kd, method = "ward.D2")

message("--- k selection (Mode A) ---")
k_result_kd <- select_k(pca_kd$scores, k_range, dist_kd, manual_k)

if (is.list(k_result_kd)) {
  best_k_kd <- k_result_kd$best_k
} else {
  best_k_kd <- k_result_kd
  k_result_kd <- list(best_k = best_k_kd)
}

clusters_kd <- cutree(hc_kd, k = best_k_kd)

# --- PAM validation ---
pam_kd <- pam(dist_kd, k = best_k_kd, diss = TRUE)
ari_kd <- {
  # Adjusted Rand Index (inline, no extra package)
  t <- table(clusters_kd, pam_kd$clustering)
  n <- sum(t)
  sum_comb_t <- sum(choose(t, 2))
  sum_comb_a <- sum(choose(rowSums(t), 2))
  sum_comb_b <- sum(choose(colSums(t), 2))
  expected <- sum_comb_a * sum_comb_b / choose(n, 2)
  max_index <- 0.5 * (sum_comb_a + sum_comb_b)
  if (max_index == expected) 1 else (sum_comb_t - expected) / (max_index - expected)
}
message(sprintf("  Ward vs PAM agreement (ARI): %.3f", ari_kd))

# --- Project drug samples into knockdown PCA space ---
drug_centered <- scale(feature_mat_drug,
                       center = pca_kd$model$center,
                       scale  = pca_kd$model$scale)
drug_pca <- drug_centered %*% pca_kd$model$rotation[, 1:pca_kd$n_pcs, drop = FALSE]

# Assign drugs to nearest cluster centroid
centroids_kd <- aggregate(pca_kd$scores,
                          by = list(cluster = clusters_kd), mean)
drug_cluster_assignments <- apply(drug_pca, 1, function(row) {
  dists <- apply(centroids_kd[, -1, drop = FALSE], 1, function(ctr) {
    sqrt(sum((row - ctr)^2))
  })
  centroids_kd$cluster[which.min(dists)]
})

# --- Mode A plot data ---
d_pca_kd <- tibble(
  PC1         = pca_kd$scores[, 1],
  PC2         = pca_kd$scores[, 2],
  EXPERIMENT  = samples_kd,
  cluster     = as.character(clusters_kd),
  sample_type = "Knockdown"
) %>%
  left_join(sample_info, by = "EXPERIMENT") %>%
  mutate(
    sample_type = if_else(is_control_label(knockdown), "Control", sample_type),
    cluster     = if_else(is_control_label(knockdown), "Control", cluster)
  )

d_pca_drug <- tibble(
  PC1         = drug_pca[, 1],
  PC2         = drug_pca[, 2],
  EXPERIMENT  = samples_drug,
  cluster     = paste0(drug_cluster_assignments, " (projected)"),
  sample_type = "Drug"
) %>%
  left_join(sample_info, by = "EXPERIMENT")

d_predictive <- bind_rows(d_pca_kd, d_pca_drug) %>%
  mutate(cluster = factor(cluster))

# ==============================================================================
# MODE A PLOTS
# ==============================================================================

# Scree plot
p_scree_kd <- {
  scree_df <- tibble(
    PC = seq_along(pca_kd$var_explained),
    Variance = pca_kd$var_explained,
    Cumulative = pca_kd$cum_var
  ) %>% filter(PC <= min(20, length(pca_kd$var_explained)))

  ggplot(scree_df) +
    geom_col(aes(x = PC, y = Variance), fill = "steelblue", alpha = 0.7) +
    geom_line(aes(x = PC, y = Cumulative), color = "red", linewidth = 1) +
    geom_point(aes(x = PC, y = Cumulative), color = "red") +
    geom_hline(yintercept = pca_variance_threshold, linetype = "dashed",
               color = "darkred") +
    geom_vline(xintercept = pca_kd$n_pcs + 0.5, linetype = "dotted",
               color = "grey40") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Principal Component", y = "Variance Explained",
         title = "Mode A: Scree plot (knockdown samples)") +
    theme_minimal()
}
print(p_scree_kd)
ggsave(file.path(fig_dir, "modeA_scree.png"), p_scree_kd,
       width = fig_width, height = fig_height, dpi = fig_dpi)

# k-selection diagnostic
if (is.list(k_result_kd) && !is.null(k_result_kd$sil_widths)) {
  p_k_kd <- {
    k_df <- tibble(k = as.integer(names(k_result_kd$sil_widths)),
                   silhouette = k_result_kd$sil_widths)
    ggplot(k_df, aes(x = k, y = silhouette)) +
      geom_line() + geom_point(size = 3) +
      geom_vline(xintercept = best_k_kd, linetype = "dashed", color = "red") +
      labs(x = "k (number of clusters)", y = "Mean silhouette width",
           title = "Mode A: Silhouette vs k") +
      theme_minimal()
  }
  print(p_k_kd)
  ggsave(file.path(fig_dir, "modeA_silhouette_vs_k.png"), p_k_kd,
         width = fig_width, height = fig_height, dpi = fig_dpi)

  # Gap statistic plot
  save_base_plot("modeA_gap_statistic.png", quote(
    plot(k_result_kd$gap, main = "Mode A: Gap statistic")
  ))
  print(plot(k_result_kd$gap, main = "Mode A: Gap statistic"))
}

# Dendrogram
save_base_plot("modeA_dendrogram.png", quote({
  plot(hc_kd, labels = s_values_analysis$knockdown[!is_drug_sample],
       main = "Mode A: Dendrogram (Ward's D2, knockdown samples)",
       xlab = "", sub = "")
  rect.hclust(hc_kd, k = best_k_kd, border = "red")
}))
plot(hc_kd, labels = s_values_analysis$knockdown[!is_drug_sample],
     main = "Mode A: Dendrogram (Ward's D2, knockdown samples)",
     xlab = "", sub = "")
rect.hclust(hc_kd, k = best_k_kd, border = "red")

# PCA biplot — coloured by cluster
pc1_var <- round(100 * pca_kd$var_explained[1], 1)
pc2_var <- round(100 * pca_kd$var_explained[2], 1)

p_pca_predictive <- ggplot(d_predictive, aes(PC1, PC2, color = cluster)) +
  geom_point(data = filter(d_predictive, sample_type != "Drug"),
             size = 3, alpha = 0.8) +
  geom_point(data = filter(d_predictive, sample_type == "Drug"),
             shape = 17, size = 4) +
  geom_text_repel(aes(label = knockdown), size = 3, max.overlaps = 20,
                  show.legend = FALSE) +
  labs(x = sprintf("PC1 (%.1f%%)", pc1_var),
       y = sprintf("PC2 (%.1f%%)", pc2_var),
       color = "Cluster",
       title = "Mode A: Drugs projected onto knockdown-trained PCA") +
  theme_minimal()

if (exists("theme_Publication")) {
  p_pca_predictive <- p_pca_predictive +
    theme_Publication() + theme(legend.position = "right")
}
print(p_pca_predictive)
ggsave(file.path(fig_dir, "modeA_pca_biplot.png"), p_pca_predictive,
       width = fig_width, height = fig_height, dpi = fig_dpi)

# PCA loadings — top features driving PC1 and PC2
loadings_df <- tibble(
  feature = rownames(pca_kd$model$rotation),
  PC1     = pca_kd$model$rotation[, 1],
  PC2     = pca_kd$model$rotation[, 2]
) %>%
  mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
  slice_max(magnitude, n = 15)

p_loadings_kd <- ggplot(loadings_df, aes(x = PC1, y = PC2)) +
  geom_segment(aes(xend = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")),
               color = "grey50") +
  geom_text_repel(aes(label = feature), size = 2.5, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = sprintf("PC1 loading (%.1f%%)", pc1_var),
       y = sprintf("PC2 loading (%.1f%%)", pc2_var),
       title = "Mode A: Top 15 feature loadings") +
  theme_minimal()
print(p_loadings_kd)
ggsave(file.path(fig_dir, "modeA_loadings.png"), p_loadings_kd,
       width = fig_width, height = fig_height, dpi = fig_dpi)

# Optional UMAP display (trained on knockdown PCA scores, drugs projected)
umap_nn_kd <- min(15L, nrow(pca_kd$scores) - 1L)
umap_kd <- uwot::umap(pca_kd$scores, n_neighbors = umap_nn_kd, min_dist = 0.2,
                       ret_model = TRUE)
umap_drug_proj <- uwot::umap_transform(drug_pca, umap_kd)

d_umap_pred <- bind_rows(
  tibble(UMAP1 = umap_kd$embedding[, 1], UMAP2 = umap_kd$embedding[, 2],
         EXPERIMENT = samples_kd, cluster = d_pca_kd$cluster,
         sample_type = d_pca_kd$sample_type) %>%
    left_join(sample_info, by = "EXPERIMENT"),
  tibble(UMAP1 = umap_drug_proj[, 1], UMAP2 = umap_drug_proj[, 2],
         EXPERIMENT = samples_drug, cluster = d_pca_drug$cluster,
         sample_type = "Drug") %>%
    left_join(sample_info, by = "EXPERIMENT")
) %>% mutate(cluster = factor(cluster))

p_umap_predictive <- ggplot(d_umap_pred, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(data = filter(d_umap_pred, sample_type != "Drug"),
             size = 3, alpha = 0.8) +
  geom_point(data = filter(d_umap_pred, sample_type == "Drug"),
             shape = 17, size = 4) +
  geom_text_repel(aes(label = knockdown), size = 3, max.overlaps = 20,
                  show.legend = FALSE) +
  labs(x = "UMAP1", y = "UMAP2", color = "Cluster",
       title = "Mode A: UMAP display (clustering from PCA, not UMAP)") +
  theme_minimal()

if (exists("theme_Publication")) {
  p_umap_predictive <- p_umap_predictive +
    theme_Publication() + theme(legend.position = "right")
}
print(p_umap_predictive)
ggsave(file.path(fig_dir, "modeA_umap.png"), p_umap_predictive,
       width = fig_width, height = fig_height, dpi = fig_dpi)

# ==============================================================================
# MODE B: Combined — PCA on all samples, drugs identified separately
# ==============================================================================

message("\n=== Mode B: Combined (all samples, drugs highlighted) ===")

# --- PCA on all samples ---
pca_all <- run_pca(feature_mat, pca_variance_threshold, label = "Mode B")

# --- k selection ---
dist_all <- dist(pca_all$scores, method = "euclidean")
hc_all   <- hclust(dist_all, method = "ward.D2")

message("--- k selection (Mode B) ---")
k_result_all <- select_k(pca_all$scores, k_range, dist_all, manual_k)

if (is.list(k_result_all)) {
  best_k_all <- k_result_all$best_k
} else {
  best_k_all <- k_result_all
  k_result_all <- list(best_k = best_k_all)
}

clusters_all <- cutree(hc_all, k = best_k_all)

# --- PAM validation ---
pam_all <- pam(dist_all, k = best_k_all, diss = TRUE)
ari_all <- {
  t <- table(clusters_all, pam_all$clustering)
  n <- sum(t)
  sum_comb_t <- sum(choose(t, 2))
  sum_comb_a <- sum(choose(rowSums(t), 2))
  sum_comb_b <- sum(choose(colSums(t), 2))
  expected <- sum_comb_a * sum_comb_b / choose(n, 2)
  max_index <- 0.5 * (sum_comb_a + sum_comb_b)
  if (max_index == expected) 1 else (sum_comb_t - expected) / (max_index - expected)
}
message(sprintf("  Ward vs PAM agreement (ARI): %.3f", ari_all))

# --- Mode B plot data ---
d_combined <- tibble(
  PC1        = pca_all$scores[, 1],
  PC2        = pca_all$scores[, 2],
  EXPERIMENT = rownames(feature_mat),
  cluster    = as.character(clusters_all)
) %>%
  left_join(sample_info, by = "EXPERIMENT") %>%
  left_join(s_values_analysis %>% dplyr::select(EXPERIMENT, is_drug, is_control),
            by = "EXPERIMENT") %>%
  mutate(
    sample_type = case_when(
      is_control ~ "Control",
      is_drug    ~ "Drug",
      TRUE       ~ "Knockdown"
    ),
    cluster = if_else(is_control, "Control", cluster),
    cluster = factor(cluster)
  )

# ==============================================================================
# MODE B PLOTS
# ==============================================================================

# Scree plot
p_scree_all <- {
  scree_df <- tibble(
    PC = seq_along(pca_all$var_explained),
    Variance = pca_all$var_explained,
    Cumulative = pca_all$cum_var
  ) %>% filter(PC <= min(20, length(pca_all$var_explained)))

  ggplot(scree_df) +
    geom_col(aes(x = PC, y = Variance), fill = "steelblue", alpha = 0.7) +
    geom_line(aes(x = PC, y = Cumulative), color = "red", linewidth = 1) +
    geom_point(aes(x = PC, y = Cumulative), color = "red") +
    geom_hline(yintercept = pca_variance_threshold, linetype = "dashed",
               color = "darkred") +
    geom_vline(xintercept = pca_all$n_pcs + 0.5, linetype = "dotted",
               color = "grey40") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Principal Component", y = "Variance Explained",
         title = "Mode B: Scree plot (all samples)") +
    theme_minimal()
}
print(p_scree_all)
ggsave(file.path(fig_dir, "modeB_scree.png"), p_scree_all,
       width = fig_width, height = fig_height, dpi = fig_dpi)

# k-selection diagnostic
if (is.list(k_result_all) && !is.null(k_result_all$sil_widths)) {
  p_k_all <- {
    k_df <- tibble(k = as.integer(names(k_result_all$sil_widths)),
                   silhouette = k_result_all$sil_widths)
    ggplot(k_df, aes(x = k, y = silhouette)) +
      geom_line() + geom_point(size = 3) +
      geom_vline(xintercept = best_k_all, linetype = "dashed", color = "red") +
      labs(x = "k (number of clusters)", y = "Mean silhouette width",
           title = "Mode B: Silhouette vs k") +
      theme_minimal()
  }
  print(p_k_all)
  ggsave(file.path(fig_dir, "modeB_silhouette_vs_k.png"), p_k_all,
         width = fig_width, height = fig_height, dpi = fig_dpi)

  save_base_plot("modeB_gap_statistic.png", quote(
    plot(k_result_all$gap, main = "Mode B: Gap statistic")
  ))
  print(plot(k_result_all$gap, main = "Mode B: Gap statistic"))
}

# Dendrogram
save_base_plot("modeB_dendrogram.png", quote({
  plot(hc_all, labels = s_values_analysis$knockdown,
       main = "Mode B: Dendrogram (Ward's D2, all samples)",
       xlab = "", sub = "")
  rect.hclust(hc_all, k = best_k_all, border = "red")
}))
all_labels <- s_values_analysis$knockdown
plot(hc_all, labels = all_labels,
     main = "Mode B: Dendrogram (Ward's D2, all samples)",
     xlab = "", sub = "")
rect.hclust(hc_all, k = best_k_all, border = "red")

# PCA biplot — coloured by cluster
pc1_var_all <- round(100 * pca_all$var_explained[1], 1)
pc2_var_all <- round(100 * pca_all$var_explained[2], 1)

p_pca_combined <- ggplot(d_combined, aes(PC1, PC2, color = cluster)) +
  geom_point(data = filter(d_combined, sample_type != "Drug"),
             size = 3, alpha = 0.8) +
  geom_point(data = filter(d_combined, sample_type == "Drug"),
             shape = 17, size = 4) +
  geom_text_repel(aes(label = knockdown), size = 3, max.overlaps = 20,
                  show.legend = FALSE) +
  labs(x = sprintf("PC1 (%.1f%%)", pc1_var_all),
       y = sprintf("PC2 (%.1f%%)", pc2_var_all),
       color = "Cluster",
       title = "Mode B: All samples (drugs = triangles)") +
  theme_minimal()

if (exists("theme_Publication")) {
  p_pca_combined <- p_pca_combined +
    theme_Publication() + theme(legend.position = "right")
}
print(p_pca_combined)
ggsave(file.path(fig_dir, "modeB_pca_biplot.png"), p_pca_combined,
       width = fig_width, height = fig_height, dpi = fig_dpi)

# PCA biplot — coloured by sample type
p_pca_combined_type <- ggplot(d_combined, aes(PC1, PC2, color = sample_type)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = knockdown), size = 3, max.overlaps = 20,
                  show.legend = FALSE) +
  labs(x = sprintf("PC1 (%.1f%%)", pc1_var_all),
       y = sprintf("PC2 (%.1f%%)", pc2_var_all),
       color = "Sample Type",
       title = "Mode B: All samples coloured by type") +
  theme_minimal()

if (exists("theme_Publication")) {
  p_pca_combined_type <- p_pca_combined_type +
    theme_Publication() + theme(legend.position = "right")
}
print(p_pca_combined_type)
ggsave(file.path(fig_dir, "modeB_pca_by_type.png"), p_pca_combined_type,
       width = fig_width, height = fig_height, dpi = fig_dpi)

# PCA loadings
loadings_all_df <- tibble(
  feature = rownames(pca_all$model$rotation),
  PC1     = pca_all$model$rotation[, 1],
  PC2     = pca_all$model$rotation[, 2]
) %>%
  mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
  slice_max(magnitude, n = 15)

p_loadings_all <- ggplot(loadings_all_df, aes(x = PC1, y = PC2)) +
  geom_segment(aes(xend = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")),
               color = "grey50") +
  geom_text_repel(aes(label = feature), size = 2.5, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = sprintf("PC1 loading (%.1f%%)", pc1_var_all),
       y = sprintf("PC2 loading (%.1f%%)", pc2_var_all),
       title = "Mode B: Top 15 feature loadings") +
  theme_minimal()
print(p_loadings_all)
ggsave(file.path(fig_dir, "modeB_loadings.png"), p_loadings_all,
       width = fig_width, height = fig_height, dpi = fig_dpi)

# Optional UMAP display
umap_nn_all <- min(15L, nrow(pca_all$scores) - 1L)
umap_all <- uwot::umap(pca_all$scores, n_neighbors = umap_nn_all, min_dist = 0.2)

d_umap_combined <- d_combined %>%
  mutate(UMAP1 = umap_all[, 1], UMAP2 = umap_all[, 2])

p_umap_combined <- ggplot(d_umap_combined, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(data = filter(d_umap_combined, sample_type != "Drug"),
             size = 3, alpha = 0.8) +
  geom_point(data = filter(d_umap_combined, sample_type == "Drug"),
             shape = 17, size = 4) +
  geom_text_repel(aes(label = knockdown), size = 3, max.overlaps = 20,
                  show.legend = FALSE) +
  labs(x = "UMAP1", y = "UMAP2", color = "Cluster",
       title = "Mode B: UMAP display (clustering from PCA, not UMAP)") +
  theme_minimal()

if (exists("theme_Publication")) {
  p_umap_combined <- p_umap_combined +
    theme_Publication() + theme(legend.position = "right")
}
print(p_umap_combined)
ggsave(file.path(fig_dir, "modeB_umap.png"), p_umap_combined,
       width = fig_width, height = fig_height, dpi = fig_dpi)

# ==============================================================================
# ROBUSTNESS ASSESSMENT
# ==============================================================================

message("\n=== Robustness assessment ===")

# --- 1. pvclust (bootstrap p-values on dendrogram nodes) ---
message("--- pvclust (Mode A: knockdowns) ---")
pv_kd <- pvclust(t(feature_mat_kd), method.hclust = "ward.D2",
                 method.dist = "euclidean", nboot = 1000, quiet = TRUE)
save_base_plot("modeA_pvclust.png", quote({
  plot(pv_kd, main = "Mode A: pvclust dendrogram (AU p-values)")
  pvrect(pv_kd, alpha = 0.95, pv = "au")
}))
plot(pv_kd, main = "Mode A: pvclust dendrogram (AU p-values)")
pvrect(pv_kd, alpha = 0.95, pv = "au")

message("--- pvclust (Mode B: all samples) ---")
pv_all <- pvclust(t(feature_mat), method.hclust = "ward.D2",
                  method.dist = "euclidean", nboot = 1000, quiet = TRUE)
save_base_plot("modeB_pvclust.png", quote({
  plot(pv_all, main = "Mode B: pvclust dendrogram (AU p-values)")
  pvrect(pv_all, alpha = 0.95, pv = "au")
}))
plot(pv_all, main = "Mode B: pvclust dendrogram (AU p-values)")
pvrect(pv_all, alpha = 0.95, pv = "au")

# --- 2. Cell-level bootstrap ---
message("--- Cell-level bootstrap (Mode A: knockdowns) ---")
cell_data_boot <- if (length(exclude_reporters) > 0) {
  cell_data %>% filter(!reporter %in% exclude_reporters)
} else {
  cell_data
}
boot_kd <- run_cell_bootstrap(
  cell_data          = cell_data_boot,
  variables          = variables_of_interest,
  fluorescence_col   = fluorescence_col,
  include_fluorescence = include_fluorescence,
  feature_meta_cols  = meta_cols,
  include_controls   = include_controls_in_analysis,
  var_threshold      = pca_variance_threshold,
  best_k             = best_k_kd,
  n_boot             = n_bootstrap,
  exclude_drug       = TRUE
)

# Co-assignment heatmap
if (boot_kd$n_valid > 0) {
  boot_labels <- s_values_analysis$knockdown[match(rownames(boot_kd$coassignment),
                                                    s_values_analysis$EXPERIMENT)]
  save_base_plot("modeA_bootstrap_heatmap.png", quote({
    heatmap(boot_kd$coassignment,
            labRow = boot_labels, labCol = boot_labels,
            main = sprintf("Mode A: Cell-level bootstrap co-assignment (%d iter)",
                           boot_kd$n_valid),
            col = colorRampPalette(c("white", "steelblue", "darkblue"))(50))
  }))
  heatmap(boot_kd$coassignment,
          labRow = boot_labels, labCol = boot_labels,
          main = sprintf("Mode A: Cell-level bootstrap co-assignment (%d iter)",
                         boot_kd$n_valid),
          col = colorRampPalette(c("white", "steelblue", "darkblue"))(50))
}

message("\n=== Pipeline complete ===")
