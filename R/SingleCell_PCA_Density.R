# ==============================================================================
# Single-cell PCA + UMAP + Diffusion Map with density heatmap
# - Runs PCA, UMAP, and diffusion map on globally z-scored single-cell features
# - Each cell plotted as a point, coloured by local density (yellow = dense,
#   blue = sparse) to reveal subpopulations
# - Faceted per condition to compare gene knockdowns vs drugs
# ==============================================================================

suppressPackageStartupMessages({
  library(MASS)       # kde2d for 2D kernel density estimation (load before tidyverse)
  library(tidyverse)  # dplyr::select masks MASS::select
  library(uwot)       # UMAP
  library(destiny)    # diffusion map
})

if (file.exists("R/Theme.R")) source("R/Theme.R")

set.seed(10)

# ----------------------------- configuration ----------------------------------

input_files <- c(
  "input_data/smeg_morphology_data.csv"  # TODO: update with actual filename
)

control_labels <- c("Plasmid")  # empty vector controls

# Helper: match control labels including numbered replicates
is_control_label <- function(x) {
  patt <- paste0("^(", paste(control_labels, collapse = "|"), ")(_.+)?$")
  str_detect(x, patt)
}

exclude_reporters <- c()
name_corrections <- list()
replicate_groups <- list()
merge_replicates <- TRUE

drug_experiment_pattern <- "drug"

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

# Outlier filtering — remove cells beyond this percentile on any feature
# (catches segmentation artefacts that distort PCA axes)
outlier_quantile <- 0.99

# UMAP parameters
umap_n_neighbors <- 30     # number of nearest neighbours
umap_min_dist    <- 0.3      # minimum distance between embedded points (0 = tighter clusters)
umap_metric      <- "euclidean"

# Diffusion map parameters
dm_k             <- 30     # number of nearest neighbours for transition kernel
dm_n_eigs        <- 10     # number of diffusion components to compute

# Density heatmap settings
kde_n           <- 200     # grid resolution for kernel density estimation
point_size      <- 3     # size of scatter points
point_alpha     <- 0.6     # transparency for scatter points
n_subsample     <- NULL    # cells per condition for plotting (NULL = all)

# Figure output
fig_dir    <- "figures"
fig_dpi    <- 300
fig_width  <- 14
fig_height <- 7

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ----------------------------- helpers ----------------------------------------

parse_sample_name <- function(name, sep = NULL) {
  stripped <- str_replace(name, "^Labelled__Drugs__", "")
  stripped <- str_replace(stripped, "^Labelled__", "")

  ctrl_match   <- str_match(stripped, "^(Plasmid)_R(\\d+)$")
  mutant_match <- str_match(stripped, "^(MSMEG_\\d+)_R(\\d+)$")
  drug_match   <- str_match(stripped, "^([A-Z][A-Z0-9]{1,4})_(\\d+X)_R(\\d+)$")

  if (!is.na(ctrl_match[1, 1])) {
    tibble(EXPERIMENT = name, experiment_type = "control",
           reporter = "", knockdown = "Plasmid")
  } else if (!is.na(mutant_match[1, 1])) {
    tibble(EXPERIMENT = name, experiment_type = "mutant",
           reporter = "", knockdown = mutant_match[1, 2])
  } else if (!is.na(drug_match[1, 1])) {
    tibble(EXPERIMENT = name, experiment_type = "drug",
           reporter = "",
           knockdown = paste0(drug_match[1, 2], "_", drug_match[1, 3]))
  } else {
    warning("Cannot parse sample name: ", name)
    tibble(EXPERIMENT = name, experiment_type = NA_character_,
           reporter = NA_character_, knockdown = NA_character_)
  }
}

# Estimate local density at each point using MASS::kde2d
# Returns a numeric vector of log-density values rescaled to [0, 1]
estimate_density <- function(x, y, n = kde_n) {
  dens <- kde2d(x, y, n = n,
                lims = c(range(x) + c(-1, 1), range(y) + c(-1, 1)))
  # Interpolate grid density back to each point
  ix <- findInterval(x, dens$x, all.inside = TRUE)
  iy <- findInterval(y, dens$y, all.inside = TRUE)
  d <- mapply(function(i, j) dens$z[i, j], ix, iy)
  # Log-transform and rescale to [0, 1] — spreads the colour scale
  # so dense vs sparse regions are distinguishable in any embedding
  d <- log1p(d)
  (d - min(d)) / (max(d) - min(d) + .Machine$double.eps)
}

# ----------------------------- load & parse -----------------------------------

cat("Loading data...\n")
raw <- bind_rows(lapply(input_files, read_csv, show_col_types = FALSE))

missing <- setdiff(variables_of_interest, names(raw))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

meta <- map_dfr(unique(raw$EXPERIMENT), parse_sample_name)

for (fix in name_corrections) {
  idx <- meta$EXPERIMENT == fix$experiment
  if (any(idx)) {
    for (field in setdiff(names(fix), "experiment")) {
      meta[[field]][idx] <- fix[[field]]
    }
  }
}

if (!merge_replicates) {
  for (grp_name in names(replicate_groups)) {
    members <- replicate_groups[[grp_name]]
    for (i in seq_along(members)) {
      idx <- meta$EXPERIMENT == members[i]
      if (any(idx)) {
        meta$knockdown[idx] <- paste0(grp_name, "_rep", i)
      }
    }
  }
}

df <- raw %>%
  dplyr::select(EXPERIMENT, all_of(variables_of_interest),
                any_of("INTENSITY.ch1.mean")) %>%
  left_join(meta, by = "EXPERIMENT") %>%
  mutate(
    is_control = is_control_label(knockdown),
    is_drug    = str_detect(experiment_type, regex(drug_experiment_pattern, ignore_case = TRUE))
  )

# Exclude reporters
if (length(exclude_reporters) > 0) {
  df <- df %>% filter(!reporter %in% exclude_reporters)
}

cat(sprintf("Loaded %d cells across %d conditions\n", nrow(df), n_distinct(df$EXPERIMENT)))

# ----------------------------- outlier filtering ------------------------------

cat("Filtering outliers...\n")

# Remove rows with NA
df <- df[complete.cases(df[, variables_of_interest]), ]

# Remove cells beyond the outlier_quantile on any feature
n_before <- nrow(df)
upper <- sapply(df[, variables_of_interest], quantile, probs = outlier_quantile, na.rm = TRUE)
lower <- sapply(df[, variables_of_interest], quantile, probs = 1 - outlier_quantile, na.rm = TRUE)
keep <- apply(df[, variables_of_interest], 1, function(row) {
  all(row <= upper & row >= lower)
})
df <- df[keep, ]
cat(sprintf("  Removed %d outlier cells (%.1f%%), %d remain\n",
            n_before - nrow(df), 100 * (n_before - nrow(df)) / n_before, nrow(df)))

# ----------------------------- PCA --------------------------------------------

cat("Running PCA on single cells...\n")

feat_mat <- df %>%
  dplyr::select(all_of(variables_of_interest))

# Drop zero/constant-variance columns (e.g. SHAPE.pole is often constant)
col_vars <- apply(feat_mat, 2, var, na.rm = TRUE)
zero_var <- names(col_vars)[col_vars == 0 | is.na(col_vars)]
if (length(zero_var) > 0) {
  cat(sprintf("  Dropping %d constant feature(s): %s\n",
              length(zero_var), paste(zero_var, collapse = ", ")))
  feat_mat <- feat_mat[, setdiff(names(feat_mat), zero_var)]
}

# Global z-score and PCA
pca_result <- prcomp(feat_mat, center = TRUE, scale. = TRUE)

# Variance explained
var_pct <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# Attach PC scores to data
df <- df %>%
  mutate(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2]
  )

# Label for condition type
df <- df %>%
  mutate(
    condition_type = case_when(
      is_control ~ "Control",
      is_drug    ~ "Drug",
      TRUE       ~ "Knockdown"
    )
  )

cat(sprintf("PC1: %.1f%% variance, PC2: %.1f%% variance\n", var_pct[1], var_pct[2]))

# ----------------------------- UMAP ------------------------------------------

cat("Running UMAP on single cells...\n")

umap_result <- umap(
  as.matrix(feat_mat),
  n_neighbors = umap_n_neighbors,
  min_dist    = umap_min_dist,
  metric      = umap_metric,
  scale       = "Z",
  n_threads   = parallel::detectCores() - 1,
  ret_model   = FALSE
)

df <- df %>%
  mutate(
    UMAP1 = umap_result[, 1],
    UMAP2 = umap_result[, 2]
  )

cat("UMAP complete.\n")

# # ----------------------------- Diffusion Map ----------------------------------
# 
# cat("Running diffusion map on single cells...\n")
# 
# # Scale features for diffusion map (same z-scoring as PCA)
# feat_scaled <- scale(as.matrix(feat_mat))
# 
# dm_result <- DiffusionMap(feat_scaled, k = dm_k, n_eigs = dm_n_eigs)
# 
# df <- df %>%
#   mutate(
#     DC1 = eigenvectors(dm_result)[, 1],
#     DC2 = eigenvectors(dm_result)[, 2]
#   )
# 
# cat("Diffusion map complete.\n")

# ----------------------------- subsample for plotting -------------------------

if (!is.null(n_subsample)) {
  set.seed(42)
  df_plot <- df %>%
    group_by(knockdown) %>%
    slice_sample(prop = 1) %>%
    slice_head(n = n_subsample) %>%
    ungroup()
  cat(sprintf("Subsampled to %d cells for plotting\n", nrow(df_plot)))
} else {
  df_plot <- df
}

# ----------------------------- shared plotting functions -----------------------

apply_theme <- function(p, legend_pos = "right") {
  if (exists("theme_Publication")) {
    p + theme_Publication() +
      theme(legend.position = legend_pos,
            strip.text = element_text(face = "bold", size = 9))
  } else {
    p + theme_minimal(base_size = 13) +
      theme(legend.position = legend_pos,
            strip.text = element_text(face = "bold", size = 9))
  }
}

# Generate all 4 density plots for a given embedding
generate_density_plots <- function(df_plot, x_col, y_col, xlab, ylab,
                                   method_label, file_prefix) {

  cat(sprintf("--- Generating %s density plots ---\n", method_label))

  # Estimate global density
  cat("  Estimating global density...\n")
  df_plot$density <- estimate_density(df_plot[[x_col]], df_plot[[y_col]])

  # Estimate per-condition density
  cat("  Estimating per-condition densities...\n")
  df_plot <- df_plot %>%
    group_by(knockdown) %>%
    mutate(density_local = estimate_density(.data[[x_col]], .data[[y_col]])) %>%
    ungroup()

  df_plot <- df_plot %>%
    mutate(
      panel_label = paste0(knockdown, " (", condition_type, ")"),
      panel_label = factor(panel_label)
    )

  axis_range_x <- range(df_plot[[x_col]])
  axis_range_y <- range(df_plot[[y_col]])

  # --- Plot A: all cells ---
  cat("  All cells plot...\n")
  p_all <- ggplot(df_plot %>% arrange(density),
                  aes(x = .data[[x_col]], y = .data[[y_col]], colour = density)) +
    geom_point(size = point_size, alpha = point_alpha, shape = 16) +
    scale_colour_viridis_c(option = "D", name = "Density") +
    labs(x = xlab, y = ylab, title = sprintf("Single-Cell %s — All Cells", method_label))
  p_all <- apply_theme(p_all)
  ggsave(file.path(fig_dir, paste0(file_prefix, "_density_all.pdf")),
         plot = p_all, width = fig_width, height = fig_height, dpi = fig_dpi)

  # --- Plot B: faceted by condition ---
  cat("  Faceted plot...\n")
  p_facet <- ggplot(df_plot %>% arrange(density_local),
                    aes(x = .data[[x_col]], y = .data[[y_col]], colour = density_local)) +
    geom_point(size = point_size * 0.7, alpha = point_alpha, shape = 16) +
    scale_colour_viridis_c(option = "D", name = "Density") +
    facet_wrap(~ panel_label, ncol = 9) +
    coord_cartesian(xlim = axis_range_x, ylim = axis_range_y) +
    labs(x = xlab, y = ylab,
         title = sprintf("Single-Cell %s — Per-Condition Density", method_label))
  n_panels <- n_distinct(df_plot$panel_label)
  facet_h  <- 2 + ceiling(n_panels / 4) * 3
  p_facet <- apply_theme(p_facet)
  ggsave(file.path(fig_dir, paste0(file_prefix, "_density_faceted.pdf")),
         plot = p_facet, width = 16, height = 8, dpi = fig_dpi)

  # --- Plot C: knockdowns + NT controls ---
  cat("  Knockdowns plot...\n")
  df_kd <- df_plot %>% filter(condition_type %in% c("Knockdown", "Control"))
  p_kd <- ggplot(df_kd %>% arrange(density_local),
                 aes(x = .data[[x_col]], y = .data[[y_col]], colour = density_local)) +
    geom_point(size = point_size, alpha = point_alpha, shape = 16) +
    scale_colour_viridis_c(option = "D", name = "Density") +
    facet_wrap(~ knockdown, ncol = 4) +
    coord_cartesian(xlim = axis_range_x, ylim = axis_range_y) +
    labs(x = xlab, y = ylab,
         title = sprintf("Single-Cell %s — Gene Knockdowns", method_label))
  kd_h <- 2 + ceiling(n_distinct(df_kd$knockdown) / 4) * 3
  p_kd <- apply_theme(p_kd)
  ggsave(file.path(fig_dir, paste0(file_prefix, "_density_knockdowns.png")),
         plot = p_kd, width = fig_width + 4, height = kd_h, dpi = fig_dpi)

  # --- Plot D: drugs only ---
  cat("  Drugs plot...\n")
  df_drug <- df_plot %>% filter(condition_type == "Drug")
  p_drug <- ggplot(df_drug %>% arrange(density_local),
                   aes(x = .data[[x_col]], y = .data[[y_col]], colour = density_local)) +
    geom_point(size = point_size, alpha = point_alpha, shape = 16) +
    scale_colour_viridis_c(option = "D", name = "Density") +
    facet_wrap(~ knockdown, ncol = 4) +
    coord_cartesian(xlim = axis_range_x, ylim = axis_range_y) +
    labs(x = xlab, y = ylab,
         title = sprintf("Single-Cell %s — Drug Treatments", method_label))
  drug_h <- 2 + ceiling(n_distinct(df_drug$knockdown) / 4) * 3
  p_drug <- apply_theme(p_drug)
  ggsave(file.path(fig_dir, paste0(file_prefix, "_density_drugs.png")),
         plot = p_drug, width = fig_width + 4, height = drug_h, dpi = fig_dpi)

  cat(sprintf("  %s plots saved.\n", method_label))
}

# ----------------------------- generate PCA plots -----------------------------

generate_density_plots(
  df_plot,
  x_col        = "PC1",
  y_col        = "PC2",
  xlab         = sprintf("PC1 (%.1f%%)", var_pct[1]),
  ylab         = sprintf("PC2 (%.1f%%)", var_pct[2]),
  method_label = "PCA",
  file_prefix  = "sc_pca"
)

# ----------------------------- generate UMAP plots ----------------------------

generate_density_plots(
  df_plot,
  x_col        = "UMAP1",
  y_col        = "UMAP2",
  xlab         = "UMAP 1",
  ylab         = "UMAP 2",
  method_label = "UMAP",
  file_prefix  = "sc_umap"
)

# ----------------------------- generate Diffusion Map plots -------------------

# generate_density_plots(
#   df_plot,
#   x_col        = "DC1",
#   y_col        = "DC2",
#   xlab         = "DC 1",
#   ylab         = "DC 2",
#   method_label = "Diffusion Map",
#   file_prefix  = "sc_dm"
# )

# ----------------------------- UMAP + fluorescence (recA) ---------------------

cat("Generating UMAP fluorescence overlay for recA strains...\n")

if ("INTENSITY.ch1.mean" %in% names(df_plot)) {

  df_recA <- df_plot %>%
    filter(reporter == "recA") %>%
    filter(!is.na(INTENSITY.ch1.mean))

  if (nrow(df_recA) > 0) {

    # Background subtract using 95th percentile of NT control intensity
    # This ensures ~95% of NT cells fall at or below zero (dark blue)
    bg_threshold <- df_recA %>%
      filter(is_control_label(knockdown), !is_drug) %>%
      pull(INTENSITY.ch1.mean) %>%
      quantile(0.95, na.rm = TRUE)
    cat(sprintf("  NT background (95th percentile): %.2f\n", bg_threshold))

    df_recA <- df_recA %>%
      mutate(
        intensity_bg_sub = INTENSITY.ch1.mean - bg_threshold,
        log_intensity    = log1p(pmax(intensity_bg_sub, 0))  # floor at 0 before log
      )

    # Shared axis range from all recA cells
    recA_x <- range(df_recA$UMAP1)
    recA_y <- range(df_recA$UMAP2)

    # --- All recA cells together ---
    p_fluor_all <- ggplot(df_recA %>% arrange(log_intensity),
                          aes(x = UMAP1, y = UMAP2, colour = log_intensity)) +
      geom_point(size = point_size, alpha = point_alpha, shape = 16) +
      scale_colour_viridis_c(option = "D", name = "log(Intensity\n− background)") +
      labs(x = "UMAP 1", y = "UMAP 2",
           title = "recA Reporter — Background-Subtracted Intensity on UMAP",
           subtitle = sprintf("Background = 95th pctl NT intensity (%.1f)", bg_threshold))
    p_fluor_all <- apply_theme(p_fluor_all)
    ggsave(file.path(fig_dir, "sc_umap_fluorescence_recA_all.png"),
           plot = p_fluor_all, width = fig_width, height = fig_height, dpi = fig_dpi)

    # --- Faceted by knockdown/drug ---
    p_fluor_facet <- ggplot(df_recA %>% arrange(log_intensity),
                            aes(x = UMAP1, y = UMAP2, colour = log_intensity)) +
      geom_point(size = point_size * 0.7, alpha = point_alpha, shape = 16) +
      scale_colour_viridis_c(option = "D", name = "log(Intensity\n− background)") +
      facet_wrap(~ knockdown, ncol = 3) +
      coord_cartesian(xlim = recA_x, ylim = recA_y) +
      labs(x = "UMAP 1", y = "UMAP 2",
           title = "recA Reporter — Background-Subtracted Intensity per Condition",
           subtitle = sprintf("Background = 95th pctl NT intensity (%.1f)", bg_threshold))
    n_recA_panels <- n_distinct(df_recA$knockdown)
    recA_h <- 2 + ceiling(n_recA_panels / 3) * 3.5
    p_fluor_facet <- apply_theme(p_fluor_facet)
    ggsave(file.path(fig_dir, "sc_umap_fluorescence_recA_faceted.png"),
           plot = p_fluor_facet, width = fig_width, height = recA_h, dpi = fig_dpi)

    cat("  recA fluorescence plots saved.\n")
  } else {
    cat("  No recA cells found after filtering — skipping.\n")
  }
} else {
  cat("  INTENSITY.ch1.mean column not found — skipping fluorescence overlay.\n")
}

# ----------------------------- done -------------------------------------------

cat("\nDone. Figures saved to ", fig_dir, "/\n")
cat("  PCA plots:            sc_pca_density_{all,faceted,knockdowns,drugs}.png\n")
cat("  UMAP plots:           sc_umap_density_{all,faceted,knockdowns,drugs}.png\n")
cat("  Fluorescence:         sc_umap_fluorescence_recA_{all,faceted}.png\n")
