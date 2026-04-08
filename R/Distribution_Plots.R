# ==============================================================================
# Single-cell ridge plot distributions of core morphological features
# - One plot per reporter, faceted by feature
# - Separate plots for ATC knockdowns vs drug treatments
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)
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

name_corrections <- list()
replicate_groups <- list()
merge_replicates <- FALSE

drug_experiment_pattern <- "drug"

# Core features to visualise
features_to_plot <- c(
  "SHAPE.length", "SHAPE.width", "SHAPE.area", "SHAPE.aspectRatio",
  "SHAPE.circularity", "SHAPE.roundness", "SHAPE.solidity",
  "INTENSITY.ch1.mean"
)

# Pretty labels for facets
feature_labels <- c(
  "SHAPE.length"      = "Length",
  "SHAPE.width"       = "Width",
  "SHAPE.area"        = "Area",
  "SHAPE.aspectRatio" = "Aspect Ratio",
  "SHAPE.circularity" = "Circularity",
  "SHAPE.roundness"   = "Roundness",
  "SHAPE.solidity"    = "Solidity",
  "INTENSITY.ch1.mean" = "Reporter Intensity"
)

# Figure output
fig_dir    <- "figures"
fig_dpi    <- 300
fig_width  <- 12
fig_height <- 10

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

# ----------------------------- load & parse -----------------------------------

cat("Loading data...\n")
raw <- bind_rows(lapply(input_files, read_csv, show_col_types = FALSE))

# Check required columns exist
missing <- setdiff(features_to_plot, names(raw))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

# Parse sample metadata
meta <- map_dfr(unique(raw$EXPERIMENT), parse_sample_name)

# Apply name corrections
for (fix in name_corrections) {
  idx <- meta$EXPERIMENT == fix$experiment
  if (any(idx)) {
    for (field in setdiff(names(fix), "experiment")) {
      meta[[field]][idx] <- fix[[field]]
    }
  }
}

# Label replicates
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

# Join metadata to raw data
df <- raw %>%
  select(EXPERIMENT, all_of(features_to_plot)) %>%
  left_join(meta, by = "EXPERIMENT")

# Classify controls and drug samples
df <- df %>%
  mutate(
    is_control = is_control_label(knockdown),
    is_drug    = str_detect(experiment_type, regex(drug_experiment_pattern, ignore_case = TRUE))
  )

cat(sprintf("Loaded %d cells across %d conditions\n", nrow(df), n_distinct(df$EXPERIMENT)))

# ----------------------------- pivot to long ----------------------------------

df_long <- df %>%
  pivot_longer(
    cols      = all_of(features_to_plot),
    names_to  = "feature",
    values_to = "value"
  ) %>%
  mutate(feature_label = feature_labels[feature]) %>%
  filter(is.finite(value)) %>%
  # Trim outliers beyond the 98th percentile per feature
  group_by(feature) %>%
  filter(value <= quantile(value, 0.98)) %>%
  ungroup()

# ----------------------------- plotting function ------------------------------

make_ridge_plot <- function(data, title_suffix) {
  # Order: controls first, then alphabetical
  ctrl   <- sort(unique(data$knockdown[data$is_control]))
  others <- sort(setdiff(unique(data$knockdown), ctrl))
  level_order <- c(others, ctrl)  # ridges plot bottom-to-top, so controls at bottom

  data <- data %>%
    mutate(knockdown = factor(knockdown, levels = level_order))

  # Subsample for jitter layer (full dataset too dense for individual points)
  max_pts_per_group <- 1000
  set.seed(42)
  data_sub <- data %>%
    group_by(knockdown, feature_label) %>%
    slice_sample(n = max_pts_per_group) %>%
    ungroup()

  p <- ggplot(data, aes(x = value, y = knockdown, fill = knockdown)) +
    geom_density_ridges(
      alpha          = 0.7,
      scale          = 0.5,
      rel_min_height = 0.01
    ) +
    geom_point(
      data     = data_sub,
      aes(x = value, y = knockdown, colour = knockdown),
      position = position_jitter(width = 0, height = 0.2),
      size     = 0.3,
      alpha    = 0.25,
      shape    = 16
    ) +
    facet_wrap(~ feature_label, scales = "free_x", ncol = 2) +
    labs(x = NULL, y = NULL, title = title_suffix) +
    theme_minimal(base_size = 13)

  if (exists("theme_Publication")) {
    p <- p + theme_Publication() +
      theme(
        legend.position  = "none",
        strip.text       = element_text(face = "bold", size = 11),
        axis.text.y      = element_text(size = 10),
        plot.title        = element_text(hjust = 0.5, size = 14, face = "bold")
      )
  } else {
    p <- p + theme(legend.position = "none")
  }

  if (exists("scale_fill_Publication")) {
    p <- p + scale_fill_Publication()
  }
  if (exists("scale_colour_Publication")) {
    p <- p + scale_colour_Publication()
  }

  p
}

# ----------------------------- generate plots ---------------------------------

reporters <- sort(unique(df_long$reporter))

for (rep in reporters) {
  # --- ATC knockdown samples ---
  atc_data <- df_long %>%
    filter(reporter == rep, !is_drug)

  if (nrow(atc_data) > 0) {
    cat(sprintf("Plotting ATC distributions for %s...\n", rep))
    p <- make_ridge_plot(atc_data, paste0("ATC Knockdowns - ", rep))
    ggsave(
      file.path(fig_dir, paste0("distributions_ATC_", rep, ".png")),
      plot = p, width = fig_width, height = fig_height, dpi = fig_dpi
    )
  }

  # --- Drug treatment samples ---
  drug_data <- df_long %>%
    filter(reporter == rep, is_drug)

  if (nrow(drug_data) > 0) {
    cat(sprintf("Plotting drug distributions for %s...\n", rep))
    p <- make_ridge_plot(drug_data, paste0("Drug Treatments - ", rep))
    ggsave(
      file.path(fig_dir, paste0("distributions_drug_", rep, ".png")),
      plot = p, width = fig_width, height = fig_height, dpi = fig_dpi
    )
  }
}

# TODO: Add M. smegmatis-specific combined figure here when conditions are known.
# The Mtb branch had a combined iniB + recA reporter figure; replace with
# an equivalent comparison relevant to M. smegmatis conditions.

cat("Done. Figures saved to ", fig_dir, "/\n")
