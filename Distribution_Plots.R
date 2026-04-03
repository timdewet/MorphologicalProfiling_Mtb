# ==============================================================================
# Single-cell ridge plot distributions of core morphological features
# - One plot per reporter, faceted by feature
# - Separate plots for ATC knockdowns vs drug treatments
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggridges)
})

if (file.exists("Theme.R")) source("Theme.R")

set.seed(10)

# ----------------------------- configuration ----------------------------------

input_files <- c(
  "input_data/data_extraction_18_03_25.csv",
  "input_data/all_morphology_combined.csv"
)

name_separator <- "__"

control_labels <- c("NT", "No_drug")

# Helper: match control labels including numbered replicates (NT_1, NT_2, ...)
is_control_label <- function(x) {
  patt <- paste0("^(", paste(control_labels, collapse = "|"), ")(_.+)?$")
  str_detect(x, patt)
}

# Name corrections (same as main pipeline)
name_corrections <- list(
  list(experiment = "WT_Reporters_+_drug__imiB__Inn",     reporter = "iniB", knockdown = "INH"),
  list(experiment = "WT_Reporters_+_drug__imiB__EMB",     reporter = "iniB"),
  list(experiment = "WT_Reporters_+_drug__imiB__No_drug", reporter = "iniB"),
  list(experiment = "WT_Reporters_+_drug__imiB__RIF",     reporter = "iniB"),
  list(experiment = "ATC_Strains__recA__dnaW2",           knockdown = "dnaN1_rep2")
)

# Replicate groups
replicate_groups <- list(
  dnaN = c("ATC_Strains__recA__dnaN1", "ATC_Strains__recA__dnaW2")
)
merge_replicates <- FALSE

# Which experiment_type values identify drug-treated samples?
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

parse_sample_name <- function(name, sep = "__") {
  parts <- str_split_fixed(name, fixed(sep), n = 3)
  tibble(
    EXPERIMENT      = name,
    experiment_type = parts[, 1],
    reporter        = parts[, 2],
    knockdown       = parts[, 3]
  )
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
meta <- map_dfr(unique(raw$EXPERIMENT), parse_sample_name, sep = name_separator)

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

# -------------------- combined iniB + recA figure -----------------------------

cat("Plotting combined iniB + recA figure...\n")

# Colour palette: reporter x experiment type
#   iniB ATC  = light green,  iniB drug = dark green
#   recA ATC  = light maroon, recA drug = dark maroon
fill_colours <- c(
  "iniB - Knockdown"  = "#5ca05c",
  "iniB - Drug"       = "#1a5e1a",
  "recA - Knockdown"  = "#c46e6e",
  "recA - Drug"       = "#7b1a1a"
)

combined <- df_long %>%
  filter(reporter %in% c("iniB", "recA")) %>%
  mutate(
    type_label  = if_else(is_drug, "Drug", "Knockdown"),
    colour_group = paste0(reporter, " - ", type_label)
  )

# Create a y-axis label combining reporter and knockdown
combined <- combined %>%
  mutate(y_label = paste0(reporter, " - ", knockdown))

# Order y-axis: within each type, controls at top, then alphabetical
ctrl_labels <- combined %>%
  filter(is_control) %>%
  distinct(y_label) %>%
  pull(y_label) %>%
  sort()
other_labels <- combined %>%
  filter(!is_control) %>%
  distinct(y_label) %>%
  pull(y_label) %>%
  sort()
level_order <- c(other_labels, ctrl_labels)
combined <- combined %>%
  mutate(y_label = factor(y_label, levels = level_order))

# Subsample for jitter
max_pts <- 1000
set.seed(42)
combined_sub <- combined %>%
  group_by(y_label, feature_label, colour_group) %>%
  slice_sample(n = max_pts) %>%
  ungroup()

p_combined <- ggplot(combined,
                     aes(x = value, y = y_label,
                         fill = colour_group)) +
  geom_density_ridges(
    alpha          = 0.7,
    scale          = 0.5,
    rel_min_height = 0.01
  ) +
  geom_point(
    data     = combined_sub,
    aes(x = value, y = y_label, colour = colour_group),
    position = position_jitter(width = 0, height = 0.2),
    size     = 0.3,
    alpha    = 0.25,
    shape    = 16
  ) +
  facet_grid(type_label ~ feature_label, scales = "free",
             space = "free_y") +
  scale_fill_manual(values = fill_colours, name = NULL) +
  scale_colour_manual(values = fill_colours, name = NULL) +
  guides(
    fill   = guide_legend(override.aes = list(alpha = 0.9)),
    colour = "none"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13)

if (exists("theme_Publication")) {
  p_combined <- p_combined + theme_Publication() +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.text      = element_text(size = 10),
      strip.text       = element_text(face = "bold", size = 11),
      axis.text.y      = element_text(size = 9),
      plot.margin      = unit(c(5, 5, 5, 5), "mm")
    )
} else {
  p_combined <- p_combined +
    theme(legend.position = "bottom", legend.direction = "horizontal")
}

ggsave(
  file.path(fig_dir, "distributions_combined_iniB_recA.png"),
  plot = p_combined, width = 16, height = 8, dpi = fig_dpi
)

cat("Done. Figures saved to ", fig_dir, "/\n")
