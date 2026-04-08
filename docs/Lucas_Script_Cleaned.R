# ==============================================================================
# UMAP + consensus clustering pipeline (cleaned + rationalised)
# - Loads/standardises datasets
# - Computes S-scores (mean + CV) vs plasmid controls for selected features
# - Runs multiple UMAP+HDBSCAN fits, builds a co-association matrix
# - Hierarchical clustering on consensus distances
# - Produces final UMAP plot + optional annotations
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, readr, stringr, ggplot2, tibble, tidyr
  library(dbscan)      # hdbscan
  library(uwot)        # umap
  library(ggrepel)     # geom_text_repel
  library(plotly)      # ggplotly (optional)
})

# Optional: custom theme (only if present)
if (file.exists("Scripts/Theme.R")) source("Scripts/Theme.R")

set.seed(10)

# ----------------------------- configuration ----------------------------------

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
  "SHAPE.width.variation",
  "MAXIMA"
)

umap_params <- list(n_neighbors = 12, min_dist = 0)
hdbscan_params <- list(minPts = 12)
n_umap_runs <- 50
consensus_k <- 3   # cutree() k

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

# Compute S-scores for mean and CV vs plasmid controls:
# S = (value - mean(plasmid_values)) / sd(plasmid_values)
compute_s_scores <- function(df, variables) {
  require_cols(df, c("Name", "EXPERIMENT"), "df")
  require_cols(df, variables, "df")
  
  gene_map <- df %>%
    distinct(Name, EXPERIMENT)
  
  gene_means <- df %>%
    group_by(Name) %>%
    summarise(across(all_of(variables), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  gene_cvs <- df %>%
    group_by(Name) %>%
    summarise(across(all_of(variables), ~ cv(.x)), .groups = "drop")
  
  long_means <- gene_means %>%
    pivot_longer(-Name, names_to = "feature", values_to = "value") %>%
    mutate(stat = "mean")
  
  long_cvs <- gene_cvs %>%
    pivot_longer(-Name, names_to = "feature", values_to = "value") %>%
    mutate(stat = "cv")
  
  long_all <- bind_rows(long_means, long_cvs)
  
  baseline <- long_all %>%
    filter(str_detect(Name, "Plasmid")) %>%
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
    select(Name, out_feature, s) %>%
    pivot_wider(names_from = out_feature, values_from = s)
  
  gene_map %>%
    left_join(scored, by = "Name")
}

make_umap_hdbscan <- function(mat, seed, n_neighbors, min_dist, minPts) {
  set.seed(seed)
  emb <- uwot::umap(mat, n_neighbors = n_neighbors, min_dist = min_dist)
  emb <- as.data.frame(emb)
  fit <- dbscan::hdbscan(emb[, 1:2, drop = FALSE], minPts = minPts)
  list(embedding = emb, cluster = as.character(fit$cluster))
}

consensus_from_clusters <- function(cluster_matrix) {
  # cluster_matrix: n_genes x n_runs (character/factor)
  n_runs <- ncol(cluster_matrix)
  n <- nrow(cluster_matrix)
  
  assoc_counts <- matrix(0, n, n)
  for (r in seq_len(n_runs)) {
    assoc_counts <- assoc_counts + outer(cluster_matrix[, r], cluster_matrix[, r], FUN = "==")
  }
  assoc_prop <- assoc_counts / n_runs
  as.dist(1 - assoc_prop)
}

# ----------------------------- load + normalise data ---------------------------

cell_data <- readr::read_csv("InputData/allClassifiedCells_final.csv", show_col_types = FALSE) %>%
  mutate(
    Name = as.character(Name),
    Name = if_else(
      Name == "Plasmid",
      paste(Name, Replica, WELL_ID),
      Name
    ),
    Name = as.factor(Name)
  ) %>%
  select(-WELL_ID) %>%
  filter(!Name %in% c("rplX", "hisI"))

drugs <- readr::read_csv("InputData/drugData_processed", show_col_types = FALSE) %>%
  mutate(
    Replica = 1L,
    EXPERIMENT = paste(Drug, Concentration, sep = "_"),
    EXPERIMENT = as.factor(EXPERIMENT),
    Name = EXPERIMENT
  ) %>%
  filter(Concentration %in% c("1X", "2X", "4X")) %>%
  select(-Drug, -Concentration, -ImageReplica)

lucas_2024_rep1 <- readr::read_csv("data_output_lucas.csv", show_col_types = FALSE) %>%
  mutate(
    Name = paste0("L_", EXPERIMENT),
    Name = if_else(Name == "L_D4", "Plasmid L_D4", Name),
    Name = as.factor(Name),
    Replica = 1L
  )

lucas_2024_rep2 <- bind_rows(
  readr::read_csv("Lucas/data_output_rep2_AB.csv", show_col_types = FALSE),
  readr::read_csv("Lucas/data_output_rep2_CD.csv", show_col_types = FALSE)
) %>%
  mutate(
    Name = paste0("LRep_", EXPERIMENT),
    Name = if_else(Name == "LRep_D4", "Plasmid LRep_D4", Name),
    Name = as.factor(Name),
    Replica = 2L
  )

lucas_2024_rep3 <- readr::read_csv("Lucas/data_output_rep3.csv", show_col_types = FALSE) %>%
  mutate(
    Name = paste0("LRep2_", EXPERIMENT),
    Name = as.factor(Name),
    Replica = 3L
  )

lucas_2023 <- readr::read_csv("InputData/data_output_12_04_23.csv", show_col_types = FALSE) %>%
  mutate(
    Name = as.factor(EXPERIMENT),
    Replica = 1L
  )
# NOTE: the original code forces factor levels to a fixed order. That is brittle.
# If you truly need this mapping, use an explicit recode() keyed by EXPERIMENT.

lucas_priorities <- readr::read_csv("Lucas/PriorityTargets/data_output_310522.csv", show_col_types = FALSE) %>%
  mutate(
    Name = as.factor(EXPERIMENT),
    Replica = 1L
  )

well_to_gene <- c(
  "A1" = "fadD32",
  "A2" = "kasA",
  "A3" = "clpP2",
  "A4" = "clpP1",
  "A5" = "nadD",
  "A6" = "secA1",
  "A7" = "topA",
  "A8" = "acpP",
  "B1" = "murA",
  "ParB" = "Plasmid_Priority"
)

lucas_priorities <- lucas_priorities %>%
  mutate(
    Name = recode(as.character(EXPERIMENT), !!!well_to_gene, .default = as.character(EXPERIMENT)),
    Name = factor(Name)
  )

# NOTE: original rbind() read the same file twice; removed.

lucas_reporters <- bind_rows(
  readr::read_csv("Lucas/ReporterStrains/data_output_cydA_atpE.csv", show_col_types = FALSE),
  readr::read_csv("Lucas/ReporterStrains/data_output_iniB_embA.csv", show_col_types = FALSE),
  readr::read_csv("Lucas/ReporterStrains/data_output_iniB_inhA.csv", show_col_types = FALSE),
  readr::read_csv("Lucas/ReporterStrains/data_output_recA_dnaE1.csv", show_col_types = FALSE)
) %>%
  mutate(
    Name = as.factor(EXPERIMENT),
    Replica = 1L
  )

both_reps <- bind_rows(
  cell_data,
  drugs,
  lucas_2024_rep1,
  lucas_2024_rep2,
  lucas_2024_rep3,
  lucas_priorities,
  lucas_reporters,
  lucas_2023
)

# Manual experiment label fixes (keep these in one place)
both_reps <- both_reps %>%
  mutate(
    EXPERIMENT = case_when(
      Name == "ftsK" ~ "02 B2",
      Name == "rpoB" ~ "02 H1",
      Name == "rpoC" ~ "01 A5",
      TRUE ~ as.character(EXPERIMENT)
    )
  )

# ----------------------------- compute S-scores --------------------------------

s_values <- compute_s_scores(both_reps, variables_of_interest)

# Build numeric matrix for UMAP
feature_mat <- s_values %>%
  select(-Name, -EXPERIMENT) %>%
  as.matrix()
rownames(feature_mat) <- as.character(s_values$Name)

# Drop columns that are all NA or contain any NA (mirrors original intent)
feature_mat <- feature_mat[, colSums(is.na(feature_mat)) != nrow(feature_mat), drop = FALSE]
feature_mat <- feature_mat[, colSums(is.na(feature_mat)) == 0, drop = FALSE]

genes <- rownames(feature_mat)

# Optional: separate drug rows (if you don’t want them included in clustering)
drug_names <- unique(as.character(drugs$Name))
is_drug <- genes %in% drug_names

feature_mat_main <- feature_mat[!is_drug, , drop = FALSE]
genes_main <- rownames(feature_mat_main)

# ----------------------- repeated UMAP + HDBSCAN runs --------------------------

runs <- map(
  1:n_umap_runs,
  ~ make_umap_hdbscan(
    mat = feature_mat_main,
    seed = .x,
    n_neighbors = umap_params$n_neighbors,
    min_dist = umap_params$min_dist,
    minPts = hdbscan_params$minPts
  )
)

cluster_mat <- do.call(cbind, lapply(runs, `[[`, "cluster"))
rownames(cluster_mat) <- genes_main
colnames(cluster_mat) <- paste0("run_", seq_len(n_umap_runs))

# Consensus distance from co-association
dist_consensus <- consensus_from_clusters(cluster_mat)
hc <- hclust(dist_consensus)
clusters_consensus <- cutree(hc, k = consensus_k)

# ----------------------------- final embedding ---------------------------------

final <- make_umap_hdbscan(
  mat = feature_mat_main,
  seed = 10,
  n_neighbors = umap_params$n_neighbors,
  min_dist = umap_params$min_dist,
  minPts = hdbscan_params$minPts
)

d_umap <- final$embedding %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(
    gene = genes_main,
    cluster = as.character(clusters_consensus),
    cluster = if_else(str_detect(gene, "Plasmid"), "Control", cluster),
    cluster = factor(cluster)
  )

p_clusters <- ggplot(d_umap, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point() +
  labs(x = "UMAP1", y = "UMAP2", color = "Cluster") +
  theme_minimal()

# If you have a custom theme:
if (exists("theme_Publication")) {
  p_clusters <- p_clusters +
    theme_Publication() +
    theme(legend.position = "right")
}

p_clusters
# plotly::ggplotly(p_clusters)

# -------------------------- Lucas annotations (vectorised) ---------------------

lucas_key <- readr::read_csv("lucasKnockdowns.csv", show_col_types = FALSE)
require_cols(lucas_key, c("Imaging.Sequence.No.", "Name_Msm"), "lucas_key")

d_umap <- d_umap %>%
  mutate(
    is_lucas = str_detect(gene, "^L(_|Rep_|Rep2_)"),
    lucas_well_raw = if_else(is_lucas, str_sub(gene, -2L), NA_character_)
  ) %>%
  # Map wells to gene names
  left_join(
    lucas_key %>%
      transmute(lucas_well_raw = as.character(`Imaging.Sequence.No.`),
                lucas_target = as.character(Name_Msm)),
    by = "lucas_well_raw"
  ) %>%
  mutate(
    lucas_well = if_else(is_lucas & !is.na(lucas_target), lucas_target, lucas_well_raw),
    cluster = if_else(is_lucas, "New Additions", as.character(cluster)),
    cluster = factor(cluster)
  ) %>%
  select(-lucas_target)

# Manual fixes kept as explicit recodes (easy to audit)
d_umap <- d_umap %>%
  mutate(
    lucas_well = recode(
      lucas_well,
      "D5" = "secA", "D6" = "secA",
      "D7" = "purB", "D8" = "purB",
      "A1" = "secA", "A2" = "rpsB", "A3" = "cdsA", "A4" = "ribD",
      "A5" = "frr",  "A6" = "MSMEG_2936",
      "B1" = "pgk",  "B2" = "coaA", "B3" = "tatA"
    )
  )

p_clusters2 <- ggplot(d_umap, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point() +
  geom_point(data = filter(d_umap, is_lucas), color = "black") +
  labs(x = "UMAP1", y = "UMAP2", color = "Cluster") +
  theme_Publication()

p_clusters2

# ----------------------- Priority targets (DRY version) ------------------------

priority_sets <- list(
  fadD32 = c("fadD32", "fabD", "inhA", "mmpL3"),
  kasA   = c("kasA",   "fabD", "inhA", "mmpL3"),
  clpP2  = c("clpP2",  "clpC1", "clpP1", "clpX"),
  clpP1  = c("clpP1",  "clpC1", "clpP2", "clpX"),
  nadD   = c("nadD",   "nadE"),
  secA1  = c("secA1",  "secD", "secF"),
  topA   = c("topA",   "gyrA", "gyrB"),
  acpP   = c("acpP"),
  murA   = c("murA", "murB", "murC", "murE")
)

priority_compare_all <- imap_dfr(priority_sets, ~ {
  tibble(priority_gene = .y, gene = .x) %>%
    left_join(d_umap, by = "gene") %>%
    mutate(priority = if_else(gene == priority_gene, "target", "similar"))
})

ggplot(d_umap, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(alpha = 0.4) +
  geom_point(data = priority_compare_all, color = "black") +
  geom_text_repel(data = priority_compare_all, aes(label = gene), color = "black") +
  facet_wrap(~ priority_gene) +
  theme_Publication()
