#!/usr/bin/env Rscript
# =============================================================================
# SupCon + OT Morphological Profiling — Publication Figures
# =============================================================================
# Generates publication-quality figures from the SupCon + OT pipeline outputs.
#
# Usage:  Rscript R/SupCon_OT_Figures.R
#
# Required packages:
#   tidyverse, dendextend, igraph, ggraph, tidygraph,
#   uwot, ggrepel, pheatmap, hdf5r, patchwork
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(dendextend)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggrepel)
  library(pheatmap)
  library(patchwork)
})

# All paths relative to project root (MorphologicalProfiling_Mtb/)
# Run from project root: Rscript R/SupCon_OT_Figures.R
if (file.exists("R/Theme.R")) source("R/Theme.R")
set.seed(42)

# =============================================================================
# CONFIGURATION
# =============================================================================

distance_file  <- "output/ot_supcon_distance_matrix.csv"
matches_file   <- "output/ot_supcon_ranked_matches.csv"
pathway_file   <- "gene_pathway_map.csv"
embedding_file <- "output/embeddings/cell_embeddings.h5"

drug_conditions <- c()  # TODO: list M. smegmatis drug conditions (e.g. "RIF_05", "INH_2")
fig_dir <- "figures/supcon_ot"
fig_dpi <- 300

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# =============================================================================
# COLOUR PALETTE
# =============================================================================

# NPG palette (ggsci) + 5 complementary muted tones for 15 pathways
npg_colors <- ggsci::pal_npg("nrc")(10)
extra_colors <- c("#6A6599", "#D6A461", "#79AF97", "#B24745", "#374E55")
all_colors <- c(npg_colors, extra_colors)

# Drug vs Gene type palette
type_palette <- c("Mutant" = "#386cb0", "Drug" = "#fdb462")

# =============================================================================
# DATA LOADING
# =============================================================================

cat("Loading data ...\n")

dist_df <- read_csv(distance_file, show_col_types = FALSE) %>%
  column_to_rownames(var = names(.)[1])
dist_mat <- as.matrix(dist_df)

pathway_map <- read_csv(pathway_file, show_col_types = FALSE)
gene_to_pw <- setNames(pathway_map$Pathway, pathway_map$Gene)

matches <- read_csv(matches_file, show_col_types = FALSE)

# Load gene name annotations for display labels
annotation_file <- "input_data/DetailedAll_UpdatedAnnotations.csv"
if (file.exists(annotation_file)) {
  ann <- readr::read_delim(annotation_file, delim = ";", show_col_types = FALSE,
                           col_select = c("Accession.no.", "geneName"))
  gene_name_map <- ann %>%
    filter(geneName != "-", !is.na(geneName)) %>%
    distinct(`Accession.no.`, .keep_all = TRUE) %>%
    {setNames(.$geneName, .$`Accession.no.`)}
} else {
  gene_name_map <- character(0)
}

# Helper: convert accession to gene name for display
display_label <- function(x) {
  ifelse(x %in% names(gene_name_map), gene_name_map[x], x)
}

# Apply display labels to distance matrix
rownames(dist_mat) <- display_label(rownames(dist_mat))
colnames(dist_mat) <- display_label(colnames(dist_mat))
rownames(dist_df) <- display_label(rownames(dist_df))
colnames(dist_df) <- display_label(colnames(dist_df))

# Gene-only distance matrix
gene_names_raw <- pathway_map$Gene[pathway_map$Gene %in% rownames(dist_df) |
                                   display_label(pathway_map$Gene) %in% rownames(dist_df)]
gene_names_display <- display_label(gene_names_raw)
gene_dist <- dist_mat[gene_names_display, gene_names_display]

# Pathway colour mapping (alphabetical)
pathways_sorted <- sort(unique(pathway_map$Pathway))
pathway_colors <- setNames(all_colors[seq_along(pathways_sorted)],
                           pathways_sorted)

# Condition type
cond_type <- ifelse(rownames(dist_mat) %in% drug_conditions, "Drug", "Mutant")
names(cond_type) <- rownames(dist_mat)

# =============================================================================
# FIGURE 1: Condition-level UMAP
# =============================================================================

cat("Figure 1: Condition-level UMAP ...\n")

# Use umap package with precomputed distances
umap_cfg <- umap::umap.defaults
umap_cfg$input <- "dist"
umap_cfg$n_neighbors <- 5
umap_cfg$min_dist <- 0.3
umap_result <- umap::umap(dist_mat, config = umap_cfg)
umap_coords <- umap_result$layout

umap_df <- tibble(
  condition = rownames(dist_mat),
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  type = cond_type[rownames(dist_mat)]
)

p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = type)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = condition), size = 2.5, max.overlaps = 30,
                  colour = "black", show.legend = FALSE) +
  scale_colour_manual(values = type_palette) +
  labs(title = "Condition-level UMAP (SupCon + Wasserstein)",
       x = "UMAP 1", y = "UMAP 2") +
  theme_Publication() +
  theme(legend.position = c(0.85, 0.15))

ggsave(file.path(fig_dir, "supcon_umap_conditions.png"), p1,
       width = 10, height = 8, dpi = fig_dpi)

# =============================================================================
# FIGURE 2: Condition UMAP with drug-gene edges
# =============================================================================

cat("Figure 2: Condition UMAP with edges ...\n")

top_k_edges <- 3
edge_data <- matches %>%
  filter(rank <= top_k_edges) %>%
  left_join(umap_df %>% select(condition, UMAP1, UMAP2),
            by = c("drug" = "condition")) %>%
  rename(x_drug = UMAP1, y_drug = UMAP2) %>%
  left_join(umap_df %>% select(condition, UMAP1, UMAP2),
            by = c("gene" = "condition")) %>%
  rename(x_gene = UMAP1, y_gene = UMAP2) %>%
  mutate(
    significant = replace_na(significant, FALSE),
    sig_label = ifelse(significant, "Significant (FDR < 0.05)", "Not significant"),
    d_scaled = (distance - min(distance)) / (max(distance) - min(distance)),
    edge_alpha = pmax(0.15, pmin(0.9, 1 - d_scaled)),
    edge_lw = ifelse(rank == 1, 1.2, 0.6)
  )

p2 <- ggplot() +
  # Edges
  geom_segment(data = edge_data,
               aes(x = x_drug, y = y_drug, xend = x_gene, yend = y_gene,
                   colour = sig_label, alpha = edge_alpha,
                   linewidth = edge_lw),
               show.legend = c(colour = TRUE, alpha = FALSE, linewidth = FALSE)) +
  scale_colour_manual(values = c("Significant (FDR < 0.05)" = "#7fc97f",
                                 "Not significant" = "#999999")) +
  scale_alpha_identity() +
  scale_linewidth_identity() +
  # Points on top
  geom_point(data = umap_df, aes(x = UMAP1, y = UMAP2, fill = type),
             shape = 21, size = 3, colour = "white", stroke = 0.5) +
  scale_fill_manual(values = type_palette) +
  geom_text_repel(data = umap_df, aes(x = UMAP1, y = UMAP2, label = condition),
                  size = 2.5, max.overlaps = 30, colour = "black") +
  labs(title = "Drug-gene matches (Top 3)",
       x = "UMAP 1", y = "UMAP 2") +
  theme_Publication() +
  theme(legend.position = "top",
        legend.direction = "vertical",
        legend.box = "horizontal",
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 4), order = 1),
         colour = guide_legend(order = 2))

ggsave(file.path(fig_dir, "supcon_umap_conditions_edges.png"), p2,
       width = 12, height = 8, dpi = fig_dpi)

# =============================================================================
# FIGURE 3: Distance heatmap with clustering
# =============================================================================

cat("Figure 3: Distance heatmap ...\n")

annotation_row <- data.frame(
  Type = factor(cond_type[rownames(dist_mat)], levels = c("Mutant", "Drug")),
  row.names = rownames(dist_mat)
)
ann_colors <- list(Type = type_palette)

pheatmap(dist_mat,
         clustering_method = "ward.D2",
         color = viridis::viridis(100, direction = -1),
         annotation_row = annotation_row,
         annotation_col = annotation_row,
         annotation_colors = ann_colors,
         fontsize = 8,
         fontsize_row = 7,
         fontsize_col = 7,
         main = "OT Distance Matrix (SupCon Embeddings)",
         filename = file.path(fig_dir, "supcon_distance_heatmap.png"),
         width = 12, height = 10)

# Capture heatmap as a grob for composite figure
p3_grob <- pheatmap(dist_mat,
                    clustering_method = "ward.D2",
                    color = viridis::viridis(100, direction = -1),
                    annotation_row = annotation_row,
                    annotation_col = annotation_row,
                    annotation_colors = ann_colors,
                    fontsize = 7,
                    fontsize_row = 6,
                    fontsize_col = 6,
                    main = "OT Distance Matrix (SupCon Embeddings)",
                    silent = TRUE)

# =============================================================================
# FIGURE 4: Top gene matches per drug (ranked bar chart)
# =============================================================================

cat("Figure 4: Ranked matches ...\n")

top_n <- 5
top_matches <- matches %>%
  filter(rank <= top_n) %>%
  mutate(
    significant = replace_na(significant, FALSE),
    bar_color = ifelse(significant, "Significant", "Not significant")
  )

# Re-order gene within each drug facet by descending distance (lowest at top)
top_matches <- top_matches %>%
  mutate(drug = factor(drug, levels = drug_conditions)) %>%
  arrange(drug, distance) %>%
  mutate(gene_drug = paste0(drug, "_", gene),
         gene_drug = fct_inorder(gene_drug))

p4 <- ggplot(top_matches, aes(x = distance, y = gene_drug, fill = bar_color)) +
  geom_col(width = 0.7) +
  geom_text(aes(x = 0, label = gene), hjust = 0, nudge_x = 0.002,
            size = 2.5, fontface = "bold") +
  facet_wrap(~ drug, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Significant" = "#7fc97f",
                               "Not significant" = "#CCCCCC")) +
  labs(title = "Top gene matches per drug (SupCon + OT)",
       x = "Wasserstein distance", y = NULL) +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave(file.path(fig_dir, "supcon_ranked_matches.png"), p4,
       width = 4 * length(drug_conditions), height = 5, dpi = fig_dpi)

# =============================================================================
# FIGURE 5: Gene dendrogram
# =============================================================================

cat("Figure 5: Gene dendrogram ...\n")

hc <- hclust(as.dist(gene_dist), method = "ward.D2")
dend <- as.dendrogram(hc)

# Colour leaves by pathway (must match dendrogram order)
leaf_labels <- labels(dend)
leaf_pw <- gene_to_pw[leaf_labels]
leaf_cols <- pathway_colors[leaf_pw]

labels_colors(dend) <- leaf_cols
labels_cex(dend) <- 0.85

png(file.path(fig_dir, "gene_dendrogram.png"),
    width = 12, height = 7, units = "in", res = fig_dpi)
par(mar = c(7, 4, 3, 14), xpd = TRUE)
plot(dend, main = "Gene Clustering by Morphological Similarity",
     ylab = "Ward Distance", xlab = "")
legend("topright", inset = c(-0.22, 0),
       legend = names(pathway_colors),
       fill = pathway_colors,
       cex = 0.55, ncol = 1, border = NA, bty = "n",
       title = "Pathway")
dev.off()

# =============================================================================
# FIGURE 6: Gene network graph
# =============================================================================

cat("Figure 6: Gene network ...\n")

# Edge threshold: 25th percentile
dist_vec <- gene_dist[lower.tri(gene_dist)]
threshold <- quantile(dist_vec, 0.25)
cat("  Network edge threshold (25th pctl):", round(threshold, 3), "\n")

# Build edge list
edge_list <- tibble()
for (i in seq_len(nrow(gene_dist))) {
  for (j in seq(i + 1, nrow(gene_dist))) {
    if (j > nrow(gene_dist)) break
    d <- gene_dist[i, j]
    if (d <= threshold) {
      edge_list <- bind_rows(edge_list, tibble(
        from = gene_names[i], to = gene_names[j],
        distance = d,
        similarity = 1 - (d / max(dist_vec))
      ))
    }
  }
}

# Node data
node_df <- tibble(
  name = gene_names,
  Pathway = gene_to_pw[gene_names]
)

# Build graph
g <- graph_from_data_frame(edge_list, directed = FALSE, vertices = node_df)
tg <- as_tbl_graph(g)

p6 <- ggraph(tg, layout = "fr") +
  geom_edge_link(aes(width = similarity, alpha = similarity),
                 colour = "grey50", show.legend = FALSE) +
  scale_edge_width(range = c(0.3, 2.5)) +
  scale_edge_alpha(range = c(0.2, 0.8)) +
  geom_node_point(aes(colour = Pathway), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_colour_manual(values = pathway_colors) +
  labs(title = "Gene Functional Network (Morphological Similarity)") +
  theme_Publication(base_size = 12) +
  theme(legend.position = "top",
        legend.text = element_text(size = 7),
        plot.title = element_text(hjust = 0.5),
        plot.title.position = "plot",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  guides(colour = guide_legend(nrow = 3, override.aes = list(size = 4)))

ggsave(file.path(fig_dir, "gene_network.png"), p6,
       width = 12, height = 8, dpi = fig_dpi)

# =============================================================================
# FIGURE 7: Cell-level UMAP (optional — requires hdf5r)
# =============================================================================

if (requireNamespace("hdf5r", quietly = TRUE) && file.exists(embedding_file)) {
  cat("Figure 7: Cell-level UMAP ...\n")
  library(hdf5r)

  h5 <- H5File$new(embedding_file, mode = "r")
  embeddings <- t(h5[["embeddings"]]$read())   # stored as 512 x N; transpose to N x 512
  labels_vec <- h5[["condition_labels"]]$read()
  is_control <- h5[["is_control"]]$read()
  h5$close_all()

  # Subsample ~5000 non-control cells
  non_ctrl <- which(!is_control)
  n_sample <- min(5000, length(non_ctrl))
  sample_idx <- sample(non_ctrl, n_sample)

  cell_umap <- uwot::umap(embeddings[sample_idx, ], n_neighbors = 15,
                           min_dist = 0.1, n_components = 2)

  cell_df <- tibble(
    UMAP1 = cell_umap[, 1],
    UMAP2 = cell_umap[, 2],
    condition = labels_vec[sample_idx]
  )

  # Use a large qualitative palette
  n_conds <- length(unique(cell_df$condition))
  cond_palette <- setNames(
    colorRampPalette(all_colors)(n_conds),
    sort(unique(cell_df$condition))
  )

  p7 <- ggplot(cell_df, aes(x = UMAP1, y = UMAP2, colour = condition)) +
    geom_point(size = 0.3, alpha = 0.4) +
    scale_colour_manual(values = cond_palette) +
    labs(title = "Cell-level Embeddings (SupCon)",
         x = "UMAP 1", y = "UMAP 2") +
    theme_Publication() +
    theme(legend.position = "top",
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.3, "cm")) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1),
                                 nrow = 2))

  ggsave(file.path(fig_dir, "supcon_umap_cells.png"), p7,
         width = 14, height = 10, dpi = fig_dpi)
} else {
  cat("Skipping cell-level UMAP (hdf5r not available or file missing)\n")
  p7 <- NULL
}

# =============================================================================
# COMPOSITE FIGURE: A4 single-page summary
# =============================================================================

cat("Composite figure: A4 summary ...\n")

# Panel A: Condition UMAP with edges (compact version)
pA <- p2 +
  labs(title = NULL, subtitle = NULL) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        plot.margin = margin(2, 2, 2, 2))

# Panel B: Heatmap (wrap pheatmap grob)
pB <- patchwork::wrap_elements(full = p3_grob$gtable)

# Panel C: Ranked matches (compact, 2 rows of facets)
p4_compact <- ggplot(top_matches, aes(x = distance, y = gene_drug, fill = bar_color)) +
  geom_col(width = 0.7) +
  geom_text(aes(x = 0, label = gene), hjust = 0, nudge_x = 0.002,
            size = 2, fontface = "bold") +
  facet_wrap(~ drug, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("Significant" = "#7fc97f",
                               "Not significant" = "#CCCCCC")) +
  labs(x = "Wasserstein distance", y = NULL) +
  theme_Publication(base_size = 9) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(2, 2, 2, 2))

# Panel D: Gene network (compact)
pD <- p6 +
  labs(title = NULL) +
  theme(legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 7),
        plot.margin = margin(2, 2, 2, 2))

# Assemble: 2x2 layout
#   A (UMAP+edges)  |  B (heatmap)
#   C (ranked bars)  |  D (network)
composite <- (pA | pB) / (p4_compact | pD) +
  plot_annotation(
    title = "SupCon + Optimal Transport Morphological Profiling",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

# A4 dimensions: 210 x 297 mm = 8.27 x 11.69 inches
ggsave(file.path(fig_dir, "composite_figure.png"), composite,
       width = 8.27, height = 11.69, dpi = fig_dpi)
ggsave(file.path(fig_dir, "composite_figure.pdf"), composite,
       width = 8.27, height = 11.69)

# =============================================================================
cat("\nAll figures saved to", fig_dir, "\n")
