# =============================================================================
# Drug-treated iniB-mScarlet cells vs the iniB CRISPRi library:
# nearest genetic (knockdown) neighbours for each drug x dose
#
# Drug data : Plate60 day 4 (EMB, INH, MOX, RIF x 0.25/0.5/1x MIC), 1:50 dilution
#             only (the recommended working dilution).
# Library   : 67-strain iniB-mScarlet CRISPRi panel (63 KDs + IFT, SSB, NT1-4),
#             canonical per-batch NT-corrected profiles (output_inib/
#             strain_profiles_corrected.csv) + raw NT cells.
#
# Method mirrors the library pipeline (inib_morphology_analysis.py /
# ribf_vs_inib.py): 36 shape features -> per-cell robust z vs a reference ->
# per-condition median profile -> standardise features across all profiles ->
# PCA (12 PCs) -> Euclidean distance.  Shape only: reporter intensity is not
# comparable across imaging runs.
#
# The drug plate is a SEPARATE imaging run with NO untreated / NT control, so
# where "untreated" sits in this run is unknown.  Two anchors bracket that:
#   * libNT : drug cells robust-z vs the library's pooled NT1-4 cells (assumes no
#             run-to-run shift); library = canonical per-batch corrected profiles.
#   * plate : drug cells robust-z vs all 1:50 drug-plate cells; library profiles
#             re-centred on the library median (assumes the average drug
#             condition ~ the average knockdown).  Removes any run-level shift.
# Neighbours that hold under BOTH anchors are the defensible ones.  Within each
# anchor, FOV bootstrap gives the technical stability of each neighbour.  One
# well per condition -> descriptive, no cell-count p-values.
#
# Sanity check built in: the library carries the direct targets of three drugs
# (EMB -> embA, INH -> inhA, MOX -> gyrA/gyrB; no rpoB for RIF), so each drug's
# expected pathway is scored by AUC against a label-permutation null.
#
# Run from the project root:
#   Rscript optimisations/dilution_factors/drug_genetic_neighbours.R \
#     [drug_features.csv] [library_dir] [library_profiles.csv]
# Outputs -> <drug_features dir>/drug_neighbours/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

if (file.exists("Theme.R")) source("Theme.R") else
  stop("Theme.R not found - run from project root")

args <- commandArgs(trailingOnly = TRUE)
drug_file <- if (length(args) >= 1) args[1] else "optimisations/dilution_factors/all_features_day2.csv"
lib_dir   <- if (length(args) >= 2) args[2] else "input_data/iniB_allstrains"
lib_prof_file <- if (length(args) >= 3) args[3] else "output_inib/strain_profiles_corrected.csv"

out_dir <- file.path(dirname(normalizePath(drug_file, mustWork = TRUE)), "drug_neighbours")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(p, name, w, h, dpi = 300) {
  # cairo is unavailable without XQuartz; quartz (pdf) + ragg (png) instead
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white",
         device = function(filename, ...) grDevices::quartz(type = "pdf", file = filename, ...))
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = dpi, bg = "white", device = ragg::agg_png)
}

# ---- parameters -------------------------------------------------------------
DILUTION   <- 50L
TOPK       <- 8L      # neighbours reported per drug x dose
TOP_STABLE <- 5L      # bootstrap stability = how often a gene is in the top 5
N_BOOT     <- 200L
N_PERM     <- 10000L
PROFILE_PCA <- 12L    # as CONFIG["profile_pca"] in the library pipeline
set.seed(42)

FEATURES <- c(
  "area_um2", "perimeter_um", "eccentricity", "solidity",
  "equivalent_diameter_um", "major_axis_length_um", "minor_axis_length_um",
  "feret_diameter_max_um", "length_um", "width_median_um", "width_mean_um",
  "width_max_um", "width_min_um", "width_std_um", "max_width_position_frac",
  "min_width_position_frac", "sinuosity", "branch_count",
  "area_um2_subpixel", "perimeter_um_subpixel", "width_amplitude_um",
  "width_variation", "aspect_ratio", "circularity_subpixel",
  "pole_width_um", "pole_taper", "angularity_mean", "angularity_max",
  "angularity_median", "angularity_std", "angularity_amplitude",
  "angularity_variation", "curvature_mean", "curvature_std",
  "curvature_max", "roundness")
NT_STRAINS <- c("NT1", "NT2", "NT3", "NT4")

# Functional groups as in inib_morphology_analysis.py (annotation only)
FUNCTIONAL_GROUPS <- list(
  "Mycolic acid / FAS-II"     = c("inhA", "kasA", "hadA", "hadB", "hadC", "fabD",
                                  "acpM", "fad32", "accA3", "pksB", "mmpL3"),
  "Arabinogalactan"           = c("embA", "dprE1", "glfT2"),
  "Peptidoglycan / D-Ala"     = c("murA", "ddl", "alr"),
  "Division / shape"          = c("ftsZ", "wag31", "pknB"),
  "DNA replication/repair"    = c("gyrA", "gyrB", "dnaN", "dnaE1", "dnaZX", "polA",
                                  "topA", "SSB", "hupB"),
  "Transcription/translation" = c("sigA", "rho", "rpsM", "fusA1", "leuS", "lysS"),
  "Respiration / ATP"         = c("atpE", "qcrB", "ndh", "cya"),
  "Central metabolism"        = c("fum", "glcB"),
  "Amino-acid biosynth"       = c("trpA", "trpB", "trpC", "aroG", "procA", "metF"),
  "Cofactor / vitamin"        = c("birA", "nadD", "ribA2", "cobQ2", "menH", "dxs1", "ppt"),
  "Proteostasis/secretion"    = c("clpP1", "clpP2", "htrA", "prcB", "secA1"),
  "Transport / efflux"        = c("efpA", "moxR1"),
  "Other / unknown"           = c("ipdC", "IFT"),
  "Control (NT)"              = NT_STRAINS)
gene_group <- setNames(rep(names(FUNCTIONAL_GROUPS), lengths(FUNCTIONAL_GROUPS)),
                       unlist(FUNCTIONAL_GROUPS))

# 8 Okabe-Ito super-groups used in the iniB A4 figure (colourblind-safe)
super_map <- c(
  "Mycolic acid / FAS-II" = "Cell wall / lipid", "Arabinogalactan" = "Cell wall / lipid",
  "Peptidoglycan / D-Ala" = "Cell wall / lipid", "Division / shape" = "Division / shape",
  "DNA replication/repair" = "DNA replication / repair",
  "Transcription/translation" = "Transcription / translation",
  "Respiration / ATP" = "Respiration / metabolism", "Central metabolism" = "Respiration / metabolism",
  "Amino-acid biosynth" = "Respiration / metabolism", "Cofactor / vitamin" = "Cofactor / vitamin",
  "Proteostasis/secretion" = "Proteostasis / transport",
  "Transport / efflux" = "Proteostasis / transport", "Other / unknown" = "Proteostasis / transport",
  "Control (NT)" = "Control (NT)")
super_levels <- c("Cell wall / lipid", "Division / shape", "DNA replication / repair",
                  "Transcription / translation", "Respiration / metabolism",
                  "Cofactor / vitamin", "Proteostasis / transport", "Control (NT)")
super_cols <- setNames(c("#D55E00", "#E69F00", "#0072B2", "#56B4E9",
                         "#009E73", "#F0E442", "#CC79A7", "#000000"), super_levels)

# Expected pathway (and direct target, where the library has it) per drug
EXPECTED <- data.table(
  drug    = c("EMB", "INH", "MOX", "RIF"),
  pathway = c("Arabinogalactan", "Mycolic acid / FAS-II",
              "DNA replication/repair", "Transcription/translation"),
  target  = c("embA", "inhA", "gyrA,gyrB", NA))

conc_levels <- c("0.25xMIC", "0.5xMIC", "1xMIC")

# ---- helpers -----------------------------------------------------------------
# Same as robust_z in ribf_vs_inib.py: centre on reference median, scale by
# MAD * 1.4826, falling back to the reference SD and then 1.
robust_z <- function(X, ref) {
  med <- apply(ref, 2, median, na.rm = TRUE)
  mad <- apply(sweep(ref, 2, med), 2, function(v) median(abs(v), na.rm = TRUE)) * 1.4826
  sdv <- apply(ref, 2, sd, na.rm = TRUE)
  sc  <- ifelse(mad > 1e-9, mad, ifelse(sdv > 1e-9, sdv, 1))
  Z <- sweep(sweep(X, 2, med), 2, sc, "/")
  Z[!is.finite(Z)] <- 0
  Z
}

col_medians <- function(Z) apply(Z, 2, median)

# Library build_network distance: standardise each feature across profiles
# (constant features -> scale 1), PCA, Euclidean on the first PROFILE_PCA PCs.
profile_dist <- function(prof) {
  sdv <- apply(prof, 2, sd)
  P <- scale(prof, center = TRUE, scale = ifelse(sdv > 1e-12, sdv, 1))
  k <- max(2L, min(PROFILE_PCA, nrow(P) - 1L, ncol(P)))
  pcs <- prcomp(P, center = FALSE)$x[, seq_len(k), drop = FALSE]
  as.matrix(dist(pcs))
}

# ---- load drug cells (1:50 only) --------------------------------------------
drug <- fread(drug_file, select = c("well", "fov", FEATURES))
drug[, c("dc", "reporter", "dil") := tstrsplit(well, "__")]
drug[, c("drug", "conc") := tstrsplit(dc, "_")]
drug[, dil := as.integer(sub("1in", "", sub("_focused", "", dil)))]
drug <- drug[dil == DILUTION]
drug[, cond := paste(drug, conc, sep = "_")]
message(sprintf("Drug cells at 1:%d: %s across %d conditions", DILUTION,
                format(nrow(drug), big.mark = ","), uniqueN(drug$cond)))
print(drug[, .N, by = cond][order(cond)])

# ---- load library -----------------------------------------------------------
lib_prof <- as.matrix(data.frame(fread(lib_prof_file), row.names = 1)[, FEATURES])
lib_files <- file.path(lib_dir, c("all_features.csv", "all_features 2.csv"))
nt_cells <- rbindlist(lapply(lib_files, function(f) {
  x <- fread(f, select = c("well", FEATURES))
  x[, strain := tstrsplit(well, "__")[[3]]]
  x[strain %in% NT_STRAINS]
}))
message(sprintf("Library: %d strain profiles; %s pooled NT cells",
                nrow(lib_prof), format(nrow(nt_cells), big.mark = ",")))
stopifnot(all(NT_STRAINS %in% rownames(lib_prof)))

Xd <- as.matrix(drug[, ..FEATURES])
Xnt <- as.matrix(nt_cells[, ..FEATURES])

# Per-cell z for each anchor + matching library profile matrix
anchors <- list(
  libNT = list(Z = robust_z(Xd, Xnt), lib = lib_prof,
               label = "library NT anchor"),
  plate = list(Z = robust_z(Xd, Xd),
               lib = sweep(lib_prof, 2, apply(lib_prof, 2, median)),
               label = "drug-plate anchor"))

conds <- sort(unique(drug$cond))
cond_rows <- split(seq_len(nrow(drug)), drug$cond)[conds]
fov_rows  <- lapply(cond_rows, function(r) split(r, drug$fov[r]))
lib_strains <- rownames(lib_prof)

drug_profiles <- function(Z, rows_by_cond) {
  t(vapply(rows_by_cond, function(r) col_medians(Z[r, , drop = FALSE]), numeric(ncol(Z))))
}

neighbour_table <- function(D) {
  rbindlist(lapply(conds, function(cnd) {
    d <- sort(D[cnd, lib_strains])
    data.table(cond = cnd, rank = seq_along(d), neighbour = names(d), distance = unname(d))
  }))
}

# ---- run each anchor --------------------------------------------------------
res <- list()
for (an in names(anchors)) {
  A <- anchors[[an]]
  dprof <- drug_profiles(A$Z, cond_rows)
  rownames(dprof) <- conds
  colnames(dprof) <- FEATURES
  D <- profile_dist(rbind(A$lib, dprof))
  nb <- neighbour_table(D)

  # phenotype strength: mean distance to NT1-4 vs typical NT-NT distance
  nt_nt <- D[NT_STRAINS, NT_STRAINS]
  nt_ref <- mean(nt_nt[upper.tri(nt_nt)])
  strength <- data.table(cond = conds,
                         dist_to_NT = rowMeans(D[conds, NT_STRAINS, drop = FALSE]),
                         typical_NT_NT = nt_ref)

  # FOV bootstrap: how often each gene lands in a condition's top TOP_STABLE
  boot_hits <- rbindlist(lapply(seq_len(N_BOOT), function(b) {
    rows_b <- lapply(fov_rows, function(fr) unlist(sample(fr, length(fr), replace = TRUE),
                                                    use.names = FALSE))
    dp <- drug_profiles(A$Z, rows_b); rownames(dp) <- conds
    Db <- profile_dist(rbind(A$lib, dp))
    rbindlist(lapply(conds, function(cnd)
      data.table(cond = cnd, neighbour = names(sort(Db[cnd, lib_strains]))[seq_len(TOP_STABLE)])))
  }))
  stab <- boot_hits[, .(boot_top5 = .N / N_BOOT), by = .(cond, neighbour)]
  nb <- merge(nb, stab, by = c("cond", "neighbour"), all.x = TRUE)
  nb[is.na(boot_top5), boot_top5 := 0]
  setorder(nb, cond, rank)
  nb[, anchor := an]

  res[[an]] <- list(prof = dprof, D = D, nb = nb, strength = strength[, anchor := an])
  fwrite(dprof, file.path(out_dir, sprintf("drug_profiles_%s.csv", an)), row.names = TRUE)
}

nb_all <- rbindlist(lapply(res, `[[`, "nb"))
nb_all[, c("drug", "conc") := tstrsplit(cond, "_")]
nb_all[, group := unname(gene_group[neighbour])]
nb_all[, super := factor(unname(super_map[group]), levels = super_levels)]
fwrite(nb_all, file.path(out_dir, "drug_neighbours_all_ranks.csv"))
str_all <- rbindlist(lapply(res, `[[`, "strength"))

# ---- cross-anchor agreement -------------------------------------------------
agree <- nb_all[rank <= TOPK, .(set = list(neighbour)), by = .(cond, anchor)]
agree <- dcast(agree, cond ~ anchor, value.var = "set")
agree[, shared := mapply(function(a, b) paste(intersect(a, b), collapse = ", "), libNT, plate)]
agree[, n_shared := mapply(function(a, b) length(intersect(a, b)), libNT, plate)]
agree[, jaccard := mapply(function(a, b) length(intersect(a, b)) / length(union(a, b)), libNT, plate)]
agree[, `:=`(libNT = NULL, plate = NULL)]

# Rank correlation of the full 67-strain distance vectors between anchors
agree[, spearman_all := vapply(cond, function(cnd)
  cor(res$libNT$D[cnd, lib_strains], res$plate$D[cnd, lib_strains], method = "spearman"), 0)]

# ---- expected-pathway check (AUC + permutation null) ------------------------
auc <- function(d, member) {
  r <- rank(-d)   # higher rank = closer
  n1 <- sum(member); n0 <- sum(!member)
  (sum(r[member]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
pw <- rbindlist(lapply(names(res), function(an) {
  D <- res[[an]]$D
  rbindlist(lapply(conds, function(cnd) {
    drg <- sub("_.*", "", cnd)
    ex <- EXPECTED[drug == drg]
    d <- D[cnd, lib_strains]
    member <- lib_strains %in% FUNCTIONAL_GROUPS[[ex$pathway]]
    obs <- auc(d, member)
    null <- replicate(N_PERM, auc(d, sample(member)))
    tg <- if (is.na(ex$target)) NA_character_ else {
      tt <- strsplit(ex$target, ",")[[1]]
      rk <- rank(d)[match(tt, lib_strains)]
      paste(sprintf("%s #%d", tt, as.integer(rk)), collapse = ", ")
    }
    data.table(anchor = an, cond = cnd, drug = drg, expected_pathway = ex$pathway,
               pathway_auc = obs, perm_p = (sum(null >= obs) + 1) / (N_PERM + 1),
               target_rank = tg)
  }))
}))
pw[, p_bh := p.adjust(perm_p, method = "BH"), by = anchor]
fwrite(pw, file.path(out_dir, "expected_pathway_check.csv"))

# ---- summary table ----------------------------------------------------------
top_str <- function(an) nb_all[anchor == an & rank <= TOPK,
  .(top = paste(sprintf("%s(%.0f%%)", neighbour, 100 * boot_top5), collapse = ", ")), by = cond]
summ <- merge(top_str("libNT"), top_str("plate"), by = "cond", suffixes = c("_libNT", "_plate"))
summ <- merge(summ, agree, by = "cond")
summ <- merge(summ, dcast(str_all, cond ~ anchor, value.var = "dist_to_NT"), by = "cond")
setnames(summ, c("libNT", "plate"), c("dist_to_NT_libNT", "dist_to_NT_plate"))
fwrite(summ, file.path(out_dir, "drug_neighbours_summary.csv"))

message("\nTop-", TOPK, " neighbours (bootstrap top-5 frequency):")
for (cnd in conds) {
  s <- summ[cond == cnd]
  message(sprintf("\n%s  [shared %d/%d, Jaccard %.2f, Spearman %.2f]",
                  cnd, s$n_shared, TOPK, s$jaccard, s$spearman_all))
  message("  libNT: ", s$top_libNT)
  message("  plate: ", s$top_plate)
}
message("\nExpected-pathway check:")
print(pw[, .(anchor, cond, expected_pathway, auc = round(pathway_auc, 2),
             p = signif(perm_p, 2), p_bh = signif(p_bh, 2), target_rank)])
message("\nPhenotype strength (distance to NT; typical NT-NT in brackets):")
print(str_all[, .(cond, anchor, dist_to_NT = round(dist_to_NT, 2), typical_NT_NT = round(typical_NT_NT, 2))])

# ---- figures ----------------------------------------------------------------
cond_lab <- function(cnd) sub("_", " ", cnd)
nb_plot <- nb_all[rank <= TOPK]
# order rows: drug blocks, doses low -> high top to bottom
row_levels <- unlist(lapply(c("EMB", "INH", "MOX", "RIF"), function(dg)
  paste(dg, conc_levels)))
nb_plot[, cond_f := factor(cond_lab(cond), levels = rev(row_levels))]
nb_plot[, txt_col := ifelse(super %in% c("Control (NT)", "DNA replication / repair"), "white", "black")]

ladder <- function(an, title, legend = TRUE) {
  dd <- nb_plot[anchor == an]
  ggplot(dd, aes(rank, cond_f)) +
    geom_tile(aes(fill = super, alpha = 0.25 + 0.75 * boot_top5),
              colour = "white", linewidth = 0.8) +
    geom_text(aes(label = neighbour, colour = txt_col), size = 2.9, fontface = "italic") +
    scale_fill_manual(values = super_cols, drop = FALSE, name = NULL) +
    scale_colour_identity() +
    scale_alpha_identity() +
    scale_x_continuous(breaks = seq_len(TOPK), expand = c(0, 0)) +
    geom_hline(yintercept = c(3.5, 6.5, 9.5), colour = "grey40", linewidth = 0.4) +
    labs(x = "Neighbour rank (1 = closest)", y = NULL, title = title,
         subtitle = "Tile opacity = FOV-bootstrap frequency in the top 5") +
    theme_Publication(base_size = 11) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.line = element_blank(), axis.ticks = element_blank(),
          plot.subtitle = element_text(size = 8)) +
    guides(fill = guide_legend(nrow = 2, override.aes = list(alpha = 1))) +
    if (!legend) theme(legend.position = "none") else NULL
}
p1 <- (ladder("libNT", "Library-NT anchor", legend = FALSE) |
       ladder("plate", "Drug-plate anchor")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")
save_fig(p1, "01_nearest_knockdowns_by_anchor", 15, 6.8)

# Expected-pathway AUC per drug x dose, both anchors
pw[, conc := factor(sub(".*_", "", cond), levels = conc_levels)]
pw[, anchor_lab := factor(ifelse(anchor == "libNT", "Library-NT", "Drug-plate"),
                          levels = c("Library-NT", "Drug-plate"))]
pw[, drug_lab := paste0(drug, " -> ", sub(" /.*|/.*", "", expected_pathway))]
p2 <- ggplot(pw, aes(conc, pathway_auc, colour = anchor_lab, group = anchor_lab)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = p_bh < 0.05), size = 3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c(`FALSE` = "BH q >= 0.05", `TRUE` = "BH q < 0.05")) +
  scale_colour_Publication() +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ drug_lab, nrow = 1) +
  labs(x = "Dose", y = "Expected-pathway AUC",
       title = "Do drugs land near their target pathway?",
       subtitle = "AUC = P(pathway knockdown closer than a non-pathway strain); 0.5 = chance; label-permutation null") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))
save_fig(p2, "02_expected_pathway_auc", 11, 4.2)

# Phenotype strength + anchor agreement
str_all[, conc := factor(sub(".*_", "", cond), levels = conc_levels)]
str_all[, drug := sub("_.*", "", cond)]
str_all[, anchor_lab := factor(ifelse(anchor == "libNT", "Library-NT", "Drug-plate"),
                               levels = c("Library-NT", "Drug-plate"))]
p3a <- ggplot(str_all, aes(conc, dist_to_NT / typical_NT_NT, colour = drug, group = drug)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  scale_colour_Publication() +
  facet_wrap(~ anchor_lab) +
  labs(x = "Dose", y = "Distance to NT / typical NT-NT", title = "Phenotype strength",
       subtitle = "1 = indistinguishable from NT-to-NT spread") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))
agree[, conc := factor(sub(".*_", "", cond), levels = conc_levels)]
agree[, drug := sub("_.*", "", cond)]
p3b <- ggplot(agree, aes(conc, n_shared, fill = drug)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_Publication() +
  scale_y_continuous(limits = c(0, TOPK), breaks = 0:TOPK) +
  labs(x = "Dose", y = sprintf("Top-%d neighbours shared", TOPK),
       title = "Anchor agreement",
       subtitle = "Same knockdown in both anchors' top 8") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))
p3 <- (p3a | p3b) + plot_layout(widths = c(2, 1.2))
save_fig(p3, "03_strength_and_anchor_agreement", 12, 4.6)

message("\nDone. Tables in ", out_dir, " | figures in ", fig_dir)
