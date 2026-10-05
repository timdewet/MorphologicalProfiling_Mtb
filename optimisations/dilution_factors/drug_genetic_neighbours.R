# =============================================================================
# Drug-treated iniB-mScarlet cells vs the iniB CRISPRi library:
# nearest genetic (knockdown) neighbours for each drug x dose
#
# Drug data : Plate60 (EMB, INH, MOX, RIF x 0.25/0.5/1x MIC). Each drug x dose
#             was plated at 5 dilutions of the same fixed resuspension; dilution
#             changes cell density but not morphology (dilution_optimisation.R,
#             fig 03), so the dilution wells are treated as replicate wells.
# Library   : 67-strain iniB-mScarlet CRISPRi panel (63 KDs + IFT, SSB, NT1-4),
#             canonical per-batch NT-corrected profiles (output_inib/
#             strain_profiles_corrected.csv) + raw NT cells.
#
# Method mirrors the library pipeline (inib_morphology_analysis.py /
# ribf_vs_inib.py): 36 shape features -> per-cell robust z vs a reference ->
# per-well median profile -> condition profile = mean of its wells (replicates
# weighted equally) -> standardise features across library + condition profiles
# -> PCA (12 PCs) -> Euclidean distance.  Individual wells are projected into
# that same space.  Shape only: reporter intensity is not comparable across runs.
#
# The drug plate is a SEPARATE imaging run with NO untreated / NT control, so
# where "untreated" sits in this run is unknown.  Two anchors bracket that:
#   * libNT : drug cells robust-z vs the library's pooled NT1-4 cells (assumes no
#             run-to-run shift); library = canonical per-batch corrected profiles.
#   * plate : drug cells robust-z vs all analysed drug-plate cells; library
#             profiles re-centred on the library median (assumes the average drug
#             condition ~ the average knockdown).  Removes any run-level shift.
# Neighbours that hold under BOTH anchors are the defensible ones.
#
# Uncertainty: two-level bootstrap (wells with replacement, then FOVs within
# each well).  With several wells per condition we also report how many
# replicate wells independently put a gene in their top 5, a per-well
# expected-pathway AUC, and whether wells cluster by condition or by dilution.
# Dilution wells share one fixed suspension, so they are technical (plating /
# imaging) replicates, not independent cultures.
#
# Sanity check built in: the library carries the direct targets of three drugs
# (EMB -> embA, INH -> inhA, MOX -> gyrA/gyrB; no rpoB for RIF), so each drug's
# expected pathway is scored by AUC against a label-permutation null.
#
# Run from the project root:
#   Rscript optimisations/dilution_factors/drug_genetic_neighbours.R \
#     [drug_features.csv] [library_dir] [library_profiles.csv] [dilutions]
# dilutions = "all" (default) or a comma list, e.g. "50".
# Outputs -> <drug_features dir>/drug_neighbours/<all_dilutions | 1in50 ...>/
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
dil_arg   <- if (length(args) >= 4) args[4] else "all"

ALL_DILUTIONS <- c(20L, 50L, 100L, 150L, 200L)
DILUTIONS <- if (dil_arg == "all") ALL_DILUTIONS else as.integer(strsplit(dil_arg, ",")[[1]])
run_tag <- if (dil_arg == "all") "all_dilutions" else paste0("1in", paste(DILUTIONS, collapse = "_"))

out_dir <- file.path(dirname(normalizePath(drug_file, mustWork = TRUE)), "drug_neighbours", run_tag)
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
TOPK       <- 8L      # neighbours reported per drug x dose
TOP_STABLE <- 5L      # stability = how often a gene is in the top 5
N_BOOT     <- 500L
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

# Library build_network space: standardise each feature across the fitting
# profiles (constant features -> scale 1), PCA, keep PROFILE_PCA PCs.  Returns
# a projector so replicate wells land in the same space as the profiles.
fit_space <- function(prof) {
  ctr <- colMeans(prof)
  sdv <- apply(prof, 2, sd); sdv <- ifelse(sdv > 1e-12, sdv, 1)
  P <- sweep(sweep(prof, 2, ctr), 2, sdv, "/")
  k <- max(2L, min(PROFILE_PCA, nrow(P) - 1L, ncol(P)))
  rot <- prcomp(P, center = FALSE)$rotation[, seq_len(k), drop = FALSE]
  function(X) sweep(sweep(X, 2, ctr), 2, sdv, "/") %*% rot
}

cross_dist <- function(A, B) {   # Euclidean distances, rows of A x rows of B
  sqrt(pmax(outer(rowSums(A^2), rowSums(B^2), "+") - 2 * A %*% t(B), 0))
}

auc <- function(d, member) {      # P(member closer than non-member)
  r <- rank(-d)
  n1 <- sum(member); n0 <- sum(!member)
  (sum(r[member]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ---- load drug cells --------------------------------------------------------
drug <- fread(drug_file, select = c("well", "fov", FEATURES))
drug[, c("dc", "reporter", "dil") := tstrsplit(well, "__")]
drug[, c("drug", "conc") := tstrsplit(dc, "_")]
drug[, dil := as.integer(sub("1in", "", sub("_focused", "", dil)))]
drug <- drug[dil %in% DILUTIONS]
drug[, cond := paste(drug, conc, sep = "_")]
drug[, rep_id := paste(cond, dil, sep = "|")]
message(sprintf("Drug cells (dilutions %s): %s, %d conditions x %d wells",
                paste(DILUTIONS, collapse = "/"), format(nrow(drug), big.mark = ","),
                uniqueN(drug$cond), uniqueN(drug$rep_id)))

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
lib_strains <- rownames(lib_prof)

Xd  <- as.matrix(drug[, ..FEATURES])
Xnt <- as.matrix(nt_cells[, ..FEATURES])

anchors <- list(
  libNT = list(Z = robust_z(Xd, Xnt), lib = lib_prof),
  plate = list(Z = robust_z(Xd, Xd), lib = sweep(lib_prof, 2, apply(lib_prof, 2, median))))

# Row indices: condition -> well -> FOV
conds <- sort(unique(drug$cond))
reps  <- unique(drug[, .(rep_id, cond, dil)])[order(cond, dil)]
rep_rows <- split(seq_len(nrow(drug)), drug$rep_id)[reps$rep_id]
rep_fov_rows <- lapply(rep_rows, function(r) split(r, drug$fov[r]))
cond_reps <- split(reps$rep_id, reps$cond)[conds]
multi_rep <- all(lengths(cond_reps) >= 2)

well_profiles <- function(Z, rows) {
  M <- t(vapply(rows, function(r) col_medians(Z[r, , drop = FALSE]), numeric(ncol(Z))))
  colnames(M) <- FEATURES; M
}
cond_profiles <- function(W, reps_by_cond) {
  M <- t(vapply(reps_by_cond, function(rr) colMeans(W[rr, , drop = FALSE]), numeric(ncol(W))))
  rownames(M) <- names(reps_by_cond); colnames(M) <- FEATURES; M
}

# ---- run each anchor --------------------------------------------------------
res <- list()
for (an in names(anchors)) {
  A <- anchors[[an]]
  W  <- well_profiles(A$Z, rep_rows)            # wells x features
  Cp <- cond_profiles(W, cond_reps)             # conditions x features
  proj <- fit_space(rbind(A$lib, Cp))
  L  <- proj(A$lib); Cs <- proj(Cp); Ws <- proj(W)

  Dcl <- cross_dist(Cs, L); dimnames(Dcl) <- list(conds, lib_strains)
  Dwl <- cross_dist(Ws, L); dimnames(Dwl) <- list(reps$rep_id, lib_strains)

  nb <- rbindlist(lapply(conds, function(cnd) {
    d <- sort(Dcl[cnd, ])
    data.table(cond = cnd, rank = seq_along(d), neighbour = names(d), distance = unname(d))
  }))

  # Replicate support: in how many wells is the gene in that well's own top 5?
  rep_top <- rbindlist(lapply(reps$rep_id, function(r)
    data.table(rep_id = r, neighbour = names(sort(Dwl[r, ]))[seq_len(TOP_STABLE)])))
  rep_top <- merge(rep_top, reps[, .(rep_id, cond)], by = "rep_id")
  support <- rep_top[, .(wells_top5 = .N), by = .(cond, neighbour)]

  # Two-level bootstrap: wells with replacement, then FOVs within each well
  boot_hits <- rbindlist(lapply(seq_len(N_BOOT), function(b) {
    Cb <- t(vapply(conds, function(cnd) {
      rr <- cond_reps[[cnd]]
      rr <- rr[sample.int(length(rr), length(rr), replace = TRUE)]
      wp <- vapply(rr, function(r) {
        fr <- rep_fov_rows[[r]]
        col_medians(A$Z[unlist(fr[sample.int(length(fr), length(fr), replace = TRUE)],
                               use.names = FALSE), , drop = FALSE])
      }, numeric(length(FEATURES)))
      rowMeans(wp)
    }, numeric(length(FEATURES))))
    colnames(Cb) <- FEATURES
    pb <- fit_space(rbind(A$lib, Cb))
    Db <- cross_dist(pb(Cb), pb(A$lib)); dimnames(Db) <- list(conds, lib_strains)
    rbindlist(lapply(conds, function(cnd)
      data.table(cond = cnd, neighbour = names(sort(Db[cnd, ]))[seq_len(TOP_STABLE)])))
  }))
  stab <- boot_hits[, .(boot_top5 = .N / N_BOOT), by = .(cond, neighbour)]

  nb <- merge(nb, stab, by = c("cond", "neighbour"), all.x = TRUE)
  nb <- merge(nb, support, by = c("cond", "neighbour"), all.x = TRUE)
  nb[is.na(boot_top5), boot_top5 := 0]
  nb[is.na(wells_top5), wells_top5 := 0L]
  nb[, n_wells := lengths(cond_reps)[cond]]
  setorder(nb, cond, rank)
  nb[, anchor := an]

  # Phenotype strength vs two noise floors: library NT-NT and replicate wells
  Dnn <- cross_dist(L[NT_STRAINS, ], L[NT_STRAINS, ])
  nt_ref <- mean(Dnn[upper.tri(Dnn)])
  rep_spread <- vapply(conds, function(cnd) {
    X <- Ws[cond_reps[[cnd]], , drop = FALSE]
    if (nrow(X) < 2) return(NA_real_)
    Dx <- cross_dist(X, X); mean(Dx[upper.tri(Dx)])
  }, 0)
  strength <- data.table(cond = conds, anchor = an,
                         dist_to_NT = rowMeans(Dcl[, NT_STRAINS, drop = FALSE]),
                         typical_NT_NT = nt_ref,
                         nn_distance = apply(Dcl, 1, min),
                         replicate_spread = rep_spread)

  res[[an]] <- list(Dcl = Dcl, Dwl = Dwl, Ws = Ws, Cs = Cs, nb = nb, strength = strength)
  fwrite(data.table(cond = conds, Cp), file.path(out_dir, sprintf("drug_profiles_%s.csv", an)))
}

nb_all <- rbindlist(lapply(res, `[[`, "nb"))
nb_all[, c("drug", "conc") := tstrsplit(cond, "_")]
nb_all[, group := unname(gene_group[neighbour])]
nb_all[, super := factor(unname(super_map[group]), levels = super_levels)]
fwrite(nb_all, file.path(out_dir, "drug_neighbours_all_ranks.csv"))
str_all <- rbindlist(lapply(res, `[[`, "strength"))
fwrite(str_all, file.path(out_dir, "phenotype_strength.csv"))

# ---- cross-anchor agreement -------------------------------------------------
agree <- nb_all[rank <= TOPK, .(set = list(neighbour)), by = .(cond, anchor)]
agree <- dcast(agree, cond ~ anchor, value.var = "set")
agree[, shared := mapply(function(a, b) paste(intersect(a, b), collapse = ", "), libNT, plate)]
agree[, n_shared := mapply(function(a, b) length(intersect(a, b)), libNT, plate)]
agree[, jaccard := mapply(function(a, b) length(intersect(a, b)) / length(union(a, b)), libNT, plate)]
agree[, `:=`(libNT = NULL, plate = NULL)]
agree[, spearman_all := vapply(cond, function(cnd)
  cor(res$libNT$Dcl[cnd, ], res$plate$Dcl[cnd, ], method = "spearman"), 0)]

# ---- expected-pathway check -------------------------------------------------
# Pooled profile: AUC vs label-permutation null.  Per well: one AUC per
# replicate, summarised by how many wells exceed chance (0.5).
pw <- rbindlist(lapply(names(res), function(an) {
  rbindlist(lapply(conds, function(cnd) {
    drg <- sub("_.*", "", cnd)
    ex <- EXPECTED[drug == drg]
    member <- lib_strains %in% FUNCTIONAL_GROUPS[[ex$pathway]]
    d <- res[[an]]$Dcl[cnd, ]
    obs <- auc(d, member)
    null <- replicate(N_PERM, auc(d, sample(member)))
    well_auc <- vapply(cond_reps[[cnd]], function(r) auc(res[[an]]$Dwl[r, ], member), 0)
    tg <- if (is.na(ex$target)) NA_character_ else {
      tt <- strsplit(ex$target, ",")[[1]]
      paste(sprintf("%s #%d", tt, as.integer(rank(d)[match(tt, lib_strains)])), collapse = ", ")
    }
    data.table(anchor = an, cond = cnd, drug = drg, expected_pathway = ex$pathway,
               pathway_auc = obs, perm_p = (sum(null >= obs) + 1) / (N_PERM + 1),
               wells_above_chance = sum(well_auc > 0.5), n_wells = length(well_auc),
               well_auc_min = min(well_auc), well_auc_max = max(well_auc),
               target_rank = tg)
  }))
}))
pw[, p_bh := p.adjust(perm_p, method = "BH"), by = anchor]
fwrite(pw, file.path(out_dir, "expected_pathway_check.csv"))

well_auc_long <- rbindlist(lapply(names(res), function(an) {
  rbindlist(lapply(reps$rep_id, function(r) {
    cnd <- reps[rep_id == r, cond]
    ex <- EXPECTED[drug == sub("_.*", "", cnd)]
    data.table(anchor = an, rep_id = r, cond = cnd, dil = reps[rep_id == r, dil],
               auc = auc(res[[an]]$Dwl[r, ], lib_strains %in% FUNCTIONAL_GROUPS[[ex$pathway]]))
  }))
}))

# ---- replicate structure: do wells group by condition or by dilution? -------
if (multi_rep) {
  rep_qc <- rbindlist(lapply(names(res), function(an) {
    Ws <- res[[an]]$Ws
    Dww <- cross_dist(Ws, Ws); diag(Dww) <- Inf
    nn <- apply(Dww, 1, which.min)
    same_cond <- mean(reps$cond[nn] == reps$cond)
    same_dil  <- mean(reps$dil[nn] == reps$dil)
    # chance levels for nearest-other-well matching
    data.table(anchor = an, nn_same_condition = same_cond,
               chance_condition = mean((lengths(cond_reps)[reps$cond] - 1) / (nrow(reps) - 1)),
               nn_same_dilution = same_dil,
               chance_dilution = mean((table(reps$dil)[as.character(reps$dil)] - 1) / (nrow(reps) - 1)))
  }))
  fwrite(rep_qc, file.path(out_dir, "replicate_structure_qc.csv"))
  message("\nReplicate structure (nearest other well shares ...):")
  print(rep_qc, digits = 2)
}

# ---- summary table ----------------------------------------------------------
top_str <- function(an) nb_all[anchor == an & rank <= TOPK,
  .(top = paste(sprintf("%s(%d/%d,%.0f%%)", neighbour, wells_top5, n_wells, 100 * boot_top5),
                collapse = ", ")), by = cond]
summ <- merge(top_str("libNT"), top_str("plate"), by = "cond", suffixes = c("_libNT", "_plate"))
summ <- merge(summ, agree, by = "cond")
summ <- merge(summ, dcast(str_all, cond ~ anchor, value.var = "dist_to_NT"), by = "cond")
setnames(summ, c("libNT", "plate"), c("dist_to_NT_libNT", "dist_to_NT_plate"))
fwrite(summ, file.path(out_dir, "drug_neighbours_summary.csv"))

message("\nTop-", TOPK, " neighbours: gene(replicate wells with it in top 5, bootstrap top-5 %):")
for (cnd in conds) {
  s <- summ[cond == cnd]
  message(sprintf("\n%s  [shared %d/%d, Spearman %.2f]", cnd, s$n_shared, TOPK, s$spearman_all))
  message("  libNT: ", s$top_libNT)
  message("  plate: ", s$top_plate)
}
message("\nExpected-pathway check:")
print(pw[, .(anchor, cond, auc = round(pathway_auc, 2), p = signif(perm_p, 2),
             q = signif(p_bh, 2), wells = paste0(wells_above_chance, "/", n_wells),
             well_range = sprintf("%.2f-%.2f", well_auc_min, well_auc_max), target_rank)])
message("\nPhenotype strength:")
print(str_all[, .(cond, anchor, to_NT = round(dist_to_NT / typical_NT_NT, 1),
                  nn_dist = round(nn_distance, 2), rep_spread = round(replicate_spread, 2))])

# ---- figures ----------------------------------------------------------------
dil_txt <- if (multi_rep) sprintf("dilutions %s pooled as replicate wells",
                                  paste0("1:", DILUTIONS, collapse = "/")) else
  sprintf("1:%s dilution only", paste(DILUTIONS, collapse = "/"))
row_levels <- unlist(lapply(c("EMB", "INH", "MOX", "RIF"), function(dg) paste(dg, conc_levels)))
anchor_labs <- c(libNT = "Library-NT", plate = "Drug-plate")

nb_plot <- nb_all[rank <= TOPK]
nb_plot[, cond_f := factor(sub("_", " ", cond), levels = rev(row_levels))]
nb_plot[, txt_col := ifelse(super %in% c("Control (NT)", "DNA replication / repair"), "white", "black")]
nb_plot[, lab := if (multi_rep) sprintf("%s\n%d/%d", neighbour, wells_top5, n_wells) else neighbour]

ladder <- function(an, title, legend = TRUE) {
  ggplot(nb_plot[anchor == an], aes(rank, cond_f)) +
    geom_tile(aes(fill = super, alpha = 0.25 + 0.75 * boot_top5),
              colour = "white", linewidth = 0.8) +
    geom_text(aes(label = lab, colour = txt_col), size = 2.6, fontface = "italic",
              lineheight = 0.85) +
    scale_fill_manual(values = super_cols, drop = FALSE, name = NULL) +
    scale_colour_identity() +
    scale_alpha_identity() +
    scale_x_continuous(breaks = seq_len(TOPK), expand = c(0, 0)) +
    geom_hline(yintercept = c(3.5, 6.5, 9.5), colour = "grey40", linewidth = 0.4) +
    labs(x = "Neighbour rank (1 = closest)", y = NULL, title = title,
         subtitle = if (multi_rep)
           "Opacity = bootstrap top-5 frequency; n/5 = replicate wells with the gene in their own top 5" else
           "Opacity = FOV-bootstrap top-5 frequency") +
    theme_Publication(base_size = 11) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.line = element_blank(), axis.ticks = element_blank(),
          plot.subtitle = element_text(size = 8)) +
    guides(fill = guide_legend(nrow = 2, override.aes = list(alpha = 1))) +
    if (!legend) theme(legend.position = "none") else NULL
}
p1 <- (ladder("libNT", "Library-NT anchor", legend = FALSE) |
       ladder("plate", "Drug-plate anchor")) +
  plot_annotation(caption = dil_txt)
save_fig(p1, "01_nearest_knockdowns_by_anchor", 15, if (multi_rep) 7.8 else 6.8)

# Expected-pathway AUC: pooled profile (line) + each replicate well (points)
pw[, conc := factor(sub(".*_", "", cond), levels = conc_levels)]
pw[, anchor_lab := factor(anchor_labs[anchor], levels = anchor_labs)]
pw[, drug_lab := paste0(drug, " -> ", sub(" /.*|/.*", "", expected_pathway))]
well_auc_long[, `:=`(conc = factor(sub(".*_", "", cond), levels = conc_levels),
                     anchor_lab = factor(anchor_labs[anchor], levels = anchor_labs),
                     drug_lab = pw$drug_lab[match(cond, pw$cond)])]
p2 <- ggplot(pw, aes(conc, pathway_auc, colour = anchor_lab, group = anchor_lab)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  { if (multi_rep) geom_point(data = well_auc_long, aes(conc, auc, colour = anchor_lab),
                              position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.5, seed = 1),
                              size = 1.1, alpha = 0.45, inherit.aes = FALSE) } +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.5)) +
  geom_point(aes(shape = p_bh < 0.05), size = 3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c(`FALSE` = "BH q >= 0.05", `TRUE` = "BH q < 0.05")) +
  scale_colour_Publication() +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ drug_lab, nrow = 1) +
  labs(x = "Dose", y = "Expected-pathway AUC",
       title = "Do drugs land near their target pathway?",
       subtitle = paste0("AUC = P(pathway knockdown closer than a non-pathway strain); 0.5 = chance. ",
                         "Large points = pooled profile (permutation q); ",
                         if (multi_rep) "small = individual replicate wells" else ""),
       caption = dil_txt) +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))
save_fig(p2, "02_expected_pathway_auc", 11, 4.4)

# Phenotype strength and anchor agreement
str_all[, `:=`(conc = factor(sub(".*_", "", cond), levels = conc_levels),
               drug = sub("_.*", "", cond),
               anchor_lab = factor(anchor_labs[anchor], levels = anchor_labs))]
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
agree[, `:=`(conc = factor(sub(".*_", "", cond), levels = conc_levels),
             drug = sub("_.*", "", cond))]
p3b <- ggplot(agree, aes(conc, n_shared, fill = drug)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single"), width = 0.75) +
  scale_fill_Publication() +
  scale_y_continuous(limits = c(0, TOPK), breaks = 0:TOPK) +
  labs(x = "Dose", y = sprintf("Top-%d neighbours shared", TOPK),
       title = "Anchor agreement", subtitle = "Same knockdown in both anchors' top 8") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))
p3 <- (p3a | p3b) + plot_layout(widths = c(2, 1.2)) + plot_annotation(caption = dil_txt)
save_fig(p3, "03_strength_and_anchor_agreement", 12, 4.8)

# Replicate wells in profile space: do they group by drug x dose?
if (multi_rep) {
  emb <- rbindlist(lapply(names(res), function(an) {
    X <- rbind(res[[an]]$Ws, res[[an]]$Cs)
    pc <- prcomp(X, center = TRUE)$x[, 1:2]
    data.table(anchor_lab = factor(anchor_labs[an], levels = anchor_labs),
               id = rownames(X), PC1 = pc[, 1], PC2 = pc[, 2],
               type = rep(c("well", "pooled"), c(nrow(res[[an]]$Ws), nrow(res[[an]]$Cs))))
  }))
  emb[, cond := sub("\\|.*", "", id)]
  emb[, `:=`(drug = sub("_.*", "", cond),
             conc = factor(sub(".*_", "", cond), levels = conc_levels),
             dil = as.integer(ifelse(type == "well", sub(".*\\|", "", id), NA)))]
  p4 <- ggplot(emb, aes(PC1, PC2, colour = drug)) +
    geom_line(data = emb[type == "well"], aes(group = cond), alpha = 0.25) +
    geom_point(data = emb[type == "well"], aes(shape = conc), size = 1.8, alpha = 0.7) +
    geom_point(data = emb[type == "pooled"], aes(shape = conc), size = 4.2, stroke = 1.1) +
    scale_colour_Publication() +
    facet_wrap(~ anchor_lab, scales = "free") +
    labs(title = "Replicate wells vs pooled profiles",
         subtitle = "Small = individual dilution wells; large = condition mean. PCA of the drug profiles in the library-defined space") +
    theme_Publication(base_size = 11) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          plot.subtitle = element_text(size = 8))
  save_fig(p4, "04_replicate_wells_pca", 11, 5.2)
}

message("\nDone. Tables in ", out_dir, " | figures in ", fig_dir)
