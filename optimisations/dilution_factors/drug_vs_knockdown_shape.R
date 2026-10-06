# =============================================================================
# Interpretable shape changes: drugs vs CRISPRi knockdowns
#
# The 36-feature neighbour distance (drug_genetic_neighbours.R) is dominated by
# the many size/length-correlated features, which hides WHICH axis drives a
# (mis)match.  Here each condition is reduced to % change in width, length and
# area plus reporter fold-change, each against its own reference:
#   * drugs      : vs RIF 0.25xMIC on the same plate and day (the condition that
#                  sits among the NT controls in the drug-plate anchor; no
#                  untreated wells exist on these plates)
#   * knockdowns : vs the NT controls of the knockdown's own library batch
# Per-well medians; drugs = mean of the 5 dilution wells.  Reporter fold-change
# between runs is descriptive only (both references sit at background).
#
# Run from the project root:
#   Rscript optimisations/dilution_factors/drug_vs_knockdown_shape.R \
#     [data_dir] [library_dir] [day2,day4]
# Outputs -> <data_dir>/shape_vs_knockdowns/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

if (file.exists("Theme.R")) source("Theme.R") else
  stop("Theme.R not found - run from project root")

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "optimisations/dilution_factors"
lib_dir  <- if (length(args) >= 2) args[2] else "input_data/iniB_allstrains"
days     <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else c("day2", "day4")

out_dir <- file.path(data_dir, "shape_vs_knockdowns")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(p, name, w, h, dpi = 300) {
  # cairo is unavailable without XQuartz; quartz (pdf) + ragg (png) instead
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white",
         device = function(filename, ...) grDevices::quartz(type = "pdf", file = filename, ...))
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = dpi, bg = "white", device = ragg::agg_png)
}

FEATS <- c(width = "width_median_um", length = "length_um", area = "area_um2",
           reporter = "intensity_median_mCher")
DRUG_REF   <- "RIF_0.25xMIC"
NT_STRAINS <- c("NT1", "NT2", "NT3", "NT4")
conc_levels <- c("0.25xMIC", "0.5xMIC", "1xMIC")

KD_GROUPS <- list(
  "Mycolic acid / FAS-II"  = c("inhA", "kasA", "hadA", "hadB", "hadC", "fabD",
                               "acpM", "fad32", "accA3", "pksB", "mmpL3"),
  "Arabinogalactan"        = c("embA", "dprE1", "glfT2"),
  "DNA replication/repair" = c("gyrA", "gyrB", "dnaN", "dnaE1", "dnaZX", "polA",
                               "topA", "SSB", "hupB"),
  "Control (NT)"           = NT_STRAINS)
kd_group <- setNames(rep(names(KD_GROUPS), lengths(KD_GROUPS)), unlist(KD_GROUPS))
group_cols <- c("Mycolic acid / FAS-II" = "#D55E00", "Arabinogalactan" = "#E69F00",
                "DNA replication/repair" = "#0072B2", "Control (NT)" = "#000000",
                "Other knockdown" = "grey75")
LABEL_KD <- c("inhA", "kasA", "hadA", "pksB", "mmpL3", "fad32",
              "embA", "dprE1", "glfT2", "gyrB", "dnaE1")

pct <- function(x, ref) 100 * (x / ref - 1)

# ---- library: per strain x batch medians vs own-batch NT --------------------
bm <- fread(file.path(lib_dir, "bulk_batch_all.csv"))
bm[, czi := basename(gsub("\\", "/", czi_path, fixed = TRUE))]
bm[, batch := as.integer(trimws(batch))]
bm <- unique(bm[, .(strain = trimws(mutant_or_drug), czi, batch)])
lib <- rbindlist(lapply(c("all_features.csv", "all_features 2.csv"), function(f)
  fread(file.path(lib_dir, f), select = c("well", "source_czi", unname(FEATS)))))
lib[, strain := tstrsplit(well, "__")[[3]]]
lib <- merge(lib, bm, by.x = c("strain", "source_czi"), by.y = c("strain", "czi"))
setnames(lib, unname(FEATS), names(FEATS))

sb <- lib[, lapply(.SD, median), by = .(strain, batch), .SDcols = names(FEATS)]
nt <- sb[strain %in% NT_STRAINS, lapply(.SD, mean), by = batch, .SDcols = names(FEATS)]
sb <- merge(sb, nt, by = "batch", suffixes = c("", "_ref"))
for (f in names(FEATS)) sb[, (f) := pct(get(f), get(paste0(f, "_ref")))]
kd <- sb[, lapply(.SD, mean), by = strain, .SDcols = names(FEATS)]
kd[, group := fifelse(strain %in% names(kd_group), kd_group[strain], "Other knockdown")]
kd[, group := factor(group, levels = names(group_cols))]
fwrite(kd, file.path(out_dir, "knockdown_shape_change_vs_NT.csv"))

# ---- drugs: per well medians -> condition mean vs RIF 0.25x -----------------
dr <- rbindlist(lapply(days, function(dy) {
  x <- fread(file.path(data_dir, sprintf("all_features_%s.csv", dy)),
             select = c("well", unname(FEATS)))
  x[, day := dy]
}))
setnames(dr, unname(FEATS), names(FEATS))
dr[, c("cond", "rep", "dil") := tstrsplit(well, "__")]
dw <- dr[, lapply(.SD, median), by = .(day, cond, dil), .SDcols = names(FEATS)]
dc <- dw[, lapply(.SD, mean), by = .(day, cond), .SDcols = names(FEATS)]
dc <- merge(dc, dc[cond == DRUG_REF, -"cond"], by = "day", suffixes = c("", "_ref"))
for (f in names(FEATS)) dc[, (f) := pct(get(f), get(paste0(f, "_ref")))]
dc <- dc[, c("day", "cond", names(FEATS)), with = FALSE]
dc[, c("drug", "conc") := tstrsplit(cond, "_")]
dc[, conc := factor(conc, levels = conc_levels)]
setorder(dc, day, drug, conc)
fwrite(dc, file.path(out_dir, "drug_shape_change_vs_RIF025.csv"))

message("% change vs reference (drugs vs ", DRUG_REF, "; knockdowns vs own-batch NT):")
print(dc[, .(day, cond, width = round(width, 1), length = round(length, 1),
             area = round(area, 1), reporter = round(reporter))])
print(kd[group != "Other knockdown", .(strain, group, width = round(width, 1),
                                       length = round(length, 1), area = round(area, 1),
                                       reporter = round(reporter))][order(group, -width)])

# ---- figures ----------------------------------------------------------------
drug_cols <- c(EMB = "#386cb0", INH = "#fdb462", MOX = "#7fc97f", RIF = "#ef3b2c")
dose_size <- c("0.25xMIC" = 2.2, "0.5xMIC" = 3.2, "1xMIC" = 4.4)
drugs_plot <- dc[cond != DRUG_REF]

# Each drug's dose path starts at the shared reference (0, 0)
path <- rbindlist(lapply(c("EMB", "INH", "MOX"), function(dg)
  rbind(data.table(day = days, drug = dg, conc = factor(NA, levels = conc_levels),
                   width = 0, length = 0),
        dc[drug == dg, .(day, drug, conc, width, length)])))
path <- rbind(path, dc[drug == "RIF", .(day, drug, conc, width, length)])
setorder(path, day, drug, conc, na.last = FALSE)

width_length_panel <- function(dy, legend) {
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +
    geom_point(data = kd, aes(width, length, colour = group), size = 1.8, alpha = 0.85) +
    geom_text_repel(data = kd[strain %in% LABEL_KD],
                    aes(width, length, label = strain, colour = group),
                    size = 2.6, fontface = "italic", min.segment.length = 0, seed = 1,
                    show.legend = FALSE, max.overlaps = 30) +
    geom_path(data = path[day == dy], aes(width, length, group = drug), colour = "grey30",
              linewidth = 0.5, arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
    geom_point(data = drugs_plot[day == dy], aes(width, length, fill = drug, size = conc),
               shape = 23, colour = "black", stroke = 0.6) +
    scale_colour_manual(values = group_cols, name = "Knockdown") +
    scale_fill_manual(values = drug_cols, name = "Drug") +
    scale_size_manual(values = dose_size, name = "Dose") +
    labs(x = "Width change (%)", y = "Length change (%)", title = dy,
         subtitle = sprintf("Drugs (diamonds, arrows = increasing dose) vs %s; knockdowns (dots) vs own-batch NT",
                            sub("_", " ", DRUG_REF))) +
    theme_Publication(base_size = 11) +
    theme(legend.position = if (legend) "right" else "none", legend.direction = "vertical",
          legend.title = element_text(size = 9), plot.subtitle = element_text(size = 8))
}
p1 <- (width_length_panel(days[1], FALSE) | width_length_panel(days[length(days)], TRUE)) +
  plot_annotation(title = "Which shape axes match: width vs length",
                  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))
save_fig(p1, "01_width_vs_length", 13, 6.2)

# Reporter vs width: does iniB induction track widening in both systems?
kd_rep <- rbindlist(lapply(days, function(dy) copy(kd)[, day := dy]))
kd_rep[, reporter_plot := pmax(reporter, -50)]
p2 <- ggplot() +
  geom_point(data = kd_rep, aes(width, reporter_plot, colour = group), size = 1.8, alpha = 0.85) +
  geom_text_repel(data = kd_rep[strain %in% LABEL_KD],
                  aes(width, reporter_plot, label = strain, colour = group),
                  size = 2.6, fontface = "italic", seed = 1, show.legend = FALSE) +
  geom_point(data = drugs_plot, aes(width, reporter, fill = drug, size = conc),
             shape = 23, colour = "black", stroke = 0.6) +
  scale_colour_manual(values = group_cols, name = "Knockdown") +
  scale_fill_manual(values = drug_cols, name = "Drug") +
  scale_size_manual(values = dose_size, name = "Dose") +
  facet_wrap(~ day, nrow = 1) +
  labs(x = "Width change (%)", y = "iniB reporter change (%)",
       title = "Reporter induction vs widening",
       subtitle = "Knockdowns repeated in each panel; reporter changes between imaging runs are descriptive") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.title = element_text(size = 9), plot.subtitle = element_text(size = 8))
save_fig(p2, "02_reporter_vs_width", 12, 5.6)

message("\nDone. Tables in ", out_dir, " | figures in ", fig_dir)
