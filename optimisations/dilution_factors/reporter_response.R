# =============================================================================
# iniB-mScarlet reporter response to drugs (Plate60, day 2 vs day 4)
#
# Reporter readout = per-cell intensity_median_mCher (never intensity_total,
# which tracks cell area).  Uninduced cells sit on a background floor (~63-69),
# so the reporter behaves as on/off at the single-cell level.
#
# No untreated wells on these plates, so "induced" is defined against an
# internal uninduced reference: MOX- and RIF-treated cells, which stay at
# background (iniBAC responds to cell-wall stress - INH, EMB - not to RIF).
# A cell is induced if it is brighter than the 99th percentile of that day's
# MOX + RIF cells (1% false-positive rate by construction).  Absolute
# intensities are compared between days only descriptively (separate plates).
#
# Each drug x dose has 5 wells (the dilution series).  Per-well FOV spread is
# shown so well-to-well differences can be told apart from imaging noise.
#
# Run from the project root:
#   Rscript optimisations/dilution_factors/reporter_response.R [data_dir] [day2,day4]
# Outputs -> <data_dir>/reporter/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggridges)
})

if (file.exists("Theme.R")) source("Theme.R") else
  stop("Theme.R not found - run from project root")

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "optimisations/dilution_factors"
days     <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("day2", "day4")

out_dir <- file.path(data_dir, "reporter")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(p, name, w, h, dpi = 300) {
  # cairo is unavailable without XQuartz; quartz (pdf) + ragg (png) instead
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white",
         device = function(filename, ...) grDevices::quartz(type = "pdf", file = filename, ...))
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = dpi, bg = "white", device = ragg::agg_png)
}

REPORTER     <- "intensity_median_mCher"
REF_DRUGS    <- c("MOX", "RIF")   # uninduced reference
REF_QUANTILE <- 0.99
MIN_CELLS    <- 20L              # per class, for within-well induced-vs-not contrasts

conc_levels <- c("0.25xMIC", "0.5xMIC", "1xMIC")
dil_levels  <- c(20L, 50L, 100L, 150L, 200L)
day_cols    <- setNames(c("#386cb0", "#ef3b2c", "#7fc97f", "#fdb462")[seq_along(days)], days)

# ---- load -------------------------------------------------------------------
d <- rbindlist(lapply(days, function(dy) {
  x <- fread(file.path(data_dir, sprintf("all_features_%s.csv", dy)),
             select = c("well", "fov", "length_um", "width_median_um", "area_um2", REPORTER))
  x[, day := dy]
}))
setnames(d, REPORTER, "reporter")
d[, c("dc", "rep", "dil") := tstrsplit(well, "__")]
d[, c("drug", "conc") := tstrsplit(dc, "_")]
d[, dil := as.integer(gsub("1in|_focused", "", dil))]
d[, conc := factor(conc, levels = conc_levels)]
d[, c("dc", "rep") := NULL]

# ---- induced threshold per day ----------------------------------------------
thr <- d[drug %in% REF_DRUGS, .(floor = median(reporter),
                                threshold = quantile(reporter, REF_QUANTILE)), by = day]
d <- merge(d, thr, by = "day")
d[, induced := reporter > threshold]
message("Uninduced reference (", paste(REF_DRUGS, collapse = " + "), "), per day:")
print(thr)

# ---- per-well and per-FOV summaries -----------------------------------------
fov <- d[, .(frac_induced = mean(induced), n = .N), by = .(day, drug, conc, dil, fov)]
wl <- d[, .(n_cells = .N, frac_induced = mean(induced),
            median_all = median(reporter),
            median_induced = if (sum(induced) >= 5) median(reporter[induced]) else NA_real_),
        by = .(day, drug, conc, dil)]
wl <- merge(wl, fov[, .(fov_q25 = quantile(frac_induced, 0.25),
                        fov_q75 = quantile(frac_induced, 0.75)),
                    by = .(day, drug, conc, dil)], by = c("day", "drug", "conc", "dil"))
setorder(wl, day, drug, conc, dil)
fwrite(wl, file.path(out_dir, "reporter_per_well.csv"))

cond <- wl[, .(frac_induced_mean = mean(frac_induced),
               frac_induced_min = min(frac_induced), frac_induced_max = max(frac_induced),
               median_induced = median(median_induced, na.rm = TRUE)),
           by = .(day, drug, conc)]
fwrite(cond, file.path(out_dir, "reporter_per_condition.csv"))
message("\nFraction of cells induced (mean of 5 wells, [min-max]) and induced-cell level:")
print(cond[, .(day, drug, conc,
               induced = sprintf("%.0f%% [%.0f-%.0f]", 100 * frac_induced_mean,
                                 100 * frac_induced_min, 100 * frac_induced_max),
               level = round(median_induced))])

# How much of the variation is between wells rather than between FOVs?
# Share of FOV-level variance in induced fraction explained by well (EMB/INH).
icc <- fov[drug %in% c("EMB", "INH"), {
  m <- lm(frac_induced ~ factor(dil))
  .(well_r2 = summary(m)$r.squared)
}, by = .(day, drug, conc)]
message("\nShare of FOV-level variance in induced fraction explained by well (R^2):")
print(icc, digits = 2)
fwrite(icc, file.path(out_dir, "well_vs_fov_variance.csv"))

# ---- reporter vs morphology, within well ------------------------------------
# Compare induced vs uninduced cells inside the same well (same drug, dose,
# culture and image set), so the contrast is not confounded by condition.
within <- d[drug %in% c("EMB", "INH"), {
  ni <- sum(induced); nu <- sum(!induced)
  if (ni >= MIN_CELLS && nu >= MIN_CELLS)
    .(n_induced = ni, n_uninduced = nu,
      d_width  = median(width_median_um[induced]) - median(width_median_um[!induced]),
      d_length = median(length_um[induced]) - median(length_um[!induced]),
      rho_width_in_induced = cor(reporter[induced], width_median_um[induced], method = "spearman"),
      rho_length_in_induced = cor(reporter[induced], length_um[induced], method = "spearman"))
}, by = .(day, drug, conc, dil)]
fwrite(within, file.path(out_dir, "reporter_vs_morphology_within_well.csv"))
message("\nInduced minus uninduced cells, same well (median um), wells with >= ", MIN_CELLS,
        " cells of each class:")
print(within[, .(wells = .N, d_width = round(median(d_width), 3),
                 width_up = sprintf("%d/%d", sum(d_width > 0), .N),
                 d_length = round(median(d_length), 2),
                 length_up = sprintf("%d/%d", sum(d_length > 0), .N)),
             by = .(day, drug, conc)])

# ---- figures ----------------------------------------------------------------
d[, cond_lab := factor(paste(drug, conc),
                       levels = rev(unlist(lapply(c("EMB", "INH", "MOX", "RIF"),
                                                  function(x) paste(x, conc_levels)))))]

# 1. Per-cell distributions (all dilution wells pooled), log scale
p1 <- ggplot(d, aes(reporter, cond_lab, fill = drug)) +
  geom_density_ridges(scale = 1.4, rel_min_height = 0.005, alpha = 0.85,
                      colour = "grey25", linewidth = 0.25) +
  geom_vline(data = thr, aes(xintercept = threshold), linetype = "dashed", colour = "grey30") +
  scale_x_log10() +
  scale_fill_Publication() +
  facet_wrap(~ day, nrow = 1) +
  labs(x = "iniB-mScarlet, per-cell median intensity (log scale)", y = NULL,
       title = "iniB reporter per cell",
       subtitle = sprintf("All dilution wells pooled. Dashed = induced threshold (%.0fth percentile of %s cells)",
                          100 * REF_QUANTILE, paste(REF_DRUGS, collapse = " + "))) +
  theme_Publication(base_size = 11) +
  theme(legend.position = "none", plot.subtitle = element_text(size = 8))
save_fig(p1, "01_reporter_distributions", 10, 7)

# 2. Fraction induced: every replicate well (with FOV IQR), dose response by day
wl[, dil_f := factor(paste0("1:", dil), levels = paste0("1:", dil_levels))]
p2 <- ggplot(wl, aes(conc, frac_induced, colour = day)) +
  geom_linerange(aes(ymin = fov_q25, ymax = fov_q75, group = interaction(day, dil)),
                 position = position_dodge(width = 0.6), alpha = 0.5) +
  geom_point(aes(group = interaction(day, dil)), position = position_dodge(width = 0.6),
             size = 1.6, alpha = 0.8) +
  stat_summary(aes(group = day), fun = mean, geom = "line", linewidth = 1) +
  scale_colour_manual(values = day_cols) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  facet_wrap(~ drug, nrow = 1) +
  labs(x = "Dose", y = "Cells induced",
       title = "Fraction of cells with iniB on",
       subtitle = "Points = the 5 replicate (dilution) wells; bars = IQR across that well's FOVs; line = mean of wells") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))

# 3. Level among induced cells
p3 <- ggplot(wl[!is.na(median_induced) & drug %in% c("EMB", "INH")],
             aes(conc, median_induced, colour = day)) +
  geom_point(position = position_dodge(width = 0.5), size = 1.6, alpha = 0.8) +
  stat_summary(aes(group = day), fun = median, geom = "line", linewidth = 1,
               position = position_dodge(width = 0.5)) +
  geom_hline(data = thr, aes(yintercept = floor, colour = day), linetype = "dotted") +
  scale_colour_manual(values = day_cols) +
  scale_y_log10() +
  facet_wrap(~ drug, nrow = 1) +
  labs(x = "Dose", y = "Median intensity of induced cells",
       title = "How bright induced cells get",
       subtitle = "Points = wells; dotted = background floor (median MOX + RIF cell)") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))
p23 <- (p2 / p3) + plot_layout(heights = c(1, 0.9))
save_fig(p23, "02_reporter_dose_response", 10, 8.5)

# 4. Replicate wells: induced fraction by dilution well (reveals well differences)
p4 <- ggplot(wl[drug %in% c("EMB", "INH")], aes(dil_f, frac_induced, colour = conc, group = conc)) +
  geom_linerange(aes(ymin = fov_q25, ymax = fov_q75), alpha = 0.6) +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  scale_colour_manual(values = c("0.25xMIC" = "#386cb0", "0.5xMIC" = "#fdb462", "1xMIC" = "#ef3b2c")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  facet_grid(day ~ drug) +
  labs(x = "Dilution well", y = "Cells induced",
       title = "Replicate-well consistency",
       subtitle = "Bars = IQR across the 16 FOVs of each well; tight bars + different wells = a real well-level difference") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p4, "03_replicate_well_consistency", 9, 7)

# 5. Reporter vs morphology within well
if (nrow(within)) {
  wm <- melt(within, id.vars = c("day", "drug", "conc", "dil"),
             measure.vars = c("d_width", "d_length"), variable.name = "feature")
  wm[, feature := factor(feature, levels = c("d_width", "d_length"),
                         labels = c("Width (µm)", "Length (µm)"))]
  p5 <- ggplot(wm, aes(conc, value, colour = day)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_point(position = position_dodge(width = 0.5), size = 1.8, alpha = 0.8) +
    stat_summary(aes(group = day), fun = median, geom = "crossbar", width = 0.35,
                 position = position_dodge(width = 0.5), linewidth = 0.4) +
    scale_colour_manual(values = day_cols) +
    facet_grid(feature ~ drug, scales = "free_y") +
    labs(x = "Dose", y = "Induced minus uninduced (median)",
         title = "Are iniB-on cells shaped differently?",
         subtitle = sprintf("Within-well contrast (same drug, dose, well); points = wells with >= %d cells in each class",
                            MIN_CELLS)) +
    theme_Publication(base_size = 11) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          plot.subtitle = element_text(size = 8))
  save_fig(p5, "04_reporter_vs_morphology", 8.5, 6.5)
}

message("\nDone. Tables in ", out_dir, " | figures in ", fig_dir)
