# =============================================================================
# Fixed-cell dilution optimisation (Plate60, day 4)
#
# Question: which dilution of the PBS-resuspended fixed pellet gives the most
# *usable* cells per FOV for downstream morphological profiling?
#
# Design: 4 drugs (EMB, INH, MOX, RIF) x 3 doses (0.25/0.5/1x MIC) x
#         5 dilutions (1:20, 1:50, 1:100, 1:150, 1:200), 16 FOVs per well,
#         IniB-mScarlet reporter strain. One well per condition x dilution.
#
# "Usable" = segmented AND not in contact with another cell. Touching cells are
# the main risk to single-cell morphology (merged / mis-split objects), so the
# yield metric is isolated cells per FOV, not raw objects per FOV.
#
# Contact is estimated by modelling each cell as a capsule (line segment along
# the major axis, radius = half the minor axis) and computing the edge-to-edge
# gap to its neighbours. skimage orientation -> major-axis direction (x, y) =
# (sin theta, cos theta); checked against bbox extents (r = 0.994).
#
# Run from the project root:
#   Rscript optimisations/dilution_factors/dilution_optimisation.R [features.csv] [od.RUC]
# Outputs (summary CSVs + figures/) are written next to features.csv.
# If no .RUC is given, the single *.RUC next to features.csv is used (OD600
# section is skipped when there is none, or more than one).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(RANN)
  library(lme4)
})

if (file.exists("Theme.R")) source("Theme.R") else
  stop("Theme.R not found - run from project root")

args    <- commandArgs(trailingOnly = TRUE)
in_file <- if (length(args)) args[1] else "optimisations/dilution_factors/all_features_day4.csv"

# All outputs (CSVs + figures/) go next to the input data
out_dir <- dirname(normalizePath(in_file, mustWork = TRUE))
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ruc_file <- if (length(args) >= 2) args[2] else
  list.files(out_dir, pattern = "\\.RUC$", full.names = TRUE)
if (length(ruc_file) > 1) {
  message("Several .RUC files found - pass the one to use as the 2nd argument; skipping OD600 section:\n  ",
          paste(basename(ruc_file), collapse = "\n  "))
  ruc_file <- character(0)
}

save_fig <- function(p, name, w, h, dpi = 300) {
  # cairo is unavailable without XQuartz; quartz (pdf) + ragg (png) render µ correctly
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white",
         device = function(filename, ...) grDevices::quartz(type = "pdf", file = filename, ...))
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = dpi, bg = "white", device = ragg::agg_png)
}

# ---- parameters -------------------------------------------------------------
N_FOV        <- 16L     # FOVs acquired per well (fov 0-15)
GAP_TOUCH_UM <- 0.3     # edge-to-edge gap below which two cells count as touching
TARGET_CELLS <- 500L    # isolated cells wanted per condition downstream
K_NN         <- 12L     # neighbours checked per cell for contact

dil_levels  <- c(20L, 50L, 100L, 150L, 200L)
dil_lab     <- function(x) factor(paste0("1:", x), levels = paste0("1:", dil_levels))
conc_levels <- c("0.25xMIC", "0.5xMIC", "1xMIC")

# ---- load -------------------------------------------------------------------
d <- fread(in_file)
d[, c("dc", "reporter", "dil") := tstrsplit(well, "__")]
d[, c("drug", "conc") := tstrsplit(dc, "_")]
d[, dil := as.integer(sub("1in", "", sub("_focused", "", dil)))]
d[, conc := factor(conc, levels = conc_levels)]
d[, condition := paste(drug, conc)]
d[, img_well := sub("^Plate\\d+_([A-H]\\d+)_.*$", "\\1", source_czi)]
d[, c("dc", "reporter") := NULL]
px_um <- sqrt(d$area_um2[1] / d$area_px[1])   # 0.0721 um / px
fov_area_um2 <- (2048 * px_um)^2               # ~21,800 um^2 per FOV

# ---- edge-to-edge gap between rods (capsule model) --------------------------
# Minimum distance between segments P1-Q1 and P2-Q2 (vectorised over pairs);
# standard clamped closest-points solution.
seg_dist <- function(p1x, p1y, q1x, q1y, p2x, p2y, q2x, q2y) {
  d1x <- q1x - p1x; d1y <- q1y - p1y
  d2x <- q2x - p2x; d2y <- q2y - p2y
  rx  <- p1x - p2x; ry  <- p1y - p2y
  a <- d1x^2 + d1y^2; e <- d2x^2 + d2y^2
  f <- d2x * rx + d2y * ry; c <- d1x * rx + d1y * ry
  b <- d1x * d2x + d1y * d2y
  den <- a * e - b^2
  s <- ifelse(den > 1e-12, pmin(pmax((b * f - c * e) / den, 0), 1), 0)
  t <- ifelse(e > 1e-12, (b * s + f) / e, 0)
  # clamp t, then recompute s
  s <- ifelse(t < 0, ifelse(a > 1e-12, pmin(pmax(-c / a, 0), 1), 0), s)
  s <- ifelse(t > 1, ifelse(a > 1e-12, pmin(pmax((b - c) / a, 0), 1), 0), s)
  t <- pmin(pmax(t, 0), 1)
  cx <- (p1x + d1x * s) - (p2x + d2x * t)
  cy <- (p1y + d1y * s) - (p2y + d2y * t)
  sqrt(cx^2 + cy^2)
}

# Capsule geometry in um
d[, `:=`(cx = centroid_x * px_um, cy = centroid_y * px_um,
         rad = minor_axis_length_um / 2)]
d[, half := pmax(0, (major_axis_length_um - minor_axis_length_um) / 2)]
d[, `:=`(ux = sin(orientation_rad), uy = cos(orientation_rad))]
d[, `:=`(px_ = cx - ux * half, py_ = cy - uy * half,
         qx_ = cx + ux * half, qy_ = cy + uy * half)]

# Nearest edge-to-edge gap + connected clusters of touching cells, per FOV
fov_contacts <- function(D) {
  n <- nrow(D)
  if (n < 2) return(list(gap = rep(Inf, n), clus = rep(1L, n)))
  k  <- min(K_NN + 1L, n)
  nn <- nn2(cbind(D$cx, D$cy), k = k)$nn.idx[, -1, drop = FALSE]
  i  <- rep(seq_len(n), ncol(nn)); j <- as.vector(nn)
  g  <- seg_dist(D$px_[i], D$py_[i], D$qx_[i], D$qy_[i],
                 D$px_[j], D$py_[j], D$qx_[j], D$qy_[j]) - D$rad[i] - D$rad[j]
  g  <- pmax(g, 0)
  gap <- tapply(g, i, min)
  # union-find over touching pairs
  par <- seq_len(n)
  find <- function(x) { while (par[x] != x) x <- par[x]; x }
  for (p in which(g < GAP_TOUCH_UM)) {
    a <- find(i[p]); b <- find(j[p]); if (a != b) par[b] <- a
  }
  root <- vapply(seq_len(n), find, 1L)
  list(gap = as.numeric(gap), clus = as.integer(ave(root, root, FUN = length)))
}

d[, c("nn_gap_um", "cluster_size") := fov_contacts(.SD), by = .(well, fov),
  .SDcols = c("cx", "cy", "rad", "px_", "py_", "qx_", "qy_")]
d[, isolated := nn_gap_um >= GAP_TOUCH_UM]

# Objects that look like segmentation failures (merged / clumped masks)
d[, suspect := branch_count > 0 | solidity < 0.75 | area_um2 > 5]

# ---- FOV-level table (fill empty FOVs with 0) -------------------------------
grid <- unique(d[, .(drug, conc, condition, dil, well)])[
  , .(fov = 0:(N_FOV - 1)), by = .(drug, conc, condition, dil, well)]
fv <- d[, .(cells = .N, iso = sum(isolated), suspect = sum(suspect),
            area_cov = sum(area_um2)), by = .(well, fov)]
fv <- merge(grid, fv, by = c("well", "fov"), all.x = TRUE)
for (v in c("cells", "iso", "suspect", "area_cov")) set(fv, which(is.na(fv[[v]])), v, 0)
fv[, `:=`(frac_touch = ifelse(cells > 0, 1 - iso / cells, NA_real_),
          coverage = 100 * area_cov / fov_area_um2,
          dil_f = dil_lab(dil))]

# ---- well-level summary ----------------------------------------------------
wl <- fv[, .(cells_fov = mean(cells), iso_fov = mean(iso),
             cv_fov = sd(cells) / mean(cells),
             coverage = mean(coverage)), by = .(drug, conc, condition, dil)]
wl2 <- d[, .(img_well = img_well[1], frac_touch = 1 - mean(isolated),
             frac_suspect = mean(suspect), frac_in_clus3 = mean(cluster_size >= 3)),
         by = .(drug, conc, condition, dil)]
wl <- merge(wl, wl2, by = c("drug", "conc", "condition", "dil"))
wl[, dil_f := dil_lab(dil)]
wl[, fov_needed := ceiling(TARGET_CELLS / iso_fov)]
setorder(wl, drug, conc, dil)

# ---- 1. Does cell number scale with dilution? -------------------------------
# Expected: halving concentration halves cells (log-log slope = -1).
fit_slope <- lmer(log(cells + 1) ~ log(dil) + (1 + log(dil) | condition), data = fv)
slope <- fixef(fit_slope)[["log(dil)"]]
slope_ci <- confint(fit_slope, parm = "log(dil)", method = "Wald")
message(sprintf("Pooled log-log slope (cells/FOV ~ dilution): %.2f [%.2f, %.2f]; linear expectation = -1",
                slope, slope_ci[1], slope_ci[2]))

conc_cols <- c("0.25xMIC" = "#386cb0", "0.5xMIC" = "#fdb462", "1xMIC" = "#ef3b2c")
dil_cols  <- setNames(c("#662506", "#ef3b2c", "#fdb462", "#7fc97f", "#386cb0"),
                      paste0("1:", dil_levels))

# Expected-if-proportional guide: anchored at each condition's 1:200 mean
exp_line <- wl[dil == 200, .(condition, drug, conc, base = cells_fov)][
  , .(dil = dil_levels, exp = base * 200 / dil_levels), by = .(condition, drug, conc)]

p1 <- ggplot(fv, aes(dil, cells + 1, colour = conc)) +
  geom_line(data = exp_line, aes(dil, exp, group = condition, colour = conc),
            linetype = "dashed", linewidth = 0.4, alpha = 0.6) +
  geom_point(position = position_jitter(width = 0.03, height = 0, seed = 1),
             size = 0.8, alpha = 0.35) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.9) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  scale_x_log10(breaks = dil_levels, labels = paste0("1:", dil_levels)) +
  scale_y_log10() +
  scale_colour_manual(values = conc_cols) +
  facet_wrap(~ drug, nrow = 1) +
  labs(x = "Dilution", y = "Cells per FOV",
       title = "Cell number vs dilution",
       subtitle = sprintf("Points = FOVs; solid = mean; dashed = proportional scaling from 1:200. Pooled log-log slope = %.2f (expected -1)",
                          slope)) +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p1, "01_cells_per_fov_vs_dilution", 11, 4.5)

# ---- 2. Crowding vs density (FOV level) -------------------------------------
fit_touch <- glmer(cbind(cells - iso, iso) ~ log(cells) + (1 | condition),
                   data = fv[cells >= 5], family = binomial)
nd <- data.frame(cells = exp(seq(log(5), log(max(fv$cells)), length.out = 200)))
nd$frac_touch <- plogis(predict(fit_touch, newdata = nd, re.form = NA))
nd$iso <- nd$cells * (1 - nd$frac_touch)
print(summary(fit_touch)$coefficients)
message(sprintf("Predicted %% touching at 25 / 100 / 250 cells per FOV: %s",
                paste(round(100 * plogis(predict(fit_touch, re.form = NA,
                  newdata = data.frame(cells = c(25, 100, 250))))), collapse = " / ")))
message(sprintf("Max area coverage of any FOV: %.1f%%", max(fv$coverage)))

p2a <- ggplot(fv[cells >= 5], aes(cells, frac_touch)) +
  geom_point(aes(colour = dil_f), size = 1, alpha = 0.6) +
  geom_line(data = nd, linewidth = 1) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_manual(values = dil_cols) +
  labs(x = "Cells per FOV", y = "Cells touching a neighbour",
       title = "Crowding",
       subtitle = sprintf("Edge-to-edge gap < %.1f \u00b5m; line = binomial GLMM", GAP_TOUCH_UM)) +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.subtitle = element_text(size = 8))

p2b <- ggplot(fv, aes(cells, iso)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
  geom_point(aes(colour = dil_f), size = 1, alpha = 0.6) +
  geom_line(data = nd, linewidth = 1) +
  scale_colour_manual(values = dil_cols) +
  labs(x = "Cells per FOV", y = "Isolated cells per FOV",
       title = "Usable yield",
       subtitle = "Dotted = every cell isolated") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.subtitle = element_text(size = 8))

p2 <- (p2a | p2b) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")
save_fig(p2, "02_crowding_vs_density", 9.5, 4.8)

# ---- 3. Does dilution bias morphology? --------------------------------------
# Within each condition, centre FOV-median features on the condition mean, then
# ask whether they trend with dilution (mixed model, condition random effect).
morph_feats <- c(length_um = "Length (\u00b5m)", width_median_um = "Width (\u00b5m)",
                 area_um2 = "Area (\u00b5m\u00b2)", solidity = "Solidity")
fm <- d[, lapply(.SD, median), by = .(condition, drug, conc, dil, well, fov),
        .SDcols = names(morph_feats)]
fm <- melt(fm, measure.vars = names(morph_feats), variable.name = "feature")
fm[, value_c := value - mean(value), by = .(condition, feature)]
fm[, value_z := value_c / sd(value), by = feature]
fm[, dil_f := dil_lab(dil)]

morph_tests <- fm[, {
  m <- lmer(value_z ~ log(dil) + (1 | condition) + (1 | well), data = .SD)
  co <- summary(m)$coefficients
  .(slope_z_per_log = co["log(dil)", "Estimate"], se = co["log(dil)", "Std. Error"],
    t = co["log(dil)", "t value"])
}, by = feature]
print(morph_tests)

fm_w <- fm[, .(value_c = median(value_c)), by = .(condition, conc, drug, feature, dil, dil_f)]
fm_w[, feature_lab := factor(morph_feats[as.character(feature)], levels = morph_feats)]
p3 <- ggplot(fm_w, aes(dil_f, value_c)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(aes(group = condition, colour = drug), alpha = 0.6, linewidth = 0.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 1.1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  scale_colour_Publication() +
  facet_wrap(~ feature_lab, scales = "free_y", nrow = 1) +
  labs(x = "Dilution", y = "Shift from condition mean",
       title = "Morphology vs dilution",
       subtitle = "Coloured = each drug x dose (centred on its own mean); black = average across conditions") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal", plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p3, "03_morphology_vs_dilution", 11, 4.2)

# ---- 4. Decision summary ----------------------------------------------------
wl[, iso_rel := iso_fov / max(iso_fov), by = condition]
wl[, best := iso_fov == max(iso_fov), by = condition]

sm <- wl[, .(cells_fov = median(cells_fov), iso_fov = median(iso_fov),
             iso_fov_min = min(iso_fov), iso_fov_max = max(iso_fov),
             iso_rel = mean(iso_rel), n_best = sum(best),
             frac_touch = median(frac_touch), frac_suspect = median(frac_suspect),
             cv_fov = median(cv_fov), fov_needed = median(fov_needed),
             fov_needed_max = max(fov_needed)),
         by = .(dil)][order(dil)]
fwrite(sm, file.path(out_dir, "dilution_summary_overall.csv"))
print(sm)

summary_panel <- function(y, ylab, title, subtitle, pct = FALSE, log_y = FALSE) {
  p <- ggplot(wl, aes(dil_f, .data[[y]])) +
    geom_line(aes(group = condition, colour = drug), alpha = 0.5) +
    geom_point(aes(colour = drug), size = 1.3, alpha = 0.7) +
    stat_summary(aes(group = 1), fun = median, geom = "line", linewidth = 1.1) +
    stat_summary(fun = median, geom = "point", size = 2.4) +
    scale_colour_Publication() +
    labs(x = "Dilution", y = ylab, title = title, subtitle = subtitle) +
    theme_Publication(base_size = 11) +
    theme(legend.position = "bottom", legend.direction = "horizontal", plot.subtitle = element_text(size = 8))
  if (pct)   p <- p + scale_y_continuous(labels = scales::percent)
  if (log_y) p <- p + scale_y_log10()
  p
}

p4a <- summary_panel("iso_fov", "Isolated cells per FOV", "Usable cells / FOV",
                     "Lines = drug x dose; black = median")
p4b <- summary_panel("fov_needed", "FOVs needed",
                     sprintf("FOVs for %d isolated cells", TARGET_CELLS),
                     sprintf("Dashed = %d FOVs acquired here", N_FOV), log_y = TRUE) +
  geom_hline(yintercept = N_FOV, linetype = "dashed", colour = "grey40")
p4c <- summary_panel("frac_suspect", "Suspect objects", "Segmentation quality",
                     "Branched, solidity < 0.75 or area > 5 \u00b5m\u00b2", pct = TRUE)
p4d <- summary_panel("cv_fov", "CV of cells per FOV", "FOV-to-FOV evenness",
                     "Lower = more predictable yield", pct = TRUE)

p4 <- (p4a | p4b) / (p4c | p4d) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")
save_fig(p4, "04_decision_summary", 10, 8.5)

# ---- 5. Per-condition heatmap: isolated cells / FOV -------------------------
p5 <- ggplot(wl, aes(dil_f, conc, fill = iso_fov)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = round(iso_fov), fontface = ifelse(best, "bold", "plain")),
            size = 3.2) +
  scale_fill_gradient(low = "#f0f0f0", high = "#386cb0", name = "Isolated\ncells/FOV") +
  facet_wrap(~ drug, ncol = 1, strip.position = "left") +
  labs(x = "Dilution", y = NULL, title = "Isolated cells per FOV",
       subtitle = "Bold = best dilution for that drug x dose") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.title = element_text(size = 9), axis.line = element_blank(),
        strip.placement = "outside", plot.subtitle = element_text(size = 8))
save_fig(p5, "05_isolated_cells_heatmap", 6, 7.5)

# ---- 6. Does culture OD600 explain cells per FOV? ---------------------------
# CLARIOstar .RUC = zip holding an "Absolute Database" file. Each well-channel
# record contains: 0xE0, well label ("A01".."H12"), ..., wavelength (double,
# +18 bytes) and reading (double, +26 bytes). Channel 600 nm is the OD; the
# 900/993 nm channels are pathlength-correction reads and are ignored.
read_clariostar_od <- function(ruc, wavelength = 600) {
  tmp <- tempfile(); on.exit(unlink(tmp, recursive = TRUE))
  utils::unzip(ruc, exdir = tmp)
  f <- list.files(tmp, pattern = "\\.abs$", full.names = TRUE)
  b <- readBin(f, "raw", file.info(f)$size)
  is_row <- b >= as.raw(0x41) & b <= as.raw(0x48)          # "A".."H"
  is_dig <- b >= as.raw(0x30) & b <= as.raw(0x39)
  n <- length(b)
  hit <- which(b[1:(n - 34)] == as.raw(0xE0))
  hit <- hit[is_row[hit + 1] & is_dig[hit + 2] & is_dig[hit + 3] & b[hit + 4] == as.raw(0)]
  rd  <- function(at) vapply(at, function(i) readBin(b[i:(i + 7)], "double", 1, 8, endian = "little"), 0)
  out <- data.table(well  = vapply(hit, function(i) rawToChar(b[(i + 1):(i + 3)]), ""),
                    wl    = rd(hit + 18),
                    value = rd(hit + 26))
  out <- unique(out[wl == wavelength, .(well, value)])
  stopifnot(nrow(out) == 96, !anyDuplicated(out$well))
  out[, `:=`(row = substr(well, 1, 1), col = as.integer(substring(well, 2)))]
  out[]
}

if (length(ruc_file) && file.exists(ruc_file)) {
  od <- read_clariostar_od(ruc_file)
  # Reader block A2-F11 = imaging wells A1-F10 with columns mirrored
  # (reader A2 -> imaging A10, reader A11 -> imaging A1).
  od[, img_well := paste0(row, 12L - col)]
  # Blank = edge wells flanking the imaging block (reader cols 1 and 12, rows A-F)
  od_blank <- median(od[col %in% c(1, 12) & row %in% LETTERS[1:6], value])
  od[, od_net := value - od_blank]
  message(sprintf("OD600 blank (reader cols 1 & 12, rows A-F): %.3f", od_blank))

  wl <- merge(wl, od[, .(img_well, reader_well = well, od600 = value, od_net)],
              by = "img_well", all.x = TRUE)
  setorder(wl, drug, conc, dil)

  # Does OD add anything once dilution is accounted for?  (well level, n = 60)
  fit_od  <- lm(log(cells_fov) ~ log(dil) + od_net, data = wl)
  fit_odc <- lm(log(cells_fov) ~ log(dil) + od_net + condition, data = wl)
  co  <- summary(fit_od)$coefficients["od_net", ]
  coc <- summary(fit_odc)$coefficients["od_net", ]
  message(sprintf("log(cells/FOV) ~ log(dil) + OD: OD coef %.2f, p = %.2f; within condition p = %.2f",
                  co[1], co[4], coc[4]))
  wl[, resid_dil := resid(lm(log(cells_fov) ~ log(dil), data = wl))]

  p6a <- ggplot(wl, aes(dil_f, conc, fill = od_net)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", od_net)), size = 3) +
    scale_fill_gradient(low = "#f0f0f0", high = "#7fc97f", name = "Net\nOD600") +
    facet_wrap(~ drug, ncol = 1, strip.position = "left") +
    labs(x = "Dilution", y = NULL, title = "Culture OD600",
         subtitle = sprintf("Blank-subtracted (blank = %.3f); same layout as isolated-cell heatmap",
                            od_blank)) +
    theme_Publication(base_size = 11) +
    theme(legend.position = "right", legend.direction = "vertical",
          legend.title = element_text(size = 9), axis.line = element_blank(),
          strip.placement = "outside", plot.subtitle = element_text(size = 8))

  p6b <- ggplot(wl, aes(od_net, exp(resid_dil))) +
    geom_hline(yintercept = 1, colour = "grey60") +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.8,
                se = TRUE, fill = "grey85") +
    geom_point(aes(colour = drug, shape = conc), size = 2.4, alpha = 0.85) +
    scale_y_log10() +
    scale_colour_Publication() +
    labs(x = "Net OD600", y = "Cells/FOV relative to dilution expectation",
         title = "OD vs cell yield",
         subtitle = sprintf("Residual after log(cells) ~ log(dilution). OD effect p = %.2f (p = %.2f within condition)",
                            co[4], coc[4])) +
    theme_Publication(base_size = 11) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.box = "vertical", plot.subtitle = element_text(size = 8))

  p6 <- (p6a | p6b) + plot_layout(widths = c(1, 1.15))
  save_fig(p6, "06_od600_vs_cell_yield", 12, 7.5)
} else {
  message("No .RUC plate-reader file found - skipping OD600 section")
}

fwrite(wl, file.path(out_dir, "dilution_summary_per_condition.csv"))
message("Done. Figures in ", fig_dir)
