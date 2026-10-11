# =============================================================================
# Culture OD600 over the drug time course (Plate60, CLARIOstar)
#
# Each imaging well has its own culture well on the reader plate (reader
# A2-F11 = imaging A1-F10, columns mirrored), so every drug x dose has 5
# replicate cultures.  Blank = reader row H, the only row that stays flat on
# every day; the other edge wells (column 1, column 12, row G) grow and are
# shown on the plate maps but not interpreted (contents not annotated).
#
# Outputs:
#   * growth curves per drug x dose (net OD vs day, 5 cultures each)
#   * whole-plate maps per day
#   * culture growth vs the imaged cells, per well: reporter induced fraction
#     and median width (needs reporter_response.R outputs).  A well whose
#     culture was inhibited but whose imaged cells look untreated flags a
#     problem between culture and imaging (fixation / dilution / plating).
#
# Run from the project root:
#   Rscript optimisations/dilution_factors/od_timecourse.R [data_dir]
# For a timepoint with several .RUC files, the one named "OD1" (else the
# first) is used.  Outputs -> <data_dir>/od_timecourse/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

if (file.exists("Theme.R")) source("Theme.R") else
  stop("Theme.R not found - run from project root")
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "plate_reader.R"))

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "optimisations/dilution_factors"

out_dir <- file.path(data_dir, "od_timecourse")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(p, name, w, h, dpi = 300) {
  # cairo is unavailable without XQuartz; quartz (pdf) + ragg (png) instead
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white",
         device = function(filename, ...) grDevices::quartz(type = "pdf", file = filename, ...))
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = dpi, bg = "white", device = ragg::agg_png)
}

conc_levels <- c("0.25xMIC", "0.5xMIC", "1xMIC")
dil_levels  <- c(20L, 50L, 100L, 150L, 200L)
conc_cols   <- c("0.25xMIC" = "#386cb0", "0.5xMIC" = "#fdb462", "1xMIC" = "#ef3b2c")

# ---- timepoints with features + a plate-reader file -------------------------
f <- list.files(data_dir, "^all_features_day[0-9]+\\.csv$")
days <- sub("^all_features_(day[0-9]+)\\.csv$", "\\1", f)
days <- days[order(as.integer(sub("day", "", days)))]
ruc <- vapply(days, function(dy) {
  r <- find_ruc(data_dir, dy)
  if (!length(r)) return(NA_character_)
  if (length(r) > 1) r <- c(r[grepl("OD1", basename(r))], r)[1]
  r
}, "")
message("Plate-reader files:\n", paste0("  ", days, ": ", basename(ruc), collapse = "\n"))
days <- days[!is.na(ruc)]; ruc <- ruc[!is.na(ruc)]
day_num <- as.integer(sub("day", "", days)); names(day_num) <- days

od <- rbindlist(lapply(days, function(dy) add_plate_layout(read_clariostar_od(ruc[[dy]]))[, day := dy]))
od[, day_n := day_num[day]]
message("\nBlank (row ", BLANK_ROW, ") per day: ",
        paste(sprintf("%s %.3f", days, od[, blank[1], by = day]$V1), collapse = ", "))

# imaging well -> condition (layout is the same on every day)
lay <- unique(fread(file.path(data_dir, sprintf("all_features_%s.csv", days[1])),
                    select = c("well", "source_czi")))
lay[, img_well := sub("^Plate\\d+_([A-H]\\d+)_.*$", "\\1", source_czi)]
lay[, c("cond", "rep", "dil") := tstrsplit(well, "__")]
lay[, c("drug", "conc") := tstrsplit(cond, "_")]
lay[, dil := as.integer(gsub("1in|_focused", "", dil))]
lay <- unique(lay[, .(img_well, drug, conc, dil)])

cw <- merge(od[!is.na(img_well)], lay, by = "img_well")
cw[, conc := factor(conc, levels = conc_levels)]
stopifnot(uniqueN(cw[, .(img_well, day)]) == 60 * length(days))
fwrite(cw[, .(day, reader_well = well, img_well, drug, conc, dil, od600 = value, blank, od_net)],
       file.path(out_dir, "od_per_culture_well.csv"))

cond <- cw[, .(od_net = mean(od_net), od_min = min(od_net), od_max = max(od_net)),
           by = .(day, day_n, drug, conc)]
fwrite(cond, file.path(out_dir, "od_per_condition.csv"))
message("\nNet OD600 (mean of 5 cultures) by day:")
print(dcast(cond, drug + conc ~ day, value.var = "od_net"), digits = 2)

# ---- 1. growth curves -------------------------------------------------------
p1 <- ggplot(cw, aes(day_n, od_net, colour = conc)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line(aes(group = img_well), alpha = 0.3) +
  stat_summary(aes(group = conc), fun = mean, geom = "line", linewidth = 1.1) +
  stat_summary(aes(group = conc), fun = mean, geom = "point", size = 2.4) +
  scale_colour_manual(values = conc_cols) +
  scale_x_continuous(breaks = day_num, labels = paste("Day", day_num)) +
  facet_wrap(~ drug, nrow = 1) +
  labs(x = NULL, y = sprintf("OD600 (minus row %s blank)", BLANK_ROW),
       title = "Culture growth under drug",
       subtitle = "Thin = each of the 5 cultures (one per imaging well); thick = mean") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))
save_fig(p1, "01_growth_curves", 11, 4.4)

# ---- 2. whole-plate maps ----------------------------------------------------
od[, `:=`(row_f = factor(row, levels = rev(LETTERS[1:8])),
          region = fifelse(!is.na(img_well), "imaging culture",
                           fifelse(row == BLANK_ROW, "blank", "other (not annotated)")))]
p2 <- ggplot(od, aes(col, row_f, fill = od_net)) +
  geom_tile(aes(colour = region), linewidth = 0.6, width = 0.94, height = 0.94) +
  geom_text(aes(label = sprintf("%.2f", od_net)), size = 2.1) +
  scale_fill_gradient(low = "#f7f7f7", high = "#386cb0", name = "Net OD600") +
  scale_colour_manual(values = c("imaging culture" = "grey35", "blank" = "#ef3b2c",
                                 "other (not annotated)" = "#fdb462"), name = NULL) +
  scale_x_continuous(breaks = 1:12, position = "top", expand = c(0, 0)) +
  facet_wrap(~ day, ncol = 1) +
  labs(x = NULL, y = NULL, title = "Reader plate, net OD600",
       subtitle = sprintf("Outline: grey = culture for an imaging well; red = blank (row %s); orange = other wells", BLANK_ROW)) +
  theme_Publication(base_size = 10) +
  theme(legend.position = "right", legend.direction = "vertical",
        axis.line = element_blank(), axis.ticks = element_blank(),
        plot.subtitle = element_text(size = 8))
save_fig(p2, "02_plate_maps", 8, 2.2 + 2.6 * length(days))

# ---- 3. culture growth vs the imaged cells ----------------------------------
rep_file <- file.path(data_dir, "reporter", "reporter_per_well.csv")
if (file.exists(rep_file)) {
  rp <- fread(rep_file)[, conc := factor(conc, levels = conc_levels)]
  # imaged-cell width per well
  wid <- rbindlist(lapply(days, function(dy) {
    x <- fread(file.path(data_dir, sprintf("all_features_%s.csv", dy)), select = c("well", "width_median_um"))
    x[, c("cond", "r", "dil") := tstrsplit(well, "__")]
    x[, c("drug", "conc") := tstrsplit(cond, "_")]
    x[, .(width = median(width_median_um)), by = .(drug, conc, dil = as.integer(gsub("1in|_focused", "", dil)))][, day := dy]
  }))
  link <- merge(cw, rp[, .(day, drug, conc, dil, frac_induced)], by = c("day", "drug", "conc", "dil"))
  link <- merge(link, wid, by = c("day", "drug", "conc", "dil"))

  # Flag: culture as inhibited as its siblings, but imaged cells far from them
  link[, `:=`(sib_induced = median(frac_induced), sib_width = median(width),
              sib_od = median(od_net)), by = .(day, drug, conc)]
  link[, flag := drug %in% c("EMB", "INH") & sib_induced > 0.5 &
         frac_induced < sib_induced - 0.4 & abs(od_net - sib_od) < 0.1]
  fwrite(link[, .(day, drug, conc, dil, img_well, reader_well = well, od_net, frac_induced, width,
                  sib_od, sib_induced, sib_width, flag)],
         file.path(out_dir, "culture_vs_imaging_per_well.csv"))
  message("\nWells whose culture grew like its siblings but whose imaged cells do not match them:")
  print(link[flag == TRUE, .(day, drug, conc, dil, img_well, od_net = round(od_net, 2),
                             sib_od = round(sib_od, 2), induced = round(frac_induced, 2),
                             sib_induced = round(sib_induced, 2), width = round(width, 3),
                             sib_width = round(sib_width, 3))])

  link[, dil_f := factor(paste0("1:", dil), levels = paste0("1:", dil_levels))]
  p3a <- ggplot(link[drug %in% c("EMB", "INH")], aes(od_net, frac_induced, colour = conc)) +
    geom_point(aes(shape = flag), size = 2.4, alpha = 0.85, stroke = 0.9) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                       labels = c(`FALSE` = "consistent", `TRUE` = "imaged cells don't match culture")) +
    scale_colour_manual(values = conc_cols) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    facet_grid(drug ~ day) +
    labs(x = "Culture net OD600", y = "Imaged cells with iniB on", shape = NULL,
         title = "Culture growth vs reporter in the imaged cells",
         subtitle = "Each point = one culture well and the imaging well made from it") +
    theme_Publication(base_size = 11) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          plot.subtitle = element_text(size = 8))
  p3b <- ggplot(link, aes(od_net, width, colour = conc)) +
    geom_point(aes(shape = flag), size = 2.2, alpha = 0.85, stroke = 0.9) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
    scale_colour_manual(values = conc_cols) +
    facet_grid(drug ~ day) +
    labs(x = "Culture net OD600", y = "Imaged cell width (\u00b5m, median)",
         title = "Culture growth vs imaged cell width") +
    theme_Publication(base_size = 11) +
    theme(legend.position = "bottom", legend.direction = "horizontal")
  p3 <- (p3a | p3b) + plot_layout(widths = c(1, 1.1))
  save_fig(p3, "03_culture_vs_imaged_cells", 15, 8.5)
} else {
  message("No reporter/reporter_per_well.csv - run reporter_response.R for the culture-vs-imaging panel")
}

message("\nDone. Tables in ", out_dir, " | figures in ", fig_dir)
