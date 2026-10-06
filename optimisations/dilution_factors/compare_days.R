# =============================================================================
# Day 2 vs day 4 of drug exposure (Plate60, iniB-mScarlet): side-by-side
# comparison of the per-day outputs of dilution_optimisation.R and
# drug_genetic_neighbours.R (all dilutions pooled as replicate wells).
#
# Reads <data_dir>/<day>/... for each day; writes <data_dir>/day_comparison/.
# Run both per-day scripts first.  From the project root:
#   Rscript optimisations/dilution_factors/compare_days.R [data_dir] [day2,day4]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

if (file.exists("Theme.R")) source("Theme.R") else
  stop("Theme.R not found - run from project root")

args     <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "optimisations/dilution_factors"
days     <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("day2", "day4")
NB_RUN   <- "all_dilutions"
TOPK     <- 8L

out_dir <- file.path(data_dir, "day_comparison")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(p, name, w, h, dpi = 300) {
  # cairo is unavailable without XQuartz; quartz (pdf) + ragg (png) instead
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white",
         device = function(filename, ...) grDevices::quartz(type = "pdf", file = filename, ...))
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = dpi, bg = "white", device = ragg::agg_png)
}

read_days <- function(...) rbindlist(lapply(days, function(dy) {
  f <- file.path(data_dir, dy, ...)
  if (!file.exists(f)) stop("Missing ", f, " - run the per-day scripts first")
  fread(f)[, day := dy]
}), fill = TRUE)

conc_levels <- c("0.25xMIC", "0.5xMIC", "1xMIC")
dil_levels  <- c(20L, 50L, 100L, 150L, 200L)
anchor_labs <- c(libNT = "Library-NT", plate = "Drug-plate")
day_cols    <- setNames(c("#386cb0", "#ef3b2c", "#7fc97f", "#fdb462")[seq_along(days)], days)

# ---- 1. Dilution: usable cells per FOV --------------------------------------
wl <- read_days("dilution_summary_per_condition.csv")
wl[, dil_f := factor(paste0("1:", dil), levels = paste0("1:", dil_levels))]
dil_sum <- wl[, .(cells_fov = median(cells_fov), iso_fov = median(iso_fov),
                  frac_touch = median(frac_touch), n_best = sum(best)), by = .(day, dil)]
setorder(dil_sum, day, dil)
fwrite(dil_sum, file.path(out_dir, "dilution_by_day.csv"))
message("Isolated cells / FOV (median across drug x dose) and # conditions where best:")
print(dcast(dil_sum, dil ~ day, value.var = c("iso_fov", "n_best")), digits = 3)

p1a <- ggplot(wl, aes(dil_f, iso_fov, colour = day)) +
  geom_line(aes(group = interaction(day, condition)), alpha = 0.2) +
  stat_summary(aes(group = day), fun = median, geom = "line", linewidth = 1.2) +
  stat_summary(aes(group = day), fun = median, geom = "point", size = 2.6) +
  scale_colour_manual(values = day_cols) +
  scale_y_log10() +
  labs(x = "Dilution", y = "Isolated cells per FOV", title = "Usable cells per FOV",
       subtitle = "Thin = each drug x dose; thick = median") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))
p1b <- ggplot(wl, aes(dil_f, frac_touch, colour = day)) +
  geom_line(aes(group = interaction(day, condition)), alpha = 0.2) +
  stat_summary(aes(group = day), fun = median, geom = "line", linewidth = 1.2) +
  stat_summary(aes(group = day), fun = median, geom = "point", size = 2.6) +
  scale_colour_manual(values = day_cols) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(x = "Dilution", y = "Cells touching a neighbour", title = "Crowding",
       subtitle = "Edge-to-edge gap < 0.3 µm") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))
p1 <- (p1a | p1b) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")
save_fig(p1, "01_dilution_by_day", 10, 4.6)

# ---- 2. Phenotype strength and expected-pathway AUC -------------------------
nbp <- function(f) file.path("drug_neighbours", NB_RUN, f)
st <- read_days(nbp("phenotype_strength.csv"))
pw <- read_days(nbp("expected_pathway_check.csv"))
add_labels <- function(x) {
  x[, `:=`(drug = sub("_.*", "", cond),
           conc = factor(sub(".*_", "", cond), levels = conc_levels),
           anchor_lab = factor(anchor_labs[anchor], levels = anchor_labs))]
}
st <- add_labels(st)
pw <- add_labels(pw)
pw[, drug_lab := paste0(drug, " -> ", sub(" /.*|/.*", "", expected_pathway))]

p2 <- ggplot(st, aes(conc, dist_to_NT / typical_NT_NT, colour = drug,
                     linetype = day, shape = day, group = interaction(drug, day))) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
  scale_colour_Publication() +
  facet_wrap(~ anchor_lab, scales = "free_y") +
  labs(x = "Dose", y = "Distance to NT / typical NT-NT", title = "Phenotype strength by day",
       subtitle = "1 = indistinguishable from NT-to-NT spread") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))
save_fig(p2, "02_phenotype_strength_by_day", 10, 4.8)

p3 <- ggplot(pw, aes(conc, pathway_auc, colour = day, linetype = anchor_lab,
                     group = interaction(day, anchor_lab))) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.3)) +
  geom_point(aes(shape = p_bh < 0.05), size = 2.8, position = position_dodge(width = 0.3)) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                     labels = c(`FALSE` = "BH q >= 0.05", `TRUE` = "BH q < 0.05")) +
  scale_colour_manual(values = day_cols) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ drug_lab, nrow = 1) +
  labs(x = "Dose", y = "Expected-pathway AUC", title = "Target-pathway proximity by day",
       subtitle = "AUC = P(pathway knockdown closer than a non-pathway strain); 0.5 = chance; solid = Library-NT, dashed = Drug-plate") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))
save_fig(p3, "03_pathway_auc_by_day", 11.5, 4.6)

pw_wide <- dcast(pw, anchor + cond ~ day, value.var = c("pathway_auc", "p_bh", "wells_above_chance"))
fwrite(pw_wide, file.path(out_dir, "pathway_auc_by_day.csv"))
message("\nExpected-pathway AUC by day:")
print(pw_wide, digits = 2)

# ---- 3. Do the nearest knockdowns persist from day to day? ------------------
nb <- read_days(nbp("drug_neighbours_all_ranks.csv"))
top <- nb[rank <= TOPK, .(set = list(neighbour), top3 = paste(neighbour[1:3], collapse = ", ")),
          by = .(day, anchor, cond)]
pairs <- combn(days, 2, simplify = FALSE)
persist <- rbindlist(lapply(pairs, function(pr) {
  a <- top[day == pr[1]]; b <- top[day == pr[2]]
  m <- merge(a, b, by = c("anchor", "cond"), suffixes = c("_a", "_b"))
  m[, .(anchor, cond, comparison = paste(pr, collapse = " vs "),
        n_shared = mapply(function(x, y) length(intersect(x, y)), set_a, set_b),
        shared = mapply(function(x, y) paste(intersect(x, y), collapse = ", "), set_a, set_b),
        top3_a = top3_a, top3_b = top3_b)]
}))
# Rank agreement over all 67 strains
rk <- dcast(nb, anchor + cond + neighbour ~ day, value.var = "rank")
persist <- merge(persist, rk[, .(spearman = cor(get(days[1]), get(days[2]), method = "spearman")),
                             by = .(anchor, cond)], by = c("anchor", "cond"))
setnames(persist, c("top3_a", "top3_b"), paste0("top3_", days[1:2]))
fwrite(persist, file.path(out_dir, "neighbour_persistence.csv"))
message("\nTop-", TOPK, " neighbours shared between days:")
print(persist[, -"comparison"], digits = 2)

persist[, `:=`(drug = sub("_.*", "", cond),
               conc = factor(sub(".*_", "", cond), levels = conc_levels),
               anchor_lab = factor(anchor_labs[anchor], levels = anchor_labs))]
p4 <- ggplot(persist, aes(conc, n_shared, fill = drug)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single"), width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", spearman), group = drug, y = n_shared + 0.35),
            position = position_dodge(width = 0.8), size = 2.4) +
  scale_fill_Publication() +
  scale_y_continuous(limits = c(0, TOPK + 0.8), breaks = 0:TOPK) +
  facet_wrap(~ anchor_lab) +
  labs(x = "Dose", y = sprintf("Top-%d neighbours shared", TOPK),
       title = sprintf("Neighbour persistence, %s vs %s", days[1], days[2]),
       subtitle = "Bars = knockdowns in both days' top 8; labels = Spearman correlation of all 67 ranks") +
  theme_Publication(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.subtitle = element_text(size = 8))
save_fig(p4, "04_neighbour_persistence", 10, 4.8)

message("\nDone. Tables in ", out_dir, " | figures in ", fig_dir)
