#!/usr/bin/env Rscript
# 21_plot_alpha_boxplots.R
#
# Fig. 3A (Chao1), 3B (Pielou's evenness), 3C (Simpson), 3D (Shannon):
# alpha-diversity boxplots by site with overlaid jittered points.
#
# Per Methods: "Boxplots were generated in R using ggplot2 v3.5.2 with
# overlaid jittered points to visualize individual sample variation.
# Sample metadata were incorporated to stratify by experimental groups."
#
# Input: alpha_diversity_per_specimen.tsv from
# 20_run_alpha_diversity_stats.sh (columns: sample-id [index], Chao1,
# Shannon, Simpson, Pielou, Site).
#
# Usage:
#   Rscript 21_plot_alpha_boxplots.R <stats_dir> <plots_dir>

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 21_plot_alpha_boxplots.R <stats_dir> <plots_dir>")
}
stats_dir <- args[1]
plots_dir <- args[2]

dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

in_path <- file.path(stats_dir, "alpha_diversity_per_specimen.tsv")
if (!file.exists(in_path)) {
  stop("Not found: ", in_path, " -- run 20_run_alpha_diversity_stats.sh first.")
}
alpha <- read.delim(in_path, check.names = FALSE)
colnames(alpha)[1] <- "sample_id"

# Site order = first-appearance order in the table (mirrors metadata.tsv,
# i.e. Table 1 order), so panels/x-axis are consistent with every other
# figure/table in the manuscript.
site_order <- unique(alpha$Site)
alpha$Site <- factor(alpha$Site, levels = site_order)

plot_metric <- function(metric_col, y_label, panel_label, out_name) {
  p <- ggplot(alpha, aes(x = Site, y = .data[[metric_col]])) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", color = "black") +
    geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.7, color = "steelblue") +
    labs(x = "Sampling site", y = y_label, title = paste0(panel_label, ": ", y_label)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  out_path <- file.path(plots_dir, out_name)
  ggsave(out_path, p, width = 6, height = 4.5, dpi = 300)
  cat("Wrote:", out_path, "\n")
}

plot_metric("Chao1",   "Chao1 richness",       "Fig. 3A", "Fig3A_chao1_boxplot.png")
plot_metric("Pielou",  "Pielou's evenness",    "Fig. 3B", "Fig3B_evenness_boxplot.png")
plot_metric("Simpson", "Simpson's diversity",  "Fig. 3C", "Fig3C_simpson_boxplot.png")
plot_metric("Shannon", "Shannon diversity",    "Fig. 3D", "Fig3D_shannon_boxplot.png")

cat("\nDone. (Faith's PD, observed features, and Good's coverage are\n")
cat("descriptive-only per Methods -- not part of Fig. 3 -- and are not\n")
cat("plotted here; they come from QIIME2's own TSV exports in\n")
cat("CORE_METRICS_DIR if you need them for internal review.)\n")
