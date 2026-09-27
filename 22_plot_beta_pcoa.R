#!/usr/bin/env Rscript
# 22_plot_beta_pcoa.R
#
# Fig. 4A (Jaccard), 4B (Bray-Curtis), 4C (weighted UniFrac), 4D
# (unweighted UniFrac): PCoA ordinations with global PERMANOVA annotated.
#
# Per Methods: "PCoA coordinates were calculated with the ape v5.8.1
# package, and plots were generated using ggplot2 v3.5.2, coloring
# samples by experimental group. PERMANOVA results were included in the
# PCoA plots to indicate statistical significance of group separation,
# along with PERMDISP results." (FIX: an earlier version of this script
# only annotated PERMANOVA -- Methods explicitly requires PERMDISP too.)
#
# PERMDISP here uses betadisper() with NO `type=` argument (vegan's
# default, type = "median" -- confirmed from the real permdisp_export.sh,
# see 17_run_permdisp.R's header note), for consistency with that script.
#
# Input: the *_distance_matrix.qza files from 16_run_beta_diversity_stats.sh,
# exported to TSV first (qiime tools export --input-path X.qza
# --output-path X_exported ; the matrix is at X_exported/distance-matrix.tsv).
#
# NOTE: this script re-runs vegan::adonis2 and betadisper()/permutest()
# itself (rather than parsing exported .qzv results) to get the
# annotation numbers, using the SAME distance matrix loaded for the plot
# (from 16_run_beta_diversity_stats.sh's output) -- so PERMANOVA and
# PERMDISP shown here are mutually consistent with each other and with
# the plotted points. They are NOT guaranteed to exactly equal Table 3's
# published PERMANOVA/PERMDISP values, though: Table 3's PERMDISP (see
# 17_run_permdisp.R) and the real pairwise scripts both draw distance
# matrices from $CORE_METRICS_DIR, a different underlying rarefaction
# draw than this script's input (see 14_run_core_diversity_metrics.sh's
# header for why -- neither rarefaction step has a fixed seed). Treat a
# mismatch against Table 3 as expected/already-explained, not a new bug.
#
# NOTE ON % VARIANCE EXPLAINED (VERIFIED): these four distance metrics are
# not all strictly Euclidean -- weighted UniFrac in particular has several
# small negative eigenvalues. ape::pcoa() is called here with NO
# Cailliez/Lingoes correction. Its `Relative_eig` column expresses each
# eigenvalue as a fraction of the sum of POSITIVE eigenvalues only (not
# all eigenvalues) -- this was checked directly against the real distance
# matrices in Python (classic double-centered PCoA, same eigendecomposition
# ape uses internally) and reproduces the manuscript's reported PC1/PC2
# percentages exactly for all four metrics (Fig. 4A: 9.9%/8.9%; 4B:
# 17.5%/13.6%; 4C: 37.5%/18.5%; 4D: 18.8%/12.1%) -- including weighted
# UniFrac, where dividing by the sum of ALL eigenvalues instead gives
# 38.3%/18.9%, which does NOT match. So: no correction, use Relative_eig
# as ape computes it by default. Do not change this without re-checking
# against the real data.
#
# Usage:
#   Rscript 22_plot_beta_pcoa.R <diversity_dir> <metadata_tsv> \
#     <group_column> <plots_dir> <metric1> [metric2 ...]

suppressPackageStartupMessages({
  library(ape)
  library(vegan)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript 22_plot_beta_pcoa.R <diversity_dir> <metadata_tsv> ",
       "<group_column> <plots_dir> <metric1> [metric2 ...]")
}
diversity_dir <- args[1]
metadata_path <- args[2]
group_column  <- args[3]
plots_dir     <- args[4]
metrics       <- args[5:length(args)]

dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read.delim(metadata_path, check.names = FALSE, row.names = 1)
if (nrow(metadata) > 0 && rownames(metadata)[1] == "#q2:types") {
  metadata <- metadata[-1, , drop = FALSE]
}

panel_label <- c(jaccard = "Fig. 4A", bray_curtis = "Fig. 4B",
                  weighted_unifrac = "Fig. 4C", unweighted_unifrac = "Fig. 4D")
display_name <- c(jaccard = "Jaccard", bray_curtis = "Bray-Curtis",
                   weighted_unifrac = "Weighted UniFrac",
                   unweighted_unifrac = "Unweighted UniFrac")

for (metric in metrics) {
  dm_path <- file.path(diversity_dir, paste0(metric, "_distance_matrix_exported"),
                        "distance-matrix.tsv")
  if (!file.exists(dm_path)) {
    stop("Not found: ", dm_path,
         "\nExport first: qiime tools export --input-path ", metric,
         "_distance_matrix.qza --output-path ", metric, "_distance_matrix_exported")
  }

  dm_df <- read.delim(dm_path, check.names = FALSE, row.names = 1)
  colnames(dm_df) <- rownames(dm_df)
  samples <- rownames(dm_df)
  groups <- factor(metadata[samples, group_column])

  d <- as.dist(as.matrix(dm_df))

  # --- PCoA ---
  pc <- pcoa(d)
  rel_eig <- pc$values$Relative_eig
  pc1_pct <- round(100 * rel_eig[1], 1)
  pc2_pct <- round(100 * rel_eig[2], 1)

  coords <- as.data.frame(pc$vectors[, 1:2])
  colnames(coords) <- c("PC1", "PC2")
  coords$Group <- groups

  # --- Global PERMANOVA (for the annotation) ---
  ad <- adonis2(d ~ groups, permutations = 999)
  F_stat <- round(ad$F[1], 3)
  R2 <- round(ad$R2[1], 3)
  p_val <- ad$`Pr(>F)`[1]
  p_label <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)

  # --- PERMDISP (for the annotation) ---
  bd <- betadisper(d, groups)  # default type = "median", see note above
  pt <- permutest(bd, permutations = 999)
  disp_F <- round(pt$tab$F[1], 3)
  disp_p <- pt$tab$`Pr(>F)`[1]
  disp_p_label <- if (disp_p < 0.001) "p < 0.001" else sprintf("p = %.3f", disp_p)

  subtitle <- sprintf(
    "PERMANOVA: pseudo-F = %.3f, R² = %.3f, %s | PERMDISP: F = %.3f, %s (999 perm.)",
    F_stat, R2, p_label, disp_F, disp_p_label
  )

  label <- if (metric %in% names(panel_label)) panel_label[metric] else metric
  name <- if (metric %in% names(display_name)) display_name[metric] else metric

  p <- ggplot(coords, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 2.5, alpha = 0.85) +
    stat_ellipse(level = 0.68, linewidth = 0.4) +
    labs(
      title = paste0(label, ": ", name, " PCoA"),
      subtitle = subtitle,
      x = paste0("PC1 (", pc1_pct, "%)"),
      y = paste0("PC2 (", pc2_pct, "%)")
    ) +
    theme_bw(base_size = 11)

  out_path <- file.path(plots_dir, paste0(gsub("Fig\\. ", "Fig", label), "_", metric, "_pcoa.png"))
  ggsave(out_path, p, width = 6.5, height = 5, dpi = 300)
  cat(sprintf("%-20s PC1=%.1f%% PC2=%.1f%%  PERMANOVA F=%.3f R2=%.3f p=%.4f  PERMDISP F=%.3f p=%.4f  -> %s\n",
              metric, pc1_pct, pc2_pct, F_stat, R2, p_val, disp_F, disp_p, out_path))
}

cat("\nDone.\n")
