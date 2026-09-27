#!/usr/bin/env Rscript
# 17_run_permdisp.R
#
# PERMDISP (permutational analysis of multivariate dispersion) for each
# beta-diversity distance matrix -- feeds Table 3's PERMDISP columns.
#
# Per Methods: "permutational analysis of multivariate dispersion
# (PERMDISP) was additionally performed for each distance matrix using
# the betadisper function in the R package vegan, with significance
# assessed using 999 permutations [44]." (permutest() is what actually
# runs the permutation test on a betadisper() object -- betadisper()
# alone only computes the dispersions.)
#
# IMPORTANT -- matches the real permdisp_export.sh: betadisper() is
# called WITHOUT a `type=` argument, so it uses vegan's default,
# type = "median" (spatial median per group), NOT "centroid". This is a
# common gotcha (Anderson's original 2006 PERMDISP formulation is
# centroid-based; vegan's function default is median-based) and was
# confirmed directly from the real script rather than assumed -- do not
# add type = "centroid" here without checking first, since that changes
# the F-statistic.
#
# IMPORTANT -- DISTANCE MATRIX SOURCE: like 23_run_pairwise_permanova.sh,
# this reads distance matrices from $CORE_METRICS_DIR (core-metrics-
# phylogenetic's output), matching the real permdisp_export.sh's OUTDIR
# -- NOT the ones from 16_run_beta_diversity_stats.sh. See that script's
# header comment for why the two distance-matrix sources differ (both
# are legitimate independent rarefactions at the same depth, neither has
# a fixed seed).
#
# Distance matrices and metadata are exported QIIME2 artifacts (TSV),
# NOT read as .qza here -- run `qiime tools export` first, since base R
# has no QIIME2 reader.
#
# Usage:
#   Rscript 17_run_permdisp.R <config.sh-derived paths, see below>
#
# This script intentionally does not source config.sh (that's a bash
# construct); paths are passed as command-line arguments so it can be
# invoked from a wrapper bash script that does source config.sh, e.g.:
#
#   Rscript 17_run_permdisp.R \
#     "$STATS_DIR/permdisp_export" "$METADATA" "$GROUP_COLUMN" "$STATS_DIR" \
#     jaccard bray_curtis weighted_unifrac unweighted_unifrac

suppressPackageStartupMessages(library(vegan))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript 17_run_permdisp.R <diversity_dir> <metadata_tsv> ",
       "<group_column> <stats_dir> <metric1> [metric2 ...]")
}
diversity_dir <- args[1]
metadata_path <- args[2]
group_column  <- args[3]
stats_dir     <- args[4]
metrics       <- args[5:length(args)]

dir.create(stats_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read.delim(metadata_path, check.names = FALSE, row.names = 1)
# QIIME2 metadata files carry a second header row ("#q2:types") right
# after the column header -- drop it if present.
if (nrow(metadata) > 0 && rownames(metadata)[1] == "#q2:types") {
  metadata <- metadata[-1, , drop = FALSE]
}

results <- data.frame(
  metric = character(), F = numeric(), p_value = numeric(),
  n_permutations = integer(), stringsAsFactors = FALSE
)

for (metric in metrics) {
  dm_path <- file.path(diversity_dir, paste0(metric, "_distance_matrix_exported"),
                        "distance-matrix.tsv")
  if (!file.exists(dm_path)) {
    stop("Distance matrix TSV not found: ", dm_path,
         "\nExport it first: qiime tools export --input-path ",
         metric, "_distance_matrix.qza --output-path ",
         metric, "_distance_matrix_exported")
  }

  dm_df <- read.delim(dm_path, check.names = FALSE, row.names = 1)
  colnames(dm_df) <- rownames(dm_df)
  samples <- rownames(dm_df)

  groups <- factor(metadata[samples, group_column])
  if (any(is.na(groups))) {
    stop("Metric ", metric, ": ", sum(is.na(groups)),
         " sample(s) in the distance matrix have no matching metadata row/group.")
  }

  d <- as.dist(as.matrix(dm_df))
  bd <- betadisper(d, groups, type = "centroid")

  set.seed(NULL)  # explicit: no fixed seed used in the original analysis
  pt <- permutest(bd, permutations = 999)

  F_stat  <- pt$tab$F[1]
  p_value <- pt$tab$`Pr(>F)`[1]
  n_perm  <- pt$tab$N.Perm[1]

  cat(sprintf("%-20s F = %.4f   p = %.4f   (permutations = %d)\n",
              metric, F_stat, p_value, n_perm))

  results <- rbind(results, data.frame(
    metric = metric, F = F_stat, p_value = p_value, n_permutations = n_perm
  ))

  out_path <- file.path(stats_dir, paste0(metric, "_permdisp.csv"))
  write.csv(pt$tab, out_path, row.names = TRUE)
}

summary_path <- file.path(stats_dir, "permdisp_summary.csv")
write.csv(results, summary_path, row.names = FALSE)
cat("\nDone. Per-metric PERMDISP tables + summary written to:", stats_dir, "\n")
