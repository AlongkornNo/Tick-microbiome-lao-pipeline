#!/bin/bash
# 16_run_beta_diversity_stats.sh
#
# Beta-diversity distance matrices + global PERMANOVA (Table 3), rewritten
# from the original per-metric scripts (run_jaccard_adonis.sh,
# bray_adonis_results.sh, run_weighted_unifrac_adonis.sh,
# run_unweighted_unifrac_adonis.sh) into one loop over $BETA_METRICS.
#
# Per Methods: "Beta diversity ... were evaluated using permutational
# multivariate analysis of variance (PERMANOVA ... 999 permutations)."
# Table 3's footnote: "PERMANOVA was performed using qiime diversity
# adonis (vegan::adonis2 in R)."
#
# IMPORTANT: distance matrices are computed directly on the already-
# rarefied table from 09_run_rarefy_normalization.sh via standalone
# `qiime diversity beta` / `beta-phylogenetic` -- NOT via
# core-metrics-phylogenetic (see the note in
# 14_run_core_diversity_metrics.sh for why that distinction matters).
#
# PERMDISP (betadisper/permutest, R package vegan) is NOT computed here --
# qiime diversity adonis does not provide it, and Methods specifies the
# R functions directly. See 17_run_permdisp.R for that step.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

table="$QIIME_DIR/${PROJECT}_table_rarefied.qza"
tree="$PHYLOGENY_DIR/${PROJECT}_rooted_tree.qza"

if [ ! -f "$table" ]; then
  echo "ERROR: rarefied table not found: $table" >&2
  echo "Run 09_run_rarefy_normalization.sh first." >&2
  exit 1
fi
if [ ! -f "$tree" ]; then
  echo "ERROR: rooted tree not found: $tree" >&2
  echo "Run 13_run_build_phylogeny.sh first." >&2
  exit 1
fi
if [ ! -f "$METADATA" ]; then
  echo "ERROR: metadata file not found: $METADATA" >&2
  exit 1
fi

mkdir -p "$DIVERSITY_DIR" "$STATS_DIR"

echo "qiime2 version: $(qiime --version 2>&1 | head -1)"

# Map this pipeline's metric names (config.sh BETA_METRICS, matching
# core-metrics-phylogenetic's file-stem convention) to the --p-metric
# values QIIME2's `beta`/`beta-phylogenetic` commands actually expect,
# and to whether the metric needs the phylogenetic tree.
metric_param() {
  case "$1" in
    bray_curtis) echo "braycurtis" ;;
    jaccard) echo "jaccard" ;;
    weighted_unifrac) echo "weighted_unifrac" ;;
    unweighted_unifrac) echo "unweighted_unifrac" ;;
    *) echo "ERROR: unknown beta metric '$1'" >&2; exit 1 ;;
  esac
}
is_phylogenetic() {
  case "$1" in
    weighted_unifrac|unweighted_unifrac) return 0 ;;
    *) return 1 ;;
  esac
}

for metric in "${BETA_METRICS[@]}"; do
  echo ""
  echo "=== $metric ==="

  dm_qza="$DIVERSITY_DIR/${metric}_distance_matrix.qza"
  adonis_qzv="$STATS_DIR/${metric}_adonis.qzv"

  if is_phylogenetic "$metric"; then
    echo "Computing $metric distance (phylogenetic)..."
    qiime diversity beta-phylogenetic \
      --i-table "$table" \
      --i-phylogeny "$tree" \
      --p-metric "$(metric_param "$metric")" \
      --o-distance-matrix "$dm_qza"
  else
    echo "Computing $metric distance..."
    qiime diversity beta \
      --i-table "$table" \
      --p-metric "$(metric_param "$metric")" \
      --o-distance-matrix "$dm_qza"
  fi

  echo "Running global PERMANOVA (adonis, $GROUP_COLUMN, 999 permutations)..."
  qiime diversity adonis \
    --i-distance-matrix "$dm_qza" \
    --m-metadata-file "$METADATA" \
    --p-formula "$GROUP_COLUMN" \
    --p-permutations 999 \
    --o-visualization "$adonis_qzv"

  echo "Exporting adonis result table..."
  qiime tools export \
    --input-path "$adonis_qzv" \
    --output-path "$STATS_DIR/${metric}_adonis_exported"
done

echo ""
echo "Done. Distance matrices: $DIVERSITY_DIR/*_distance_matrix.qza"
echo "      Global PERMANOVA:  $STATS_DIR/*_adonis.qzv (+ exported/ CSVs)"
echo "Pairwise post hoc PERMANOVA per site pair (Tables S9-S12) and PERMDISP"
echo "(Table 3) use a DIFFERENT distance-matrix source (core-metrics-"
echo "phylogenetic's own output, \$CORE_METRICS_DIR) -- see"
echo "23_run_pairwise_permanova.sh and 17_run_permdisp.R, and the note in"
echo "14_run_core_diversity_metrics.sh explaining why the sources differ."
