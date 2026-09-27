#!/bin/bash
# 24_run_alpha_group_significance.sh
#
# QIIME2-native alpha-diversity group significance testing, rewritten
# from the real group_stats.sh (alpha section) -- hardcoded macOS paths
# replaced with config.sh variables, driven by $ALPHA_METRICS.
#
# This is the QIIME2-native route (`qiime diversity alpha-group-
# significance`, skbio-based Kruskal-Wallis) to the same result Tables
# S5-S8 report. It is a DIFFERENT tool than alpha_diversity_stats.py
# (20_run_alpha_diversity_stats.sh), the from-scratch Python
# reimplementation already verified against the published Tables S5-S8
# (60/60 exact match: H, p, q for all 4 formally-tested metrics x 15 site
# pairs). Keep both: this one is the more literal reproduction of the
# real script; that one is the independently-verified cross-check.
#
# All 7 vector.qza inputs are assumed already built by
# 14_run_core_diversity_metrics.sh (shannon, observed_features, evenness,
# faith_pd) and 15_run_additional_alpha_metrics.sh (chao1, simpson,
# goods_coverage) -- both write into $CORE_METRICS_DIR, even though (see
# 14's header note) they draw on two DIFFERENT underlying rarefactions.
# This script does not recompute them.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

if [ ! -f "$METADATA" ]; then
  echo "ERROR: metadata file not found: $METADATA" >&2
  exit 1
fi

group_stats_dir="$STATS_DIR/alpha_group_significance"
mkdir -p "$group_stats_dir"

echo "qiime2 version: $(qiime --version 2>&1 | head -1)"

for metric in "${ALPHA_METRICS[@]}"; do
  vector_qza="$CORE_METRICS_DIR/${metric}_vector.qza"
  if [ ! -f "$vector_qza" ]; then
    echo "WARNING: not found, skipping $metric: $vector_qza" >&2
    echo "  (run 14_run_core_diversity_metrics.sh / 15_run_additional_alpha_metrics.sh first)" >&2
    continue
  fi

  echo "Processing alpha metric: $metric ..."

  qiime diversity alpha-group-significance \
    --i-alpha-diversity "$vector_qza" \
    --m-metadata-file "$METADATA" \
    --o-visualization "$group_stats_dir/${metric}-group-significance.qzv"

  export_tmp="$group_stats_dir/${metric}_tsv"
  mkdir -p "$export_tmp"
  qiime tools export \
    --input-path "$vector_qza" \
    --output-path "$export_tmp"
  mv "$export_tmp/alpha-diversity.tsv" "$group_stats_dir/${metric}_vector.tsv"
  rm -r "$export_tmp"
done

echo ""
echo "Done. Results saved in: $group_stats_dir"
