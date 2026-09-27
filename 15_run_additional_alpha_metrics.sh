#!/bin/bash
# 15_run_additional_alpha_metrics.sh
#
# Chao1, Simpson, and Good's coverage (all reported in the manuscript's
# alpha-diversity results) are NOT part of core-metrics-phylogenetic's
# fixed output set, so compute them separately here. Good's coverage is
# included since it is listed among the alpha metrics imported in the
# Methods ("Shannon, Simpson, Chao1, Faith's Phylogenetic Diversity,
# observed features, Pielou's evenness, and Good's coverage").
#
# CORRECTED input table, per the real group_stats.sh: these three metrics
# are computed on the STANDALONE rarefied table from
# 09_run_rarefy_normalization.sh ($QIIME_DIR/${PROJECT}_table_rarefied.qza,
# named dada2_table_rm_organelles_rarefied.qza in that script), NOT on
# $CORE_METRICS_DIR/rarefied_table.qza (core-metrics-phylogenetic's own
# internal resampling, step 14) as an earlier version of this script
# assumed. This is a genuine, confirmed inconsistency worth flagging: it
# means Chao1/Simpson/Good's-coverage and Shannon/Faith's-PD/Observed-
# Features/Evenness come from TWO DIFFERENT random rarefaction draws at
# the same depth (neither `qiime feature-table rarefy` nor
# core-metrics-phylogenetic's internal step has a fixed seed) -- not all
# 7 alpha-diversity metrics in the manuscript share one resampling, despite
# what the header comment here previously claimed.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

rarefied_table="$QIIME_DIR/${PROJECT}_table_rarefied.qza"
if [ ! -f "$rarefied_table" ]; then
  echo "ERROR: rarefied table not found: $rarefied_table" >&2
  echo "Run 09_run_rarefy_normalization.sh first." >&2
  exit 1
fi

echo "qiime2 version: $(qiime --version 2>&1 | head -1)"

for metric in chao1 simpson goods_coverage; do
  echo "Computing $metric..."
  qiime diversity alpha \
    --i-table "$rarefied_table" \
    --p-metric "$metric" \
    --o-alpha-diversity "$CORE_METRICS_DIR/${metric}_vector.qza"
done

echo "Done. Additional alpha metrics written to: $CORE_METRICS_DIR"
