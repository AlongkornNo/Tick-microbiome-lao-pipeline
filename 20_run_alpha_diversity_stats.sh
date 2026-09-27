#!/bin/bash
# 20_run_alpha_diversity_stats.sh
#
# Per-specimen Chao1, Shannon, Simpson, and Pielou's evenness, plus
# global + pairwise Kruskal-Wallis (Benjamini-Hochberg FDR) by site --
# feeds Fig. 3A-D and Tables S5-S8.
#
# See alpha_diversity_stats.py's own docstring for verification details:
# this is a from-scratch reimplementation (not a copy of an original
# lab script -- none was available to review for this specific step),
# but it was checked against the manuscript's already-published Tables
# S5-S8 on the real rarefied feature table and matched exactly (all 60
# pairwise H/p/q values, 4 metrics x 15 site pairs). Report it to
# reviewers as independently-verified, not as the literal original code.
#
# NOTE: this does not cover Faith's PD, observed features, or Good's
# coverage (Methods: examined descriptively via boxplots only, not
# formally tested) -- those require the phylogeny-aware QIIME2 outputs
# from 14_run_core_diversity_metrics.sh / 15_run_additional_alpha_metrics.sh
# directly, not this script.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

feature_table_tsv="$DIVERSITY_DIR/dada2_table_rm_organelles_rarefied.tsv"

if [ ! -f "$feature_table_tsv" ]; then
  echo "ERROR: rarefied feature table (TSV) not found: $feature_table_tsv" >&2
  echo "Run 09_run_rarefy_normalization.sh first." >&2
  exit 1
fi
if [ ! -f "$METADATA" ]; then
  echo "ERROR: metadata file not found: $METADATA" >&2
  exit 1
fi

mkdir -p "$STATS_DIR"

python3 "$SCRIPT_DIR/alpha_diversity_stats.py" \
  --feature-table "$feature_table_tsv" \
  --metadata "$METADATA" \
  --group-column "$GROUP_COLUMN" \
  --outdir "$STATS_DIR"

echo ""
echo "Done. Output: $STATS_DIR"
echo "  alpha_diversity_per_specimen.tsv  (input to 21_plot_alpha_boxplots.R / Fig. 3)"
echo "  <metric>_pairwise_stats.tsv       (Tables S5-S8)"
