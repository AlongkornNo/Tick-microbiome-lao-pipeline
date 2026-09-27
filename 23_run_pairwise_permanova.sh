#!/bin/bash
# 23_run_pairwise_permanova.sh
#
# Pairwise (site-vs-site) PERMANOVA for all four beta-diversity metrics --
# feeds Tables S9-S12. Rewritten from the original run_permanova.sh:
# hardcoded macOS paths replaced with config.sh variables, loop driven by
# $BETA_METRICS, and (see FIX note below) --p-pairwise added.
#
# Uses `qiime diversity beta-group-significance --p-method permanova`
# (scikit-bio's permanova), which is a DIFFERENT QIIME2/tool path than
# Table 3's global PERMANOVA (`qiime diversity adonis`, vegan::adonis2 --
# see 16_run_beta_diversity_stats.sh). Both implement the same Anderson
# (2001) pseudo-F statistic, so that tool difference alone would not
# cause any numeric mismatch between the two for the SAME input distance
# matrix.
#
# FIX: the original run_permanova.sh did not pass --p-pairwise, so as
# written it would only produce the omnibus (all 6 groups at once) test,
# not the site-by-site pairwise table Tables S9-S12 actually need. Added
# --p-pairwise below -- confirmed correct by a second, independent real
# script (group_stats.sh), whose beta-group-significance loop does
# include --p-pairwise for all four metrics.
#
# IMPORTANT -- DISTANCE MATRIX SOURCE: this script's original OUTDIR
# pointed at core-metrics-results (i.e. the output of
# `qiime diversity core-metrics-phylogenetic`, see
# 14_run_core_diversity_metrics.sh), NOT the distance matrices from
# 16_run_beta_diversity_stats.sh (standalone `qiime diversity beta` /
# `beta-phylogenetic` calls on the separately-rarefied table). Those are
# two independent random subsamples at the same depth (neither
# `qiime feature-table rarefy` nor core-metrics-phylogenetic's internal
# rarefaction step has a fixed seed) -- confirmed, by directly comparing
# the two real script sets, to be the reason a from-scratch PERMANOVA on
# one set of matrices reproduces Tables S9-S12 (pairwise, this script's
# source) exactly but does not exactly reproduce Table 3 (global, script
# 16's source). This script intentionally keeps using
# $CORE_METRICS_DIR, matching the original -- do not "fix" it to use
# script 16's matrices without checking first, since that would break
# the exact match to Tables S9-S12 that we've already verified.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

export_base="$STATS_DIR/beta_group_significance_export"
mkdir -p "$export_base"

if [ ! -f "$METADATA" ]; then
  echo "ERROR: metadata file not found: $METADATA" >&2
  exit 1
fi

echo "=== SAMPLE SIZE PER GROUP ($GROUP_COLUMN) ==="
awk -v col="$GROUP_COLUMN" '
  NR==1 { for (i=1;i<=NF;i++) if ($i==col) c=i; next }
  $1=="#q2:types" { next }
  { count[$c]++ }
  END { for (g in count) print g, count[g] }
' FS='\t' "$METADATA" | sort
echo "=============================="
echo ""

for metric in "${BETA_METRICS[@]}"; do
  dist_matrix="$CORE_METRICS_DIR/${metric}_distance_matrix.qza"
  export_dir="$export_base/${metric}_permanova"
  mkdir -p "$export_dir"

  if [[ ! -f "$dist_matrix" ]]; then
    echo "WARNING: distance matrix not found: $dist_matrix -- skipping $metric" >&2
    continue
  fi

  echo "Running pairwise PERMANOVA for: $metric ..."
  echo "Using file: $dist_matrix"

  qiime diversity beta-group-significance \
    --i-distance-matrix "$dist_matrix" \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$GROUP_COLUMN" \
    --p-method permanova \
    --p-pairwise \
    --p-permutations 999 \
    --o-visualization "$export_dir/${metric}_permanova.qzv"

  qiime tools export \
    --input-path "$export_dir/${metric}_permanova.qzv" \
    --output-path "$export_dir"

  echo "Done: $export_dir"
  echo ""
done

echo "=== ALL PAIRWISE PERMANOVA COMPLETED ==="
echo "Results exported to: $export_base"
