#!/bin/bash
# 14_run_core_diversity_metrics.sh
#
# Compute the standard phylogeny-aware alpha/beta diversity suite in one
# QIIME2 call -- per Methods: "Phylogeny-aware diversity analyses were
# performed using the core-metrics-phylogenetic pipeline with a sampling
# depth of 1,800 reads per sample matching the rarefaction threshold."
#   alpha: Shannon, Observed Features, Pielou's Evenness, Faith's PD
#   beta:  unweighted/weighted UniFrac, Jaccard, Bray-Curtis
#
# This command performs ITS OWN internal rarefaction at --p-sampling-depth,
# so it takes the pre-rarefaction, organelle-filtered table (step 08's
# output) as input -- NOT the already-rarefied table from step 09/10.
#
# IMPORTANT -- NOT ALL ALPHA/BETA RESULTS IN THE MANUSCRIPT SHARE THIS ONE
# RESAMPLING. Confirmed directly from the real analysis scripts (not
# assumed), there are actually THREE separate rarefactions in play, all at
# the same depth but each an independent random draw (neither
# `qiime feature-table rarefy` nor core-metrics-phylogenetic's internal
# step has a fixed seed):
#   1. THIS command's internal resampling -- source for Shannon, Observed
#      Features, Evenness, Faith's PD (alpha), AND for the pairwise
#      PERMANOVA (Tables S9-S12, see 23_run_pairwise_permanova.sh) and
#      PERMDISP (Table 3's PERMDISP columns, see 17_run_permdisp.R) --
#      both of those read distance matrices from $CORE_METRICS_DIR,
#      i.e. from here.
#   2. The standalone rarefied table from 09_run_rarefy_normalization.sh
#      ($QIIME_DIR/${PROJECT}_table_rarefied.qza) -- source for Chao1,
#      Simpson, Good's coverage (see 15_run_additional_alpha_metrics.sh,
#      corrected to use this table) AND for Table 3's GLOBAL PERMANOVA
#      pseudo-F/R²/p (see 16_run_beta_diversity_stats.sh, which computes
#      its own distance matrices from this same table via standalone
#      `qiime diversity beta`/`beta-phylogenetic`).
#   3. (Not yet encountered as a separate draw beyond the above two.)
# Net effect: Table 3's global PERMANOVA (source 2) is NOT computed from
# the same distance matrices as its own PERMDISP column or as Tables
# S9-S12 (both source 1) -- this is confirmed, not speculative, and is
# the resolved explanation for why a from-scratch recomputation of the
# global PERMANOVA using source-1 matrices doesn't exactly match Table 3.
# Similarly, Chao1/Simpson/Good's-coverage (source 2) don't share a
# resampling with Shannon/Faith's-PD/Observed-Features/Evenness (source 1),
# despite reflecting "the same" rarefaction depth and being reported
# together in Fig. 3 as if from one pass.
#
# This command's own *_distance_matrix.qza outputs ARE used (source 1,
# above) for Tables S9-S12 and PERMDISP -- just not for Table 3's global
# PERMANOVA or the PCoA plots (Fig. 4), which need
# 16_run_beta_diversity_stats.sh's separately-computed matrices instead.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

table="$QIIME_DIR/${PROJECT}_table_bacteria_only.qza"
tree="$PHYLOGENY_DIR/${PROJECT}_rooted_tree.qza"

if [[ ! -f "$table" || ! -f "$tree" ]]; then
  echo "ERROR: input file(s) not found:" >&2
  echo "  $table" >&2
  echo "  $tree" >&2
  echo "Run 08_run_remove_organelle.sh and 13_run_build_phylogeny.sh first." >&2
  exit 1
fi
if [ ! -f "$METADATA" ]; then
  echo "ERROR: metadata file not found: $METADATA" >&2
  exit 1
fi

if [ -d "$CORE_METRICS_DIR" ]; then
  echo "ERROR: output directory already exists: $CORE_METRICS_DIR" >&2
  echo "QIIME2 refuses to write into an existing directory -- remove or rename it first." >&2
  exit 1
fi

echo "qiime2 version: $(qiime --version 2>&1 | head -1)"
echo "Running core-metrics-phylogenetic (sampling depth: $RAREFACTION_DEPTH)..."

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$tree" \
  --i-table "$table" \
  --p-sampling-depth "$RAREFACTION_DEPTH" \
  --m-metadata-file "$METADATA" \
  --output-dir "$CORE_METRICS_DIR"

echo "Done. Core diversity metrics: $CORE_METRICS_DIR"
echo "  Alpha: shannon_vector.qza, observed_features_vector.qza, evenness_vector.qza, faith_pd_vector.qza"
echo "  Beta:  unweighted_unifrac_distance_matrix.qza, weighted_unifrac_distance_matrix.qza, jaccard_distance_matrix.qza, bray_curtis_distance_matrix.qza"
