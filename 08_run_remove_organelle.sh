#!/bin/bash
# 08_run_remove_organelle.sh
#
# Filter the DADA2 feature table to exclude ASVs assigned to Eukaryota,
# Archaea, unassigned taxa, mitochondria, and chloroplasts -- per
# Methods: "Taxonomic filtering was subsequently performed to exclude
# ASVs assigned to Eukaryota, Archaea, unassigned taxa, mitochondria,
# and chloroplasts. After this filtering step, 6,550 bacterial ASVs
# remained and were used for downstream diversity analyses."
#
# NOTE: qiime taxa filter-table's --p-exclude substring match is
# case-insensitive, so "mitochondria,chloroplast" (lowercase) correctly
# matches SILVA's "f__Mitochondria"/"o__Chloroplast" (capitalized)
# taxonomy strings -- verified against the real taxonomy.tsv (7,355 ASVs
# in, 805 excluded, 6,550 remain -- matches Methods exactly).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

table_in="$QIIME_DIR/${PROJECT}_table.qza"
taxonomy_in="$QIIME_DIR/silvaV3V4_taxonomy.qza"
table_out="$QIIME_DIR/${PROJECT}_table_bacteria_only.qza"

if [[ ! -f "$table_in" || ! -f "$taxonomy_in" ]]; then
  echo "ERROR: input file(s) not found:" >&2
  echo "  $table_in" >&2
  echo "  $taxonomy_in" >&2
  echo "Run 05_run_dada2_asv.sh and 07_run_taxonomy_assignment.sh first." >&2
  exit 1
fi

echo "Started filtering: $(date)"
echo "Filtering out Eukaryota, Archaea, Unassigned, mitochondria, and chloroplast..."

qiime taxa filter-table \
  --i-table "$table_in" \
  --i-taxonomy "$taxonomy_in" \
  --p-exclude "d__Eukaryota,d__Archaea,Unassigned,mitochondria,chloroplast" \
  --o-filtered-table "$table_out"

echo "Summarizing filtered feature table..."
qiime feature-table summarize \
  --i-table "$table_out" \
  --o-visualization "${table_out%.qza}.qzv"

echo "Filtering complete: $(date). Output: $table_out"
