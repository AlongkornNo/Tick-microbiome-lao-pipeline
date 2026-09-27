#!/bin/bash
# 18_export_taxa_top10.sh
#
# Export top-10 (+ "Other") true relative-abundance tables at Phylum and
# Family level, for the stacked bar charts in Fig. 2A/2B.
#
# Per Methods: "Microbial composition at the phylum and family levels was
# visualized using stacked bar charts of relative abundance, with the top
# 10 taxa per level exported as CSV files via Python v3.11. All remaining,
# lower-abundance taxa at each level were grouped into a single 'Other'
# category; the top 10 taxa and the 'Other' category together reflect the
# true relative abundance of each specimen and were not renormalized to
# sum to 100% among the top 10 taxa alone."
#
# Uses summarize_taxonomy.py (this directory), which computes relative
# abundance against each specimen's TRUE total read count (all taxa at
# that rank), not just the sum of the top 10 -- this was verified against
# the manuscript's Table S3 (exact match, all 30 specimens, Phylum and
# Family) before being adopted into this pipeline. See that script's
# docstring for what bug this avoids.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

feature_table_tsv="$QIIME_DIR/barplot_export/feature-table.tsv"
taxonomy_tsv="$QIIME_DIR/barplot_export/taxonomy_export/taxonomy.tsv"

if [[ ! -f "$feature_table_tsv" || ! -f "$taxonomy_tsv" ]]; then
  echo "ERROR: expected exports not found:" >&2
  echo "  $feature_table_tsv" >&2
  echo "  $taxonomy_tsv" >&2
  echo "Run 12_export_barplot_csv.sh first (it produces both)." >&2
  exit 1
fi

mkdir -p "$TAXA_SUMMARY_DIR"

python3 "$SCRIPT_DIR/summarize_taxonomy.py" \
  --feature-table "$feature_table_tsv" \
  --taxonomy "$taxonomy_tsv" \
  --rank Phylum,Family \
  --top-n "$TOP_N_BARPLOT" \
  --outdir "$TAXA_SUMMARY_DIR"

echo ""
echo "Done. Top-$TOP_N_BARPLOT (+Other) relative-abundance tables written to: $TAXA_SUMMARY_DIR"
echo "  top${TOP_N_BARPLOT}_phylum_relative_abundance.csv  (Fig. 2A)"
echo "  top${TOP_N_BARPLOT}_family_relative_abundance.csv  (Fig. 2B)"
