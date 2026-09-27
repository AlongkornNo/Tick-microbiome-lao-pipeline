#!/bin/bash
# 25_export_genus_top30.sh
#
# Genus-level top-30 relative abundance -- NOT the same output as the
# Phylum/Family top-10 bar charts (18/19). Per Methods: "At the genus
# level, relative-abundance values for the top 30 taxa were exported via
# Python v3.11 and summarized as mean relative abundance per sampling
# site using pandas, with per-specimen values retained in a supplementary
# table."
#
# This was a step described in Methods but not yet covered by any script
# in this pipeline until now (18/19 only handle Phylum/Family, top 10,
# per-specimen). Two outputs, matching the two things Methods says exist:
#   1. top30_genus_relative_abundance.csv (summarize_taxonomy.py, reused
#      from 18_export_taxa_top10.sh, just with --rank Genus --top-n 30)
#      -- per-specimen values, "retained in a supplementary table"
#   2. genus_mean_relative_abundance_by_site.csv (new: aggregate_genus_by_site.py)
#      -- mean relative abundance per site, the main reported genus table

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

echo "Step 1: per-specimen top-30 genus relative abundance..."
python3 "$SCRIPT_DIR/summarize_taxonomy.py" \
  --feature-table "$feature_table_tsv" \
  --taxonomy "$taxonomy_tsv" \
  --rank Genus \
  --top-n 30 \
  --outdir "$TAXA_SUMMARY_DIR"

echo ""
echo "Step 2: mean relative abundance per sampling site..."
python3 "$SCRIPT_DIR/aggregate_genus_by_site.py" \
  --top-table "$TAXA_SUMMARY_DIR/top30_genus_relative_abundance.csv" \
  --metadata "$METADATA" \
  --group-column "$GROUP_COLUMN" \
  --outdir "$TAXA_SUMMARY_DIR"

echo ""
echo "Done. Output: $TAXA_SUMMARY_DIR"
echo "  top30_genus_relative_abundance.csv          (per-specimen, supplementary table)"
echo "  genus_mean_relative_abundance_by_site.csv    (mean per site, main genus table)"
