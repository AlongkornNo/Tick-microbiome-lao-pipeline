#!/bin/bash
# 09_run_rarefy_normalization.sh
#
# Rarefy the bacteria-only feature table to a fixed sampling depth
# (RAREFACTION_DEPTH, 1800 reads/specimen by default -- config.sh), then
# export it to a QIIME2-format TSV (dada2_table_rm_organelles_rarefied.tsv)
# -- the exact file used for all downstream alpha-diversity and taxonomy
# verification in this reproducibility package (confirmed with the
# corresponding author: it is exported from this step's output,
# ${PROJECT}_table_rarefied.qza).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

min_depth="${1:-$RAREFACTION_DEPTH}"

echo "Starting rarefaction with depth: $min_depth"

input_table="$QIIME_DIR/${PROJECT}_table_bacteria_only.qza"
if [ ! -f "$input_table" ]; then
  echo "ERROR: input table not found: $input_table" >&2
  echo "Run 08_run_remove_organelle.sh first." >&2
  exit 1
fi
if [ ! -f "$METADATA" ]; then
  echo "ERROR: metadata file not found: $METADATA" >&2
  exit 1
fi

output_table="$QIIME_DIR/${PROJECT}_table_rarefied.qza"
output_viz="$QIIME_DIR/${PROJECT}_table_rarefied.qzv"

echo "qiime2 version: $(qiime --version 2>&1 | head -1)"

qiime feature-table rarefy \
  --i-table "$input_table" \
  --p-sampling-depth "$min_depth" \
  --o-rarefied-table "$output_table"

echo "Summarizing rarefied table..."
qiime feature-table summarize \
  --i-table "$output_table" \
  --m-sample-metadata-file "$METADATA" \
  --o-visualization "$output_viz"

echo "Exporting rarefied table to TSV (dada2_table_rm_organelles_rarefied.tsv)..."
export_dir="$DIVERSITY_DIR/rarefied_table_export"
mkdir -p "$export_dir"
qiime tools export --input-path "$output_table" --output-path "$export_dir"
biom convert \
  -i "$export_dir/feature-table.biom" \
  -o "$DIVERSITY_DIR/dada2_table_rm_organelles_rarefied.tsv" \
  --to-tsv

echo "Rarefaction and export complete."
echo "  Rarefied table (qza): $output_table"
echo "  Visualization: $output_viz"
echo "  Rarefied table (tsv): $DIVERSITY_DIR/dada2_table_rm_organelles_rarefied.tsv"
