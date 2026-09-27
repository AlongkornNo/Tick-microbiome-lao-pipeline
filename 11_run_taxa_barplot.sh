#!/bin/bash
# 11_run_taxa_barplot.sh
#
# Generate the interactive QIIME2 taxa barplot from the rarefied,
# organelle-filtered feature table + assigned taxonomy -- the basis for
# the manuscript's Figure 2. Uses metadata.tsv (has the Group column),
# not manifest.tsv, so the barplot can be interactively grouped/colored
# by sampling site.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

table="$QIIME_DIR/${PROJECT}_table_rarefied.qza"
taxonomy="$QIIME_DIR/silvaV3V4_taxonomy.qza"
output_viz="$QIIME_DIR/${PROJECT}_taxa_barplot.qzv"

if [[ ! -f "$table" || ! -f "$taxonomy" ]]; then
  echo "ERROR: input file(s) not found:" >&2
  echo "  $table" >&2
  echo "  $taxonomy" >&2
  echo "Run 09_run_rarefy_normalization.sh and 07_run_taxonomy_assignment.sh first." >&2
  exit 1
fi
if [ ! -f "$METADATA" ]; then
  echo "ERROR: metadata file not found: $METADATA" >&2
  exit 1
fi

echo "qiime2 version: $(qiime --version 2>&1 | head -1)"
echo "Generating taxa barplot..."

qiime taxa barplot \
  --i-table "$table" \
  --i-taxonomy "$taxonomy" \
  --m-metadata-file "$METADATA" \
  --o-visualization "$output_viz"

echo "Done. Visualization: $output_viz"
