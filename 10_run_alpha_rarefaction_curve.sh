#!/bin/bash
# 10_run_alpha_rarefaction_curve.sh
#
# Generate the alpha-rarefaction curve visualization, confirming that the
# chosen sampling depth (RAREFACTION_DEPTH, config.sh) adequately
# captures each specimen's diversity.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

input_table="$QIIME_DIR/${PROJECT}_table_rarefied.qza"
if [ ! -f "$input_table" ]; then
  echo "ERROR: rarefied table not found: $input_table" >&2
  echo "Run 09_run_rarefy_normalization.sh first." >&2
  exit 1
fi
if [ ! -f "$METADATA" ]; then
  echo "ERROR: metadata file not found: $METADATA" >&2
  exit 1
fi

output_viz="$QIIME_DIR/dada2_table_rm_organelles_rarefied_rarefaction.qzv"

echo "qiime2 version: $(qiime --version 2>&1 | head -1)"
echo "Calculating rarefaction curve (min-depth=10, max-depth=$RAREFACTION_DEPTH)..."

qiime diversity alpha-rarefaction \
  --i-table "$input_table" \
  --p-max-depth "$RAREFACTION_DEPTH" \
  --p-min-depth 10 \
  --m-metadata-file "$METADATA" \
  --o-visualization "$output_viz"

echo "Done. Visualization: $output_viz"
