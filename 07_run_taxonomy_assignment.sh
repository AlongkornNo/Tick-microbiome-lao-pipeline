#!/bin/bash
# 07_run_taxonomy_assignment.sh
#
# Assign taxonomy to the DADA2 representative sequences using the
# in-house-trained classifier (06_run_classifier_training.sh) -- per
# Methods: "Representative sequences were taxonomically classified using
# the QIIME 2 classify-sklearn method with a naive Bayes classifier
# trained on the SILVA 138 reference database..."
#
# Reads the classifier path from $CLASSIFIER_QZA (config.sh) rather than
# a hardcoded filename, so a filename typo here can't silently point at
# a nonexistent classifier.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

rep_seqs="$QIIME_DIR/${PROJECT}_rep_seqs.qza"

if [ ! -f "$CLASSIFIER_QZA" ]; then
  echo "ERROR: classifier not found: $CLASSIFIER_QZA" >&2
  echo "Run 06_run_classifier_training.sh first." >&2
  exit 1
fi
if [ ! -f "$rep_seqs" ]; then
  echo "ERROR: representative sequences not found: $rep_seqs" >&2
  echo "Run 05_run_dada2_asv.sh first." >&2
  exit 1
fi

echo "qiime2 version: $(qiime --version 2>&1 | head -1)"
echo "Classifier: $CLASSIFIER_QZA"
echo "Assigning taxonomy (classify-sklearn, confidence=$CLASSIFIER_CONFIDENCE)..."

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER_QZA" \
  --i-reads "$rep_seqs" \
  --p-confidence "$CLASSIFIER_CONFIDENCE" \
  --o-classification "$QIIME_DIR/silvaV3V4_taxonomy.qza"

echo "Generating taxonomy summary visualization..."
qiime metadata tabulate \
  --m-input-file "$QIIME_DIR/silvaV3V4_taxonomy.qza" \
  --o-visualization "$QIIME_DIR/silvaV3V4_taxonomy.qzv"

echo "Done. Output: $QIIME_DIR/silvaV3V4_taxonomy.qza(.qzv)"
