#!/bin/bash
# 05_run_dada2_asv.sh
#
# Denoise paired-end reads, merge, remove chimeras, and build the ASV
# feature table with DADA2, per Methods: "DADA2, as implemented in QIIME
# 2 2024.10, was used for denoising, paired-end read merging, chimera
# removal, and feature-table construction. DADA2 was run with
# trim-left-f = 0, trim-left-r = 0, trunc-len-f = 0, trunc-len-r = 0,
# max-ee-f = 2, max-ee-r = 2, and trunc-q = 2."
#
# Uses metadata.tsv (has the Group column) -- NOT manifest.tsv (only
# has file paths, no per-sample metadata) -- for the feature-table
# summary visualization.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

input_seqs="$QIIME_DIR/${PROJECT}_demux_seqs.qza"
if [ ! -f "$input_seqs" ]; then
  echo "ERROR: demultiplexed sequences not found: $input_seqs" >&2
  echo "Run 04_import_fastq_data.sh first." >&2
  exit 1
fi
if [ ! -f "$METADATA" ]; then
  echo "ERROR: metadata file not found: $METADATA" >&2
  exit 1
fi

echo "qiime2 version: $(qiime --version 2>&1 | head -1)"
echo "Running DADA2 denoising (trim-left-f/r=0, trunc-len-f/r=0, max-ee-f/r=2, trunc-q=2)..."

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$input_seqs" \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-trunc-q 2 \
  --p-n-threads 0 \
  --o-table "$QIIME_DIR/${PROJECT}_table.qza" \
  --o-representative-sequences "$QIIME_DIR/${PROJECT}_rep_seqs.qza" \
  --o-denoising-stats "$QIIME_DIR/${PROJECT}_stats.qza" \
  --verbose

echo "Summarizing feature table (against $METADATA)..."
qiime feature-table summarize \
  --i-table "$QIIME_DIR/${PROJECT}_table.qza" \
  --m-sample-metadata-file "$METADATA" \
  --o-visualization "$QIIME_DIR/${PROJECT}_table.qzv"

qiime feature-table tabulate-seqs \
  --i-data "$QIIME_DIR/${PROJECT}_rep_seqs.qza" \
  --o-visualization "$QIIME_DIR/${PROJECT}_rep_seqs.qzv"

qiime metadata tabulate \
  --m-input-file "$QIIME_DIR/${PROJECT}_stats.qza" \
  --o-visualization "$QIIME_DIR/${PROJECT}_stats.qzv"

echo "DADA2 step complete. Output: $QIIME_DIR/${PROJECT}_table.qza, ${PROJECT}_rep_seqs.qza, ${PROJECT}_stats.qza"
