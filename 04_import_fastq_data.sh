#!/bin/bash
# 04_import_fastq_data.sh
#
# Import the primer-trimmed paired-end fastq files into QIIME2 using a
# manifest of absolute file paths -- per Methods: "Paired-end reads were
# imported into QIIME2 (qiime2-amplicon-2024.10) using the
# PairedEndFastqManifestPhred33V2 format."
#
# Requires metadata/manifest.tsv (see the manifest-generation step, which
# builds it from sample_index_map.tsv + the PRIMER_TRIMMED_DIR fastq
# files -- run that first if $MANIFEST does not exist yet).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

if [ ! -f "$MANIFEST" ]; then
  echo "ERROR: manifest file not found: $MANIFEST" >&2
  echo "Generate it first (sample-id,forward-absolute-filepath,reverse-absolute-filepath)." >&2
  exit 1
fi

mkdir -p "$QIIME_DIR"

echo "qiime2 version: $(qiime --version 2>&1 | head -1)"
echo "Importing paired-end fastq data from $MANIFEST ..."

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path "$MANIFEST" \
  --output-path "$QIIME_DIR/${PROJECT}_demux_seqs.qza"

qiime demux summarize \
  --i-data "$QIIME_DIR/${PROJECT}_demux_seqs.qza" \
  --o-visualization "$QIIME_DIR/${PROJECT}_demux_seqs.qzv"

echo "Import complete. Output: $QIIME_DIR/${PROJECT}_demux_seqs.qza(.qzv)"
