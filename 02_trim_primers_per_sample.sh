#!/bin/bash
# 02_trim_primers_per_sample.sh
#
# Remove the 16S V3-V4 PCR primers (341F/806R) from each index-trimmed
# read pair (output of 01_run_cutadapt_index_trim.sh), trimming both the
# expected 5' primer and any 3' read-through of the opposite primer, then
# discard reads shorter than 200 bp -- per Methods: "Low-quality bases,
# adapters, and primer sequences were removed using Cutadapt v5.5 ... and
# reads shorter than 200 bps were discarded."
#
# Primers (from config.sh -- must match everywhere they're used):
#   341F = CCTAYGGGRBGCASCAG
#   806R = GGACTACNNGGGTATCTAAT

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

if [ ! -f "$SAMPLE_INDEX_MAP" ]; then
  echo "ERROR: sample index map not found: $SAMPLE_INDEX_MAP" >&2
  exit 1
fi

echo "cutadapt version: $(cutadapt --version)"
echo "PRIMER_F=$PRIMER_F"
echo "PRIMER_R=$PRIMER_R"

mkdir -p "$PRIMER_TRIMMED_DIR"

tail -n +2 "$SAMPLE_INDEX_MAP" | cut -f1 | while read -r sample; do
  [ -z "$sample" ] && continue
  r1_file="$INDEX_TRIMMED_DIR/${sample}_R1_trimmed.fastq.gz"
  r2_file="$INDEX_TRIMMED_DIR/${sample}_R2_trimmed.fastq.gz"

  if [[ ! -f "$r1_file" || ! -f "$r2_file" ]]; then
    echo "WARNING: skipping $sample -- missing index-trimmed R1/R2 ($r1_file / $r2_file)" >&2
    continue
  fi

  output_r1="$PRIMER_TRIMMED_DIR/${sample}_R1_primertrimmed.fastq.gz"
  output_r2="$PRIMER_TRIMMED_DIR/${sample}_R2_primertrimmed.fastq.gz"

  echo "Trimming primers for sample $sample..."
  cutadapt \
    -g "$PRIMER_F" -G "$PRIMER_R" \
    -a "$PRIMER_R" -A "$PRIMER_F" \
    --minimum-length=200 \
    -o "$output_r1" -p "$output_r2" \
    "$r1_file" "$r2_file"
done

echo
echo "Primer trimming complete. Output: $PRIMER_TRIMMED_DIR"
