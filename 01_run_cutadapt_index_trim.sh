#!/bin/bash
# 01_run_cutadapt_index_trim.sh
#
# Trim the per-sample sequencing index (barcode) from each raw paired-end
# read, using the sample -> index-pair mapping in samples/sample_index_map.tsv
# (one row per sample, columns: sample_id, fwd_index, rev_index) instead of
# hardcoding each sample's cutadapt call individually.
#
# NOTE ON COVERAGE: any sample listed in metadata.tsv but missing from
# sample_index_map.tsv is reported below and skipped, rather than
# silently dropped or given a guessed index.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

if [ ! -f "$SAMPLE_INDEX_MAP" ]; then
  echo "ERROR: sample index map not found: $SAMPLE_INDEX_MAP" >&2
  exit 1
fi

echo "cutadapt version: $(cutadapt --version)"

mkdir -p "$INDEX_TRIMMED_DIR"

tail -n +2 "$SAMPLE_INDEX_MAP" | while IFS=$'\t' read -r sample fwd_index rev_index; do
  [ -z "$sample" ] && continue
  r1="$RAW_DIR/${sample}.raw_1.fastq.gz"
  r2="$RAW_DIR/${sample}.raw_2.fastq.gz"

  if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "WARNING: skipping $sample -- raw fastq files not found ($r1 / $r2)" >&2
    continue
  fi

  echo "Processing $sample (fwd_index=$fwd_index, rev_index=$rev_index)..."
  cutadapt \
    -g "$fwd_index" -G "$rev_index" \
    -o "$INDEX_TRIMMED_DIR/${sample}_R1_trimmed.fastq.gz" \
    -p "$INDEX_TRIMMED_DIR/${sample}_R2_trimmed.fastq.gz" \
    "$r1" "$r2"
done

# Cross-check against metadata.tsv: flag any sample expected by the study
# design that never got an index-trim (e.g. because it has no row in
# sample_index_map.tsv yet).
if [ -f "$METADATA" ]; then
  echo
  echo "Checking sample_index_map.tsv coverage against $METADATA ..."
  comm -23 \
    <(tail -n +2 "$METADATA" | cut -f1 | grep -v '^#q2:types$' | sort -u) \
    <(tail -n +2 "$SAMPLE_INDEX_MAP" | cut -f1 | sort -u) \
    > "$INDEX_TRIMMED_DIR/samples_missing_index.txt" || true
  if [ -s "$INDEX_TRIMMED_DIR/samples_missing_index.txt" ]; then
    echo "WARNING: the following samples are in $METADATA but have NO row in" >&2
    echo "$SAMPLE_INDEX_MAP and were NOT index-trimmed:" >&2
    cat "$INDEX_TRIMMED_DIR/samples_missing_index.txt" >&2
  else
    echo "All samples in $METADATA have a matching index-map entry."
  fi
fi

echo
echo "Index trimming complete. Output: $INDEX_TRIMMED_DIR"
