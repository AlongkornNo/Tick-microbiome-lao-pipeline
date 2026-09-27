#!/bin/bash
# 00_run_fastqc_multiqc_rawdata.sh
#
# Quality-check the raw paired-end fastq.gz reads before any trimming,
# per Methods: "Initial read quality was assessed with FastQC v0.12.1 and
# summarized using MultiQC v1.28."
#
# Usage: source config.sh from the same pipeline directory, then run this
# script; RAW_DIR (set in config.sh) must contain the raw *.fastq.gz files.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

if ! compgen -G "$RAW_DIR"/*.fastq.gz > /dev/null; then
  echo "ERROR: no .fastq.gz files found in $RAW_DIR" >&2
  exit 1
fi

echo "FastQC version: $(fastqc --version)"
echo "MultiQC version: $(multiqc --version)"

FASTQC_OUT="$RAW_DIR/fastqc_results"
mkdir -p "$FASTQC_OUT"

echo "Running FastQC on raw reads in $RAW_DIR ..."
fastqc "$RAW_DIR"/*.fastq.gz -o "$FASTQC_OUT"

echo "Running MultiQC on FastQC output ..."
multiqc "$FASTQC_OUT" --outdir "$RAW_DIR" --filename multiqc_rawdata_qc --force

echo "Done. Report: $RAW_DIR/multiqc_rawdata_qc.html"
