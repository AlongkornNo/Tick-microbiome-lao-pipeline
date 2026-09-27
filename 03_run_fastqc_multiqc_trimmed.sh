#!/bin/bash
# 03_run_fastqc_multiqc_trimmed.sh
#
# Quality-check the primer-trimmed reads (output of
# 02_trim_primers_per_sample.sh) to verify successful adapter and primer
# removal -- per Methods: "Trimmed reads were then re-evaluated to verify
# successful adapter and primer removal."
#
# Reads from PRIMER_TRIMMED_DIR (config.sh) -- the actual output
# directory of step 02 -- not a separately hardcoded path, which is what
# caused the original script to silently find zero files.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

if ! compgen -G "$PRIMER_TRIMMED_DIR"/*.fastq.gz > /dev/null; then
  echo "ERROR: no .fastq.gz files found in $PRIMER_TRIMMED_DIR" >&2
  echo "(did 02_trim_primers_per_sample.sh run first?)" >&2
  exit 1
fi

echo "FastQC version: $(fastqc --version)"
echo "MultiQC version: $(multiqc --version)"

FASTQC_DIR="$QC_TRIMMED_DIR/fastqc_reports"
mkdir -p "$FASTQC_DIR"

echo "Running FastQC on primer-trimmed reads in $PRIMER_TRIMMED_DIR ..."
fastqc "$PRIMER_TRIMMED_DIR"/*.fastq.gz -o "$FASTQC_DIR"

echo "Running MultiQC ..."
multiqc "$FASTQC_DIR" --outdir "$QC_TRIMMED_DIR" --filename multiqc_trimmed_qc --force

echo "Done. Report: $QC_TRIMMED_DIR/multiqc_trimmed_qc.html"
