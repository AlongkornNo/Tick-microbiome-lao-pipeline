#!/bin/bash
# config.sh — central configuration for the 16S microbiome pipeline.
#
# Every script under scripts/ sources this file, so this is the ONLY place
# you should need to edit paths/parameters for your own environment.
#
# Usage inside a script:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#   source "$SCRIPT_DIR/../config.sh"

# --- Project name (used as a filename prefix throughout) ---
export PROJECT="demo"

# --- Base directory: EDIT THIS to point at your own data location ---
export BASE_DIR="${BASE_DIR:-$HOME/microbiome}"

# --- Pipeline stage directories (all derived from BASE_DIR) ---
export RAW_DIR="$BASE_DIR/00_rawdata"
export INDEX_TRIMMED_DIR="$BASE_DIR/01_cutadapt_index_trimmed"
export PRIMER_TRIMMED_DIR="$BASE_DIR/02_primer_trimmed"
export QC_TRIMMED_DIR="$BASE_DIR/03_qc_trimmed"
export QIIME_DIR="$BASE_DIR/04_qiime2"
export CLASSIFIER_DIR="$BASE_DIR/classifier"

# --- Metadata / manifest files ---
# manifest.tsv  = sample-id,forward-absolute-filepath,reverse-absolute-filepath
#                 (fed to `qiime tools import`)
# metadata.tsv  = sample-id, Group, ...   (fed to --m-metadata-file everywhere else)
export METADATA_DIR="$BASE_DIR/metadata"
export MANIFEST="$METADATA_DIR/manifest.tsv"
export METADATA="$METADATA_DIR/metadata.tsv"

# --- Per-sample index map used by 01_run_cutadapt_index_trim.sh ---
# columns: sample_id  fwd_index  rev_index
export SAMPLE_INDEX_MAP="$BASE_DIR/samples/sample_index_map.tsv"

# --- Primers: V3-V4 region, 341F/806R (must match everywhere they're used) ---
export PRIMER_F="CCTAYGGGRBGCASCAG"
export PRIMER_R="GGACTACNNGGGTATCTAAT"

# --- Sampling depth for rarefaction (run 09 first to help you pick this) ---
export RAREFACTION_DEPTH="${RAREFACTION_DEPTH:-1800}"

# --- Taxonomy classification confidence (QIIME2 classify-sklearn) ---
export CLASSIFIER_CONFIDENCE="${CLASSIFIER_CONFIDENCE:-0.7}"

# --- In-house SILVA 138 V3-V4 classifier ---
# The classifier used for taxonomy assignment was trained IN-HOUSE, not
# downloaded pretrained. The workflow is:
#   1. qiime feature-classifier extract-reads on the raw SILVA 138 99%
#      reference (CLASSIFIER_RAW_SEQS/CLASSIFIER_RAW_TAX below), using
#      PRIMER_F/PRIMER_R (the SAME 341F/806R primers used to trim the
#      actual sequenced reads elsewhere in this pipeline -- there is no
#      separate "classifier extraction primer" pair; the two are identical).
#   2. qiime feature-classifier fit-classifier-naive-bayes on the
#      extracted reference reads + their taxonomy -> CLASSIFIER_QZA.
# See scripts/06_run_classifier_training.sh.
#
# Raw SILVA 138 99% OTU reference files (input to extract-reads):
export CLASSIFIER_RAW_SEQS="$CLASSIFIER_DIR/silva-138-99-seqs.qza"
export CLASSIFIER_RAW_TAX="$CLASSIFIER_DIR/silva-138-99-tax.qza"
# Output of the in-house training (extract-reads + fit-classifier-naive-bayes):
export CLASSIFIER_QZA="$CLASSIFIER_DIR/silva-138-99-seqs-341F-806R_classifier.qza"

# --- Downstream analysis directories ---
export PHYLOGENY_DIR="$BASE_DIR/05_phylogeny"
export CORE_METRICS_DIR="$QIIME_DIR/${PROJECT}_core_metrics"
export DIVERSITY_DIR="$BASE_DIR/06_diversity_tables"
export TAXA_SUMMARY_DIR="$BASE_DIR/07_taxonomy_summary"
export STATS_DIR="$BASE_DIR/08_statistics"
export PLOTS_DIR="$BASE_DIR/09_plots"

# --- Grouping column in metadata.tsv used for stats/plots/PERMANOVA ---
export GROUP_COLUMN="${GROUP_COLUMN:-Group}"

# --- How many top taxa to summarize / plot ---
export TOP_N_BARPLOT="${TOP_N_BARPLOT:-10}"
# TOP_N_HEATMAP is unused by any script in this repo (the taxonomy heatmap
# was dropped from the manuscript); kept only in case a heatmap is added
# back later.
export TOP_N_HEATMAP="${TOP_N_HEATMAP:-30}"

# --- Beta-diversity distance matrices to test/plot (must match the file
#     stems produced by core-metrics-phylogenetic) ---
export BETA_METRICS=(unweighted_unifrac weighted_unifrac jaccard bray_curtis)

# --- Alpha-diversity metrics assembled into one summary table ---
export ALPHA_METRICS=(shannon observed_features evenness faith_pd chao1 simpson goods_coverage)
