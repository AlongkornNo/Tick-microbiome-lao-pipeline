#!/bin/bash
# 13_run_build_phylogeny.sh
#
# Build a rooted phylogenetic tree from the DADA2 representative
# sequences (MAFFT alignment -> masking -> FastTree -> midpoint
# rooting). Required before core-metrics-phylogenetic, since UniFrac
# and Faith's Phylogenetic Diversity are phylogeny-aware metrics -- per
# Methods: "The rarefied feature table, rooted phylogenetic tree, and
# sample metadata were loaded as QIIME2 artifacts."

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

rep_seqs="$QIIME_DIR/${PROJECT}_rep_seqs.qza"
if [ ! -f "$rep_seqs" ]; then
  echo "ERROR: representative sequences not found: $rep_seqs" >&2
  echo "Run 05_run_dada2_asv.sh first." >&2
  exit 1
fi

mkdir -p "$PHYLOGENY_DIR"

echo "qiime2 version: $(qiime --version 2>&1 | head -1)"
echo "Aligning representative sequences and building a rooted tree (MAFFT + FastTree)..."

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$rep_seqs" \
  --o-alignment "$PHYLOGENY_DIR/${PROJECT}_aligned_rep_seqs.qza" \
  --o-masked-alignment "$PHYLOGENY_DIR/${PROJECT}_masked_aligned_rep_seqs.qza" \
  --o-tree "$PHYLOGENY_DIR/${PROJECT}_unrooted_tree.qza" \
  --o-rooted-tree "$PHYLOGENY_DIR/${PROJECT}_rooted_tree.qza"

echo "Done. Rooted tree: $PHYLOGENY_DIR/${PROJECT}_rooted_tree.qza"
