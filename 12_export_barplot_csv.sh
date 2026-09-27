#!/bin/bash
# 12_export_barplot_csv.sh
#
# Export the rarefied feature table + taxonomy, merge them, and
# summarize raw read counts by taxonomic rank (Phylum..Species).
#
# BUG FIX vs. the original version of this script: taxonomy strings
# split on "; " have index 0 = Domain (d__), not Phylum -- e.g.
# "d__Bacteria; p__Proteobacteria; c__...; o__...; f__...; g__...; s__..."
# The original script's `levels = ["Phylum", "Class", "Order", "Family",
# "Genus", "Species"]` list was paired with positional indices 0-5,
# silently extracting each rank ONE LEVEL TOO SHALLOW (the file named
# "phylum" actually contained Domain data, "family" actually contained
# Order data, etc.), and true Species-level data (index 6) was never
# extracted at all. This version matches ranks by their "x__" prefix
# instead of by position, which is also robust to taxonomy strings that
# are missing trailing ranks (no genus/species assigned).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../config.sh"

feature_table="$QIIME_DIR/${PROJECT}_table_rarefied.qza"
taxonomy="$QIIME_DIR/silvaV3V4_taxonomy.qza"
export_dir="$QIIME_DIR/barplot_export"
output_csv="$export_dir/taxa_table_with_names.csv"

if [ ! -f "$feature_table" ]; then
  echo "ERROR: feature table not found: $feature_table" >&2
  exit 1
fi
if [ ! -f "$taxonomy" ]; then
  echo "ERROR: taxonomy file not found: $taxonomy" >&2
  exit 1
fi

mkdir -p "$export_dir"

echo "1. Exporting feature table (BIOM)..."
qiime tools export --input-path "$feature_table" --output-path "$export_dir"

echo "2. Converting BIOM -> TSV..."
biom convert -i "$export_dir/feature-table.biom" -o "$export_dir/feature-table.tsv" --to-tsv

echo "3. Exporting taxonomy..."
qiime tools export --input-path "$taxonomy" --output-path "$export_dir/taxonomy_export"

echo "4. Merging + summarizing by taxonomic rank (Python)..."
python3 << EOF
import re
import pandas as pd
from pathlib import Path

export_dir = Path("$export_dir")
ft_path = export_dir / "feature-table.tsv"
tax_path = export_dir / "taxonomy_export" / "taxonomy.tsv"
out_csv = Path("$output_csv")

# --- Load feature table ---
ft = pd.read_csv(ft_path, sep='\t', skiprows=1)
ft.rename(columns={ft.columns[0]: "Feature ID"}, inplace=True)
ft.set_index("Feature ID", inplace=True)

# --- Load taxonomy ---
tax = pd.read_csv(tax_path, sep='\t')
tax = tax[tax["Feature ID"] != "#q2:types"]
tax.set_index("Feature ID", inplace=True)

# --- Merge ---
merged = ft.join(tax["Taxon"], how="left")
cols = ["Taxon"] + [c for c in merged.columns if c != "Taxon"]
merged = merged[cols]
merged.to_csv(out_csv)
print(f"Saved merged feature+taxonomy table -> {out_csv}")

# ------------------------------------------------
# Summarize by taxonomic rank (Phylum -> Species), matched by "x__"
# prefix rather than fixed position (see bug-fix note above).
# ------------------------------------------------
RANK_PREFIX = {
    "Phylum": "p__", "Class": "c__", "Order": "o__",
    "Family": "f__", "Genus": "g__", "Species": "s__",
}

def extract_rank(taxon_str, prefix):
    parts = [p.strip() for p in str(taxon_str).split(";")]
    for part in parts:
        if part.startswith(prefix):
            name = part.split("__", 1)[-1].strip()
            return name if name else "Unassigned"
    return "Unassigned"

for lvl, prefix in RANK_PREFIX.items():
    df = merged.copy()
    df[lvl] = df["Taxon"].apply(lambda x: extract_rank(x, prefix))
    summarized = df.groupby(lvl).sum(numeric_only=True)
    out_path = export_dir / f"taxa_barplot_{lvl.lower()}.csv"
    summarized.to_csv(out_path)
    print(f"Summarized by {lvl}: {out_path}")
EOF

echo "All exports complete."
echo "Output directory: $export_dir"
