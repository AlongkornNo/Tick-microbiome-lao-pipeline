"""
summarize_taxonomy.py — per-specimen relative abundance by taxonomic rank,
from a rarefied ASV feature table and a QIIME2 taxonomy assignment table.

WHY THIS SCRIPT EXISTS / WHAT IT FIXES
---------------------------------------
An earlier version of Table S3 / Figure 2 (top 10 phyla, top 10 families)
was built by summing counts across ALL ASVs at a rank, but then dividing
each specimen's top-10 taxa by the SUM OF THE TOP 10 ONLY, instead of by
the specimen's true total read count. That silently renormalized the top
10 taxa to 100% and hid the fact that, especially at family level, a large
and variable fraction of each specimen's reads belong to taxa outside the
top 10 (up to ~98% for some specimens).

This script avoids that bug structurally: for every rank, it first sums
counts across ALL taxa observed at that rank (every ASV, not just the top
N) to get each specimen's TRUE total, and only then computes relative
abundance = (reads assigned to taxon X) / (specimen's true total reads).
The top-N table it writes always includes an explicit "Other" row equal to
1 - sum(top N), so the columns for every specimen still sum to exactly 1.0.

INPUT
-----
--feature-table : rarefied ASV feature table, QIIME2 TSV export
                  (comment line "# Constructed from biom file", then a
                  header "#OTU ID<TAB>sample1<TAB>...", one row per ASV).
--taxonomy      : QIIME2 taxonomy.tsv (from `qiime feature-classifier
                  classify-sklearn` / exported taxonomy.qza), columns
                  Feature ID, Taxon, Confidence. A "#q2:types" metadata
                  row, if present, is dropped automatically.
--rank          : which rank(s) to summarize (default: Phylum,Family --
                  the two ranks reported in Fig. 2 / Table S3). Any of
                  Domain,Phylum,Class,Order,Family,Genus,Species.
--top-n         : how many top taxa to keep explicitly before folding the
                  remainder into "Other" (default 10, matching
                  TOP_N_BARPLOT in the pipeline config).

OUTPUT (per rank, in --outdir)
-------------------------------
relative_abundance_<rank>.csv        long format, ALL taxa, every specimen
                                      (no top-N cutoff -- the full, true
                                      breakdown; useful for auditing "Other").
top<N>_<rank>_relative_abundance.csv wide format: rows = top N taxa (by
                                      mean relative abundance across
                                      specimens) + an explicit "Other" row;
                                      columns = specimens; every column
                                      sums to 1.0.

VERIFICATION
------------
Run on the real data (dada2_table_rm_organelles_rarefied.tsv, 1800
reads/specimen, and the SILVA V3-V4 classifier's taxonomy.tsv export),
this script's top10_Phylum and top10_Family tables were compared
line-by-line against the corrected Table S3
(Table_S3_top10_Phylum-Family_relative_abundance_with_Other.xlsx, produced
independently by compute_relative_abundance_with_other.py from Table S2)
and matched exactly for all 30 specimens at both ranks (see
verify_summarize_taxonomy.py). This confirms the two independent
computation paths -- one from the merged Table S2, one from the raw
feature-table + taxonomy.tsv pair, as the actual pipeline produces them --
agree.

USAGE
-----
    python3 summarize_taxonomy.py \\
        --feature-table dada2_table_rm_organelles_rarefied.tsv \\
        --taxonomy taxonomy.tsv \\
        --rank Phylum,Family --top-n 10 \\
        --outdir ./out
"""

import argparse
from pathlib import Path

import pandas as pd

RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def load_feature_table(path):
    with open(path) as f:
        first_line = f.readline()
    skip = 1 if first_line.startswith("#") and "Constructed from biom" in first_line else 0
    df = pd.read_csv(path, sep="\t", skiprows=skip)
    df = df.rename(columns={df.columns[0]: "Feature ID"})
    return df.set_index("Feature ID")


def load_taxonomy(path):
    tax = pd.read_csv(path, sep="\t")
    tax = tax.rename(columns={tax.columns[0]: "Feature ID"})
    tax = tax[tax["Feature ID"] != "#q2:types"]  # drop QIIME2 metadata-type row, if present
    tax = tax.set_index("Feature ID")

    def split_taxon(taxon):
        parts = [p.strip() for p in str(taxon).split(";")]
        cleaned = [p.split("__", 1)[-1] if "__" in p else p for p in parts]
        cleaned = [c if c else "Unassigned" for c in cleaned]
        cleaned = (cleaned + ["Unassigned"] * len(RANKS))[: len(RANKS)]
        return cleaned

    split_cols = tax["Taxon"].apply(split_taxon).apply(pd.Series)
    split_cols.columns = RANKS
    return split_cols


def true_relative_abundance(counts_df, taxa_df, rank):
    """Sum counts across ALL taxa at `rank` (not just a top-N subset), then
    divide by each specimen's TRUE total (sum across every taxon at that
    rank, which equals the specimen's total read count). This is what
    prevents the renormalization bug described in the module docstring."""
    merged = counts_df.join(taxa_df[[rank]], how="left")
    merged[rank] = merged[rank].fillna("Unassigned")
    grouped = merged.groupby(rank)[counts_df.columns].sum()  # taxon x specimen, TRUE counts
    totals = grouped.sum(axis=0)  # specimen's true total reads (sum across ALL taxa)
    rel_abund = grouped.div(totals, axis=1)
    return rel_abund  # every column already sums to 1.0 (all taxa included)


def top_n_with_other(rel_abund, n):
    top_taxa = rel_abund.mean(axis=1).sort_values(ascending=False).head(n).index
    top_table = rel_abund.loc[top_taxa].copy()
    other = 1.0 - top_table.sum(axis=0)
    top_table.loc["Other"] = other
    return top_table


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--feature-table", required=True, help="Rarefied ASV feature table (QIIME2 TSV export)")
    ap.add_argument("--taxonomy", required=True, help="QIIME2 taxonomy.tsv (Feature ID, Taxon, Confidence)")
    ap.add_argument("--rank", default="Phylum,Family", help="Comma-separated ranks to summarize (default: Phylum,Family)")
    ap.add_argument("--top-n", type=int, default=10, help="Number of top taxa to keep before folding into Other (default 10)")
    ap.add_argument("--outdir", default=".", help="Output directory")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    counts_df = load_feature_table(args.feature_table)
    taxa_df = load_taxonomy(args.taxonomy)

    ranks = [r.strip() for r in args.rank.split(",")]
    for rank in ranks:
        if rank not in RANKS:
            raise SystemExit(f"Unknown rank '{rank}'; choose from {RANKS}")

        rel_abund = true_relative_abundance(counts_df, taxa_df, rank)

        # sanity check: every specimen's full (all-taxa) breakdown must sum to 1.0
        col_sums = rel_abund.sum(axis=0)
        bad = col_sums[(col_sums - 1.0).abs() > 1e-9]
        if len(bad):
            print(f"  WARNING ({rank}): {len(bad)} specimens do not sum to 1.0: {list(bad.index)}")

        rank_lower = rank.lower()
        long_path = outdir / f"relative_abundance_{rank_lower}.csv"
        rel_abund.reset_index().melt(
            id_vars=rank, var_name="sample-id", value_name="relative_abundance"
        ).to_csv(long_path, index=False)

        top_table = top_n_with_other(rel_abund, args.top_n)
        top_path = outdir / f"top{args.top_n}_{rank_lower}_relative_abundance.csv"
        top_table.to_csv(top_path)

        print(f"{rank}: {len(rel_abund)} taxa observed across {len(rel_abund.columns)} specimens")
        print(f"  wrote {long_path} (all taxa, long format)")
        print(f"  wrote {top_path} (top {args.top_n} + Other, wide format)")

    print("\nDone.")


if __name__ == "__main__":
    main()
