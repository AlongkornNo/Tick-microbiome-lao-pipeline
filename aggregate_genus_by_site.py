"""
aggregate_genus_by_site.py — mean relative abundance per sampling site,
from the per-specimen top-N genus table.

Per Methods: "At the genus level, relative-abundance values for the top
30 taxa were exported via Python v3.11 and summarized as mean relative
abundance per sampling site using pandas, with per-specimen values
retained in a supplementary table."

INPUT
-----
--top-table : top30_genus_relative_abundance.csv from summarize_taxonomy.py
              (wide format: rows = top 30 genera + "Other", columns =
              specimens). This file IS the "per-specimen values retained
              in a supplementary table" -- this script does not modify it,
              only reads it.
--metadata  : TSV with columns sample-id, Group (site name).

OUTPUT
------
genus_mean_relative_abundance_by_site.csv : rows = top 30 genera + "Other",
columns = sampling sites, values = mean relative abundance across
specimens at that site.
"""

import argparse
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top-table", required=True, help="top30_genus_relative_abundance.csv (wide: taxa x specimen)")
    ap.add_argument("--metadata", required=True, help="TSV with sample-id, Group columns")
    ap.add_argument("--group-column", default="Group")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    top_table = pd.read_csv(args.top_table, index_col=0)  # taxon x specimen

    meta = pd.read_csv(args.metadata, sep="\t")
    meta = meta[meta.iloc[:, 0] != "#q2:types"]
    sample_col = meta.columns[0]
    site_of = dict(zip(meta[sample_col], meta[args.group_column]))

    missing = [s for s in top_table.columns if s not in site_of]
    if missing:
        print(f"WARNING: {len(missing)} specimen(s) in the top table have no metadata "
              f"entry and will be dropped from the site means: {missing[:5]}...")
    specimens = [s for s in top_table.columns if s in site_of]

    site_series = pd.Series({s: site_of[s] for s in specimens})
    # preserve first-appearance order in metadata (matches Table 1 order)
    site_order = list(dict.fromkeys(meta[args.group_column]))

    # groupby-mean across specimens per site (transpose so specimens are
    # rows -- avoids relying on groupby's deprecated axis=1 on columns)
    means = top_table[specimens].transpose().groupby(site_series).mean().transpose()
    means = means[[s for s in site_order if s in means.columns]]

    out_path = outdir / "genus_mean_relative_abundance_by_site.csv"
    means.to_csv(out_path)
    print(f"Wrote {out_path} ({len(means)} taxa x {len(means.columns)} sites)")
    print("(Per-specimen values are the input --top-table itself -- unchanged, "
          "serves as the supplementary table per Methods.)")


if __name__ == "__main__":
    main()
