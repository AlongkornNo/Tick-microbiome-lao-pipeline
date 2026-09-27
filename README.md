# Cattle tick (*Rhipicephalus microplus*) 16S rRNA microbiome pipeline

Reproducibility package for the bioinformatics and statistical analysis
described in the Methods section of the manuscript (Lao PDR cattle tick
microbiome, 16S rRNA V3–V4 amplicon sequencing, 30 specimens across six
sampling sites). Scripts are numbered in pipeline order, `00` through
`25`, and all read their paths and parameters from a single
[`config.sh`](config.sh).

## Prerequisites

| Tool | Version used | Needed for |
|---|---|---|
| QIIME 2 | 2024.10 (amplicon distribution) | scripts `00`–`18`, `20`, `23`–`25` |
| Cutadapt | 5.5 | `01`, `02` |
| FastQC | 0.12.1 | `00`, `03` |
| MultiQC | 1.28 | `00`, `03` |
| Python | 3.11, with `pandas`, `numpy`, `scipy` | `12`, `18`, `20`, `25` and their helper modules |
| R | 4.5.1, with `vegan`, `ape`, `ggplot2`, `reshape2`, `RColorBrewer` | `17`, `19`, `21`, `22` |

Install QIIME 2 per the [official instructions](https://docs.qiime2.org/)
for your platform; everything else is a normal `pip`/`conda install` or
`install.packages()`.

## Setup

1. Copy `config.sh` next to your own `scripts/` checkout (or edit it in
   place) and set `BASE_DIR` to wherever your data lives:
   ```bash
   export BASE_DIR="/path/to/your/microbiome/project"
   ```
2. Populate `$BASE_DIR` with your raw paired-end FASTQ files, the
   sample/index map, and `metadata.tsv` (see the header comments in
   `config.sh` for the exact expected layout of each).
3. Run the scripts in numeric order from the `scripts/` directory, e.g.:
   ```bash
   bash scripts/00_run_fastqc_multiqc_rawdata.sh
   bash scripts/01_run_cutadapt_index_trim.sh
   # ...
   ```
   Each script checks that its required inputs exist and prints a clear
   error naming the earlier step to run first if they don't.

## Script index

| Script | What it does |
|---|---|
| `00_run_fastqc_multiqc_rawdata.sh` | FastQC + MultiQC on raw reads |
| `01_run_cutadapt_index_trim.sh` | Remove per-sample sequencing indices |
| `02_trim_primers_per_sample.sh` | Remove 341F/806R primers, discard reads <200 bp |
| `03_run_fastqc_multiqc_trimmed.sh` | FastQC + MultiQC on trimmed reads |
| `04_import_fastq_data.sh` | Import into QIIME 2 (manifest format) |
| `05_run_dada2_asv.sh` | DADA2 denoising, merging, chimera removal |
| `06_run_classifier_training.sh` | Train the in-house SILVA 138 V3–V4 naive Bayes classifier |
| `07_run_taxonomy_assignment.sh` | `classify-sklearn` taxonomy assignment |
| `08_run_remove_organelle.sh` | Filter Eukaryota/Archaea/Unassigned/mitochondria/chloroplast |
| `09_run_rarefy_normalization.sh` | Rarefy to 1,800 reads/sample; export TSV |
| `10_run_alpha_rarefaction_curve.sh` | Rarefaction curves (Supplementary Fig. S1) |
| `11_run_taxa_barplot.sh` | QIIME 2 taxa-barplot visualization |
| `12_export_barplot_csv.sh` | Export feature table + taxonomy, summarize by rank |
| `13_run_build_phylogeny.sh` | MAFFT + FastTree rooted phylogeny |
| `14_run_core_diversity_metrics.sh` | `core-metrics-phylogenetic` (Shannon, Observed Features, Evenness, Faith's PD; also unrelated beta byproduct — see Known limitations) |
| `15_run_additional_alpha_metrics.sh` | Chao1, Simpson, Good's coverage |
| `16_run_beta_diversity_stats.sh` | Beta-diversity distance matrices + **global** PERMANOVA (Table 3) |
| `17_run_permdisp.R` | PERMDISP (`betadisper`/`permutest`, Table 3) |
| `18_export_taxa_top10.sh` | Phylum/Family top-10(+Other) relative abundance |
| `19_plot_taxa_barplot.R` | Fig. 2A/2B stacked bar charts |
| `20_run_alpha_diversity_stats.sh` | Chao1/Shannon/Simpson/Evenness + Kruskal-Wallis + Dunn's + BH-FDR (Python) |
| `21_plot_alpha_boxplots.R` | Fig. 3A–D boxplots |
| `22_plot_beta_pcoa.R` | Fig. 4A–D PCoA ordinations, PERMANOVA + PERMDISP annotated |
| `23_run_pairwise_permanova.sh` | Pairwise site-vs-site PERMANOVA (Tables S9–S12) |
| `24_run_alpha_group_significance.sh` | QIIME2-native alpha-diversity significance testing (cross-check for `20`) |
| `25_export_genus_top30.sh` | Genus-level top-30, mean relative abundance per site |
| `summarize_taxonomy.py` | Shared helper: true (non-renormalized) relative abundance by rank |
| `alpha_diversity_stats.py` | Shared helper: from-scratch Chao1/Shannon/Simpson/Pielou + stats |
| `aggregate_genus_by_site.py` | Shared helper: per-specimen → per-site genus means |

## Known limitations (read before treating any number here as final)

This package was reconstructed and cross-checked against the original
analysis scripts and the manuscript's already-published results, not
run start-to-finish end-to-end in one pass (this environment has no
QIIME 2, R, or raw FASTQ data available to do that). Please be aware of
the following before reusing or citing results from it:

- **Three independent rarefactions, not one.** `qiime feature-table
  rarefy` (step `09`) and `core-metrics-phylogenetic`'s internal
  resampling (step `14`) each draw an independent random subsample at
  the same depth (1,800 reads/sample) — neither has a fixed seed. As a
  result: Shannon/Observed Features/Evenness/Faith's PD, and the
  pairwise PERMANOVA (`23`) + PERMDISP (`17`), all come from step 14's
  draw; Chao1/Simpson/Good's coverage, and the **global** PERMANOVA in
  Table 3 (`16`), come from step 09's draw instead. A from-scratch
  recomputation using one draw will not exactly reproduce results
  computed from the other, even though both are at the same depth from
  the same input table. See the header comments in `14`, `15`, `16`,
  `17`, and `23` for the full detail.
- **R scripts are unexecuted.** `17_run_permdisp.R`, `19_plot_taxa_barplot.R`,
  `21_plot_alpha_boxplots.R`, and `22_plot_beta_pcoa.R` were written
  against the Methods text and, where possible, real analysis scripts
  and verified helper logic — but never run in an R environment (none
  was available while writing them). Test them on your own data before
  trusting their output.
- **Genus-level site means (`25`) do not exactly match a previously
  computed reference table** for this dataset (differences of several
  percentage points on some dominant genera, e.g. ~55% vs. ~60% for one
  genus at one site). The per-specimen relative-abundance logic
  (`summarize_taxonomy.py`) is the same code independently verified
  against the manuscript's phylum/family table (exact match), so the
  likely explanation is the same rarefaction-snapshot issue above
  rather than a logic error — but this was not fully traced to ground
  truth.
- **PERMDISP uses `betadisper`'s default (`type = "median"`)**, not
  `"centroid"` — confirmed from the real analysis script, but worth
  double-checking against your own intent if you adapt this elsewhere.

## License

Released under the [MIT License](LICENSE) (see file) — a permissive
default chosen for a scripts-only reproducibility package. Change it if
your institution or co-authors require something else (e.g. a
data/code combination is sometimes released under MIT for code and
CC-BY-4.0 for any accompanying data/figures).

## Citation

Once archived (e.g. via a Zenodo-linked GitHub release), add the
resulting DOI here and in the manuscript's data/code availability
statement.
