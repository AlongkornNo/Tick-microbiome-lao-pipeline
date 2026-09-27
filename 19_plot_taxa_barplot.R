#!/usr/bin/env Rscript
# 19_plot_taxa_barplot.R
#
# Fig. 2A (Phylum) / Fig. 2B (Family): stacked bar charts of relative
# abundance, top 10 taxa + "Other" per specimen.
#
# Per Methods: "Data were reshaped in R v4.5.1 using the reshape2 v1.4.4
# package. Stacked plots were generated with ggplot2 v3.5.2 with colors
# assigned using the RColorBrewer v1.1.3 palettes."
#
# Input: the wide-format top-N(+Other) CSVs from 18_export_taxa_top10.sh
# (rows = taxa, columns = specimens, already true relative abundance --
# do NOT renormalize here).
#
# NOTE ON LAYOUT: Methods does not specify bar order/grouping explicitly.
# This script orders specimens by sampling site (metadata Group column,
# in the site order given in Table 1) since that is how every other
# figure/table in the manuscript is organized -- confirm this matches the
# actual published Fig. 2 before treating this as final.
#
# Usage:
#   Rscript 19_plot_taxa_barplot.R <taxa_summary_dir> <metadata_tsv> \
#     <group_column> <plots_dir> <top_n>

suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript 19_plot_taxa_barplot.R <taxa_summary_dir> ",
       "<metadata_tsv> <group_column> <plots_dir> <top_n>")
}
taxa_summary_dir <- args[1]
metadata_path    <- args[2]
group_column     <- args[3]
plots_dir        <- args[4]
top_n            <- as.integer(args[5])

dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read.delim(metadata_path, check.names = FALSE, row.names = 1)
if (nrow(metadata) > 0 && rownames(metadata)[1] == "#q2:types") {
  metadata <- metadata[-1, , drop = FALSE]
}
# Site order = first-appearance order in metadata.tsv (matches Table 1 order)
site_order <- unique(metadata[[group_column]])

plot_rank <- function(rank_label, csv_name, out_name) {
  csv_path <- file.path(taxa_summary_dir, csv_name)
  if (!file.exists(csv_path)) {
    stop("Not found: ", csv_path, " -- run 18_export_taxa_top10.sh first.")
  }

  wide <- read.csv(csv_path, check.names = FALSE, row.names = 1)
  # Sanity check: every specimen column should already sum to 1.0
  col_sums <- colSums(wide)
  bad <- abs(col_sums - 1.0) > 1e-6
  if (any(bad)) {
    warning(rank_label, ": ", sum(bad), " specimen(s) do not sum to 1.0 -- ",
            "input CSV may not be the true-relative-abundance export.")
  }

  wide$Taxon <- rownames(wide)
  # "Other" always last in the legend/stack, remaining taxa by descending
  # mean abundance (already the row order summarize_taxonomy.py wrote).
  taxon_order <- c(setdiff(wide$Taxon, "Other"), "Other")
  wide$Taxon <- factor(wide$Taxon, levels = rev(taxon_order))

  long <- melt(wide, id.vars = "Taxon", variable.name = "sample_id",
               value.name = "relative_abundance")

  long$Group <- metadata[as.character(long$sample_id), group_column]
  long$sample_id <- factor(long$sample_id,
                            levels = rownames(metadata)[order(match(metadata[[group_column]], site_order))])

  n_colors <- top_n + 1  # + "Other"
  palette <- if (n_colors <= 12) {
    colorRampPalette(brewer.pal(min(n_colors, 9), "Set3"))(n_colors)
  } else {
    colorRampPalette(brewer.pal(9, "Set3"))(n_colors)
  }
  names(palette) <- rev(taxon_order)
  palette["Other"] <- "grey70"

  p <- ggplot(long, aes(x = sample_id, y = relative_abundance, fill = Taxon)) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = palette, breaks = taxon_order) +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    labs(x = "Specimen", y = "Relative abundance",
         fill = rank_label, title = paste0(rank_label, "-level composition")) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          panel.spacing = unit(0.2, "lines"),
          strip.text = element_text(size = 8))

  out_path <- file.path(plots_dir, out_name)
  ggsave(out_path, p, width = 10, height = 5, dpi = 300)
  cat("Wrote:", out_path, "\n")
}

plot_rank("Phylum", paste0("top", top_n, "_phylum_relative_abundance.csv"), "Fig2A_phylum_barplot.png")
plot_rank("Family", paste0("top", top_n, "_family_relative_abundance.csv"), "Fig2B_family_barplot.png")

cat("\nDone.\n")
