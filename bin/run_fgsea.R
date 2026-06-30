#!/usr/bin/env Rscript
# Driver script for modules/gsea_fgsea.nf.
#
# Runs ranked-list GSEA (fgsea) against a user-supplied MSigDB .gmt file,
# using the DESeq2 Wald-statistic ranking produced by bin/run_deseq2.R.

suppressPackageStartupMessages({
    library(optparse)
    library(fgsea)
    library(dplyr)
    library(readr)
    library(ggplot2)
})

option_list <- list(
    make_option("--ranked_genes", type = "character", help = "ranked_genes.rnk from run_deseq2.R (gene_id<TAB>stat, no header)"),
    make_option("--gmt",          type = "character", help = "MSigDB .gmt gene set file (e.g. Hallmark h.all.*.symbols.gmt)"),
    make_option("--outdir",       type = "character", default = "gsea_results", help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$ranked_genes), !is.null(opt$gmt))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

ranks_df <- read_tsv(opt$ranked_genes, col_names = c("gene_id", "stat"), show_col_types = FALSE)
ranks <- setNames(ranks_df$stat, ranks_df$gene_id)
ranks <- ranks[!is.na(ranks)]

pathways <- gmtPathways(opt$gmt)

fgsea_res <- fgsea(pathways = pathways, stats = ranks, eps = 0.0) %>%
    as.data.frame() %>%
    arrange(padj)

# leadingEdge is a list column -- flatten for a plain TSV
fgsea_out <- fgsea_res %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ","))
write_tsv(fgsea_out, file.path(opt$outdir, "gsea_results.tsv"))

top_pathways <- fgsea_res %>% arrange(padj) %>% head(20)
if (nrow(top_pathways) > 0) {
    plot <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = padj < 0.05)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("grey60", "firebrick"), na.value = "grey60") +
        theme_bw() +
        labs(title = "Top 20 pathways by padj", x = NULL, y = "Normalized Enrichment Score", fill = "padj < 0.05")
    ggsave(file.path(opt$outdir, "gsea_top_pathways.pdf"), plot, width = 8, height = 6)
}

cat("fgsea analysis complete:", opt$outdir, "\n")
