#!/usr/bin/env Rscript
# Driver script for modules/deseq2.nf.
#
# Takes a featureCounts gene-count matrix + the pipeline samplesheet
# (condition/batch columns) and runs a standard DESeq2 differential
# expression analysis for a single two-group contrast.

suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(dplyr)
    library(readr)
    library(tibble)
    library(ggplot2)
})

option_list <- list(
    make_option("--counts",      type = "character", help = "featureCounts gene_counts.txt matrix"),
    make_option("--samplesheet", type = "character", help = "pipeline samplesheet CSV (sample,fastq_1,fastq_2,condition,batch)"),
    make_option("--contrast",    type = "character", help = "e.g. condition:treated_vs_control"),
    make_option("--design",      type = "character", default = "~ condition", help = "DESeq2 design formula"),
    make_option("--outdir",      type = "character", default = "deseq2_results", help = "output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$counts), !is.null(opt$samplesheet), !is.null(opt$contrast))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# --- Parse --contrast "column:levelA_vs_levelB" ---
contrast_parts <- strsplit(opt$contrast, ":")[[1]]
contrast_column <- contrast_parts[1]
levels_parts    <- strsplit(contrast_parts[2], "_vs_")[[1]]
level_numerator   <- levels_parts[1]
level_denominator <- levels_parts[2]

# --- Load featureCounts matrix ---
# featureCounts output has a leading "# Program:featureCounts ..." comment
# line, then a header: Geneid,Chr,Start,End,Strand,Length,<bam1>,<bam2>,...
# Sample columns are full BAM paths -- strip directory + the
# ".sorted.bam" suffix added by modules/samtools.nf to recover sample IDs.
counts_raw <- read.delim(opt$counts, header = TRUE, comment.char = "#", check.names = FALSE)
gene_ids <- counts_raw$Geneid
count_cols <- setdiff(colnames(counts_raw), c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
count_matrix <- as.matrix(counts_raw[, count_cols, drop = FALSE])
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- basename(colnames(count_matrix))
colnames(count_matrix) <- sub("\\.sorted\\.bam$", "", colnames(count_matrix))

# --- Load samplesheet metadata ---
samplesheet <- read_csv(opt$samplesheet, show_col_types = FALSE)
stopifnot("condition" %in% colnames(samplesheet))
col_data <- samplesheet %>%
    distinct(sample, .keep_all = TRUE) %>%
    column_to_rownames("sample")

common_samples <- intersect(colnames(count_matrix), rownames(col_data))
if (length(common_samples) == 0) {
    stop("No overlap between count matrix columns and samplesheet 'sample' values -- check that featureCounts BAM names match samplesheet sample IDs.")
}
count_matrix <- count_matrix[, common_samples, drop = FALSE]
col_data <- col_data[common_samples, , drop = FALSE]

for (col in all.vars(as.formula(opt$design))) {
    col_data[[col]] <- as.factor(col_data[[col]])
}

# --- Run DESeq2 ---
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = as.formula(opt$design))
dds <- DESeq(dds)

res <- results(dds, contrast = c(contrast_column, level_numerator, level_denominator))
res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)

write_tsv(res_df, file.path(opt$outdir, "deseq2_results.tsv"))
write_tsv(
    as.data.frame(counts(dds, normalized = TRUE)) %>% rownames_to_column("gene_id"),
    file.path(opt$outdir, "normalized_counts.tsv")
)
saveRDS(dds, file.path(opt$outdir, "dds.rds"))

# Ranked gene list for fgsea: DESeq2's Wald statistic is a simple, defensible
# ranking metric (sign and magnitude both carry information fgsea needs).
res_df %>%
    filter(!is.na(stat)) %>%
    select(gene_id, stat) %>%
    arrange(desc(stat)) %>%
    write_tsv(file.path(opt$outdir, "ranked_genes.rnk"), col_names = FALSE)

# --- Plots ---
# vst() subsamples ~1000 rows to fit the dispersion trend and errors out below
# that; varianceStabilizingTransformation() fits on all rows instead, which is
# DESeq2's own documented fallback for small gene counts (e.g. a chr21-only
# test run or a small targeted panel).
vsd <- tryCatch(
    vst(dds, blind = TRUE),
    error = function(e) varianceStabilizingTransformation(dds, blind = TRUE)
)
pca_data <- plotPCA(vsd, intgroup = contrast_column, returnData = TRUE)
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = .data[[contrast_column]])) +
    geom_point(size = 3) +
    theme_bw() +
    labs(title = "PCA (VST-transformed counts)")
ggsave(file.path(opt$outdir, "pca_plot.pdf"), pca_plot, width = 6, height = 5)

volcano_plot <- ggplot(res_df, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(aes(color = padj < 0.05), alpha = 0.6, na.rm = TRUE) +
    scale_color_manual(values = c("grey60", "firebrick"), na.value = "grey60") +
    theme_bw() +
    labs(title = paste0(level_numerator, " vs ", level_denominator), color = "padj < 0.05")
ggsave(file.path(opt$outdir, "volcano_plot.pdf"), volcano_plot, width = 6, height = 5)

cat("DESeq2 analysis complete:", opt$outdir, "\n")
