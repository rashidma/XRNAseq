process GSEA_FGSEA {
    tag "gsea"
    publishDir "${params.outdir}/gsea", mode: 'copy'

    // Same custom image as modules/deseq2.nf -- see that file's container note.
    container 'xrnaseq-r-analysis:1.0.0'
    label 'process_low'

    input:
    path ranked_genes
    path msigdb_gmt

    output:
    path "gsea_results.tsv", emit: results
    path "*.pdf", optional: true, emit: plots

    script:
    """
    run_fgsea.R \\
        --ranked_genes ${ranked_genes} \\
        --gmt ${msigdb_gmt} \\
        --outdir .
    """
}
