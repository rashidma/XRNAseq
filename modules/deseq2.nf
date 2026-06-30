process DESEQ2 {
    tag "deseq2"
    publishDir "${params.outdir}/deseq2", mode: 'copy'

    // Custom image -- build containers/r_analysis/Dockerfile and push it
    // somewhere you control, then set this to your own tag, e.g.:
    //   container 'ghcr.io/<you>/xrnaseq-r-analysis:1.0.0'
    container 'PLACEHOLDER_VERIFY_AT_IMPLEMENTATION:xrnaseq-r-analysis'
    label 'process_medium'

    input:
    path counts_matrix
    path samplesheet
    val  contrast
    val  design_formula

    output:
    path "deseq2_results.tsv",    emit: results
    path "normalized_counts.tsv", emit: normalized_counts
    path "ranked_genes.rnk",      emit: ranked_genes
    path "dds.rds",               emit: dds
    path "*.pdf",                 emit: plots

    script:
    """
    run_deseq2.R \\
        --counts ${counts_matrix} \\
        --samplesheet ${samplesheet} \\
        --contrast "${contrast}" \\
        --design "${design_formula}" \\
        --outdir .
    """
}
