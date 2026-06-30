process FEATURECOUNTS {
    tag "featurecounts"
    publishDir "${params.outdir}/counts", mode: 'copy'

    container 'quay.io/biocontainers/subread:2.0.3--h7132678_2'
    label 'process_high'

    input:
    path bam_files
    path gtf

    output:
    path "gene_counts.txt", emit: counts
    path "gene_counts.txt.summary", emit: summary

    script:
    """
    featureCounts \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -o gene_counts.txt \\
        -p --countReadPairs --primary \\
        ${bam_files.join(' ')}
    """
}
