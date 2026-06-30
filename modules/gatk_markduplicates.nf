process GATK_MARKDUPLICATES {
    tag "$sample_id"
    publishDir "${params.outdir}/variant_calling/markduplicates", mode: 'copy', pattern: "*.metrics.txt"

    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1'
    label 'process_high_memory'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.md.bam"), path("${sample_id}.md.bai"), emit: bam
    path "${sample_id}.md.metrics.txt", emit: metrics

    script:
    """
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${sample_id}.md.bam \\
        -M ${sample_id}.md.metrics.txt \\
        --CREATE_INDEX true
    """
}
