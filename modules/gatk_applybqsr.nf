process GATK_APPLYBQSR {
    tag "$sample_id"
    publishDir "${params.outdir}/variant_calling/bqsr", mode: 'copy'

    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1'
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam), path(bai), path(recal_table)
    path genome_fasta
    path fai
    path dict

    output:
    tuple val(sample_id), path("${sample_id}.bqsr.bam"), path("${sample_id}.bqsr.bai"), emit: bam

    script:
    """
    gatk ApplyBQSR \\
        -R ${genome_fasta} \\
        -I ${bam} \\
        --bqsr-recal-file ${recal_table} \\
        -O ${sample_id}.bqsr.bam
    """
}
