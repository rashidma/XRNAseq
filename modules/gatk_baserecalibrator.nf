process GATK_BASERECALIBRATOR {
    tag "$sample_id"

    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1'
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome_fasta
    path fai
    path dict
    path known_sites
    path known_sites_index

    output:
    tuple val(sample_id), path("${sample_id}.recal.table"), emit: table

    script:
    """
    gatk BaseRecalibrator \\
        -R ${genome_fasta} \\
        -I ${bam} \\
        --known-sites ${known_sites} \\
        -O ${sample_id}.recal.table
    """
}
