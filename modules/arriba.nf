process ARRIBA {
    tag "$sample_id"
    publishDir "${params.outdir}/arriba", mode: 'copy'

    container 'quay.io/biocontainers/arriba:2.5.1--h87b9561_0'
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome_fasta
    path gtf
    path blacklist
    path known_fusions

    output:
    tuple val(sample_id), path("${sample_id}.fusions.tsv"),          emit: fusions
    tuple val(sample_id), path("${sample_id}.fusions.discarded.tsv"), emit: discarded

    script:
    """
    arriba \\
        -x ${bam} \\
        -a ${genome_fasta} \\
        -g ${gtf} \\
        -b ${blacklist} \\
        -k ${known_fusions} \\
        -o ${sample_id}.fusions.tsv \\
        -O ${sample_id}.fusions.discarded.tsv
    """
}
