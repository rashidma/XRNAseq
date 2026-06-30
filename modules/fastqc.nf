process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    container 'biocontainers/fastqc:v0.11.9_cv8'
    label 'process_low'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.zip"), emit: zip
    tuple val(sample_id), path("*.html"), emit: html

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}
