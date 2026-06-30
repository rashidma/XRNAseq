process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    container 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'
    label 'process_low'

    input:
    path '*'

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .
    """
}
