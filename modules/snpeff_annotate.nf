process SNPEFF_ANNOTATE {
    tag "$sample_id"
    publishDir "${params.outdir}/variant_calling/annotated", mode: 'copy'

    // snpEff/SnpSift chosen over VEP: self-contained DB auto-download,
    // simpler to keep working identically across local/Slurm/cloud than
    // VEP's offline cache directory. Verify this tag is current before
    // relying on it: https://quay.io/repository/biocontainers/snpeff?tab=tags
    container 'quay.io/biocontainers/snpeff:5.2--hdfd78af_1'
    label 'process_medium'

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    val snpeff_db // e.g. 'GRCh38.105' -- genome build string, not hardcoded

    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf"), emit: vcf
    path "${sample_id}.snpeff_summary.html", emit: summary

    script:
    """
    snpEff -v -stats ${sample_id}.snpeff_summary.html ${snpeff_db} ${vcf} > ${sample_id}.annotated.vcf
    """
}
