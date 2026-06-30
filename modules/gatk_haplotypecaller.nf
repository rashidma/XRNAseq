process GATK_HAPLOTYPECALLER {
    tag "$sample_id"
    publishDir "${params.outdir}/variant_calling/vcf", mode: 'copy'

    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1'
    label 'process_high'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome_fasta
    path fai
    path dict

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf

    script:
    // Per-sample calling, not joint genotyping -- a deliberate choice for
    // this cohort style. For cohort-level joint calling later, swap this for
    // -ERC GVCF + CombineGVCFs/GenotypeGVCFs.
    """
    gatk HaplotypeCaller \\
        -R ${genome_fasta} \\
        -I ${bam} \\
        -O ${sample_id}.vcf.gz \\
        --dont-use-soft-clipped-bases \\
        -stand-call-conf 20
    """
}
