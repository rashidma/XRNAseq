process GATK_VARIANTFILTRATION {
    tag "$sample_id"
    publishDir "${params.outdir}/variant_calling/filtered", mode: 'copy'

    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1'
    label 'process_medium'

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path genome_fasta
    path fai
    path dict

    output:
    tuple val(sample_id), path("${sample_id}.filtered.vcf.gz"), path("${sample_id}.filtered.vcf.gz.tbi"), emit: vcf

    script:
    // GATK's documented RNA-seq-specific simplified hard filter (from
    // gatk-workflows/gatk4-rnaseq-germline-snps-indels) -- deliberately
    // lighter than the full DNA germline filter set, since there isn't
    // enough RNA-seq validation data for a VQSR-style model. --window
    // 35 --cluster 3 flags clusters of >=3 SNPs in a 35bp window
    // (SnpCluster), targeting RNA-editing-like artifacts near splice
    // junctions. FILTER is annotated, not subsetted -- PASS/FAIL stays
    // auditable in the output VCF.
    """
    gatk VariantFiltration \\
        -R ${genome_fasta} \\
        -V ${vcf} \\
        --window 35 \\
        --cluster 3 \\
        --filter-name FS --filter-expression "FS > 30.0" \\
        --filter-name QD --filter-expression "QD < 2.0" \\
        -O ${sample_id}.filtered.vcf.gz
    """
}
