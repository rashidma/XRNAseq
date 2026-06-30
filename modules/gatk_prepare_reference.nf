// Reference indexing for the GATK variant-calling branch -- run once per
// pipeline run (not per sample) and cached/reused across all samples.

process SAMTOOLS_FAIDX {
    tag "faidx"
    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
    label 'process_low'

    input:
    path genome_fasta

    output:
    path "${genome_fasta}.fai", emit: fai

    script:
    """
    samtools faidx ${genome_fasta}
    """
}

process GATK_CREATESEQUENCEDICTIONARY {
    tag "dict"
    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1'
    label 'process_low'

    input:
    path genome_fasta

    output:
    path "${genome_fasta.baseName}.dict", emit: dict

    script:
    """
    gatk CreateSequenceDictionary -R ${genome_fasta} -O ${genome_fasta.baseName}.dict
    """
}
