process GATK_SPLITNCIGARREADS {
    tag "$sample_id"

    container 'quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1'
    label 'process_high_memory'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome_fasta
    path fai
    path dict

    output:
    tuple val(sample_id), path("${sample_id}.split.bam"), path("${sample_id}.split.bai"), emit: bam

    script:
    // RNA-seq-specific step with no DNA-seq equivalent: splits reads with N
    // CIGAR operations (spliced reads) into multiple reads, hard-clipping
    // into the splice junction. Required before HaplotypeCaller, otherwise
    // spliced reads produce spurious "indels" at exon boundaries.
    """
    gatk SplitNCigarReads \\
        -R ${genome_fasta} \\
        -I ${bam} \\
        -O ${sample_id}.split.bam
    """
}
