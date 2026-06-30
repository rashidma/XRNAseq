process SAMTOOLS_SORT_INDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/bam", mode: 'copy'

    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
    label 'process_medium'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: sorted_bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    """
}
