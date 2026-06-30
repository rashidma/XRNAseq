process STAR_ALIGN_2PASS {
    tag "$sample_id"
    publishDir "${params.outdir}/star_align", mode: 'copy'

    container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'
    label 'process_high_memory'

    input:
    tuple val(sample_id), path(reads)
    path star_index
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.Log.final.out"),                emit: log
    tuple val(sample_id), path("${sample_id}.ReadsPerGene.out.tab"),         emit: star_gene_counts
    tuple val(sample_id), path("${sample_id}.SJ.out.tab"),                  emit: sj_tab

    script:
    /*
     * One alignment pass per sample feeds every downstream consumer:
     * featureCounts/DESeq2, the GATK variant-calling chain, Arriba, and
     * rMATS-turbo. --twopassMode Basic makes "2-pass" a single STAR
     * invocation (STAR reruns pass 2 internally with pass-1 splice
     * junctions), so this is still exactly one process execution per sample.
     *
     * --outSAMattrRGline is required by GATK SplitNCigarReads/HaplotypeCaller
     * downstream, even though samtools/featureCounts never needed it before.
     *
     * --chim* flags are Arriba's own documented recommended STAR parameters;
     * --chimOutType WithinBAM keeps chimeric reads inside this one sorted BAM
     * so Arriba reads them directly with no second STAR output stream.
     */
    """
    STAR \\
        --runMode alignReads \\
        --genomeDir ${star_index} \\
        --readFilesIn ${reads.join(' ')} \\
        --readFilesCommand zcat \\
        --twopassMode Basic \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --sjdbGTFfile ${gtf} \\
        --outFileNamePrefix ${sample_id}. \\
        --runThreadN ${task.cpus} \\
        --quantMode GeneCounts \\
        --outSAMattrRGline ID:${sample_id} SM:${sample_id} PL:ILLUMINA \\
        --chimSegmentMin 10 \\
        --chimOutType WithinBAM SoftClip \\
        --chimJunctionOverhangMin 10 \\
        --chimScoreMin 1 \\
        --chimScoreDropMax 30 \\
        --chimScoreJunctionNonGTAG 0 \\
        --chimScoreSeparation 1 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --chimSegmentReadGapMax 3
    """
}
