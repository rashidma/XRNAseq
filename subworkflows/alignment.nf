include { STAR_INDEX }          from '../modules/star_index.nf'
include { STAR_ALIGN_2PASS }    from '../modules/star_align_2pass.nf'
include { SAMTOOLS_SORT_INDEX } from '../modules/samtools.nf'

// Builds (or reuses) the STAR index, runs the single shared 2-pass +
// chimeric-aware alignment per sample, and emits one indexed BAM channel
// that every downstream modality (DGE, variant calling, fusion, splicing)
// consumes independently -- STAR runs exactly once per sample regardless of
// how many --run_* toggles are on.
workflow ALIGNMENT {
    take:
    trimmed_reads      // tuple(sample_id, [reads])
    genome_fasta       // path or null
    gtf                // path
    star_index         // path or null
    sa_index_n_bases   // val: STAR --genomeSAindexNbases (only used if star_index is null)

    main:
    if (star_index) {
        // A user-supplied path needs .first() to become a reusable value
        // channel (one STAR_ALIGN_2PASS call per sample, same index each time).
        index_ch = Channel.fromPath(star_index, checkIfExists: true).first()
    } else {
        // A single process call's output is already a value channel.
        STAR_INDEX(genome_fasta, gtf, sa_index_n_bases)
        index_ch = STAR_INDEX.out.index
    }

    STAR_ALIGN_2PASS(trimmed_reads, index_ch, gtf)
    SAMTOOLS_SORT_INDEX(STAR_ALIGN_2PASS.out.bam)

    emit:
    bam_indexed      = SAMTOOLS_SORT_INDEX.out.sorted_bam // tuple(sample_id, bam, bai) -- THE shared channel
    star_log         = STAR_ALIGN_2PASS.out.log
    star_gene_counts = STAR_ALIGN_2PASS.out.star_gene_counts
}
