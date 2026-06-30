include { ARRIBA } from '../modules/arriba.nf'

// Consumes the shared bam_indexed channel directly (Arriba reads chimeric
// records out of the same coordinate-sorted BAM produced by the one shared
// STAR_ALIGN_2PASS run) -- no re-alignment.
workflow FUSION {
    take:
    bam_indexed   // tuple(sample_id, bam, bai)
    genome_fasta
    gtf
    blacklist
    known_fusions

    main:
    ARRIBA(bam_indexed, genome_fasta, gtf, blacklist, known_fusions)

    emit:
    fusions   = ARRIBA.out.fusions
    discarded = ARRIBA.out.discarded
}
