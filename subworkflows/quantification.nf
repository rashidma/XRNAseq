include { FEATURECOUNTS } from '../modules/featurecounts.nf'

workflow QUANTIFICATION {
    take:
    bam_indexed // tuple(sample_id, bam, bai)
    gtf

    main:
    FEATURECOUNTS(bam_indexed.map { it[1] }.collect(), gtf)

    emit:
    counts  = FEATURECOUNTS.out.counts
    summary = FEATURECOUNTS.out.summary
}
