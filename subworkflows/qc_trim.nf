include { FASTQC as FASTQC_RAW } from '../modules/fastqc.nf'
include { TRIM_GALORE }           from '../modules/trim_galore.nf'

workflow QC_TRIM {
    take:
    read_pairs_ch // tuple(sample_id, [reads])

    main:
    FASTQC_RAW(read_pairs_ch)
    TRIM_GALORE(read_pairs_ch)

    emit:
    trimmed_reads = TRIM_GALORE.out.trimmed_reads
    fastqc_zip    = FASTQC_RAW.out.zip
    trim_log      = TRIM_GALORE.out.log
}
