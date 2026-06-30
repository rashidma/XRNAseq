include { RMATS_TURBO } from '../modules/rmats_turbo.nf'

// Groups the shared bam_indexed channel into the two factor levels named by
// --contrast, writes rMATS-turbo's b1.txt/b2.txt BAM-list files, and runs
// once per contrast (not per sample) -- structurally different from the
// per-sample fan-out used by the variant-calling and fusion branches.
workflow SPLICING {
    take:
    bam_indexed       // tuple(sample_id, bam, bai)
    samplesheet_rows  // channel of samplesheet row maps (sample, condition, batch, ...)
    contrast          // val: "column:levelA_vs_levelB"
    gtf
    read_length
    single_end
    rmats_fdr

    main:
    def (column, level_pair) = contrast.split(':', 2)
    def (level1, level2)     = level_pair.split('_vs_', 2)

    condition_ch = samplesheet_rows.map { row -> tuple(row.sample, row[column]) }

    grouped = bam_indexed
        .map { sample_id, bam, bai -> tuple(sample_id, bam) }
        .join(condition_ch)
        .branch {
            grp1: it[2] == level1
            grp2: it[2] == level2
        }

    b1_txt = grouped.grp1
        .map { sample_id, bam, group -> bam.toString() }
        .collect()
        .map { paths -> paths.join(',') }
        .collectFile(name: 'b1.txt', newLine: false)

    b2_txt = grouped.grp2
        .map { sample_id, bam, group -> bam.toString() }
        .collect()
        .map { paths -> paths.join(',') }
        .collectFile(name: 'b2.txt', newLine: false)

    reads_type = single_end ? 'single' : 'paired'

    RMATS_TURBO(b1_txt, b2_txt, gtf, read_length, reads_type, rmats_fdr)

    emit:
    events      = RMATS_TURBO.out.events
    significant = RMATS_TURBO.out.significant
}
