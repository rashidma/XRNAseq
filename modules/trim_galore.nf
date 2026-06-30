process TRIM_GALORE {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*trimming_report.txt"

    container 'quay.io/biocontainers/trim-galore:0.6.11--hdfd78af_0'
    label 'process_medium'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_val_{1,2}.fq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("*trimming_report.txt"), emit: log

    script:
    def is_paired = reads instanceof List && reads.size() == 2
    if (is_paired) {
        """
        trim_galore --paired --gzip --cores ${task.cpus} ${reads[0]} ${reads[1]}
        """
    } else {
        """
        trim_galore --gzip --cores ${task.cpus} ${reads[0]}
        """
    }
}
