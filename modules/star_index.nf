process STAR_INDEX {
    tag "star_index"
    publishDir "${params.outdir}/star_index", mode: 'copy'

    container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'
    label 'process_high_memory'

    input:
    path genome_fasta
    path gtf
    val  sa_index_n_bases // STAR recommends scaling this down from the default 14 for small genomes (e.g. a single chromosome test subset) to avoid wasting index memory

    output:
    path "star_index", emit: index

    script:
    """
    mkdir star_index
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang 100 \\
        --genomeSAindexNbases ${sa_index_n_bases} \\
        --runThreadN ${task.cpus}
    """
}
