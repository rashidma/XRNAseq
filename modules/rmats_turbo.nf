process RMATS_TURBO {
    tag "rmats"
    publishDir "${params.outdir}/rmats", mode: 'copy'

    // No usable rmats-turbo tag exists on quay.io/biocontainers as of
    // implementation time -- the upstream-maintained image is the right
    // call here. Verify this tag is still current before relying on it:
    // https://hub.docker.com/r/xinglab/rmats/tags
    container 'xinglab/rmats:v4.3.0'
    label 'process_high'

    input:
    path b1_txt
    path b2_txt
    path gtf
    val  read_length
    val  reads_type // 'paired' or 'single'
    val  fdr_threshold

    output:
    path "rmats_out/*.MATS.JC.txt",        emit: events
    path "rmats_out/*_significant.txt",    emit: significant
    path "rmats_out",                       emit: outdir

    script:
    """
    mkdir -p rmats_tmp
    rmats.py \\
        --b1 ${b1_txt} --b2 ${b2_txt} \\
        --gtf ${gtf} \\
        --od rmats_out --tmp rmats_tmp \\
        -t ${reads_type} \\
        --readLength ${read_length} \\
        --nthread ${task.cpus}

    for f in rmats_out/*.MATS.JC.txt; do
        out="\${f%.MATS.JC.txt}_significant.txt"
        python3 -c "
import csv, sys
with open('\$f') as inp, open('\$out', 'w') as outp:
    reader = csv.reader(inp, delimiter='\t')
    writer = csv.writer(outp, delimiter='\t')
    header = next(reader)
    writer.writerow(header)
    fdr_idx = header.index('FDR')
    for row in reader:
        if row[fdr_idx] != 'NA' and float(row[fdr_idx]) < ${fdr_threshold}:
            writer.writerow(row)
"
    done
    """
}
