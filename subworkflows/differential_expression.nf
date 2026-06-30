include { DESEQ2 }     from '../modules/deseq2.nf'
include { GSEA_FGSEA } from '../modules/gsea_fgsea.nf'

workflow DIFFERENTIAL_EXPRESSION {
    take:
    counts_matrix  // path: featureCounts gene_counts.txt
    samplesheet    // path
    contrast       // val
    design_formula // val
    run_gsea       // val: bool
    msigdb_gmt     // path or null

    main:
    DESEQ2(counts_matrix, samplesheet, contrast, design_formula)

    // run_gsea defaults to true (it follows run_dge), but msigdb_gmt has no
    // default -- skip GSEA with a warning rather than failing the whole
    // pipeline on a zero-config run, so --run_dge alone still behaves like
    // the original pipeline plus DESeq2 results.
    gsea_results = Channel.empty()
    if (run_gsea && !msigdb_gmt) {
        log.warn "--run_gsea is true but --msigdb_gmt was not provided -- skipping GSEA. Download a gene set file (e.g. MSigDB Hallmark) and pass it with --msigdb_gmt to enable it."
    } else if (run_gsea) {
        GSEA_FGSEA(DESEQ2.out.ranked_genes, file(msigdb_gmt))
        gsea_results = GSEA_FGSEA.out.results
    }

    emit:
    deseq2_results = DESEQ2.out.results
    gsea_results    = gsea_results
}
