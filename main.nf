#!/usr/bin/env nextflow
/*
 * Custom multi-modal RNA-seq pipeline
 * Author: Dr. Mamoon Rashid, KAIMRC
 *
 * Always runs: FastQC -> Trim Galore -> STAR (2-pass, chimeric-aware) -> samtools sort/index
 * Then, gated by --run_* flags, fans the one shared aligned BAM per sample out to:
 *   - featureCounts -> DESeq2 -> fgsea/GSEA           (--run_dge / --run_gsea, default on)
 *   - GATK RNA-seq variant calling                     (--run_variant_calling)
 *   - Arriba fusion detection                            (--run_fusion)
 *   - rMATS-turbo alternative splicing                    (--run_splicing)
 * -> MultiQC aggregates QC across whichever branches ran.
 *
 * Usage:
 *   nextflow run main.nf --samplesheet samplesheet.csv --genome_fasta ref.fa \
 *     --gtf annotation.gtf --outdir results -profile docker
 *
 * Params are defined centrally in nextflow.config; this file only validates
 * and wires subworkflows together.
 */

nextflow.enable.dsl=2

include { QC_TRIM }                  from './subworkflows/qc_trim.nf'
include { ALIGNMENT }                 from './subworkflows/alignment.nf'
include { QUANTIFICATION }            from './subworkflows/quantification.nf'
include { DIFFERENTIAL_EXPRESSION }   from './subworkflows/differential_expression.nf'
include { VARIANT_CALLING }           from './subworkflows/variant_calling.nf'
include { FUSION }                    from './subworkflows/fusion.nf'
include { SPLICING }                  from './subworkflows/splicing.nf'
include { MULTIQC }                   from './modules/multiqc.nf'

def helpMessage() {
    log.info """
    Custom multi-modal RNA-seq pipeline
    ------------------------------------
    Required:
      --samplesheet   CSV with columns: sample,fastq_1,fastq_2,condition,batch
      --gtf           Gene annotation GTF file
      --genome_fasta  Reference genome FASTA (required if --star_index not given;
                       ALSO required by --run_variant_calling / --run_fusion even
                       if a prebuilt --star_index is supplied)
      OR
      --star_index    Pre-built STAR genome index directory

    Optional (core):
      --outdir        Output directory (default: results)
      --single_end    Set true for single-end reads (default: false)

    Modality toggles (see README for the reference files each one needs):
      --run_dge=true/false              gene-level featureCounts -> DESeq2 (default: true)
      --run_gsea=true/false             fgsea pathway enrichment on DESeq2 results (default: true)
      --run_variant_calling=true/false  GATK RNA-seq short variant discovery (default: false)
      --run_annotation=true/false       snpEff annotation of filtered VCFs (default: true, only if run_variant_calling)
      --run_fusion=true/false           Arriba fusion detection (default: false)
      --run_splicing=true/false         rMATS-turbo differential splicing (default: false)

    DGE/splicing/GSEA params:
      --contrast "condition:treated_vs_control"
      --design_formula "~ condition"   (or "~ batch + condition")
      --msigdb_gmt /path/to/h.all.symbols.gmt   (required if --run_gsea)

    Variant calling params:
      --known_sites /path/or/gs://.../dbsnp.vcf.gz   (.tbi must sit alongside it)
      --snpeff_db GRCh38.105                          (required if --run_annotation)

    Fusion params:
      --arriba_blacklist / --arriba_known_fusions     (from Arriba's GitHub release tarball)

    Splicing params:
      --read_length 100   (required if --run_splicing; must match real FASTQ read length)
      --rmats_fdr 0.05

    Example:
      nextflow run main.nf \\
        --samplesheet samplesheet.csv \\
        --genome_fasta genome.fa \\
        --gtf annotation.gtf \\
        --outdir results -profile docker
    """
}

// Nextflow (modern versions) does not allow bare statements at script top
// level alongside process/workflow declarations, so validation lives in a
// function called from inside the workflow block below.
def checkParams() {
    if (params.help) {
        helpMessage()
        exit 0
    }
    if (!params.samplesheet || !params.gtf || (!params.genome_fasta && !params.star_index)) {
        log.error "Missing required parameters. Run with --help for usage."
        exit 1
    }
    if ((params.run_variant_calling || params.run_fusion) && !params.genome_fasta) {
        log.error "--run_variant_calling and --run_fusion both need --genome_fasta (the raw reference FASTA), even when a prebuilt --star_index is supplied."
        exit 1
    }
    if (params.run_variant_calling && !params.known_sites) {
        log.error "--run_variant_calling requires --known_sites (a dbSNP-style known-sites VCF for BQSR)."
        exit 1
    }
    if (params.run_fusion && (!params.arriba_blacklist || !params.arriba_known_fusions)) {
        log.error "--run_fusion requires --arriba_blacklist and --arriba_known_fusions (download from Arriba's GitHub release tarball for your genome build)."
        exit 1
    }
    if (params.run_splicing && !params.read_length) {
        log.error "--run_splicing requires --read_length (must match your FASTQ read length)."
        exit 1
    }
    // Note: --run_gsea defaults to true but is skipped with a warning (not a
    // hard failure) if --msigdb_gmt isn't supplied -- see subworkflows/differential_expression.nf.
}

workflow {

    checkParams()

    def genome_fasta_ch = params.genome_fasta ? file(params.genome_fasta) : null
    def gtf_ch           = file(params.gtf)

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def reads = params.single_end \
                ? [ file(row.fastq_1) ] \
                : [ file(row.fastq_1), file(row.fastq_2) ]
            tuple(row.sample, reads)
        }
        .set { read_pairs_ch }

    // Re-read the samplesheet separately for the splicing branch, which needs
    // the condition column per-sample rather than just the FASTQ paths.
    samplesheet_rows_ch = Channel.fromPath(params.samplesheet).splitCsv(header: true)

    // ---- Always-on: QC, trimming, one shared 2-pass+chimeric alignment per sample ----
    QC_TRIM(read_pairs_ch)
    ALIGNMENT(QC_TRIM.out.trimmed_reads, genome_fasta_ch, gtf_ch, params.star_index, params.star_genome_sa_index_n_bases)

    multiqc_inputs = QC_TRIM.out.fastqc_zip.map { it[1] }
        .mix(QC_TRIM.out.trim_log.map { it[1] })
        .mix(ALIGNMENT.out.star_log.map { it[1] })

    // ---- Differential expression + GSEA ----
    if (params.run_dge) {
        QUANTIFICATION(ALIGNMENT.out.bam_indexed, gtf_ch)
        DIFFERENTIAL_EXPRESSION(
            QUANTIFICATION.out.counts,
            file(params.samplesheet),
            params.contrast,
            params.design_formula,
            params.run_gsea,
            params.msigdb_gmt
        )
        multiqc_inputs = multiqc_inputs.mix(QUANTIFICATION.out.summary)
    }

    // ---- GATK RNA-seq variant calling ----
    if (params.run_variant_calling) {
        VARIANT_CALLING(
            ALIGNMENT.out.bam_indexed,
            genome_fasta_ch,
            file(params.known_sites),
            file("${params.known_sites}.tbi"),
            params.snpeff_db,
            params.run_annotation
        )
    }

    // ---- Arriba fusion detection ----
    if (params.run_fusion) {
        FUSION(
            ALIGNMENT.out.bam_indexed,
            genome_fasta_ch,
            gtf_ch,
            file(params.arriba_blacklist),
            file(params.arriba_known_fusions)
        )
    }

    // ---- rMATS-turbo alternative splicing ----
    if (params.run_splicing) {
        SPLICING(
            ALIGNMENT.out.bam_indexed,
            samplesheet_rows_ch,
            params.contrast,
            gtf_ch,
            params.read_length,
            params.single_end,
            params.rmats_fdr
        )
    }

    MULTIQC(multiqc_inputs.collect())

    workflow.onComplete = {
        log.info "Pipeline completed at: $workflow.complete"
        log.info "Execution status: ${ workflow.success ? 'OK' : 'FAILED' }"
        log.info "Results in: ${params.outdir}"
    }
}
