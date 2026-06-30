include { SAMTOOLS_FAIDX; GATK_CREATESEQUENCEDICTIONARY } from '../modules/gatk_prepare_reference.nf'
include { GATK_MARKDUPLICATES }    from '../modules/gatk_markduplicates.nf'
include { GATK_SPLITNCIGARREADS } from '../modules/gatk_splitncigarreads.nf'
include { GATK_BASERECALIBRATOR } from '../modules/gatk_baserecalibrator.nf'
include { GATK_APPLYBQSR }        from '../modules/gatk_applybqsr.nf'
include { GATK_HAPLOTYPECALLER }  from '../modules/gatk_haplotypecaller.nf'
include { GATK_VARIANTFILTRATION } from '../modules/gatk_variantfiltration.nf'
include { SNPEFF_ANNOTATE }        from '../modules/snpeff_annotate.nf'

// GATK Best Practices RNA-seq short variant discovery, consuming the shared
// bam_indexed channel per-sample (no re-alignment): MarkDuplicates ->
// SplitNCigarReads -> BQSR -> HaplotypeCaller -> VariantFiltration ->
// optional snpEff annotation. Per-sample VCFs, not joint genotyping.
workflow VARIANT_CALLING {
    take:
    bam_indexed     // tuple(sample_id, bam, bai)
    genome_fasta
    known_sites
    known_sites_index
    snpeff_db
    run_annotation  // val: bool

    main:
    SAMTOOLS_FAIDX(genome_fasta)
    GATK_CREATESEQUENCEDICTIONARY(genome_fasta)
    fai  = SAMTOOLS_FAIDX.out.fai
    dict = GATK_CREATESEQUENCEDICTIONARY.out.dict

    GATK_MARKDUPLICATES(bam_indexed)
    GATK_SPLITNCIGARREADS(GATK_MARKDUPLICATES.out.bam, genome_fasta, fai, dict)
    GATK_BASERECALIBRATOR(GATK_SPLITNCIGARREADS.out.bam, genome_fasta, fai, dict, known_sites, known_sites_index)

    applybqsr_input = GATK_SPLITNCIGARREADS.out.bam.join(GATK_BASERECALIBRATOR.out.table)
    GATK_APPLYBQSR(applybqsr_input, genome_fasta, fai, dict)

    GATK_HAPLOTYPECALLER(GATK_APPLYBQSR.out.bam, genome_fasta, fai, dict)
    GATK_VARIANTFILTRATION(GATK_HAPLOTYPECALLER.out.vcf, genome_fasta, fai, dict)

    filtered_vcf  = GATK_VARIANTFILTRATION.out.vcf
    annotated_vcf = Channel.empty()
    // run_annotation defaults to true, but snpeff_db has no default -- skip
    // annotation with a warning rather than failing the whole run if it's
    // not set, leaving the unannotated filtered VCFs as the usable output.
    if (run_annotation && !snpeff_db) {
        log.warn "--run_annotation is true but --snpeff_db was not provided -- skipping snpEff annotation. Pass e.g. --snpeff_db GRCh38.105 to enable it."
    } else if (run_annotation) {
        SNPEFF_ANNOTATE(filtered_vcf, snpeff_db)
        annotated_vcf = SNPEFF_ANNOTATE.out.vcf
    }

    emit:
    filtered_vcf  = filtered_vcf
    annotated_vcf = annotated_vcf
}
