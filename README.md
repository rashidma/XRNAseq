# XRNAseq — Custom Multi-Modal RNA-seq Pipeline (Nextflow DSL2)

A hand-rolled (no nf-core) Nextflow DSL2 pipeline that runs unmodified on a local Ubuntu
workstation, KAIMRC's Slurm HPC, and Google Cloud Batch, and performs five RNA-seq analysis
modalities from one shared alignment per sample:

- **Differential gene expression** — featureCounts → DESeq2
- **RNA-seq variant calling** — GATK Best Practices short variant discovery
- **Fusion detection** — Arriba
- **Alternative splicing** — rMATS-turbo
- **Pathway / GSEA enrichment** — fgsea against MSigDB gene sets, automatic on DESeq2 results

Always-on backbone: **FastQC → Trim Galore → STAR (2-pass, chimeric-aware) → samtools sort/index
→ MultiQC**. Everything else fans out from that one alignment per sample — see
[Architecture](#architecture) below.

## Repository structure

```
XRNAseq/
├── main.nf                       # entrypoint: validates params, wires subworkflows together
├── nextflow.config                # params + profile wiring + manifest + reporting
├── conf/
│   ├── base.config                 # resource labels (process_low/medium/high/high_memory)
│   ├── local.config                 # local workstation executor tuning
│   ├── slurm.config                 # KAIMRC Slurm settings -- TODO placeholders, fill in before use
│   ├── google_batch.config           # Google Cloud Batch settings -- TODO placeholders
│   └── test.config                   # chr21-subset smoke test wiring
├── modules/                          # one process per file
│   ├── fastqc.nf, trim_galore.nf, star_index.nf, star_align_2pass.nf, samtools.nf,
│   │   featurecounts.nf, multiqc.nf            # core QC/align/count backbone
│   ├── gatk_*.nf, snpeff_annotate.nf            # variant calling chain
│   ├── arriba.nf                                  # fusion detection
│   ├── rmats_turbo.nf                             # alternative splicing
│   └── deseq2.nf, gsea_fgsea.nf                    # DGE + GSEA (R wrappers)
├── subworkflows/
│   ├── qc_trim.nf, alignment.nf, quantification.nf      # backbone
│   ├── variant_calling.nf, fusion.nf, splicing.nf       # opt-in modalities
│   └── differential_expression.nf                       # DGE + GSEA
├── assets/
│   ├── samplesheet_example.csv, test_samplesheet.csv
│   └── test_genome.fa, test_annotation.gtf, test_known_sites.vcf.gz  # NOT included -- you generate these, see below
├── bin/
│   ├── run_deseq2.R, run_fgsea.R                # R driver scripts called by deseq2.nf / gsea_fgsea.nf
└── containers/r_analysis/Dockerfile             # custom DESeq2+fgsea image -- you build & host this
```

## Architecture: one alignment per sample, four independent consumers

`STAR_ALIGN_2PASS` runs **once per sample**, with `--twopassMode Basic` (STAR's own 2-pass
re-run, still a single process execution) and Arriba's documented `--chim*` flags always on, so
chimeric reads end up inside the same coordinate-sorted BAM. The resulting indexed BAM channel
(`bam_indexed`) is then consumed independently by up to four branches, with **no re-alignment**:

1. `featureCounts → DESeq2 → fgsea` (`--run_dge` / `--run_gsea`, default **on**)
2. GATK `MarkDuplicates → SplitNCigarReads → BQSR → HaplotypeCaller → VariantFiltration → snpEff` (`--run_variant_calling`, default **off**)
3. `Arriba` (`--run_fusion`, default **off**)
4. `rMATS-turbo` (`--run_splicing`, default **off**)

DGE/GSEA stay on by default so a zero-config run behaves like the original pipeline plus DE/pathway
results for free. Variant calling, fusion, and splicing are opt-in because each needs extra
reference files you won't have staged on a first run, and GATK's chain is by far the most
expensive branch.

## Prerequisites

1. **Nextflow** ≥22.10 (tested against 26.04): `curl -s https://get.nextflow.io | bash`, put the
   binary on your `$PATH`.
2. **Java 17+**.
3. **Docker** (local workstation, Google Batch) or **Singularity** (KAIMRC Slurm — Docker is
   typically not permitted on shared HPC nodes).
4. For the DGE/GSEA branch: build and push the custom R container (see
   [Custom R container](#custom-r-container) below) before your first `--run_dge` run.

## Samplesheet

CSV with `sample,fastq_1,fastq_2,condition,batch` (see `assets/samplesheet_example.csv`).
`condition` is the biological grouping used by DESeq2, rMATS-turbo, and the `--contrast` you pass;
required whenever `--run_dge` or `--run_splicing` is on. `batch` only feeds the DESeq2 design
formula (e.g. `~ batch + condition`) — **rMATS-turbo has no batch-covariate mechanism**, so if your
cohort spans multiple collection sites/timepoints (e.g. a multi-site cohort with real batch
effects), the splicing results have no built-in defense against that confounding the way DGE does.
For single-end data, leave `fastq_2` blank and pass `--single_end true`.

`--contrast "condition:treated_vs_control"` drives both DESeq2's contrast and which samples become
rMATS-turbo's two BAM groups.

## Running it

### Local Ubuntu workstation

```bash
nextflow run main.nf \
  --samplesheet samplesheet.csv \
  --genome_fasta /ref/GRCh38.primary_assembly.fa \
  --gtf /ref/gencode.v44.annotation.gtf \
  --outdir results \
  -profile docker,local
```

### KAIMRC HPC (Slurm)

Fill in the `TODO` placeholders in `conf/slurm.config` first (real partition/queue name, account
string, a shared `singularity.cacheDir` path) — the pipeline will fail to submit jobs otherwise.

```bash
nextflow run main.nf \
  --samplesheet samplesheet.csv \
  --genome_fasta /ref/GRCh38.primary_assembly.fa \
  --gtf /ref/gencode.v44.annotation.gtf \
  --outdir results \
  -profile singularity,slurm
```

(`-profile docker,slurm` is not a useful combination — `slurm.config` forces Singularity on
regardless, since Docker generally isn't permitted on shared HPC nodes.)

### Google Cloud Batch

Fill in the `TODO` placeholders in `conf/google_batch.config` (GCP project, region, GCS bucket for
`workDir`) first, and ensure you're authenticated (`gcloud auth application-default login`).

```bash
nextflow run main.nf \
  --samplesheet samplesheet.csv \
  --genome_fasta gs://your-bucket/ref/GRCh38.primary_assembly.fa \
  --gtf gs://your-bucket/ref/gencode.v44.annotation.gtf \
  --outdir gs://your-bucket/results \
  -profile docker,google_batch
```
Google Batch intermediate files (`workDir`) accrue GCS storage cost until cleaned up — run
`nextflow clean` after each run, or set a bucket lifecycle rule, especially given how much bigger
the GATK/rMATS intermediates are than the original pipeline's.

Reusing an existing STAR index (any environment, faster for repeat runs): swap `--genome_fasta` for
`--star_index /path/to/star_index`. **Note**: `--genome_fasta` is still required even with
`--star_index` if you enable `--run_variant_calling` or `--run_fusion` — GATK and Arriba both need
the raw reference FASTA, not just the STAR index.

Resume after a crash: add `-resume` (Nextflow caches completed tasks).

## Enabling each modality

### Differential expression + GSEA (on by default)
```bash
--contrast "condition:treated_vs_control" \
--design_formula "~ batch + condition" \
--msigdb_gmt /path/to/h.all.v2024.1.Hs.symbols.gmt   # required for GSEA; download from MSigDB (registration required)
```
If `--msigdb_gmt` isn't supplied, GSEA is skipped with a warning rather than failing the run — you
still get DESeq2 results.

### RNA-seq variant calling
```bash
--run_variant_calling true \
--known_sites /path/or/gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
--snpeff_db GRCh38.105
```
The `.tbi` index must sit alongside `--known_sites`. The `gs://genomics-public-data/...` path
resolves directly only when running `-profile google_batch`; for Slurm/local, download a local
copy instead. **Caveat**: RNA-seq variant calls are noisier than DNA-seq calls — GATK's own
position is that there isn't enough RNA-seq validation data for a VQSR-style filter, hence the
lighter `FS > 30 / QD < 2 / SnpCluster` hard filter used here instead of the full DNA germline
filter set. Treat these calls as hypothesis-generating, not clinical-grade, without orthogonal
DNA-seq confirmation. RNA-seq variant calling also structurally can't see variants in genes that
aren't expressed in the sampled tissue/timepoint.

### Fusion detection (Arriba)
```bash
--run_fusion true \
--arriba_blacklist /path/to/blacklist_hg38_GRCh38_vX.X.X.tsv.gz \
--arriba_known_fusions /path/to/known_fusions_hg38_GRCh38_vX.X.X.tsv.gz
```
Download both from [Arriba's GitHub releases](https://github.com/suhrig/arriba/releases) tarball
(`database/`) for your genome build — they ship only in the release tarball, not via `git clone`,
and must be kept in lockstep with the pinned Arriba container version (bump both together
deliberately).

### Alternative splicing (rMATS-turbo)
```bash
--run_splicing true \
--read_length 100 \
--rmats_fdr 0.05
```
`--read_length` must match your actual FASTQ read length — rMATS-turbo's statistical model assumes
uniform read length across all input BAMs; mixed-length libraries within the same `condition` group
will degrade results.

## Outputs

```
results/
├── fastqc/, trimmed/, star_align/, bam/      # backbone QC + alignment
├── counts/                                    # gene_counts.txt featureCounts matrix
├── deseq2/                                     # deseq2_results.tsv, normalized_counts.tsv, dds.rds, PCA/volcano plots
├── gsea/                                        # gsea_results.tsv, top-pathway plot
├── variant_calling/                              # markduplicates/, bqsr/, vcf/, filtered/, annotated/
├── arriba/                                        # <sample>.fusions.tsv, <sample>.fusions.discarded.tsv
├── rmats/                                         # SE/A5SS/A3SS/MXE/RI MATS.JC.txt + *_significant.txt at FDR < rmats_fdr
├── multiqc/                                       # aggregated HTML report
└── pipeline_info/                                  # execution report, timeline, trace, DAG
```

## Custom R container

No existing public image was confirmed to contain the exact DESeq2 + fgsea + apeglm +
annotation-package combination this pipeline needs, so `containers/r_analysis/Dockerfile` is
something you build and host yourself once:

```bash
docker build -t ghcr.io/<you>/xrnaseq-r-analysis:1.0.0 -f containers/r_analysis/Dockerfile .
docker push ghcr.io/<you>/xrnaseq-r-analysis:1.0.0
```

Then set that tag in `modules/deseq2.nf` and `modules/gsea_fgsea.nf` (currently
`PLACEHOLDER_VERIFY_AT_IMPLEMENTATION:xrnaseq-r-analysis`).

## Container tags — verify before relying on them

Most container tags below were checked live against their registries while building this
pipeline; a couple have no clean public option and need a decision before first use:

| Tool | Tag used | Confidence |
|---|---|---|
| GATK4 | `quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1` | verified live, Sep 2025 |
| Arriba | `quay.io/biocontainers/arriba:2.5.1--h87b9561_0` | verified live, Aug 2025 |
| rMATS-turbo | `xinglab/rmats:v4.3.0` (Docker Hub, **not** biocontainers — no working biocontainers tag exists for rMATS-turbo) | verified live |
| Bioconductor base | `bioconductor/bioconductor_docker:3.23-r-4.6.0` | verified live |
| snpEff | `quay.io/biocontainers/snpeff:5.2--hdfd78af_1` | not independently verified — check before relying on it |

## chr21-subset test assets

These are generated and present in the repo:

- `assets/test_genome.fa` — chr21 only, from UCSC hg38 (`hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz`), 46,709,983bp, `chr21`-prefixed naming.
- `assets/test_annotation.gtf` — chr21 features filtered from GENCODE v44's full annotation (38,256 lines).
- `assets/test_known_sites.vcf.gz` (+ `.tbi`) — chr21 variants from Ensembl release 111's variation VCF, downsampled 1:150 (96,211 variants) and re-headered to `chr21` naming, for the BQSR test.
- `test_data/*.fastq.gz` — synthetic paired-end 100bp reads, generated by `assets/generate_test_fastqs.py` (deterministic, `seed=42`) by sampling fragments from spliced chr21 transcript sequences (so fragments genuinely span splice junctions, exercising STAR's splice-aware alignment and rMATS-turbo, unlike naive genomic-window sampling). Three signals are deliberately injected so each modality has something real to find rather than an empty result:
  - **DE**: 3 transcripts get extra reads only in `treated` samples (8,000 pairs/sample vs 5,600 for `control`).
  - **Splicing**: one 4-exon transcript has its middle exon skipped in `control` samples only — a textbook skipped-exon (SE) event for rMATS-turbo.
  - **Variant**: a SNP is injected into ~50% of reads from one transcript, `treated` samples only — gives HaplotypeCaller something to actually call.

Regenerate from scratch any time with:
```bash
curl -sL -o /tmp/chr21.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz
gunzip -c /tmp/chr21.fa.gz > assets/test_genome.fa

curl -sL -o /tmp/gencode.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip -c /tmp/gencode.gtf.gz | awk -F'\t' '$1=="chr21"' > assets/test_annotation.gtf

curl -sL -o /tmp/chr21_var.vcf.gz "https://ftp.ensembl.org/pub/release-111/variation/vcf/homo_sapiens/homo_sapiens-chr21.vcf.gz"
awk -v OFS='\t' '
  /^##fileformat/ { print; print "##contig=<ID=chr21,length=46709983>"; next }
  /^##contig/ { next }
  /^#/ { print; next }
  (NR%150)==0 { $1="chr21"; print }
' <(gunzip -c /tmp/chr21_var.vcf.gz) > /tmp/known_sites.vcf
# sort the variant body by position (the awk pass above preserves input order)
grep "^#" /tmp/known_sites.vcf > /tmp/known_sites.header.vcf
grep -v "^#" /tmp/known_sites.vcf | sort -k1,1 -k2,2n > /tmp/known_sites.body.vcf
cat /tmp/known_sites.header.vcf /tmp/known_sites.body.vcf > assets/test_known_sites.vcf
bgzip -f assets/test_known_sites.vcf && tabix -p vcf assets/test_known_sites.vcf.gz

python3 assets/generate_test_fastqs.py   # writes test_data/*.fastq.gz
```
Then: `nextflow run main.nf -profile docker,test` (or `singularity,test` / add the appropriate
profile for your environment).

**This has only been validated with `nextflow run main.nf -profile docker,test -preview`** (full
DAG/channel-wiring dry run, no Docker daemon was available in the environment these assets were
built in) — no container has actually been pulled and no tool has actually run against this data
yet. Run it for real (without `-preview`) on your Ubuntu workstation as the first real validation
step.

**Limitation**: a single-chromosome test set structurally cannot validate fusion detection (real
fusions are frequently inter-chromosomal) — a green test run proves pipeline mechanics for Arriba,
not its actual sensitivity/specificity. `conf/test.config` leaves `--run_fusion` off for this
reason; test Arriba separately against one of its own published test datasets if you want deeper
local validation.

## Note on testing

This pipeline's syntax, config-profile merging, and full channel wiring (including the GATK
7-process chain, the rMATS-turbo group-splitting logic, and the multi-consumer fan-out from one
shared alignment) were validated with `nextflow run main.nf -preview` across all five modalities —
that confirms the DAG is structurally sound, but **no actual bioinformatics tool has been run** (no
containers pulled, no real reference data available in the environment this was built in). Before
trusting results on a real cohort: run the chr21 smoke test above on your own workstation, then a
small real subset on KAIMRC Slurm and Google Batch, checking outputs at each new modality before
scaling up. Paste any Nextflow/tool error output and I'll help debug the relevant module.
