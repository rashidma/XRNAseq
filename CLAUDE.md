# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A hand-rolled (no nf-core) Nextflow DSL2 RNA-seq pipeline that runs unmodified on a local Ubuntu
workstation, KAIMRC's Slurm HPC, and Google Cloud Batch. One shared alignment per sample fans out
to five analysis modalities: differential expression (featureCounts → DESeq2 → fgsea), RNA-seq
variant calling (GATK), fusion detection (Arriba), and alternative splicing (rMATS-turbo).

**Status: the `docker,test` profile (chr21 synthetic data) has been run for real end-to-end** —
DGE/DESeq2, GATK variant calling, and rMATS-turbo splicing all completed and produced genuine
results (significant DE genes, significant SE/MXE/etc. splicing events, filtered VCFs). Fusion
(Arriba) and GSEA are off in `conf/test.config` and remain unvalidated. Slurm and Google Batch
profiles remain unvalidated on real infrastructure — only `-preview`'d there.

## Commands

```bash
# Full test run (chr21 synthetic data) -- the primary way to validate any change
nextflow run main.nf -profile docker,test          # or singularity,test on Slurm-like systems

# Dry-run: validate DAG/channel wiring without executing anything or needing Docker
nextflow run main.nf -profile docker,test -preview

# Resume after a crash/edit (Nextflow caches completed tasks)
nextflow run main.nf -profile docker,test -resume

# Real run
nextflow run main.nf --samplesheet s.csv --genome_fasta ref.fa --gtf ann.gtf --outdir results -profile docker,local
```

There is no unit test suite, linter, or build step — this is a workflow pipeline, not an
application. "Testing a change" means running (or at minimum `-preview`-ing) the `test` profile
against the chr21 synthetic assets in `assets/` and `test_data/`, and checking the relevant
`results/<stage>/` output. If those assets are missing, regenerate them with the commands in
README.md ("chr21-subset test assets" section) or `WORKSTATION_SETUP.txt`.

Before the DGE/GSEA branch can run for real, `containers/r_analysis/Dockerfile` must be built and
its tag substituted for the `PLACEHOLDER_VERIFY_AT_IMPLEMENTATION:xrnaseq-r-analysis` placeholder
in `modules/deseq2.nf` and `modules/gsea_fgsea.nf`.

## Architecture: one alignment, four independent consumers

`STAR_ALIGN_2PASS` (`modules/star_align_2pass.nf`) runs **once per sample** with
`--twopassMode Basic` and Arriba's `--chim*` flags always on, so chimeric reads land in the same
coordinate-sorted BAM. That BAM channel (`ALIGNMENT.out.bam_indexed` in `main.nf`) is then consumed
independently, with no re-alignment, by up to four branches gated by `--run_*` params:

1. `QUANTIFICATION` → `DIFFERENTIAL_EXPRESSION` — featureCounts → DESeq2 → fgsea (`--run_dge` /
   `--run_gsea`, **default on**)
2. `VARIANT_CALLING` — GATK MarkDuplicates → SplitNCigarReads → BQSR → HaplotypeCaller →
   VariantFiltration → snpEff (`--run_variant_calling`, default off — most expensive branch)
3. `FUSION` — Arriba (`--run_fusion`, default off)
4. `SPLICING` — rMATS-turbo (`--run_splicing`, default off)

`main.nf` is intentionally thin: it validates params (`checkParams()`) and wires subworkflows
together with plain `if (params.run_x) { ... }` blocks — no logic beyond wiring belongs there. Each
subworkflow in `subworkflows/` composes one or more single-process modules from `modules/` (one
Nextflow process per file, named after the tool/step). When adding a new modality, follow this same
three-layer split (module → subworkflow → wiring in `main.nf`) rather than inlining a process
directly into a subworkflow file.

## Configuration layering

- `nextflow.config` — all `params` defaults live here centrally (not scattered in modules) plus
  profile definitions. Should not need to change when moving between environments.
- `conf/base.config` — resource labels (`process_low/medium/high/high_memory`) applied via
  `label 'process_x'` on each module process. To retune CPU/memory/time for one environment, edit
  the environment's config (`conf/local.config`, `conf/slurm.config`, `conf/google_batch.config`,
  `conf/test.config`), not `base.config` or the module itself.
- `conf/slurm.config` / `conf/google_batch.config` — contain `TODO` placeholders (partition,
  account, cache dirs / GCP project, region, bucket) that must be filled in before first use on
  those backends.
- `conf/test.config` — chr21-subset smoke test wiring; scales resource labels down and points at
  the synthetic assets. `run_fusion` is deliberately left off here — a single chromosome can't
  meaningfully exercise fusion detection.
- Profiles compose via comma-separation, e.g. `-profile docker,slurm` — but `docker,slurm` is not
  useful since `slurm.config` forces Singularity regardless (Docker generally isn't permitted on
  shared HPC nodes).

## Conventions to follow when editing modules/subworkflows

- One process per module file; module filename matches the tool (`arriba.nf`, `featurecounts.nf`).
- Every process declares `container '<pinned tag>'` and `label 'process_x'` — never hardcode
  cpus/memory in a module. See the README's "Container tags" table before changing or adding a
  container tag; several were verified live against registries and one (snpEff) is flagged as
  unverified.
- `--outSAMattrRGline` on STAR and similar cross-cutting flags exist because a *later*, seemingly
  unrelated branch (GATK) needs them — read the comment above a script block before removing a flag
  that looks unused by the module's own immediate output.
- `differential_expression.nf`'s pattern for an optional sub-step that needs an unset param (GSEA
  needing `--msigdb_gmt`) is to `log.warn` and skip via `Channel.empty()`, not to hard-fail the
  whole run — follow this pattern for any new optional dependency rather than adding a new required
  param.
- R driver scripts (`bin/run_deseq2.R`, `bin/run_fgsea.R`) are called by their corresponding module
  (`modules/deseq2.nf`, `modules/gsea_fgsea.nf`) inside the custom `xrnaseq-r-analysis` container —
  keep script CLI args in sync with the `script:` block in those modules if you change either side.
  **Gotcha**: `containers/r_analysis/Dockerfile` also `COPY`s these same scripts to
  `/usr/local/bin/` at image build time, and Nextflow appends (not prepends) the repo's `bin/` dir
  to `PATH` inside the container — so the container's baked-in copy wins over an edited repo script
  until you rebuild the image (`docker build -t xrnaseq-r-analysis:1.0.0 -f containers/r_analysis/Dockerfile .`).
  Editing `bin/run_deseq2.R` or `bin/run_fgsea.R` without rebuilding is a silent no-op.
- `rmats.py` is not on `PATH` in `xinglab/rmats:v4.3.0` — `modules/rmats_turbo.nf` invokes it as
  `python3 /rmats/rmats.py` (confirmed via `docker run --entrypoint sh xinglab/rmats:v4.3.0 -c 'find / -iname rmats.py'`).
- `bin/run_deseq2.R`'s `vst(dds, blind = TRUE)` call falls back to
  `varianceStabilizingTransformation()` on error, since `vst()`'s default subsampling needs ~1000
  genes and fails on small inputs (the chr21 test set, or any small gene panel) — this is DESeq2's
  own documented remedy, not a workaround.

## Known caveats worth preserving in behavior/docs

- RNA-seq variant calls use a lighter hard filter (`FS > 30 / QD < 2 / SnpCluster`) than DNA-seq
  GATK germline pipelines, per GATK's own guidance that RNA-seq lacks VQSR-quality validation data.
  Don't silently upgrade this to a stricter DNA-seq-style filter without discussing — it's a
  deliberate choice, not an oversight.
- rMATS-turbo has no batch-covariate mechanism, unlike DESeq2's design formula — multi-site/batch
  cohorts have no built-in confounding defense on the splicing branch.
- `--read_length` must match actual FASTQ read length exactly (rMATS-turbo assumes uniform length).
- A single-chromosome test run structurally cannot validate Arriba's real sensitivity (fusions are
  often inter-chromosomal) — a green `test` profile run only proves pipeline mechanics for that
  branch.
