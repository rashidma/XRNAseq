# assets/

Reference files used for local/test runs are **not tracked in git** (too large).

## Files needed for `-profile test`

| File | Source |
|------|--------|
| `test_genome.fa` | Generate via `python generate_test_fastqs.py` or download a chr21 slice from UCSC/Ensembl |
| `test_genome.fa.fai` | `samtools faidx test_genome.fa` |
| `test_annotation.gtf` | Matching Ensembl GTF (chr21 slice or full GRCh38) |
| `test_known_sites.vcf.gz` | dbSNP / Mills subset for chr21; `bgzip` + `tabix -p vcf` |
| `test_known_sites.vcf.gz.tbi` | Created by `tabix` above |

For real cohort runs, point `params.genome`, `params.gtf`, and `params.known_sites` in
`conf/base.config` (or on the CLI) to the full reference paths on your system.
