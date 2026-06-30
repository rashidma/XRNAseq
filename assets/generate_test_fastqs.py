#!/usr/bin/env python3
"""
Generates synthetic paired-end RNA-seq FASTQs for the chr21 smoke test
(assets/test_samplesheet.csv). Not a substitute for real reads -- this exists
to exercise pipeline mechanics with a real, interpretable signal in each
branch, not to be biologically realistic.

Reads are sampled from spliced transcript sequences built from
assets/test_annotation.gtf exons over assets/test_genome.fa, so fragments
genuinely span splice junctions (exercising STAR's splice-aware alignment and
rMATS-turbo). Three deliberate signals are injected so each modality has
something real to find, not just an empty/trivial result:

  - DE:       a handful of transcripts get extra reads only in "treated" samples
  - Splicing: one transcript has its middle exon skipped in "control" samples
              only (a textbook skipped-exon event for rMATS-turbo)
  - Variant:  a heterozygous-style SNP is injected into ~50% of reads from one
              transcript, "treated" samples only (for GATK HaplotypeCaller)

Run from the repo root: python3 assets/generate_test_fastqs.py
"""
import gzip
import os
import random
import re
from collections import defaultdict

random.seed(42)

GENOME_FA = "assets/test_genome.fa"
GTF = "assets/test_annotation.gtf"
OUT_DIR = "test_data"
READ_LEN = 100
FRAG_LEN = 300

SAMPLES = {
    "test_treated_1": "treated",
    "test_treated_2": "treated",
    "test_control_1": "control",
    "test_control_2": "control",
}
N_BACKGROUND_PAIRS = 4000
N_SIGNAL_PAIRS = 800

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def load_chr21(path):
    chunks = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                chunks.append(line.rstrip())
    return "".join(chunks).upper()


def parse_exons(path):
    transcripts = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            start, end, strand, attrs = int(cols[3]), int(cols[4]), cols[6], cols[8]
            m_t = re.search(r'transcript_id "([^"]+)"', attrs)
            if not m_t:
                continue
            transcripts[m_t.group(1)].append([start, end, strand])
    for tid in transcripts:
        transcripts[tid].sort(key=lambda e: e[0])
    return transcripts


def build_transcript_seq(chr21, exons, skip_exon_index=None):
    strand = exons[0][2]
    parts = []
    for i, (start, end, _strand) in enumerate(exons):
        if skip_exon_index is not None and i == skip_exon_index:
            continue
        parts.append(chr21[start - 1:end])  # GTF is 1-based inclusive
    seq = "".join(parts)
    return revcomp(seq) if strand == "-" else seq


def inject_snp(seq, pos):
    alt = {"A": "T", "T": "A", "C": "G", "G": "C"}.get(seq[pos], "N")
    return seq[:pos] + alt + seq[pos + 1:]


def sample_fragment(seq):
    if len(seq) < FRAG_LEN:
        return None
    start = random.randint(0, len(seq) - FRAG_LEN)
    frag = seq[start:start + FRAG_LEN]
    return frag[:READ_LEN], revcomp(frag[-READ_LEN:])


def write_fastq_pair(sample, r1_records, r2_records):
    os.makedirs(OUT_DIR, exist_ok=True)
    for suffix, records in (("R1", r1_records), ("R2", r2_records)):
        path = f"{OUT_DIR}/{sample}_{suffix}.fastq.gz"
        with gzip.open(path, "wt") as fh:
            for i, seq in enumerate(records):
                qual = "I" * len(seq)
                mate = 1 if suffix == "R1" else 2
                fh.write(f"@{sample}:{i}/{mate}\n{seq}\n+\n{qual}\n")


def main():
    chr21 = load_chr21(GENOME_FA)
    transcripts = parse_exons(GTF)

    candidates = []
    for tid, exons in transcripts.items():
        if len(exons) < 3:
            continue
        total_len = sum(e[1] - e[0] + 1 for e in exons)
        if total_len < FRAG_LEN + 50:
            continue
        candidates.append(tid)
    candidates.sort()  # deterministic ordering before shuffling
    random.shuffle(candidates)
    if len(candidates) < 65:
        raise SystemExit(f"Only found {len(candidates)} usable multi-exon transcripts on chr21 -- expected >=65.")

    expr_transcripts = candidates[:60]
    splicing_tid = candidates[60]
    de_up_tids = candidates[61:64]
    variant_tid = candidates[64]

    print(f"{len(transcripts)} transcripts parsed, {len(candidates)} usable for sampling")
    print(f"splicing transcript: {splicing_tid} ({len(transcripts[splicing_tid])} exons)")
    print(f"DE-upregulated transcripts (treated only): {de_up_tids}")
    print(f"variant transcript (treated only, ~50% of reads): {variant_tid}")

    splicing_exons = transcripts[splicing_tid]
    skip_idx = len(splicing_exons) // 2
    seq_with_exon = build_transcript_seq(chr21, splicing_exons)
    seq_skip_exon = build_transcript_seq(chr21, splicing_exons, skip_exon_index=skip_idx)

    variant_seq_ref = build_transcript_seq(chr21, transcripts[variant_tid])
    variant_pos = len(variant_seq_ref) // 2
    variant_seq_alt = inject_snp(variant_seq_ref, variant_pos)

    for sample, condition in SAMPLES.items():
        r1s, r2s = [], []

        # Background expression: same transcript pool/rate across all samples.
        for _ in range(N_BACKGROUND_PAIRS):
            tid = random.choice(expr_transcripts)
            seq = build_transcript_seq(chr21, transcripts[tid])
            frag = sample_fragment(seq)
            if frag:
                r1s.append(frag[0]); r2s.append(frag[1])

        # DE signal: extra reads from de_up_tids, treated only.
        if condition == "treated":
            for tid in de_up_tids:
                seq = build_transcript_seq(chr21, transcripts[tid])
                for _ in range(N_SIGNAL_PAIRS):
                    frag = sample_fragment(seq)
                    if frag:
                        r1s.append(frag[0]); r2s.append(frag[1])

        # Splicing signal: control skips the middle exon, treated includes it.
        splice_seq = seq_with_exon if condition == "treated" else seq_skip_exon
        for _ in range(N_SIGNAL_PAIRS):
            frag = sample_fragment(splice_seq)
            if frag:
                r1s.append(frag[0]); r2s.append(frag[1])

        # Variant signal: heterozygous-style SNP in ~50% of treated reads only.
        for _ in range(N_SIGNAL_PAIRS):
            use_alt = condition == "treated" and random.random() < 0.5
            seq = variant_seq_alt if use_alt else variant_seq_ref
            frag = sample_fragment(seq)
            if frag:
                r1s.append(frag[0]); r2s.append(frag[1])

        combined = list(zip(r1s, r2s))
        random.shuffle(combined)
        r1s, r2s = zip(*combined) if combined else ([], [])

        write_fastq_pair(sample, list(r1s), list(r2s))
        print(f"{sample} ({condition}): {len(r1s)} pairs")


if __name__ == "__main__":
    main()
