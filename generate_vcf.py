#!/usr/bin/env python3
"""
Synthetic Genomic Pipeline: FASTQ -> BAM -> VCF
================================================
Generates synthetic paired-end reads from a mock reference genome,
aligns them, and calls variants to produce a valid VCF file for
one simulated patient.

Required external tools:
    - bwa        (alignment)
    - samtools   (SAM/BAM manipulation)
    - bcftools   (variant calling)

Install on macOS:
    brew install bwa samtools bcftools

Install on Ubuntu/Debian:
    sudo apt-get install bwa samtools bcftools

Usage:
    python scripts/generate_vcf.py                     # defaults
    python scripts/generate_vcf.py --outdir output      # custom output dir
    python scripts/generate_vcf.py --num-variants 30    # more variants
"""

from __future__ import annotations

import argparse
import os
import random
import shutil
import subprocess
import sys
import tempfile
import textwrap


# ---------------------------------------------------------------------------
# 1. Reference Genome Generation
# ---------------------------------------------------------------------------

def generate_reference_genome(length: int, seed: int = 42) -> str:
    """Return a random DNA sequence of *length* bases (deterministic via seed)."""
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=length))


def write_reference_fasta(sequence: str, path: str, chrom: str = "chr_synth") -> None:
    """Write a single-contig FASTA file."""
    with open(path, "w") as fh:
        fh.write(f">{chrom}\n")
        # Wrap to 80 characters per line (standard FASTA convention)
        for i in range(0, len(sequence), 80):
            fh.write(sequence[i : i + 80] + "\n")
    print(f"[ref]  Wrote reference ({len(sequence)} bp) -> {path}")


# ---------------------------------------------------------------------------
# 2. Variant Introduction & Read Simulation
# ---------------------------------------------------------------------------

def introduce_variants(
    reference: str,
    num_variants: int,
    seed: int = 99,
) -> tuple[str, str, list[dict]]:
    """
    Introduce random SNPs and small indels into *reference*.

    Returns two haplotype sequences (hap1, hap2) and a variant list.
    For 0/1 (het) variants: hap1 = reference allele, hap2 = alt allele.
    For 1/1 (hom-alt) variants: both haplotypes carry the alt allele.

    Each variant record has: pos (0-based), ref_allele, alt_allele, type, genotype
    """
    rng = random.Random(seed)
    variants: list[dict] = []

    # Eligible positions (leave a margin at the edges)
    margin = 200
    eligible = list(range(margin, len(reference) - margin))
    rng.shuffle(eligible)

    placed = 0
    used_positions: set[int] = set()

    for pos in eligible:
        if placed >= num_variants:
            break
        # Skip if too close to an already-placed variant
        if any(abs(pos - u) < 10 for u in used_positions):
            continue

        ref_base = reference[pos]
        variant_type = rng.choices(["snp", "ins", "del"], weights=[0.55, 0.25, 0.20])[0]
        genotype = "0/1"  # all variants are heterozygous

        if variant_type == "snp":
            alt_base = rng.choice([b for b in "ACGT" if b != ref_base])
            variants.append(dict(
                pos=pos, ref_allele=ref_base, alt_allele=alt_base,
                type="snp", genotype=genotype,
            ))

        elif variant_type == "ins":
            insert_seq = "".join(rng.choices("ACGT", k=rng.randint(1, 3)))
            variants.append(dict(
                pos=pos, ref_allele=ref_base, alt_allele=ref_base + insert_seq,
                type="ins", genotype=genotype,
            ))

        elif variant_type == "del":
            del_len = rng.randint(1, 3)
            if pos + del_len + 1 >= len(reference):
                continue
            deleted = reference[pos : pos + del_len + 1]
            variants.append(dict(
                pos=pos, ref_allele=deleted, alt_allele=ref_base,
                type="del", genotype=genotype,
            ))

        used_positions.add(pos)
        placed += 1

    # Sort by position
    variants.sort(key=lambda v: v["pos"])

    # Build haplotype sequences by applying variants from right to left
    # (reverse order avoids coordinate shifts from indels)
    hap1_list = list(reference)
    hap2_list = list(reference)

    for v in reversed(variants):
        p = v["pos"]
        ref_len = len(v["ref_allele"])
        alt_seq = v["alt_allele"]

        if v["genotype"] == "1/1":
            # Both haplotypes get the alt allele
            hap1_list[p : p + ref_len] = list(alt_seq)
            hap2_list[p : p + ref_len] = list(alt_seq)
        else:
            # 0/1: only hap2 gets the alt allele
            hap2_list[p : p + ref_len] = list(alt_seq)

    hap1 = "".join(hap1_list)
    hap2 = "".join(hap2_list)

    print(f"[var]  Placed {len(variants)} variants "
          f"(SNPs: {sum(1 for v in variants if v['type']=='snp')}, "
          f"INS: {sum(1 for v in variants if v['type']=='ins')}, "
          f"DEL: {sum(1 for v in variants if v['type']=='del')})")
    gt_summary = {gt: sum(1 for v in variants if v["genotype"] == gt)
                  for gt in ("0/1", "1/1")}
    print(f"[var]  Genotypes: {gt_summary}")
    return hap1, hap2, variants


def _phred_char(q: int) -> str:
    """Convert a Phred quality score (int) to the corresponding ASCII character."""
    return chr(q + 33)


def simulate_paired_reads(
    hap1: str,
    hap2: str,
    num_pairs: int,
    read_length: int,
    fragment_size: int,
    error_rate: float,
    seed: int = 7,
) -> tuple[list[tuple[str, str, str]], list[tuple[str, str, str]]]:
    """
    Simulate paired-end reads from a diploid genome (hap1 + hap2).

    Each read pair is drawn from one of the two haplotypes at random (50/50),
    mimicking diploid sequencing. This naturally produces het (0/1) evidence
    when only hap2 carries the alt allele, and hom-alt (1/1) evidence when
    both haplotypes carry it.

    Returns two lists of (name, sequence, quality) tuples for R1 and R2.
    """
    rng = random.Random(seed)
    haplotypes = (hap1, hap2)

    r1_reads: list[tuple[str, str, str]] = []
    r2_reads: list[tuple[str, str, str]] = []

    # Fragment size standard deviation (~15% of mean) mimics real
    # Illumina library prep variation.  This ensures bwa computes a wide
    # enough insert-size distribution to flag indel-spanning reads as
    # properly paired (critical for bcftools to count them).
    frag_sd = max(int(fragment_size * 0.15), 20)

    for i in range(num_pairs):
        # Pick a random haplotype (diploid: 50/50)
        source = haplotypes[rng.randint(0, 1)]

        # Sample a fragment size from a normal distribution
        frag_len = int(rng.gauss(fragment_size, frag_sd))
        frag_len = max(read_length * 2 + 10, min(frag_len, len(source) - 1))

        # Pick a random starting position for the fragment
        max_start = len(source) - frag_len
        if max_start < 1:
            max_start = 1
        start = rng.randint(0, max_start)

        fragment = source[start : start + frag_len]
        if len(fragment) < read_length * 2:
            continue

        # Forward read (R1) and reverse-complement read (R2)
        r1_seq_list = list(fragment[:read_length])
        r2_seq_list = list(_reverse_complement(fragment[-read_length:]))

        # Introduce sequencing errors and generate quality scores
        r1_seq, r1_qual = _add_errors_and_quals(r1_seq_list, error_rate, rng)
        r2_seq, r2_qual = _add_errors_and_quals(r2_seq_list, error_rate, rng)

        name = f"synth_read_{i+1}"
        r1_reads.append((name, r1_seq, r1_qual))
        r2_reads.append((name, r2_seq, r2_qual))

    print(f"[sim]  Simulated {len(r1_reads)} read pairs "
          f"(read_len={read_length}, frag={fragment_size}, err={error_rate})")
    return r1_reads, r2_reads


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def _add_errors_and_quals(
    seq_list: list[str],
    error_rate: float,
    rng: random.Random,
) -> tuple[str, str]:
    """
    Mutate bases at *error_rate* probability and produce a realistic quality
    string. Bases near the 3' end (higher index) get slightly lower qualities
    to mimic Illumina degradation.
    """
    quals = []
    for idx in range(len(seq_list)):
        # Quality degrades toward the end of the read
        base_q = rng.gauss(35 - idx * 0.05, 3)
        base_q = max(2, min(40, int(base_q)))

        if rng.random() < error_rate:
            original = seq_list[idx]
            seq_list[idx] = rng.choice([b for b in "ACGT" if b != original])
            base_q = min(base_q, 15)  # errors get low quality

        quals.append(_phred_char(base_q))

    return "".join(seq_list), "".join(quals)


def write_fastq(reads: list[tuple[str, str, str]], path: str, read_num: int) -> None:
    """Write reads to a FASTQ file."""
    with open(path, "w") as fh:
        for name, seq, qual in reads:
            fh.write(f"@{name}/{read_num}\n{seq}\n+\n{qual}\n")
    print(f"[fq]   Wrote {len(reads)} reads -> {path}")


# ---------------------------------------------------------------------------
# 3. Alignment: FASTQ -> sorted, indexed BAM
# ---------------------------------------------------------------------------

def check_tools() -> None:
    """Verify that all required external tools are on PATH."""
    missing = []
    for tool in ("bwa", "samtools", "bcftools"):
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        print(
            f"ERROR: The following tools are not installed or not on PATH: "
            f"{', '.join(missing)}\n"
            f"Install with:\n"
            f"  macOS:  brew install {' '.join(missing)}\n"
            f"  Linux:  sudo apt-get install {' '.join(missing)}",
            file=sys.stderr,
        )
        sys.exit(1)


def index_reference(ref_path: str) -> None:
    """Create bwa index and samtools faidx for the reference genome."""
    subprocess.run(["bwa", "index", ref_path], check=True, capture_output=True)
    subprocess.run(["samtools", "faidx", ref_path], check=True, capture_output=True)
    print(f"[idx]  Indexed reference: {ref_path}")


def align_reads(
    ref_path: str,
    fq1_path: str,
    fq2_path: str,
    bam_path: str,
    sample_name: str = "SYNTH_PATIENT_01",
    threads: int = 2,
) -> None:
    """
    Align paired-end reads with bwa mem, pipe through samtools to produce
    a coordinate-sorted, indexed BAM file.

    Steps:
        bwa mem  ->  samtools sort  ->  .bam
        samtools index .bam          ->  .bam.bai
    """
    rg_tag = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tLB:lib1"

    # bwa mem | samtools sort -> BAM
    bwa_cmd = [
        "bwa", "mem",
        "-t", str(threads),
        "-R", rg_tag,
        ref_path, fq1_path, fq2_path,
    ]
    sort_cmd = [
        "samtools", "sort",
        "-@", str(threads),
        "-o", bam_path,
        "-",  # read from stdin
    ]

    print(f"[aln]  Running: bwa mem | samtools sort -> {bam_path}")
    bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    sort_proc = subprocess.run(sort_cmd, stdin=bwa_proc.stdout, capture_output=True)
    bwa_proc.stdout.close()
    bwa_proc.wait()

    if bwa_proc.returncode != 0 or sort_proc.returncode != 0:
        print("ERROR: Alignment or sorting failed.", file=sys.stderr)
        sys.exit(1)

    # Index the BAM
    subprocess.run(["samtools", "index", bam_path], check=True, capture_output=True)
    print(f"[aln]  Sorted & indexed BAM: {bam_path}")


# ---------------------------------------------------------------------------
# 4. Variant Calling: BAM -> VCF
# ---------------------------------------------------------------------------

def call_variants(ref_path: str, bam_path: str, vcf_path: str) -> None:
    """
    Call variants using bcftools mpileup + bcftools call.

    Produces a VCF with genotype (GT) fields including 0/0, 0/1, 1/1.
    """
    mpileup_cmd = [
        "bcftools", "mpileup",
        "-f", ref_path,
        "--max-depth", "1000",
        "-a", "FORMAT/AD,FORMAT/DP",   # annotate with allelic depth & depth
        "-q", "1",                      # min mapping quality
        "-Q", "10",                     # min base quality
        bam_path,
    ]
    call_cmd = [
        "bcftools", "call",
        "-mv",                          # multiallelic caller, output variants only
        "--ploidy", "2",
        "-Oz",                          # compressed VCF output
        "-o", vcf_path + ".gz",
    ]

    print(f"[vc]   Running: bcftools mpileup | bcftools call -> {vcf_path}")
    mpileup_proc = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    call_proc = subprocess.run(call_cmd, stdin=mpileup_proc.stdout, capture_output=True)
    mpileup_proc.stdout.close()
    mpileup_proc.wait()

    if mpileup_proc.returncode != 0 or call_proc.returncode != 0:
        print("ERROR: Variant calling failed.", file=sys.stderr)
        sys.exit(1)

    # Index, then decompress to plain VCF for readability
    subprocess.run(
        ["bcftools", "index", vcf_path + ".gz"],
        check=True, capture_output=True,
    )
    subprocess.run(
        ["bcftools", "view", vcf_path + ".gz", "-o", vcf_path],
        check=True, capture_output=True,
    )
    print(f"[vc]   VCF written: {vcf_path}")


# ---------------------------------------------------------------------------
# 5. Summary / QC
# ---------------------------------------------------------------------------

def summarize_vcf(vcf_path: str) -> None:
    """Print a quick summary of the VCF contents."""
    result = subprocess.run(
        ["bcftools", "stats", vcf_path],
        capture_output=True, text=True,
    )
    snps = indels = total = 0
    for line in result.stdout.splitlines():
        if line.startswith("SN") and "number of records:" in line:
            total = int(line.strip().split("\t")[-1])
        elif line.startswith("SN") and "number of SNPs:" in line:
            snps = int(line.strip().split("\t")[-1])
        elif line.startswith("SN") and "number of indels:" in line:
            indels = int(line.strip().split("\t")[-1])

    # Count genotypes
    gt_counts = {"0/0": 0, "0/1": 0, "1/1": 0, "other": 0}
    result_query = subprocess.run(
        ["bcftools", "query", "-f", "[%GT]\n", vcf_path],
        capture_output=True, text=True,
    )
    for gt in result_query.stdout.strip().splitlines():
        gt = gt.strip()
        if gt in gt_counts:
            gt_counts[gt] += 1
        else:
            gt_counts["other"] += 1

    print("\n" + "=" * 50)
    print("  PIPELINE SUMMARY")
    print("=" * 50)
    print(f"  Total variant records : {total}")
    print(f"  SNPs                  : {snps}")
    print(f"  Indels                : {indels}")
    print(f"  Genotypes             : {dict(gt_counts)}")
    print("=" * 50)


# ---------------------------------------------------------------------------
# Main Pipeline
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Synthetic genomic pipeline: FASTQ -> BAM -> VCF",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Example:
                python scripts/generate_vcf.py --outdir output --num-variants 20
        """),
    )
    p.add_argument("--outdir", default="pipeline_output",
                    help="Directory for all output files (default: pipeline_output)")
    p.add_argument("--ref-length", type=int, default=50_000,
                    help="Length of the synthetic reference genome in bp (default: 50000)")
    p.add_argument("--num-variants", type=int, default=15,
                    help="Number of variants to introduce (default: 15)")
    p.add_argument("--num-reads", type=int, default=10_000,
                    help="Number of read pairs to simulate (default: 10000)")
    p.add_argument("--read-length", type=int, default=150,
                    help="Read length in bp (default: 150)")
    p.add_argument("--fragment-size", type=int, default=400,
                    help="Mean fragment size in bp (default: 400)")
    p.add_argument("--error-rate", type=float, default=0.005,
                    help="Per-base sequencing error rate (default: 0.005)")
    p.add_argument("--seed", type=int, default=42,
                    help="Random seed for reproducibility (default: 42)")
    p.add_argument("--sample-name", default="SYNTH_PATIENT_01",
                    help="Sample name in BAM/VCF (default: SYNTH_PATIENT_01)")
    p.add_argument("--auto", action="store_true",
                    help="Generate 10 patients each for klinik1, klinik2, klinik3 (30 total)")
    p.add_argument("--joint", action="store_true",
                    help="Same as --auto, but also merge each clinic's VCFs into a joint "
                         "multi-sample VCF (klinikN.vcf) alongside the per-patient files")
    return p.parse_args()


def run_single_patient(
    args: argparse.Namespace,
    outdir: str,
    sample_name: str,
    seed: int,
    shared_ref: str | None = None,
) -> str:
    """
    Run the full pipeline for one patient, outputting all files to *outdir*.

    If *shared_ref* is provided, it is used as the reference sequence instead
    of generating a new one.  This is required when multiple patients must
    share the same reference (e.g. for bcftools merge in --joint mode).

    Returns the path to the output VCF.
    """
    os.makedirs(outdir, exist_ok=True)

    ref_path = os.path.join(outdir, "reference.fa")
    fq1_path = os.path.join(outdir, "reads_R1.fastq")
    fq2_path = os.path.join(outdir, "reads_R2.fastq")
    bam_path = os.path.join(outdir, "aligned.sorted.bam")
    vcf_path = os.path.join(outdir, "variants.vcf")

    if shared_ref is not None:
        ref_seq = shared_ref
    else:
        ref_seq = generate_reference_genome(args.ref_length, seed=seed)
    write_reference_fasta(ref_seq, ref_path)

    hap1, hap2, variants = introduce_variants(
        ref_seq, args.num_variants, seed=seed + 1,
    )

    r1, r2 = simulate_paired_reads(
        hap1=hap1, hap2=hap2,
        num_pairs=args.num_reads,
        read_length=args.read_length,
        fragment_size=args.fragment_size,
        error_rate=args.error_rate,
        seed=seed + 2,
    )
    write_fastq(r1, fq1_path, read_num=1)
    write_fastq(r2, fq2_path, read_num=2)

    index_reference(ref_path)
    align_reads(ref_path, fq1_path, fq2_path, bam_path, sample_name=sample_name)
    call_variants(ref_path, bam_path, vcf_path)

    return vcf_path


def merge_clinic_vcfs(clinic: str, clinic_num: int) -> None:
    """
    Merge all per-patient VCFs in a clinic into a single multi-sample VCF.

    Each per-patient VCF is compressed and indexed, then bcftools merge
    combines them into klinikN.vcf (placed next to raw_variants/).
    """
    vcf_dir = os.path.join(clinic, "raw_variants")
    patient_vcfs = sorted(
        [os.path.join(vcf_dir, f) for f in os.listdir(vcf_dir) if f.endswith(".vcf")],
    )

    # Compress and index each per-patient VCF for bcftools merge
    gz_paths = []
    for vcf in patient_vcfs:
        gz = vcf + ".gz"
        subprocess.run(["bcftools", "view", vcf, "-Oz", "-o", gz],
                        check=True, capture_output=True)
        subprocess.run(["bcftools", "index", gz],
                        check=True, capture_output=True)
        gz_paths.append(gz)

    # Merge into a joint multi-sample VCF
    joint_vcf = os.path.join(clinic, f"klinik{clinic_num}.vcf")
    merge_cmd = ["bcftools", "merge"] + gz_paths + ["-o", joint_vcf]
    subprocess.run(merge_cmd, check=True, capture_output=True)

    # Clean up per-patient .gz and .csi files (the plain .vcf files are kept)
    for gz in gz_paths:
        os.remove(gz)
        csi = gz + ".csi"
        if os.path.exists(csi):
            os.remove(csi)

    print(f"[merge] {joint_vcf}  ({len(patient_vcfs)} samples)")


def auto_generate(args: argparse.Namespace, joint: bool = False) -> None:
    """Generate 10 patients x 3 clinics. VCFs go to klinikN/raw_variants/."""
    clinics = ["klinik1", "klinik2", "klinik3"]
    patients_per_clinic = 10

    mode_label = "JOINT" if joint else "AUTO"
    print(f"\n{'='*50}")
    print(f"  {mode_label} MODE: {len(clinics)} clinics x {patients_per_clinic} patients")
    print(f"{'='*50}\n")

    for clinic_idx, clinic in enumerate(clinics):
        vcf_dir = os.path.join(clinic, "raw_variants")
        os.makedirs(vcf_dir, exist_ok=True)

        # In joint mode all patients in a clinic must share the same reference
        # so that bcftools merge can combine their VCFs.
        clinic_seed = args.seed + (clinic_idx * 1000)
        if joint:
            shared_ref = generate_reference_genome(args.ref_length, seed=clinic_seed)
        else:
            shared_ref = None

        for patient_num in range(1, patients_per_clinic + 1):
            sample_name = f"patient{patient_num}"
            seed = clinic_seed + (patient_num * 100)

            print(f"\n>>> [{clinic}] Generating {sample_name} (seed={seed})")

            with tempfile.TemporaryDirectory() as tmpdir:
                run_single_patient(args, tmpdir, sample_name, seed,
                                   shared_ref=shared_ref)

                src_vcf = os.path.join(tmpdir, "variants.vcf")
                dst_vcf = os.path.join(vcf_dir, f"{sample_name}.vcf")
                shutil.copy2(src_vcf, dst_vcf)

            print(f"    -> {dst_vcf}")

        # After all patients in this clinic are done, merge if joint mode
        if joint:
            clinic_num = clinic_idx + 1
            merge_clinic_vcfs(clinic, clinic_num)

    # Final summary
    print(f"\n{'='*50}")
    print(f"  {mode_label} MODE COMPLETE")
    print(f"{'='*50}")
    for clinic_idx, clinic in enumerate(clinics):
        vcf_dir = os.path.join(clinic, "raw_variants")
        count = len([f for f in os.listdir(vcf_dir) if f.endswith(".vcf")])
        print(f"  {vcf_dir}: {count} VCFs")
        if joint:
            joint_vcf = os.path.join(clinic, f"klinik{clinic_idx + 1}.vcf")
            print(f"  {joint_vcf} (joint multi-sample)")
    print(f"{'='*50}")


def main() -> None:
    args = parse_args()

    # Verify external tools are available
    check_tools()

    if args.joint:
        auto_generate(args, joint=True)
        return

    if args.auto:
        auto_generate(args)
        return

    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    ref_path = os.path.join(outdir, "reference.fa")
    fq1_path = os.path.join(outdir, "reads_R1.fastq")
    fq2_path = os.path.join(outdir, "reads_R2.fastq")
    bam_path = os.path.join(outdir, "aligned.sorted.bam")
    vcf_path = os.path.join(outdir, "variants.vcf")

    print("\n>>> Step 1: Generating synthetic reference genome")
    ref_seq = generate_reference_genome(args.ref_length, seed=args.seed)
    write_reference_fasta(ref_seq, ref_path)

    print("\n>>> Step 2: Introducing variants")
    hap1, hap2, variants = introduce_variants(
        ref_seq, args.num_variants, seed=args.seed + 1,
    )

    print("\n>>> Step 3: Simulating paired-end reads")
    r1, r2 = simulate_paired_reads(
        hap1=hap1, hap2=hap2,
        num_pairs=args.num_reads,
        read_length=args.read_length,
        fragment_size=args.fragment_size,
        error_rate=args.error_rate,
        seed=args.seed + 2,
    )
    write_fastq(r1, fq1_path, read_num=1)
    write_fastq(r2, fq2_path, read_num=2)

    print("\n>>> Step 4: Aligning reads (bwa mem -> samtools sort)")
    index_reference(ref_path)
    align_reads(ref_path, fq1_path, fq2_path, bam_path, sample_name=args.sample_name)

    print("\n>>> Step 5: Calling variants (bcftools mpileup | call)")
    call_variants(ref_path, bam_path, vcf_path)

    summarize_vcf(vcf_path)

    print(f"\nAll outputs in: {os.path.abspath(outdir)}/")
    print("Done.")


if __name__ == "__main__":
    main()
