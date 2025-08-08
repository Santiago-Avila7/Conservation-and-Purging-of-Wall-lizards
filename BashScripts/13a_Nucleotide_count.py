#!/usr/bin/env python
"""
VCF Nucleotide Counter for Polarization Analysis

This script counts nucleotides per variant across specified sample groups and outputs counts in the format:
<ingroup_counts>\t<outgroup1_counts> <outgroup2_counts> <outgroup3_counts>

Key features:
- First column (ingroup) is separated by tab
- Outgroups are space-separated
- Output format: "A,C,G,T" counts
- Handles phased/unphased genotypes
- Skips non-SNPs and non-ACGT bases

Usage:
    python count_nucleotides.py input.vcf.gz ingroup.txt outgroup1.txt outgroup2.txt outgroup3.txt > output_counts.tsv
"""
import sys
import pysam  # type: ignore
from collections import defaultdict

def load_sample_groups(group_files):
    """Load sample names from group files with validation

    Args:
        group_files: List of 4 filenames [ingroup, outgroup1, outgroup2, outgroup3]
    Returns:
        List of sample lists in order: [ingroup, outgroup1, outgroup2, outgroup3]
    """
    groups = []
    group_names = ["ingroup", "outgroup1", "outgroup2", "outgroup3"]
    for idx, filename in enumerate(group_files):
        try:
            with open(filename) as f:
                samples = [l.strip() for l in f if l.strip()]
            if not samples:
                sys.stderr.write(f"Error: Empty group file - {filename}\n")
                sys.exit(1)
            groups.append(samples)
            sys.stderr.write(f"Loaded {len(samples)} samples for {group_names[idx]}\n")
        except FileNotFoundError:
            sys.stderr.write(f"Error: Group file not found - {filename}\n")
            sys.exit(1)
    return groups

def count_nucleotides(vcf_path, groups):
    """Count nucleotides per group per variant and output in requested format

    Args:
        vcf_path: Path to bgzipped VCF file
        groups: List of sample lists [ingroup, outgroup1, outgroup2, outgroup3]
    """
    # Open VCF and validate samples
    try:
        vcf = pysam.VariantFile(vcf_path)
    except (IOError, ValueError) as e:
        sys.stderr.write(f"Error: Could not open VCF file {vcf_path}: {e}\n")
        sys.exit(1)

    # Verify all samples exist in VCF
    all_samples = [s for grp in groups for s in grp]
    vcf_samples = set(vcf.header.samples)
    missing = [s for s in all_samples if s not in vcf_samples]
    if missing:
        sys.stderr.write(f"Error: {len(missing)} samples not found in VCF\n")
        for s in missing[:5]:
            sys.stderr.write(f"  - {s}\n")
        if len(missing) > 5:
            sys.stderr.write(f"  ... and {len(missing) - 5} more\n")
        sys.exit(1)

    sys.stderr.write(f"Processing {vcf_path} with {len(all_samples)} samples\n")

    variant_count = 0
    for record in vcf:
        variant_count += 1
        if variant_count % 10000 == 0:
            sys.stderr.write(f"Processed {variant_count} variants\n")

        # STEP 1: Filter non-SNP and non-ACGT variants
        ref = record.ref.upper()
        if len(ref) != 1 or ref not in 'ACGT':
            continue
        alts = tuple(alt.upper() for alt in record.alts or ())
        if any(len(a) != 1 or a not in 'ACGT' for a in alts):
            continue

        # STEP 2: Initialize nucleotide counters for all groups
        group_counts = [defaultdict(int) for _ in groups]

        # STEP 3: Process genotypes for each group
        for gi, sample_list in enumerate(groups):
            for sample in sample_list:
                gt = record.samples[sample].get('GT')
                if not gt or any(allele is None for allele in gt):
                    continue
                for allele in gt:
                    if allele == 0:
                        base = ref
                    elif allele == 1 and alts:
                        base = alts[0]
                    else:
                        continue  # skip missing or higher-order ALT
                    if base in 'ACGT':
                        group_counts[gi][base] += 1

        # STEP 4: Format output
        # ingroup: comma-separated, then tab
        ing = group_counts[0]
        ing_str = f"{ing['A']},{ing['C']},{ing['G']},{ing['T']}"
        # outgroups: space-separated
        out_strs = []
        for cnts in group_counts[1:]:
            out_strs.append(f"{cnts['A']},{cnts['C']},{cnts['G']},{cnts['T']}")
        print(f"{ing_str}\t{' '.join(out_strs)}")

    sys.stderr.write(f"Completed. Processed {variant_count} total variants\n")

if __name__ == '__main__':
    if len(sys.argv) != 6:
        sys.stderr.write(
            "Usage: python count_nucleotides.py <input.vcf.gz> "
            "<ingroup.txt> <outgroup1.txt> <outgroup2.txt> <outgroup3.txt>\n"
        )
        sys.stderr.write(
            "Output format: <ingroup_counts>\t<outgroup1_counts> "
            "<outgroup2_counts> <outgroup3_counts>\n"
        )
        sys.exit(1)

    group_files = sys.argv[2:]
    groups = load_sample_groups(group_files)
    count_nucleotides(sys.argv[1], groups)
