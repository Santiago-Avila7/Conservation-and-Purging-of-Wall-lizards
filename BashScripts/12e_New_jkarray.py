#!/usr/bin/env python3

import sys
import gzip
import os

# Configuration
JACKKNIFE_ID = int(sys.argv[1])  # SLURM_ARRAY_TASK_ID passed as argument
VCF_PATH = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/All_Origins_annotated.vcf.gz"
BLOCK_FILE = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/renamed_blocks.bed"
ITALIAN_SAMPLES = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Italian"
INTRODUCED_SAMPLES = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Italian_Introduced"
NATIVE_SAMPLES = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Italian_Native"
OUTPUT_DIR = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife2"

# Load genomic blocks
with open(BLOCK_FILE) as f:
    blocks = [line.strip().split('\t') for line in f]
current_block = blocks[JACKKNIFE_ID - 1]
exclude_region = f"{current_block[0]}:{int(current_block[1]) + 1}-{current_block[2]}"

# Load sample lists
def load_samples(filename):
    with open(filename) as f:
        return set(s.strip() for s in f)

italian_samples = load_samples(ITALIAN_SAMPLES)
introduced = load_samples(INTRODUCED_SAMPLES)
native = load_samples(NATIVE_SAMPLES)

# Initialize allele count dictionary
counters = {
    'Introduced': {'HIGH': 0, 'MODERATE': 0, 'LOW': 0, 'MODIFIER': 0, 'Total_Alleles': 0},
    'Native': {'HIGH': 0, 'MODERATE': 0, 'LOW': 0, 'MODIFIER': 0, 'Total_Alleles': 0}
}

def get_impact(info):
    """Parse SnpEff impact from INFO field."""
    for field in info.split(';'):
        if field.startswith('ANN='):
            ann = field.split('|')[2]  # Extract impact type
            return ann.split('&')[0].upper()  # Use the first impact if multiple
    return 'MODIFIER'  # Default to MODIFIER if unknown

def process_variant(line, samples, origin):
    """Count alleles for a variant in specified samples."""
    fields = line.strip().split('\t')
    chrom, pos, _, ref, alt, _, filt, info = fields[:8]
    genotypes = fields[9:]

    if filt != 'PASS':  # Skip non-PASS variants
        return

    impact = get_impact(info)
    total_alleles = 0

    for i, gt_info in enumerate(genotypes):
        sample = header[i + 9]  # Extract sample name from VCF header
        genotype = gt_info.split(":")[0]

        if sample not in samples:
            continue

        if genotype.startswith("."):  # Skip missing genotypes
            continue

        alleles = genotype.replace("|", "/").split("/")
        for a in alleles:
            if a == "1":
                total_alleles += 1

    # Update allele counts
    counters[origin][impact] += total_alleles
    counters[origin]['Total_Alleles'] += total_alleles

# Process VCF while skipping the excluded region
with gzip.open(VCF_PATH, 'rt') as vcf:
    for line in vcf:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')  # Extract sample headers
            continue

        chrom, pos = line.split('\t')[:2]
        if f"{chrom}:{pos}" in exclude_region:
            continue

        process_variant(line, introduced, 'Introduced')
        process_variant(line, native, 'Native')

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Save results
for origin in ['Introduced', 'Native']:
    output_file = f"{OUTPUT_DIR}/jackknife_{JACKKNIFE_ID}_{origin}.tsv"
    with open(output_file, 'w') as f:
        f.write("\t".join([
            origin,
            str(counters[origin]['HIGH']),
            str(counters[origin]['MODERATE']),
            str(counters[origin]['LOW']),
            str(counters[origin]['MODIFIER']),
            str(counters[origin]['Total_Alleles']),
        ]) + "\n")
