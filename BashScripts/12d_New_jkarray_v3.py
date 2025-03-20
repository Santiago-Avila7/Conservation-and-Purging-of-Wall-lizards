#!/usr/bin/env python3
"""
Jackknife Analysis Script for Italian-Only VCF

This script performs jackknife analysis using a pre-generated set of 100k intergenic sites.
It calculates allele counts for different variant impacts (HIGH, MODERATE, LOW, MODIFIER)
and intergenic sites across genomic blocks.
"""

import sys
import pysam
from collections import defaultdict

# Configuration
JACKKNIFE_ID = int(sys.argv[1])  # SLURM_ARRAY_TASK_ID
VCF_PATH = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Origin/All_Italian_annotated.vcf.gz"  # Updated VCF path
FIXED_INTERGENIC_FILE = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/4th/Italian_intergenic_sites.txt"
BLOCK_FILE = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/renamed_blocks.bed"
SAMPLE_LISTS = {
    'introduced': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Introduced_Purging_Test",
    'native': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Italian_Native"
}
OUTPUT_DIR = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/4th"

# ========================
# INITIALIZATION SECTION
# ========================
# Load genomic blocks
with open(BLOCK_FILE) as f:
    blocks = [line.strip().split('\t') for line in f]
current_block = blocks[JACKKNIFE_ID-1]
exclude_chrom, exclude_start = current_block[0], int(current_block[1])+1
exclude_end = int(current_block[2])

# Load fixed intergenic sites
with open(FIXED_INTERGENIC_FILE) as f:
    fixed_intergenic = set(tuple(line.strip().split('\t')) for line in f)

# Load sample lists
def load_samples(path):
    with open(path) as f:
        return set(s.strip() for s in f)

samples = {k: load_samples(v) for k, v in SAMPLE_LISTS.items()}

# ========================
# COUNTER INITIALIZATION
# ========================
counters = {
    'introduced': defaultdict(int),
    'native': defaultdict(int)
}

# Initialize counter structure
for group in ['introduced', 'native']:
    counters[group].update({
        'HIGH': 0, 'MODERATE': 0, 'LOW': 0, 'MODIFIER': 0,
        'Total_Alleles': 0, 'Intergenic_Alternate': 0, 'Intergenic_Total': 0
    })

# ========================
# CORE PROCESSING FUNCTIONS
# ========================
def parse_ann(ann_entry):
    """
    Parse SnpEff ANN field to extract impact and intergenic status
    Example ANN format: 'C|intergenic_region|MODIFIER|...'
    """
    for effect in ann_entry.split(','):
        parts = effect.split('|')
        if len(parts) > 2:
            effect_type = parts[1].lower()
            impact = parts[2].upper()
            is_intergenic = 'intergenic' in effect_type
            return impact, is_intergenic
    return 'MODIFIER', False

def process_genotype(gt):
    """
    Process genotype string to count alternate alleles
    Handles both phased and unphased genotypes:
    - 0/0 → 0 alts
    - 0/1 → 1 alt
    - 1/1 → 2 alts
    - ./. → Missing (returns 0)
    """
    if None in gt:  # Missing genotype (./.)
        return 0
    return sum(1 for allele in gt if allele == 1)

# ========================
# MAIN PROCESSING LOGIC
# ========================
with pysam.VariantFile(VCF_PATH) as vcf:
    for record in vcf:
        chrom, pos = record.chrom, record.start + 1  # VCF is 1-based
        
        # Skip variants in excluded block
        if chrom == exclude_chrom and exclude_start <= pos <= exclude_end:
            continue
        
        # Determine if variant is in fixed intergenic set
        is_fixed_intergenic = (chrom, str(pos)) in fixed_intergenic
        
        # Parse variant information
        ann_entry = record.info.get('ANN', [''])[0]
        impact, is_intergenic = parse_ann(ann_entry)
        
        # Process both groups in a single pass
        for group in ['introduced', 'native']:
            group_alts = 0
            intergenic_alts = 0
            intergenic_total = 0
            
            for sample in samples[group]:
                if sample not in record.samples:
                    continue
                
                # Get genotype and count alts
                gt = record.samples[sample].get('GT', (None, None))
                alt_count = process_genotype(gt)
                
                # Update counters
                group_alts += alt_count
                
                # Track intergenic counts for fixed sites
                if is_fixed_intergenic and is_intergenic:
                    intergenic_alts += alt_count
                    intergenic_total += 2  # Diploid total
                    
            # Update impact counters
            counters[group][impact] += group_alts
            counters[group]['Total_Alleles'] += group_alts
            
            # Update intergenic counters
            if is_fixed_intergenic and is_intergenic:
                counters[group]['Intergenic_Alternate'] += intergenic_alts
                counters[group]['Intergenic_Total'] += intergenic_total

# ========================
# OUTPUT GENERATION
# ========================
output_path = f"{OUTPUT_DIR}/jackknife_{JACKKNIFE_ID}_results.tsv"
with open(output_path, 'w') as f:
    # Write header
    header = "Origin\tHIGH\tMODERATE\tLOW\tMODIFIER\tTotal_Alleles\tIntergenic_Alternate\tIntergenic_Total\n"
    f.write(header)
    
    # Write data
    for group in ['introduced', 'native']:
        line = "\t".join([
            group,
            str(counters[group]['HIGH']),
            str(counters[group]['MODERATE']),
            str(counters[group]['LOW']),
            str(counters[group]['MODIFIER']),
            str(counters[group]['Total_Alleles']),
            str(counters[group]['Intergenic_Alternate']),
            str(counters[group]['Intergenic_Total'])
        ]) + "\n"
        f.write(line)