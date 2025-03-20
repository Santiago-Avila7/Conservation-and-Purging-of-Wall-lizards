#!/usr/bin/env python3

import sys
import random
import pysam

# Configuration
JACKKNIFE_ID = int(sys.argv[1])  # SLURM_ARRAY_TASK_ID passed as argument
VCF_PATH = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/All_Origins_annotated.vcf.gz"
BLOCK_FILE = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/renamed_blocks.bed"
ITALIAN_SAMPLES = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Italian"
INTRODUCED_SAMPLES = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Introduced_Purging_Test"
NATIVE_SAMPLES = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Italian_Native"
OUTPUT_DIR = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife2"

# Load genomic blocks
with open(BLOCK_FILE) as f:
    blocks = [line.strip().split('\t') for line in f]
current_block = blocks[JACKKNIFE_ID-1]
exclude_chrom, exclude_start = current_block[0], int(current_block[1])+1
exclude_end = int(current_block[2])

# Load sample lists
with open(ITALIAN_SAMPLES) as f:
    italian_samples = set(s.strip() for s in f)

with open(INTRODUCED_SAMPLES) as f:
    introduced = set(s.strip() for s in f)

with open(NATIVE_SAMPLES) as f:
    native = set(s.strip() for s in f)

# Initialize counters
counters = {
    'Introduced': {'HIGH': 0, 'MODERATE': 0, 'LOW': 0, 
                   'MODIFIER': 0, 'Total_Alleles': 0,
                   'Intergenic_Alternate': 0, 'Intergenic_Total': 0},
    'Native': {'HIGH': 0, 'MODERATE': 0, 'LOW': 0,
               'MODIFIER': 0, 'Total_Alleles': 0,
               'Intergenic_Alternate': 0, 'Intergenic_Total': 0}
}

# Store intergenic variants for sampling
intergenic_variants = {
    'Introduced': {'alts': [], 'totals': []},
    'Native': {'alts': [], 'totals': []}
}

def parse_ann(ann_entry):
    """Parse SnpEff ANN field to get effect and impact"""
    effects = ann_entry.split(',')
    for effect in effects:
        parts = effect.split('|')
        if len(parts) > 2:
            effect_type = parts[1].lower()
            impact = parts[2].upper()
            is_intergenic = 'intergenic' in effect_type
            return is_intergenic, impact
    return False, 'MODIFIER'

def process_variant(record, samples, origin):
    """Count alleles for a variant in specified samples"""
    # Parse INFO field
    info = record.info
    ann_entry = info.get('ANN', [''])[0]
    is_intergenic, impact = parse_ann(ann_entry)
    
    # Count alleles
    total_alleles = 0
    intergenic_alts = 0
    intergenic_total = 0
    
    for sample in samples:
        if sample not in record.samples:
            continue
            
        gt = record.samples[sample]['GT']
        if None in gt:  # Skip missing genotypes
            continue
            
        # Count alternate alleles
        alleles = sum(1 for a in gt if a == 1)
        total_alleles += alleles * 2  # Diploid count
        
        # Track intergenic counts
        if is_intergenic:
            intergenic_alts += alleles
            intergenic_total += 2  # Diploid count
    
    # Update main counters
    counters[origin][impact] += total_alleles
    counters[origin]['Total_Alleles'] += total_alleles
    
    # Store intergenic data for sampling
    if is_intergenic:
        intergenic_variants[origin]['alts'].append(intergenic_alts)
        intergenic_variants[origin]['totals'].append(intergenic_total)

# Process VCF
with pysam.VariantFile(VCF_PATH) as vcf:
    # Filter samples
    vcf.subset_samples(list(italian_samples))
    
    # Iterate over variants, excluding the current block
    for record in vcf:
        chrom, pos = record.chrom, record.pos
        
        # Skip excluded region
        if chrom == exclude_chrom and exclude_start <= pos <= exclude_end:
            continue
        
        # Process Introduced samples
        process_variant(record, introduced, 'Introduced')
        
        # Process Native samples
        process_variant(record, native, 'Native')

# Process intergenic counts
for origin in ['Introduced', 'Native']:
    alts = intergenic_variants[origin]['alts']
    totals = intergenic_variants[origin]['totals']
    
    if not alts:
        continue
    
    # Sample 100k intergenic variants (with replacement if needed)
    n_samples = min(100000, len(alts))
    sampled_indices = random.choices(range(len(alts)), k=n_samples)
    
    sampled_alts = sum(alts[i] for i in sampled_indices)
    sampled_totals = sum(totals[i] for i in sampled_indices)
    
    counters[origin]['Intergenic_Alternate'] = sampled_alts
    counters[origin]['Intergenic_Total'] = sampled_totals

# Write results with headers
with open(f"{OUTPUT_DIR}/jackknife_{JACKKNIFE_ID}_results.tsv", 'w') as f:
    f.write("Origin\tHIGH\tMODERATE\tLOW\tMODIFIER\tTotal_Alleles\tIntergenic_Alternate\tIntergenic_Total\n")
    for origin in ['Introduced', 'Native']:
        f.write("\t".join([
            origin,
            str(counters[origin]['HIGH']),
            str(counters[origin]['MODERATE']),
            str(counters[origin]['LOW']),
            str(counters[origin]['MODIFIER']),
            str(counters[origin]['Total_Alleles']),
            str(counters[origin]['Intergenic_Alternate']),
            str(counters[origin]['Intergenic_Total'])
        ]) + "\n")
