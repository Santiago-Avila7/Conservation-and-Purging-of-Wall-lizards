#!/usr/bin/env python3
"""
Jackknife Analysis Script for Categorized VCFs (HIGH, MODERATE, LOW, MODIFIER, INTERGENIC)
Counts alternate alleles and total called alleles per category, excluding missing genotypes.
"""

import sys
import pysam
from collections import defaultdict

# Configuration
JACKKNIFE_ID = int(sys.argv[1])  # SLURM_ARRAY_TASK_ID
VCF_PATHS = {
    'HIGH': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts/All_Italian_high.vcf.gz",
    'MODERATE': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts/All_Italian_moderate.vcf.gz",
    'LOW': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts/All_Italian_low.vcf.gz",
    'MODIFIER': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts/All_Italian_modifier.vcf.gz",
    'INTERGENIC': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts/All_Italian_intergenic.vcf.gz"
}
BLOCK_FILE = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/renamed_blocks.bed"
SAMPLE_LISTS = {
    'introduced': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Introduced_Purging_Test",
    'native': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Italian_Native"
}
OUTPUT_DIR = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/Italian_SizeControl_ItalianVCF"

# ========================
# INITIALIZATION SECTION
# ========================
# Initialize counters with nested defaultdicts
# Structure: counters[population][impact_category]['Alternate'|'Total']
counters = {
    'introduced': defaultdict(lambda: {'Alternate': 0, 'Total': 0}),
    'native': defaultdict(lambda: {'Alternate': 0, 'Total': 0})
}

def load_samples(path):
    """Load sample names from file, preserving original case and whitespace"""
    with open(path) as f:
        return [s.strip() for s in f if s.strip()]  # Keep empty line checking

# Load sample lists exactly as they appear in files
sample_sets = {k: load_samples(v) for k, v in SAMPLE_LISTS.items()}

# Debug: Verify sample loading
print(f"\n=== SAMPLE LOADING DEBUG ===", file=sys.stdout)
print(f"Loaded {len(sample_sets['introduced'])} introduced samples", file=sys.stdout)
print(f"Loaded {len(sample_sets['native'])} native samples", file=sys.stdout)
print(f"First 3 introduced: {sample_sets['introduced'][:3]}", file=sys.stdout)
print(f"First 3 native: {sample_sets['native'][:3]}", file=sys.stdout)

# Load genomic block to exclude for this jackknife iteration
with open(BLOCK_FILE) as f:
    blocks = [line.strip().split('\t') for line in f]
current_block = blocks[JACKKNIFE_ID - 1]
exclude_chrom, exclude_start = current_block[0], int(current_block[1]) + 1  # 1-based
exclude_end = int(current_block[2])

# ========================
# GENOTYPE PROCESSING FUNCTIONS
# ========================
def process_genotype(gt):
    """
    Robust genotype processor that handles:
    - Phased (0|1) and unphased (0/1) genotypes
    - Missing data (./., None)
    - Multiple alternate alleles
    
    Returns:
    - alternate_count: Number of non-reference alleles (1, 2, etc.)
    - is_called: Boolean indicating if genotype was valid
    """
    if gt is None or gt == (None, None) or gt == "./.":
        return 0, False
    
    # Convert all genotypes to uniform format
    gt_str = "|".join(str(a) for a in gt if a is not None)
    alleles = []
    for part in gt_str.replace("/", "|").split("|"):
        if part.isdigit():
            alleles.append(int(part))
    
    # Count alternate alleles (anything > 0)
    alternate_count = sum(1 for a in alleles if a > 0)
    is_called = len(alleles) > 0
    
    return alternate_count, is_called

# ========================
# MAIN PROCESSING LOOP
# ========================
print("\n=== BEGIN PROCESSING VCFs ===", file=sys.stdout)

for impact_category, vcf_path in VCF_PATHS.items():
    print(f"\nProcessing {impact_category} VCF...", file=sys.stdout)
    
    try:
        with pysam.VariantFile(vcf_path) as vcf:
            # Get actual sample names in VCF (case-sensitive)
            vcf_samples = set(vcf.header.samples)
            print(f"VCF contains {len(vcf_samples)} samples", file=sys.stdout)
            
            # Process each population group
            for group in ['introduced', 'native']:
                print(f"\nProcessing {group} samples...", file=sys.stdout)
                
                # Find intersection between our samples and VCF samples
                valid_samples = [s for s in sample_sets[group] if s in vcf_samples]
                print(f"Found {len(valid_samples)} valid samples", file=sys.stdout)
                
                if not valid_samples:
                    print(f"WARNING: No valid samples for {group} in {impact_category} VCF!", file=sys.stdout)
                    continue
                
                # Reset VCF reader for each group to ensure complete processing
                vcf.reset()
                records_processed = 0
                
                for record in vcf:
                    # Skip variants in excluded block
                    chrom, pos = record.chrom, record.start + 1
                    if chrom == exclude_chrom and exclude_start <= pos <= exclude_end:
                        continue
                    
                    # Initialize counters for this variant
                    alt_count = total_alleles = 0
                    
                    # Process all samples in this group
                    for sample in valid_samples:
                        gt = record.samples[sample].get('GT', None)
                        alts, is_called = process_genotype(gt)
                        if is_called:
                            alt_count += alts
                            total_alleles += 2  # Diploid count
                    
                    # Update global counters
                    counters[group][impact_category]['Alternate'] += alt_count
                    counters[group][impact_category]['Total'] += total_alleles
                    records_processed += 1
                
                print(f"Processed {records_processed} variants for {group}", file=sys.stdout)
    
    except Exception as e:
        print(f"ERROR processing {vcf_path}: {str(e)}", file=sys.stderr)
        continue

# ========================
# OUTPUT GENERATION
# ========================
output_path = f"{OUTPUT_DIR}/jackknife_{JACKKNIFE_ID}_results.tsv"
print(f"\n=== WRITING RESULTS TO {output_path} ===", file=sys.stdout)

with open(output_path, 'w') as f:
    # Write header
    header_columns = []
    for cat in ['HIGH', 'MODERATE', 'LOW', 'MODIFIER', 'INTERGENIC']:
        header_columns.extend([f"Alternate_{cat}", f"Total_{cat}"])
    header = "Origin\t" + "\t".join(header_columns) + "\n"
    f.write(header)
    
    # Write data for each group
    for group in ['introduced', 'native']:
        # Collect all counts in order
        count_values = []
        for cat in ['HIGH', 'MODERATE', 'LOW', 'MODIFIER', 'INTERGENIC']:
            count_values.append(str(counters[group][cat]['Alternate']))
            count_values.append(str(counters[group][cat]['Total']))
        
        line = f"{group}\t" + "\t".join(count_values) + "\n"
        f.write(line)
        
        # Debug output
        print(f"\nFinal counts for {group}:", file=sys.stdout)
        for cat in ['HIGH', 'MODERATE', 'LOW', 'MODIFIER', 'INTERGENIC']:
            print(f"{cat}: {counters[group][cat]['Alternate']}/{counters[group][cat]['Total']}", file=sys.stdout)

print("\n=== ANALYSIS COMPLETE ===", file=sys.stdout)