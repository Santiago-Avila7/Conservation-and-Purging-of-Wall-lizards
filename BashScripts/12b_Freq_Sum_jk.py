#!/usr/bin/env python3
"""
Jackknife Analysis for Relative Frequency Sums (fxy/fyx)
Calculates sum of per-site relative frequencies between introduced/native groups
"""

import sys
import pysam # type: ignore
from collections import defaultdict

# Configuration
JACKKNIFE_ID = int(sys.argv[1])  # SLURM_ARRAY_TASK_ID
VCF_PATHS = {
    'HIGH': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts/All_high.vcf.gz",
    'MODERATE': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts/All_moderate.vcf.gz",
    'LOW': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts/All_low.vcf.gz",
    'MODIFIER': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts/All_modifier.vcf.gz",
    'INTERGENIC': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Impacts/All_intergenic.vcf.gz"
}
BLOCK_FILE = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/renamed_blocks.bed"
SAMPLE_LISTS = {
    'introduced': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Italian_Introduced",
    'native': "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Italian_Native"
}
OUTPUT_DIR = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/IntroducedVsNative/Italian"

# ========================
# INITIALIZATION SECTION
# ========================
# Initialize relative frequency counters
sum_relative = {
    'fxy': defaultdict(float),  # Σ(f_introduced * (1 - f_native))
    'fyx': defaultdict(float)   # Σ(f_native * (1 - f_introduced))
}

def load_samples(path):
    """Load sample names while preserving case and order"""
    with open(path) as f:
        return [s.strip() for s in f if s.strip()]

# Load sample lists
sample_sets = {k: load_samples(v) for k, v in SAMPLE_LISTS.items()}

# Debug: Verify sample loading with names
print("\n=== SAMPLE LOADING DEBUG ===")
print(f"Introduced samples ({len(sample_sets['introduced'])}):")
print("---------------------------")
for i, sample in enumerate(sample_sets['introduced'][:5]):  # First 5 samples
    print(f"{i+1}. {sample}")
if len(sample_sets['introduced']) > 5:
    print(f"... plus {len(sample_sets['introduced']) - 5} more")

print(f"\nNative samples ({len(sample_sets['native'])}):")
print("-----------------------")
for i, sample in enumerate(sample_sets['native'][:5]):  # First 5 samples
    print(f"{i+1}. {sample}")
if len(sample_sets['native']) > 5:
    print(f"... plus {len(sample_sets['native']) - 5} more")
    
# Load genomic block to exclude
with open(BLOCK_FILE) as f:
    blocks = [line.strip().split('\t') for line in f]
current_block = blocks[JACKKNIFE_ID - 1]
exclude_chrom, exclude_start = current_block[0], int(current_block[1]) + 1
exclude_end = int(current_block[2])

# ========================
# GENOTYPE PROCESSING FUNCTIONS
# ========================
def process_genotype(gt):
    """Process genotype and return alternate count and validity"""
    if gt in [(None, None), './.']:
        return 0, False
    return sum(1 for a in gt if a > 0), True

# ========================
# MAIN PROCESSING LOOP
# ========================
print("\n=== BEGIN PROCESSING VCFs ===")

for impact_category, vcf_path in VCF_PATHS.items():
    print(f"\nProcessing {impact_category} category...")
    
    try:
        with pysam.VariantFile(vcf_path) as vcf:
            # Get valid samples for both groups
            vcf_samples = set(vcf.header.samples)
            intro_samples = [s for s in sample_sets['introduced'] if s in vcf_samples]
            native_samples = [s for s in sample_sets['native'] if s in vcf_samples]
            
            # Process each variant record
            for record in vcf:
                # Skip variants in excluded block
                chrom, pos = record.chrom, record.start + 1
                if chrom == exclude_chrom and exclude_start <= pos <= exclude_end:
                    continue
                
                # Initialize per-site counters
                intro_alt, intro_total = 0, 0
                native_alt, native_total = 0, 0
                
                # Process introduced samples
                for sample in intro_samples:
                    gt = record.samples[sample].get('GT', (None, None))
                    alts, is_called = process_genotype(gt)
                    if is_called:
                        intro_alt += alts
                        intro_total += 2  # Diploid count
                
                # Process native samples
                for sample in native_samples:
                    gt = record.samples[sample].get('GT', (None, None))
                    alts, is_called = process_genotype(gt)
                    if is_called:
                        native_alt += alts
                        native_total += 2
                
                # Calculate relative frequencies if both groups have data
                if intro_total > 0 and native_total > 0:
                    f_intro = intro_alt / intro_total
                    f_native = native_alt / native_total
                    
                    # Calculate and sum relative frequencies
                    sum_relative['fxy'][impact_category] += f_intro * (1 - f_native)
                    sum_relative['fyx'][impact_category] += f_native * (1 - f_intro)
                    
    except Exception as e:
        print(f"ERROR processing {vcf_path}: {str(e)}", file=sys.stderr)
        continue

# ========================
# OUTPUT GENERATION
# ========================
output_path = f"{OUTPUT_DIR}/jackknife_{JACKKNIFE_ID}_relative.tsv"
print(f"\n=== WRITING RESULTS TO {output_path} ===")

# Create header
categories = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER', 'INTERGENIC']
header = "\t".join([f"Sum_fxy_{cat}\tSum_fyx_{cat}" for cat in categories]) + "\n"

with open(output_path, 'w') as f:
    f.write(header)
    
    # Build data line
    line_parts = []
    for cat in categories:
        line_parts.append(f"{sum_relative['fxy'][cat]:.6f}")
        line_parts.append(f"{sum_relative['fyx'][cat]:.6f}")
    
    f.write("\t".join(line_parts) + "\n")

print("\n=== ANALYSIS COMPLETE ===")