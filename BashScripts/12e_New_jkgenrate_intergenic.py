#!/usr/bin/env python3
"""
Generate Fixed Intergenic Sites Script

This script generates a fixed set of 100k intergenic sites from the Italian samples in a VCF file.
The output is saved to a file for use in the jackknife analysis.
"""

import random
import pysam

# Configuration
VCF_PATH = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/All_Origins_annotated.vcf.gz"
ITALIAN_SAMPLES = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Italian"
OUTPUT_FILE = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/3rd/fixed_intergenic_sites.txt"

def load_samples(path):
    """Load sample names from a file into a set."""
    with open(path) as f:
        return set(s.strip() for s in f)

def main():
    # Load Italian samples
    italian_samples = load_samples(ITALIAN_SAMPLES)
    
    # Collect intergenic sites
    intergenic_sites = []
    
    with pysam.VariantFile(VCF_PATH) as vcf:
        # Filter to Italian samples
        vcf.subset_samples(list(italian_samples))
        
        for record in vcf:
            # Parse ANN field to check if variant is intergenic
            ann_entry = record.info.get('ANN', [''])[0]
            if 'intergenic' in ann_entry.lower():
                intergenic_sites.append((record.chrom, record.pos))
    
    # Randomly select 100k intergenic sites (with replacement if needed)
    selected = random.choices(intergenic_sites, k=100000)
    
    # Save to file
    with open(OUTPUT_FILE, 'w') as f:
        for chrom, pos in selected:
            f.write(f"{chrom}\t{pos}\n")
    
    print(f"Saved {len(selected)} intergenic sites to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()