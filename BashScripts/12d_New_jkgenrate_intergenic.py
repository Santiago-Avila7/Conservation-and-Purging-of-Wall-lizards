#!/usr/bin/env python3
"""
Generate Fixed Intergenic Sites Script

This script generates a fixed set of 100k intergenic sites from an Italian-only VCF file.
The output is saved to a file for use in jackknife analysis.
"""

import random
import pysam

# Configuration
VCF_PATH = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Origin/All_Italian_annotated.vcf.gz"  # Updated to reflect Italian-only VCF
OUTPUT_FILE = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging/Jackknife/4th/italian_intergenic_sites.txt"

def main():
    # Collect intergenic sites
    intergenic_sites = []
    
    with pysam.VariantFile(VCF_PATH) as vcf:
        for record in vcf:
            # Extract ANN field and check if it's intergenic
            ann_entries = record.info.get('ANN', []) 
            if any('intergenic' in entry.lower() for entry in ann_entries):
                intergenic_sites.append((record.chrom, record.pos))
    
    # Randomly select 100k intergenic sites (with replacement if needed)
    selected = random.choices(intergenic_sites, k=100000) if len(intergenic_sites) >= 100000 else intergenic_sites
    
    # Save to file
    with open(OUTPUT_FILE, 'w') as f:
        for chrom, pos in selected:
            f.write(f"{chrom}\t{pos}\n")
    
    print(f"Saved {len(selected)} intergenic sites to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
