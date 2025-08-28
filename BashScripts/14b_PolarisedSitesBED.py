#!/usr/bin/env python3
"""
Using the initial VCF and the est-sfs output
Identifies sites where reference allele is ancestral with high confidence threshold min 0.95
Outputs only BED file for bcftools filtering
"""

import pysam # type: ignore
import time
from collections import defaultdict

def main():
    # Configuration
    VCF_PATH = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization/Pol_VCFs/merged_sites_in_ingroup.vcf.gz"
    OUTGROUPS = {"Psi_CAL_N84", "Pva_OU_E88", "SO08"}
    PVAL_FILE = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization/Pmuralis/PmuralisKimura-pval.txt"
    OUTPUT_DIR = "/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarization/Bedtest"
    THRESHOLD = 0.95
    
    # Start timing
    start_time = time.time()
    
    # Read est-sfs probabilities
    print("Reading est-sfs probabilities...")
    est_sfs_probs = {}
    site_counter = 0

    try:
        with open(PVAL_FILE, 'rt') as f:
            for line in f:
                if line.startswith('0'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 3:
                    est_sfs_probs[site_counter] = float(parts[2])
                    site_counter += 1
    except IOError as e:
        print(f"Error reading est-sfs file: {e}")
        return
        
    print(f"Loaded {len(est_sfs_probs)} probabilities from est-sfs")
    
    # Open VCF file
    try:
        vcf = pysam.VariantFile(VCF_PATH)
    except IOError as e:
        print(f"Error opening VCF file: {e}")
        return
    
    # Identify ingroup samples
    all_samples = list(vcf.header.samples)
    ingroup_samples = [s for s in all_samples if s not in OUTGROUPS]
    
    # Precompute sample indices for faster access
    sample_indices = {sample: i for i, sample in enumerate(all_samples)}
    ingroup_indices = [sample_indices[s] for s in ingroup_samples]
    
    print(f"Processing {len(ingroup_samples)} ingroup samples")
    
    # Prepare output files
    high_conf_sites = set()  # Use set to avoid duplicates
    output_file = f"{OUTPUT_DIR}/high_confidence_ancestral_ref.txt"
    
    # Process each variant
    print("Processing variants...")
    processed_count = 0
    skipped_non_snp = 0
    skipped_no_prob = 0
    
    with open(output_file, 'wt') as out_f:
        out_f.write("CHROM\tPOS\tREF\tALT\tREF_FREQ\tALT_FREQ\tP_ANCIENTRAL\tANCESTRAL_ALLELE\n")
        
        for record in vcf:
            processed_count += 1
            if processed_count % 100000 == 0:
                elapsed = time.time() - start_time
                print(f"Processed {processed_count} variants in {elapsed:.1f}s...")
                
            # Skip if not a SNP
            if len(record.ref) != 1 or any(len(alt) != 1 for alt in record.alts):
                skipped_non_snp += 1
                continue
                
            # Get est-sfs probability for this site
            site_idx = processed_count - 1
            if site_idx not in est_sfs_probs:
                skipped_no_prob += 1
                continue
                
            p_ancestral = est_sfs_probs[site_idx]
            
            # Calculate allele frequencies in ingroup only
            ref_count = 0
            alt_count = 0
            
            for sample_idx in ingroup_indices:
                sample_name = all_samples[sample_idx]
                gt = record.samples[sample_name].get('GT', None)
                if gt is None or None in gt:
                    continue
                    
                for allele in gt:
                    if allele == 0:
                        ref_count += 1
                    elif allele == 1:
                        alt_count += 1
            
            total = ref_count + alt_count
            if total == 0:
                continue
                
            ref_freq = ref_count / total
            alt_freq = alt_count / total
            
            # Determine which allele is ancestral based on est-sfs probability
            ancestral_allele = "unknown"
            
            if p_ancestral > THRESHOLD:
                # Major allele is ancestral (could be REF or ALT)
                if ref_freq >= 0.5:
                    ancestral_allele = "REF"
                    high_conf_sites.add((record.chrom, record.pos))
                else:
                    ancestral_allele = "ALT"
            
            elif p_ancestral < (1 - THRESHOLD):
                # Major allele is derived (minor allele is ancestral)
                if ref_freq >= 0.5:
                    ancestral_allele = "ALT"  # REF is major but derived
                else:
                    ancestral_allele = "REF"  # ALT is major but derived, so REF is ancestral
                    high_conf_sites.add((record.chrom, record.pos))
            
            # Write to output for all sites with confident calls
            if ancestral_allele != "unknown":
                out_f.write(f"{record.chrom}\t{record.pos}\t{record.ref}\t{record.alts[0]}\t{ref_freq:.4f}\t{alt_freq:.4f}\t{p_ancestral:.6f}\t{ancestral_allele}\n")
    
    # Create BED file for sites where REF is ancestral with high confidence
    bed_file = f"{OUTPUT_DIR}/high_confidence_ref_ancestral.bed"
    with open(bed_file, 'wt') as bed_f:
        for chrom, pos in sorted(high_conf_sites):
            bed_f.write(f"{chrom}\t{pos-1}\t{pos}\n")
    
    # Report results
    elapsed = time.time() - start_time
    print(f"Analysis complete in {elapsed:.2f} seconds")
    print(f"Processed {processed_count} variants")
    print(f"Skipped {skipped_non_snp} non-SNP variants")
    print(f"Skipped {skipped_no_prob} variants without est-sfs probabilities")
    print(f"Found {len(high_conf_sites)} high-confidence sites where REF is ancestral")
    print(f"BED file created: {bed_file}")
    print(f"Detailed results saved to: {output_file}")

if __name__ == "__main__":
    main()