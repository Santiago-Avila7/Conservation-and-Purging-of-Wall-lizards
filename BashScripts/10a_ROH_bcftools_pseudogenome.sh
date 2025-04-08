#!/bin/bash -l
#SBATCH -A naiss2025-22-189
#SBATCH -J ROH_Pseudogenome
#SBATCH --output=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/ROHs/%x_%j.out
#SBATCH --error=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/logs/ROHs/%x_%j.err
#SBATCH -p shared
#SBATCH -t 4:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sa5674av-s@student.lu.se

# =============================================
# 1. Load required modules
# =============================================
module load bcftools/1.20

# =============================================
# 2. Define input VCFs and directories
# =============================================
datadir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/Diversity_Data"
savedir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Gene_Diversity/ROHs"

FRENCH_VCF="$datadir/All_French_final.vcf.gz"
ITALIAN_VCF="$datadir/All_Italian_final.vcf.gz"

# =============================================
# 3. Function to create pseudogenome and compute total ROH length
# =============================================
run_pseudogenome_roh() {
  local INPUT_VCF="$1"
  local LABEL="$2"

  local PSEUDO_VCF="${savedir}/${LABEL}_pseudohom.vcf.gz"
  local ROH_OUTPUT="${savedir}/${LABEL}_pseudohom_ROH.txt.gz"
  local LENGTH_FILE="${savedir}/${LABEL}_pseudohom_ROH_length.tsv"
  local HEADER_TMP="${savedir}/${LABEL}_header.tmp"

  echo "[${LABEL}] Creating pseudogenome VCF with all 0/0 genotypes..."

  mkdir -p "$savedir"

  # Build header
  zcat "$INPUT_VCF" | grep "^##" > "$HEADER_TMP"
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpseudohom" >> "$HEADER_TMP"

  # Create records with 0/0 genotypes
  zcat "$INPUT_VCF" | grep -v "^#" | \
    sed -e 's/\\$/\t0\/0/' | \
    cat "$HEADER_TMP" - | bgzip > "$PSEUDO_VCF"

  bcftools index "$PSEUDO_VCF"

  echo "[${LABEL}] Running bcftools roh..."
  bcftools roh -G 30 -Orz -o "$ROH_OUTPUT" "$PSEUDO_VCF"

  echo "[${LABEL}] Summarizing ROH length..."
  zcat "$ROH_OUTPUT" | \
    awk -v s="pseudohom" 'BEGIN{sum=0} $2==s {sum+=$6} END{print s"\t"sum}' > "$LENGTH_FILE"

  echo "[${LABEL}] Done. Length saved to $LENGTH_FILE"

  # Cleanup temp files
  rm -f "$HEADER_TMP" "$ROH_OUTPUT"
}

# =============================================
# 4. Run for French and Italian VCFs
# =============================================
run_pseudogenome_roh "$FRENCH_VCF" "French"
run_pseudogenome_roh "$ITALIAN_VCF" "Italian"

# =============================================
# 5. Merge summary outputs
# =============================================
echo -e "Label\tROH_Total_Length" > "$savedir/pseudogenome_roh_summary.tsv"
cat "$savedir/French_pseudohom_ROH_length.tsv" >> "$savedir/pseudogenome_roh_summary.tsv"
cat "$savedir/Italian_pseudohom_ROH_length.tsv" >> "$savedir/pseudogenome_roh_summary.tsv"

echo "All pseudogenome ROH analyses completed. Summary saved to pseudogenome_roh_summary.tsv"
