#!/bin/bash
#SBATCH --job-name=make_pseudogenome
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/make_pseudogenome_%j.out
#SBATCH --error=logs/make_pseudogenome_%j.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

module load bcftools/1.15.1
module load samtools/1.16.1
module load bwa/0.7.17

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
VCF_DIR=${PROJECT_DIR}/gatk/trimmomatic/vcf_filtered
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
OUT_DIR=${PROJECT_DIR}/reference/pseudogenome

mkdir -p "$OUT_DIR"

TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR
trap "rm -rf $TMPDIR" EXIT

CONTROLS="Control_MNP_I_S54_L002,Control_MNP_II_S55_L002,Control_MNP_III_S56_L002"
AF_THRESH=0.667

echo "=== Step 1: Filter SNPs to Control samples at majority frequency ==="
bcftools view -f PASS -s "$CONTROLS" "${VCF_DIR}/snps_filtered.vcf.gz" | \
  bcftools view -g ^miss | \
  bcftools +fill-tags -Oz -o "${TMPDIR}/ctrl_snps_filled.vcf.gz" -- -t AF,AC,AN
bcftools index -t "${TMPDIR}/ctrl_snps_filled.vcf.gz"

bcftools view -e "AF<${AF_THRESH}" -Oz \
  -o "${TMPDIR}/ctrl_snps_hq.vcf.gz" \
  "${TMPDIR}/ctrl_snps_filled.vcf.gz"
bcftools index -t "${TMPDIR}/ctrl_snps_hq.vcf.gz"
echo "Control SNPs (AF >= ${AF_THRESH}): $(bcftools view -H ${TMPDIR}/ctrl_snps_hq.vcf.gz | wc -l)"

echo "=== Step 2: Filter INDELs to Control samples at majority frequency ==="
bcftools view -f PASS -s "$CONTROLS" "${VCF_DIR}/indels_filtered.vcf.gz" | \
  bcftools view -g ^miss | \
  bcftools +fill-tags -Oz -o "${TMPDIR}/ctrl_indels_filled.vcf.gz" -- -t AF,AC,AN
bcftools index -t "${TMPDIR}/ctrl_indels_filled.vcf.gz"

bcftools view -e "AF<${AF_THRESH}" -Oz \
  -o "${TMPDIR}/ctrl_indels_hq.vcf.gz" \
  "${TMPDIR}/ctrl_indels_filled.vcf.gz"
bcftools index -t "${TMPDIR}/ctrl_indels_hq.vcf.gz"
echo "Control INDELs (AF >= ${AF_THRESH}): $(bcftools view -H ${TMPDIR}/ctrl_indels_hq.vcf.gz | wc -l)"

echo "=== Step 3: Merge SNPs + INDELs ==="
bcftools concat -a -D \
  "${TMPDIR}/ctrl_snps_hq.vcf.gz" \
  "${TMPDIR}/ctrl_indels_hq.vcf.gz" | \
  bcftools sort -Oz -o "${TMPDIR}/ctrl_variants_hq.vcf.gz"
bcftools index -t "${TMPDIR}/ctrl_variants_hq.vcf.gz"

echo "=== Variant count audit ==="
bcftools stats "${TMPDIR}/ctrl_variants_hq.vcf.gz" | grep "^SN"

echo "=== Step 4: Apply variants to reference (bcftools consensus) ==="
bcftools consensus \
  -f "$REF" \
  -s Control_MNP_I_S54_L002 \
  -H A \
  -c "${OUT_DIR}/colombian_pseudogenome.chain" \
  "${TMPDIR}/ctrl_variants_hq.vcf.gz" \
  -o "${OUT_DIR}/colombian_pseudogenome.fna"

echo "=== Step 5: Index pseudogenome ==="
samtools faidx "${OUT_DIR}/colombian_pseudogenome.fna"
bwa index -a bwtsw "${OUT_DIR}/colombian_pseudogenome.fna"

echo "=== Verification ==="
echo "Contig count: $(grep -c '^>' ${OUT_DIR}/colombian_pseudogenome.fna)"
echo "BWA index files:"
ls -lh "${OUT_DIR}/colombian_pseudogenome.fna."{amb,ann,bwt,pac,sa}
echo "Done. Pseudogenome → ${OUT_DIR}/colombian_pseudogenome.fna"
