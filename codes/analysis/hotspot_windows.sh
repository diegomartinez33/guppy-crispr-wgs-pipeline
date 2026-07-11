#!/bin/bash
#SBATCH --job-name=hotspot_windows
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=logs/hotspot_windows_%j.out
#SBATCH --error=logs/hotspot_windows_%j.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

module load bcftools/1.15.1
module load bedtools/2.30.0

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
VCF_DIR=${PROJECT_DIR}/gatk/trimmomatic/vcf_filtered
OUT_DIR=${PROJECT_DIR}/gatk/trimmomatic/hotspots

mkdir -p "$OUT_DIR"

TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR
trap "rm -rf $TMPDIR" EXIT

echo "=== Step 1: Build genome.txt from VCF header ==="
bcftools view -h "${VCF_DIR}/snps_filtered.vcf.gz" | \
  grep "^##contig" | \
  sed 's/.*ID=\([^,]*\),length=\([0-9]*\).*/\1\t\2/' \
  > "${TMPDIR}/genome.txt"

echo "Chromosomes found: $(wc -l < ${TMPDIR}/genome.txt)"

echo "=== Step 2: Generate 10 kb sliding windows (2 kb step) ==="
bedtools makewindows \
  -g "${TMPDIR}/genome.txt" \
  -w 10000 \
  -s 2000 \
  > "${TMPDIR}/windows_10kb_2kb.bed"

echo "Windows generated: $(wc -l < ${TMPDIR}/windows_10kb_2kb.bed)"

echo "=== Step 3: Extract PASS SNP positions as 0-based BED ==="
bcftools view -f PASS "${VCF_DIR}/snps_filtered.vcf.gz" | \
  bcftools query -f '%CHROM\t%POS\t%END\n' | \
  awk '{print $1"\t"$2-1"\t"$3}' \
  > "${TMPDIR}/snps.bed"

echo "PASS SNP records: $(wc -l < ${TMPDIR}/snps.bed)"

echo "=== Step 4: Extract PASS INDEL positions as 0-based BED ==="
bcftools view -f PASS "${VCF_DIR}/indels_filtered.vcf.gz" | \
  bcftools query -f '%CHROM\t%POS\t%END\n' | \
  awk '{print $1"\t"$2-1"\t"$3}' \
  > "${TMPDIR}/indels.bed"

echo "PASS INDEL records: $(wc -l < ${TMPDIR}/indels.bed)"

echo "=== Step 5: Merge SNPs + INDELs into combined BED ==="
cat "${TMPDIR}/snps.bed" "${TMPDIR}/indels.bed" | \
  sort -k1,1 -k2,2n \
  > "${TMPDIR}/combined_variants.bed"

echo "=== Step 6: Count variants per window ==="
bedtools coverage -counts \
  -a "${TMPDIR}/windows_10kb_2kb.bed" \
  -b "${TMPDIR}/snps.bed" \
  > "${TMPDIR}/window_counts_snp.bed"

bedtools coverage -counts \
  -a "${TMPDIR}/windows_10kb_2kb.bed" \
  -b "${TMPDIR}/indels.bed" \
  > "${TMPDIR}/window_counts_indel.bed"

bedtools coverage -counts \
  -a "${TMPDIR}/windows_10kb_2kb.bed" \
  -b "${TMPDIR}/combined_variants.bed" \
  > "${TMPDIR}/window_counts_combined.bed"

echo "=== Step 7: Combine into 6-column TSV ==="
paste \
  <(awk '{print $1"\t"$2"\t"$3"\t"$4}' "${TMPDIR}/window_counts_snp.bed") \
  <(awk '{print $4}' "${TMPDIR}/window_counts_indel.bed") \
  <(awk '{print $4}' "${TMPDIR}/window_counts_combined.bed") \
  > "${OUT_DIR}/window_counts.tsv"

echo "Windows written: $(wc -l < ${OUT_DIR}/window_counts.tsv)"
echo "Total SNP counts in windows:   $(awk '{s+=$4} END{print s}' ${OUT_DIR}/window_counts.tsv)"
echo "Total INDEL counts in windows: $(awk '{s+=$5} END{print s}' ${OUT_DIR}/window_counts.tsv)"
echo "Done: ${OUT_DIR}/window_counts.tsv"
