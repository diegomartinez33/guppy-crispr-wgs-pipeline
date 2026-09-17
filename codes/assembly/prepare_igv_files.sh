#!/bin/bash
#SBATCH --job-name=prepare_igv
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/prepare_igv_%j.out
#SBATCH --error=logs/prepare_igv_%j.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

module load samtools/1.16.1

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
source "${PROJECT_DIR}/codes/genome_versions.sh"
PSEUDO_DIR=${PROJECT_DIR}/reference/pseudogenome${OUT_SUFFIX}
MERGED_DIR=${PROJECT_DIR}/mapping/trimmomatic${OUT_SUFFIX}/merged
IGV_DIR=${PROJECT_DIR}/igv_files${OUT_SUFFIX}
# features_of_interest.bed is hand-curated per version (bdnf/sgRNA/off-target
# coordinates are genome-specific) and lives directly in IGV_DIR, same as
# v1 - this script never generated it and doesn't touch it.

mkdir -p "$IGV_DIR"
echo "Start: $(date)"

# ── 1. GFF3 — sort, bgzip, tabix ─────────────────────────────────────────────
echo ""
echo "[ 1. Preparing GFF3 annotation ]"
GFF_IN=${PSEUDO_DIR}/colombian_pseudogenome.liftoff.gff3
GFF_SORTED=${IGV_DIR}/colombian_pseudogenome.gff3.gz

echo "  Sorting..."
(grep "^#" "$GFF_IN"; grep -v "^#" "$GFF_IN" | sort -k1,1 -k4,4n) \
    | bgzip -c > "$GFF_SORTED"

echo "  Indexing with tabix..."
tabix -p gff "$GFF_SORTED"

echo "  ✅ GFF3 ready: $(basename $GFF_SORTED)"

# ── 2. Genome FASTA — symlink (already has .fai) ──────────────────────────────
echo ""
echo "[ 2. Genome FASTA ]"
ln -sf "${PSEUDO_DIR}/colombian_pseudogenome.fna"     "${IGV_DIR}/colombian_pseudogenome.fna"
ln -sf "${PSEUDO_DIR}/colombian_pseudogenome.fna.fai" "${IGV_DIR}/colombian_pseudogenome.fna.fai"
echo "  ✅ FASTA + .fai linked"

# ── 3. BAMs — symlink merged group BAMs (already sorted and indexed) ──────────
echo ""
echo "[ 3. Merged group BAMs ]"
for GROUP in Control Only_MNP Plasmid_Ko RNP_Cas; do
    BAM="${MERGED_DIR}/${GROUP}_merged.sorted.bam"
    BAI="${BAM}.bai"
    if [ -f "$BAM" ] && [ -f "$BAI" ]; then
        ln -sf "$BAM" "${IGV_DIR}/${GROUP}_merged.sorted.bam"
        ln -sf "$BAI" "${IGV_DIR}/${GROUP}_merged.sorted.bam.bai"
        echo "  ✅ ${GROUP}"
    else
        echo "  ⚠️  ${GROUP} BAM or index not found — skipping"
    fi
done

# ── 4. Print file manifest and download instructions ──────────────────────────
echo ""
echo "========================================================"
echo " FILES TO DOWNLOAD TO YOUR LOCAL MACHINE"
echo "========================================================"
echo ""
echo "  scp -r da.martinez33@hypatia:${IGV_DIR} ~/igv_guppy/"
echo ""
echo " Or download individual files:"
for f in "$IGV_DIR"/*; do
    echo "    $(du -sh "$f" 2>/dev/null | awk '{print $1}')  $(basename "$f")"
done

echo ""
echo "========================================================"
echo " HOW TO OPEN IN IGV"
echo "========================================================"
echo "  1. Genomes → Load Genome from File"
echo "     → colombian_pseudogenome.fna"
echo ""
echo "  2. File → Load from File"
echo "     → colombian_pseudogenome.gff3.gz  (annotation)"
echo "     → Control_merged.sorted.bam"
echo "     → Only_MNP_merged.sorted.bam"
echo "     → Plasmid_Ko_merged.sorted.bam"
echo "     → RNP_Cas_merged.sorted.bam"
echo ""
echo "  3. Search bar:"
echo "     → bdnf                           (by gene name)"
if [ "$REF_VERSION" = "v1" ]; then
    echo "     → NC_024333.1:15923726-15938393  (bdnf, by coordinates)"
    echo "     → NC_024334.1:13468145-13847119  (nlgn1)"
else
    echo "     → NC_088832.1:15851100-15865643  (bdnf, by coordinates)"
    echo "     → NC_088833.1:31182502-31568166  (nlgn1)"
fi
echo "========================================================"

echo ""
echo "End: $(date)"
echo "Done. Files in: $IGV_DIR"
