#!/bin/bash
#SBATCH --job-name=crispor_offtarget_scan
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --output=logs/crispor_offtarget_scan_%j.out
#SBATCH --error=logs/crispor_offtarget_scan_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL
set -euo pipefail

# Regenerate the "CRISPOR side" of off-target discovery for the bdnf sgRNA
# (TGAGAGACGCCCCGGGCATG + NGG) by running crispor.py directly against a
# registered reference genome, instead of the old workflow of manually
# downloading a .xls export from crispor.org (data/crispor_offtargets.xls -
# a one-time artifact for the OLD genome only, with no way to regenerate it
# for a new one without repeating that manual step). User-approved
# modernization, see CLAUDE.md, "Migration to GCF_904066995.2 (v2)".
#
# Output columns match crisporWebsite's documented offtarget TSV format:
# seqId, guideId, guideSeq, offtargetSeq, mismatchPos, mismatchCount,
# mitOfftargetScore, cfdOfftargetScore, chrom, start, end, strand, locusDesc
# - convert_crispor_offtargets.py reads this directly (adapted from reading
# the .xls), matching by guideSeq instead of the old website's guideId
# (which won't exist/match in a fresh container run).

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
source "${PROJECT_DIR}/codes/genome_versions.sh"

SIF=${PROJECT_DIR}/codes/analysis/crispor_singularity/crispor_v5.2c_amd64.sif
GENOMES_DIR=${PROJECT_DIR}/codes/analysis/crispor_singularity/genomes

if [ "$REF_VERSION" = "v1" ]; then
    GENOME_ID="guppyRefTrinidad"
    AMPLICON=${PROJECT_DIR}/reference/amplicon_bdnf_100bp.fa
elif [ "$REF_VERSION" = "v2" ]; then
    GENOME_ID="guppyRefMaleV2"
    AMPLICON=${PROJECT_DIR}/reference/amplicon_bdnf_100bp_v2.fa
fi

OUT_DIR=${PROJECT_DIR}/crispresso${OUT_SUFFIX}/offtargets/crispor_container
mkdir -p "$OUT_DIR" logs/

module load singularity/3.7.1

GUIDE_OUT=${OUT_DIR}/bdnf_guides.tsv
OFFT_OUT=${OUT_DIR}/bdnf_offtargets.tsv

echo "Start time: $(date)"
echo "Genome: ${GENOME_ID}  Amplicon: ${AMPLICON}"

singularity exec -B "${GENOMES_DIR}:/data/genomes" "$SIF" \
    /data/www/crispor/crispor.py \
    -g /data/genomes \
    "$GENOME_ID" "$AMPLICON" "$GUIDE_OUT" \
    -o "$OFFT_OUT"

echo ""
echo "=== Guides found ==="
wc -l "$GUIDE_OUT"
echo "=== Off-targets found ==="
wc -l "$OFFT_OUT"

echo ""
echo "=== bdnf guide row (spacer TGAGAGACGCCCCGGGCATG + NGG) ==="
awk -F'\t' 'NR==1 || $3 ~ /^TGAGAGACGCCCCGGGCATG/' "$GUIDE_OUT"

echo "End time: $(date)"
