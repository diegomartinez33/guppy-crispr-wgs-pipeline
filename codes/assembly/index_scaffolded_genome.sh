#!/bin/bash
#SBATCH --job-name=index_scaffolded
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/index_scaffolded_%j.out
#SBATCH --error=logs/index_scaffolded_%j.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

module load samtools/1.16.1
module load bwa/0.7.17
module load blast/2.14.1+

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
GENOME=${PROJECT_DIR}/reference/colombian_scaffolded_genome/colombian_scaffolded.fna
BLAST_DB=${PROJECT_DIR}/reference/colombian_scaffolded_genome/blast_db/colombian_scaffolded

mkdir -p ${PROJECT_DIR}/reference/colombian_scaffolded_genome/blast_db
echo "Start: $(date)"

# ── 1. samtools faidx ─────────────────────────────────────────────────────────
echo ""
echo "[ 1. samtools faidx ]"
samtools faidx "$GENOME"
echo "  ✅ $(basename ${GENOME}).fai — $(wc -l < ${GENOME}.fai) contigs indexed"

# ── 2. BWA index ──────────────────────────────────────────────────────────────
echo ""
echo "[ 2. BWA index ]"
bwa index "$GENOME"
BWA_FILES=$(ls ${GENOME}.{amb,ann,bwt,pac,sa} 2>/dev/null | wc -l)
echo "  ✅ BWA index complete ($BWA_FILES/5 files)"

# ── 3. BLAST database ─────────────────────────────────────────────────────────
echo ""
echo "[ 3. BLAST database ]"
makeblastdb \
    -in "$GENOME" \
    -dbtype nucl \
    -out "$BLAST_DB" \
    -title "Colombian_scaffolded_genome" \
    -parse_seqids
echo "  ✅ BLAST database → $(dirname $BLAST_DB)/"

# ── 4. GATK sequence dictionary ───────────────────────────────────────────────
echo ""
echo "[ 4. GATK sequence dictionary ]"
module load gatk4/4.4.0.0
gatk CreateSequenceDictionary -R "$GENOME"
echo "  ✅ $(basename ${GENOME%.fna}).dict"

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "========================================================"
echo " INDEX SUMMARY"
echo "========================================================"
ls -lh ${PROJECT_DIR}/reference/colombian_scaffolded_genome/
echo "========================================================"
echo "End: $(date)"
