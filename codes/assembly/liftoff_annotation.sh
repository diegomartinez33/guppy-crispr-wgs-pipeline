#!/bin/bash
#SBATCH --job-name=liftoff_annotation
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --output=logs/liftoff_annotation_%j.out
#SBATCH --error=logs/liftoff_annotation_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Phase 4: Annotation transfer (Liftoff) - maps Trinidad reference gene
# models, including bdnf, onto the new Colombian scaffolded genome.
# Requires 00_download_gff.sh and 00_setup_liftoff_env.sh to have been run.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
TRINIDAD_GFF=${PROJECT_DIR}/reference/GCF_000633615.1_annotation.gff
NEW_FASTA=${PROJECT_DIR}/assembly/ragtag_output/ragtag.scaffold.fasta
OUT_DIR=${PROJECT_DIR}/reference/colombian_scaffolded_genome
NEW_GFF=${OUT_DIR}/colombian_scaffolded.liftoff.gff3

mkdir -p "$OUT_DIR" logs/

if [ ! -f "$TRINIDAD_GFF" ]; then
    echo "ERROR: Trinidad GFF not found at $TRINIDAD_GFF - run 00_download_gff.sh first"
    exit 1
fi
if [ ! -f "$NEW_FASTA" ]; then
    echo "ERROR: ragtag.scaffold.fasta not found at $NEW_FASTA - did Phase 2 (ragtag_scaffold.sh) finish?"
    exit 1
fi

CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate liftoff_env

cp "$NEW_FASTA" "${OUT_DIR}/colombian_scaffolded.fna"

echo "Start time: $(date)"
liftoff -g "$TRINIDAD_GFF" -o "$NEW_GFF" -p "${SLURM_CPUS_PER_TASK}" "$NEW_FASTA" "$REF"
echo "End time: $(date)"

echo "=== bdnf locus check ==="
if grep -q "ID=gene-bdnf" "$NEW_GFF"; then
    echo "bdnf gene successfully lifted over:"
    grep "ID=gene-bdnf" "$NEW_GFF"
else
    echo "WARNING: bdnf gene not found in lifted annotation"
fi
