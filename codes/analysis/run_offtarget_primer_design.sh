#!/bin/bash
#SBATCH --job-name=offtarget_primers
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --output=logs/offtarget_primers_%j.out
#SBATCH --error=logs/offtarget_primers_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Design PCR primer pairs around CRISPR on-target/off-target cut sites -
# see codes/analysis/design_offtarget_primers.py for the full pipeline
# (eprimer3 design + primersearch specificity + pseudogenome population-
# variant check). Requires codes/analysis/setup_primer3.sh to have been
# run once already (installs the primer3_core binary eprimer3 needs).
#
# ── MODIFICAR AQUÍ para otro gen/CSV de sitios ──────────────────────────
GENE=${GENE:-bdnf}
SITES_CSV=${SITES_CSV:-${PROJECT_DIR:-/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data}/crispresso/offtargets/combined/combined_offtargets.csv}
REF_VERSION=${REF_VERSION:-v1}
# ─────────────────────────────────────────────────────────────────────────

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
mkdir -p logs/ analysis/offtarget_primers

module load emboss/6.6.0
module load minimap2
module load samtools/1.16.1

export EMBOSS_PRIMER3_CORE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso/envs/primer3_env/bin/primer3_core

echo "Start time: $(date)"
echo "Gene: $GENE  Sites: $SITES_CSV  REF_VERSION: $REF_VERSION"

python3 "${PROJECT_DIR}/codes/analysis/design_offtarget_primers.py" \
    --gene "$GENE" \
    --sites-csv "$SITES_CSV" \
    --ref-version "$REF_VERSION"

echo "End time: $(date)"
