#!/bin/bash
#SBATCH --job-name=rtqpcr_primer_verify
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --output=logs/rtqpcr_primer_verify_%j.out
#SBATCH --error=logs/rtqpcr_primer_verify_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Verify RT-qPCR primers against v1/v2 and the Colombian pseudogenome -
# see codes/analysis/verify_rtqpcr_primers.py for the full checklist
# (BLAST locate -> GFF exon overlap -> EMBOSS water realign -> liftover
# population check). Needs a BLAST-short primer CSV as input.
#
# ── EDIT HERE for a different primer set ──────────────────────────────────
PRIMERS_CSV=${PRIMERS_CSV:-${PROJECT_DIR:-/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data}/analysis/rtqpcr_verification/rtqpcr_primers.csv}
REF_VERSIONS=${REF_VERSIONS:-v1,v2}
POPULATION_CHECK=${POPULATION_CHECK:-1}
# ─────────────────────────────────────────────────────────────────────────

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
mkdir -p logs/ analysis/rtqpcr_verification

module load blast/2.14.1+
module load emboss/6.6.0
module load samtools/1.16.1

POP_FLAG=""
if [ "$POPULATION_CHECK" = "1" ]; then
    POP_FLAG="--population-check"
fi

echo "Start time: $(date)"
echo "Primers CSV: $PRIMERS_CSV  Ref versions: $REF_VERSIONS"

python3 "${PROJECT_DIR}/codes/analysis/verify_rtqpcr_primers.py" \
    --primers-csv "$PRIMERS_CSV" \
    --ref-versions "$REF_VERSIONS" \
    $POP_FLAG

echo "End time: $(date)"
