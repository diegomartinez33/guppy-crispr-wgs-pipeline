#!/bin/bash
#SBATCH --job-name=crossmap_pseudo
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/crossmap_pseudogenome_%j.out
#SBATCH --error=logs/crossmap_pseudogenome_%j.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate crossmap_env
export PATH="${CONDA_BASE}/envs/crossmap_env/bin:$PATH"

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
CHAIN=${PROJECT_DIR}/reference/pseudogenome/colombian_pseudogenome.chain
GFF_IN=${PROJECT_DIR}/reference/GCF_000633615.1_annotation.gff
OUT_DIR=${PROJECT_DIR}/reference/pseudogenome
GFF_OUT=${OUT_DIR}/colombian_pseudogenome.gff3

echo "Start time: $(date)"
echo "CrossMap version: $(CrossMap --version 2>&1)"
echo "Chain file: $CHAIN"
echo "Input GFF:  $GFF_IN"
echo "Output GFF: $GFF_OUT"

echo ""
echo "=== Running CrossMap GFF liftover ==="
CrossMap gff \
    "$CHAIN" \
    "$GFF_IN" \
    "$GFF_OUT"

echo ""
echo "=== Verification ==="
TOTAL_IN=$(grep -v "^#" "$GFF_IN" | wc -l)
TOTAL_OUT=$(grep -v "^#" "$GFF_OUT" | wc -l)
echo "Input features:  $TOTAL_IN"
echo "Output features: $TOTAL_OUT"
echo "Transfer rate:   $(awk "BEGIN {printf \"%.1f%%\", ${TOTAL_OUT}/${TOTAL_IN}*100}")"

echo ""
echo "=== bdnf gene check ==="
grep -i "bdnf\|gene-bdnf" "$GFF_OUT" | head -3 || echo "bdnf not found by name — checking by region"
grep "NC_024333.1" "$GFF_OUT" | awk '$3=="gene"' | \
    awk '$4 > 15800000 && $5 < 16100000' | head -3

echo ""
echo "=== Feature type counts ==="
grep -v "^#" "$GFF_OUT" | awk '{print $3}' | sort | uniq -c | sort -rn | head -15

echo ""
echo "End time: $(date)"
echo "Done. Annotation → $GFF_OUT"
