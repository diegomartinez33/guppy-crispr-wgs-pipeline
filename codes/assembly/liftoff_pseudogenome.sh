#!/bin/bash
#SBATCH --job-name=liftoff_pseudo
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/liftoff_pseudogenome_%j.out
#SBATCH --error=logs/liftoff_pseudogenome_%j.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate liftoff_env

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
TRINIDAD_GFF=${PROJECT_DIR}/reference/GCF_000633615.1_annotation.gff
PSEUDO=${PROJECT_DIR}/reference/pseudogenome/colombian_pseudogenome.fna
OUT_GFF=${PROJECT_DIR}/reference/pseudogenome/colombian_pseudogenome.liftoff.gff3

echo "Start: $(date)"
echo "Target:     $PSEUDO"
echo "Reference:  $REF"
echo "Input GFF:  $TRINIDAD_GFF"
echo "Output GFF: $OUT_GFF"

liftoff \
    -g "$TRINIDAD_GFF" \
    -o "$OUT_GFF" \
    -p "${SLURM_CPUS_PER_TASK}" \
    "$PSEUDO" \
    "$REF"

echo ""
echo "=== Transfer statistics ==="
TOTAL_IN=$(grep -v "^#" "$TRINIDAD_GFF" | wc -l)
TOTAL_OUT=$(grep -v "^#" "$OUT_GFF" | wc -l)
GENES_OUT=$(grep -v "^#" "$OUT_GFF" | awk '$3=="gene"' | wc -l)
echo "Input features:  $TOTAL_IN"
echo "Output features: $TOTAL_OUT"
echo "Genes transferred: $GENES_OUT"
echo "Transfer rate: $(awk "BEGIN {printf \"%.1f%%\", ${TOTAL_OUT}/${TOTAL_IN}*100}")"

echo ""
echo "=== Feature type counts ==="
grep -v "^#" "$OUT_GFF" | awk '{print $3}' | sort | uniq -c | sort -rn | head -10

echo ""
echo "=== bdnf check ==="
if grep -q "ID=gene-bdnf" "$OUT_GFF"; then
    echo "bdnf gene transferred successfully:"
    grep "ID=gene-bdnf" "$OUT_GFF" | awk '{print $1, $4, $5, $7}'
else
    echo "WARNING: bdnf gene not found"
fi

echo ""
echo "=== Hierarchy integrity check ==="
echo "  Genes without children (orphan gene records):"
GENE_IDS=$(grep -v "^#" "$OUT_GFF" | awk '$3=="gene"' | grep -oP 'ID=[^;]+' | sort)
PARENT_IDS=$(grep -v "^#" "$OUT_GFF" | grep -oP 'Parent=[^;]+' | sed 's/Parent=//' | sort -u)
ORPHAN_GENES=$(comm -23 <(echo "$GENE_IDS" | sed 's/ID=//') <(echo "$PARENT_IDS") | wc -l)
echo "  $ORPHAN_GENES gene records with no children"

echo ""
echo "End: $(date)"
echo "Done. Output: $OUT_GFF"
