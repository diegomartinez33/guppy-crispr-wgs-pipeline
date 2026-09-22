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

CONDA_BASE=${CONDA_BASE:-${HOME}/miniconda3}
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate liftoff_env

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.." >/dev/null 2>&1 && pwd)}}"
source "${PROJECT_DIR}/codes/genome_versions.sh"
TRINIDAD_GFF=${REF_GFF}
PSEUDO=${PROJECT_DIR}/reference/pseudogenome${OUT_SUFFIX}/colombian_pseudogenome.fna
OUT_GFF=${PROJECT_DIR}/reference/pseudogenome${OUT_SUFFIX}/colombian_pseudogenome.liftoff.gff3

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

# liftoff writes this into the CURRENT WORKING DIRECTORY, not next to
# $OUT_GFF - undocumented (found 2026-09-15 on the v2 run). Relocate it
# next to the annotation it describes, matching v1's naming convention.
if [ -f "unmapped_features.txt" ]; then
    mv "unmapped_features.txt" "$(dirname "$OUT_GFF")/liftoff_unmapped_genes.txt"
fi

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
# Match ID=/Parent= only at the start of the (tab-isolated) attributes
# column, not anywhere in the line - a naive `grep -oP 'ID=[^;]+'` over
# the whole line also matches substrings like `sequence_ID=` and
# `copy_num_ID=` (both real Liftoff attributes on partial/low-identity
# transfers), inflating the "gene ID" count 3x and producing a nonsense
# orphan count exceeding the total gene count. Found 2026-09-15 on the v2
# pseudogenome run (93,678 false "gene IDs" vs the real 31,226; reported
# 62,452 "orphans" when the true count, verified independently in Python,
# is 0) - v1's run didn't trigger this as visibly, but the bug was latent
# there too.
GENE_IDS=$(grep -v "^#" "$OUT_GFF" | awk -F'\t' '$3=="gene"{print $9}' | grep -oP '(?:^|;)ID=\K[^;]+' | sort)
PARENT_IDS=$(grep -v "^#" "$OUT_GFF" | awk -F'\t' '{print $9}' | grep -oP '(?:^|;)Parent=\K[^;]+' | sort -u)
ORPHAN_GENES=$(comm -23 <(echo "$GENE_IDS") <(echo "$PARENT_IDS") | wc -l)
echo "  $ORPHAN_GENES gene records with no children"

echo ""
echo "End: $(date)"
echo "Done. Output: $OUT_GFF"
