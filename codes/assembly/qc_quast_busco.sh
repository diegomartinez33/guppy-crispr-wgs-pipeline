#!/bin/bash
#SBATCH --job-name=assembly_qc
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=short
#SBATCH --output=logs/assembly_qc_%j.out
#SBATCH --error=logs/assembly_qc_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Phase 3: Assembly quality assessment (QUAST + BUSCO).
# Uses actinopterygii_odb12.2 (actinopterygii_odb10 no longer exists in the
# BUSCO 5.7.1 catalog). Requires the lineage dataset pre-staged by
# 00_download_busco_lineage.sh so this can run fully --offline.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
SCAFFOLD=${PROJECT_DIR}/assembly/ragtag_output/ragtag.scaffold.fasta
QC_DIR=${PROJECT_DIR}/assembly/qc_results
BUSCO_DL=${QC_DIR}/busco_downloads

mkdir -p "${QC_DIR}/quast" "${QC_DIR}/busco" logs/

if [ ! -f "$SCAFFOLD" ]; then
    echo "ERROR: ragtag.scaffold.fasta not found at $SCAFFOLD - did Phase 2 (ragtag_scaffold.sh) finish?"
    exit 1
fi

echo "=== QUAST ==="
module load quast/5.0.2
echo "Start time: $(date)"
quast.py "$SCAFFOLD" -r "$REF" -o "${QC_DIR}/quast" -t "${SLURM_CPUS_PER_TASK}"
echo "QUAST end time: $(date)"

echo "=== BUSCO ==="
module load busco/5.7.1
cd "${QC_DIR}/busco"
busco -i "$SCAFFOLD" \
      -l actinopterygii_odb12.2 \
      -m genome \
      -c "${SLURM_CPUS_PER_TASK}" \
      --download_path "$BUSCO_DL" \
      --offline \
      -o busco_colombian_scaffold
echo "BUSCO end time: $(date)"

echo "=== Summary ==="
cat "${QC_DIR}/quast/report.txt" 2>/dev/null | head -30
find "${QC_DIR}/busco/busco_colombian_scaffold" -name "short_summary*.txt" -exec cat {} \;
