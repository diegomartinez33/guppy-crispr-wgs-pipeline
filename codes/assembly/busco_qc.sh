#!/bin/bash
#SBATCH --job-name=busco_qc
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --output=logs/busco_qc_%j.out
#SBATCH --error=logs/busco_qc_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Phase 3b: Assembly completeness assessment (BUSCO).
#
# Split out from the original combined qc_quast_busco.sh so a slow QUAST
# run (see quast_qc.sh) can't block or eat into BUSCO's time budget - both
# now depend only on RagTag (Phase 2) and run independently/in parallel.
#
# Uses actinopterygii_odb10 (2024-01-08 archived dataset). BUSCO 5.7.1's
# own code (BuscoConfig.py) HARDCODES a requirement for datasets_version ==
# "odb10" - it fatally errors on any other version (e.g. odb12.2, the
# current interactive-catalog default) regardless of the --datasets_version
# flag; the version check at BuscoConfig.py:673 is unconditional. odb10 is
# delisted from `busco --list-datasets` but the tarball is still hosted -
# see 00_download_busco_lineage.sh and "BUSCO — hardcoded odb10 version
# check" in CLAUDE.md Known Issues for the full debugging history.
# Requires the lineage dataset pre-staged by 00_download_busco_lineage.sh
# so this can run fully --offline.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SCAFFOLD=${PROJECT_DIR}/assembly/ragtag_output/ragtag.scaffold.fasta
QC_DIR=${PROJECT_DIR}/assembly/qc_results
BUSCO_DL=${QC_DIR}/busco_downloads

mkdir -p "${QC_DIR}/busco" logs/

if [ ! -f "$SCAFFOLD" ]; then
    echo "ERROR: ragtag.scaffold.fasta not found at $SCAFFOLD - did Phase 2 (ragtag_scaffold.sh) finish?"
    exit 1
fi

module load busco/5.7.1
cd "${QC_DIR}/busco"

echo "Start time: $(date)"
busco -i "$SCAFFOLD" \
      -l actinopterygii_odb10 \
      -m genome \
      -c "${SLURM_CPUS_PER_TASK}" \
      --download_path "$BUSCO_DL" \
      --offline \
      -o busco_colombian_scaffold
echo "End time: $(date)"

echo "=== Summary ==="
find "${QC_DIR}/busco/busco_colombian_scaffold" -name "short_summary*.txt" -exec cat {} \;
