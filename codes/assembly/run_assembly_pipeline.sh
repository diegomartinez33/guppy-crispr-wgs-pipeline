#!/bin/bash
# Colombian guppy de novo assembly pipeline
# Chain: SPAdes co-assembly → RagTag scaffold → {QUAST+BUSCO QC, Liftoff annotation}
#
# PREREQUISITES - run these once on the login node before this driver:
#   ./00_download_gff.sh
#   ./00_download_busco_lineage.sh
#   ./00_setup_liftoff_env.sh

mkdir -p logs/

cd "$(dirname "$0")"

SPADES_JID=$(sbatch --parsable spades_coassembly.sh)
echo "SPAdes co-assembly job ID: ${SPADES_JID}  (bigmem partition, up to 10 days)"

RAGTAG_JID=$(sbatch --parsable \
    --dependency=afterok:${SPADES_JID} \
    ragtag_scaffold.sh)
echo "RagTag scaffold job ID: ${RAGTAG_JID}"

QC_JID=$(sbatch --parsable \
    --dependency=afterok:${RAGTAG_JID} \
    qc_quast_busco.sh)
echo "QUAST+BUSCO QC job ID: ${QC_JID}"

LIFTOFF_JID=$(sbatch --parsable \
    --dependency=afterok:${RAGTAG_JID} \
    liftoff_annotation.sh)
echo "Liftoff annotation job ID: ${LIFTOFF_JID}"

echo ""
echo "Pipeline: ${SPADES_JID} → ${RAGTAG_JID} → ${QC_JID} + ${LIFTOFF_JID}"
