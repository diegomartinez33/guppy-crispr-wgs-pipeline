#!/bin/bash
# Colombian guppy de novo assembly pipeline
# Chain: SPAdes co-assembly → RagTag scaffold → {QUAST QC, BUSCO QC, Liftoff annotation}
# (QUAST and BUSCO run as independent parallel jobs - see quast_qc.sh /
# busco_qc.sh, split from an original combined script after a slow QUAST
# run ate into BUSCO's shared time budget - see CLAUDE.md Known Issues)
#
# NOTE: contigs.fasta must be filtered to >=500bp (contigs.min500.fasta)
# before ragtag_scaffold.sh runs - RagTag cannot handle the raw ~3.4M
# unfiltered SPAdes contigs (see "RagTag Hang on 3.4M Unfiltered Contigs" in
# CLAUDE.md). This driver does not do that filtering step automatically;
# run it manually after spades_coassembly.sh completes:
#   awk '/^>/{split($0,p,"_"); keep=(p[4]+0>=500)} keep' contigs.fasta > contigs.min500.fasta
#
# PREREQUISITES - run these once on the login node before this driver:
#   ./00_download_gff.sh
#   ./00_download_busco_lineage.sh
#   ./00_setup_liftoff_env.sh

mkdir -p logs/

cd "$(dirname "$0")"

SPADES_JID=$(sbatch --parsable spades_coassembly.sh)
echo "SPAdes co-assembly job ID: ${SPADES_JID}  (bigmem partition, up to 10 days)"
echo "NOTE: after this completes, manually filter contigs.fasta to >=500bp"
echo "      (see comment above) before submitting ragtag_scaffold.sh"

RAGTAG_JID=$(sbatch --parsable \
    --dependency=afterok:${SPADES_JID} \
    ragtag_scaffold.sh)
echo "RagTag scaffold job ID: ${RAGTAG_JID}"

QUAST_JID=$(sbatch --parsable \
    --dependency=afterok:${RAGTAG_JID} \
    quast_qc.sh)
echo "QUAST QC job ID: ${QUAST_JID}"

BUSCO_JID=$(sbatch --parsable \
    --dependency=afterok:${RAGTAG_JID} \
    busco_qc.sh)
echo "BUSCO QC job ID: ${BUSCO_JID}"

LIFTOFF_JID=$(sbatch --parsable \
    --dependency=afterok:${RAGTAG_JID} \
    liftoff_annotation.sh)
echo "Liftoff annotation job ID: ${LIFTOFF_JID}"

echo ""
echo "Pipeline: ${SPADES_JID} → ${RAGTAG_JID} → ${QUAST_JID} + ${BUSCO_JID} + ${LIFTOFF_JID}"
