#!/bin/bash
#SBATCH --job-name=spades_coassembly_resume
#SBATCH --cpus-per-task=16
#SBATCH --mem=470G
#SBATCH --time=10-00:00:00
#SBATCH --partition=bigmem
#SBATCH --output=logs/spades_coassembly_resume_%j.out
#SBATCH --error=logs/spades_coassembly_resume_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# One-time recovery script for job 652298, which ran ~24h and reached the
# K77 Distance Estimation stage before a genuine OOM (mimalloc ENOMEM,
# peak RSS ~249GB against the old 250G limit). The last checkpoint saved
# was "late_pair_info_count" within K77, 10 minutes before the crash.
#
# --restart-from last resumes from that checkpoint with updated options
# (here: a higher -m), skipping error correction + K21/K33/K55 + most of
# K77 (~23 of the 24 hours already spent) and only re-running Distance
# Estimation onward. See "SPAdes — mmap ENOMEM" Known Issue in CLAUDE.md.
#
# This script is NOT part of the standing pipeline (spades_coassembly.sh
# already has -m 450 for any future from-scratch run) - it exists only to
# recover this specific interrupted output directory.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
OUTPUT_DIR=${PROJECT_DIR}/assembly/spades_control_coassembly

mkdir -p logs/

if [ ! -d "${OUTPUT_DIR}/K77/saves" ]; then
    echo "ERROR: no checkpoint saves found at ${OUTPUT_DIR}/K77/saves — nothing to resume"
    exit 1
fi

module load spades/4.0.0

echo "Start time: $(date)"
echo "Available space: $(df -h "$PROJECT_DIR" | awk 'NR==2 {print $4}')"

spades.py \
    --restart-from last \
    -m 450 \
    -o "$OUTPUT_DIR"

echo "End time: $(date)"

CONTIGS="${OUTPUT_DIR}/contigs.fasta"
if [ -f "$CONTIGS" ]; then
    echo "Contigs: $(grep -c '^>' "$CONTIGS") sequences"
    ls -lh "$CONTIGS"
else
    echo "ERROR: contigs.fasta not produced"
    exit 1
fi
