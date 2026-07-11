#!/bin/bash
#SBATCH --job-name=ragtag_scaffold
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --output=logs/ragtag_scaffold_%j.out
#SBATCH --error=logs/ragtag_scaffold_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Phase 2: Reference-guided scaffolding of SPAdes contigs (RagTag).
#
# Uses contigs.min500.fasta (>=500bp), not the raw contigs.fasta: job 676655
# fed RagTag the full 3.4M-contig SPAdes output and it hung indefinitely at
# the "Writing scaffolds" step (never produced ragtag.scaffold.agp even after
# 4h, despite alignment+ordering finishing in under 4 minutes) - RagTag has
# no built-in filter for query sequence length, only for alignment length
# (-f), so pre-filtering the input is the fix. See "RagTag hang on 3.4M
# contigs" in Known Issues. Filtered set keeps 501,663 contigs (>=500bp) -
# the vast majority of real assembled sequence; the full unfiltered
# contigs.fasta is untouched on disk for any future use.
# --remove-small placement: filtering is applied to the INPUT fasta, not via
# -f/--remove-small (which only governs alignment trust, not which query
# sequences are processed).
# -C concatenates any remaining unplaced contigs into one chr0 record
# instead of writing each individually, as a second safeguard against the
# same class of hang even at 500K contigs.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
CONTIGS=${PROJECT_DIR}/assembly/spades_control_coassembly/contigs.min500.fasta
OUTPUT_DIR=${PROJECT_DIR}/assembly/ragtag_output

mkdir -p logs/

if [ ! -f "$CONTIGS" ]; then
    echo "ERROR: contigs.min500.fasta not found at $CONTIGS - did Phase 1 (spades_coassembly.sh) finish and get filtered?"
    exit 1
fi

module load ragtag/2.1.0

echo "Start time: $(date)"
ragtag.py scaffold "$REF" "$CONTIGS" -o "$OUTPUT_DIR" -t "${SLURM_CPUS_PER_TASK}" -C
echo "End time: $(date)"

SCAFFOLD="${OUTPUT_DIR}/ragtag.scaffold.fasta"
if [ -f "$SCAFFOLD" ]; then
    ls -lh "$SCAFFOLD"
else
    echo "ERROR: ragtag.scaffold.fasta not produced"
    exit 1
fi
