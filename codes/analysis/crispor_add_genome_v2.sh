#!/bin/bash
#SBATCH --job-name=crispor_add_genome_v2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=short
#SBATCH --output=logs/crispor_add_genome_v2_%j.out
#SBATCH --error=logs/crispor_add_genome_v2_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL
set -euo pipefail

# Register the new reference genome (GCF_904066995.2, male, PacBio+Hi-C) as
# a custom CRISPOR genome ("guppyRefMaleV2"), so crispor.py can be used for
# real off-target discovery/scoring against it - same pattern and same
# container as crispor_add_genomes.sh (v1 genomes), see that script's
# comments for the full history of bugs found while setting this up
# (Docker unusable on this cluster -> Singularity; "latest" tag broken
# arm64-only -> use v5.2c; --gff crashes on missing libpng12.so.0 in the
# container's UCSC binaries -> register fasta-only).
#
# Needed now (ahead of the rest of the v2 pseudogenome/ko_guide_scan work)
# because the off-target discovery step for the bdnf sgRNA (previously done
# by manually downloading a .xls from crispor.org against the OLD genome)
# is being modernized to run crispor.py directly against this registered
# genome instead - see convert_crispor_offtargets.py and CLAUDE.md,
# "Migration to GCF_904066995.2 (v2)".

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SIF=${PROJECT_DIR}/codes/analysis/crispor_singularity/crispor_v5.2c_amd64.sif
GENOMES_DIR=${PROJECT_DIR}/codes/analysis/crispor_singularity/genomes

mkdir -p "$GENOMES_DIR" logs/

module load singularity/3.7.1

# /tmp is node-local and persists across job submissions on the same node -
# clear any orphaned staging dir from a crashed prior run before starting.
rm -rf /tmp/guppyRefMaleV2

echo "Start time: $(date)"

echo "=== Adding new reference genome (GCF_904066995.2, male, PacBio+Hi-C) ==="
singularity exec -B "${GENOMES_DIR}:/data/genomes" "$SIF" \
    /data/www/crispor/tools/crisporAddGenome fasta \
    "${PROJECT_DIR}/reference/GCF_904066995.2_P_reticulata-male-v2_genomic.fna" \
    --desc 'guppyRefMaleV2|Poecilia reticulata|Guppy (male reference, PacBio+Hi-C)|GCF_904066995.2' \
    --baseDir /data/genomes \
    -f

echo ""
echo "End time: $(date)"
echo "=== Genomes directory contents ==="
find "$GENOMES_DIR" -maxdepth 2

echo ""
echo "=== Verifying genome actually persisted (2bit + BWA index) ==="
if [ ! -f "${GENOMES_DIR}/guppyRefMaleV2/guppyRefMaleV2.2bit" ] || [ ! -f "${GENOMES_DIR}/guppyRefMaleV2/guppyRefMaleV2.fa.bwt" ]; then
    echo "ERROR: guppyRefMaleV2 is missing .2bit or .fa.bwt - registration failed" >&2
    exit 1
fi
echo "OK: guppyRefMaleV2"
