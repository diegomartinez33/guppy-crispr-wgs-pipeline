#!/bin/bash
#SBATCH --job-name=crispor_add_genomes
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=short
#SBATCH --output=logs/crispor_add_genomes_%j.out
#SBATCH --error=logs/crispor_add_genomes_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL
set -euo pipefail

# Register the reference (Trinidad) and Colombian pseudogenome as custom
# CRISPOR genomes, so crispor.py can be used for real efficiency (Doench/
# Azimuth) and specificity (MIT/CFD) guide scoring - complementing (not
# replacing) the manual PAM scan + variant classification in
# ko_guide_scan.py, which stays as the primary variant-vs-guide-window
# comparison tool.
#
# Uses the CRISPOR docker image via Singularity (Docker itself needs the
# docker group / root-equivalent daemon access, not available to regular
# users on this cluster - confirmed via `docker ps` -> permission denied).
# IMPORTANT: the "latest" tag on Docker Hub is a broken multi-arch build
# (only has an arm64 variant, pushed 2026-03-01 - probably from an Apple
# Silicon Mac) and will not run on this x86_64 cluster. Use v5.2c instead,
# confirmed to have both amd64 and arm64 variants.
#
# CRISPOR's internal genome directory (/data/www/crispor/genomes, and
# /data/www/genomes, both symlink to /data/genomes) is empty inside the
# read-only container image - bind-mount a persistent external directory
# there so genomes added here survive across container invocations
# (this same bind mount must be used again whenever crispor.py is run
# later to actually score guides).
#
# NO --gff here (deliberately, first attempt 2026-09-05 included --gff and
# looked COMPLETED in squeue/sacct, but was a silent failure): the
# container's bundled UCSC binaries (bedSort/bedToExons/genePredToBed) are
# missing libpng12.so.0 (stale Ubuntu base image dependency) and crash
# with an AssertionError while building the gene-locus annotation track.
# crisporAddGenome builds everything under /tmp/<genomeId>/ first and only
# copies it into --baseDir at the very end, so that crash meant NOTHING
# was ever written to the bind-mounted genomes dir, even though the
# outer script (no `set -e` in the first version) exited 0 and SLURM
# reported COMPLETED. Registering fasta-only avoids the broken GFF/exon
# conversion step entirely - it's a fine tradeoff here since ko_guide_scan.py
# only reads targetSeq/mitSpecScore/offtargetCount/efficiency scores from
# crispor.py's output, never targetGenomeGeneLocus (which needs the gene
# track). `set -euo pipefail` above now makes this script actually fail
# (non-zero exit, SLURM state FAILED) if crisporAddGenome errors again.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SIF=${PROJECT_DIR}/codes/analysis/crispor_singularity/crispor_v5.2c_amd64.sif
GENOMES_DIR=${PROJECT_DIR}/codes/analysis/crispor_singularity/genomes

mkdir -p "$GENOMES_DIR" logs/

module load singularity/3.7.1

# crisporAddGenome stages everything under /tmp/<genomeId>/ before copying
# the finished result into --baseDir, and refuses to run if that staging
# dir already exists (guards against concurrent runs). /tmp is node-local
# on this cluster, so a crashed prior run (e.g. the first attempt here,
# which got through BWA indexing before crashing on the GFF step) leaves
# an orphaned staging dir that persists on whichever node this job lands
# on next, regardless of anything cleaned up from the login node. Clear it
# unconditionally at the start so this script is idempotent no matter
# which node picks it up.
rm -rf /tmp/guppyRefTrinidad /tmp/guppyColPseudogenome

echo "Start time: $(date)"

echo "=== Adding reference genome (GCF_000633615.1, Trinidad/Guanapo) ==="
singularity exec -B "${GENOMES_DIR}:/data/genomes" "$SIF" \
    /data/www/crispor/tools/crisporAddGenome fasta \
    "${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna" \
    --desc 'guppyRefTrinidad|Poecilia reticulata|Guppy (Trinidad/Guanapo reference)|GCF_000633615.1' \
    --baseDir /data/genomes \
    -f

echo ""
echo "=== Adding Colombian pseudogenome ==="
singularity exec -B "${GENOMES_DIR}:/data/genomes" "$SIF" \
    /data/www/crispor/tools/crisporAddGenome fasta \
    "${PROJECT_DIR}/reference/pseudogenome/colombian_pseudogenome.fna" \
    --desc 'guppyColPseudogenome|Poecilia reticulata|Guppy (Colombian pseudogenome)|bcftools consensus 2026' \
    --baseDir /data/genomes \
    -f

echo ""
echo "End time: $(date)"
echo "=== Genomes directory contents ==="
find "$GENOMES_DIR" -maxdepth 2

echo ""
echo "=== Verifying both genomes actually persisted (2bit + BWA index) ==="
for g in guppyRefTrinidad guppyColPseudogenome; do
    if [ ! -f "${GENOMES_DIR}/${g}/${g}.2bit" ] || [ ! -f "${GENOMES_DIR}/${g}/${g}.fa.bwt" ]; then
        echo "ERROR: ${g} is missing .2bit or .fa.bwt - registration failed" >&2
        exit 1
    fi
    echo "OK: ${g}"
done
