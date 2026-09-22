#!/bin/bash
#SBATCH --job-name=liftoff_annotation
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --output=logs/liftoff_annotation_%j.out
#SBATCH --error=logs/liftoff_annotation_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Phase 4: Annotation transfer (Liftoff) - maps reference gene models,
# including bdnf, onto the new Colombian scaffolded genome.
# Requires 00_download_gff.sh (or _v2 variant) and 00_setup_liftoff_env.sh
# to have been run.
#
# GENOME_STAGE selects which assembly stage to annotate + promote to the
# canonical colombian_scaffolded.fna:
#   scaffold (default) - the raw RagTag output (assembly/ragtag_output${OUT_SUFFIX}/
#             ragtag.scaffold.fasta) - the original 3-step SPAdes->RagTag->
#             Liftoff pipeline described in reference/colombian_scaffolded_genome/README.md.
#   final   - the gap-filled + polished genome (assembly/nextpolish_output_gapfilled${OUT_SUFFIX}/
#             genome.nextpolish.fasta) - the actually-validated, adopted
#             assembly (see CLAUDE.md, "Targeted Post-Gap-Fill Polishing").
# v1's canonical colombian_scaffolded.fna was stuck at the "scaffold" stage
# from Jul 11 despite the "final" stage being validated back in
# early September - the "final" genome was computed but never promoted
# here, so its Liftoff annotation/indices never existed. Found + fixed
# 2026-09-17 while scoping the v2 de novo assembly work - the old
# scaffold-stage files were archived to pre_gapfill_archive/, not deleted.

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.." >/dev/null 2>&1 && pwd)}}"
source "${PROJECT_DIR}/codes/genome_versions.sh"
GENOME_STAGE=${GENOME_STAGE:-scaffold}
case "$GENOME_STAGE" in
    scaffold) NEW_FASTA=${PROJECT_DIR}/assembly/ragtag_output${OUT_SUFFIX}/ragtag.scaffold.fasta ;;
    final)    NEW_FASTA=${PROJECT_DIR}/assembly/nextpolish_output_gapfilled${OUT_SUFFIX}/genome.nextpolish.fasta ;;
    *) echo "ERROR: GENOME_STAGE must be 'scaffold' or 'final' (got '${GENOME_STAGE}')" >&2; exit 1 ;;
esac
OUT_DIR=${PROJECT_DIR}/reference/colombian_scaffolded_genome${OUT_SUFFIX}
NEW_GFF=${OUT_DIR}/colombian_scaffolded.liftoff.gff3

mkdir -p "$OUT_DIR" logs/

if [ ! -f "$REF_GFF" ]; then
    echo "ERROR: reference GFF not found at $REF_GFF - run 00_download_gff.sh (or _v2) first"
    exit 1
fi
if [ ! -f "$NEW_FASTA" ]; then
    echo "ERROR: genome not found at $NEW_FASTA (stage=${GENOME_STAGE}) - did the prior stage finish?"
    exit 1
fi

CONDA_BASE=${CONDA_BASE:-${HOME}/miniconda3}
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate liftoff_env

cp "$NEW_FASTA" "${OUT_DIR}/colombian_scaffolded.fna"

echo "Start time: $(date)"
echo "REF_VERSION=${REF_VERSION}  GENOME_STAGE=${GENOME_STAGE}  NEW_FASTA=${NEW_FASTA}"
liftoff -g "$REF_GFF" -o "$NEW_GFF" -p "${SLURM_CPUS_PER_TASK}" "$NEW_FASTA" "$REF"
echo "End time: $(date)"

# liftoff writes this into the CURRENT WORKING DIRECTORY, not next to
# $NEW_GFF - same undocumented quirk already fixed in liftoff_pseudogenome.sh.
if [ -f "unmapped_features.txt" ]; then
    mv "unmapped_features.txt" "${OUT_DIR}/liftoff_unmapped_genes.txt"
fi

echo "=== bdnf locus check ==="
if grep -q "ID=gene-bdnf" "$NEW_GFF"; then
    echo "bdnf gene successfully lifted over:"
    grep "ID=gene-bdnf" "$NEW_GFF"
else
    echo "WARNING: bdnf gene not found in lifted annotation"
fi
