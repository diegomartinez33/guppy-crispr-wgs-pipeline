#!/bin/bash
# Shared reference-genome version config, sourced by every pipeline script
# that needs a genome/annotation path or an output-directory suffix.
#
# Why this exists: GCF_000633615.1 (the reference used by this whole project
# since its start) was suppressed by NCBI; GCF_904066995.2 (male,
# PacBio+Hi-C, 2025) is now the official RefSeq reference genome for
# P. reticulata, with much better contiguity/completeness (see
# CLAUDE.md, "Migration to GCF_904066995.2 (v2)"). Rather than duplicating
# ~20 scripts, each one sources this file and gets REF/REF_GFF/INTERVALS/
# OUT_SUFFIX for whichever REF_VERSION it's run with.
#
# Usage: `source "${PROJECT_DIR}/codes/genome_versions.sh"` after PROJECT_DIR
# is set. Select the version via `REF_VERSION=v2 sbatch some_script.sh`
# (SLURM propagates exported env vars to the job by default). Leaving
# REF_VERSION unset defaults to v1 - byte-for-byte identical paths to
# before this file existed, so nothing about the existing v1 pipeline or
# its outputs changes.

REF_VERSION=${REF_VERSION:-v1}

case "$REF_VERSION" in
    v1)
        # Trinidad/Guanapo, female, short-read (2014) - suppressed by NCBI
        # as of 2026, kept as-is for continuity/comparison with prior results.
        REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
        REF_GFF=${PROJECT_DIR}/reference/GCF_000633615.1_annotation.gff
        INTERVALS=${PROJECT_DIR}/reference/intervals.list
        OUT_SUFFIX=""
        ;;
    v2)
        # P_reticulata-male-v2, male, PacBio+Hi-C (2025) - current NCBI
        # RefSeq reference genome for the species.
        REF=${PROJECT_DIR}/reference/GCF_904066995.2_P_reticulata-male-v2_genomic.fna
        REF_GFF=${PROJECT_DIR}/reference/GCF_904066995.2_annotation.gff
        INTERVALS=${PROJECT_DIR}/reference/intervals_v2.list
        OUT_SUFFIX="_v2"
        ;;
    *)
        echo "ERROR: REF_VERSION must be 'v1' or 'v2' (got '${REF_VERSION}')" >&2
        exit 1
        ;;
esac

echo "genome_versions.sh: REF_VERSION=${REF_VERSION} REF=${REF}"
