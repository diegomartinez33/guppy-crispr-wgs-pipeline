#!/bin/bash
# Copies the population-genome resources + the files needed to run the
# reusable CRISPR guide-design and primer-design/verification pipelines
# (ko_guide_scan.py, crispri_tss_scan.py, design_offtarget_primers.py,
# verify_rtqpcr_primers.py) to the shared guppy-genome@hypatia account, so
# colleagues can use them without needing the author's personal account.
#
# Deliberately excludes: raw_fastq/, trimmed_*/, mapping/, gatk/, crispresso*/,
# assembly/ intermediate work, and CRISPResso/GATK off-target analysis output
# - none of the 4 scripts above touch them, and replicating the WGS
# off-target pipeline itself needs the colleague's own samples anyway (see
# docs/PIPELINE.md / docs/TUTORIAL.md).
#
# Run this yourself with your own guppy-genome credentials/keys - see
# docs/CLUSTER_ACCESS.md. Safe to re-run: rsync skips files that are already
# up to date on the destination, so an interrupted run just resumes.

set -e

# No SLURM_SUBMIT_DIR tier here, unlike the sbatch-submitted scripts
# elsewhere in codes/ - this script is always run directly (`bash
# codes/copy_to_shared_account.sh`), never via sbatch, so trusting
# SLURM_SUBMIT_DIR would only pick up a stale value inherited from some
# unrelated earlier interactive SLURM session, silently overriding the
# correct self-location below (found 2026-09-22: a user with a lingering
# SLURM_SUBMIT_DIR from an interactive session launched inside codes/ hit
# exactly this, and PROJECT_DIR resolved one level too deep).
PROJECT_DIR="${PROJECT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")/.." >/dev/null 2>&1 && pwd)}"
cd "$PROJECT_DIR"

DEST_USER="guppy-genome"
DEST_HOST="hypatia.uniandes.edu.co"
DEST_BASE="~/CRISPRGuppy/guppy-crispr-wgs-pipeline"   # adjust if the shared account uses a different layout

# 11 separate rsync calls below - without connection reuse, that's 11
# separate password prompts (the shared account has no SSH key set up from
# this session, so it's password-only). SSH ControlMaster/ControlPersist
# multiplexes them: the first call authenticates and becomes the "master"
# connection, every later call to the same host reuses it automatically -
# one password prompt total, not eleven. ControlPersist=600 keeps the
# connection open for 10 minutes after the last rsync exits, then SSH
# closes it on its own (the explicit `-O exit` at the end also closes it
# immediately once this script finishes, whichever comes first).
SSH_CONTROL_PATH="/tmp/ssh-cm-${DEST_USER}-${DEST_HOST}-$$"
SSH_OPTS="-o ControlMaster=auto -o ControlPath=${SSH_CONTROL_PATH} -o ControlPersist=600"
RSYNC=(rsync -avP --info=progress2 -e "ssh ${SSH_OPTS}")

cleanup() {
    ssh -o ControlPath="${SSH_CONTROL_PATH}" -O exit "${DEST_USER}@${DEST_HOST}" 2>/dev/null || true
}
trap cleanup EXIT

copy_dir() {
    local src="$1" dst="$2"
    echo ""
    echo "=== ${src} -> ${DEST_USER}@${DEST_HOST}:${dst} ==="
    "${RSYNC[@]}" "${src}/" "${DEST_USER}@${DEST_HOST}:${dst}/"
}

copy_file() {
    local src="$1" dst_dir="$2"
    echo ""
    echo "=== ${src} -> ${DEST_USER}@${DEST_HOST}:${dst_dir}/ ==="
    "${RSYNC[@]}" "${src}" "${DEST_USER}@${DEST_HOST}:${dst_dir}/"
}

START=$(date +%s)

# ── Population genome resources ────────────────────────────────────────────
copy_dir "reference/pseudogenome"                    "${DEST_BASE}/reference/pseudogenome"
copy_dir "reference/pseudogenome_v2"                  "${DEST_BASE}/reference/pseudogenome_v2"
copy_dir "reference/colombian_scaffolded_genome"       "${DEST_BASE}/reference/colombian_scaffolded_genome"
copy_dir "reference/colombian_scaffolded_genome_v2"    "${DEST_BASE}/reference/colombian_scaffolded_genome_v2"

# ── Base reference genomes + annotation (needed by ko_guide_scan.py /
#    crispri_tss_scan.py / design_offtarget_primers.py / verify_rtqpcr_primers.py) ──
copy_file "reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna"      "${DEST_BASE}/reference"
copy_file "reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna.fai"  "${DEST_BASE}/reference"
copy_file "reference/GCF_000633615.1_annotation.gff"                       "${DEST_BASE}/reference"
copy_file "reference/GCF_904066995.2_P_reticulata-male-v2_genomic.fna"     "${DEST_BASE}/reference"
copy_file "reference/GCF_904066995.2_P_reticulata-male-v2_genomic.fna.fai" "${DEST_BASE}/reference"
copy_file "reference/GCF_904066995.2_annotation.gff"                       "${DEST_BASE}/reference"

# ── CRISPOR Singularity container + registered genomes (required for all
#    off-target scoring in guide design) ───────────────────────────────────
copy_dir "codes/analysis/crispor_singularity"          "${DEST_BASE}/codes/analysis/crispor_singularity"

END=$(date +%s)
echo ""
echo "=== Done in $(( (END-START)/60 ))m $(( (END-START)%60 ))s ==="
echo "Remember: colleagues still need to build their own primer3_env locally"
echo "via codes/analysis/setup_primer3.sh - it's a conda environment, not a"
echo "file, so it isn't copied here."
