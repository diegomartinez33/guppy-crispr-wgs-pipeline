#!/bin/bash
#SBATCH --job-name=nextpolish_gapfilled
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --partition=medium
#SBATCH --output=logs/nextpolish_gapfilled_%j.out
#SBATCH --error=logs/nextpolish_gapfilled_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Targeted post-gap-fill polishing (see "Targeted Post-Gap-Fill Polishing"
# Known Issue in CLAUDE.md, and the matching section in
# reference/colombian_scaffolded_genome/README.md).
#
# Rationale: TGS-GapCloser's gap-filling (2026-08-31) gave a huge
# completeness gain (genome fraction 82.8%->92.1%, BUSCO 87.1%->95.4%) but
# also a real precision cost (misassemblies 7,589->23,532, NA50
# 198K->115K) - because the newly-filled sequence came from Nanopore
# reads, error-corrected only with racon, and was NEVER reconciled against
# the high-precision Illumina reads.
#
# The original whole-genome NextPolish experiment (see "NextPolish — No
# Improvement" Known Issue) found nothing to fix because it reprocessed
# reads the assembly was ALREADY built from. This run is different in
# effect even though it's the same whole-genome procedure: run on the
# GAP-FILLED genome, regions that were already in the pre-gap-fill genome
# should see little-to-no change (as before), while the Nanopore-derived
# filled regions - genuinely never compared against Illumina data - are
# where real, correctable disagreements should exist this time.
#
# Identical tool/config/2-round recipe as nextpolish_genome.sh (see that
# script for the full -N / forkserver-fix / resource-sizing rationale) -
# only GENOME and WORKDIR differ. Kept as a separate script (not a
# parameterized shared one) so both runs remain independently reproducible
# and their logs/configs don't collide.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
INPUT_DIR=${PROJECT_DIR}/trimmed_trimmomatic
GENOME=${PROJECT_DIR}/assembly/tgsgapcloser_output/colombian_gapfilled.fasta
WORKDIR=${PROJECT_DIR}/assembly/nextpolish_output_gapfilled

mkdir -p "$WORKDIR" logs/

if [ ! -f "$GENOME" ]; then
    echo "ERROR: genome not found at $GENOME"
    exit 1
fi

SAMPLES="Control_MNP_I_S54_L002 Control_MNP_II_S55_L002 Control_MNP_III_S56_L002"
SGS_FOFN=${WORKDIR}/sgs.fofn
: > "$SGS_FOFN"
for S in $SAMPLES; do
    for R in R1 R2; do
        F="${INPUT_DIR}/${S}_${R}_paired.fastq.gz"
        if [ ! -f "$F" ]; then
            echo "ERROR: missing input file: $F"
            exit 1
        fi
        echo "$F" >> "$SGS_FOFN"
    done
done
echo "sgs_fofn written: ${SGS_FOFN}"
cat "$SGS_FOFN"

RUN_CFG=${WORKDIR}/run.cfg
cat > "$RUN_CFG" <<EOF
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 4
multithread_jobs = 8
genome = ${GENOME}
workdir = ${WORKDIR}
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ${SGS_FOFN}
sgs_options = -max_depth 100 -N
EOF
echo "run.cfg written: ${RUN_CFG}"
cat "$RUN_CFG"

CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate nextpolish_env

echo "Start time: $(date)"
nextPolish "$RUN_CFG"
echo "End time: $(date)"

POLISHED="${WORKDIR}/genome.nextpolish.fasta"
if [ -f "$POLISHED" ]; then
    echo "Polished genome: $(grep -c '^>' "$POLISHED") sequences"
    ls -lh "$POLISHED" "${POLISHED}.stat" 2>/dev/null
else
    echo "ERROR: genome.nextpolish.fasta not produced"
    exit 1
fi
