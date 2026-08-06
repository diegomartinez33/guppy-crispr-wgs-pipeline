#!/bin/bash
#SBATCH --job-name=nextpolish_genome
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --partition=medium
#SBATCH --output=logs/nextpolish_genome_%j.out
#SBATCH --error=logs/nextpolish_genome_%j.err
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Phase 5: Short-read polishing of the scaffolded Colombian genome (NextPolish).
#
# Corrects small-scale errors (indels/SNPs) introduced during assembly by
# realigning the same Illumina short reads used to build the genome back
# onto it. Targets the ~6.9% of "complete" BUSCO genes found with internal
# stop codons (likely frameshift artifacts) - see README.md in
# reference/colombian_scaffolded_genome/ for the quality baseline this is
# meant to improve.
#
# Uses NextPolish 1.4.1 (conda env nextpolish_env, see 00_setup_nextpolish_env.sh).
# Default task=best with only sgs_fofn provided resolves to 2 rounds of
# short-read polishing (task 1,2,1,2) - NextPolish's own recommended recipe
# for Illumina-only polishing.
#
# Resource sizing: the project's own single-sample BWA alignment against a
# similarly-sized reference (bwa_trimmomatic_array.sh) used 8 CPUs/32GB/20h
# PER SAMPLE. NextPolish aligns all 3 Control samples, twice (2 rounds), so
# sized generously here (32 CPUs, 64GB, 3 days, medium partition) rather
# than risk a repeat of this project's history of first-attempt
# under-provisioning (see SPAdes/RagTag/QUAST Known Issues in CLAUDE.md).

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
INPUT_DIR=${PROJECT_DIR}/trimmed_trimmomatic
GENOME=${PROJECT_DIR}/reference/colombian_scaffolded_genome/colombian_scaffolded.fna
WORKDIR=${PROJECT_DIR}/assembly/nextpolish_output

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
sgs_options = -max_depth 100
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
