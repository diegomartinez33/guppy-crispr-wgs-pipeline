#!/bin/bash
#SBATCH --job-name=crispresso_aggregate_ontarget
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/crispresso_aggregate_ontarget.out
#SBATCH --error=logs/crispresso_aggregate_ontarget.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# Pending Analyses #1 Track C: CRISPRessoAggregate — all 15 samples summary
# (individual on-target results, not merged-by-group).
#
# Same two fixes as crispresso_wgs_aggregate.sh (see that script + CLAUDE.md
# Known Issues for full context):
#   1. CRISPRessoAggregate has no output-directory flag - must cd into
#      OUTPUT_DIR first, it always writes into the current working directory.
#   2. --suppress_report: this CRISPResso2 install (2.3.1) has a bug in
#      make_multi_report() (missing positional args) that crashes ONLY the
#      final HTML report step - the actual data tables/plots generate fine
#      before that point.
#
# -p prefix: each sample's on-target run is a single exact-match folder
# (crispresso/ontarget/trimmomatic/<SAMPLE>/CRISPResso_on_<SAMPLE>), not a
# glob over multiple per-site subfolders like the WGS case - so one -p per
# sample, pointing directly at that exact path.

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
ONTARGET_DIR=${PROJECT_DIR}/crispresso/ontarget/trimmomatic
OUTPUT_DIR=${PROJECT_DIR}/crispresso/aggregate

mkdir -p "$OUTPUT_DIR" logs/

CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate crispresso2_env
export PATH="${CONDA_BASE}/envs/crispresso2_env/bin:$PATH"

echo "Start time: $(date)"

cd "$OUTPUT_DIR"

SAMPLES="Control_MNP_I_S54_L002 Control_MNP_II_S55_L002 Control_MNP_III_S56_L002 Only_MNP_C1_S57_L002 Only_MNP_C2_S58_L002 Only_MNP_C3_S59_L002 Only_MNP_C4_S60_L002 Plasmid_Ko_P1_S61_L002 Plasmid_Ko_P2_S62_L002 Plasmid_Ko_P3_S63_L002 Plasmid_Ko_P4_S64_L002 RNP_Cas1_S65_L002 RNP_Cas2_S66_L002 RNP_Cas3_S67_L002 RNP_Cas4_S68_L002"

PREFIX_ARGS=""
for SAMPLE in $SAMPLES; do
    RUN_DIR="${ONTARGET_DIR}/${SAMPLE}/CRISPResso_on_${SAMPLE}"
    if [ ! -d "$RUN_DIR" ]; then
        echo "ERROR: missing on-target run: $RUN_DIR"
        exit 1
    fi
    PREFIX_ARGS="$PREFIX_ARGS -p ${RUN_DIR}"
done

CRISPRessoAggregate \
  $PREFIX_ARGS \
  --name "all_samples_ontarget_aggregate" \
  --min_reads_for_inclusion 10 \
  --place_report_in_output_folder \
  --suppress_report \
  --n_processes "${SLURM_CPUS_PER_TASK}"

echo "✅ all_samples_ontarget aggregated"
echo "End time: $(date)"
