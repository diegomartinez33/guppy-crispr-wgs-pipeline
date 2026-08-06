#!/bin/bash
#SBATCH --job-name=crispresso_pooled_groups
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/crispresso_pooled_groups.out
#SBATCH --error=logs/crispresso_pooled_groups.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
ONTARGET_DIR=${PROJECT_DIR}/crispresso/ontarget/trimmomatic
OUTPUT_DIR=${PROJECT_DIR}/crispresso/pooled/trimmomatic
mkdir -p "$OUTPUT_DIR" logs/

# mamba shell hook fails in SLURM — use conda directly (see "CRISPResso2
# Activation in SLURM" Known Issue in CLAUDE.md)
CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate crispresso2_env
export PATH="${CONDA_BASE}/envs/crispresso2_env/bin:$PATH"

echo "Start time: $(date)"

# ── Pool Control ──────────────────────────────────────────────────────────────
CRISPRessoPooled \
  --crispresso_output_folders \
    "${ONTARGET_DIR}/Control_MNP_I_S54_L002/CRISPResso_on_Control_MNP_I_S54_L002" \
    "${ONTARGET_DIR}/Control_MNP_II_S55_L002/CRISPResso_on_Control_MNP_II_S55_L002" \
    "${ONTARGET_DIR}/Control_MNP_III_S56_L002/CRISPResso_on_Control_MNP_III_S56_L002" \
  --output_folder "${OUTPUT_DIR}/Control_pooled" \
  --name "Control_pooled" \
  --n_processes "${SLURM_CPUS_PER_TASK}"

echo "✅ Control pooled done"

# ── Pool RNP_Cas ──────────────────────────────────────────────────────────────
CRISPRessoPooled \
  --crispresso_output_folders \
    "${ONTARGET_DIR}/RNP_Cas1_S65_L002/CRISPResso_on_RNP_Cas1_S65_L002" \
    "${ONTARGET_DIR}/RNP_Cas2_S66_L002/CRISPResso_on_RNP_Cas2_S66_L002" \
    "${ONTARGET_DIR}/RNP_Cas3_S67_L002/CRISPResso_on_RNP_Cas3_S67_L002" \
    "${ONTARGET_DIR}/RNP_Cas4_S68_L002/CRISPResso_on_RNP_Cas4_S68_L002" \
  --output_folder "${OUTPUT_DIR}/RNP_Cas_pooled" \
  --name "RNP_Cas_pooled" \
  --n_processes "${SLURM_CPUS_PER_TASK}"

echo "✅ RNP_Cas pooled done"

# ── Pool Plasmid_Ko ───────────────────────────────────────────────────────────
CRISPRessoPooled \
  --crispresso_output_folders \
    "${ONTARGET_DIR}/Plasmid_Ko_P1_S61_L002/CRISPResso_on_Plasmid_Ko_P1_S61_L002" \
    "${ONTARGET_DIR}/Plasmid_Ko_P2_S62_L002/CRISPResso_on_Plasmid_Ko_P2_S62_L002" \
    "${ONTARGET_DIR}/Plasmid_Ko_P3_S63_L002/CRISPResso_on_Plasmid_Ko_P3_S63_L002" \
    "${ONTARGET_DIR}/Plasmid_Ko_P4_S64_L002/CRISPResso_on_Plasmid_Ko_P4_S64_L002" \
  --output_folder "${OUTPUT_DIR}/Plasmid_Ko_pooled" \
  --name "Plasmid_Ko_pooled" \
  --n_processes "${SLURM_CPUS_PER_TASK}"

echo "✅ Plasmid_Ko pooled done"

# ── Pool Only_MNP ─────────────────────────────────────────────────────────────
CRISPRessoPooled \
  --crispresso_output_folders \
    "${ONTARGET_DIR}/Only_MNP_C1_S57_L002/CRISPResso_on_Only_MNP_C1_S57_L002" \
    "${ONTARGET_DIR}/Only_MNP_C2_S58_L002/CRISPResso_on_Only_MNP_C2_S58_L002" \
    "${ONTARGET_DIR}/Only_MNP_C3_S59_L002/CRISPResso_on_Only_MNP_C3_S59_L002" \
    "${ONTARGET_DIR}/Only_MNP_C4_S60_L002/CRISPResso_on_Only_MNP_C4_S60_L002" \
  --output_folder "${OUTPUT_DIR}/Only_MNP_pooled" \
  --name "Only_MNP_pooled" \
  --n_processes "${SLURM_CPUS_PER_TASK}"

echo "✅ Only_MNP pooled done"
echo "End time: $(date)"