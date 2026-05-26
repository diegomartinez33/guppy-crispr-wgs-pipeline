#!/bin/bash
#SBATCH --job-name=crispresso_compare_groups
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/crispresso_compare_groups.out
#SBATCH --error=logs/crispresso_compare_groups.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
POOLED_DIR=${PROJECT_DIR}/crispresso/pooled/trimmomatic
OUTPUT_DIR=${PROJECT_DIR}/crispresso/compare/trimmomatic

mkdir -p "$OUTPUT_DIR" logs/

#source /hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso/etc/profile.d/conda.sh
#eval "$(mamba shell hook --shell bash)"
#mamba activate crispresso2_env
source /hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso/etc/profile.d/conda.sh
conda activate crispresso2_env

echo "Start time: $(date)"

# ── Control vs RNP_Cas ────────────────────────────────────────────────────────
CRISPRessoCompare \
  "${POOLED_DIR}/Control_pooled" \
  "${POOLED_DIR}/RNP_Cas_pooled" \
  --sample_1_name "Control" \
  --sample_2_name "RNP_Cas" \
  --output_folder "${OUTPUT_DIR}/RNP_Cas_vs_Control" \
  --name "RNP_Cas_vs_Control"

echo "✅ RNP_Cas vs Control done"

# ── Control vs Plasmid_Ko ─────────────────────────────────────────────────────
CRISPRessoCompare \
  "${POOLED_DIR}/Control_pooled" \
  "${POOLED_DIR}/Plasmid_Ko_pooled" \
  --sample_1_name "Control" \
  --sample_2_name "Plasmid_Ko" \
  --output_folder "${OUTPUT_DIR}/Plasmid_Ko_vs_Control" \
  --name "Plasmid_Ko_vs_Control"

echo "✅ Plasmid_Ko vs Control done"

# ── Control vs Only_MNP ───────────────────────────────────────────────────────
CRISPRessoCompare \
  "${POOLED_DIR}/Control_pooled" \
  "${POOLED_DIR}/Only_MNP_pooled" \
  --sample_1_name "Control" \
  --sample_2_name "Only_MNP" \
  --output_folder "${OUTPUT_DIR}/Only_MNP_vs_Control" \
  --name "Only_MNP_vs_Control"

echo "✅ Only_MNP vs Control done"

# ── RNP_Cas vs Plasmid_Ko (comparación entre tratamientos) ───────────────────
CRISPRessoCompare \
  "${POOLED_DIR}/RNP_Cas_pooled" \
  "${POOLED_DIR}/Plasmid_Ko_pooled" \
  --sample_1_name "RNP_Cas" \
  --sample_2_name "Plasmid_Ko" \
  --output_folder "${OUTPUT_DIR}/RNP_Cas_vs_Plasmid_Ko" \
  --name "RNP_Cas_vs_Plasmid_Ko"

echo "✅ RNP_Cas vs Plasmid_Ko done"

# ── RNP_Cas vs Only_MNP (comparación entre tratamientos) ───────────────────
CRISPRessoCompare \
  "${POOLED_DIR}/RNP_Cas_pooled" \
  "${POOLED_DIR}/Only_MNP_pooled" \
  --sample_1_name "RNP_Cas" \
  --sample_2_name "Only_MNP" \
  --output_folder "${OUTPUT_DIR}/RNP_Cas_vs_Only_MNP" \
  --name "RNP_Cas_vs_Only_MNP"

echo "✅ RNP_Cas vs Only_MNP done"

# ── Plasmid_Ko vs Only_MNP (comparación entre tratamientos) ───────────────────
CRISPRessoCompare \
  "${POOLED_DIR}/Plasmid_Ko_pooled" \
  "${POOLED_DIR}/Only_MNP_pooled" \
  --sample_1_name "Plasmid_Ko" \
  --sample_2_name "Only_MNP" \
  --output_folder "${OUTPUT_DIR}/Plasmid_Ko_vs_Only_MNP" \
  --name "Plasmid_Ko_vs_Only_MNP"

echo "✅ Plasmid_Ko vs Only_MNP done"
echo "End time: $(date)"