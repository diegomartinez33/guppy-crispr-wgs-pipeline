#!/bin/bash
#SBATCH --job-name=hotspot_analysis
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/hotspot_analysis_%j.out
#SBATCH --error=logs/hotspot_analysis_%j.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

module load bcftools/1.15.1
module load bedtools/2.30.0

source $(conda info --base)/etc/profile.d/conda.sh
conda activate fastp_env

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data

python "${PROJECT_DIR}/codes/analysis/hotspot_analysis.py"
