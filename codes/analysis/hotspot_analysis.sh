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

# This script's original `conda activate fastp_env` never actually worked:
# `conda info --base` resolves to whatever ~/.bashrc's conda-init block
# points at (miniconda3_crispresso, which has no fastp_env at all), and even
# the real fastp_env (under anaconda3/) is missing statsmodels/scipy, which
# hotspot_analysis.py needs for its BH-FDR hotspot calling. Both v1's
# original runs (jobs 683994/684016) and this script therefore always
# silently fell through to whatever plain `python3`/`module load` left on
# PATH - which happened to be the cluster's anaconda/conda4.12.0 module
# (confirmed: it has pandas/matplotlib/statsmodels/scipy all present, and
# reproduces v1's exact "1,780 windows, BH-corrected FDR<0.05, 403 merged
# regions" result). Load that module explicitly instead of depending on
# fragile ambient PATH/conda state - found 2026-09-15 on the v2 run.
module load anaconda/conda4.12.0

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data

python "${PROJECT_DIR}/codes/analysis/hotspot_analysis.py"
