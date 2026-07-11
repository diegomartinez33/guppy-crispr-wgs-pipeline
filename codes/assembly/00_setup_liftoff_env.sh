#!/bin/bash
# One-time setup: create a dedicated conda environment for Liftoff, under the
# existing miniconda3_crispresso base (mamba confirmed available there).
# Run directly on the login node (not via sbatch).
set -euo pipefail

CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh

echo "Creating liftoff_env..."
mamba create -y -n liftoff_env -c bioconda -c conda-forge --channel-priority flexible liftoff

conda activate liftoff_env
echo "Verifying installation:"
liftoff --version

echo "Done. Activate with: conda activate liftoff_env (after sourcing ${CONDA_BASE}/etc/profile.d/conda.sh)"
