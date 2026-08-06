#!/bin/bash
# One-time setup: create a dedicated conda environment for NextPolish, under
# the existing miniconda3_crispresso base (mamba confirmed available there -
# same pattern as liftoff_env). Run directly on the login node (not sbatch).
set -euo pipefail

CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh

echo "Creating nextpolish_env..."
mamba create -y -n nextpolish_env -c bioconda -c conda-forge --channel-priority flexible nextpolish

conda activate nextpolish_env
echo "Verifying installation:"
nextPolish --version
which nextPolish

echo "Done. Activate with: conda activate nextpolish_env (after sourcing ${CONDA_BASE}/etc/profile.d/conda.sh)"
