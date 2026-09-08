#!/bin/bash
# One-time setup: create a dedicated conda environment for TGS-GapCloser,
# under the existing miniconda3_crispresso base (mamba confirmed available
# there - same pattern as liftoff_env, nextpolish_env). Run directly on the
# login node (not via sbatch).
set -euo pipefail

CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh

echo "Creating tgsgapcloser_env..."
mamba create -y -n tgsgapcloser_env -c bioconda -c conda-forge --channel-priority flexible tgsgapcloser

conda activate tgsgapcloser_env
echo "Verifying installation:"
tgsgapcloser --version 2>&1 || tgsgapcloser -v 2>&1 || echo "check binary name/help manually"
which tgsgapcloser

echo "Done. Activate with: conda activate tgsgapcloser_env (after sourcing ${CONDA_BASE}/etc/profile.d/conda.sh)"
