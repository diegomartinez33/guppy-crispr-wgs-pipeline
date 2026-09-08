#!/bin/bash
# One-time setup: install the classic (boulder-IO) Primer3 engine that
# EMBOSS's `eprimer3` wrapper needs but does not ship itself.
#
# `module load emboss/6.6.0` only provides the wrapper (`eprimer3`), which
# shells out to an external binary literally named `primer3_core` using the
# old boulder-IO protocol (Primer3 <=1.x) - confirmed by running eprimer3
# on this cluster: "Died: eprimer3 uses external program 'primer3_core'
# which is not in the PATH or defined as EMBOSS_PRIMER3_CORE". The
# default/current primer3 build (2.6.1, e.g. from the base conda channel)
# uses a different, incompatible tag-value format and will NOT work here -
# must be exactly the legacy 1.1.4 build.
#
# Run this once, directly on the login node (internet access confirmed
# there, same as the NCBI GFF downloads elsewhere in this project). Not
# needed again unless the primer3_env conda environment is removed.

CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh

mamba create -y -n primer3_env -c bioconda -c conda-forge primer3=1.1.4

echo ""
echo "Done. Every script that calls eprimer3 must export:"
echo "  export EMBOSS_PRIMER3_CORE=${CONDA_BASE}/envs/primer3_env/bin/primer3_core"
echo "(already done by run_offtarget_primer_design.sh)"
