#!/bin/bash

mkdir -p logs/

# Esperar a que el mapeo termine primero
# Reemplaza XXXXX con los job IDs de bwa_fastp y bwa_trimmomatic

# ── Pipeline fastp ────────────────────────────────────────────────────────────
#MARKDUP_F=$(sbatch --parsable mark_duplicates.sh)
#echo "Mark duplicates job ID: $MARKDUP_F"

# --dependency=afterok:${MARKDUP_F} \
# HAPLO_F=$(sbatch --parsable \
#   haplotype_caller.sh)
# echo "HaplotypeCaller job ID: $HAPLO_F"

#  --dependency=afterok:${HAPLO_F} \
GENOMICSDB_F=$(sbatch --parsable \
  genomics_db_import.sh)
echo "GenomicsDBImport job ID: $GENOMICSDB_F"

GENOTYPE_F=$(sbatch --parsable \
  --dependency=afterok:${GENOMICSDB_F} \
  genotype_gvcf.sh)
echo "GenotypeGVCFs job ID: $GENOTYPE_F"

VARFILT_F=$(sbatch --parsable \
  --dependency=afterok:${GENOTYPE_F} \
  variant_filtration.sh)
echo "VariantFiltration job ID: $VARFILT_F"

OFFT_JOB=$(sbatch --parsable \
  --dependency=afterok:${GENOTYPE_F} \
  select_offtargets.sh)
echo "Off-target variants: $OFFT_JOB"

#echo "Pipeline fastp: $MARKDUP_F → $HAPLO_F → $GENOMICSDB_F → $GENOTYPE_F → $VARFILT_F"
echo "Pipeline trimmomatic: $GENOMICSDB_F → $GENOTYPE_F → $VARFILT_F → $OFFT_JOB"