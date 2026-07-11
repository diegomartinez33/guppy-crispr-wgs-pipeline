#!/bin/bash
# Full GATK variant-calling pipeline — trimmomatic track
# Chain: HaplotypeCaller (15 samples, medium partition) → GenomicsDB → GenotypeGVCFs → VariantFiltration

mkdir -p logs/

cd "$(dirname "$0")"

HAPLO_JID=$(sbatch --parsable haplotype_caller.sh)
echo "HaplotypeCaller job ID: ${HAPLO_JID}  (15 samples, medium partition)"

GENOMICSDB_JID=$(sbatch --parsable \
    --dependency=afterok:${HAPLO_JID} \
    genomics_db_import.sh)
echo "GenomicsDBImport job ID: ${GENOMICSDB_JID}"

GENOTYPE_JID=$(sbatch --parsable \
    --dependency=afterok:${GENOMICSDB_JID} \
    genotype_gvcf.sh)
echo "GenotypeGVCFs job ID: ${GENOTYPE_JID}"

VARFILT_JID=$(sbatch --parsable \
    --dependency=afterok:${GENOTYPE_JID} \
    variant_filtration.sh)
echo "VariantFiltration job ID: ${VARFILT_JID}"

OFFT_JID=$(sbatch --parsable \
    --dependency=afterok:${GENOTYPE_JID} \
    select_offtargets.sh)
echo "Off-target variants job ID: ${OFFT_JID}"

echo ""
echo "Pipeline: ${GENOMICSDB_JID} → ${GENOTYPE_JID} → ${VARFILT_JID} + ${OFFT_JID}"
