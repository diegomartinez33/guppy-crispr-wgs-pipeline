#!/bin/bash
# Full GATK variant-calling pipeline — trimmomatic track
# Chain: HaplotypeCaller (15 samples, medium partition) → GenomicsDB → GenotypeGVCFs → VariantFiltration
#
# REF_VERSION=v2 ./run_variant_calling.sh runs the whole chain against the
# new reference genome (GCF_904066995.2) instead of v1 - see
# codes/genome_versions.sh and CLAUDE.md, "Migration to GCF_904066995.2 (v2)".
# All v1 paths/outputs are untouched either way (OUT_SUFFIX is empty for v1).

REF_VERSION=${REF_VERSION:-v1}

mkdir -p logs/

cd "$(dirname "$0")"

HAPLO_JID=$(sbatch --parsable --export=ALL,REF_VERSION="${REF_VERSION}" haplotype_caller.sh)
echo "HaplotypeCaller job ID: ${HAPLO_JID}  (15 samples, medium partition, REF_VERSION=${REF_VERSION})"

GENOMICSDB_JID=$(sbatch --parsable \
    --dependency=afterok:${HAPLO_JID} \
    --export=ALL,REF_VERSION="${REF_VERSION}" \
    genomics_db_import.sh)
echo "GenomicsDBImport job ID: ${GENOMICSDB_JID}"

GENOTYPE_JID=$(sbatch --parsable \
    --dependency=afterok:${GENOMICSDB_JID} \
    --export=ALL,REF_VERSION="${REF_VERSION}" \
    genotype_gvcf.sh)
echo "GenotypeGVCFs job ID: ${GENOTYPE_JID}"

VARFILT_JID=$(sbatch --parsable \
    --dependency=afterok:${GENOTYPE_JID} \
    --export=ALL,REF_VERSION="${REF_VERSION}" \
    variant_filtration.sh)
echo "VariantFiltration job ID: ${VARFILT_JID}"

echo ""
echo "Pipeline: ${HAPLO_JID} → ${GENOMICSDB_JID} → ${GENOTYPE_JID} → ${VARFILT_JID}"
echo "(select_offtargets.sh is run separately, once CRISPOR/Cas-OFFinder off-target sites for this REF_VERSION are known)"
