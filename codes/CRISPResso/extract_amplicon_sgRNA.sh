#!/bin/bash
module load samtools/1.16.1

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
source "${PROJECT_DIR}/codes/genome_versions.sh"

# REGION is genome-specific and can't be derived automatically each run -
# v2's value was found by aligning the v1 61bp amplicon window (minimap2
# --cs, sr preset) against the new bdnf gene span (NC_088832.1:15848607-
# 15863146, from GCF_904066995.2_annotation.gff): 100% identity (NM:i:0),
# confirming this sgRNA/cut-site region is fully conserved between the two
# assemblies. See CLAUDE.md, "Migration to GCF_904066995.2 (v2)".
if [ "$REF_VERSION" = "v1" ]; then
    REGION="NC_024333.1:15922011-15922071"
    OUT_SUFFIX_FILE=""
elif [ "$REF_VERSION" = "v2" ]; then
    REGION="NC_088832.1:15849666-15849726"
    OUT_SUFFIX_FILE="_v2"
fi

# Extraer secuencia del amplicon (±150bp alrededor del sitio de corte)
samtools faidx "$REF" \
  $REGION \
  > "${PROJECT_DIR}/reference/amplicon_bdnf_60bp${OUT_SUFFIX_FILE}.fa"

cat "${PROJECT_DIR}/reference/amplicon_bdnf_60bp${OUT_SUFFIX_FILE}.fa"

# Extraer complemento reverso (hebra negativa)
samtools faidx "$REF" \
  -i \
  --mark-strand sign \
  $REGION \
  > "${PROJECT_DIR}/reference/amplicon_bdnf_60bp${OUT_SUFFIX_FILE}_rc.fa"

cat "${PROJECT_DIR}/reference/amplicon_bdnf_60bp${OUT_SUFFIX_FILE}_rc.fa"