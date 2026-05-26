#!/bin/bash
module load samtools/1.16.1

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna
REGION="NC_024333.1:15922011-15922071"

# Extraer secuencia del amplicon (±150bp alrededor del sitio de corte)
# Ajusta CHROMOSOME al nombre exacto que aparezca
samtools faidx "$REF" \
  $REGION \
  > ${PROJECT_DIR}/reference/amplicon_bdnf_60bp.fa

cat ${PROJECT_DIR}/reference/amplicon_bdnf_60bp.fa

# Extraer complemento reverso (hebra negativa)
samtools faidx "$REF" \
  -i \
  --mark-strand sign \
  $REGION \
  > ${PROJECT_DIR}/reference/amplicon_bdnf_60bp_rc.fa

cat ${PROJECT_DIR}/reference/amplicon_bdnf_60bp_rc.fa