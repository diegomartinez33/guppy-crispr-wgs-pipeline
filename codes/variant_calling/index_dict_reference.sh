# Cargar módulos
module load samtools/1.16.1
module load gatk4/4.4.0.0

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF=${PROJECT_DIR}/reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna

# 1. Crear .fai (samtools — segundos)
samtools faidx "$REF"

# 2. Crear .dict (GATK — 1-2 minutos)
gatk CreateSequenceDictionary -R "$REF"

# Verificar que se crearon correctamente
ls -lh "${REF}.fai"
ls -lh "${REF%.fna}.dict"