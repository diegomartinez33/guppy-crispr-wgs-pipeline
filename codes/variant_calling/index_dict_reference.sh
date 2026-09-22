# Cargar módulos
module load samtools/1.16.1
module load gatk4/4.4.0.0

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.." >/dev/null 2>&1 && pwd)}}"
source "${PROJECT_DIR}/codes/genome_versions.sh"

# 1. Crear .fai (samtools — segundos)
samtools faidx "$REF"

# 2. Crear .dict (GATK — 1-2 minutos)
gatk CreateSequenceDictionary -R "$REF"

# Verificar que se crearon correctamente
ls -lh "${REF}.fai"
ls -lh "${REF%.fna}.dict"