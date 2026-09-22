#!/bin/bash

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.." >/dev/null 2>&1 && pwd)}}"

# Limpiar si existe una lista previa
rm -f fastqc_trimmed_filelist.txt

TRIMMOMATIC_DIR=${PROJECT_DIR}/trimmed_trimmomatic
FASTP_DIR=${PROJECT_DIR}/trimmed_fastp

FASTQC_TRIMMOMATIC=${TRIMMOMATIC_DIR}/fastqc_results
FASTQC_FASTP=${FASTP_DIR}/fastqc_results

for f in "${TRIMMOMATIC_DIR}"/*.fastq.gz; do
    echo -e "$f\t$FASTQC_TRIMMOMATIC"
done >> fastqc_trimmed_filelist.txt

for f in "${FASTP_DIR}"/*.fastq.gz; do
    echo -e "$f\t$FASTQC_FASTP"
done >> fastqc_trimmed_filelist.txt

# Verificar total exacto
wc -l fastqc_trimmed_filelist.txt    # ajusta --array=1-N con este número

# Crear directorios de salida si no existen
mkdir -p "$FASTQC_TRIMMOMATIC" "$FASTQC_FASTP" logs/