#!/bin/bash
#SBATCH --job-name=fastqc_all
#SBATCH --array=1-90%16         # 90 archivos, máx 16 simultáneos
#SBATCH --cpus-per-task=2        # FastQC usa pocos recursos
#SBATCH --mem=4G
#SBATCH --time=10:00:00
#SBATCH --output=logs/fastqc_%A_%a.out
#SBATCH --error=logs/fastqc_%A_%a.err
#SBATCH --partition=short
#SBATCH --mail-user=diegomartinez3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
TRIMMOMATIC_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data/trimmed_trimmomatic
FASTP_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data/trimmed_fastp

FASTQC_TRIMMOMATIC=${TRIMMOMATIC_DIR}/fastqc_results
FASTQC_FASTP=${FASTP_DIR}/fastqc_results

# ── Cargar módulo ─────────────────────────────────────────────────────────────
module load fastqc    # ajusta al nombre exacto

# ── Construir lista de todos los archivos con su output dir ───────────────────
# Genera una lista combinada: "archivo\tdirectorio_output"
FILE_LIST=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data/codes/filtering/fastqc_trimmed_filelist.txt

if [ ! -f "$FILE_LIST" ]; then

    # Trimmomatic files
    for f in "${TRIMMOMATIC_DIR}"/*.fastq.gz; do
        echo -e "$f\t$FASTQC_TRIMMOMATIC"
    done >> "$FILE_LIST"

    # fastp files
    for f in "${FASTP_DIR}"/*.fastq.gz; do
        echo -e "$f\t$FASTQC_FASTP"
    done >> "$FILE_LIST"
fi

# Verificar cantidad de archivos
TOTAL=$(wc -l < "$FILE_LIST")
echo "Total files in list: $TOTAL"

# ── Obtener archivo para esta tarea ───────────────────────────────────────────
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST" | cut -f1)
OUTDIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST" | cut -f2)

echo "Analyzing: $FILE"
echo "Output to: $OUTDIR"
echo "Start time: $(date)"

# ── Correr FastQC ─────────────────────────────────────────────────────────────
fastqc \
  "$FILE" \
  --outdir "$OUTDIR" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --noextract

echo "End time: $(date)"
echo "Done: $(basename $FILE)"