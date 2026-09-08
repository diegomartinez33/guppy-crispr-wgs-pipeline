#!/bin/bash
#SBATCH --job-name=crispresso_ontarget
#SBATCH --array=1-15%8
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/crispresso_%A_%a.out
#SBATCH --error=logs/crispresso_%A_%a.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
source "${PROJECT_DIR}/codes/genome_versions.sh"
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
INPUT_DIR=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/markdup
OUTPUT_DIR=${PROJECT_DIR}/crispresso${OUT_SUFFIX}/ontarget/trimmomatic
# See extract_amplicon_sgRNA.sh for how the v2 coordinate was relocated
# (minimap2 --cs alignment of the v1 window against the new bdnf gene span
# - 100% identity, NM:i:0).
if [ "$REF_VERSION" = "v1" ]; then
    REGION="NC_024333.1:15922000-15922100"
    AMPLICON_SUFFIX=""
elif [ "$REF_VERSION" = "v2" ]; then
    REGION="NC_088832.1:15849655-15849755"
    AMPLICON_SUFFIX="_v2"
fi

mkdir -p "$OUTPUT_DIR" logs/

# ── Activar ambiente CRISPResso2 ──────────────────────────────────────────────
source /hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso/etc/profile.d/conda.sh
eval "$(mamba shell hook --shell bash)"
mamba activate crispresso2_env

# ── Samtools del ambiente conda (NO module load) ──────────────────────────────
SAMTOOLS=$(which samtools)
echo "Using samtools: $SAMTOOLS"

# ── Secuencias del experimento ────────────────────────────────────────────────
SGRNA="TGAGAGACGCCCCGGGCATG"
# Leer el amplicon RC directamente del archivo
AMPLICON_FWD=$(grep -v "^>" \
  "${PROJECT_DIR}/reference/amplicon_bdnf_100bp${AMPLICON_SUFFIX}.fa" | tr -d '\n')
AMPLICON_RC=$(grep -v "^>" \
  "${PROJECT_DIR}/reference/amplicon_bdnf_100bp${AMPLICON_SUFFIX}_rc.fa" | tr -d '\n')
# ── Muestra para esta tarea ───────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
BAM="${INPUT_DIR}/${SAMPLE}.markdup.bam"

echo "Sample:        $SAMPLE"
echo "Input BAM:     $BAM"
echo "Region:        $REGION"
echo "Start time:    $(date)"

# ── Verificar que el input existe ─────────────────────────────────────────────
if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM not found: $BAM"
    exit 1
fi

# ── Filtrar duplicados y región de interés ────────────────────────────────────
FILTERED_BAM="${OUTPUT_DIR}/${SAMPLE}_filtered.bam"
SORTED_BAM="${OUTPUT_DIR}/${SAMPLE}_filtered.sorted.bam"

# ✅ Todo en una línea — sin \ ni comentarios intermedios
# ✅ La región va como argumento posicional al final, sin -r
$SAMTOOLS view -b -F 0x400 -h "$BAM" "$REGION" -o "$FILTERED_BAM"

# Verificar que el BAM filtrado no está vacío
if [ ! -s "$FILTERED_BAM" ]; then
    echo "ERROR: Filtered BAM is empty — no reads in region $REGION"
    exit 1
fi

echo "Reads before sort: $($SAMTOOLS view -c $FILTERED_BAM)"

# Ordenar e indexar
$SAMTOOLS sort -o "$SORTED_BAM" "$FILTERED_BAM"
$SAMTOOLS index "$SORTED_BAM"

echo "Reads after sort:  $($SAMTOOLS view -c $SORTED_BAM)"

# ── CRISPResso2 ───────────────────────────────────────────────────────────────
CRISPResso \
  --bam_input "$SORTED_BAM" \
  --amplicon_seq "$AMPLICON_FWD,$AMPLICON_RC" \
  --amplicon_name "bdnf_fwd,bdnf_rc" \
  --guide_seq "$SGRNA" \
  --bam_chr_loc "$REGION" \
  --output_folder "${OUTPUT_DIR}/${SAMPLE}" \
  --name "${SAMPLE}" \
  --n_processes "${SLURM_CPUS_PER_TASK}" \
  --min_frequency_alleles_around_cut_to_plot 0.05 \
  --expand_ambiguous_alignments \
  --place_report_in_output_folder \
  --write_cleaned_report

# ── Limpiar BAMs temporales ───────────────────────────────────────────────────
rm -f "$FILTERED_BAM" "$SORTED_BAM" "${SORTED_BAM}.bai"

echo "End time: $(date)"
echo "Done: $SAMPLE"