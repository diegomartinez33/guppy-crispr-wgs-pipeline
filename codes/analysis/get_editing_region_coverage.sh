#!/bin/bash
#SBATCH --job-name=coverage_bdnf
#SBATCH --array=1-15%8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/coverage_%A_%a.out
#SBATCH --error=logs/coverage_%A_%a.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
BAM_DIR=${PROJECT_DIR}/mapping/trimmomatic
OUTPUT_DIR=${PROJECT_DIR}/coverage/bdnf_site_v2

mkdir -p "$OUTPUT_DIR" logs/

module load samtools/1.16.1

# ── Coordenadas del sitio CRISPR ──────────────────────────────────────────────
CHROMOSOME="NC_024333.1"
SGRNA_START=15922039
SGRNA_END=15922058
WINDOW=600

# ── Verificar tamaños de zonas ────────────────────────────────────────────────
echo "=== Zone size verification ==="
echo "upstream_100bp:   $(( (SGRNA_START-1) - (SGRNA_START-100) + 1 ))bp"
echo "upstream_500bp:   $(( (SGRNA_START-101) - (SGRNA_START-WINDOW) + 1 ))bp"
echo "sgRNA_site:       $(( SGRNA_END - SGRNA_START + 1 ))bp"
echo "downstream_100bp: $(( (SGRNA_END+100) - (SGRNA_END+1) + 1 ))bp"
echo "downstream_500bp: $(( (SGRNA_END+WINDOW) - (SGRNA_END+101) + 1 ))bp"
echo ""

# ── Definir zonas ─────────────────────────────────────────────────────────────
REGION_FULL="${CHROMOSOME}:$((SGRNA_START - WINDOW))-$((SGRNA_END + WINDOW))"
REGION_SGRNA="${CHROMOSOME}:${SGRNA_START}-${SGRNA_END}"
REGION_UP100="${CHROMOSOME}:$((SGRNA_START - 100))-$((SGRNA_START - 1))"
REGION_UP500="${CHROMOSOME}:$((SGRNA_START - WINDOW))-$((SGRNA_START - 101))"
REGION_DN100="${CHROMOSOME}:$((SGRNA_END + 1))-$((SGRNA_END + 100))"
REGION_DN500="${CHROMOSOME}:$((SGRNA_END + 101))-$((SGRNA_END + WINDOW))"

# ── Muestra para esta tarea ───────────────────────────────────────────────────
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
BAM="${BAM_DIR}/${SAMPLE}.sorted.bam"

echo "Sample:      $SAMPLE"
echo "Full region: $REGION_FULL"
echo "Start time:  $(date)"
echo ""

# ── Verificar BAM ─────────────────────────────────────────────────────────────
if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM not found: $BAM"
    exit 1
fi

# ── Cobertura por posición (región completa) ──────────────────────────────────
samtools depth \
  -a \
  -r "$REGION_FULL" \
  "$BAM" \
  > "${OUTPUT_DIR}/${SAMPLE}_depth.txt"

# ── Coverage por zona — llamadas directas ─────────────────────────────────────
ZONE_COV="${OUTPUT_DIR}/${SAMPLE}_coverage_by_zone.txt"

# Header
echo -e "Zone\tRegion\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" \
  > "$ZONE_COV"

# full_region
samtools coverage -r "$REGION_FULL" "$BAM" \
  | awk -v z="full_region" -v r="$REGION_FULL" \
    'NR>1 {printf z"\t"r"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\n"}' \
  >> "$ZONE_COV"

# sgRNA_site
samtools coverage -r "$REGION_SGRNA" "$BAM" \
  | awk -v z="sgRNA_site" -v r="$REGION_SGRNA" \
    'NR>1 {printf z"\t"r"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\n"}' \
  >> "$ZONE_COV"

# upstream_100bp
samtools coverage -r "$REGION_UP100" "$BAM" \
  | awk -v z="upstream_100bp" -v r="$REGION_UP100" \
    'NR>1 {printf z"\t"r"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\n"}' \
  >> "$ZONE_COV"

# upstream_500bp
samtools coverage -r "$REGION_UP500" "$BAM" \
  | awk -v z="upstream_500bp" -v r="$REGION_UP500" \
    'NR>1 {printf z"\t"r"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\n"}' \
  >> "$ZONE_COV"

# downstream_100bp
samtools coverage -r "$REGION_DN100" "$BAM" \
  | awk -v z="downstream_100bp" -v r="$REGION_DN100" \
    'NR>1 {printf z"\t"r"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\n"}' \
  >> "$ZONE_COV"

# downstream_500bp
samtools coverage -r "$REGION_DN500" "$BAM" \
  | awk -v z="downstream_500bp" -v r="$REGION_DN500" \
    'NR>1 {printf z"\t"r"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\n"}' \
  >> "$ZONE_COV"

# ── Coverage global (región completa) ─────────────────────────────────────────
samtools coverage \
  -r "$REGION_FULL" \
  "$BAM" \
  > "${OUTPUT_DIR}/${SAMPLE}_coverage.txt"

# ── Resumen en el log ─────────────────────────────────────────────────────────
echo "=== Coverage by zone: $SAMPLE ==="
column -t "$ZONE_COV"

echo ""
echo "=== Depth stats by zone ==="
for ZONE_NAME in sgRNA_site upstream_100bp upstream_500bp downstream_100bp downstream_500bp; do
    case $ZONE_NAME in
        sgRNA_site)       REGION="$REGION_SGRNA" ;;
        upstream_100bp)   REGION="$REGION_UP100" ;;
        upstream_500bp)   REGION="$REGION_UP500" ;;
        downstream_100bp) REGION="$REGION_DN100" ;;
        downstream_500bp) REGION="$REGION_DN500" ;;
    esac

    echo "--- $ZONE_NAME ($REGION) ---"
    samtools depth -a -r "$REGION" "$BAM" \
      | awk '{
          sum+=$3; count++
          if($3>max) max=$3
          if(min=="" || $3<min) min=$3
        }
        END {
          printf "  Mean: %.2f×  Min: %s×  Max: %s×  Positions: %s\n",
            sum/count, min, max, count
        }'
done

echo ""
echo "End time: $(date)"
echo "Done: $SAMPLE"