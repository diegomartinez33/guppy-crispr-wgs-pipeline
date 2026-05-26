#!/bin/bash
#SBATCH --job-name=parse_coverage
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/parse_coverage.out
#SBATCH --error=logs/parse_coverage.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
COVERAGE_DIR=${PROJECT_DIR}/coverage/bdnf_site_v2
OUTPUT_DIR=${PROJECT_DIR}/coverage/csv
SAMPLE_LIST=${PROJECT_DIR}/samples.txt

mkdir -p "$OUTPUT_DIR"

# ── Coordenadas del sitio CRISPR ──────────────────────────────────────────────
SGRNA_START=15922039
SGRNA_END=15922058
CUT_SITE=$(( (SGRNA_START + SGRNA_END) / 2 ))   # centro del sgRNA ~15922048

echo "Start time: $(date)"
echo "Cut site position: $CUT_SITE"

# ─────────────────────────────────────────────────────────────────────────────
# CSV 1: Resumen de cobertura por muestra (desde *_coverage.txt)
# ─────────────────────────────────────────────────────────────────────────────
SUMMARY_CSV="${OUTPUT_DIR}/coverage_summary_all_samples.csv"

echo "Sample,Group,Chromosome,Region_Start,Region_End,\
Total_Reads,Covered_Bases,Coverage_Pct,Mean_Depth,\
Mean_BaseQ,Mean_MapQ" > "$SUMMARY_CSV"

while IFS= read -r SAMPLE; do
    COV_FILE="${COVERAGE_DIR}/${SAMPLE}_coverage.txt"

    if [ ! -f "$COV_FILE" ]; then
        echo "WARNING: $COV_FILE not found, skipping"
        continue
    fi

    # Extraer grupo del nombre de muestra
    if [[ $SAMPLE == Control* ]];   then GROUP="Control"
    elif [[ $SAMPLE == Only_MNP* ]]; then GROUP="Only_MNP"
    elif [[ $SAMPLE == Plasmid* ]];  then GROUP="Plasmid_Ko"
    elif [[ $SAMPLE == RNP* ]];      then GROUP="RNP_Cas"
    else                                  GROUP="Unknown"
    fi

    # Parsear el archivo coverage (saltar header con NR>1)
    awk -v sample="$SAMPLE" -v group="$GROUP" 'NR>1 {
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
            sample, group,
            $1, $2, $3,     # chr, start, end
            $4, $5, $6,     # numreads, covbases, coverage%
            $7, $8, $9      # meandepth, meanbaseq, meanmapq
    }' "$COV_FILE" >> "$SUMMARY_CSV"

done < "$SAMPLE_LIST"

echo "✅ Coverage summary CSV: $SUMMARY_CSV"

# ─────────────────────────────────────────────────────────────────────────────
# CSV 2: Profundidad por posición — todas las muestras combinadas
# ─────────────────────────────────────────────────────────────────────────────
DEPTH_CSV="${OUTPUT_DIR}/depth_per_position_all_samples.csv"

# Header dinámico con nombre de cada muestra
HEADER="Position"
while IFS= read -r SAMPLE; do
    HEADER="${HEADER},${SAMPLE}"
done < "$SAMPLE_LIST"
echo "$HEADER" > "$DEPTH_CSV"

# Crear archivos temporales de depth por muestra (solo columnas pos + depth)
TMPDIR_LOCAL="/tmp/${USER}_coverage_$$"
mkdir -p "$TMPDIR_LOCAL"

while IFS= read -r SAMPLE; do
    DEPTH_FILE="${COVERAGE_DIR}/${SAMPLE}_depth.txt"
    if [ -f "$DEPTH_FILE" ]; then
        awk '{print $2"\t"$3}' "$DEPTH_FILE" \
          > "${TMPDIR_LOCAL}/${SAMPLE}_pos_depth.txt"
    fi
done < "$SAMPLE_LIST"

# Combinar todos los archivos por posición usando paste
FIRST_SAMPLE=$(head -1 "$SAMPLE_LIST")
POSITIONS_FILE="${TMPDIR_LOCAL}/${FIRST_SAMPLE}_pos_depth.txt"

# Extraer posiciones del primer archivo
awk '{print $1}' "$POSITIONS_FILE" > "${TMPDIR_LOCAL}/positions.txt"

# Construir matriz: posición + depth de cada muestra
PASTE_FILES="${TMPDIR_LOCAL}/positions.txt"
while IFS= read -r SAMPLE; do
    PASTE_FILES="$PASTE_FILES ${TMPDIR_LOCAL}/${SAMPLE}_pos_depth.txt"
done < "$SAMPLE_LIST"

paste $PASTE_FILES \
  | awk 'BEGIN{OFS=","} {
        # posición en col 1, luego depth en col 3,5,7... (skip pos repetidas)
        printf $1
        for(i=2; i<=NF; i+=2) printf "," $i
        printf "\n"
    }' >> "$DEPTH_CSV"

rm -rf "$TMPDIR_LOCAL"

echo "✅ Depth per position CSV: $DEPTH_CSV"

# ─────────────────────────────────────────────────────────────────────────────
# CSV 3: Estadísticas de depth por muestra y zona relativa al sitio de corte
# ─────────────────────────────────────────────────────────────────────────────
ZONE_CSV="${OUTPUT_DIR}/depth_by_zone.csv"

echo "Sample,Group,Zone,Zone_Start,Zone_End,\
Mean_Depth,Min_Depth,Max_Depth,Positions" > "$ZONE_CSV"

while IFS= read -r SAMPLE; do
    DEPTH_FILE="${COVERAGE_DIR}/${SAMPLE}_depth.txt"

    if [ ! -f "$DEPTH_FILE" ]; then continue; fi

    if [[ $SAMPLE == Control* ]];    then GROUP="Control"
    elif [[ $SAMPLE == Only_MNP* ]]; then GROUP="Only_MNP"
    elif [[ $SAMPLE == Plasmid* ]];  then GROUP="Plasmid_Ko"
    elif [[ $SAMPLE == RNP* ]];      then GROUP="RNP_Cas"
    else                                  GROUP="Unknown"
    fi

    # Código awk guardado en variable separada
    AWK_CODE='
    {
        pos   = $2
        depth = $3

        if (pos >= sg_start && pos <= sg_end)
            zone = "sgRNA_site"
        else if (pos < sg_start && pos >= sg_start - 100)
            zone = "upstream_100bp"
        else if (pos > sg_end && pos <= sg_end + 100)
            zone = "downstream_100bp"
        else if (pos < sg_start - 100)
            zone = "upstream_500bp"
        else
            zone = "downstream_500bp"

        sum[zone]   += depth
        count[zone] += 1

        if (depth > max[zone] || max[zone] == "") max[zone] = depth
        if (depth < min[zone] || min[zone] == "") min[zone] = depth

        start[zone] = (start[zone] == "" || pos < start[zone]) ? pos : start[zone]
        end[zone]   = (end[zone]   == "" || pos > end[zone])   ? pos : end[zone]
    }
    END {
        n        = split("upstream_500bp,upstream_100bp,sgRNA_site,downstream_100bp,downstream_500bp", zones, ",")
        for (i = 1; i <= n; i++) {
            z = zones[i]
            if (count[z] > 0) {
                printf "%s,%s,%s,%s,%s,%.2f,%s,%s,%s\n",
                    sample, group, z,
                    start[z], end[z],
                    sum[z] / count[z],
                    min[z], max[z], count[z]
            }
        }
    }'

    # Llamada a awk limpia — variables bash pasadas con -v
    awk \
        -v sample="$SAMPLE" \
        -v group="$GROUP" \
        -v cut="$CUT_SITE" \
        -v sg_start="$SGRNA_START" \
        -v sg_end="$SGRNA_END" \
        "$AWK_CODE" \
        "$DEPTH_FILE" >> "$ZONE_CSV"

done < "$SAMPLE_LIST"

echo "✅ Depth by zone CSV: $ZONE_CSV"

# ─────────────────────────────────────────────────────────────────────────────
# CSV 4: Coverage por zona (num reads, cov bases, coverage%) para cada muestra
# ─────────────────────────────────────────────────────────────────────────────

# CSV de coverage por zona
ZONE_COV_CSV="${OUTPUT_DIR}/coverage_by_zone_all_samples.csv"
echo "Sample,Group,Zone,Region,NumReads,CovBases,Coverage_Pct,Mean_Depth,Mean_BaseQ,Mean_MapQ" \
  > "$ZONE_COV_CSV"

while IFS= read -r SAMPLE; do
    ZONE_FILE="${COVERAGE_DIR}/${SAMPLE}_coverage_by_zone.txt"

    if [ ! -f "$ZONE_FILE" ]; then continue; fi

    if [[ $SAMPLE == Control* ]];    then GROUP="Control"
    elif [[ $SAMPLE == Only_MNP* ]]; then GROUP="Only_MNP"
    elif [[ $SAMPLE == Plasmid* ]];  then GROUP="Plasmid_Ko"
    elif [[ $SAMPLE == RNP* ]];      then GROUP="RNP_Cas"
    fi

    awk -F'\t' -v sample="$SAMPLE" -v group="$GROUP" \
    'NR>1 {print sample","group","$1","$2","$3","$4","$5","$6","$7","$8}' \
    "$ZONE_FILE" >> "$ZONE_COV_CSV"

done < "$SAMPLE_LIST"

echo "✅ Coverage by zone CSV: $ZONE_COV_CSV"

# ─────────────────────────────────────────────────────────────────────────────
# Resumen final
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "=== CSVs generados ==="
echo "1. $(basename $SUMMARY_CSV)     → resumen por muestra"
echo "2. $(basename $DEPTH_CSV) → profundidad por posición"
echo "3. $(basename $ZONE_CSV)        → stats por zona relativa al corte"
echo "4. $(basename $ZONE_COV_CSV)   → coverage por zona"
echo ""
echo "=== Preview: Coverage Summary ==="
column -t -s',' "$SUMMARY_CSV"

echo ""
echo "=== Preview: Depth by Zone ==="
column -t -s',' "$ZONE_CSV"

echo "End time: $(date)"