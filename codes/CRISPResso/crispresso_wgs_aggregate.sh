#!/bin/bash
#SBATCH --job-name=crispresso_wgs_aggregate
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/crispresso_wgs_aggregate.out
#SBATCH --error=logs/crispresso_wgs_aggregate.err
#SBATCH --partition=short

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
WGS_DIR=${PROJECT_DIR}/crispresso/wgs/trimmomatic
OUTPUT_DIR=${PROJECT_DIR}/crispresso/wgs/trimmomatic/aggregate

mkdir -p "$OUTPUT_DIR" logs/

CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate crispresso2_env
export PATH="${CONDA_BASE}/envs/crispresso2_env/bin:$PATH"

echo "Start time: $(date)"

# CRISPRessoAggregate has no output-directory flag - it always creates its
# folder in the current working directory, so we must cd into OUTPUT_DIR
# first (first attempt, May 2026, missed this and wrote output into
# codes/CRISPResso/ instead - see "CRISPRessoAggregate" Known Issue).
cd "$OUTPUT_DIR"

# ── Aggregate por grupo ───────────────────────────────────────────────────────
for GROUP in Control Only_MNP Plasmid_Ko RNP_Cas; do
    echo "=== Aggregating: $GROUP ==="

    case $GROUP in
        Control)    SAMPLES="Control_MNP_I_S54_L002 Control_MNP_II_S55_L002 Control_MNP_III_S56_L002" ;;
        Only_MNP)   SAMPLES="Only_MNP_C1_S57_L002 Only_MNP_C2_S58_L002 Only_MNP_C3_S59_L002 Only_MNP_C4_S60_L002" ;;
        Plasmid_Ko) SAMPLES="Plasmid_Ko_P1_S61_L002 Plasmid_Ko_P2_S62_L002 Plasmid_Ko_P3_S63_L002 Plasmid_Ko_P4_S64_L002" ;;
        RNP_Cas)    SAMPLES="RNP_Cas1_S65_L002 RNP_Cas2_S66_L002 RNP_Cas3_S67_L002 RNP_Cas4_S68_L002" ;;
    esac

    # Construir prefijos: -p debe apuntar DENTRO de CRISPRessoWGS_on_<sample>/
    # a "CRISPResso_on_" (prefijo de las 8 corridas por sitio), no a la
    # carpeta de la muestra en sí - CRISPRessoAggregate hace glob-match sobre
    # el prefijo dado, no escanea subdirectorios recursivamente.
    PREFIX_ARGS=""
    for SAMPLE in $SAMPLES; do
        PREFIX_ARGS="$PREFIX_ARGS -p ${WGS_DIR}/${SAMPLE}/CRISPRessoWGS_on_${SAMPLE}/CRISPResso_on_"
    done

    # --suppress_report: esta versión de CRISPResso2 (2.3.1) tiene un bug en
    # make_multi_report() (faltan argumentos posicionales) que hace fallar
    # SOLO la generación del reporte HTML final - las tablas/plots de datos
    # (lo que realmente importa) se generan correctamente antes de ese punto.
    CRISPRessoAggregate \
      $PREFIX_ARGS \
      --name "${GROUP}_wgs_aggregate" \
      --min_reads_for_inclusion 10 \
      --place_report_in_output_folder \
      --suppress_report \
      --n_processes 4

    echo "✅ $GROUP aggregated"
done

echo "End time: $(date)"