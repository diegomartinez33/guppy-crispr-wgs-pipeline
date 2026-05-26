#!/bin/bash

source /hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso/etc/profile.d/conda.sh
eval "$(mamba shell hook --shell bash)"
mamba activate crispresso2_env

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
#REGION="NC_024333.1:15922039-15922058"
REGION="NC_024333.1:15921941-15922141"
#REGION="NC_024333.1:15922000-15922100"
#REGION="NC_024333.1"
SGRNA="TGAGAGACGCCCCGGGCATG"

SAMPLE=$(sed -n "1p" "$SAMPLE_LIST")

echo "Sample: $SAMPLE"

# Leer el amplicon RC directamente del archivo
AMPLICON_FWD=$(grep -v "^>" \
  ${PROJECT_DIR}/reference/amplicon_bdnf_200bp.fa | tr -d '\n')
AMPLICON_RC=$(grep -v "^>" \
  ${PROJECT_DIR}/reference/amplicon_bdnf_200bp_rc.fa | tr -d '\n')

echo "Amplicon primeros 50bp: ${AMPLICON_FWD:0:50}"
echo "Amplicon longitud: ${#AMPLICON_FWD}bp"
echo "sgRNA en amplicon: $(echo $AMPLICON_RC | grep -o $SGRNA || echo 'NO ENCONTRADO')"

# Usar un BAM de prueba
BAM="${PROJECT_DIR}/gatk/trimmomatic/markdup/Control_MNP_I_S54_L002.markdup.bam"
SAMTOOLS=$(which samtools)

# Filtrar región
$SAMTOOLS view -b -F 0x400 -h "$BAM" "NC_024333.1:15921941-15922141" -o /tmp/test_200.bam
$SAMTOOLS sort -o /tmp/test_200.sorted.bam /tmp/test_200.bam
$SAMTOOLS index /tmp/test_200.sorted.bam

# Ver exactamente qué reads está recibiendo CRISPResso2
$SAMTOOLS view /tmp/test_200.sorted.bam | head -3 | awk '{print $4, $10}' 

$SAMTOOLS view /tmp/test_200.sorted.bam > ./test_reads.sam

cp /tmp/test_200.sorted.bam ./test_200.sorted.bam
cp /tmp/test_200.sorted.bam.bai ./test_200.sorted.bam.bai

# Ver si el amplicon RC coincide con los reads
READ=$(samtools view /tmp/test_200.sorted.bam | head -1 | awk '{print $10}')
echo "Read:    ${READ:0:40}"
echo "Amplicon:${AMPLICON_FWD:0:40}"

# Buscar el read dentro del amplicon
echo "$AMPLICON_FWD" | grep -o "${READ:0:20}" \
  && echo "✅ Read encontrado en amplicon" \
  || echo "❌ Read NO encontrado en amplicon"

echo "Reads en región: $($SAMTOOLS view -c /tmp/test_200.sorted.bam)"


# -r1 ${PROJECT_DIR}/raw_fastq/${SAMPLE}_R1_001.fastq.gz \
# -r2 ${PROJECT_DIR}/raw_fastq/${SAMPLE}_R2_001.fastq.gz \
# Correr CRISPResso2
CRISPResso \
  --bam_input /tmp/test_200.sorted.bam \
  --guide_seq "$SGRNA" \
  --bam_chr_loc "$REGION" \
  --amplicon_seq "$AMPLICON_FWD,$AMPLICON_RC" \
  --amplicon_name "bdnf_fwd,bdnf_rc" \
  --output_folder /tmp/crispresso_200bp_test \
  --name "test_200bp" \
  --n_processes 4 \
  --expand_ambiguous_alignments \
  --place_report_in_output_folder \
  --write_cleaned_report

# Ver resultado
grep "N_TOT_READS\|N_COMPUTED_ALN\|No alignments" \
  /tmp/crispresso_200bp_test/CRISPResso_on_test/CRISPResso_RUNNING_LOG.txt 2>/dev/null \
  || echo "Ver log manualmente"
