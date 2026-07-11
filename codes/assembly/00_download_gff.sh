#!/bin/bash
# One-time setup: download the Trinidad/Guanapo (GCF_000633615.1) RefSeq GFF3
# annotation from NCBI. Run directly on the login node (not via sbatch) -
# it's a 14.8MB download that completes in seconds, and compute-node
# internet access on hypatia is unverified.
set -euo pipefail

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF_DIR=${PROJECT_DIR}/reference
BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/633/615/GCF_000633615.1_Guppy_female_1.0_MT"
GFF_NAME="GCF_000633615.1_Guppy_female_1.0_MT_genomic.gff.gz"
URL="${BASE_URL}/${GFF_NAME}"
MD5_URL="${BASE_URL}/md5checksums.txt"
GZ="${REF_DIR}/${GFF_NAME}"
OUT="${REF_DIR}/GCF_000633615.1_annotation.gff"

echo "Downloading Trinidad reference GFF annotation..."
wget -O "$GZ" "$URL"

echo "Verifying checksum against NCBI md5checksums.txt..."
EXPECTED_MD5=$(curl -s "$MD5_URL" | grep "\./${GFF_NAME}$" | awk '{print $1}')
if [ -z "$EXPECTED_MD5" ]; then
    echo "ERROR: could not find ${GFF_NAME} entry in md5checksums.txt"
    exit 1
fi
ACTUAL_MD5=$(md5sum "$GZ" | awk '{print $1}')
if [ "$EXPECTED_MD5" != "$ACTUAL_MD5" ]; then
    echo "ERROR: md5 mismatch for ${GFF_NAME}"
    echo "  expected: $EXPECTED_MD5"
    echo "  actual:   $ACTUAL_MD5"
    exit 1
fi
echo "Checksum OK: $ACTUAL_MD5"

gunzip -k "$GZ"
mv "${REF_DIR}/GCF_000633615.1_Guppy_female_1.0_MT_genomic.gff" "$OUT"

echo "Sanity check: bdnf gene feature must be present"
BDNF_COUNT=$(grep -c "ID=gene-bdnf" "$OUT")
echo "bdnf gene features found: ${BDNF_COUNT}"
if [ "$BDNF_COUNT" -lt 1 ]; then
    echo "ERROR: bdnf gene not found in downloaded GFF"
    exit 1
fi

echo "Done. Annotation saved to: ${OUT}"
