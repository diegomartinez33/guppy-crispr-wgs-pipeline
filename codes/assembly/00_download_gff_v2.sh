#!/bin/bash
# One-time setup: download the RefSeq GFF3 annotation for the new reference
# genome, GCF_904066995.2 (P_reticulata-male-v2, University of Exeter,
# PacBio+Hi-C) - the current NCBI "reference genome" for P. reticulata,
# replacing the now-suppressed GCF_000633615.1. Run directly on the login
# node (not via sbatch) - small download, compute-node internet access on
# hypatia is unverified.
set -euo pipefail

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
REF_DIR=${PROJECT_DIR}/reference
BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/904/066/995/GCF_904066995.2_P_reticulata-male-v2"
GFF_NAME="GCF_904066995.2_P_reticulata-male-v2_genomic.gff.gz"
URL="${BASE_URL}/${GFF_NAME}"
MD5_URL="${BASE_URL}/md5checksums.txt"
GZ="${REF_DIR}/${GFF_NAME}"
OUT="${REF_DIR}/GCF_904066995.2_annotation.gff"

echo "Downloading new reference (GCF_904066995.2, male, PacBio+Hi-C) GFF annotation..."
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
mv "${REF_DIR}/GCF_904066995.2_P_reticulata-male-v2_genomic.gff" "$OUT"

echo "Sanity check: bdnf gene feature must be present"
BDNF_COUNT=$(grep -c "ID=gene-bdnf" "$OUT")
echo "bdnf gene features found: ${BDNF_COUNT}"
if [ "$BDNF_COUNT" -lt 1 ]; then
    echo "ERROR: bdnf gene not found in downloaded GFF"
    exit 1
fi

echo "Done. Annotation saved to: ${OUT}"
