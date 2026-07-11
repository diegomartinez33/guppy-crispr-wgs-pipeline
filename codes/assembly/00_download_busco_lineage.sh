#!/bin/bash
# One-time setup: pre-stage the BUSCO actinopterygii_odb12.2 lineage dataset
# so Phase 3 (qc_quast_busco.sh) can run fully --offline on a compute node.
# actinopterygii_odb10 no longer exists in the BUSCO 5.7.1 catalog - odb12.2
# is the current replacement for this clade.
# Run directly on the login node (not via sbatch).
#
# NOTE: `busco --download actinopterygii_odb12.2` is broken in this module
# (busco/5.7.1) - it mis-parses the ".2" suffix and requests
# ".../actinopterygii_odb12.2026-05-13.tar.gz" instead of the real
# ".../actinopterygii_odb12.2.2026-05-13.tar.gz", which 404s every time
# (confirmed reproducibly). Downloading and staging the tarball manually
# works around this - the file itself is fine, only busco's own URL
# construction for this specific lineage name is buggy.
set -euo pipefail

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
DOWNLOAD_PATH=${PROJECT_DIR}/assembly/qc_results/busco_downloads
LINEAGES_DIR=${DOWNLOAD_PATH}/lineages
URL="https://busco-data.ezlab.org/v5/data/lineages/actinopterygii_odb12.2.2026-05-13.tar.gz"
TARBALL=${DOWNLOAD_PATH}/actinopterygii_odb12.2.tar.gz

mkdir -p "$LINEAGES_DIR"

echo "Downloading BUSCO lineage dataset actinopterygii_odb12.2 (manual staging, see script comment)..."
wget -O "$TARBALL" "$URL"

echo "Extracting to ${LINEAGES_DIR}/..."
tar -xzf "$TARBALL" -C "$LINEAGES_DIR"

if [ ! -d "${LINEAGES_DIR}/actinopterygii_odb12.2" ]; then
    echo "ERROR: expected extracted folder ${LINEAGES_DIR}/actinopterygii_odb12.2 not found"
    exit 1
fi

rm -f "$TARBALL"
echo "Done. Lineage dataset staged at: ${LINEAGES_DIR}/actinopterygii_odb12.2"
