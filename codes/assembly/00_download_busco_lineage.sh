#!/bin/bash
# One-time setup: pre-stage the BUSCO actinopterygii_odb10 lineage dataset
# so Phase 3 (busco_qc.sh) can run fully --offline on a compute node.
# Run directly on the login node (not via sbatch).
#
# NOTE: actinopterygii_odb10 no longer appears in `busco --list-datasets`
# (the interactive catalog now shows odb12.2 as current for this clade),
# but BUSCO 5.7.1's own code (BuscoConfig.py) HARDCODES a requirement for
# datasets_version == "odb10" - it fatally errors on any other version
# regardless of the --datasets_version flag (confirmed by reading the
# source: the final check at BuscoConfig.py:673 is unconditional). So for
# this BUSCO module version, odb10 is not just the old default - it's the
# only version that will actually run. The archived odb10 tarball
# (2024-01-08) is still hosted and downloadable, just delisted from the
# interactive catalog, similar to the odb12.2 URL-construction bug found
# earlier. See "BUSCO — hardcoded odb10 version check" in CLAUDE.md Known
# Issues for the full debugging history (including the abandoned odb12.2
# attempt).
set -euo pipefail

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
DOWNLOAD_PATH=${PROJECT_DIR}/assembly/qc_results/busco_downloads
LINEAGES_DIR=${DOWNLOAD_PATH}/lineages
URL="https://busco-data.ezlab.org/v5/data/lineages/actinopterygii_odb10.2024-01-08.tar.gz"
TARBALL=${DOWNLOAD_PATH}/actinopterygii_odb10.tar.gz

mkdir -p "$LINEAGES_DIR"

echo "Downloading BUSCO lineage dataset actinopterygii_odb10 (manual staging, see script comment)..."
wget -O "$TARBALL" "$URL"

echo "Extracting to ${LINEAGES_DIR}/..."
tar -xzf "$TARBALL" -C "$LINEAGES_DIR"

if [ ! -d "${LINEAGES_DIR}/actinopterygii_odb10" ]; then
    echo "ERROR: expected extracted folder ${LINEAGES_DIR}/actinopterygii_odb10 not found"
    exit 1
fi

rm -f "$TARBALL"
echo "Done. Lineage dataset staged at: ${LINEAGES_DIR}/actinopterygii_odb10"
