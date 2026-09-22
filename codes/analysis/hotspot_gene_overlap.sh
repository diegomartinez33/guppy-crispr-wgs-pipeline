#!/bin/bash
# Genes overlapping the merged hotspot regions (hotspots.bed), intersected
# against the population-specific Liftoff annotation (not the raw NCBI GFF -
# the pseudogenome's own liftoff.gff3 carries this population's actual gene
# set/coordinates). For v1 this step was originally done ad hoc/by hand
# (hotspot_gene_overlaps.tsv/hotspot_gene_list.txt/hotspot_geneIDs_all.txt
# have no script behind them in this repo) - reconstructed here as a proper,
# reproducible script by reverse-engineering the exact command from the
# existing v1 output format, so it can be re-run identically for v2 and any
# future reference version. 2026-09-15.
#
# Usage:
#   bash codes/analysis/hotspot_gene_overlap.sh
#   REF_VERSION=v2 bash codes/analysis/hotspot_gene_overlap.sh
#
# Requires: module load bedtools/2.30.0 (bcftools not needed here)

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.." >/dev/null 2>&1 && pwd)}}"
source "${PROJECT_DIR}/codes/genome_versions.sh"

HOT_DIR=${PROJECT_DIR}/gatk/trimmomatic${OUT_SUFFIX}/hotspots
PSEUDO_GFF=${PROJECT_DIR}/reference/pseudogenome${OUT_SUFFIX}/colombian_pseudogenome.liftoff.gff3
HOTSPOTS_BED=${HOT_DIR}/hotspots.bed
OUT_TSV=${HOT_DIR}/hotspot_gene_overlaps.tsv
OUT_GENE_LIST=${HOT_DIR}/hotspot_gene_list.txt
OUT_GENEIDS_ALL=${HOT_DIR}/hotspot_geneIDs_all.txt

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

echo "REF_VERSION=${REF_VERSION}  HOTSPOTS_BED=${HOTSPOTS_BED}  PSEUDO_GFF=${PSEUDO_GFF}"

echo "=== Step 1: strip header from hotspots.bed ==="
tail -n +2 "$HOTSPOTS_BED" > "${TMPDIR}/hotspots_noheader.bed"
echo "Hotspot regions: $(wc -l < ${TMPDIR}/hotspots_noheader.bed)"

echo "=== Step 2: extract gene features from the Liftoff annotation ==="
awk -F'\t' '!/^#/ && $3=="gene"' "$PSEUDO_GFF" > "${TMPDIR}/genes.gff3"
echo "Genes in Liftoff annotation: $(wc -l < ${TMPDIR}/genes.gff3)"

echo "=== Step 3: intersect hotspots with genes ==="
bedtools intersect -a "${TMPDIR}/hotspots_noheader.bed" -b "${TMPDIR}/genes.gff3" -wa -wb \
  > "$OUT_TSV"
echo "Overlap records: $(wc -l < ${OUT_TSV})"
echo "✅ ${OUT_TSV}"

echo "=== Step 4: unique gene symbol list ==="
grep -oP 'Name=\K[^;]+' "$OUT_TSV" | sort -u > "$OUT_GENE_LIST"
echo "Unique gene symbols: $(wc -l < ${OUT_GENE_LIST})"
echo "✅ ${OUT_GENE_LIST}"

echo "=== Step 5: unique NCBI GeneIDs ==="
grep -oP 'Dbxref=GeneID:\K[0-9]+' "$OUT_TSV" | sort -u -n > "$OUT_GENEIDS_ALL"
echo "Unique GeneIDs: $(wc -l < ${OUT_GENEIDS_ALL})"
echo "✅ ${OUT_GENEIDS_ALL}"

echo ""
echo "Done. Next (manual, external tool): upload ${OUT_GENE_LIST} to gProfiler"
echo "(https://biit.cs.ut.ee/gprofiler/) against Danio rerio for ortholog"
echo "mapping + enrichment, matching v1's hotspot_genes_zebrafish.txt/"
echo "hotspot_zebrafish_ENSDARG.txt/hotspot_named_genes_gProfiler.txt/gProfiler_*.csv."
