# Results Map

Where each already-generated result lives, organized by objective. ✅ = complete, 🔄 = in
progress, ⏳ = pending. For detail on how each thing was generated see [PIPELINE.md](PIPELINE.md).

## 1. Off-target WGS Analysis (GATK + CRISPResso2) — ✅ complete (v1)

| What | Path |
|---|---|
| Per-sample variant summary (counts, Ti/Tv) | `codes/analysis/gatk_summary/gatk_variant_summary.csv` + `gatk_summary_barplot.png`, `gatk_titv_boxplot.png` |
| Genotypes at the 8 off-target sites (only OT4 has variants, pre-existing in Control) | `codes/analysis/gatk_summary/offtarget_genotypes.csv` + `offtarget_genotype_heatmap.png` |
| % on-target vs off-target editing per sample/group | `codes/analysis/editing_comparison/editing_summary.csv` + `editing_heatmap.png`, `ontarget_barplot.png`, `offtarget_dotplot.png` |
| Final list of the 8 off-targets (coordinates, MIT/CFD, locus) | `crispresso/offtargets/combined/combined_offtargets.csv` + `combined_offtargets_igv.bed` |
| Final filtered VCFs (SNP/INDEL, whole genome) | `gatk/trimmomatic/vcf_filtered/{snps,indels}_filtered.vcf.gz` |
| VCF restricted to the 8 off-target loci | `gatk/trimmomatic/vcf_offtargets/offtarget_variants.vcf.gz` |
| Per-sample CRISPResso2 reports (on-target, 15 folders + `merged/`) | `crispresso/ontarget/trimmomatic/<SAMPLE>/CRISPResso_on_<SAMPLE>/` |
| Per-sample CRISPRessoWGS reports (8 off-target sites, 15 folders + `aggregate/`) | `crispresso/wgs/trimmomatic/<SAMPLE>/` |
| **Consolidated visual report** | [`analysis/reports/offtarget_wgs_report.html`](../analysis/reports/offtarget_wgs_report.html) |

**Main finding:** no CRISPR-induced indel detected at any of the 8 off-target sites by GATK; the
only site with variants (OT4) already had them in the Control group — a pre-existing population
polymorphism, not an editing effect.

## 2. Variant Hotspots — ✅ complete, v1 only

| What | Path |
|---|---|
| Variant density per window (10kb/2kb) | `gatk/trimmomatic/hotspots/window_counts_annotated.csv`, `window_counts.tsv` |
| Final hotspot regions (403 merged, FDR<0.05) | `gatk/trimmomatic/hotspots/hotspots.bed` |
| Hotspot gene annotation (incl. zebrafish orthologs, gProfiler enrichment) | `hotspot_gene_summary.tsv`, `hotspot_gene_overlaps.tsv`, `hotspot_genes_zebrafish.txt`, `gProfiler_*.csv` |
| Genome-wide Manhattan plot + per-chromosome density | `hotspot_manhattan_genome.png`, `hotspot_plots/` |
| 5 final summary figures | `gatk/trimmomatic/hotspots/summary_plots/hotspot_summary_0{1..5}_*.png` |
| **Consolidated visual report** | [`analysis/reports/hotspots_report.html`](../analysis/reports/hotspots_report.html) |

**Not run under v2 yet** — `gatk/trimmomatic_v2/` has no `hotspots/` subfolder (see objective 7).

## 3. Colombian Pseudogenome — ✅ complete (v1)

`reference/pseudogenome/` — genome + all indices (`.fai`, `.dict`, BWA, minimap2) + Liftoff
annotation (99.5% transfer, 26,264 genes, 0 orphans) + `.chain` for exact liftover. Full detail,
method, and limitations: **[reference/pseudogenome/README.md](../reference/pseudogenome/README.md)**.

## 4. Colombian De Novo Assembly — ✅ complete (v1)

`reference/colombian_scaffolded_genome/` — final genome (gap-filled + polished) + indices +
transferred Liftoff annotation (bdnf: coverage=0.945, sequence_ID=0.923). Comparative QC across
the 4 stages (raw → polished → gapfilled → gapfilled+polished) in `assembly/qc_results/`. Full
detail, method, and limitations:
**[reference/colombian_scaffolded_genome/README.md](../reference/colombian_scaffolded_genome/README.md)**.
Comparative visual report for both population genomes (pseudogenome + de novo assembly):
[`analysis/reports/genome_resources_report.html`](../analysis/reports/genome_resources_report.html).

## 5. CRISPR KO/CRISPRi Guide Design (8 genes) — ✅ complete

`analysis/ko_guide_scan/` — one file set per gene (bdnf, agap3, grin1a, grin1b, gria1a, gria1b,
gria2b, nlgn1): CRISPRko candidate comparison (`*_guide_comparison.csv`), CRISPRi candidates
(`*_crispri_candidates.csv`), full gene-body variants (`*_gene_body_variants.csv`), and CRISPOR's
raw TSVs where available.

**Consolidated report (Guppy CRISPR Atlas):**
[`analysis/ko_guide_scan/report/guppy_crispr_atlas.html`](../analysis/ko_guide_scan/report/guppy_crispr_atlas.html)
— also published as an Artifact (shareable link on request).

**Limitation:** agap3, grin1a, gria1a have no CRISPOR scores (ambiguous IUPAC codes in the v1
reference crash `crispor.py`) — will be repeated against v2 once available. Detail in
[PIPELINE.md §8](PIPELINE.md#8-crispr-guide-design-ko--crispri-per-gene).

## 6. PCR Primer Design — 🔄 partial (bdnf/v1 only)

| What | Path |
|---|---|
| Final primer table (9 sites: on-target + 8 off-target) | `analysis/offtarget_primers/bdnf_v1_primers.csv` |
| Raw `eprimer3` output per site | `analysis/offtarget_primers/raw/bdnf_v1/` |
| **Visual report** | [`analysis/reports/primer_design_report.html`](../analysis/reports/primer_design_report.html) |

**Pending:** the other 7 candidate genes (agap3, grin1a, grin1b, gria1a, gria1b, gria2b, nlgn1)
and the v2 version — the script is already parameterized (`--gene`/`--ref-version`), it just
needs to be run (see [TUTORIAL.md §3](TUTORIAL.md#3-primer-design-for-a-new-gene)).

**Verification of existing RT-qPCR primers (2026-09-08)** — a one-off query, not part of the
automated pipeline: 4 primer pairs already in use in the lab (bdnf + housekeeping
myosin/beta-actin/rpl13a) verified against v1, v2, and the Colombian pseudogenome. Main finding:
**a real Colombian SNP in `rpl_13a_F`** (also confirmed against v2, which matches the primer's
original allele) and **`miosina_guppy_F/R` with no identifiable binding site** in the guppy genome
(possibly a primer designed for a different species). Full detail in the "RT-qPCR" section of the
[primer report](../analysis/reports/primer_design_report.html) and in `CLAUDE.md`
("RT-qPCR Primer Verification").

## 7. Migration to the v2 Reference Genome — 🔄 in progress

| Stage | Status |
|---|---|
| Off-target discovery (Cas-OFFinder + CRISPOR) for bdnf | ✅ done, cross-validated (8/8 match) — `crispresso_v2/offtargets/combined/combined_offtargets.csv` |
| Mapping (BWA) + MarkDuplicates + HaplotypeCaller (15 samples) | 🔄 in progress — `gatk/trimmomatic_v2/{markdup,gvcf}/` exist; HaplotypeCaller running |
| GenomicsDBImport / GenotypeGVCFs / VariantFiltration | ⏳ pending (depends on the above) |
| v2 pseudogenome | ⏳ pending (needs the filtered VCF from above) |
| De novo assembly re-scaffolded against v2 (RagTag Phase 2) | ⏳ pending (deliberately postponed until the above finishes) |
| CRISPResso on-target/off-target under v2 | ⏳ pending — neither `crispresso_v2/ontarget/` nor `wgs/` exist yet |
| Hotspots under v2 | ⏳ pending |
| bdnf guide/primer design already supports `--ref-version v2`/`--population pseudogenome_v2` | ✅ code ready, waiting on the v2 pseudogenome |

## 8. IGV Files — ✅ complete, v1 only

`igv_files/` — pseudogenome genome + annotation (bgzip+tabix) + per-group merged BAMs (4) +
`features_of_interest.bed` (bdnf, guide site, cut site, 8 off-targets). Ready to load directly
into IGV Desktop. No v2 equivalent yet (consistent with objective 7).

---

## Visual Reports — Summary

Each report exists in two forms: published as an Artifact (shareable link, private by default —
share it from the page's own menu) and as a self-contained HTML copy in the repository (to send
directly to a colleague without needing a claude.ai account).

| Report | Objective | Artifact | Local copy |
|---|---|---|---|
| Guppy CRISPR Atlas | KO/CRISPRi guide design, 8 genes | (share on request) | [`analysis/ko_guide_scan/report/guppy_crispr_atlas.html`](../analysis/ko_guide_scan/report/guppy_crispr_atlas.html) |
| Off-Target WGS Report | Off-target WGS (GATK + CRISPResso) | [link](https://claude.ai/code/artifact/3290263b-0f56-4d76-84fe-825d4d98110c) | [`analysis/reports/offtarget_wgs_report.html`](../analysis/reports/offtarget_wgs_report.html) |
| Variant Hotspots Report | Variant hotspots | [link](https://claude.ai/code/artifact/bd831a9e-276a-4aef-a8ef-d35be58f539c) | [`analysis/reports/hotspots_report.html`](../analysis/reports/hotspots_report.html) |
| Colombian Genome Resources | Pseudogenome + de novo assembly | [link](https://claude.ai/code/artifact/beb5bc85-ac67-4f3f-912f-ab4dc82a0d5d) | [`analysis/reports/genome_resources_report.html`](../analysis/reports/genome_resources_report.html) |
| Primer Design Report | Primer design (bdnf/v1) | [link](https://claude.ai/code/artifact/eb740fdb-261b-41a5-b44b-e8530a82c215) | [`analysis/reports/primer_design_report.html`](../analysis/reports/primer_design_report.html) |
