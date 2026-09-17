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

## 2. Variant Hotspots — ✅ complete, both v1 and v2

| What | v1 path | v2 path |
|---|---|---|
| Variant density per window (10kb/2kb) | `gatk/trimmomatic/hotspots/window_counts_annotated.csv`, `window_counts.tsv` | `gatk/trimmomatic_v2/hotspots/window_counts_annotated.csv`, `window_counts.tsv` |
| Final hotspot regions (403 merged, FDR<0.05) | `gatk/trimmomatic/hotspots/hotspots.bed` | `gatk/trimmomatic_v2/hotspots/hotspots.bed` (459 merged, FDR<0.05) |
| Hotspot gene overlap (reconstructed as `codes/analysis/hotspot_gene_overlap.sh`, 2026-09-15) | `hotspot_gene_overlaps.tsv`, `hotspot_gene_list.txt` (705 genes), `hotspot_geneIDs_all.txt` | `hotspot_gene_overlaps.tsv`, `hotspot_gene_list.txt` (882 genes), `hotspot_geneIDs_all.txt` |
| Zebrafish orthologs (g:Orth, manual/external tool) + gene summary (reconstructed as `codes/analysis/hotspot_gene_summary.py`, 2026-09-15) | `hotspot_gene_summary.tsv` (349/403 regions with ≥1 gene), `hotspot_zebrafish_ENSDARG.txt` (315 orthologs), `gProfiler_*.csv` | `hotspot_gene_summary.tsv` (434/459 regions with ≥1 gene), `hotspot_zebrafish_ENSDARG.txt` (314 orthologs), `gProfiler_*.csv` — ✅ done |
| Genome-wide Manhattan plot + per-chromosome density | `hotspot_manhattan_genome.png`, `hotspot_plots/` (573 PNGs) | `hotspot_manhattan_genome.png`, `hotspot_plots/` (182 PNGs) |
| 5 final summary figures | `gatk/trimmomatic/hotspots/summary_plots/hotspot_summary_0{1..5}_*.png` | `gatk/trimmomatic_v2/hotspots/summary_plots/hotspot_summary_0{1..5}_*.png` |
| **Consolidated visual report** | [`analysis/reports/hotspots_report.html`](../analysis/reports/hotspots_report.html) (v1 numbers unchanged, see note) | not yet built |

**2026-09-15**: found and fixed two real bugs while running this for v2 (see CLAUDE.md §8) — a
broken `conda activate fastp_env` that silently fell through to a system Python, and a
pandas/matplotlib version incompatibility that crashed `hotspot_analysis.py`'s per-chromosome
plots and 2 of `plot_hotspot_summary.py`'s 5 plots. v1's `hotspot_plots/` had actually been empty
since 2026-06-03 (never worked); re-running the fixed scripts for v1 reproduced the documented
numbers exactly (403 regions, byte-identical CSV/BED) while finally populating those plots — no
factual numbers changed, so the existing v1 report/Artifact needed no update.

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
— also published as an Artifact: [link](https://claude.ai/code/artifact/71120b70-10ea-4817-b4d6-883a4a1b8856).

**Limitation:** agap3, grin1a, gria1a have no CRISPOR scores (ambiguous IUPAC codes in the v1
reference crash `crispor.py`) — will be repeated against v2 once available. Detail in
[PIPELINE.md §8](PIPELINE.md#8-crispr-guide-design-ko--crispri-per-gene).

**Isoform constitutivity (2026-09-12):** CRISPRko candidates are now checked for whether their
genomic footprint is present in EVERY annotated isoform of the gene, not just the one
`ko_guide_scan.py` designs against — surfaced as `isoform_coverage_n`/`isoform_coverage_total`/
`is_constitutive` in each `*_guide_comparison.csv`, and used as the top ranking criterion for the
Atlas's recommended candidates. Mostly a minor correction (6/8 genes have an 80-96% constitutive
core), but critical for **agap3**: only 19% of its guide candidates are constitutive across its 7
very differently-structured isoforms — the Atlas's top-5 picks for agap3 now guarantee full
isoform coverage, which the previous position-only ranking did not. Full methodology in
`CLAUDE.md`, item #6, "Isoform constitutivity check".

**CRISPRi TSS coverage (2026-09-14):** same question, extended to CRISPRi — does the chosen
-50/+300bp window around the TSS actually reach every isoform's own transcription start? Most
genes cluster tightly (one window covers all), but **bdnf** (3.8kb TSS spread) and **agap3**
(63.8kb spread) have real, independently RNA-seq-supported alternative promoters, not annotation
noise. The CRISPRi representative transcript is now chosen by strongest RNA-seq support for its
own TSS (not CDS length, which is meaningless for a promoter question) — this changed agap3's
CRISPRi representative to a transcript 63.8kb upstream of the one CRISPRko uses, since it has the
strongest support of all 7 isoforms (44 RNA-seq samples). The Atlas now shows an explicit
coverage note per gene (e.g. "covers only 1/7 isoforms — real alternative promoter, not
addressed"). Full methodology in `CLAUDE.md`, item #6, "Same process extended to CRISPRi".

## 6. PCR Primer Design — 🔄 partial (bdnf only, v1 + v2)

| What | v1 | v2 |
|---|---|---|
| Final primer table (9 sites: on-target + 8 off-target) | `analysis/offtarget_primers/bdnf_v1_primers.csv` | `analysis/offtarget_primers/bdnf_v2_primers.csv` (45 pairs, 36/45 genome-wide specific, all 45 free of population variants) |
| Raw `eprimer3` output per site | `analysis/offtarget_primers/raw/bdnf_v1/` | `analysis/offtarget_primers/raw/bdnf_v2/` |
| **Visual report** | [`analysis/reports/primer_design_report.html`](../analysis/reports/primer_design_report.html) (v1 numbers only) | — |

**Pending:** the other 7 candidate genes (agap3, grin1a, grin1b, gria1a, gria1b, gria2b, nlgn1),
any version — the script is already parameterized (`--gene`/`--ref-version`), it just
needs to be run (see [TUTORIAL.md §3](TUTORIAL.md#3-primer-design-for-a-new-gene)).

## 6b. RT-qPCR Primer Verification — ✅ complete (formalized script)

Verification of existing lab RT-qPCR primers is a separate objective from #6 above — the primers
are already in use (not designed by this project), the target is spliced mRNA rather than genomic
DNA, and the gene is confirmed rather than known up front. Originally a one-off query
(2026-09-08); formalized into a reusable, parameterized script (2026-09-09):

| What | Path |
|---|---|
| Verification script | `codes/analysis/verify_rtqpcr_primers.py` (+ `run_rtqpcr_primer_verification.sh`) |
| Result table (8 pairs × v1/v2) | `analysis/rtqpcr_verification/rtqpcr_primer_verification.csv` |
| **Visual report** | [`analysis/reports/rtqpcr_primer_verification_report.html`](../analysis/reports/rtqpcr_primer_verification_report.html) |

**New housekeeping primer candidates (2026-09-11, not yet in lab use):** since myosin's actual role
is as the reference/housekeeping gene contrasted against `bdnf` expression after KO, 3 new
candidates were designed and verified: `myosin_conserved` (targets a region shared by ≥6 tandem
myosin paralog copies, found via `muscle` MSA + `eprimer3` + `primersearch`, rather than one
specific transcript), and standard single-copy alternatives `gapdh_new`/`ef1a_new`. All 3: 100%
identity in v1 and v2, IDENTICAL in the Colombian pseudogenome. Full design methodology in
`CLAUDE.md` ("RT-qPCR Primer Verification", follow-up 2026-09-11).

5 primer pairs checked (bdnf, beta-actin, rpl13a original, rpl13a replacement, myosin) against v1,
v2, and the Colombian pseudogenome. Main findings: **a real Colombian SNP in the original
`rpl_13a_F`** (also confirmed against v2, which matches the primer's original allele — since
replaced by a verified SNP-free primer) and **`miosina_guppy_F/R` matching a real tandem myosin
heavy chain gene family at 90-100% identity** (the original 2026-09-08 conclusion of "no
identifiable binding site" was a search-coverage gap — guppy's myosin genes carry no short gene
symbol — corrected once BLAST/GFF-overlap-based search replaced the symbol-regex search). Full
detail in the visual report above, in `CLAUDE.md` ("RT-qPCR Primer Verification"), and in method
form in [PIPELINE.md §10](PIPELINE.md#10-rt-qpcr-primer-verification).

## 7. Migration to the v2 Reference Genome — 🔄 in progress

| Stage | Status |
|---|---|
| Off-target discovery (Cas-OFFinder + CRISPOR) for bdnf | ✅ done, cross-validated (8/8 match) — `crispresso_v2/offtargets/combined/combined_offtargets.csv` |
| Mapping (BWA) + MarkDuplicates + HaplotypeCaller (15 samples) | ✅ done — all 15 samples completed |
| GenomicsDBImport / GenotypeGVCFs / VariantFiltration | ✅ done — `gatk/trimmomatic_v2/vcf_filtered/` (12.9M SNPs, 3.4M indels, all 23 chromosomes) |
| v2 pseudogenome | ✅ done, 12/12 verification checks passed — `reference/pseudogenome_v2/` (fna + chain + BWA/GATK/minimap2/BLAST indices + Liftoff annotation, 31,226 genes, 99.6% transfer) |
| De novo assembly re-scaffolded against v2 (RagTag Phase 2) | ✅ done — `assembly/ragtag_output_v2/ragtag.scaffold.fasta` (474,480 sequences placed, 611.6Mb, vs v1's 468,519/606.2Mb — v2 places slightly more, consistent with its better contiguity) |
| CRISPResso on-target/off-target under v2 | ✅ done — individual (15), WGS (15×8 sites), and merged (4 groups, authoritative) all completed; merged track exactly reproduces v1's historical editing numbers (Control/RNP_Cas/Plasmid_Ko 0%, Only_MNP 1.47% noise) — `crispresso_v2/ontarget/`, `crispresso_v2/wgs/` |
| Hotspots under v2 (incl. zebrafish orthologs) | ✅ done — 459 merged hotspot regions (FDR<0.05), gene overlap (882 genes), 314 zebrafish orthologs — `gatk/trimmomatic_v2/hotspots/`. Found + fixed 2 real bugs along the way (broken conda env, pandas/matplotlib plotting crash — also affected v1, now fixed there too) — see objective 2 |
| `select_offtargets` genotyping for the 8 v2 off-target sites | ✅ done — 2 SNPs at off_target_4, present in Control (pre-existing, not CRISPR-induced), reproducing v1's exact finding — `gatk/trimmomatic_v2/vcf_offtargets/` |
| GATK genotype/summary plots for v2 (`gatk_offtarget_genotypes.py`, `gatk_variant_summary.py`, `plot_editing_comparison.py`) | ✅ done — all 3 scripts parameterized (`REF_VERSION`) and run for v2; off-target genotypes exactly reproduce v1 (same 2 background SNPs at OT4, remapped coordinates); on-target editing numbers match v1 exactly per sample; genome-wide summary (~8-9M SNPs/sample, Ti/Tv 1.36-1.37) tracks v1's corrected numbers — `codes/analysis/gatk_summary_v2/`, `codes/analysis/editing_comparison_v2/`. Side finding: this run also caught that the *committed v1* `gatk_summary/gatk_variant_summary.csv` was stale (computed before the final hard-filtered VCF existed, 5.5× undercounted) and has now been refreshed — see CLAUDE.md §8 |
| `guppyColPseudogenomeV2` registered with CRISPOR | ✅ done — `codes/analysis/crispor_add_genome_v2.sh` extended to also register the v2 pseudogenome (mirrors v1's `crispor_add_genomes.sh`); verified (2bit + BWA index present) — unblocks CRISPOR scoring for `ko_guide_scan.py`/`crispri_tss_scan.py --population pseudogenome_v2` |
| bdnf primers for v2 | ✅ done — 45 candidate pairs (5/site × 9 sites), 36/45 confirmed genome-wide specific, all 45 free of population variants in the primer footprint — `analysis/offtarget_primers/bdnf_v2_primers.csv` |
| IGV files for v2 | ✅ done — `codes/assembly/prepare_igv_files.sh` parameterized (`REF_VERSION`); `igv_files_v2/` (genome, sorted+bgzipped+tabixed annotation, 4 merged BAMs, `features_of_interest.bed`). The off-target/sgRNA coordinates in the BED were lifted from reference to pseudogenome coordinates via CrossMap (`colombian_pseudogenome.chain`) and spot-verified against the actual sequence — this also surfaced that v1's equivalent BED mixes reference- and pseudogenome-native coordinates (pre-existing imprecision, not fixed retroactively) |

## 8. IGV Files — ✅ complete, both v1 and v2

`igv_files{,_v2}/` — pseudogenome genome + annotation (bgzip+tabix) + per-group merged BAMs (4) +
`features_of_interest.bed` (bdnf, guide site, cut site, 8 off-targets). Ready to load directly
into IGV Desktop.

---

## Visual Reports — Summary

Each report exists in two forms: published as an Artifact (shareable link, private by default —
share it from the page's own menu) and as a self-contained HTML copy in the repository (to send
directly to a colleague without needing a claude.ai account).

| Report | Objective | Artifact | Local copy |
|---|---|---|---|
| Guppy CRISPR Atlas | KO/CRISPRi guide design, 8 genes | [link](https://claude.ai/code/artifact/71120b70-10ea-4817-b4d6-883a4a1b8856) | [`analysis/ko_guide_scan/report/guppy_crispr_atlas.html`](../analysis/ko_guide_scan/report/guppy_crispr_atlas.html) |
| Off-Target WGS Report | Off-target WGS (GATK + CRISPResso) | [link](https://claude.ai/code/artifact/3290263b-0f56-4d76-84fe-825d4d98110c) | [`analysis/reports/offtarget_wgs_report.html`](../analysis/reports/offtarget_wgs_report.html) |
| Variant Hotspots Report | Variant hotspots, v1/v2 toggle | [link](https://claude.ai/code/artifact/bd831a9e-276a-4aef-a8ef-d35be58f539c) | [`analysis/reports/hotspots_report.html`](../analysis/reports/hotspots_report.html) |
| Colombian Genome Resources | Pseudogenome + de novo assembly | [link](https://claude.ai/code/artifact/beb5bc85-ac67-4f3f-912f-ab4dc82a0d5d) | [`analysis/reports/genome_resources_report.html`](../analysis/reports/genome_resources_report.html) |
| Primer Design Report | Primer design (bdnf/v1) | [link](https://claude.ai/code/artifact/eb740fdb-261b-41a5-b44b-e8530a82c215) | [`analysis/reports/primer_design_report.html`](../analysis/reports/primer_design_report.html) |
| RT-qPCR Primer Verification Report | RT-qPCR primer verification (5 in-use pairs + 3 candidates) | [link](https://claude.ai/code/artifact/7b447aed-e947-42b6-96e5-085bfdaea0e9) | [`analysis/reports/rtqpcr_primer_verification_report.html`](../analysis/reports/rtqpcr_primer_verification_report.html) |
