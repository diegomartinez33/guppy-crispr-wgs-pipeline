# Results Map

Where each already-generated result lives, organized by objective. ✅ = complete, 🔄 = in
progress, ⏳ = pending. For detail on how each thing was generated see [PIPELINE.md](PIPELINE.md).

## 1. Off-target WGS Analysis (GATK + CRISPResso2) — ✅ complete, both v1 and v2

| What | v1 path | v2 path |
|---|---|---|
| Per-sample variant summary (counts, Ti/Tv) | `codes/analysis/gatk_summary/gatk_variant_summary.csv` + `gatk_summary_barplot.png`, `gatk_titv_boxplot.png` | `codes/analysis/gatk_summary_v2/` (same files) |
| Genotypes at the 8 off-target sites (only OT4 has variants, pre-existing in Control) | `codes/analysis/gatk_summary/offtarget_genotypes.csv` + `offtarget_genotype_heatmap.png` | `codes/analysis/gatk_summary_v2/offtarget_genotypes.csv` + heatmap — same finding (2 SNPs at OT4, pre-existing in Control) |
| % on-target vs off-target editing per sample/group | `codes/analysis/editing_comparison/editing_summary.csv` + `editing_heatmap.png`, `ontarget_barplot.png`, `offtarget_dotplot.png` | `codes/analysis/editing_comparison_v2/` (same files) — matches v1 exactly per sample |
| Final list of the 8 off-targets (coordinates, MIT/CFD, locus) | `crispresso/offtargets/combined/combined_offtargets.csv` + `combined_offtargets_igv.bed` | `crispresso_v2/offtargets/combined/combined_offtargets.csv` — 8/8 coordinates + scores match v1 (remapped) |
| Final filtered VCFs (SNP/INDEL, whole genome) | `gatk/trimmomatic/vcf_filtered/{snps,indels}_filtered.vcf.gz` | `gatk/trimmomatic_v2/vcf_filtered/{snps,indels}_filtered.vcf.gz` (12.9M SNPs, 3.4M indels) |
| VCF restricted to the 8 off-target loci | `gatk/trimmomatic/vcf_offtargets/offtarget_variants.vcf.gz` | `gatk/trimmomatic_v2/vcf_offtargets/offtarget_variants.vcf.gz` |
| Per-sample CRISPResso2 reports (on-target, 15 folders + `merged/`) | `crispresso/ontarget/trimmomatic/<SAMPLE>/CRISPResso_on_<SAMPLE>/` | `crispresso_v2/ontarget/trimmomatic/<SAMPLE>/` (15 + `merged/`) |
| Per-sample CRISPRessoWGS reports (8 off-target sites, 15 folders + `aggregate/`) | `crispresso/wgs/trimmomatic/<SAMPLE>/` | `crispresso_v2/wgs/trimmomatic/<SAMPLE>/` (15 + `aggregate/`) |
| **Consolidated visual report** | [`analysis/reports/offtarget_wgs_report.html`](../analysis/reports/offtarget_wgs_report.html) — v1/v2 toggle | (same report, `v2` tab) |

**Main finding:** no CRISPR-induced indel detected at any of the 8 off-target sites by GATK; the
only site with variants (OT4) already had them in the Control group — a pre-existing population
polymorphism, not an editing effect.

**Cross-validation (2026-09-20):** this 8-site off-target list is independently reproduced by an
entirely separate pipeline — `codes/analysis/ko_guide_scan.py` (built for objective 5's guide
*design*, not off-target *discovery*) scans every candidate guide in bdnf's CDS and necessarily
includes the guide actually used in the lab experiment among its candidates. Diffing its
CRISPOR off-target output for that exact guide against `combined_offtargets.csv` here: 8/8
coordinates and MIT/CFD scores match exactly, for both v1 and v2. Detail in CLAUDE.md §8.

## 2. Variant Hotspots — ✅ complete, both v1 and v2

| What | v1 path | v2 path |
|---|---|---|
| Variant density per window (10kb/2kb) | `gatk/trimmomatic/hotspots/window_counts_annotated.csv`, `window_counts.tsv` | `gatk/trimmomatic_v2/hotspots/window_counts_annotated.csv`, `window_counts.tsv` |
| Final hotspot regions (403 merged, FDR<0.05) | `gatk/trimmomatic/hotspots/hotspots.bed` | `gatk/trimmomatic_v2/hotspots/hotspots.bed` (459 merged, FDR<0.05) |
| Hotspot gene overlap (reconstructed as `codes/analysis/hotspot_gene_overlap.sh`, 2026-09-15) | `hotspot_gene_overlaps.tsv`, `hotspot_gene_list.txt` (705 genes), `hotspot_geneIDs_all.txt` | `hotspot_gene_overlaps.tsv`, `hotspot_gene_list.txt` (882 genes), `hotspot_geneIDs_all.txt` |
| Zebrafish orthologs (g:Orth, manual/external tool) + gene summary (reconstructed as `codes/analysis/hotspot_gene_summary.py`, 2026-09-15) | `hotspot_gene_summary.tsv` (349/403 regions with ≥1 gene), `hotspot_zebrafish_ENSDARG.txt` (315 orthologs), `gProfiler_*.csv` | `hotspot_gene_summary.tsv` (434/459 regions with ≥1 gene), `hotspot_zebrafish_ENSDARG.txt` (314 orthologs), `gProfiler_*.csv` — ✅ done |
| Genome-wide Manhattan plot + per-chromosome density | `hotspot_manhattan_genome.png`, `hotspot_plots/` (573 PNGs) | `hotspot_manhattan_genome.png`, `hotspot_plots/` (182 PNGs) |
| 5 final summary figures | `gatk/trimmomatic/hotspots/summary_plots/hotspot_summary_0{1..5}_*.png` | `gatk/trimmomatic_v2/hotspots/summary_plots/hotspot_summary_0{1..5}_*.png` |
| **Consolidated visual report** | [`analysis/reports/hotspots_report.html`](../analysis/reports/hotspots_report.html) — v1/v2 toggle, both tabs have real data | (same report, `v2` tab) |

**2026-09-15**: found and fixed two real bugs while running this for v2 (see CLAUDE.md §8) — a
broken `conda activate fastp_env` that silently fell through to a system Python, and a
pandas/matplotlib version incompatibility that crashed `hotspot_analysis.py`'s per-chromosome
plots and 2 of `plot_hotspot_summary.py`'s 5 plots. v1's `hotspot_plots/` had actually been empty
since 2026-06-03 (never worked); re-running the fixed scripts for v1 reproduced the documented
numbers exactly (403 regions, byte-identical CSV/BED) while finally populating those plots — no
factual numbers changed, so the existing v1 report/Artifact needed no update.

## 3. Colombian Pseudogenome — ✅ complete, both v1 and v2

`reference/pseudogenome/` (v1) — genome + all indices (`.fai`, `.dict`, BWA, minimap2) + Liftoff
annotation (99.5% transfer, 26,264 genes, 0 orphans) + `.chain` for exact liftover.
`reference/pseudogenome_v2/` (v2) — same file set, 12/12 verification checks passed, Liftoff
annotation 99.6% transfer (31,226 genes). Full detail, method, and limitations:
**[reference/pseudogenome/README.md](../reference/pseudogenome/README.md)**.

## 4. Colombian De Novo Assembly — ✅ complete, both v1 and v2

`reference/colombian_scaffolded_genome/` (v1) and `reference/colombian_scaffolded_genome_v2/`
(v2) — final genome (gap-filled + polished) + indices + transferred Liftoff annotation
(bdnf: coverage=0.960/sequence_ID=0.957 v1, 0.980/0.964 v2). Comparative QC across the 4 stages
(raw → polished → gapfilled → gapfilled+polished) in `assembly/qc_results/` for v1; v2 was QC'd
only at the final stage (see objective 7's table for the direct v1-vs-v2 BUSCO/QUAST comparison).
Full detail, method, and limitations:
**[reference/colombian_scaffolded_genome/README.md](../reference/colombian_scaffolded_genome/README.md)**.
Comparative visual report for both population genomes (pseudogenome + de novo assembly), now
including the full v1-vs-v2 final-assembly comparison:
[`analysis/reports/genome_resources_report.html`](../analysis/reports/genome_resources_report.html).

**Corrected 2026-09-17:** the previous bdnf figures (coverage=0.945, sequence_ID=0.923) were
computed against the raw RagTag scaffold, not the actually-adopted gap-filled+polished assembly —
that final genome was validated back in September but never promoted to the canonical
`colombian_scaffolded.fna`, so its annotation/indices never existed until now. See CLAUDE.md §8,
"v1 de novo assembly promotion bug".

## 5. CRISPR KO/CRISPRi Guide Design (8 genes) — ✅ complete

`analysis/ko_guide_scan/` — one file set per gene (bdnf, agap3, grin1a, grin1b, gria1a, gria1b,
gria2b, nlgn1): CRISPRko candidate comparison (`*_guide_comparison.csv`), CRISPRi candidates
(`*_crispri_candidates.csv`), full gene-body variants (`*_gene_body_variants.csv`), and CRISPOR's
raw TSVs where available.

**Consolidated report (Guppy CRISPR Atlas):**
[`analysis/ko_guide_scan/report/guppy_crispr_atlas.html`](../analysis/ko_guide_scan/report/guppy_crispr_atlas.html)
— also published as an Artifact: [link](https://claude.ai/code/artifact/71120b70-10ea-4817-b4d6-883a4a1b8856).
**Now has a v1/v2 toggle** (added 2026-09-20, same masthead button pattern as the other 3 toggled
reports) — `report_data.json` (v1) and `report_data_pseudogenome_v2.json` (v2) both embedded.

**Limitation (v1 only):** agap3, grin1a, gria1a have no CRISPOR scores (ambiguous IUPAC codes in
the v1 reference crash `crispor.py`). **Resolved under v2** (2026-09-17): all 3 genes now have
100% CRISPOR coverage (559/559, 442/442, 387/387 guide-comparison rows respectively) —
`analysis/ko_guide_scan/{agap3,grin1a,gria1a}_pseudogenome_v2_guide_comparison.csv`. The Atlas's
v2 toggle now shows this directly (`crispor_available: true` for all 8 genes under v2, vs 5/8
under v1). Detail in [PIPELINE.md §8](PIPELINE.md#8-crispr-guide-design-ko--crispri-per-gene).

**Cross-validation with objective 1 (2026-09-20):** bdnf's own candidate list necessarily
includes the guide already used in the lab experiment (spacer `TGAGAGACGCCCCGGGCATG`, PAM
`CGG` — v1 CDS pos 137, v2 CDS pos 167). Its CRISPOR off-target output here matches objective
1's dedicated `combined_offtargets.csv` exactly (8/8 coordinates + MIT/CFD scores, both
versions) — see objective 1's note and CLAUDE.md §8 for detail.

**Bug found + fixed while rebuilding the Atlas (2026-09-20):** `ko_guide_scan.py`/
`crispri_tss_scan.py` wrote their reference-side CRISPOR TSVs to a filename keyed only on gene
name (e.g. `bdnf_reference_crispor_guides.tsv`), with no reference-version tag — so the 2026-09-17
v2 run had silently overwritten v1's reference-side files with v2 data (population-side files and
the already-published v1 Atlas were unaffected, since both predate the overwrite). Fixed by keying
the filename suffix on `REF_FASTA` identity (`_v2` when the v2 genome is in use); the v1 guide +
CRISPRi scan was then re-run cleanly to regenerate proper v1 reference-side files. See CLAUDE.md
§8 for the full root-cause writeup.

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
| Final primer table (9 sites: on-target + 8 off-target) | `analysis/offtarget_primers/bdnf_v1_primers.csv` (7/9 sites — 2 blocked by IUPAC ambiguity) | `analysis/offtarget_primers/bdnf_v2_primers.csv` (**9/9 sites** — v2's cleaner assembly has no ambiguous bases at those 2 windows; 45 pairs, 36/45 genome-wide specific, all 45 free of population variants) |
| Raw `eprimer3` output per site | `analysis/offtarget_primers/raw/bdnf_v1/` | `analysis/offtarget_primers/raw/bdnf_v2/` |
| **Visual report** | [`analysis/reports/primer_design_report.html`](../analysis/reports/primer_design_report.html) — v1/v2 toggle | |

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

**Population-variant check closed for v2 (2026-09-21):** `population_check()` was hardcoded to
`v1` only (the v2 pseudogenome didn't exist when it was written) even though
`REF_BY_VERSION`/`PSEUDOGENOME_BY_VERSION`/`CHAIN_BY_VERSION` — imported from
`design_offtarget_primers.py` — already had real v2 entries once `reference/pseudogenome_v2/`
was built (2026-09-15); the restriction was just never lifted, so every v2 row silently got
`population_status=not_checked`. Generalized to use `row["ref_version"]` and re-ran: all 16 v2
rows now have a real result (15 `IDENTICAL`, 1 `VARIANT_FOUND` — the same `rpl_13a_original_F`
SNP as v1, at the same position, confirming the finding below).

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

## 7. Migration to the v2 Reference Genome — ✅ complete

| Stage | Status |
|---|---|
| Off-target discovery (Cas-OFFinder + CRISPOR) for bdnf | ✅ done, cross-validated (8/8 match) — `crispresso_v2/offtargets/combined/combined_offtargets.csv` |
| Mapping (BWA) + MarkDuplicates + HaplotypeCaller (15 samples) | ✅ done — all 15 samples completed |
| GenomicsDBImport / GenotypeGVCFs / VariantFiltration | ✅ done — `gatk/trimmomatic_v2/vcf_filtered/` (12.9M SNPs, 3.4M indels, all 23 chromosomes) |
| v2 pseudogenome | ✅ done, 12/12 verification checks passed — `reference/pseudogenome_v2/` (fna + chain + BWA/GATK/minimap2/BLAST indices + Liftoff annotation, 31,226 genes, 99.6% transfer) |
| De novo assembly re-scaffolded against v2 (RagTag Phase 2) | ✅ done — `assembly/ragtag_output_v2/ragtag.scaffold.fasta` (474,480 sequences placed, 611.6Mb, vs v1's 468,519/606.2Mb — v2 places slightly more, consistent with its better contiguity) |
| CRISPResso on-target/off-target under v2 | ✅ done — individual (15), WGS (15×8 sites), and merged (4 groups, authoritative) all completed; merged track exactly reproduces v1's historical editing numbers (Control/RNP_Cas/Plasmid_Ko 0%, Only_MNP 1.47% noise) — `crispresso_v2/ontarget/`, `crispresso_v2/wgs/` |
| CRISPResso post-processing (aggregate + group-compare) under v2 | ✅ done (2026-09-20) — `crispresso_aggregate_ontarget.sh`, `crispresso_compare_merged.sh`, `crispresso_wgs_aggregate.sh` were real, done-for-v1 steps still hardcoded to v1 paths (flagged "v1-only" in PIPELINE.md but never actually closed) — parameterized for `REF_VERSION` and run for v2; v1 re-run first as a regression check (byte-identical folder structure/counts). v2's `RNP_Cas_vs_Control` comparison shows 0% substitutions in both groups, consistent with the no-CRISPR-induced-editing finding — `crispresso_v2/aggregate/`, `crispresso_v2/compare/trimmomatic/merged/`, `crispresso_v2/wgs/trimmomatic/aggregate/` |
| Hotspots under v2 (incl. zebrafish orthologs) | ✅ done — 459 merged hotspot regions (FDR<0.05), gene overlap (882 genes), 314 zebrafish orthologs — `gatk/trimmomatic_v2/hotspots/`. Found + fixed 2 real bugs along the way (broken conda env, pandas/matplotlib plotting crash — also affected v1, now fixed there too) — see objective 2 |
| `select_offtargets` genotyping for the 8 v2 off-target sites | ✅ done — 2 SNPs at off_target_4, present in Control (pre-existing, not CRISPR-induced), reproducing v1's exact finding — `gatk/trimmomatic_v2/vcf_offtargets/` |
| GATK genotype/summary plots for v2 (`gatk_offtarget_genotypes.py`, `gatk_variant_summary.py`, `plot_editing_comparison.py`) | ✅ done — all 3 scripts parameterized (`REF_VERSION`) and run for v2; off-target genotypes exactly reproduce v1 (same 2 background SNPs at OT4, remapped coordinates); on-target editing numbers match v1 exactly per sample; genome-wide summary (~8-9M SNPs/sample, Ti/Tv 1.36-1.37) tracks v1's corrected numbers — `codes/analysis/gatk_summary_v2/`, `codes/analysis/editing_comparison_v2/`. Side finding: this run also caught that the *committed v1* `gatk_summary/gatk_variant_summary.csv` was stale (computed before the final hard-filtered VCF existed, 5.5× undercounted) and has now been refreshed — see CLAUDE.md §8 |
| `guppyColPseudogenomeV2` registered with CRISPOR | ✅ done — `codes/analysis/crispor_add_genome_v2.sh` extended to also register the v2 pseudogenome (mirrors v1's `crispor_add_genomes.sh`); verified (2bit + BWA index present) — unblocks CRISPOR scoring for `ko_guide_scan.py`/`crispri_tss_scan.py --population pseudogenome_v2` |
| bdnf primers for v2 | ✅ done — 45 candidate pairs (5/site × 9 sites), 36/45 confirmed genome-wide specific, all 45 free of population variants in the primer footprint — `analysis/offtarget_primers/bdnf_v2_primers.csv` |
| IGV files for v2 | ✅ done — `codes/assembly/prepare_igv_files.sh` parameterized (`REF_VERSION`); `igv_files_v2/` (genome, sorted+bgzipped+tabixed annotation, 4 merged BAMs, `features_of_interest.bed`). The off-target/sgRNA coordinates in the BED were lifted from reference to pseudogenome coordinates via CrossMap (`colombian_pseudogenome.chain`) and spot-verified against the actual sequence — this also surfaced that v1's equivalent BED mixes reference- and pseudogenome-native coordinates (pre-existing imprecision, not fixed retroactively) |
| 8-gene KO/CRISPRi guide comparison against pseudogenome_v2 | ✅ done — all 8 genes, both CRISPRko (`*_pseudogenome_v2_guide_comparison.csv`) and CRISPRi (`*_pseudogenome_v2_crispri_candidates.csv`). **agap3/grin1a/gria1a's v1 CRISPOR-scoring limitation is fully resolved**: 0/557, 0/441, 0/333 scored rows under v1 → 559/559, 442/442, 387/387 under v2 (100% coverage) — v2's cleaner assembly has no ambiguous IUPAC bases at these genes' guide windows. Atlas report rebuilt with a v1/v2 toggle (2026-09-20) — see objective 5 |
| bdnf editing-site coverage-by-zone analysis for v2 | ✅ done (2026-09-20) — `get_editing_region_coverage.sh` parameterized (`BAM_DIR`/`OUTPUT_DIR`/`CHROMOSOME`/`SGRNA_START`/`SGRNA_END`), re-run for both versions — see objective 9 |
| De novo assembly gap-filling/polishing/annotation for v2 | ✅ done (2026-09-20) — TGS-GapCloser → NextPolish (8h21m) → Liftoff → indexing → BUSCO → QUAST QC all complete (job 729122, 7h04m). bdnf Liftoff coverage=0.980/sequence_ID=0.964 (v1: 0.960/0.957); BUSCO C:96.3% (S:95.9%,D:0.4%) F:2.0% M:1.7% (v1: C:95.5% F:2.4% M:2.1%) — `reference/colombian_scaffolded_genome_v2/`. QUAST: N50=30.7Mb, NA50=159Kb, misassemblies=13,768, duplication=1.029 — all better than v1 (29.4Mb, 114Kb, 23,914, 1.042); genome fraction is the one metric v1 leads on (92.2% vs 85.5%), explained by a reference-composition effect, not an assembly regression — see the "Reading the genome fraction gap" note in `genome_resources_report.html` and CLAUDE.md §8. Found + fixed 2 real bugs along the way: `tgsgapcloser_genome.sh` never copied its `.scaff_seqs` output to the `.fasta` name downstream tools expect (v1 had this done by hand, never scripted); QUAST itself was initially skipped as "not warranted" — wrong, BUSCO alone can't give the structural metrics the comparison report needs |

**Every pipeline stage above is now done for v2**, including the de novo assembly's QUAST QC
(2026-09-20). The Guppy CRISPR Atlas (objective 5) has its v1/v2 toggle. The Colombian Genome
Resources report (objective 4) now has a full v1-vs-v2 final-assembly comparison (BUSCO + QUAST
side by side, with an explanation of the one metric — genome fraction — where v1 numerically
leads and why that isn't a real quality regression in v2).

## 8. IGV Files — ✅ complete, both v1 and v2

`igv_files{,_v2}/` — pseudogenome genome + annotation (bgzip+tabix) + per-group merged BAMs (4) +
`features_of_interest.bed` (bdnf, guide site, cut site, 8 off-targets). Ready to load directly
into IGV Desktop.

## 9. bdnf Editing-Site Coverage-by-Zone Analysis — ✅ complete, both v1 and v2

A supplementary per-sample depth/coverage breakdown around the bdnf CRISPR cut site (sgRNA site,
±100bp, ±500bp zones), independent of the GATK/CRISPResso objective-1 analysis above. Predates
this project's `REF_VERSION` convention — its scripts were fully hardcoded to v1 until
parameterized on 2026-09-20 (`get_editing_region_coverage.sh`: `BAM_DIR`/`OUTPUT_DIR`/
`CHROMOSOME`/`SGRNA_START`/`SGRNA_END`; `summary_coverage.sh` and its 4 plotting scripts in
`coverage/csv/`: `COVERAGE_DIR`/`OUTPUT_DIR`/`INPUT_CSV`/`CHROMOSOME`/`SGRNA_START`/`SGRNA_END`/
`SGRNA_SEQ`/`REGION_LABEL`), then run once per version.

| What | v1 path | v2 path |
|---|---|---|
| Per-sample raw depth/coverage (15 samples × 3 files) | `coverage/bdnf_site/` | `coverage/bdnf_site_v2/` |
| Aggregated CSVs (4: summary, depth-by-position, depth-by-zone, coverage-by-zone) | `coverage/csv/` | `coverage/csv_v2/` |
| Plots (5 per script × 4 scripts) | `coverage/csv/{coverage_plots,depth_zone_plots,depth_position_plots,coverage_zone_plots}/` | `coverage/csv_v2/{coverage_plots_v2,depth_zone_plots_v2,depth_position_plots_v2,coverage_zone_plots_v2}/` |

**Naming note:** `coverage/bdnf_site_v2/` predates this project's `REF_VERSION` convention and
used to contain v1-coordinate data despite its name (a leftover from before "_v2" meant "second
reference genome" in this project) — `coverage/bdnf_site/` (no suffix) also existed with stale
data from an older script revision (`WINDOW=500` vs the current `WINDOW=600`). Both were
regenerated 2026-09-20 with the corrected, version-appropriate data; the naming collision is now
resolved (`bdnf_site/` = genuine v1, `bdnf_site_v2/` = genuine v2).

**Bug found + fixed (2026-09-20):** `summary_coverage.sh`'s per-position depth CSV builder had an
off-by-one in its `paste`-output field extraction (`for(i=2;i<=NF;i+=2)` instead of `i=3`), which
silently extracted each sample's duplicated *position* column instead of its *depth* column —
`depth_per_position_all_samples.csv` had values in the millions (genomic coordinates) instead of
real depth (~40-70×). Fixed; also fixed a separate plotting bug in `plot_depth_by_position.py`
(a text annotation placed at `y=-1` in axes-fraction coordinates, which made `bbox_inches="tight"`
try to render an ~970-inch-tall canvas).

v2's cut-site coordinates (`NC_088832.1:15849694-15849713`) were already established in an earlier
session via minimap2 liftover from the v1 site — see CLAUDE.md §8.

---

## Visual Reports — Summary

Each report exists in two forms: published as an Artifact (shareable link, private by default —
share it from the page's own menu) and as a self-contained HTML copy in the repository (to send
directly to a colleague without needing a claude.ai account).

| Report | Objective | Artifact | Local copy |
|---|---|---|---|
| Guppy CRISPR Atlas | KO/CRISPRi guide design, 8 genes, v1/v2 toggle | [link](https://claude.ai/code/artifact/71120b70-10ea-4817-b4d6-883a4a1b8856) | [`analysis/ko_guide_scan/report/guppy_crispr_atlas.html`](../analysis/ko_guide_scan/report/guppy_crispr_atlas.html) |
| Off-Target WGS Report | Off-target WGS (GATK + CRISPResso), v1/v2 toggle | [link](https://claude.ai/code/artifact/3290263b-0f56-4d76-84fe-825d4d98110c) | [`analysis/reports/offtarget_wgs_report.html`](../analysis/reports/offtarget_wgs_report.html) |
| Variant Hotspots Report | Variant hotspots, v1/v2 toggle | [link](https://claude.ai/code/artifact/bd831a9e-276a-4aef-a8ef-d35be58f539c) | [`analysis/reports/hotspots_report.html`](../analysis/reports/hotspots_report.html) |
| Colombian Genome Resources | Pseudogenome + de novo assembly, v1-vs-v2 final-assembly comparison | [link](https://claude.ai/code/artifact/beb5bc85-ac67-4f3f-912f-ab4dc82a0d5d) | [`analysis/reports/genome_resources_report.html`](../analysis/reports/genome_resources_report.html) |
| Primer Design Report | Primer design (bdnf), v1/v2 toggle | [link](https://claude.ai/code/artifact/eb740fdb-261b-41a5-b44b-e8530a82c215) | [`analysis/reports/primer_design_report.html`](../analysis/reports/primer_design_report.html) |
| RT-qPCR Primer Verification Report | RT-qPCR primer verification (5 in-use pairs + 3 candidates) | [link](https://claude.ai/code/artifact/7b447aed-e947-42b6-96e5-085bfdaea0e9) | [`analysis/reports/rtqpcr_primer_verification_report.html`](../analysis/reports/rtqpcr_primer_verification_report.html) |
