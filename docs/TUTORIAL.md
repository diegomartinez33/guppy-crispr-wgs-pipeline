# Tutorial: Running the Parameterized Pipelines for Your Own Gene/Guide

This guide is for colleagues who want to use the two pipelines built for direct reuse on
**new genes**, without touching any code: **CRISPR guide design** (KO/CRISPRi) and **PCR primer
design**. Both already accept the gene and genome version as command-line parameters. For the
technical detail of each step see [PIPELINE.md](PIPELINE.md); for where results already live,
[RESULTS.md](RESULTS.md); for cluster access and the required base files,
[CLUSTER_ACCESS.md](CLUSTER_ACCESS.md).

## 0. Prerequisites (once per user)

1. Access to the shared cluster account — see [CLUSTER_ACCESS.md](CLUSTER_ACCESS.md).
2. The base genomes (`reference/`) and the CRISPOR container
   (`codes/analysis/crispor_singularity/`) must exist in your copy of the repository — they are
   the only heavy files copied to the shared account, precisely so this tutorial works without
   needing to regenerate them.
3. Only for primer design: install `primer3_core` once (not shipped with the EMBOSS module):
   ```bash
   bash codes/analysis/setup_primer3.sh
   ```
   This creates the `primer3_env` conda environment. No need to repeat this unless it's deleted.

## 1. Choosing a Reference Genome Version

The whole pipeline supports two genomes: `v1` (Trinidad/Guanapo, female, 2014 — used in all of
this project's historical results) and `v2` (male, PacBio+Hi-C, 2025 — the current RefSeq genome,
migration in progress). Select it with the `REF_VERSION` environment variable:

```bash
# Default (unspecified) uses v1 - identical paths to how the project has always worked
sbatch codes/mapping/bwa_index.sh

# To explicitly use v2:
REF_VERSION=v2 sbatch codes/mapping/bwa_index.sh
```

`REF_VERSION` is automatically propagated to every script that does
`source codes/genome_versions.sh` (most of the pipeline — see the "Dual-genome" column in
[PIPELINE.md](PIPELINE.md) for which ones). No script needs editing to switch versions.

## 2. Guide Design for a New Gene

The guide scan (`ko_guide_scan.py` for CRISPRko, `crispri_tss_scan.py` for CRISPRi) only needs
the gene symbol — it looks up the gene by name (`gene=NAME;`) in both the reference's and the
pseudogenome's GFF, so it must exist in both annotations.

```bash
module load minimap2
module load samtools/1.16.1   # order matters: minimap2 first, then samtools
module load singularity/3.7.1

# CRISPRko - cutting candidates in the CDS
python3 codes/analysis/ko_guide_scan.py --gene my_gene --population pseudogenome

# CRISPRi - silencing candidates near the TSS
python3 codes/analysis/crispri_tss_scan.py --gene my_gene --population pseudogenome
```

Parameters:
- `--gene` (required): gene symbol exactly as it appears in the GFF (e.g. `bdnf`, `grin1a`).
- `--population`: `pseudogenome` (recommended, preserves exon/intron structure — see
  [PIPELINE.md §7](PIPELINE.md#7-colombian-population-genomes)), `scaffolded` (de novo assembly),
  or `pseudogenome_v2` (once it exists, see [RESULTS.md](RESULTS.md)).
- `--no-crispor`: skips the additional CRISPOR scoring (the manual PAM scan and population
  variant classification always run, with or without this flag).

To run several genes in a row, edit the array at the top of
`codes/analysis/run_ko_guide_scan.sh`:
```bash
GENES=(bdnf agap3 grin1a grin1b gria1a gria1b gria2b nlgn1 my_new_gene)
```
and run it with `bash codes/analysis/run_ko_guide_scan.sh`.

**Output:** `analysis/ko_guide_scan/<gene>_<population>_guide_comparison.csv` (CRISPRko) and
`_crispri_candidates.csv` (CRISPRi), plus CRISPOR's raw TSVs where available. To regenerate the
consolidated visual report (like the Guppy CRISPR Atlas) after adding a new gene:
`python3 codes/analysis/build_guide_report.py` (edit the `GENES` list at the top of the script
first).

**If your new gene needs a genome not yet registered with CRISPOR** (uncommon — the project's 3
genomes, `guppyRefTrinidad`/`guppyColPseudogenome`/`guppyRefMaleV2`, already cover
v1/pseudogenome/v2), register it first with `codes/analysis/crispor_add_genomes.sh` (or
`crispor_add_genome_v2.sh` as a pattern reference).

## 3. Primer Design for a New Gene

Requires a sites CSV in the same format `combine_offtargets.py` produces (columns:
`chromosome,start,end,strand,offtarget_seq,mismatches,mit_score,cfd_score,locus,source,found_by_both`).
If your gene has already gone through off-target discovery
([PIPELINE.md §5](PIPELINE.md#5-off-target-discovery-and-analysis)), that CSV already exists at
`crispresso${OUT_SUFFIX}/offtargets/combined/combined_offtargets.csv`.

```bash
sbatch --export=ALL,GENE=my_gene,SITES_CSV=/path/to/my_gene_offtargets.csv,REF_VERSION=v1 \
  codes/analysis/run_offtarget_primer_design.sh
```

Or directly in Python (useful for fine-tuning parameters without editing the wrapper):
```bash
module load emboss/6.6.0 minimap2 samtools/1.16.1
export EMBOSS_PRIMER3_CORE=<path_to>/envs/primer3_env/bin/primer3_core

python3 codes/analysis/design_offtarget_primers.py \
  --gene my_gene --sites-csv my_gene_offtargets.csv --ref-version v1 \
  --window-halfsize 500 --exclude-halfsize 75 \
  --product-min 200 --product-max 400 \
  --opt-tm 60 --min-tm 58 --max-tm 62 --min-gc 40 --max-gc 60 \
  --num-return 5 --specificity-mismatch-pct 10
```

Most useful parameters to adjust:
- `--window-halfsize`: how much sequence around the cut site is given to `eprimer3` to search
  (default 500bp on each side).
- `--exclude-halfsize`: region around the cut where NO primer may fall (default 75bp) — avoids
  designing a primer right over the site where the indel is expected.
- `--product-min`/`--product-max`: desired amplicon size (default 200-400bp, a good range for gel
  + Illumina/Sanger sequencing).
- `--no-population-check`: skips the check against the pseudogenome (faster, but won't detect if
  a primer falls on a real population variant).

**Output:** `analysis/offtarget_primers/<gene>_<ref_version>_primers.csv` — one row per candidate
pair, with specificity (`primersearch`) and population-variant columns at the primer's binding
site. Always check the `design_status` column: `OK` = a real candidate,
`INPUT_ERROR_AMBIGUOUS_BASES` = the site has ambiguous IUPAC codes in the reference and nothing
could be designed (see the limitation in [PIPELINE.md §9](PIPELINE.md#9-pcr-primer-design)),
`NO_CANDIDATES_FOUND` = a search was run but no candidate met the parameters (try relaxing
`--min-tm`/`--max-tm`/`--min-gc`/`--max-gc`).

## 4. Verifying Existing RT-qPCR Primers

If you already have RT-qPCR primers (from a paper, a colleague, or a primer-design tool like
Primer3Plus) and want to check they still work correctly for the Colombian population — and
that they actually target real mRNA, not just genomic DNA — use
`codes/analysis/verify_rtqpcr_primers.py`. Full method and rationale:
[PIPELINE.md §10](PIPELINE.md#10-rt-qpcr-primer-verification).

```bash
module load blast/2.14.1+ emboss/6.6.0 samtools/1.16.1

# 1. Create a CSV with your primer pairs
cat > my_primers.csv << 'EOF'
pair_name,forward_seq,reverse_seq
my_gene,ACGT...,ACGT...
EOF

# 2. Run the checklist against both genome versions + the pseudogenome
python3 codes/analysis/verify_rtqpcr_primers.py \
  --primers-csv my_primers.csv \
  --ref-versions v1,v2 \
  --population-check
```

**Output:** `analysis/rtqpcr_verification/rtqpcr_primer_verification.csv`, one row per
primer×genome-version. Read the `status` column: `OK` = exact match inside a real exon in this
genome version; `WARN_PARTIAL_MATCH` = the closest real gene match isn't 100% identical (check
`mrna_identity_pct`/`mrna_gaps` — a handful of bp off, especially away from the 3' end, may still
amplify but is worth flagging); `GENOMIC_MATCH_NO_GENE` = the primer matches genomic DNA but not
inside any annotated gene (will not amplify from cDNA); `NO_MATCH` = no usable hit at all in that
genome version. For any `OK` row from v1, check `population_status`: `IDENTICAL` = safe;
`VARIANT_FOUND` = a Colombian population variant sits inside this primer's binding site (see
`population_note` for the exact position/base change).

**If you get `GENOMIC_MATCH_NO_GENE` or a low `mrna_identity_pct`, don't assume the primer is
wrong before checking the gene it landed on** — genes can be missing or fragmented in one genome
version and fine in the other (re-run with `--ref-versions v2` alone, or vice versa), and a
primer designed from a slightly different transcript build can still be a real, usable match to
the intended gene despite a few bp of difference (see the myosin heavy chain finding in
`CLAUDE.md` for a worked example).

## 5. Frequently Asked Questions

**Can I run this on a CRISPR guide I already picked, instead of a whole gene?** Both pipelines
work at the gene level (they search for candidates themselves within the CDS/TSS window). If you
already have a specific guide selected and only want its off-targets/primers, use
`casoffinder.sh` + `crispor_offtarget_scan.sh` directly with your guide sequence (see
[PIPELINE.md §5](PIPELINE.md#5-off-target-discovery-and-analysis)) to generate your own
`combined_offtargets.csv`, then follow step 3 of this tutorial from there.

**What happens if my gene isn't in the reference GFF?** `ko_guide_scan.py`/
`crispri_tss_scan.py` fail with a clear error (`ERROR: gene 'X' not found in reference GFF`).
Verify the exact symbol with `grep "gene=" reference/GCF_*_annotation.gff | grep -i my_gene`.

**Do I need to run the whole WGS pipeline (mapping, GATK) to use these two scripts?** No — both
only need the files in `reference/` (genomes + annotations + pseudogenome) and, for guides, the
CRISPOR container. They don't depend on any sample's BAMs or VCFs.
