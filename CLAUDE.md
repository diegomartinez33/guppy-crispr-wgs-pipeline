# CLAUDE.md — Project Context for Claude Code
# Guppy CRISPR-Cas9 WGS Analysis Pipeline
# Author: da.martinez33 | Universidad de los Andes (Colombia)
# Last updated: August 6 2026

---

## Project Overview

Whole genome sequencing (WGS) analysis of CRISPR-Cas9 editing efficiency
and off-target effects in the guppy fish (*Poecilia reticulata*) from a
Colombian population. The target gene is **bdnf** (brain-derived neurotrophic
factor), edited using SpCas9 with a single guide RNA.

---

## Biological Context

| Parameter | Value |
|---|---|
| Species | *Poecilia reticulata* (guppy) |
| Population | Colombia (divergent from reference genome population) |
| Reference genome | GCF_000633615.1 (Guanapo strain, Trinidad) |
| Sequencing platform | Illumina NovaSeq X, 25B flowcell |
| Read type | 2×150bp paired-end |
| Library prep | Illumina DNA Prep (Nextera tagmentation) |
| Adapters R1 | CTGTCTCTTATACACATCT |
| Adapters R2 | ATGTGTATAAGAGACA |

---

## CRISPR Target Information

| Parameter | Value |
|---|---|
| Target gene | bdnf |
| Cas9 variant | SpCas9 NLS (NEB) |
| sgRNA sequence | TGAGAGACGCCCCGGGCATG |
| sgRNA strand | Negative (-) |
| Chromosome | NC_024333.1 (LG3) |
| sgRNA coordinates | 15922039–15922058 |
| Cut site | ~15922041 |
| PAM sequence | CGG (on negative strand) |
| Amplicon region | NC_024333.1:15922011–15922071 (60bp) |
| Amplicon FWD file | reference/amplicon_bdnf_60bp_fwd.fa |
| Amplicon RC file | reference/amplicon_bdnf_60bp_rc.fa |

> **Important:** The sgRNA maps to the **negative strand**. CRISPResso2
> requires both forward and reverse complement amplicons:
> `--amplicon_seq "FWD_SEQ,RC_SEQ" --amplicon_name "bdnf_fwd,bdnf_rc"`

---

## Experimental Groups (15 samples total)

| Group | Samples | n |
|---|---|---|
| Control | Control_MNP_I_S54_L002, Control_MNP_II_S55_L002, Control_MNP_III_S56_L002 | 3 |
| Only_MNP | Only_MNP_C1_S57_L002, Only_MNP_C2_S58_L002, Only_MNP_C3_S59_L002, Only_MNP_C4_S60_L002 | 4 |
| Plasmid_Ko | Plasmid_Ko_P1_S61_L002, Plasmid_Ko_P2_S62_L002, Plasmid_Ko_P3_S63_L002, Plasmid_Ko_P4_S64_L002 | 4 |
| RNP_Cas | RNP_Cas1_S65_L002, RNP_Cas2_S66_L002, RNP_Cas3_S67_L002, RNP_Cas4_S68_L002 | 4 |

Sample list file: `samples.txt` (one sample per line, without extensions)

---

## Cluster Configuration

| Parameter | Value |
|---|---|
| Cluster | hypatia (SLURM) |
| Login node | hypatia |
| Compute nodes | nodei-1 through nodei-10 |
| Default partition | short |
| Storage limit | 2TB per user |
| /tmp size | ~322GB available per node |

### Modules Available
```bash
module load samtools/1.16.1
module load bwa/0.7.17
module load gatk4/4.4.0.0
module load trimmomatic        # version available as module
module load fastqc             # version available as module
```

### Conda Installations
```bash
# System conda (fastp, multiqc, etc.) — use dynamic base, not hardcoded path
source $(conda info --base)/etc/profile.d/conda.sh
conda activate fastp_env       # fastp, pandas, matplotlib, multiqc

# Personal conda (CRISPResso2, CrossMap, Liftoff, NextPolish)
CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate crispresso2_env # CRISPResso2, cas-offinder, pandas
conda activate liftoff_env     # Liftoff v1.5.1 (annotation transfer for scaffolded assembly)
conda activate crossmap_env    # CrossMap v0.7.3 (coordinate liftover using chain file)
conda activate nextpolish_env  # NextPolish v1.4.1 (short-read polishing of scaffolded assembly)
```

### CRISPResso2 Activation (use in all CRISPResso scripts)
```bash
CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate crispresso2_env
export PATH="${CONDA_BASE}/envs/crispresso2_env/bin:$PATH"
```

---

## Project Directory Structure

```
PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
│
├── samples.txt                    → 15 sample names (one per line)
├── CLAUDE.md                      → this file
│
├── reference/
│   ├── GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna  → reference genome
│   ├── GCF_000633615.1_*.fna.{amb,ann,bwt,pac,sa}       → BWA index
│   ├── GCF_000633615.1_*.fna.fai                         → samtools index
│   ├── GCF_000633615.1_*.dict                            → GATK dict
│   ├── intervals.list                                     → 24 chromosomes
│   ├── amplicon_bdnf_60bp_fwd.fa                         → amplicon FWD
│   ├── amplicon_bdnf_60bp_rc.fa                          → amplicon RC
│   ├── amplicon_bdnf_200bp_rc.fa                         → 200bp amplicon RC
│   ├── GCF_000633615.1_annotation.gff                    → Trinidad RefSeq annotation (NCBI)
│   ├── pseudogenome/                                      → bcftools-consensus pseudogenome
│   │   ├── colombian_pseudogenome.fna                     → pseudogenome sequence
│   │   ├── colombian_pseudogenome.fna.{fai,amb,ann,bwt,pac,sa,mmi} → samtools + BWA + minimap2 index
│   │   ├── colombian_pseudogenome.dict                    → GATK sequence dictionary
│   │   ├── colombian_pseudogenome.chain                   → coordinate liftover chain (bcftools)
│   │   ├── colombian_pseudogenome.liftoff.gff3            → ✅ PRIMARY annotation (Liftoff, Aug 6 2026)
│   │   ├── liftoff_unmapped_genes.txt                     → 4 genes not transferred by Liftoff (tRNAs + LOC103465184)
│   │   ├── README.md                                      → build method, QC, limitations, quick-start
│   │   ├── blast_db/                                      → BLAST nucleotide database
│   │   └── crossmap_superseded/                           → superseded CrossMap output (64% broken hierarchy)
│   └── colombian_scaffolded_genome/                       → SPAdes+RagTag+Liftoff assembly
│       ├── colombian_scaffolded.fna                       → scaffolded genome sequence
│       ├── colombian_scaffolded.fna.{fai,amb,ann,bwt,pac,sa} → samtools + BWA index
│       ├── colombian_scaffolded.dict                      → GATK sequence dictionary
│       ├── colombian_scaffolded.liftoff.gff3              → annotation (Liftoff)
│       ├── README.md                                      → build method, QC, limitations, quick-start
│       └── blast_db/                                      → BLAST nucleotide database
│
├── raw_fastq/                     → original *_R1_001.fastq.gz files
│
├── trimmed_fastp/                 → fastp filtered reads
│   └── fastqc_results/            → FastQC results
│
├── trimmed_trimmomatic/           → trimmomatic paired reads
│   └── fastqc_results/            → FastQC results
│
├── mapping/
│   ├── fastp/                     → sorted BAMs from fastp
│   │   └── fastqc_results/
│   └── trimmomatic/               → sorted BAMs from trimmomatic
│       └── merged/                → merged BAMs by group
│           ├── Control_merged.sorted.bam
│           ├── RNP_Cas_merged.sorted.bam
│           ├── Plasmid_Ko_merged.sorted.bam
│           └── Only_MNP_merged.sorted.bam
│
├── gatk/
│   ├── fastp/
│   │   ├── markdup/               → markdup BAMs + .bai
│   │   ├── metrics/               → duplication metrics
│   │   ├── gvcf/                  → per-sample GVCFs (.g.vcf.gz + .tbi)
│   │   ├── genomicsdb/            → GenomicsDB workspace
│   │   ├── vcf/                   → joint genotyped VCF
│   │   └── vcf_filtered/          → filtered SNPs and INDELs
│   └── trimmomatic/               → same structure as fastp
│
├── crispresso/
│   ├── ontarget/
│   │   └── trimmomatic/           → individual sample results (per markdup BAM)
│   │       └── merged/            → merged by group results (per merged BAM)
│   ├── pooled/
│   │   └── trimmomatic/           → CRISPRessoPooled results by group
│   ├── compare/
│   │   └── trimmomatic/           → CRISPRessoCompare results (pooled groups)
│   │       └── merged/            → CRISPRessoCompare results (merged groups)
│   ├── aggregate/                 → CRISPRessoAggregate results
│   ├── wgs/
│   │   └── trimmomatic/           → CRISPRessoWGS per sample
│   │       └── aggregate/         → CRISPRessoWGS aggregate by group
│   └── offtargets/
│       ├── offtargets.txt         → Cas-OFFinder raw results
│       └── combined/              → merged Cas-OFFinder + CRISPOR
│           ├── combined_offtargets.csv
│           ├── offtargets_crispresso_wgs.bed
│           ├── combined_offtargets_igv.bed
│           └── offtargets_intervals.list
│
├── coverage/
│   ├── bdnf_site_v2/              → raw depth/coverage text files (per sample)
│   └── csv/                       → parsed CSVs + Python plot scripts
│       ├── coverage_summary_all_samples.csv
│       ├── depth_per_position_all_samples.csv
│       ├── depth_by_zone.csv
│       ├── coverage_by_zone_all_samples.csv
│       ├── plot_coverage.py
│       ├── plot_depth_by_zone.py
│       ├── plot_depth_by_position.py
│       ├── plot_coverage_by_zone.py
│       ├── coverage_plots/         → PNGs from plot_coverage.py
│       ├── depth_zone_plots/       → PNGs from plot_depth_by_zone.py
│       ├── depth_position_plots/   → PNGs from plot_depth_by_position.py
│       └── coverage_zone_plots/    → PNGs from plot_coverage_by_zone.py
│
├── multiqc_report/
│   └── multiqc_trimming_comparison.html  → MultiQC comparison report
│
├── data/
│   ├── crispor_offtargets.xls     → raw CRISPOR output (guide ID: 326forw)
│   └── offtargets/
│       └── crispor/               → CRISPOR converted outputs
│           └── crispor_326forw_casoffinder_format.txt
│
├── assembly/                      → de novo assembly pipeline outputs
│   ├── spades_control_coassembly/ → SPAdes contigs.fasta (3 Control replicates)
│   ├── ragtag_output/             → RagTag scaffolded genome
│   └── qc_results/
│       ├── quast/                 → QUAST report vs reference
│       ├── busco/                 → BUSCO completeness results
│       └── busco_downloads/       → pre-staged actinopterygii_odb12.2 lineage
│
├── igv_files/                     → IGV-ready files for visualization
│   ├── features_of_interest.bed  → CRISPR features BED (bdnf, sgRNA, cut site, 8 off-targets)
│   ├── colombian_pseudogenome.*  → FASTA symlinks + sorted/indexed Liftoff GFF3 (.gz + .tbi)
│   └── *_merged.sorted.bam(.bai) → BAM symlinks for 4 groups (large, download separately)
│
└── codes/                         → all pipeline scripts
    ├── filtering/                 → trimming + FastQC scripts (6 scripts)
    ├── mapping/                   → BWA scripts (5 scripts)
    ├── variant_calling/           → GATK pipeline scripts (8 scripts)
    ├── CRISPResso/                → CRISPResso2 + off-target scripts (21 scripts)
    ├── analysis/                  → coverage computation scripts (2 scripts)
    └── assembly/                  → de novo + pseudogenome pipeline scripts (14 scripts)
```

---

## Pipeline Status

| Step | Tool | Status | Notes |
|---|---|---|---|
| Trimming | fastp | ✅ Done | Both fastp and trimmomatic completed |
| Trimming | trimmomatic | ✅ Done | Using NexteraPE-PE.fa adapters |
| QC | FastQC + MultiQC | ✅ Done | Poly-G artifacts from NovaSeq X confirmed normal |
| Mapping | BWA MEM | ✅ Done | bwa/0.7.17, both fastp and trimmomatic |
| Coverage | samtools | ✅ Done | ±500bp around sgRNA site, by zone |
| Mark duplicates | GATK | ✅ Done | --REMOVE_DUPLICATES false |
| HaplotypeCaller | GATK | ✅ Done | 15 samples, job 627832+669462, medium partition 7-day limit. ~3.9h/chrom × 24 chroms. |
| GVCF reindex | tabix | ✅ Done | All 15 GVCFs bgzip+tbi, completed Jul 12 2026 |
| GenomicsDBImport | GATK | ✅ Done | job 669470, completed Jul 12 23:47 |
| GenotypeGVCFs | GATK | ✅ Done | job 669471, all_samples.vcf.gz 3.5G, all 24 chroms, completed Jul 13 04:47 |
| VariantFiltration | GATK | ✅ Done | job 669472, snps_filtered + indels_filtered VCFs, completed Jul 13 05:03 |
| SelectVariants (off-target sites) | GATK | ✅ Done | job 669473, offtarget_variants.vcf.gz + offtarget_indels.vcf.gz, completed Jul 13 04:47 |
| Hotspot analysis | bedtools + Python | ✅ Done | job 683993-684016, 403 merged hotspot regions, 1780 windows FDR<0.05, Jul 14 2026 |
| Colombian pseudogenome | bcftools consensus | ✅ Done | job 683995, Control AF≥0.667 SNPs+INDELs applied, BWA-indexed, Jul 14 2026 |
| Pseudogenome annotation | CrossMap v0.7.3 | ⚠️ Superseded | job 684246, 95% transfer, but 64% of genes had broken GFF3 hierarchy — replaced by Liftoff |
| Pseudogenome annotation | Liftoff v1.5.1 | ✅ Done | job 692709, 99.5% transfer, 26,264 genes, 0 orphaned records, bdnf ✅, Aug 6 2026 |
| Pseudogenome verification | verify_pseudogenome.sh | ✅ Done | 12/12 checks passed (job 685661), Jul 15 2026 |
| Scaffolded genome indexing | samtools+BWA+GATK+BLAST | ✅ Done | job 692704, all indices generated, Aug 6 2026 |
| Pseudogenome indexing | GATK CreateSequenceDictionary | ✅ Done | .dict generated Aug 6 2026 |
| Merge BAMs | samtools | ✅ Done | markdup BAMs merged by group (4 merged BAMs) |
| CRISPResso2 on-target (individual) | CRISPResso2 | ✅ Done | 15 samples completed |
| CRISPResso2 on-target (merged) | CRISPResso2 | ✅ Done | 4 groups completed May 18 |
| CRISPRessoPooled | CRISPResso2 | ⚠️ Superseded | Never worked (invalid flag); would exactly duplicate CRISPResso2 on-target (merged) results — see Known Issues |
| CRISPRessoCompare (pooled) | CRISPResso2 | ⚠️ Superseded | CLAUDE.md previously said "done" — WRONG, its input (pooled/) was always empty. Superseded along with CRISPRessoPooled |
| CRISPRessoCompare (merged) | CRISPResso2 | ✅ Done | 6 pairwise merged group comparisons |
| CRISPRessoAggregate (on-target, all samples) | CRISPResso2 | ✅ Done | job 692758, 15 samples, Aug 6 2026 — see Known Issues for -p prefix + report-crash fixes |
| Off-target prediction | Cas-OFFinder | ✅ Done | 4 mismatches, 9 sites (1 on-target) |
| Off-target prediction | CRISPOR | ✅ Done | 8 sites, 4 mismatches, guide ID: 326forw |
| Combine off-targets | Python | ✅ Done | 8 unique off-target sites |
| CRISPRessoWGS (individual) | CRISPResso2 | ✅ Done | 15 samples × 8 sites completed May 13 |
| CRISPRessoWGS aggregate | CRISPResso2 | ✅ Done | job 692749, 4 groups, Aug 6 2026 — see Known Issues for -p prefix + report-crash fixes |
| De novo co-assembly | SPAdes v4.0.0 | ✅ Done | 3,379,459 contigs, N50 669bp — highly fragmented, see Known Issues |
| Reference-guided scaffolding | RagTag v2.1.0 | ✅ Done | N50 28.4Mb, 692Mb total, 24 chromosome-scale scaffolds + Chr0 |
| Assembly QC | QUAST v5.0.2 + BUSCO v5.7.1 | ✅ Done | Genome fraction 82.8%, BUSCO C:87.1% [S:86.2%,D:0.9%] — see Known Issues |
| Annotation transfer | Liftoff | ✅ Done | bdnf transferred at 94.5% coverage / 92.3% identity |

---

## Key Parameters and Decisions

### Trimming
```
fastp:       --qualified_quality_phred 12 --length_required 50
             --adapter_sequence CTGTCTCTTATACACATCT
             --adapter_sequence_r2 ATGTGTATAAGAGACA
trimmomatic: ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:8:keepBothReads
             LEADING:2 TRAILING:2 SLIDINGWINDOW:4:12 MINLEN:50
Note:        NovaSeq X uses binned quality scores (Q2/Q12/Q23/Q37)
             Aggressive quality trimming is not appropriate
Unpaired reads: EXCLUDED from all downstream analysis (mapping, GATK, CRISPResso2)
             Reason: lower quality (partner discarded), minimal coverage gain,
             risk of false positive CRISPR signal. Only paired R1+R2 used for BWA.
Adapter file path (hypatia cluster):
             /hpcfs/apps/conda4.12.0/envs/trimmomatic-0.39/share/trimmomatic/adapters/NexteraPE-PE.fa
```

### Mapping
```
Tool:    bwa mem (bwa/0.7.17)
Command: bwa mem -t THREADS -R "@RG\tID:SAMPLE\tSM:SAMPLE\tPL:ILLUMINA\tLB:lib1\tPU:unit1"
         Output piped directly to samtools sort (no intermediate SAM file)
Index:   samtools faidx + gatk CreateSequenceDictionary required for GATK
Merge:   samtools merge of markdup BAMs (gatk/trimmomatic/markdup/*.markdup.bam)
         → mapping/trimmomatic/merged/${GROUP}_merged.sorted.bam
         Merged BAMs are sourced from POST-GATK markdup BAMs, not from raw sorted BAMs
FastQC:  90 files processed (60 trimmomatic paired + 30 fastp filtered)
         File list generated by generate_fastqc_files.sh → fastqc_trimmed_filelist.txt
         Flag: --noextract (reports kept zipped)
```

### GATK
```
MarkDuplicates: --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY SILENT
HaplotypeCaller: -ERC GVCF -ploidy 2 -L intervals.list --native-pair-hmm-threads 4
  Runtime: ~3.9h/chromosome × 24 chromosomes = ~93h total per sample
  Partition: MUST use medium (7-day limit). short (2 days) and 3-day limits both failed.
  Array: --array=1-15%8 (8 concurrent × 4 CPUs = 32 CPUs; medium MaxCPUPU=48)
  Note: progress goes to .err log, not .out — check .err to see current chromosome
GVCF format: must be bgzip compressed (not regular gzip) for tabix indexing
  Truncated GVCFs (from time-limit cancellation) cause "Invalid GZIP header" in GenomicsDB
GenomicsDBImport: uses intervals.list file (NOT inline intervals)
  Cannot write to existing directory — always rm -rf genomicsdb/ before re-running
VariantFiltration: hard filtering only (no VQSR - no known variant DB for guppy)
SNP filters: QD<2, FS>60, MQ<40, MQRankSum<-12.5, ReadPosRankSum<-8
INDEL filters: QD<2, FS>200, ReadPosRankSum<-20
TMPDIR: /tmp/${USER}_${SLURM_JOB_ID} (cleanup with trap)
```

### CRISPResso2
```
sgRNA:  TGAGAGACGCCCCGGGCATG (on negative strand)
BAM mode: --bam_input + --bam_chr_loc (no -r1/-r2 needed)
Filtering: samtools view -b -F 0x400 -h BAM REGION -o filtered.bam
           (F 0x400 excludes duplicates; no -r flag — region as positional arg)
Sort+index: required before passing to CRISPResso2 (both samtools sort AND index)
Invalid params: --genome and --exclude_duplicates do NOT exist in --bam_input mode
           (genome is implicit in BAM; filter duplicates with samtools -F 0x400 first)
```

On-target analyses, by BAM source (region/amplicon corrected 2026-08-06 —
previously documented as 60bp individual vs 200bp merged, which was WRONG;
both scripts actually use the same 100bp amplicon and same region):

| Analysis | Script | BAM source | Region filtered | Amplicon |
|---|---|---|---|---|
| Individual (15 samples) | crispresso_ontarget.sh | gatk/trimmomatic/markdup/*.markdup.bam | NC_024333.1:15922000-15922100 (100bp) | amplicon_bdnf_100bp(_rc).fa |
| Merged groups (4 groups) | crispresso_ontarget_merged.sh | mapping/trimmomatic/merged/*.sorted.bam | NC_024333.1:15922000-15922100 (100bp) | amplicon_bdnf_100bp(_rc).fa |
| Off-target WGS (15 samples) | crispresso_wgs.sh | gatk/trimmomatic/markdup/*.markdup.bam | from offtargets_crispresso_wgs.bed | — (WGS mode) |

Because Individual and Merged use identical region/amplicon parameters, a
group-level "pooled" analysis (combining the 3-4 individual-track samples
per group) would be numerically identical to the already-completed Merged
track — this is why CRISPRessoPooled was marked Superseded rather than
fixed, see "CRISPRessoPooled — Wrong Tool + Redundant with Merged Track"
in Known Issues.

```
Amplicon: reference/amplicon_bdnf_100bp.fa + amplicon_bdnf_100bp_rc.fa
          used for BOTH individual and merged on-target analysis
  (60bp/200bp amplicon files also exist in reference/ from earlier
  iterations but are not what the current scripts use)

CRISPRessoWGS key parameters:
  --min_reads_to_use_region 10
  --min_frequency_alleles_around_cut_to_plot 0.05
  --expand_ambiguous_alignments
  --skip_failed

CRISPResso on-target key parameters:
  --min_frequency_alleles_around_cut_to_plot 0.05
  --expand_ambiguous_alignments
  --place_report_in_output_folder
  --write_cleaned_report
```

### Assembly Pipeline (SPAdes / RagTag / QUAST / BUSCO / Liftoff)
```
Goal: population-specific Colombian guppy genome via TRUE de novo co-assembly
      (recovers structural variants / novel sequence), as opposed to the
      existing reference/pseudogenome/ (bcftools consensus - majority-allele
      SNP/indel substitution onto a single linear copy of the reference, no
      SV recovery). Both genomes are complementary, not redundant.

Modules: spades/4.0.0, ragtag/2.1.0, quast/5.0.2, busco/5.7.1
Liftoff: not a module — installed via conda in a dedicated env:
  CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
  source ${CONDA_BASE}/etc/profile.d/conda.sh
  conda activate liftoff_env
  (created with: mamba create -n liftoff_env -c bioconda -c conda-forge
   --channel-priority flexible liftoff — plain `mamba create ... liftoff`
   without --channel-priority flexible fails to solve, biopython version
   conflict between bioconda/conda-forge)

Phase 1 - SPAdes co-assembly (codes/assembly/spades_coassembly.sh):
  Input: 3 Control replicates, trimmomatic-paired reads (--pe1..--pe3)
  Input size: ~47GB compressed (confirmed via du) — sized generously
  Resources: --mem=250G --time=10-00:00:00 --partition=bigmem
    (draft's original 64G/24h/short was undersized given the confirmed
    input size — same class of problem as the HaplotypeCaller history above)
  --checkpoints last enables resume via `spades.py --continue -o <dir>`
    if the job times out, instead of restarting from scratch
  --isolate required for this combined depth (~160-270x) — see "SPAdes —
    mmap ENOMEM Without --isolate" in Known Issues below
  Output: assembly/spades_control_coassembly/contigs.fasta

Phase 2 - RagTag scaffolding (codes/assembly/ragtag_scaffold.sh):
  Input MUST be filtered to >=500bp first (contigs.min500.fasta, not the
  raw 3.4M-contig contigs.fasta) — see "RagTag Hang on 3.4M Unfiltered
  Contigs" in Known Issues.
  ragtag.py scaffold REF CONTIGS.min500 -o assembly/ragtag_output -t 8 -C
  Output: assembly/ragtag_output/ragtag.scaffold.fasta
  Result: N50 28.4Mb, total 692Mb, 24 chromosome-scale scaffolds + ~2039
    small unplaced fragments + 1 Chr0 (39Mb, concatenated unplaced via -C)

Phase 3 - QUAST + BUSCO QC (codes/assembly/quast_qc.sh, busco_qc.sh):
  Split into two independent parallel jobs (both depend only on RagTag) —
  see "QUAST/BUSCO — Slow Contig Analyzer on Fragmented Ref+Scaffold" in
  Known Issues for why.
  QUAST: quast.py SCAFFOLD -r REF -o assembly/qc_results/quast --fragmented
    (--fragmented because both REF and the RagTag scaffold are fragmented)
    Input is ragtag.scaffold.noChr0.fasta (Chr0 excluded) with a 48h limit —
    see "QUAST/BUSCO — Slow Contig Analyzer" Known Issue for why (2 failed
    attempts first: an 8h combined job, then a 24h --fragmented-only job).
    Result (2026-07-15): genome fraction 82.807%, N50 28.3Mb (Chr0-excluded
    set), duplication ratio 1.078, 7589 misassemblies.
  BUSCO: actinopterygii_odb10, NOT odb12.2 — BUSCO 5.7.1's own code
    hardcodes datasets_version=="odb10" as a fatal requirement, no flag
    bypasses it (odb12.2 is what `busco --list-datasets` shows as current,
    but this BUSCO version literally cannot run with it). See "BUSCO —
    Hardcoded odb10 Version Check" in Known Issues.
    wget https://busco-data.ezlab.org/v5/data/lineages/actinopterygii_odb10.2024-01-08.tar.gz
    extracted to assembly/qc_results/busco_downloads/lineages/actinopterygii_odb10/
    then run with --offline --download_path assembly/qc_results/busco_downloads
    Result (2026-07-13): C:87.1% [S:86.2%,D:0.9%], F:6.2%, M:6.7%,
    n:3640 (3169 complete BUSCOs)
  See codes/assembly/00_download_busco_lineage.sh for the workaround.

Phase 3b - Liftoff annotation for pseudogenome (codes/assembly/liftoff_pseudogenome.sh):
  Annotates the bcftools-consensus pseudogenome using Liftoff v1.5.1 (liftoff_env).
  CrossMap was attempted first (crossmap_pseudogenome.sh, job 684246) but produced
  broken GFF3 hierarchy: 16,637 of 26,000 genes had gene/mRNA records dropped to
  .unmap while exons transferred correctly — making 64% of genes unusable by tools
  that require parent-child hierarchy. Root cause: CrossMap cannot remap large gene
  spans (>10 kb) when multiple INDELs create coordinate ambiguity at the boundaries.
  Liftoff uses minimap2 sequence alignment per gene, handles the hierarchy correctly
  by design, and never produces orphaned records.
  liftoff -g TRINIDAD_GFF -o OUT_GFF -p 8 PSEUDO_FASTA REF_FASTA
  Transfer rate: 99.5% (1,287,616 / 1,293,974 features), 26,264 genes, 0 orphaned.
  bdnf: NC_024333.1:15923726-15938393 (- strand), 8 mRNA isoforms, 22 exons.
  PRIMARY annotation: reference/pseudogenome/colombian_pseudogenome.liftoff.gff3
  CrossMap output retained as reference: colombian_pseudogenome.gff3 (superseded).
  Verification: codes/assembly/verify_pseudogenome.sh — 12/12 checks passed (Jul 15).
  IGV annotation: igv_files/colombian_pseudogenome.gff3.gz rebuilt from Liftoff output.

Phase 4 - Liftoff annotation transfer (codes/assembly/liftoff_annotation.sh):
  Transfers Trinidad gene models (incl. bdnf) onto the new scaffold coords.
  Requires reference/GCF_000633615.1_annotation.gff (downloaded from NCBI —
  did not exist anywhere in the project before this pipeline; see
  codes/assembly/00_download_gff.sh, md5-checksum verified against NCBI).
  liftoff -g TRINIDAD_GFF -o NEW_GFF -p 8 NEW_FASTA REF
  Verification: grep "ID=gene-bdnf" on the output GFF3 — bdnf feature is at
  NC_024333.1:15920888-15935548 in the Trinidad GFF, must reappear (same ID)
  at new coordinates after liftover.
  Result (2026-07-11): bdnf lifted to NC_024333.1_RagTag:14113623-14129072,
  coverage=0.945, sequence_ID=0.923 — high-confidence transfer.

Phase 5 - NextPolish short-read polishing (codes/assembly/nextpolish_genome.sh):
  Corrects small-scale errors (indels/SNPs) by realigning the same Illumina
  short reads used to build the genome back onto colombian_scaffolded.fna.
  Targets the ~6.9% of "complete" BUSCO genes with internal stop codons
  (likely frameshift artifacts) found during QC — see README.md in
  reference/colombian_scaffolded_genome/ for the pre-polish quality baseline.
  NextPolish v1.4.1 (nextpolish_env, installed via
  codes/assembly/00_setup_nextpolish_env.sh — mamba create ... nextpolish,
  bundles bwa/samtools/minimap2 as dependencies).
  Config format (run.cfg, generated dynamically by nextpolish_genome.sh, not
  a static file): required keys are `genome` and one of
  `sgs_fofn`/`lgs_fofn`/`hifi_fofn`. `sgs_fofn` is a plain list of fastq
  file paths, ONE PER LINE (not space-separated R1/R2 pairs) — NextPolish
  pairs them automatically. `task = best` with only sgs_fofn provided
  auto-resolves to `[1, 2, 1, 2]` (2 rounds of short-read polishing) — this
  is NextPolish's own recommended recipe for Illumina-only data, confirmed
  by reading config_parser.py directly (the C-compiled binary doesn't
  expose this in --help). `polish_options = -p {multithread_jobs}` uses
  literal template substitution — do not hardcode a thread count there.
  Command: nextPolish run.cfg
  Output: ${workdir}/genome.nextpolish.fasta (+ .stat)
  Resource sizing rationale: the project's own single-sample BWA alignment
  against a similarly-sized reference (bwa_trimmomatic_array.sh) actually
  took only ~2-3h at 8 threads in practice (per sacct, job 437545 array —
  the requested 20h was a safety margin, not the real runtime). NextPolish
  aligns all 3 Control samples (parallel_jobs=4, so effectively
  concurrent), twice. Sized generously anyway (32 CPUs, 64GB, 3 days,
  medium partition) per this project's established pattern of avoiding
  first-attempt under-provisioning.
  Status: job 692710 failed after ~1h (N-content rejection, see "NextPolish
  — N-content Rejection" Known Issue) — fixed with -N in sgs_options.
  job 692760 (resubmit) failed after ~2h44m at the polish_genome step
  (forkserver/fork multiprocessing bug in NextPolish 1.4.1 itself, see
  "NextPolish — forkserver vs fork Multiprocessing" Known Issue) — fixed by
  patching nextpolish1.py/nextpolish2.py in the conda env directly, fix
  verified against the exact failed command before resubmitting. Third
  submission (job 699678) COMPLETED 2026-08-10, ~9h, produced
  genome.nextpolish.fasta (2064 seqs, 691,921,577bp — Chr0 renamed
  Chr0_RagTag_np1212, NextPolish appends the task sequence to every
  sequence name).

  QC comparison (job 700892 QUAST, job 700893 BUSCO, both 2026-08-11) —
  RESULT: NO MEANINGFUL IMPROVEMENT, several structural metrics slightly
  WORSE. Polished genome NOT adopted; colombian_scaffolded.fna remains
  authoritative. See "NextPolish — No Improvement" Known Issue below and
  reference/colombian_scaffolded_genome/README.md's "Polishing experiment"
  section for the full comparison table and explanation (short-read
  polishing has a low ceiling when using the exact same reads that built
  the assembly — nothing to correct where the assembly already matches
  what those reads say, even if it differs from the BUSCO ortholog set).

Driver: codes/assembly/run_assembly_pipeline.sh
  Chain: SPAdes → RagTag → {QUAST, BUSCO, Liftoff} (all three fan out from
  RagTag independently). Driver does NOT auto-filter contigs to >=500bp
  between SPAdes and RagTag — that must be done manually (see driver's own
  comment header for the exact command).
  Prerequisites (run once on login node before the driver):
    00_download_gff.sh, 00_download_busco_lineage.sh, 00_setup_liftoff_env.sh
```

### Coverage Analysis
```
Region: NC_024333.1:15921439-15922658 (±600bp around sgRNA)
Zones:
  upstream_500bp:   15921439-15921938 (500bp)
  upstream_100bp:   15921939-15922038 (100bp)
  sgRNA_site:       15922039-15922058 (20bp)
  downstream_100bp: 15922059-15922158 (100bp)
  downstream_500bp: 15922159-15922658 (500bp)
Note: Coverage_Pct = 100% expected for WGS with ~50x depth
      Use Mean_Depth for meaningful comparisons between zones
```

---

## Off-Target Sites for WGS Analysis

| Chromosome | Position | Strand | MM | MIT Score | CFD Score | Source | Locus |
|---|---|---|---|---|---|---|---|
| NC_024331.1 | 5708724 | + | 3 | 0.312 | 0.071 | CRISPOR | exon:XM_008429791.2 |
| NC_024331.1 | 13951199 | - | 4 | 0.092 | 0.141 | CRISPOR | exon:XM_017306200.1 |
| NC_024331.1 | 26228796 | + | 4 | 0.020 | 0.144 | CRISPOR | intergenic |
| NC_024332.1 | 5810651 | + | 4 | 0.113 | 0.012 | CRISPOR | intron:XM_008425835.2 |
| NC_024338.1 | 20512932 | - | 4 | 0.017 | 0.109 | CRISPOR | intergenic |
| NC_024339.1 | 7034820 | - | 4 | 0.118 | 0.065 | CRISPOR | intron:XM_017306598.1 |
| NC_024340.1 | 12200655 | - | 4 | 0.200 | 0.044 | CRISPOR | intergenic |
| NC_024349.1 | 24882037 | + | 4 | 0.199 | 0.031 | CRISPOR | intron:XM_008437971.2 |

BED file: `crispresso/offtargets/combined/offtargets_crispresso_wgs.bed`

### Off-Target Discovery Workflow

**Cas-OFFinder:**
```
Input:    genome FASTA path list + PAM (NNNNNNNNNNNNNNNNNNNNNGG) + sgRNA + max mismatches
Command:  cas-offinder input.txt C offtargets.txt
MM used:  4 (returns all sites with 0–4 mismatches)
Output:   tab-separated (guideSeq, chrom, pos, sequence, strand, MM_count)
Result:   9 sites total (1 on-target NC_024333.1:15922039 + 8 off-targets)
```

**CRISPOR (web tool, converted output):**
```
Input:    data/crispor_offtargets.xls (raw CRISPOR Excel export)
Guide ID: "326forw" (filter applied in convert_crispor_offtargets.py)
Script:   codes/CRISPResso/convert_crispor_offtargets.py
Result:   8 off-target sites with MIT score, CFD score, locus annotation
Output:   data/offtargets/crispor/crispor_326forw_casoffinder_format.txt
          (Cas-OFFinder-compatible format, used as input to combine_offtargets.py)
```

**Deduplication (±2bp tolerance):**
```
All 8 CRISPOR sites found by CasOFFinder within ±2bp
Priority: CRISPOR entry kept (has scores + annotation) — see deduplication Known Issue above
Final:    8 unique off-target sites in combined_offtargets.csv
```

---

## Common SLURM Script Template

```bash
#!/bin/bash
#SBATCH --job-name=job_name
#SBATCH --array=1-15%8          # 15 samples, max 8 concurrent
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --partition=short
#SBATCH --mail-user=diegoandres3322@gmail.com
#SBATCH --mail-type=ALL

PROJECT_DIR=/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data
SAMPLE_LIST=${PROJECT_DIR}/samples.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

# TMPDIR setup (always include for GATK jobs)
TMPDIR="/tmp/${USER}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
export TMPDIR
trap "rm -rf $TMPDIR" EXIT
```

---

## Known Issues and Solutions

### samtools view — Region Specification
```bash
# CORRECT: region as positional argument at end, no -r flag
samtools view -b -F 0x400 -h "$BAM" "NC_024333.1:15922011-15922071" -o out.bam

# WRONG: -r flag does not exist in samtools view
samtools view -r "NC_024333.1:..." "$BAM"  # ❌
```

### Inline Comments After Backslash
```bash
# WRONG: comment after \ breaks line continuation
samtools view -b \
  -F 0x400 \          # this breaks the command ❌

# CORRECT: no inline comments after \
samtools view -b \
  -F 0x400 \
  "$BAM"              # comment on separate line is ok
```

### GenomicsDB Workspace
```bash
# GenomicsDB cannot write to existing directory
# Always remove before re-running:
rm -rf ${PROJECT_DIR}/gatk/trimmomatic/genomicsdb/
# Do NOT mkdir after rm — GATK creates it internally
```

### GVCF Indexing
```bash
# GVCFs must be bgzip compressed (not gzip) for tabix
# Re-compress if needed:
zcat sample.g.vcf.gz | bgzip -c > sample.g.vcf.gz.tmp
mv sample.g.vcf.gz.tmp sample.g.vcf.gz
tabix -p vcf sample.g.vcf.gz
```

### CRISPResso2 — No Alignments Found
```bash
# Cause: amplicon too large (>290bp for 150bp reads)
# Rule: max amplicon = (read_length × 2) - 10 = 290bp
# Solution: use 60bp amplicon centered on cut site
# Also: provide both FWD and RC strands as amplicon
```

### CRISPResso2 Activation in SLURM
```bash
# mamba shell hook fails in SLURM — use conda directly
CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate crispresso2_env
export PATH="${CONDA_BASE}/envs/crispresso2_env/bin:$PATH"
```

### CRISPResso2 — Invalid Parameters with `--bam_input`
```bash
# WRONG: --genome and --exclude_duplicates don't exist in --bam_input mode
CRISPResso2 --bam_input input.bam --genome ref.fa --exclude_duplicates ...

# CORRECT: genome is implicit in the BAM header; exclude duplicates via samtools
samtools view -b -F 0x400 -h "$BAM" "NC_024333.1:15921600-15922400" -o filtered.bam
samtools sort -o filtered.sorted.bam filtered.bam
samtools index filtered.sorted.bam
CRISPResso2 --bam_input filtered.sorted.bam --amplicon_seq "FWD,RC" ...
```

### CRISPResso2 — BAM Must Be Sorted AND Indexed
```bash
# Both steps are required before passing to CRISPResso2
samtools sort -o output.sorted.bam input.bam
samtools index output.sorted.bam
# Skipping sort or index causes CRISPResso2 to fail with an empty alignment error
```

### FastQC — GC Content Warnings (NovaSeq X Artifact)
```
Small peak at 95-100% GC in FastQC "Per Sequence GC Content" is EXPECTED
for NovaSeq X patterned flowcell data. Caused by:
  - Poly-G artifacts (NovaSeq X empty clusters report G)
  - Mitochondrial / rRNA sequences
  - Telomeric repeats
This is NOT contamination. Cross-reference "Overrepresented Sequences" module
to confirm. ~10 warning flags per 120 FASTQ files is normal for this platform.
```

### fastp — `--detect_adapter_for_pe` Version Compatibility
```bash
# --detect_adapter_for_pe requires fastp >= v0.20.0
# Error: "undefined option: --detect_adapter_for_pe"
# Fix: explicitly specify adapter sequences (preferred since adapters are known):
fastp --adapter_sequence CTGTCTCTTATACACATCT \
      --adapter_sequence_r2 ATGTGTATAAGAGACA ...
```

### HaplotypeCaller — Time Limit History (Critical)
```
Run 1 (May 2026): --array=1-8, --time=12:00:00, no -L intervals.list
  Result: only 8/15 samples, GVCFs truncated at LG4 (GATK processes chroms sequentially)
  Root cause of LG6–LG23 zero variants in downstream VCF: truncated GVCFs

Run 2 (Jun 2026): --array=1-15%8, --time=20:00:00, with -L intervals.list
  Result: all 15 samples started; cancelled at 20h (SLURM TIME LIMIT)
  GVCFs truncated around LG14 → "Invalid GZIP header" in GenomicsDB

Run 3 (Jul 2026): --array=1-15%8, --time=7-00:00:00, partition=medium
  Tasks 1-8 cancelled (could not extend running jobs with scontrol) and resubmitted as job 627832
  Tasks 9-15 extended to 7 days via: scontrol update JobId=JOBID TimeLimit=7-00:00:00
  Note: scontrol can extend PENDING jobs but returns "Access/permission denied" for RUNNING jobs
  Task 627832_2 (Control_MNP_II) running 3.5x slower than others (~10h/chrom vs ~3h/chrom)
  — may time out again; resubmit with: sbatch --array=2 haplotype_caller.sh when it cancels

Lesson: WGS HaplotypeCaller needs 7-day limit on medium partition. Never use short partition.
Diagnosis: tail -30 logs/haplotype_trimmomatic_JOBID_TASKID.err | grep "NC_0"
           Progress meter shows elapsed minutes and current chromosome position.
```

### SPAdes — mmap ENOMEM Without `--isolate` on High-Coverage Data
```
Job 646213 (spades_coassembly, --mem=250G, bigmem/nodei-4) crashed after 46min:
  "mmap(2) failed. Reason: Cannot allocate memory. Error code: 12"
  at the k-mer counting step, with MaxRSS only ~15GB (far under the 250GB
  requested — confirmed via sacct; node had 562GB free; ulimits checked via
  srun test on the same node were correct: -m 250GB, -v ~275GB).

Root cause: combined depth across the 3 Control replicates is ~160-270x
(47GB compressed reads vs a ~700Mb genome) — well beyond SPAdes's default
"standard" mode assumptions (tuned for single-cell/uneven-coverage data).
The default code path tried to build an oversized k-mer index for this much
redundant high-coverage data and hit a single oversized mmap request.

Fix: add --isolate (SPAdes's own params.txt output explicitly recommends
this for "high-coverage isolate and multi-cell data" — exactly this
dataset — but it was missed in the initial draft/implementation).
  spades.py ... --isolate --checkpoints last -o $OUTPUT_DIR

Lesson: always pass --isolate for standard (non-single-cell, non-metagenomic)
high-coverage WGS co-assembly with SPAdes. Diagnosis path: sacct MaxRSS vs
--mem requested (rules out true OOM), srun ulimit -a / vm.max_map_count on
the same node (rules out cgroup/ulimit misconfiguration) — isolates the
cause to SPAdes's own mode-dependent memory behavior, not the SLURM request.

Round 2 (job 652298, --isolate, -m 250): ran ~24h, reached the K77 Distance
Estimation stage, then hit a GENUINE OOM (mimalloc ENOMEM, peak RSS ~249GB
against the 250G ceiling — confirmed via sacct MaxRSS, not a config bug this
time). Fix: SPAdes supports --restart-from last, which resumes from the
most recent saved checkpoint WITH UPDATED OPTIONS (e.g. a higher -m),
avoiding a full restart:
  spades.py --restart-from last -m 450 -o $OUTPUT_DIR
This resumed from the "late_pair_info_count" checkpoint (saved 10 min
before the crash) and only re-ran Distance Estimation onward — recovered
in ~16h instead of repeating the full ~24h. See
codes/assembly/spades_coassembly_resume.sh (one-off recovery script; the
standing spades_coassembly.sh was updated to -m 450 / --mem=470G for any
future from-scratch run, with a 20G buffer between SLURM's cgroup limit and
SPAdes's own self-limit so SPAdes self-terminates gracefully before a hard
cgroup kill).

Result (job 669505, completed 2026-07-10): contigs.fasta produced, but
HIGHLY FRAGMENTED — 3,379,459 contigs, 1.12 Gbp total (~1.5x the ~700-750Mb
expected genome), median length 155bp, N50 669bp, only 213,030 contigs
(6.3%) ≥1000bp. Root cause is NOT primarily the 150bp Illumina read length
(a contributing but secondary factor) — it's pooling 3 distinct wild-caught
Colombian guppies (not clonal replicates) at ~160-270x combined depth: every
heterozygous SNP/indel between the 3 individuals' genomes creates a graph
"bubble" that gets output as a separate short contig instead of collapsing.
Decision (2026-07-10): proceed to Phase 2 (RagTag) rather than attempting a
hybrid long-read reassembly — see [[nanopore_epigenome_data]] memory for a
~3.2x-depth Nanopore dataset (pooled telencephalon, 10 fish, same
population) that was considered but deferred as a future post-hoc
scaffolding/gap-filling supplement rather than a full hybrid rerun, given
the shallow depth. Initially tried Phase 2 with the fragmented contigs
as-is (no length filter) — see next Known Issue, this had to be revisited
once RagTag itself couldn't handle 3.4M contigs.
```

### RagTag Hang on 3.4M Unfiltered Contigs
```
Job 676655 (ragtag_scaffold.sh, unfiltered contigs.fasta, 4h limit) hung
indefinitely: minimap2 alignment + reading + filtering + ordering all
finished in under 4 minutes, then it hung at "Writing scaffolds" for the
remaining ~3h56m until the 4h wall-clock killed it — ragtag.scaffold.agp
was never produced.

Root cause: RagTag has no option to filter/exclude short QUERY sequences.
-f/--remove-small only govern which ALIGNMENTS are trusted, not which input
contigs get processed - RagTag always attempts to place or individually
write out every query sequence given to it. Feeding it 3.4M contigs (85%
of them <500bp noise fragments — see previous Known Issue) made the output-
writing step for the unplaced majority impractically slow.

Fix: pre-filter the query fasta before invoking RagTag (RagTag itself has
no equivalent option, this must be done upstream):
  awk '/^>/{split($0,p,"_"); keep=(p[4]+0>=500)} keep' contigs.fasta \
    > contigs.min500.fasta
  (works because SPAdes headers embed length directly, e.g.
  NODE_1_length_145897_cov_994.365889 — p[4] is the length field)
This dropped 3.4M → 501,663 contigs (>=500bp), covering the vast majority
of real assembled sequence. Also added RagTag's -C flag (concatenate any
remaining unplaced contigs into one chr0 record instead of writing each
individually) as a second safeguard against the same class of hang even at
500K contigs. ragtag_scaffold.sh time limit raised 4h → 12h for margin.
Filtering was a live decision revisited mid-run, not planned upfront — see
"SPAdes — mmap ENOMEM" Known Issue above for the original 2026-07-10
rationale for NOT filtering, superseded once this hang occurred.

Lesson: population-level de novo co-assemblies of pooled wild individuals
can produce contig counts (millions) that downstream tools like RagTag were
never designed to handle directly — check a tool's expected input scale,
not just its correctness, when chaining it after a highly fragmented
assembly.
```

### QUAST/BUSCO — Slow Contig Analyzer on Fragmented Ref+Scaffold
```
Job 683049 (qc_quast_busco.sh, QUAST+BUSCO sequential in one 8h job) timed
out mid-QUAST — BUSCO never started. minimap2 alignment itself finished in
~7min, but QUAST's single-threaded Python misassembly classification was
still actively writing output (ragtag-scaffold.coords.filtered still being
appended) at the moment of the timeout kill, not actually hung. Root cause:
BOTH the reference (2768 fragments) and the RagTag scaffold (~7220 N's per
100kbp from join gaps) are fragmented, generating a very large number of
individual alignment segments for QUAST to classify one at a time — exactly
the condition QUAST's own log recommends `--fragmented` for.

Fix: split into two independent parallel jobs (quast_qc.sh, busco_qc.sh),
both depending only on RagTag, each with a generous 24h limit — so a slow
QUAST run can't block or share a time budget with BUSCO. Added
`--fragmented` to the QUAST call per its own recommendation.

Lesson: sequentially chaining two independently-slow QC tools inside one
SLURM job compounds their time budgets and lets an early one starve a later
one of runtime — split into parallel jobs when steps don't depend on each
other's output.

Round 2 (job 683419, --fragmented, 24h limit): NOT a hang — filtered
alignment output grew ~4x vs. round 1 (15MB → 60MB of
ragtag-scaffold.coords.filtered), confirming steady real progress — but
still didn't finish in 24h. Checked whether the artificial Chr0_RagTag
sequence (RagTag's -C flag concatenates ~2039 unplaced fragments into one
39Mb pseudo-chromosome with no biological continuity between components)
was the dominant cause: it accounted for 53,735 of 525,459 total alignment
records (~10%) — a real but not dominant contributor, since the bulk of
the volume comes from genuinely comparing 692Mb of assembled sequence
against the reference's 2768 fragments.

Fix: two changes together — (1) time limit raised 24h → 48h given the
confirmed slow-but-real progress rate, (2) QUAST input switched to a
Chr0-excluded copy of the scaffold (samtools faidx -r <keep_list.txt>,
naming convention ragtag.scaffold.noChr0.fasta), since "misassemblies"
inside an artificial concatenation of unrelated fragments aren't
scientifically meaningful to detect anyway. The full ragtag.scaffold.fasta
(with Chr0) is untouched and still used by BUSCO/Liftoff, which don't hit
this bottleneck (BUSCO's gene-completeness scan is unaffected by scaffold
boundaries — it completed in 3 minutes).

Lesson: when a tool is slow-but-progressing (not hung) on a specific input,
check whether a piece of that input is scientifically meaningless for the
analysis at hand (like an artificial concatenation) before just throwing
more walltime at it — removing dead weight can compound with a larger time
budget rather than substituting for one.
```

### BUSCO — Hardcoded odb10 Version Check
```
Even after the odb12.2 lineage was successfully manually staged (working
around the module's download URL bug — see "Assembly Pipeline" section
above), busco_qc.sh (jobs 683420 and 683422) failed instantly (~16s) with:
  "ERROR: BUSCO v5 only works with datasets from OrthoDB v10 (with the
  suffix '_odb10')"
This happened even after passing --datasets_version odb12.2 explicitly.

Root cause (confirmed by reading BUSCO 5.7.1 source,
busco/BuscoConfig.py:648-680, check_lineage_present()): the final
validation is UNCONDITIONAL —
  if datasets_version != "odb10": raise BatchFatalError(...)
--datasets_version only controls what gets appended to a bare lineage name
that lacks a "_odb..." suffix; once a version is resolved (from either the
lineage name or the flag), BUSCO 5.7.1's code refuses to proceed unless
that resolved version is literally "odb10" — no CLI flag bypasses this.
This BUSCO module version was released before the OrthoDB catalog moved to
v12, and its dataset-version validation was never updated.

Fix: use actinopterygii_odb10 instead of odb12.2. It no longer appears in
`busco --list-datasets` (interactive catalog shows only odb12.2 for this
clade now), but the archived tarball is still hosted and downloadable:
  https://busco-data.ezlab.org/v5/data/lineages/actinopterygii_odb10.2024-01-08.tar.gz
Staged the same way as odb12.2 was (manual wget + extract to
busco_downloads/lineages/, since `busco --download` has its own separate
URL-construction bug for versioned names — see 00_download_busco_lineage.sh,
now downloads odb10). busco_qc.sh updated to `-l actinopterygii_odb10` (no
--datasets_version override needed — odb10 is what the name resolves to by
default). This was tried BEFORE getting the mmap/hang issues below sorted,
so ended up being the third distinct BUSCO/QUAST-stage issue in this
pipeline, after the RagTag hang.

Lesson: a "current" dataset version from a tool's live catalog is not the
same as what an older pinned software version actually supports — check
the tool's own source/validation logic (not just the catalog) when a
dataset-version mismatch error appears, especially for slow-moving cluster
modules that lag behind fast-moving companion data catalogs.
```

### CrossMap — bdnf Gene-Level Record Dropped to UNMAP
```
CrossMap v0.7.3 (job 684246) transferred 95% of features from Trinidad GFF
to the pseudogenome using colombian_pseudogenome.chain. However, the bdnf
gene-level record (NC_024333.1:15920888-15935548, 14.6 kb span) and its
mRNA container records ended up in the .unmap file. CrossMap could not
confidently remap the full gene boundary because multiple INDELs within
that 14.6 kb span produced coordinate ambiguity at the gene level.

Importantly: the individual exon and CDS records for bdnf DID transfer
correctly — CrossMap handles short, specific intervals better than large
gene-spanning records.

Fix: reconstructed gene and mRNA records from the coordinates of the
transferred exons (min/max of each mRNA's exon set), then inserted them
into the GFF3 before the first bdnf exon record using a Python script.
Result: gene-bdnf at NC_024333.1:15923726-15938393 (- strand), 8 mRNA
isoforms restored, 21 exons intact.

Script used: inline Python in verify session (Jul 15 2026) — see git history
or reproduce from the .unmap file + exon coordinates in the fixed GFF3.
Backup of original CrossMap output: colombian_pseudogenome.gff3.bak
```

### Pseudogenome Verification — Spot Check False Negatives
```
verify_pseudogenome.sh section 5 (variant spot check) shows ❌ symbols
but reports 0 failures in the SUMMARY. This is expected behavior:
  - The spot check uses snps_filtered.vcf.gz (all 15 samples, PASS)
  - The pseudogenome only contains Control-specific variants at AF ≥ 0.667
  - Positions 76, 102, 275 on NC_024333.1 from the all-sample VCF were
    not necessarily in the Control-filtered subset that was applied
  - The ❌ lines are printed by a while-read loop that does NOT call the
    check() function, so they do not increment the FAIL counter

The 12 formal checks (using check() and check_range()) all passed.
The sgRNA site sequence difference (section 4) is also expected: INDELs
applied upstream shift coordinates so the same position in the pseudogenome
corresponds to a different genomic region than in the Trinidad reference.
```

### CRISPRessoAggregate — `-p` Semantics, Output Directory, and a Report Bug
```
crispresso_wgs_aggregate.sh (original) reported "Analysis Complete! 100%"
for all 4 groups but produced NOTHING useful — "Read 0 folders (0 not
imported)". Two separate bugs, found 2026-08-06:

1. -p is a PREFIX for glob-matching folder paths, not a parent directory to
   scan. The script passed `-p ${WGS_DIR}/${SAMPLE}/` (the sample's WGS
   output dir itself), but the actual per-site run folders CRISPRessoAggregate
   looks for live one level deeper: `${WGS_DIR}/${SAMPLE}/CRISPRessoWGS_on_${SAMPLE}/CRISPResso_on_<site>`.
   Fix: point -p at `.../CRISPRessoWGS_on_${SAMPLE}/CRISPResso_on_` (trailing
   partial-name prefix, glob-matches all per-site subfolders). Confirmed via
   CRISPRessoAggregate --help ("-p PREFIX ... Prefix for CRISPResso folders
   to aggregate") and by testing interactively before fixing the SLURM script.
   Same fix applied to the new crispresso_aggregate_ontarget.sh (Pending
   Analyses #1 Track C), except there each sample's on-target run is a single
   exact-match folder (no per-site glob needed): one -p per sample pointing
   directly at `crispresso/ontarget/trimmomatic/<SAMPLE>/CRISPResso_on_<SAMPLE>`.

2. CRISPRessoAggregate has NO output-directory flag - it always creates its
   `CRISPRessoAggregate_on_<name>` folder in the current working directory.
   The original script defined OUTPUT_DIR but never `cd`ed into it, so the
   (empty, broken) output landed in codes/CRISPResso/ instead. Fix: `cd
   "$OUTPUT_DIR"` before invoking CRISPRessoAggregate.

3. Separately, this CRISPResso2 install (2.3.1) has a real bug:
   `make_multi_report() missing 2 required positional arguments:
   'crispresso_tool' and 'logger'` - crashes ONLY the final HTML report
   step. All the actual data (txt tables, pdf/png plots) generate correctly
   before that point. Fix: add `--suppress_report`. A second, non-fatal
   bug from the same version also appears in the logs and can be ignored:
   `Error in plot pool: plot_nucleotide_quilt() missing 1 required
   positional argument: 'custom_colors'` (breaks one specific plot type,
   doesn't block completion or the rest of the output).

Lesson: "Analysis Complete!" / exit code 0 does not mean useful output was
produced for this tool - always check "Read N folders (M not imported)" in
the log, and spot-check that expected output files actually contain data,
not just that the job didn't crash.
```

### CRISPRessoPooled — Wrong Tool + Redundant with Merged Track
```
crispresso_pooled_groups.sh (Pending Analyses #1 Track A) never worked:
`CRISPRessoPooled: error: unrecognized arguments: --crispresso_output_folders
...` for all 4 groups (the script had no error checking, so it printed "✅
done" after each failed call anyway - a second, independent bug). Also used
the `mamba shell hook` activation pattern the "CRISPResso2 Activation in
SLURM" Known Issue already warns against (fixed in passing, but moot given
the finding below).

Root cause: CRISPRessoPooled's actual purpose is to take RAW reads (-r1/-r2)
plus multiple amplicons and do its OWN demultiplexing + per-amplicon
analysis in one step ("pooled amplicon sequencing"). It has no flag to
combine pre-existing CRISPResso_on_X output folders - that is
CRISPRessoAggregate's job, not CRISPRessoPooled's. `--crispresso_output_folders`
is not a real flag in CRISPResso2 v2.3.1 (confirmed via --help).

The correct way to build a genuine group-level "pooled" analysis (not just
a summary) would be to merge BAMs by group and run CRISPResso directly via
--bam_input, exactly mirroring crispresso_ontarget_merged.sh (Track B).
But: reading both scripts confirmed Track A (individual, crispresso_ontarget.sh)
and Track B (merged) use IDENTICAL region and amplicon parameters (see
"Assembly Pipeline" — actually "CRISPResso2" section above, corrected
2026-08-06 from a previously-wrong 60bp-vs-200bp distinction). Since Track
B already merges the same BAMs with the same parameters, a Track A "pooled"
step would exactly reproduce Track B's already-completed results.

Decision (2026-08-06): mark CRISPRessoPooled + CRISPRessoCompare (pooled)
as Superseded by the Merged track rather than fixing them. CLAUDE.md
previously claimed CRISPRessoCompare (pooled) was "✅ Done" — this was
WRONG (its only possible input, crispresso/pooled/trimmomatic/, was
confirmed empty); that claim was never validated against actual output
before being written down.

Lesson: don't mark a step "done" because a dependent step "ran" - verify
the actual output directory has real content, especially for multi-step
chains where step N's success message doesn't guarantee step N-1 actually
produced usable input.
```

### NextPolish — N-content Rejection
```
Job 692710 (nextpolish_genome.sh, first attempt) ran ~1h then failed at the
very first step (db_split, before any real alignment/polishing work):
  "Too many[0.109336] reads contains N base, please do QC first."
NextPolish's own seq_split tool refuses to proceed if too large a fraction
of input reads contain an N base anywhere - 10.9% of reads did here,
likely related to NovaSeq X's lenient quality-trimming thresholds already
documented in this project (aggressive trimming isn't appropriate given
its binned quality scores - see "fastp" Known Issue above) - permissive
trimming leaves more low-confidence/N-called bases in than a platform with
continuous quality scores would.

Fix: add -N to sgs_options in run.cfg (tells seq_split to keep N-containing
reads instead of refusing to run) - confirmed via reading
config_parser.py: `if '-N' in self.cfg['sgs_options']:
self.cfg['sgs_rm_nread'] = 0`, which is passed straight through to
seq_split's own -N flag ("don't discard a read/pair if the read contains N
base"). Chosen over pre-filtering N-containing reads out entirely (would
discard whole read pairs rather than just the ambiguous bases within them,
and adds an extra preprocessing step) - see codes/assembly/nextpolish_genome.sh.

Lesson: NextPolish's read-splitting step does its own data-quality gate
independent of whatever upstream QC/trimming already happened - a tool can
have a stricter internal threshold for a specific metric (N-content here)
even when the input already passed this project's established trimming
pipeline for a different, unrelated purpose (quality-score-based trimming,
not N-content filtering).
```

### NextPolish — `forkserver` vs `fork` Multiprocessing (Real Package Bug)
```
Job 692760 (nextpolish_genome.sh, with the -N fix applied) got much further
- past db_split, through alignment, ~2h44m in - then failed at the actual
polish_genome step, all 4 parallel workers:
  NameError: name 'FUN' is not defined
  File ".../nextpolish-1.4.1/lib/nextpolish1.py", line 182, in worker
    c_seq = FUN(seed_name, CFG)

Root cause: a genuine bug in NextPolish 1.4.1 itself (not our config/data),
triggered by a Python version mismatch. nextpolish1.py's main() sets `global
CFG, FUN` right before creating a `Pool(...)` (nextpolish1.py:218-223),
relying on true `fork` semantics - worker processes inheriting a
copy-on-write snapshot of the live parent process's globals at the moment
Pool() runs. Confirmed via
`/hpcfs/.../nextpolish_env/bin/python -c "import multiprocessing;
print(multiprocessing.get_start_method())"` → this Python 3.14 build
defaults to `forkserver`, not `fork`. Under forkserver, worker processes
fork from a long-lived server process whose module state predates main()'s
later `FUN = {...}[args.task]` assignment - so workers never see it. This
is a compatibility gap between NextPolish 1.4.1's code (written assuming
fork-everywhere, true of most Python versions historically on Linux) and
whatever changed this specific build's default start method.
nextpolish2.py has the identical `global ... ; Pool(...)` pattern
(confirmed via grep) and would hit the same bug if a later task reaches it
- patched defensively even though only nextpolish1.py was actually
exercised so far.

Fix: patched both vendored files directly (NOT a config change - there's
no run.cfg option for this) -
/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso/envs/nextpolish_env/share/nextpolish-1.4.1/lib/{nextpolish1,nextpolish2}.py
- added `import multiprocessing` +
`multiprocessing.set_start_method('fork', force=True)` under
`if __name__ == '__main__':`, right after the existing imports. Verified
by re-running the exact failed command
(`python nextpolish1.py -p 8 -g ... -t 1 -s sgs.sort.bam -l lgs.sort.bam -o
/tmp/test_polish_out.fasta`) directly - produced 28MB+ of real polished
sequence within 2 minutes with no error, vs. instant NameError before the
patch. Full job resubmitted (workdir wiped first, rewrite=yes in our
config means no partial-state resume anyway) - job 692760's ~2h44m of
alignment work had to be redone from scratch.

Lesson: a third-party tool's parallelism can silently break across Python
version upgrades if it relies on implicit fork-inherited global state -
this class of bug won't show up in the tool's own test suite unless tested
against the same Python version/build users actually run it with. When a
NameError points at a variable that IS assigned earlier in the same
function, suspect a multiprocessing start-method mismatch, not a typo -
check `multiprocessing.get_start_method()` before assuming the vendored
code is simply broken.
```

### NextPolish — No Improvement (Both Bugs Fixed, Result Still Negative)
```
After both NextPolish bugs above were fixed, job 699678 completed cleanly
(2026-08-10, ~9h) and produced genome.nextpolish.fasta. QC comparison (job
700892 QUAST, job 700893 BUSCO, both 2026-08-11) against the pre-polish
baseline:

Metric                    Original      Polished      Change
N50                       28.31 Mb      28.30 Mb      ~same
Genome fraction           82.807%       82.038%       WORSE (-0.77pp)
Duplication ratio         1.078         1.087         WORSE
# misassemblies           7,589         8,043         WORSE (+6%)
NA50 (aligned N50)        198,159       193,546       WORSE
Mismatches/100kbp         567.95        557.72        better (~1.8%)
Indels/100kbp             137.54        136.60        better (~0.7%)
BUSCO Complete            87.1% (3169)  87.1% (3173)  ~same
BUSCO w/ internal stops   220           221           ~same (no fix)

Decision: polished genome NOT adopted. colombian_scaffolded.fna (original)
remains authoritative. No file replacement, no re-run of Liftoff.

Why polishing didn't help: NextPolish corrects disagreements between the
assembly and aligned reads - but polishing used the SAME reads already
used to build the assembly with SPAdes. Wherever the assembly already
matched what those reads say, there was nothing to correct, even where
that sequence differs from the BUSCO ortholog gene models. The small
increase in misassemblies suggests polishing even introduced a few new
small-scale artifacts (plausibly from indel edits shifting local alignment
behavior against the already-fragmented reference). This reframes the 220
internal-stop-codon BUSCO genes as more likely genuine population
divergence or gene-prediction quirks than assembly errors.

Lesson: short-read polishing has a low ceiling when polishing with the
exact same reads that built the assembly - there's no new information for
the polisher to act on. Meaningful further improvement needs genuinely new
data (long reads) or a different assembly strategy (e.g. per-individual
assembly to avoid cross-individual heterozygosity), not reprocessing the
same short reads through a different tool. See "Genome Assembly — Further
Improvement Options" in Pending Analyses for what's next.
```

### TGS-GapCloser — Multi-stream Gzip Crashes the FASTQ Reader
```
Job 705937 (tgsgapcloser_genome.sh, roadmap option 2) "completed" in
2m50s - too fast to be real (compare: QUAST/BUSCO alignment steps on this
genome take minutes-to-hours). Checked the log rather than trusting the
COMPLETED status (same lesson as the CRISPRessoAggregate "0 folders"
Known Issue): the actual TGSGapCandidate binary aborted -

  tgsgapcandidate: ../biocommon/fasta/fasta.cpp:63: void
  BGIQD::FASTA::Id_Desc_Head::Init(const string&): Assertion
  'line.size() > 1' failed.
  Aborted (core dumped)

- while "LoadONTReads", producing a 0-byte .ont.fasta.

Root cause: the combined Nanopore reads file was built with
`cat run1/*.fastq.gz run2/*.fastq.gz > combined.fastq.gz` - valid gzip,
but a MULTI-STREAM file (~490 concatenated gzip members, one per source
fastq_pass file). TGS-GapCloser's internal FASTA/FASTQ reader (BGIQD
library, hand-rolled C++) doesn't handle multi-stream gzip boundaries
correctly - it almost certainly hit a blank/malformed line right at a
stream transition and its own internal assertion caught it (crashed
instead of silently producing wrong output, which is at least honest).

Fix: build the combined file via `zcat run1/*.fastq.gz run2/*.fastq.gz |
gzip > combined.fastq.gz` instead - zcat correctly decompresses
multi-member gzip (standard, robust), and re-piping through a fresh
`gzip` produces a single-member output stream that naive downstream
parsers can't misread. Also added a line-count sanity check
(divisible-by-4, >=1M reads) right after combining, so a similar problem
would fail loudly and immediately next time instead of silently producing
a small/malformed file that only surfaces as a confusing crash deep
inside a third-party binary later.

Lesson: `cat file1.gz file2.gz > combined.gz` is technically valid gzip
but is a known landmine for naive/custom decompression code (common in
older bioinformatics C/C++ tools that hand-roll their own gzip reading
rather than using zlib's stream-aware APIs correctly). Prefer
`zcat ... | gzip > combined.gz` whenever concatenating compressed FASTQ/FASTA
for a tool whose gzip-handling robustness is unknown - the cost is
negligible (one extra decompress/recompress pass) and it eliminates this
entire class of failure. Also: a fast "COMPLETED" SLURM status is not
proof of success - check the tool's own log for what actually happened,
especially for tools with multi-step internal pipelines where later steps
can fail invisibly to SLURM's own exit-code tracking if error propagation
is imperfect.

CORRECTION (2026-08-20): this fix was necessary but NOT sufficient. Job
710309, run with the clean single-stream gzip fix above, hit the EXACT
SAME assertion crash again - proving multi-stream gzip was not the (whole)
root cause. See next Known Issue for what it actually was. Keeping the
zcat|gzip rebuild anyway since it's still correct practice and does no
harm, but it alone did not fix this.
```

### TGS-GapCloser — Wrapper Always Uses `--ont_reads_a` (Wants FASTA, Not FASTQ)
```
After the multi-stream gzip fix (previous Known Issue) still didn't
resolve the crash, and after also fixing a second, unrelated bug (job
707123 silently resumed from a PREVIOUS failed run's leftover
done_step1_tag/done_step2.1_tag marker files - TGS-GapCloser writes these
into the CURRENT WORKING DIRECTORY, not under --output, completely
undocumented; fixed with `rm -f done_step*_tag` at the top of the script),
job 710309 still hit the identical assertion crash with a CORRECTLY BUILT
single-stream gzip FASTQ file and fresh (non-resumed) state. This ruled
out both prior hypotheses and pointed at the actual root cause:

The tgsgapcloser wrapper script unconditionally invokes its internal
binaries with `--ont_reads_a` (the FASTA-format flag) - confirmed via
`grep -n "ont_reads_a\|ont_reads_q" .../bin/tgsgapcloser`: `--ont_reads_a`
appears 3 times, `--ont_reads_q` appears ZERO times anywhere in the
wrapper. The internal tgsgapcandidate binary's own --help (run directly:
`.../tgsgapcloserbin/tgsgapcandidate` with no args) reveals two separate,
mutually exclusive flags: `--ont_reads_q "the ont reads in fastq format"`
vs `--ont_reads_a "the ont reads in fasta format"` - but the wrapper only
ever uses the FASTA one, regardless of what format the user's `--reads`
file actually is. The top-level `tgsgapcloser --help` never mentions this
requirement ("--reads <tgs_reads_file> input TGS read file" - no format
specified), so feeding it raw Nanopore FASTQ (the only format ONT
basecalling produces) silently sets up a crash: the FASTA-mode parser
chokes on '@'/'+' FASTQ header lines it isn't expecting.

Fix: convert reads to FASTA before calling tgsgapcloser -
  zcat combined.fastq.gz | awk 'NR%4==1 {print ">"substr($0,2)}
  NR%4==2 {print}' | gzip > combined.fasta.gz
Added a sanity check (FASTA sequence count via `grep -c "^>"` must equal
the original FASTQ read count) to catch a bad conversion immediately
rather than downstream. Verified the fix directly against the actual
failing binary before resubmitting the full pipeline: reused the
already-computed (and valid - produced by minimap2, which handles FASTQ
input fine) .sub.filter.paf and .orignial_scaff_infos from the failed
run, ran tgsgapcandidate directly against the new FASTA reads - no
assertion crash within 5 minutes (vs. an almost-instant crash before),
strong evidence the fix works even though the direct test itself timed
out before finishing (large file, no rush to complete it - SLURM has no
such timeout).

Lesson: a tool's top-level CLI help can describe an input flag generically
("input TGS read file") while the wrapper's actual invocation of internal
binaries hardcodes a specific, undocumented format requirement - when a
parser crashes on what looks like valid, well-formed input in an
apparently-supported format, check what the wrapper script ACTUALLY passes
to its internal binaries (grep the wrapper itself), not just what the
top-level --help claims to accept. This is the third distinct real bug
found in third-party genome-improvement tools this month (NextPolish had
two - N-content rejection needing user awareness, and a genuine
forkserver/fork multiprocessing incompatibility) - treat "the tool crashed
on our real data" as a signal to read the tool's source/wrapper, not
necessarily a sign our data or config is wrong.
```

### TGS-GapCloser — Result: Genuine Gap-Filling Success (Third Attempt)
```
Job 710348 (all three fixes applied: zcat|gzip rebuild, stale-tag cleanup,
FASTQ->FASTA conversion) ran 42m51s and printed "ALL DONE !!!" - the first
of three attempts to actually complete the real pipeline (attempts 1 and 2
crashed within minutes each).

Output: assembly/tgsgapcloser_output/colombian_gapfilled.scaff_seqs (copied
to colombian_gapfilled.fasta for downstream tools - TGS-GapCloser's native
output extension isn't a standard FASTA suffix). 2064 sequences (matches
original), 694MB (vs 692MB original - grew, consistent with real sequence
replacing some N gaps). Chr0 still present, still named "Chr0_RagTag" -
unlike NextPolish, TGS-GapCloser does NOT rename sequences.

Gap-fill rate (from colombian_gapfilled.gap_fill_detail, type field F=filled
vs N=still-a-gap): 179,359 / 318,572 gap regions filled = 56.3%. A real,
substantial result given the shallow ~3.2x Nanopore depth - most gaps that
got filled likely had at least one spanning read; the 43.7% that didn't
likely had zero reads crossing them, consistent with Poisson coverage
gaps at 3.2x mean depth.

QC comparison (jobs 714184 QUAST, 714185 BUSCO, 2026-08-20) submitted
against the ORIGINAL baseline (not the rejected NextPolish version).

RESULT (2026-08-31, job 714184 completed after 18h58m — job 714185 BUSCO
completed earlier in 3m31s): a genuine trade-off, not a clean win.

Metric                    Original      Gap-filled    Change
N50                       28.31 Mb      29.43 Mb      better (+4%)
Genome fraction           82.807%       92.088%       MUCH better (+9.3pp)
Duplication ratio         1.078         1.042         better
N's per 100kbp            7,145         1,564         much better (-78%)
# misassemblies           7,589         23,532        MUCH worse (+210%)
NA50 (aligned N50)        198,159       115,382       worse (-42%)
Mismatches/100kbp         567.95        700.33        worse
Indels/100kbp             137.54        210.08        worse
BUSCO Complete            87.1% (3169)  95.4% (3476)  MUCH better
BUSCO Missing             6.7% (245)    2.1% (73)     MUCH better
BUSCO Fragmented          6.2% (226)    2.5% (91)     MUCH better

Interpretation: filling a gap replaces an N placeholder with real sequence
derived from Nanopore reads, corrected only with racon (long-read
consensus) - never reconciled against the high-precision Illumina reads.
That new sequence is where essentially all the completeness gain comes
from (genome fraction, BUSOC completeness) - but also where the new
mismatches/indels and most new misassemblies almost certainly concentrate:
a previously-invisible N gap can't register as "misassembled", but once
it's real (imperfect) sequence, QUAST can now detect local disagreement
there against the fragmented reference.

Decision (2026-08-31): NOT adopted as-is. colombian_scaffolded.fna
(unfilled) remains authoritative pending the targeted post-gap-fill
polishing follow-up - see "Targeted Post-Gap-Fill Polishing" Known Issue
below. Documented in full in
reference/colombian_scaffolded_genome/README.md ("Gap-filling experiment"
section).

Lesson: a QC comparison against a single baseline can show a real,
substantial improvement on one axis (completeness) simultaneously with a
real regression on another (structural precision) - report both rather
than picking whichever framing looks better, and let the actual downstream
use case (here: does the 8.4x-worse misassembly count actually matter for
what this genome is used for, vs. does the +9.3pp genome fraction matter
more) drive the adoption decision instead of a single composite "better/
worse" verdict.
```

### Targeted Post-Gap-Fill Polishing
```
Rationale: the gap-filling trade-off above (huge completeness gain, real
structural/precision cost) is explained by the newly-filled sequence never
having been reconciled against the high-precision Illumina reads (it came
from Nanopore + racon only). The original whole-genome NextPolish
experiment (see "NextPolish — No Improvement" Known Issue) found nothing
to correct because it reprocessed reads the assembly was ALREADY built
from - but the newly-filled regions are different: they were never built
from or compared against Illumina reads at all, so genuine, correctable
disagreements should exist there this time.

Method: re-ran nextpolish_genome.sh (identical tool/config/2-rounds recipe
as the original polishing experiment - see codes/assembly/nextpolish_genome.sh,
no code changes needed) pointed at colombian_gapfilled.fasta instead of
colombian_scaffolded.fna. This is technically a whole-genome NextPolish
run, not literally restricted to the filled coordinates - but the expected
effect is a de facto targeted correction: regions already consistent with
Illumina reads (everything that was already in the pre-gap-fill genome)
should see little-to-no change (exactly as observed in the original
polishing experiment), while the Nanopore-derived filled regions - never
previously reconciled with Illumina data - are where real corrections
should land. Simpler and lower-risk than surgically extracting/patching
just the filled coordinates (which would need careful re-anchoring if
NextPolish's own indel corrections shift local coordinates).

Output dir: assembly/nextpolish_output_gapfilled/ (kept separate from the
original assembly/nextpolish_output/ run).

Status: job 716452 COMPLETED 2026-09-06 (10h04m, started once the
medium-partition node backlog cleared) -> assembly/nextpolish_output_gapfilled/
genome.nextpolish.fasta (2064 sequences). QC added this same session
(quast_qc_gapfilled_polished.sh, busco_qc_gapfilled_polished.sh - new
scripts, same pattern as the prior QUAST/BUSCO stages; Chr0 excluded for
QUAST via genome.nextpolish.noChr0.fasta, renamed Chr0_RagTag_np1212 by
NextPolish same as the original polishing run).

BUSCO result - genuine, real improvement, exactly the kind hoped for:
| | Pre-polish gap-filled | Post-polish (this run) |
|---|---|---|
| Complete | 95.4% (3476) | 95.5% (3479) |
| **Internal stop codons** | **184** | **136** |
| Fragmented | 91 | 87 |
| Missing | 73 | 74 |

Completeness barely moved (as expected - gap-filling, not polishing, is
what fixes Missing/Fragmented), but internal stop codons - the artifact
this whole experiment targeted (Nanopore-derived indel errors causing
frameshifts) - dropped by 48 genes (184->136, -26%). This is Illumina
short-read correction doing exactly its job on sequence that was
genuinely never checked against it before.

QUAST result (2026-09-07) - structurally flat, exactly as expected for a
short-read polish (it fixes point-like indel errors, not large-scale
placement):
| | Pre-polish gap-filled | Post-polish (this run) |
|---|---|---|
| Genome fraction | 92.088% | 92.209% |
| Misassemblies | 23,532 | 23,914 |
| Duplication ratio | 1.042 | 1.042 |

Net verdict for the whole gap-fill+polish experiment: the polish is a
clean, low-risk addition on top of gap-filling - it recovers 48 genes'
worth of internal-stop-codon artifacts (real correctness win, see BUSCO
above) at essentially zero structural cost (misassembly count within
noise, genome fraction and duplication ratio unchanged). It does NOT
address the structural precision lost during gap-filling itself (the
23,532 vs. 7,589 misassemblies gap against the original pre-gap-fill
assembly remains) - that trade-off (completeness for structural risk) is
inherent to TGS-GapCloser's long-read gap-filling and was already the
known, accepted cost documented in "TGS-GapCloser - Result" above.

### Off-Target Deduplication — CasOFFinder vs CRISPOR Priority
```python
# When CasOFFinder and CRISPOR find the same off-target site within ±2bp,
# keep the CRISPOR entry (it has MIT/CFD scores and locus annotation).
# WRONG: alphabetic sort keeps CasOFFinder ('Cas' < 'CRI' alphabetically)
# CORRECT: explicitly assign priority (CRISPOR=0, CasOFFinder=1),
#          then sort ascending and drop_duplicates(keep='first')
df['priority'] = df['Source'].map({'CRISPOR': 0, 'CasOFFinder': 1})
df = df.sort_values('priority').drop_duplicates(subset='site_key', keep='first')
```

### PCR Primer Design for On-/Off-Target Validation (2026-09-08)

New pipeline (`codes/analysis/design_offtarget_primers.py` +
`run_offtarget_primer_design.sh`, one call per gene) designs PCR primers
around each site in a gene's `combined_offtargets.csv` (on-target + all
known off-targets), for gel-based indel checks and deep targeted
sequencing of the same sites already covered by the WGS off-target
analysis. Parameterized via `--gene`/`--sites-csv`/`--ref-version` —
`bdnf` is the first gene run (real edited samples exist); reusable as-is
for the other 7 candidate genes once their edits/off-target lists exist.
Design is done against the reference genome first (matching
`ko_guide_scan.py`'s pattern), then candidate primer footprints are
checked against the Colombian pseudogenome for population variants.

**`eprimer3` needs an external `primer3_core` binary not shipped by the
`emboss/6.6.0` module** — confirmed via `Died: eprimer3 uses external
program 'primer3_core' which is not in the PATH or defined as
EMBOSS_PRIMER3_CORE`. It uses the legacy boulder-IO protocol (Primer3
≤1.x), incompatible with the modern default build (2.6.1). Fix: one-time
`mamba create -n primer3_env -c bioconda -c conda-forge primer3=1.1.4`
(`codes/analysis/setup_primer3.sh`), then every caller exports
`EMBOSS_PRIMER3_CORE=<...>/envs/primer3_env/bin/primer3_core`.

**Pseudogenome coordinate drift breaks naive fixed-buffer alignment
between reference and pseudogenome windows.** Initial assumption of a
small (~100bp) buffer to locate the pseudogenome-equivalent window failed
for 6/9 bdnf sites — actual indel-driven drift is much larger and
non-uniform: the bdnf-locus chromosome (`NC_024333.1`) is 5,314bp longer
in the pseudogenome overall, with ~2,838bp already accumulated by the
~15.9Mb bdnf locus specifically. Fix: reuse the existing
`reference/pseudogenome/colombian_pseudogenome.chain` file (a byproduct
of `make_pseudogenome.sh`'s `bcftools consensus -c`) via
`CrossMap bed <chain> <bed>` (from the `crossmap_env` conda env) to
liftover the exact reference window to pseudogenome coordinates before
extraction. Fixed population-check coverage from 3/9 to 7/9 sites (the
remaining 2 are the IUPAC-blocked sites below, expected).

**IUPAC ambiguity codes in the reference FASTA silently break
`eprimer3`/`primer3_core` — and `eprimer3` still exits 0.** bdnf's
`off_target_1` and `off_target_7` returned 0 primer candidates; root
cause was `primer3_core: Error: Unrecognized base in input sequence` from
leftover IUPAC ambiguity codes (K/M/R/S/W/Y) at unresolved-heterozygous
positions in the 2014 short-read reference assembly
(`GCF_000633615.1`) — same root cause already documented for the
CRISPOR failures on agap3/grin1a/gria1a. `primer3_core` hard-rejects
**any** such code anywhere in its input window (lowercase soft-masking
alone is fine). The bug: `eprimer3` returns exit code 0 even when
`primer3_core` fails this way internally, so a bare `returncode != 0`
check silently swallows it. Fix: `design_offtarget_primers.py` now
pre-flight-checks the extraction window for non-ACGTN bases before
calling `eprimer3` (skips with an explicit `INPUT_ERROR_AMBIGUOUS_BASES`
status + exact per-code counts), and also scans stdout/stderr text for
`Error:`/`Died:` as a defensive fallback. The output CSV carries a
`design_status`/`design_note` column pair so a hard input failure is
never confused with a genuine "searched cleanly, found nothing" result
(`NO_CANDIDATES_FOUND`).

**Does this affect the already-completed GATK+CRISPResso2 WGS off-target
results for these same 2 sites? No, for both, based on direct checks:**
the ±500bp window used for *primer design* is far wider than the window
the original WGS analysis actually used (the exact `start-end` interval
for GATK's `select_offtargets.sh`, or that ±`WGS_PADDING=40` from
`combine_offtargets.py` for CRISPRessoWGS's amplicon, ~103bp). Checked
ambiguous-base presence specifically in those narrower windows:
- `off_target_1` (`NC_024331.1:5708724-5708747` core;
  `5708684-5708787` padded): **zero** ambiguous bases in either window.
  All 143 ambiguous bases found were artifacts of the much wider
  primer-design window only.
- `off_target_7` (`NC_024340.1:12200655-12200678` core;
  `12200615-12200718` padded): core interval is clean; the padded
  window has exactly **one** ambiguous base (`R`), at genomic position
  12200714 — 36bp past the core interval's end, ~4bp from the padded
  window's outer edge (i.e., in the padding margin, not near the actual
  cut site).
- Read depth/MAPQ spot-check (one markdup BAM) at both sites: normal
  coverage (~45-51x), MAPQ predominantly 60 — no sign BWA mapping itself
  was disrupted by the ambiguous bases nearby.

**Follow-up integrity sweep (2026-09-08)**: extended the same core-vs-
padded-window ambiguous-base check to all 9 bdnf sites x both reference
versions (18 site x version combinations) — the only gene with real
off-target data on disk today (the other 7 candidate genes have no
edited samples / `combined_offtargets.csv` yet, so nothing to sweep for
them). Result: **0/18 core intervals** (the exact GATK `SelectVariants`
window) have any ambiguous base; **17/18 padded windows** (the
CRISPRessoWGS amplicon window) are clean, the sole exception being the
already-identified `off_target_7`/v1 case above. v2 (the new long-read
assembly) is clean in all 18 windows, consistent with it having far less
unresolved heterozygosity than the 2014 short-read v1 assembly. No new
cases found — this confirms off_target_7/v1 remains the only actual
point of contact between this issue and the already-reported WGS
results, for everything currently on disk.

---

## Pending Analyses

### 1. CRISPResso2 On-Target — RESOLVED 2026-08-06
```
Track A — Individual (15 samples, markdup BAMs, 100bp amplicon/region):
  ✅ crispresso_ontarget.sh (done)
  ⚠️ crispresso_pooled_groups.sh — SUPERSEDED, not fixed. CRISPRessoPooled
     was the wrong tool (see Known Issues) and would exactly duplicate
     Track B's results since both tracks use identical parameters.
  ⚠️ crispresso_compare_groups.sh — SUPERSEDED along with the above.
     Previously marked "✅ Done" in this file — that was WRONG, its input
     was always empty. See "CRISPRessoPooled — Wrong Tool + Redundant"
     Known Issue for the full correction.

Track B — Merged groups (4 groups, merged BAMs, 100bp amplicon/region):
  ✅ crispresso_ontarget_merged.sh (done May 18)
  ✅ crispresso_compare_merged.sh (6 pairwise comparisons, done)
  This is now the sole/authoritative group-level on-target comparison.

Track C — Aggregate:
  ✅ crispresso_aggregate_ontarget.sh (job 692758, 15 samples, Aug 6 2026)
     New script — none existed before. See "CRISPRessoAggregate" Known
     Issue for the -p prefix / output-dir / --suppress_report fixes.
     Output: crispresso/aggregate/CRISPRessoAggregate_on_all_samples_ontarget_aggregate/
```

### 2. CRISPRessoWGS Aggregate — DONE 2026-08-06
```
✅ crispresso_wgs.sh — per-sample WGS done (15 samples × 8 sites, May 13)
✅ crispresso_wgs_aggregate.sh — job 692749, 4 groups (20-28 folders each)
   Fixed the same two bugs as crispresso_aggregate_ontarget.sh (-p prefix,
   cd to OUTPUT_DIR, --suppress_report) — see "CRISPRessoAggregate" Known
   Issue for full detail.
   Output: crispresso/wgs/trimmomatic/aggregate/CRISPRessoAggregate_on_${GROUP}_wgs_aggregate/
```

### 3. De Novo Assembly + Scaffolding + Annotation (Controls)
```
Tools: SPAdes → RagTag → QUAST/BUSCO → Liftoff (see "Assembly Pipeline" above
       for full parameter details)
Input: Control_MNP_I, II, III trimmed reads (~47GB compressed)
Purpose: Generate population-specific genome via true de novo co-assembly
         (Colombian guppy vs Guanapo/Trinidad reference), complementary to
         the existing reference/pseudogenome/ (bcftools consensus method)

Setup completed (2026-07-06):
  ✅ reference/GCF_000633615.1_annotation.gff downloaded + md5-verified
  ✅ BUSCO actinopterygii_odb10 lineage staged (odb12.2 abandoned — BUSCO
     5.7.1 hardcodes an odb10 requirement, see Known Issues)
  ✅ liftoff_env conda environment created (liftoff v1.5.1)

Phase 1 completed (2026-07-10, job 669505, after 2 failed attempts — see
  "SPAdes — mmap ENOMEM" Known Issue above for the full debugging history):
  ✅ contigs.fasta produced but highly fragmented (3.4M contigs, N50 669bp)
     — Nanopore-based gap-filling deferred, see [[nanopore_epigenome_data]]
     memory

Phase 2 completed (2026-07-11, job 683048, after 1 failed attempt — see
  "RagTag Hang on 3.4M Unfiltered Contigs" Known Issue):
  ✅ ragtag.scaffold.fasta: N50 28.4Mb, 692Mb total, 24 chromosome-scale
     scaffolds + ~2039 small unplaced fragments + 1 Chr0 (39Mb, -C flag)

Phase 3 completed (BUSCO: 2026-07-13 job 683423 after 2 failed attempts;
  QUAST: 2026-07-15 job 683925 after 1 failed attempt — see "QUAST/BUSCO"
  and "BUSCO — Hardcoded odb10" Known Issues):
  ✅ BUSCO: C:87.1% [S:86.2%,D:0.9%], F:6.2%, M:6.7% (3169/3640 complete)
  ✅ QUAST (Chr0-excluded scaffold): genome fraction 82.807%, N50 28.3Mb,
     duplication ratio 1.078, 7589 misassemblies (largely explained by
     reference fragmentation + RagTag join gaps + genuine population
     divergence from the Trinidad reference, not assembly error)

Phase 4 completed (2026-07-11, job 683050):
  ✅ colombian_scaffolded.liftoff.gff3: bdnf transferred to
     NC_024333.1_RagTag:14113623-14129072, coverage=0.945, sequence_ID=0.923

Core pipeline (Phases 1-4) complete. Final outputs:
  reference/colombian_scaffolded_genome/ (colombian_scaffolded.fna +
  colombian_scaffolded.liftoff.gff3 + README.md for external sharing) and
  assembly/qc_results/ (QUAST + BUSCO reports).

Phase 5 DONE (2026-08-11, jobs 699678/700892/700893 — see "Assembly
  Pipeline" above for full parameter details, "NextPolish — No
  Improvement" Known Issue for the full result): NextPolish short-read
  polishing tested, QC compared against baseline.
  ✅ Result: no meaningful improvement, several structural metrics
     slightly worse (genome fraction, duplication ratio, misassemblies,
     NA50). colombian_scaffolded.fna (original) remains authoritative -
     genome NOT replaced, Liftoff NOT re-run.
  → Next: see "Genome Assembly — Further Improvement Options" below for
    what's left to try (options 2-5; option 1, polishing, is now closed).
```

### 4. Coverage Plots — DONE (confirmed 2026-08-06, originally run ~May 6)
```
Scripts in coverage/csv/ (run from that directory). All 4 confirmed to have
real output (5 plots + CSV summary each, 20 plots total) — this item's
status was ambiguous in earlier versions of this file (no ✅/⏳ marker like
other sections), verified against actual output directories, not assumed:
- plot_coverage.py           → 5 plots: depth per sample, boxplot by group,
                               coverage %, heatmap, depth vs MapQ
                               → output: coverage_plots/ ✅
- plot_depth_by_zone.py      → 5 plots: depth per zone/sample/group,
                               sgRNA vs flanking scatter
                               → output: depth_zone_plots/ ✅
- plot_depth_by_position.py  → 5 plots: position-level depth with 10bp rolling avg,
                               sgRNA site annotated, ratio vs flanking
                               → output: depth_position_plots/ ✅
- plot_coverage_by_zone.py   → 5 plots: coverage metrics by zone + group summary
                               → output: coverage_zone_plots/ ✅
Cut site midpoint used: 15922048 (=(15922039+15922058)//2)
```

### 5. Genome Assembly — Further Improvement Options (Future)
```
Discussed 2026-08-06, to revisit after seeing Phase 5 (NextPolish) results.
Applies to reference/colombian_scaffolded_genome/ — see its README.md for
the quality baseline (QUAST/BUSCO) these options aim to improve, and for
required Limitations disclosure if any of these change what's shared
externally.

Ranked roughly by cost/effort vs. expected payoff:

1. Short-read polishing (NextPolish) — DONE 2026-08-11, CLOSED. No
   meaningful improvement (see Phase 5 above / "NextPolish — No
   Improvement" Known Issue). Confirmed this option's low ceiling: reusing
   the same reads that built the assembly leaves nothing new to correct
   with. Do not retry this option without genuinely new short-read data.

2. Gap-filling with existing Nanopore data — IN PROGRESS, 2026-08-20.
   TGS-GapCloser v1.2.1 (tgsgapcloser_env), 3 attempts before it actually
   worked (multi-stream gzip, stale resume tags, and a wrapper-always-uses-
   FASTA bug — all three real, see the three "TGS-GapCloser" Known Issues
   above for the debugging history). Job 710348: "ALL DONE !!!",
   colombian_gapfilled.fasta produced, 56.3% of gap regions filled
   (179,359/318,572) — a real, substantial result despite the shallow
   ~3.2x Nanopore depth. QC comparison (QUAST job 714184, BUSCO job
   714185) submitted, result pending.

3. Assemble each of the 3 Control individuals separately, then compare —
   moderate effort, no new data needed. Root cause of SPAdes's fragmentation
   (N50 669bp before scaffolding) was co-assembling 3 genetically distinct
   wild individuals at once — every heterozygous difference became a graph
   bubble. Per-individual assembly avoids cross-individual heterozygosity
   fragmenting the graph; pick the best single assembly as a new base, or
   reconcile the three. More compute than Phase 1 (3x SPAdes runs) but
   works with existing trimmed_trimmomatic/ reads.

4. New deeper long-read sequencing (Nanopore or PacBio HiFi, single
   individual, 20-30x+) — highest cost (new sequencing run), highest
   payoff. Would enable genuine long-read or hybrid assembly, likely
   resolving most of what's currently fragmented/missing (currently: 82.8%
   genome fraction, 87.1% BUSCO completeness). The real fix if this genome
   needs to go beyond "solid draft" quality (e.g. for a genome paper or a
   resource other labs build on).

5. Targeted investigation of what's missing before investing further —
   characterize the ~13% missing/fragmented BUSCO genes and ~17% of
   reference genome fraction not represented: are they concentrated in
   repetitive regions, specific chromosomes, or near regions relevant to
   the CRISPR work (bdnf locus, the 8 off-target sites)? Would clarify
   whether targeted effort beats a genome-wide fix.
```

### 6. Population-Specific Knockout Guide Comparison — DONE 2026-08-31
```
New need (not from the original pipeline): compare CRISPR knockout guide
candidates (SpCas9, NGG PAM) designed against the NCBI Trinidad reference
vs. the Colombian population genome, per gene - i.e. would a guide designed
from the reference actually work (same PAM, same seed sequence) in the
real Colombian fish, or does a population variant break/weaken it? Genes
requested: bdnf, agap3, grin1, gria1, gria2 (grin1/gria1/gria2 each have
teleost-specific "a"/"b" paralogs in this genome - both processed where
both exist; gria2a is not annotated in this assembly, only gria2b).

Scripts (parameterized, gene list configurable at the top - as requested):
  codes/analysis/run_ko_guide_scan.sh - bash wrapper, GENES=(...) array at
    the top is the thing to edit to analyze different genes; also selects
    which Colombian genome to compare against (POPULATION=pseudogenome by
    default - see rationale in the script's own header comment: preserves
    exon/intron structure almost exactly since it's reference + SNP/indel
    substitution, so gene-by-gene comparison is reliable; the alternative
    "scaffolded" de novo genome carries assembly-artifact risk, e.g. the
    internal-stop-codon BUSCO genes documented in its own README).
  codes/analysis/ko_guide_scan.py - does the actual work, one gene per
    invocation (--gene NAME --population pseudogenome|scaffolded).

No pysam/biopython dependency - crispresso2_env's samtools is broken
(missing libcrypto.so.1.0.0, pre-existing, not something this session
caused) and the script only needs Python stdlib, so it uses samtools
(module) + minimap2 (module) via subprocess instead. IMPORTANT module
load order: `module load minimap2` THEN `module load samtools/1.16.1` -
loading minimap2 second reloads an older anaconda base that shadows
samtools/1.16.1's own libs and silently swaps in the broken conda
samtools. Also: shell state does not persist between separate Bash tool
calls in this environment - module loads and the python3 invocation must
happen in the SAME call.

Pipeline per gene:
  1. Look up the gene in both GFF3s by `gene=NAME;` (Liftoff preserves the
     reference's gene naming).
  2. Extract the full gene span (both genomes), align with minimap2 --cs,
     report all variants across the gene body (includes introns).
  3. Pick a transcript ID present in BOTH genomes (longest ref CDS as the
     representative-isoform proxy) - CRITICAL: Liftoff preserves transcript
     IDs but NOT their listing order in the GFF, so naively taking "the
     first mRNA listed" in each file independently silently compares TWO
     DIFFERENT ISOFORMS between reference and population (caught during
     testing: bdnf's first-listed reference transcript was XM_008405157.2,
     but the pseudogenome's first-listed was XM_008405147.2 - same gene,
     wrong comparison, would have produced meaningless CDS-length mismatches
     and zero real variant detection). Fixed by matching on transcript ID
     (stripping the "rna-" ID prefix, which IS shared) rather than file
     order. Verified the fix by re-running bdnf: CDS length now matches
     exactly (885bp/2 exons) between reference and pseudogenome.
  4. Align CDS-to-CDS (minimap2 --cs) for base-precise coordinate mapping.
  5. Enumerate every NGG-PAM candidate (20bp spacer + NGG, both strands) in
     the reference CDS, classify each via the CDS alignment: IDENTICAL /
     PAM_BROKEN (variant in the 3bp PAM) / SEED_VARIANT (variant in the
     10bp proximal to PAM - most critical for Cas9 binding) /
     DISTAL_VARIANT (variant in the distal 10bp of the spacer) /
     NO_ALIGNMENT (window not covered by the CDS alignment).
  6. Also report population-only novel PAM sites (a population variant
     creating a new NGG absent from the reference) - invisible if you only
     ever design against the reference.
Sanity-checked before trusting the "0 CDS variants" result for
bdnf/agap3/grin1a/grin1b (biologically plausible - purifying selection on
coding sequence): cross-referenced agap3's 85 gene-body variants against
its 18 CDS exon coordinates directly - confirmed all 85 fall in
introns/UTRs, zero in CDS. Not a bug; real biology.

Results (pseudogenome, 2026-08-31):
| Gene | CDS len | Gene-body variants | CDS variants | Guides affected / total |
|---|---|---|---|---|
| bdnf | 885bp | 6 | 0 | 0/150 |
| agap3 | 3966bp | 85 | 0 | 0/557 |
| grin1a | 2919bp | 60 | 0 | 0/441 |
| grin1b | 2817bp | 0 | 0 | 0/484 |
| gria1a | 2451bp | 51 | 1 | 1/333 (1 PAM_BROKEN) |
| gria1b | 2877bp | 49 | 2 | 0/392 (variants missed all guide windows) |
| gria2b | 2691bp | 118 | 5 | 8/341 (3 PAM_BROKEN, 4 SEED_VARIANT, 1 DISTAL_VARIANT) |
| nlgn1 | 2649bp | 241 | 1 | 4/418 (4 DISTAL_VARIANT) |

nlgn1 (neuroligin-1) added 2026-09-06, no a/b paralog for this gene (unlike
nlgn2/nlgn3/nlgn4x, which do) - single unambiguous `gene=nlgn1` symbol in
both GFFs. Large gene (NC_024334.1:13466999-13845833, ~379kb) with a lot of
intronic variation (241 gene-body variants) but only 1 landing in the CDS,
affecting 4 overlapping NGG windows, all DISTAL_VARIANT (none PAM_BROKEN or
SEED_VARIANT) - lowest-severity category in this project's classification,
consistent with a reference-designed guide there likely still working in
the Colombian population.

**Correction 2026-09-06** (prompted by a user question about CRISPOR guide
design that generalizes to this tool too): find_ngg_candidates() scans the
spliced/concatenated CDS, so a 23bp window straddling an exon-exon junction
in that concatenated string isn't a real contiguous stretch of genomic DNA
- Cas9 can't actually target it. Audited all 7 genes' already-reported
candidates (codes/analysis/check_exon_junction_candidates.py): 353/2699
candidates crossed a junction, but 352 were already classified IDENTICAL
(harmless to the conclusions either way) - the ONE exception was a real
false positive, one of gria2b's two reported DISTAL_VARIANT calls (originally
9 affected guides, corrected to 8 above). Fixed going forward: ko_guide_scan.py
now classifies any junction-crossing candidate as EXON_JUNCTION_ARTIFACT
(both for the reference-CDS scan and the population-only novel-site scan)
instead of whatever it would otherwise have been classified as.

gria2b is the standout finding: a cluster of variants around CDS position
~2139 breaks or weakens candidate guides across both strands (single
variant, multiple overlapping NGG windows affected) - a reference-designed
guide there would likely fail or behave unpredictably in the real
Colombian population. Most other genes show no CDS-level difference,
consistent with strong functional conservation in these
signaling/receptor genes.

Output: analysis/ko_guide_scan/<gene>_<population>_guide_comparison.csv
(every candidate + classification) and
analysis/ko_guide_scan/<gene>_<population>_gene_body_variants.csv (full
variant list, includes introns).

Next steps (not yet done): (a) for gria2b and gria1a's flagged guides,
consider running the surviving/IDENTICAL candidates through
casoffinder.sh's existing pattern for an off-target check against both
genomes before final guide selection; (b) if more genes are needed, edit
the GENES=(...) array in run_ko_guide_scan.sh - no code changes needed.

**nlgn1 added 2026-09-06** (no a/b paralog, unambiguous `gene=nlgn1` in
both GFFs) - see results table above.

**Consolidated report - DONE 2026-09-06**: codes/analysis/build_guide_report.py
reuses ko_guide_scan.py's own functions (find_gene_features,
exon_junction_boundaries) to build a ranked-candidate JSON
(analysis/ko_guide_scan/report_data.json) across all 8 genes - top 5
CRISPRko candidates per gene (IDENTICAL classification, ranked: not in
last exon > CRISPOR mitSpecScore > offtargetCount > Doench'16, with a
position-only fallback rank for agap3/grin1a/gria1a, which lack CRISPOR
scores - see Methods below) plus the full list of variant-affected guides
to avoid. Published as an HTML artifact ("Guppy CRISPR Atlas") with the
CRISPRko vs CRISPRi guide-selection criteria the user asked about:
CRISPRko wants an early/constitutive exon (NOT simply "anywhere in the
CDS" - the last exon is specifically a bad target, since a premature stop
there often escapes NMD and yields a partially-functional truncated
protein instead of a true null), no population variant in PAM/seed, high
CRISPOR specificity, decent Doench'16 efficiency (weakest of the four
signals). CRISPRi (dCas9-KRAB) wants a narrow window around the TSS
(-50/+300bp per Gilbert 2014/Horlbeck 2016 - not "the promoter" loosely,
and strand matters less than proximity to the TSS), and essentially never
inside the CDS.

**CRISPRi TSS-window scan - DONE 2026-09-07** (user asked to add it right
after reading the criteria above): new script codes/analysis/crispri_tss_scan.py,
reuses ko_guide_scan.py's alignment/classification functions directly
(find_gene_features, faidx_seq, align_cs, parse_cs_variants,
find_ngg_candidates, classify_candidate, run_crispor) rather than
reimplementing them - only the window definition and coordinate math are
new. TSS = the representative transcript's mRNA start (not CDS start,
which is downstream of any 5'UTR) - for a "-" strand gene this is the
mRNA's END coordinate (GFF stores start<end regardless of strand). Window
extracted in genomic (+strand) coordinates then reverse-complemented for
"-" strand genes so it reads 5'->3' in the gene's own transcriptional
direction, with position 0 = TSS-50 and the last position = TSS+300 for
BOTH strands after this correction - i.e. tss_relative_position =
-50 + local_0based_index uniformly. No exon-junction check needed here
(genomic DNA is always contiguous, unlike the spliced CDS). Classified
against the pseudogenome exactly like the CDS scan (IDENTICAL/PAM_BROKEN/
SEED_VARIANT/DISTAL_VARIANT/NO_ALIGNMENT) and CRISPOR-scored the same way,
except Doench'16 is deliberately NOT surfaced for CRISPRi candidates in
the report - it's trained to predict cutting/indel efficiency, not
dCas9-KRAB silencing, so it isn't a valid signal there; ranking uses only
specificity (mitSpecScore/offtargetCount) then proximity to the TSS.

Run for all 8 genes (pseudogenome population): 36-71 NGG candidates per
gene in the 351bp window (much smaller n than the CDS scan, as expected).
Only **agap3** hit the IUPAC-ambiguity-code CRISPOR crash here (see #7) -
grin1a and gria1a, which DO hit it in their CDS, turned out clean in their
(much smaller, different genomic location) TSS windows, confirming the
ambiguity codes are scattered through the reference assembly rather than
gene-wide. **gria2b** has a real population variant in its TSS window (2
PAM_BROKEN + 1 SEED_VARIANT) - the only gene with a CRISPRi-relevant
finding; all other genes are fully IDENTICAL in this window. Wired into
build_guide_report.py (report[gene]["crispri"]: total_guides,
classification_counts, variant_affected_guides, top_candidates,
crispor_available) and the "Guppy CRISPR Atlas" artifact - each gene card
now shows CRISPRko and CRISPRi candidate tables side by side, plus a
second executive-summary table for CRISPRi. Output CSVs:
analysis/ko_guide_scan/<gene>_pseudogenome_crispri_candidates.csv.

**Guide source switched to CRISPOR's own predictions - DONE 2026-09-08**
(user request: show CRISPOR-predicted guides instead of the manual NGG
scan's own ranking, with their real metrics and real off-targets).
build_guide_report.py rewritten: the candidate list and its displayed
metrics/off-targets now come directly from CRISPOR's own guide/offtarget
TSVs (`<gene>_reference_crispor_guides.tsv` / `_crispor_offs.tsv` for KO,
`_crispri_crispor_guides.tsv` / `_crispri_crispor_offs.tsv` for CRISPRi -
these were already being written to disk by run_crispor(), just never
read back beyond 3 scalar fields before). The manual scan's CSVs
(guide_comparison.csv / crispri_candidates.csv) are still required and
still used - but now only as the join key for population-safety
classification and gene-structure position (CDS%/last-exon or TSS
offset), since CRISPOR's own guide list has no concept of exon structure
or the pseudogenome at all. Join key: exact 23bp targetSeq (spacer+PAM) -
both scans use the same NGG/20bp rule on the same input, confirmed to
match essentially 1:1 in practice.

Now surfaced per guide: CRISPOR's own guideId, mitSpecScore, cfdSpecScore
(not used before), offtargetCount, Doench'16 + Moreno-Mateos (KO only -
Moreno-Mateos/CRISPRscan is specifically the model recommended for guides
made by in-vitro T7 transcription, per crisporWebsite's own README, which
is how this project's guides are actually made for embryo microinjection -
worth weighting alongside Doench'16, not just as a footnote), and the top
3 REAL off-targets (genomic position, strand, mismatch count, CFD score)
parsed from the offtarget TSV and sorted mismatches-ascending/CFD-
descending (worst-case-first). Fallback (agap3 both contexts; grin1a/
gria1a KO context only - the known IUPAC-ambiguity CRISPOR crash, see #7)
still shows manual-scan candidates ranked by position only, clearly
labeled per-row ("sin datos de CRISPOR para esta guía").

### 7. CRISPOR Integration into ko_guide_scan.py — DONE 2026-09-05
```
Goal: complement (not replace) the manual PAM-scan/variant-classification
pipeline above with real efficiency (Doench'16/Rule Set 2) and specificity
(MIT score, real BWA-based off-target count) scoring from CRISPOR, per the
user's explicit instruction to keep the manual scan as the primary tool.

Docker vs Singularity: Docker needs the `docker` group (root-equivalent
daemon access) - confirmed unusable for this user (`docker ps` -> permission
denied). Singularity (module singularity/3.7.1, already on this cluster)
pulls the same Docker Hub image directly (`singularity pull docker://...`)
without needing Docker or root - the correct path here.

Docker Hub multi-arch bug: the `maximilianh/crispor` "latest" tag (pushed
2026-03-01) is ARM64-ONLY (broken build, likely pushed from Apple Silicon) -
`singularity exec ... uname -m` on this cluster returned "the image's
architecture (arm64) could not run on the host's (amd64)". Confirmed via
`curl https://hub.docker.com/v2/repositories/maximilianh/crispor/tags` that
`v5.2`/`v5.2c` (pushed 2026-01-21, earlier but correctly multi-arch) have
both amd64/arm64. Used v5.2c. Image kept at
codes/analysis/crispor_singularity/crispor_v5.2c_amd64.sif.

Custom genomes registered via crisporAddGenome (see crispor_add_genomes.sh,
crispor_add_genome_v2.sh), bind-mounted at
codes/analysis/crispor_singularity/genomes/ -> /data/genomes inside the
container (the container's own genome dirs are empty/read-only, and a
Singularity container is otherwise ephemeral, so this bind mount is what
makes registered genomes persist across invocations). Registered IDs:
guppyRefTrinidad, guppyColPseudogenome, guppyRefMaleV2 (v2 reference; v2
pseudogenome pending, see #8 below).

ko_guide_scan.py's run_crispor() feeds each gene's already-extracted CDS to
`crispor.py <genomeId> <cds.fa> <guides.tsv> -o <offs.tsv>` and matches
CRISPOR's guides back to the manual scan's own candidates by EXACT 23bp
spacer+PAM sequence (both use the same NGG/20bp definition, reported 5'->3'
on the guide's own strand) - confirmed via `crispor.py noGenome` test that
this container version's actual TSV header includes targetSeq,
mitSpecScore, offtargetCount, "Doench '16-Score" (all consumed here), even
though it otherwise differs from the older sample TSVs bundled in the
crisporWebsite repo docs.

**Follow-up 2026-09-08 - why agap3/grin1a/gria1a have NO CRISPOR TSVs at
all** (`crispor_available: false` in report_data.json; the manual PAM-scan
CSVs still exist and are unaffected, per #6). Re-ran ko_guide_scan.py for
all 3 to capture the actual crispor.py traceback (previously only the
truncated last-800-char stderr tail was ever seen, not logged to disk).
Root cause confirmed as the same recurring IUPAC ambiguity-code issue
(K/M/R/S/W/Y in the 2014 short-read reference, see primer-design and
combine_offtargets Known Issues above) - but here it's a HARD, UNCAUGHT
crash inside crispor.py itself, in two different internal code paths that
don't guard against non-ACGTN input at all (unlike the Azimuth soft-fail
already noted for individual guides overlapping an ambiguous site -
gria1a's case is actually this SAME Azimuth code path, just crashing the
whole run instead of one guide):
  - agap3, grin1a: `revComp()` (called from `writePamFlank`/
    `flankSeqIter` during off-target/flank sequence extraction) ->
    `KeyError: 'S'` - its complement-base dict has no 'S' entry.
  - gria1a: Azimuth-2.0's `nucleotide_features()` (one-hot encoding for
    the Doench'16 efficiency model) -> `ValueError: 'R' is not in list` -
    its fixed alphabet list has no 'R' entry.
Both failures happen against BOTH guppyRefTrinidad and
guppyColPseudogenome for all 3 genes - confirms the ambiguous base is
physically present in the v1 reference assembly near these genes, not
introduced only by the pseudogenome's `bcftools consensus` step.
`run_crispor()`'s existing try/except-equivalent (checking `returncode !=
0`) already handles this gracefully at the ko_guide_scan.py level (prints
a WARNING, returns `{}`, the manual scan's classification is entirely
unaffected) - this is a real gap in CRISPOR's own code, not a bug in this
project's wrapper, and not fixable without patching the container's
bundled crispor.py/Azimuth-2.0 to skip or N-mask ambiguous bases before
scoring (not done - out of scope, efficiency/off-target scoring is a
supplement to the manual scan, never the primary classifier, per this
item's original goal statement above).
```

### 8. Migration to GCF_904066995.2 (v2) — IN PROGRESS 2026-09-05
```
Why: GCF_000633615.1 (Trinidad/Guanapo, female, short-read, 2014 - the
reference this whole project has used) is now marked "suppressed" by NCBI
(confirmed via the NCBI datasets API `assembly_info.assembly_status`).
GCF_904066995.2 (P_reticulata-male-v2, University of Exeter, PacBio+Hi-C,
released 2025-08-03) is the current RefSeq "reference genome" for the
species, with drastically better contiguity/completeness: contig N50
41.9kb -> 7.94Mb, scaffold N50 5.27Mb -> 32.9Mb, BUSCO complete 96.84% ->
98.69% (cyprinodontiformes_odb10, both via the NCBI API, not re-run here).

Of the 7 knockout-guide genes (#6 above): grin1b sits on an UNPLACED
scaffold in the old assembly (NW_007615041.1, confirmed via grep of both
GFFs) - the new assembly may finally anchor it to a real chromosome.
grin1a sits on LG12 (NC_024342.1 in v1; confirmed via NCBI
`sequence_reports` that LG12 is the guppy XY sex chromosome) - since every
experiment in this project uses only females (XX, no Y), any off-target or
Y-specific sequence CRISPOR/GATK finds there against the v2 (male) genome
needs an interpretation caveat: it may not exist in the real animals. The
other 5 genes are autosomal - a clean improvement, no caveat needed.

User explicitly asked for a PARALLEL replica: same scripts, parametrized by
genome version, writing to sibling paths, WITHOUT touching any v1 path or
output. Implementation: a single new file, codes/genome_versions.sh, is
`source`d by every reference-dependent script; it sets REF/REF_GFF/
INTERVALS/OUT_SUFFIX based on `REF_VERSION` (env var, defaults to "v1" =
byte-identical to pre-migration paths). Scripts append `${OUT_SUFFIX}` to
their output directory names (e.g. `gatk/trimmomatic${OUT_SUFFIX}/...`,
`crispresso${OUT_SUFFIX}/...`). Full architecture/plan document:
/hpcfs/home/ing_civil/da.martinez33/.claude/plans/hi-claude-while-i-fluttering-graham.md

Setup done so far:
  - reference/GCF_904066995.2_annotation.gff downloaded + md5-verified
    (adapted codes/assembly/00_download_gff_v2.sh from the v1 version).
  - reference/intervals_v2.list: 23 real chromosome accessions (NC_088830.1
    - NC_088852.1) from NCBI `sequence_reports`, NOT the fragile
    `grep '>' REF | head -24` fallback in genomics_db_import.sh (which
    v1's own intervals.list quietly depends on being pre-built anyway - the
    NW_007615013.1 entry in v1's list is just whatever sequence happened to
    be 24th in file order, NOT the mitochondrion; confirmed the real MT
    accession NC_024238.1/NC_024238.1 is absent from v1's intervals.list
    and is intentionally excluded here too for v2).
  - reference/GCF_904066995.2_P_reticulata-male-v2_genomic.fna indexed:
    samtools faidx, gatk CreateSequenceDictionary, bwa index (~11.5 min).
  - Parametrized with genome_versions.sh (source + ${OUT_SUFFIX} renames):
    codes/mapping/{bwa_index.sh,bwa_trimmomatic_array.sh},
    codes/variant_calling/{index_dict_reference.sh,mark_duplicates.sh,
    haplotype_caller_scatter.sh,merge_sample_gvcfs.sh,reindex_gvcf.sh,
    genomics_db_import.sh,genotype_gvcf.sh,variant_filtration.sh,
    select_offtargets.sh}, codes/CRISPResso/{casoffinder.sh,
    extract_amplicon_sgRNA.sh,crispresso_ontarget.sh,
    crispresso_ontarget_merged.sh,crispresso_wgs.sh},
    codes/analysis/{make_pseudogenome.sh,hotspot_windows.sh,
    hotspot_analysis.py,plot_hotspot_summary.py},
    codes/assembly/{liftoff_pseudogenome.sh,verify_pseudogenome.sh}.
    haplotype_caller_scatter.sh's `#SBATCH --array=1-360` is v1-specific
    (15 samples x 24 intervals.list lines); v2 needs `sbatch
    --array=1-345%40 --export=ALL,REF_VERSION=v2 ...` (15 x 23) since
    SBATCH directives can't read the interval file at parse time.

  - Bug fixes found/applied while porting (not new bugs introduced by the
    migration itself - pre-existing issues this work surfaced):
    - crispor_add_genomes.sh: looked COMPLETED (exit 0) on its first real
      run but was a SILENT FAILURE - crisporAddGenome stages everything
      under /tmp/<genomeId>/ and only copies the finished result into
      --baseDir at the very end; the GFF->gene-locus conversion step
      crashed (container's bundled UCSC binaries bedSort/bedToExons/
      genePredToBed are missing libpng12.so.0, a stale Ubuntu base-image
      dependency) with an AssertionError, so nothing was ever copied out -
      but the outer script had no `set -e`, so it still printed its final
      echo/find lines and exited 0. Fixed: (1) register genomes fasta-only
      (no --gff) - ko_guide_scan.py's run_crispor() never reads
      targetGenomeGeneLocus anyway, so this loses nothing it needs; (2)
      added `set -euo pipefail` so a real failure now actually shows as
      SLURM state FAILED; (3) added a post-run check that .2bit/.fa.bwt
      actually exist before declaring success. Also found and fixed: /tmp
      is NODE-LOCAL and persists across job submissions on the same node
      (confirmed via `sacct --format=NodeList` - three consecutive
      submissions all landed on nodea-1), so a crashed run's orphaned /tmp
      staging dir survives and blocks the next attempt regardless of
      cleanup done from the login node; fixed by having the script
      unconditionally `rm -rf` its own known staging dirs at the very
      start (idempotent regardless of which node picks it up).
    - codes/CRISPResso/casoffinder.sh: referenced `${CONDA_BASE}` before
      defining it (real pre-existing bug, harmless only because the
      cluster happened to have a stale CONDA_BASE in the environment
      already) - fixed to match the correct pattern already used in
      crispresso_wgs.sh (define CONDA_BASE explicitly before sourcing its
      conda.sh).
    - codes/assembly/verify_pseudogenome.sh: was checking against
      `colombian_pseudogenome.gff3` (the CrossMap output, superseded - see
      #pseudogenome README, 64% broken parent/child hierarchy) instead of
      `.liftoff.gff3` (the primary annotation, 99.5% transfer, 0 orphans) -
      fixed while parametrizing.

  - On-target site relocation (bdnf sgRNA, guide TGAGAGACGCCCCGGGCATG+NGG):
    the CRISPResso2 on-target scripts use a hardcoded coordinate window
    (not a sequence search), so it must be relocated by hand per genome.
    Method: minimap2 --cs (-x sr preset) aligning the v1 window sequence
    against the new bdnf gene span (NC_088832.1:15848607-15863146, from
    the newly-downloaded v2 GFF). Both the 61bp and 101bp windows aligned
    with 100% identity (NM:i:0) - the sgRNA/cut-site region is fully
    conserved between assemblies. Relocated coordinates:
      v1 60bp window  NC_024333.1:15922011-15922071 -> v2 NC_088832.1:15849666-15849726
      v1 100bp window NC_024333.1:15922000-15922100 -> v2 NC_088832.1:15849655-15849755
      v1 20bp cut site NC_024333.1:15922039-15922058 -> v2 NC_088832.1:15849694-15849713
    New FASTAs written: reference/amplicon_bdnf_{60,100}bp_v2{_rc,}.fa.
    Cas-OFFinder itself needed NO relocation (it searches the new genome's
    whole sequence directly for the same guide+PAM string - reference-
    agnostic by construction); only the coordinate-based CRISPResso2
    on-target scripts needed this treatment.

  - Off-target discovery redesign: the existing pipeline's "CRISPOR side"
    of off-target discovery (convert_crispor_offtargets.py) reads a
    `data/crispor_offtargets.xls` that was manually downloaded from the
    crispor.org website for the OLD genome one time - there is no
    equivalent file for v2 and none can be "pointed at" a new path, it has
    to be regenerated. Decision (user-approved): use the crispor.py
    container set up in #7 instead of the website - run it against
    guppyRefMaleV2 with the bdnf amplicon sequence, then adapt
    convert_crispor_offtargets.py to parse crispor.py's own offtarget TSV
    (columns already characterized in #7) instead of the .xls, matching by
    guideSeq instead of the old CRISPOR-website guideId ("326forw", which
    won't exist in a fresh container run). DONE 2026-09-06.

Off-target discovery for v2 - DONE 2026-09-06 and cross-validated:
Cas-OFFinder (genome-only, no mapping/GATK needed) found 9 sites for the
bdnf guide against the new genome (8 off-target + the on-target itself, 0
mismatches, at NC_088832.1:15849690 - exactly where relocation predicted).
crispor_offtarget_scan.sh (new script) ran crispor.py in the container
against guppyRefMaleV2 with the 100bp genomic amplicon window as input
(NOT a spliced transcript - see below) and independently found the same 8
off-targets with MIT/CFD scores. combine_offtargets.py (parametrized,
REF_VERSION-aware BDNF_CHROM for the on-target-exclusion check) merged
both sources: all 8 CRISPOR off-targets matched a Cas-OFFinder hit
(found_by_both=8), confirming both tools agree completely for this guide
on the new genome. Outputs at crispresso_v2/offtargets/combined/.
NOTE: guppyRefMaleV2/guppyColPseudogenome were registered fasta-only (see
#7's libpng12 bug), so crispor.py's offtarget TSV "guideSeq" column holds
the full 23bp spacer+PAM on-target sequence (not just the 20bp spacer) -
convert_crispor_offtargets.py filters on that full sequence, not guideId.

User question 2026-09-06, investigated before running any of this: does
the bdnf transcript re-annotation between genomes (RefSeq's Gnomon pipeline
relabeled XM_008405147 from "variant X1" to "variant X2" because a brand
new isoform, XM_081605136.1, was predicted from the better assembly and
took the X1 slot) affect the previously-selected guide's validity, and
what should be fed to CRISPOR (full gene / exons / one exon)? Findings:
all 8 pre-existing bdnf transcripts persist unchanged (same root accession,
version bump only, same exon count) in v2, just renumbered by one; ALL 9
v2 transcripts (old lineage + the new X1) share the IDENTICAL terminal exon
NC_088832.1:15848607-15849795, and the guide site sits safely inside it
(899bp/82bp from either boundary) - so the guide is isoform-independent,
not affected by the relabeling. Read crispor.py's source
(extendAndGetSeq/getExtSeq): efficiency scores are computed from genomic
flanking sequence re-fetched via the aligned position (twoBitToFa), NOT
from the input file's own sequence - so the input should always be
genomic DNA (gene body or exon-of-interest + generous real flanking
sequence), never a spliced transcript/isoform sequence, which risks
placing a candidate PAM across an exon-exon junction (not real, contiguous
genomic DNA - see the ko_guide_scan.py finding immediately below, which is
the same risk in a different tool).

**ko_guide_scan.py correctness fix, same root cause, found while answering
the question above:** find_ngg_candidates() scans the CONCATENATED CDS
(exons spliced together) - a candidate window straddling an exon-exon
junction there isn't real contiguous genomic DNA either, so it's not an
actually editable Cas9 target no matter how it gets classified. Audited
all 7 previously-reported genes (codes/analysis/check_exon_junction_candidates.py):
353/2699 candidates crossed a junction; 352 were already classified
IDENTICAL (harmless), but ONE was a genuine false positive - one of
gria2b's two DISTAL_VARIANT calls (originally reported as 9 affected
guides, corrected to 8: 3 PAM_BROKEN + 4 SEED_VARIANT + 1 DISTAL_VARIANT -
see #6's results table, now corrected). Fixed: ko_guide_scan.py now has
exon_junction_boundaries()/crosses_exon_junction() and classifies any
junction-crossing candidate as EXON_JUNCTION_ARTIFACT (both the reference
scan and the population-only novel-site scan) instead of whatever it
would otherwise have been classified as. All 7 genes re-run against v1
pseudogenome with the fix - CSVs on disk now reflect the correction.
Separately noticed while re-running with CRISPOR scoring: the pseudogenome
consensus (bcftools consensus -H A) leaves IUPAC ambiguity codes at some
sites (e.g. 114 R/Y/S/W/K/M codes in gria1a's CDS region) - CRISPOR's
Azimuth efficiency model crashes on non-ACGT characters (soft-fails per
guide, already handled gracefully by run_crispor()'s existing try/except -
not a blocker, just a known gap in efficiency scoring for guides
overlapping an ambiguous site).

Registered/done as of this entry (check `sacct -j <id>` for current state
before relying on any SLURM job below):
  - guppyRefMaleV2 CRISPOR genome: registered (crispor_add_genome_v2.sh).
  - Cas-OFFinder + CRISPOR off-target discovery for the bdnf guide vs v2:
    done and cross-validated (see above).
  - Mapping+GATK v2 chain SUBMITTED (dependency-chained, not yet complete -
    this is the multi-day bottleneck, historically the longest part of
    this whole project): bwa_trimmomatic_array.sh -> mark_duplicates.sh ->
    haplotype_caller.sh -> genomics_db_import.sh -> genotype_gvcf.sh ->
    variant_filtration.sh, all with REF_VERSION=v2 exported through the
    chain (run_variant_calling.sh itself also now accepts REF_VERSION).

Still pending: pseudogenome v2 (make_pseudogenome.sh needs the v2 filtered
VCF from the GATK chain above, once it completes); ko_guide_scan.py's
"pseudogenome_v2" GENOME_CHOICES entry is already wired in (fasta/gff
paths + CRISPOR genome ids) but the pseudogenome files don't exist yet;
gatk_offtarget_genotypes.py/plot_editing_comparison.py's 8 hardcoded
OT-site coordinates need updating to the new v2 off-target list (now
available: crispresso_v2/offtargets/combined/combined_offtargets.csv, 8
sites); Phase 2 (RagTag re-scaffolding of the ALREADY
polished/gap-filled Colombian contigs against the new reference - SPAdes/
NextPolish/TGS-GapCloser do NOT need to be re-run) deliberately deferred
until Phase 1 is done, per the user-approved plan.

Status check-in 2026-09-08: HaplotypeCaller (job 716490, array 1-15%8) at
~56% by genomic position across the 8 currently-running tasks (samples.txt
lines 1-8, started 2026-09-06 11:45); lines 9-15 still queued behind the
8-task concurrency limit. Downstream GenomicsDBImport/GenotypeGVCFs/
VariantFiltration (716491-716493) remain PENDING (Dependency). Rough
estimate ~4-5 more days for all 15 to clear HaplotypeCaller.

Considered and explicitly rejected: branching off a Control-only (3-sample)
joint-genotyping side-chain to unblock pseudogenome v2/assembly Phase 2
early, since Control_MNP_I/II/III happen to be samples.txt lines 1-3 (in
the first, already-running batch) and would finish well before the other
12. Rejected because (1) genomics_db_import.sh's SAMPLE_LIST is currently
the full 15-sample samples.txt with no subset option (a real code change,
not just a config flip), and (2) more importantly, it would be a
methodological deviation from how the v1 pseudogenome was built (full
15-sample joint genotyping, THEN `bcftools view -s $CONTROLS` + fill-tags
recomputes AF from just the controls) - VariantFiltration's site-level hard
filters (QD/FS/MQ/etc.) are computed from pooled read evidence across
however many samples are in the joint call, so a 3-sample-only cohort would
compute them differently than a 15-sample cohort, breaking v1/v2
comparability. No actual bias risk either way for the pseudogenome itself:
AF_THRESH=0.667 is recomputed from the Controls' own GT values AFTER
subsetting, regardless of cohort size, so edited/treatment samples' unique
CRISPR-induced variants cannot leak into the Control consensus through this
path. Decision: keep joint-genotyping all 15 together, subset to Controls
only at the very end, same as v1.
```

### 9. PCR Primer Design for On-/Off-Target Validation — bdnf v1 DONE 2026-09-08
```
codes/analysis/design_offtarget_primers.py (+ run_offtarget_primer_design.sh),
parameterized by --gene/--sites-csv/--ref-version. See "PCR Primer Design
for On-/Off-Target Validation" in Known Issues for the primer3_core and
CrossMap-liftover fixes this required.

bdnf v1 (crispresso/offtargets/combined/combined_offtargets.csv, 9 sites):
DONE. analysis/offtarget_primers/bdnf_v1_primers.csv - 37 rows: 35 real
candidate pairs (7 sites x 5) + 2 clean placeholder rows for off_target_1/
off_target_7 (design_status=INPUT_ERROR_AMBIGUOUS_BASES, see Known Issues).
on_target's top candidate: SPECIFIC (1 amplimer), primer_variant_flag=NONE.
Population-variant check (vs colombian_pseudogenome.fna) available for 7/9
sites (the 2 IUPAC-blocked sites have no primers to check by definition).

Pending:
  - bdnf v2: blocked on reference/pseudogenome_v2/ (item #8's GATK v2
    chain must finish first; the script's --ref-version v2 path is
    already wired but pseudogenome_v2 doesn't exist yet).
  - Other 7 candidate genes (agap3, grin1a/b, gria1a/b, gria2b, nlgn1):
    blocked on those genes having real edited samples + their own
    combined_offtargets.csv - no code changes needed, same command with
    --gene/--sites-csv swapped.
```

---

## Contact and Resources

- Email: diegoandres3322@gmail.com
- GATK docs: https://gatk.broadinstitute.org
- CRISPResso2 docs: http://crispresso.pinellolab.org
- Cas-OFFinder: http://www.rgenome.net/cas-offinder/
- CRISPOR: http://crispor.tefor.net
