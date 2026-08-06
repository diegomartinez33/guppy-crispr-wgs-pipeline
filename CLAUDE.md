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

# Personal conda (CRISPResso2, CrossMap, Liftoff)
CONDA_BASE=/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate crispresso2_env # CRISPResso2, cas-offinder, pandas
conda activate liftoff_env     # Liftoff v1.5.1 (annotation transfer for scaffolded assembly)
conda activate crossmap_env    # CrossMap v0.7.3 (coordinate liftover using chain file)
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
| CRISPRessoPooled | CRISPResso2 | ⏳ Pending | Pool individual results by group |
| CRISPRessoCompare (pooled) | CRISPResso2 | ✅ Done | 6 pairwise group comparisons |
| CRISPRessoCompare (merged) | CRISPResso2 | ✅ Done | 6 pairwise merged group comparisons |
| CRISPRessoAggregate | CRISPResso2 | ⏳ Pending | All samples summary |
| Off-target prediction | Cas-OFFinder | ✅ Done | 4 mismatches, 9 sites (1 on-target) |
| Off-target prediction | CRISPOR | ✅ Done | 8 sites, 4 mismatches, guide ID: 326forw |
| Combine off-targets | Python | ✅ Done | 8 unique off-target sites |
| CRISPRessoWGS (individual) | CRISPResso2 | ✅ Done | 15 samples × 8 sites completed May 13 |
| CRISPRessoWGS aggregate | CRISPResso2 | ⏳ Pending | Aggregate per group after WGS |
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

Three distinct on-target analyses with different regions and amplicons:

| Analysis | Script | BAM source | Region filtered | Amplicon |
|---|---|---|---|---|
| Individual (15 samples) | crispresso_ontarget.sh | gatk/trimmomatic/markdup/*.markdup.bam | NC_024333.1:15922000-15922100 (100bp) | 60bp FWD + RC |
| Merged groups (4 groups) | crispresso_ontarget_merged.sh | mapping/trimmomatic/merged/*.sorted.bam | NC_024333.1:15921941-15922141 (200bp) | amplicon_bdnf_200bp_rc.fa |
| Off-target WGS (15 samples) | crispresso_wgs.sh | gatk/trimmomatic/markdup/*.markdup.bam | from offtargets_crispresso_wgs.bed | — (WGS mode) |

```
Amplicons:
  60bp:  reference/amplicon_bdnf_60bp_fwd.fa + amplicon_bdnf_60bp_rc.fa
         used for individual on-target analysis
  200bp: reference/amplicon_bdnf_200bp_rc.fa
         used for merged group on-target analysis (larger window, more context)

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

---

## Pending Analyses

### 1. CRISPResso2 On-Target (remaining steps)
```
Track A — Individual (15 samples, markdup BAMs, 60bp amplicon, 100bp region):
  ✅ crispresso_ontarget.sh (done)
  → crispresso_pooled_groups.sh (CRISPRessoPooled per group, pending)
  ✅ crispresso_compare_groups.sh (6 pairwise comparisons of pooled, done)

Track B — Merged groups (4 groups, merged BAMs, 200bp amplicon, 200bp region):
  ✅ crispresso_ontarget_merged.sh (done May 18)
  ✅ crispresso_compare_merged.sh (6 pairwise comparisons, done)

Track C — Aggregate:
  → CRISPRessoAggregate: all samples summary (pending)
```

### 2. CRISPRessoWGS Aggregate
```
✅ crispresso_wgs.sh — per-sample WGS done (15 samples × 8 sites, May 13)
→ crispresso_wgs_aggregate.sh (pending — aggregate directory empty)
   Output: crispresso/wgs/trimmomatic/aggregate/${GROUP}_wgs_aggregate/
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

Pipeline complete. Final outputs: reference/colombian_scaffolded_genome/
  (colombian_scaffolded.fna + colombian_scaffolded.liftoff.gff3) and
  assembly/qc_results/ (QUAST + BUSCO reports).
```

### 4. Coverage Plots
```
Scripts in coverage/csv/ (run from that directory):
- plot_coverage.py           → 5 plots: depth per sample, boxplot by group,
                               coverage %, heatmap, depth vs MapQ
                               → output: coverage_plots/
- plot_depth_by_zone.py      → 5 plots: depth per zone/sample/group
                               → output: depth_zone_plots/
- plot_depth_by_position.py  → 5 plots: position-level depth with 10bp rolling avg,
                               sgRNA site annotated, ratio vs flanking
                               → output: depth_position_plots/
- plot_coverage_by_zone.py   → 5 plots: coverage metrics by zone + group summary
                               → output: coverage_zone_plots/
Cut site midpoint used: 15922048 (=(15922039+15922058)//2)
```

---

## Contact and Resources

- Email: diegoandres3322@gmail.com
- GATK docs: https://gatk.broadinstitute.org
- CRISPResso2 docs: http://crispresso.pinellolab.org
- Cas-OFFinder: http://www.rgenome.net/cas-offinder/
- CRISPOR: http://crispor.tefor.net
