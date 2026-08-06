# Guppy CRISPR-Cas9 WGS Analysis Pipeline

Bioinformatics pipeline for whole-genome sequencing (WGS) analysis of
CRISPR-Cas9 editing efficiency and off-target effects in *Poecilia reticulata*
(guppy fish) from a Colombian population.

## Overview

This repository contains all scripts and code used to analyze Illumina NovaSeq X
whole-genome sequencing data from guppy fish subjected to CRISPR-Cas9 editing of
the **bdnf** (brain-derived neurotrophic factor) gene. The pipeline spans read
trimming and quality control, genome alignment, variant calling, CRISPR editing
quantification, off-target analysis, population genomics, and construction of
two population-specific reference genomes for the Colombian guppy population.

The experimental fish are from a **Colombian population** of *P. reticulata*,
sequenced and analyzed against the publicly available **Guanapo (Trinidad)
reference genome** (GCF_000633615.1). A key component of this work is addressing
the genomic divergence between these two populations in the context of CRISPR
specificity assessment.

## Experimental Design

| Group | Description | n |
|---|---|---|
| Control | No CRISPR components delivered | 3 |
| Only_MNP | Nanoparticles only, no Cas9 | 4 |
| Plasmid_Ko | Cas9 delivered via plasmid | 4 |
| RNP_Cas | Cas9 delivered as ribonucleoprotein | 4 |

- **Target gene:** bdnf (NC_024333.1, LG3)
- **sgRNA:** TGAGAGACGCCCCGGGCATG (negative strand)
- **Cut site:** ~NC_024333.1:15922041
- **Cas9:** SpCas9 NLS (NEB)
- **Sequencing:** Illumina NovaSeq X, 25B flowcell, 2×150bp paired-end
- **Library prep:** Illumina DNA Prep (Nextera tagmentation)

## Pipeline Steps

### 1. Read quality control and trimming (`codes/filtering/`)
- FastQC for raw read quality assessment
- fastp and Trimmomatic for adapter trimming and quality filtering
- MultiQC for aggregated QC reports

### 2. Genome alignment (`codes/mapping/`)
- BWA-MEM alignment to GCF_000633615.1
- samtools for sorting, indexing, and statistics
- Group-level BAM merging (4 merged BAMs by experimental group)

### 3. Coverage analysis (`codes/analysis/`)
- samtools depth over ±600bp around sgRNA site
- Zone-based coverage: upstream 500bp, upstream 100bp, sgRNA site, downstream 100bp/500bp

### 4. Variant calling — GATK best practices (`codes/variant_calling/`)
- MarkDuplicates — PCR duplicate flagging (--REMOVE_DUPLICATES false)
- HaplotypeCaller — per-sample GVCF (-ERC GVCF, -ploidy 2, medium partition 7-day limit)
- GenomicsDBImport — joint sample database across 24 chromosomes
- GenotypeGVCFs — joint genotyping across all 15 samples
- VariantFiltration — hard filtering (SNPs: QD<2/FS>60/MQ<40; INDELs: QD<2/FS>200)
- Hotspot analysis — sliding window variant density, BH-corrected FDR < 0.05

### 5. CRISPR on-target analysis (`codes/CRISPResso/`)
- CRISPResso2 — per-sample editing quantification (60bp amplicon, both strands)
- CRISPResso2 — per-group analysis on merged BAMs (200bp amplicon)
- CRISPRessoCompare — 6 pairwise group comparisons
- CRISPRessoAggregate — multi-sample summary reports

### 6. Off-target analysis (`codes/CRISPResso/`)
- Cas-OFFinder — genome-wide off-target prediction (≤4 mismatches, 9 sites)
- CRISPOR — MIT/CFD scores and locus annotation (8 sites, guide ID: 326forw)
- Deduplication — combined 8 unique off-target sites (±2bp tolerance, CRISPOR priority)
- CRISPRessoWGS — editing quantification at 8 predicted off-target sites per sample

### 7. Population-specific genomes (`codes/assembly/`)

Two complementary Colombian-population reference genomes were built:

**A) Pseudogenome** (bcftools consensus):
- Applies Control-group variants (AF ≥ 0.667) to the Trinidad reference
- Near-identical coordinate system to reference (use chain file for exact liftover)
- Annotated with Liftoff v1.5.1: 99.5% transfer rate, 26,264 genes, 0 orphaned records
- See `reference/pseudogenome/README.md` for full details

**B) De novo scaffolded genome** (SPAdes + RagTag):
- True de novo co-assembly of 3 Control replicates → reference-guided scaffolding
- N50 = 28.3 Mb, BUSCO completeness 87.1%, genome fraction 82.8%
- Annotated with Liftoff v1.5.1: bdnf coverage=0.945, sequence_ID=0.923
- See `reference/colombian_scaffolded_genome/README.md` for full details

### 8. Visualization (`codes/assembly/`, `igv_files/`)
- IGV-ready files: pseudogenome FASTA + sorted/indexed Liftoff GFF3 + merged BAMs
- BED track with key CRISPR features: bdnf gene, sgRNA site, cut site, 8 off-target sites

## Repository Structure

```
codes/
├── filtering/          → trimming + FastQC scripts (fastp, trimmomatic, multiqc)
├── mapping/            → BWA-MEM alignment + BAM merging scripts
├── variant_calling/    → GATK pipeline (markdup → HaplotypeCaller → GenomicsDB
│                         → GenotypeGVCFs → VariantFiltration + hotspot analysis)
├── CRISPResso/         → CRISPResso2 on-target, off-target, WGS, compare scripts
│                         + Cas-OFFinder + CRISPOR conversion + combine_offtargets.py
├── analysis/           → coverage computation + hotspot density plots
└── assembly/           → Colombian genome pipeline:
                          pseudogenome (bcftools consensus + CrossMap + Liftoff)
                          de novo (SPAdes → RagTag → QUAST/BUSCO → Liftoff)
                          indexing (samtools, BWA, GATK, BLAST)
                          IGV preparation

igv_files/
├── features_of_interest.bed   → CRISPR features BED (bdnf, sgRNA, cut site, 8 OTs)
└── [genome + annotation]      → FASTA + Liftoff GFF3 (download separately — large files)

reference/
├── pseudogenome/              → Colombian pseudogenome + all indices + Liftoff annotation
│   └── README.md              → build method, QC, limitations, quick-start
└── colombian_scaffolded_genome/ → de novo scaffolded genome + all indices + annotation
    └── README.md              → build method, QC, limitations, quick-start
```

## Dependencies

| Tool | Version | Use |
|---|---|---|
| FastQC | 0.11.x | Read QC |
| fastp | 0.23.x | Trimming |
| Trimmomatic | 0.39 | Trimming |
| MultiQC | 1.x | QC aggregation |
| BWA | 0.7.17 | Alignment |
| samtools | 1.16.1 | BAM/FASTA processing |
| GATK | 4.4.0.0 | Variant calling |
| bedtools | 2.x | Interval operations |
| CRISPResso2 | 2.x | CRISPR editing analysis |
| Cas-OFFinder | 2.4 | Off-target prediction |
| bcftools | 1.15.1 | VCF processing + consensus |
| CrossMap | 0.7.3 | Coordinate liftover (chain file) |
| Liftoff | 1.5.1 | GFF3 annotation transfer |
| SPAdes | 4.0.0 | De novo genome assembly |
| RagTag | 2.1.0 | Reference-guided scaffolding |
| QUAST | 5.0.2 | Assembly QC |
| BUSCO | 5.7.1 | Assembly completeness (odb10) |
| BLAST | 2.14.1 | Local sequence search |
| Python | 3.9+ | Data processing |
| pandas / matplotlib / scipy | — | Analysis and visualization |

### Conda Environments

```bash
# System conda — trimming, QC, plotting
conda activate fastp_env        # fastp, multiqc, pandas, matplotlib

# Personal conda (miniconda3_crispresso)
conda activate crispresso2_env  # CRISPResso2, cas-offinder, pandas
conda activate liftoff_env      # Liftoff v1.5.1
conda activate crossmap_env     # CrossMap v0.7.3
```

## Reference Data

- **Reference genome:** GCF_000633615.1 (*Poecilia reticulata* Guanapo female 1.0 + MT)
- **NCBI:** https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000633615.1/
- **sgRNA target:** NC_024333.1 (LG3), position 15922039–15922058, negative strand
- **bdnf locus (pseudogenome):** NC_024333.1:15923726–15938393 (−)

## Key Genomic Landmarks (pseudogenome coordinates)

| Feature | Chromosome | Start | End | Strand |
|---|---|---|---|---|
| bdnf gene | NC_024333.1 | 15,923,726 | 15,938,393 | − |
| sgRNA target site | NC_024333.1 | 15,922,039 | 15,922,058 | − |
| Cas9 cut site | NC_024333.1 | ~15,922,041 | — | − |
| NLGN3-like (LOC103471143) | NC_024340.1 | 10,168,263 | 10,248,649 | − |

## Cluster

SLURM scheduler on the **hypatia** cluster (Universidad de los Andes).
Scripts use `--partition=short` (2 days) or `--partition=medium` (7 days)
depending on runtime requirements. GATK HaplotypeCaller requires medium partition.

## Key Results

- ~1.4M SNPs and ~370K INDELs per sample genome-wide, reflecting substantial
  divergence between the Colombian population and the Guanapo reference
- Consistent Ti/Tv ratio of 1.37 across all groups, confirming high-quality variant calls
- 403 merged hotspot regions of elevated variant density (FDR < 0.05, 1,780 windows)
- No CRISPR-induced INDELs detected at any of the 8 predicted off-target sites by GATK
- Colombian pseudogenome: 99.5% annotation transfer (Liftoff), bdnf confirmed present
- De novo scaffolded genome: N50 28.3 Mb, BUSCO 87.1% complete

## Citation

If you use this pipeline, please cite the relevant tools:

- **GATK:** Van der Auwera & O'Connor (2020). *Genomics in the Cloud*. O'Reilly.
- **CRISPResso2:** Clement et al. (2019). *Nature Biotechnology*.
- **Cas-OFFinder:** Bae et al. (2014). *Bioinformatics*.
- **BWA:** Li & Durbin (2009). *Bioinformatics*.
- **samtools:** Danecek et al. (2021). *GigaScience*.
- **SPAdes:** Bankevich et al. (2012). *Journal of Computational Biology*.
- **RagTag:** Alonge et al. (2022). *Genome Biology*.
- **Liftoff:** Shumate & Salzberg (2021). *Bioinformatics*.

## Author

Diego Andrés Martínez
Biomedical Engineer and Biologist
Universidad de los Andes, Colombia
da.martinez33@uniandes.edu.co | diegoandres3322@gmail.com

## License

MIT License — see LICENSE file for details.
