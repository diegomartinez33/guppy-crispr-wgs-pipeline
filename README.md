# Guppy CRISPR-Cas9 WGS Analysis Pipeline

Bioinformatics pipeline for whole-genome sequencing (WGS) analysis 
of CRISPR-Cas9 editing efficiency and off-target effects in 
*Poecilia reticulata* (guppy fish) from a Colombian population.

## Overview

This repository contains all scripts and code used to analyze 
Illumina NovaSeq X whole-genome sequencing data from guppy fish 
subjected to CRISPR-Cas9 editing of the **bdnf** (brain-derived 
neurotrophic factor) gene. The pipeline spans read trimming and 
quality control, genome alignment, variant calling, CRISPR editing 
quantification, and population genomics analysis.

The experimental fish are from a **Colombian population** of 
*P. reticulata*, sequenced and analyzed against the publicly 
available **Guanapo (Trinidad) reference genome** 
(GCF_000633615.1). A key component of this work is addressing the 
genomic divergence between these two populations (~1.4M SNPs 
genome-wide) in the context of CRISPR specificity assessment.

## Experimental Design

| Group | Description | n |
|---|---|---|
| Control | No CRISPR components delivered | 3 |
| Only_MNP | Nanoparticles only, no Cas9 | 4 |
| Plasmid_Ko | Cas9 delivered via plasmid | 4 |
| RNP_Cas | Cas9 delivered as ribonucleoprotein | 4 |

- **Target gene:** bdnf
- **sgRNA:** TGAGAGACGCCCCGGGCATG (negative strand)
- **Cas9:** SpCas9 NLS (NEB)
- **Sequencing:** Illumina NovaSeq X, 25B flowcell, 2×150bp PE
- **Library prep:** Illumina DNA Prep (Nextera tagmentation)

## Pipeline Steps

### 1. Read quality control and trimming
- FastQC for raw read quality assessment
- fastp and Trimmomatic for adapter trimming and quality filtering
- MultiQC for aggregated QC reports

### 2. Genome alignment
- BWA-MEM alignment to GCF_000633615.1
- samtools for sorting, indexing, and flagstat statistics

### 3. Variant calling (GATK best practices)
- MarkDuplicates — PCR duplicate flagging
- HaplotypeCaller — per-sample GVCF generation
- GenomicsDBImport — joint sample database
- GenotypeGVCFs — joint genotyping across all 15 samples
- VariantFiltration — hard filtering (SNPs and INDELs separately)

### 4. CRISPR on-target analysis (CRISPResso2)
- Merged BAM files by experimental group
- CRISPResso2 — editing efficiency quantification at bdnf cut site
- CRISPRessoCompare — pairwise group comparisons
- CRISPRessoAggregate — multi-sample summary reports

### 5. Off-target analysis
- Cas-OFFinder — genome-wide off-target site prediction (≤4 mismatches)
- CRISPOR — MIT score and locus annotation for predicted sites
- CRISPRessoWGS — editing quantification at 8 predicted off-target sites
- GATK variant inspection at predicted off-target coordinates

### 6. Population genomics
- SNP density analysis — divergence hotspots between Colombian 
  population and Guanapo reference
- Nucleotide diversity (π) and Tajima's D
- Colombian pseudogenome construction via 
  FastaAlternateReferenceMaker
- Re-evaluation of off-target predictions on population-specific 
  reference

## Repository Structure
```
├── trimming/
│   ├── run_fastp.sh
│   └── run_trimmomatic.sh
├── mapping/
│   ├── run_bwa.sh
│   └── merge_bams_by_group.sh
├── gatk/
│   ├── markduplicates.sh
│   ├── haplotypecaller.sh
│   ├── reindex_gvcf.sh
│   ├── genomicsdb_import.sh
│   ├── genotype_gvcfs.sh
│   ├── select_variants.sh
│   ├── variant_filtration.sh
│   └── select_offtargets.sh
├── crispresso/
│   ├── crispresso_ontarget_merged.sh
│   ├── crispresso_compare.sh
│   ├── crispresso_aggregate.sh
│   ├── crispresso_wgs_offtarget.sh
│   └── combine_offtargets.py
├── casoffinder/
│   ├── run_casoffinder.sh
│   └── convert_crispor.py
├── coverage/
│   ├── run_coverage.sh
│   └── plot_coverage.py
├── population_genomics/
│   ├── snp_density.sh
│   ├── diversity_metrics.sh
│   ├── fst_analysis.sh
│   ├── pseudogenome.sh
│   └── manhattan_plot.py
└── CLAUDE.md
```
## Dependencies

### Tools
```
| Tool | Version | Use |
|---|---|---|
| FastQC | 0.11.9 | Read QC |
| fastp | 0.23.x | Trimming |
| Trimmomatic | 0.39 | Trimming |
| MultiQC | 1.x | QC aggregation |
| BWA | 0.7.17 | Alignment |
| samtools | 1.16.1 | BAM processing |
| GATK | 4.4.0.0 | Variant calling |
| CRISPResso2 | 2.x | CRISPR analysis |
| Cas-OFFinder | 2.4 | Off-target prediction |
| vcftools | 0.1.16 | Population genetics |
| bcftools | 1.x | VCF processing |
| Python | 3.9 | Data processing |
| pandas | 1.x | Data analysis |
| matplotlib | 3.x | Visualization |
```

### Conda Environments
```bash
# System environment (trimming, alignment, QC)
conda create -n fastp_env fastp multiqc pandas matplotlib

# CRISPResso2 environment
conda create -n crispresso2_env -c bioconda crispresso2 cas-offinder
```

## Reference Data

- **Reference genome:** GCF_000633615.1 
  (*Poecilia reticulata* Guanapo female 1.0 + MT)
- **Download:** 
  https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000633615.1/
- **sgRNA target:** NC_024333.1 (LG3), position 15922039–15922058, 
  negative strand

## Cluster Configuration

All scripts are written for SLURM job scheduler. Key parameters:

```bash
#SBATCH --partition=short
#SBATCH --mail-type=ALL
```

Adjust `--mem`, `--cpus-per-task`, and `--time` according to your 
cluster's available resources and policies.

## Key Results

- ~1.4M SNPs and ~370K INDELs per sample genome-wide, reflecting 
  substantial divergence between the Colombian experimental 
  population and the Guanapo reference genome
- Consistent Ti/Tv ratio of 1.37 across all experimental groups, 
  confirming high-quality variant calling
- No CRISPR-induced INDELs detected at any of the 8 predicted 
  off-target sites by GATK; two background SNPs at OT4 
  (NC_024332.1) consistent with natural population polymorphism
- On-target editing analysis performed with CRISPResso2 using 
  merged BAMs per group (60bp amplicon, both strands)

## Citation

If you use this pipeline, please cite the tools it depends on:

- **GATK:** Van der Auwera & O'Connor (2020). Genomics in the Cloud. O'Reilly.
- **CRISPResso2:** Clement et al. (2019). Nature Biotechnology.
- **Cas-OFFinder:** Bae et al. (2014). Bioinformatics.
- **BWA:** Li & Durbin (2009). Bioinformatics.
- **samtools:** Danecek et al. (2021). GigaScience.

## Author

Diego Andrés Martínez  
Biomedical engineer and Biologist  
Universidad de los Andes, Colombia  
da.martinez33@uniandes.edu.co

## License

MIT License — see LICENSE file for details.