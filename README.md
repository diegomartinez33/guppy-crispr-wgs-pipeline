# Guppy CRISPR-Cas9 WGS Analysis Pipeline

Bioinformatics pipeline for whole-genome sequencing (WGS) analysis of CRISPR-Cas9 gene editing
and its off-target effects in *Poecilia reticulata* (guppy) from a Colombian population, with
proprietary population genomic resources and reusable tools for CRISPR guide and PCR primer
design in future genes.

## Overview

This repository contains all the code used to analyze Illumina NovaSeq X sequencing data from
guppy fish edited with CRISPR-Cas9 in the **bdnf** gene (brain-derived neurotrophic factor). The
pipeline covers: read trimming and QC, genomic alignment, variant calling, CRISPR editing
quantification, off-target analysis, population genomics, construction of two Colombian
population-specific reference genomes, and two reusable tools (CRISPR guide design and PCR
primer design) parameterized for use on any future gene.

The experimental fish are from a **Colombian population** of *P. reticulata*, sequenced and
analyzed against two public reference genomes:

- **v1 — Guanapo (Trinidad)**, female, short-read (GCF_000633615.1, 2014) — the genome used in
  all of this project's historical results. Marked as *suppressed* by NCBI in 2026.
- **v2 — male, PacBio+Hi-C** (GCF_904066995.2, 2025) — the species' current RefSeq genome, with
  far superior contiguity. Migration in progress (see [docs/RESULTS.md](docs/RESULTS.md#7-migration-to-the-v2-reference-genome--in-progress)).

A core component of this work is addressing the genomic divergence between the Colombian
population and these references in the context of CRISPR specificity — hence the two proprietary
population genomes (pseudogenome + de novo assembly, see pipeline step 7).

## Experimental Design

| Group | Description | n |
|---|---|---|
| Control | No CRISPR components | 3 |
| Only_MNP | Nanoparticles only, no Cas9 | 4 |
| Plasmid_Ko | Cas9 delivered via plasmid | 4 |
| RNP_Cas | Cas9 delivered as ribonucleoprotein | 4 |

- **Target gene:** bdnf (v1: `NC_024333.1`; v2: `NC_088832.1`)
- **sgRNA:** `TGAGAGACGCCCCGGGCATG` (negative strand)
- **Cas9:** SpCas9 NLS (NEB)
- **Sequencing:** Illumina NovaSeq X, 25B flow cell, 2×150bp paired-end
- **Library prep:** Illumina DNA Prep (Nextera tagmentation)

## Documentation

| Document | Contents |
|---|---|
| **[docs/PIPELINE.md](docs/PIPELINE.md)** | Each pipeline stage in detail: purpose, scripts, exact tool parameters, flow diagrams |
| **[docs/TUTORIAL.md](docs/TUTORIAL.md)** | How to run CRISPR guide design and primer design for a new gene (for colleagues) |
| **[docs/RESULTS.md](docs/RESULTS.md)** | Where each generated result lives, by objective, with status (complete/in progress/pending) |
| **[docs/CLUSTER_ACCESS.md](docs/CLUSTER_ACCESS.md)** | Shared cluster account + which heavy files were copied and where |
| [reference/pseudogenome/README.md](reference/pseudogenome/README.md) | Colombian pseudogenome: method, QC, limitations |
| [reference/colombian_scaffolded_genome/README.md](reference/colombian_scaffolded_genome/README.md) | Colombian de novo assembly: method, QC, limitations |

**Visual result reports** (Artifacts + downloadable HTML, see
[docs/RESULTS.md](docs/RESULTS.md#visual-reports--summary) for the full list): Guppy CRISPR
Atlas (guide design), plus off-target WGS, hotspots, population genomes, and primer design
reports under `analysis/reports/`.

## Quick Start

```bash
# 1. Clone and enter the repo
git clone <remote> off-target_data && cd off-target_data

# 2. Base genomes (reference/) and the CRISPOR container are NOT in git (too heavy) —
#    see docs/CLUSTER_ACCESS.md for where to get them (shared cluster account).

# 3. Choose a reference version (v1 by default, or v2) via REF_VERSION — see docs/TUTORIAL.md
REF_VERSION=v1 sbatch codes/mapping/bwa_index.sh

# 4. Design CRISPR guides for a new gene (needs no sample BAMs/VCFs):
module load minimap2 samtools/1.16.1 singularity/3.7.1
python3 codes/analysis/ko_guide_scan.py --gene my_gene --population pseudogenome

# 5. Design PCR primers for a gene's off-targets (needs its combined_offtargets.csv):
bash codes/analysis/setup_primer3.sh   # one-time setup
sbatch --export=ALL,GENE=my_gene,SITES_CSV=<path>,REF_VERSION=v1 \
  codes/analysis/run_offtarget_primer_design.sh
```
Full guide with all parameters: [docs/TUTORIAL.md](docs/TUTORIAL.md).

## Repository Structure

```
codes/
├── filtering/          → trimming + FastQC (fastp, trimmomatic, multiqc)
├── mapping/             → BWA-MEM alignment + per-group BAM merging
├── variant_calling/     → GATK pipeline (markdup → HaplotypeCaller → GenomicsDB
│                          → GenotypeGVCFs → VariantFiltration)
├── CRISPResso/          → CRISPResso2 on-target/off-target/WGS/compare + off-target
│                          discovery (Cas-OFFinder + CRISPOR)
├── analysis/             → hotspots, CRISPR KO/CRISPRi guide design (ko_guide_scan.py,
│                          crispri_tss_scan.py), primer design, consolidated reports
├── assembly/             → pseudogenome (bcftools consensus + CrossMap + Liftoff) +
│                          de novo assembly (SPAdes → RagTag → gap-filling/polishing → Liftoff)
└── genome_versions.sh    → shared v1/v2 config (REF_VERSION)

docs/                     → detailed documentation (see table above)

analysis/
├── ko_guide_scan/        → guide design results (8 genes) + report/ (Atlas HTML)
├── offtarget_primers/    → primer design results
└── reports/              → visual reports for the other sub-pipelines

reference/
├── pseudogenome/          → Colombian pseudogenome + its own README
└── colombian_scaffolded_genome/ → Colombian de novo assembly + its own README

igv_files/                → ready-to-load IGV Desktop bundle (v1)
```
See [docs/PIPELINE.md](docs/PIPELINE.md) for the full detail of each stage, and
[docs/CLUSTER_ACCESS.md](docs/CLUSTER_ACCESS.md) for which heavy directories live outside git
and where to find them.

## Dependencies

| Tool | Version | Use |
|---|---|---|
| FastQC | 0.11.x | Read QC |
| fastp | 0.23.x | Trimming (comparison) |
| Trimmomatic | 0.39 | Trimming (primary) |
| MultiQC | 1.x | QC aggregation |
| BWA | 0.7.17 | Alignment |
| samtools | 1.16.1 | BAM/FASTA processing |
| GATK | 4.4.0.0 | Variant calling |
| bedtools | 2.30.0 | Interval operations (hotspots) |
| CRISPResso2 | 2.x | CRISPR editing quantification |
| Cas-OFFinder | 2.4 | Off-target prediction |
| CRISPOR (Singularity) | v5.2c | MIT/CFD/Doench'16 scores, real off-targets |
| EMBOSS (`eprimer3`/`primersearch`) | 6.6.0 | Primer design and validation |
| primer3_core | 1.1.4 (legacy boulder-IO) | Engine required by `eprimer3` |
| bcftools | 1.15.1 | VCF processing + consensus |
| CrossMap | 0.7.3 | Coordinate liftover (chain file) |
| Liftoff | 1.5.1 | GFF3 annotation transfer |
| minimap2 | 2.24 | CDS/window alignment for variant classification |
| SPAdes | 4.0.0 | De novo assembly |
| RagTag | 2.1.0 | Reference-guided scaffolding |
| TGS-GapCloser | 1.2.1 | Gap-filling with Nanopore reads |
| NextPolish | 1.4.1 | Polishing with Illumina reads |
| QUAST | 5.0.2 | Structural assembly QC |
| BUSCO | 5.7.1 | Assembly completeness (`actinopterygii_odb10`) |
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
conda activate primer3_env      # primer3_core v1.1.4 (for eprimer3)
conda activate nextpolish_env   # NextPolish v1.4.1
conda activate tgsgapcloser_env # TGS-GapCloser v1.2.1
```

## Reference Data

- **v1:** GCF_000633615.1 — https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000633615.1/
- **v2:** GCF_904066995.2 — https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_904066995.2/
- **sgRNA site (v1):** `NC_024333.1:15922039-15922058`, negative strand
- **bdnf locus (v1 pseudogenome):** `NC_024333.1:15923726-15938393` (−)

## Cluster

SLURM scheduler on the **hypatia** cluster (Universidad de los Andes). Scripts use
`--partition=short` (2 days), `--partition=medium` (7 days), or `--partition=bigmem` depending on
the requirement — GATK HaplotypeCaller needs medium, SPAdes co-assembly needs bigmem. See
[docs/CLUSTER_ACCESS.md](docs/CLUSTER_ACCESS.md) for the shared account available to colleagues.

## Key Results

- ~1.4M SNPs and ~370K INDELs per sample genome-wide, reflecting substantial divergence between
  the Colombian population and the Guanapo reference
- No CRISPR-induced indel detected at any of the 8 off-target sites predicted by GATK; the only
  site with variants is a pre-existing polymorphism in the Control group
- 403 elevated variant-density hotspot regions (FDR<0.05)
- Colombian pseudogenome: 99.5% annotation transfer (Liftoff)
- Colombian de novo assembly: BUSCO 95.5% complete after gap-filling + polishing
- CRISPR guide design completed for 8 candidate genes (bdnf, agap3, grin1a/b, gria1a/b, gria2b,
  nlgn1) — see [Guppy CRISPR Atlas](analysis/ko_guide_scan/report/guppy_crispr_atlas.html)
- PCR primers designed and validated in silico for bdnf (9 on-/off-target sites)

See [docs/RESULTS.md](docs/RESULTS.md) for the full map with exact paths.

## Citation

If you use this pipeline, please cite the relevant tools:

- **GATK:** Van der Auwera & O'Connor (2020). *Genomics in the Cloud*. O'Reilly.
- **CRISPResso2:** Clement et al. (2019). *Nature Biotechnology*.
- **CRISPOR:** Haeussler et al. (2016). *Genome Biology*.
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

MIT License — see LICENSE file.
