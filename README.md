# Guppy CRISPR-Cas9 WGS Analysis Pipeline

Pipeline bioinformático para el análisis de secuenciación de genoma completo (WGS) de edición
génica CRISPR-Cas9 y sus efectos off-target en *Poecilia reticulata* (guppy) de una población
colombiana, con recursos genómicos poblacionales propios y herramientas reutilizables para el
diseño de guías CRISPR y primers PCR en genes futuros.

## Overview

Este repositorio contiene todo el código usado para analizar datos de secuenciación Illumina
NovaSeq X de peces guppy editados con CRISPR-Cas9 en el gen **bdnf** (brain-derived neurotrophic
factor). El pipeline cubre: recorte y QC de lecturas, alineamiento genómico, llamado de
variantes, cuantificación de edición CRISPR, análisis de off-targets, genómica poblacional,
construcción de dos genomas de referencia específicos para la población colombiana, y dos
herramientas reutilizables (diseño de guías CRISPR y diseño de primers PCR) parametrizadas para
aplicarse a cualquier gen futuro.

Los peces experimentales son de una **población colombiana** de *P. reticulata*, secuenciados y
analizados contra dos genomas de referencia públicos:

- **v1 — Guanapo (Trinidad)**, hembra, short-read (GCF_000633615.1, 2014) — el genoma usado en
  todos los resultados históricos de este proyecto. Marcado como *suppressed* por NCBI en 2026.
- **v2 — macho, PacBio+Hi-C** (GCF_904066995.2, 2025) — el genoma RefSeq actual de la especie,
  con contigüidad muy superior. Migración en progreso (ver [docs/RESULTS.md](docs/RESULTS.md#7-migración-al-genoma-de-referencia-v2--en-progreso)).

Un componente central de este trabajo es abordar la divergencia genómica entre la población
colombiana y estas referencias en el contexto de la especificidad de CRISPR — de ahí los dos
genomas poblacionales propios (pseudogenoma + ensamblaje de novo, ver punto 7 del pipeline).

## Diseño experimental

| Grupo | Descripción | n |
|---|---|---|
| Control | Sin componentes CRISPR | 3 |
| Only_MNP | Solo nanopartículas, sin Cas9 | 4 |
| Plasmid_Ko | Cas9 entregado vía plásmido | 4 |
| RNP_Cas | Cas9 entregado como ribonucleoproteína | 4 |

- **Gen objetivo:** bdnf (v1: `NC_024333.1`; v2: `NC_088832.1`)
- **sgRNA:** `TGAGAGACGCCCCGGGCATG` (hebra negativa)
- **Cas9:** SpCas9 NLS (NEB)
- **Secuenciación:** Illumina NovaSeq X, celda 25B, 2×150pb pareado
- **Preparación de librería:** Illumina DNA Prep (tagmentación Nextera)

## Documentación

| Documento | Contenido |
|---|---|
| **[docs/PIPELINE.md](docs/PIPELINE.md)** | Cada etapa del pipeline en detalle: objetivo, scripts, parámetros exactos de cada herramienta, diagramas de flujo |
| **[docs/TUTORIAL.md](docs/TUTORIAL.md)** | Cómo correr el diseño de guías CRISPR y el diseño de primers para un gen nuevo (colegas) |
| **[docs/RESULTS.md](docs/RESULTS.md)** | Dónde está cada resultado ya generado, por objetivo, con estado (completo/en progreso/pendiente) |
| **[docs/CLUSTER_ACCESS.md](docs/CLUSTER_ACCESS.md)** | Cuenta compartida del clúster + qué archivos pesados se copiaron y dónde |
| [reference/pseudogenome/README.md](reference/pseudogenome/README.md) | Pseudogenoma colombiano: método, QC, limitaciones |
| [reference/colombian_scaffolded_genome/README.md](reference/colombian_scaffolded_genome/README.md) | Ensamblaje de novo colombiano: método, QC, limitaciones |

**Reportes visuales de resultados** (Artifacts + HTML descargable, ver
[docs/RESULTS.md](docs/RESULTS.md#reportes-visuales--resumen) para la lista completa): Guppy
CRISPR Atlas (diseño de guías), y reportes de off-targets WGS, hotspots, genomas poblacionales,
y diseño de primers en `analysis/reports/`.

## Quick start

```bash
# 1. Clonar y ubicarse en el repo
git clone <remote> off-target_data && cd off-target_data

# 2. Los genomas base (reference/) y el contenedor CRISPOR NO están en git (son pesados) —
#    ver docs/CLUSTER_ACCESS.md para dónde obtenerlos (cuenta compartida del clúster).

# 3. Elegir versión de referencia (v1 por defecto, o v2) con REF_VERSION — ver docs/TUTORIAL.md
REF_VERSION=v1 sbatch codes/mapping/bwa_index.sh

# 4. Diseñar guías CRISPR para un gen nuevo (no requiere BAMs/VCFs de ninguna muestra):
module load minimap2 samtools/1.16.1 singularity/3.7.1
python3 codes/analysis/ko_guide_scan.py --gene mi_gen --population pseudogenome

# 5. Diseñar primers PCR para los off-targets de un gen (requiere su combined_offtargets.csv):
bash codes/analysis/setup_primer3.sh   # una sola vez
sbatch --export=ALL,GENE=mi_gen,SITES_CSV=<ruta>,REF_VERSION=v1 \
  codes/analysis/run_offtarget_primer_design.sh
```
Guía completa con todos los parámetros: [docs/TUTORIAL.md](docs/TUTORIAL.md).

## Estructura del repositorio

```
codes/
├── filtering/          → recorte + FastQC (fastp, trimmomatic, multiqc)
├── mapping/             → alineamiento BWA-MEM + fusión de BAMs por grupo
├── variant_calling/     → pipeline GATK (markdup → HaplotypeCaller → GenomicsDB
│                          → GenotypeGVCFs → VariantFiltration)
├── CRISPResso/          → CRISPResso2 on-target/off-target/WGS/compare + descubrimiento
│                          de off-targets (Cas-OFFinder + CRISPOR)
├── analysis/             → hotspots, guías CRISPR KO/CRISPRi (ko_guide_scan.py,
│                          crispri_tss_scan.py), diseño de primers, reportes consolidados
├── assembly/             → pseudogenoma (bcftools consensus + CrossMap + Liftoff) +
│                          ensamblaje de novo (SPAdes → RagTag → gap-filling/pulido → Liftoff)
└── genome_versions.sh    → config compartida v1/v2 (REF_VERSION)

docs/                     → documentación detallada (ver tabla arriba)

analysis/
├── ko_guide_scan/        → resultados de diseño de guías (8 genes) + report/ (Atlas HTML)
├── offtarget_primers/    → resultados de diseño de primers
└── reports/              → reportes visuales de las demás sub-pipelines

reference/
├── pseudogenome/          → genoma pseudogenoma colombiano + README propio
└── colombian_scaffolded_genome/ → ensamblaje de novo colombiano + README propio

igv_files/                → paquete listo para IGV Desktop (v1)
```
Ver [docs/PIPELINE.md](docs/PIPELINE.md) para el detalle completo de cada etapa, y
[docs/CLUSTER_ACCESS.md](docs/CLUSTER_ACCESS.md) para qué directorios pesados existen fuera de
git y dónde encontrarlos.

## Dependencias

| Herramienta | Versión | Uso |
|---|---|---|
| FastQC | 0.11.x | QC de lecturas |
| fastp | 0.23.x | Recorte (comparación) |
| Trimmomatic | 0.39 | Recorte (principal) |
| MultiQC | 1.x | Agregación de QC |
| BWA | 0.7.17 | Alineamiento |
| samtools | 1.16.1 | Procesamiento BAM/FASTA |
| GATK | 4.4.0.0 | Llamado de variantes |
| bedtools | 2.30.0 | Operaciones de intervalos (hotspots) |
| CRISPResso2 | 2.x | Cuantificación de edición CRISPR |
| Cas-OFFinder | 2.4 | Predicción de off-targets |
| CRISPOR (Singularity) | v5.2c | Puntajes MIT/CFD/Doench'16, off-targets reales |
| EMBOSS (`eprimer3`/`primersearch`) | 6.6.0 | Diseño y validación de primers |
| primer3_core | 1.1.4 (legado boulder-IO) | Motor requerido por `eprimer3` |
| bcftools | 1.15.1 | Procesamiento VCF + consenso |
| CrossMap | 0.7.3 | Liftover de coordenadas (chain file) |
| Liftoff | 1.5.1 | Transferencia de anotación GFF3 |
| minimap2 | 2.24 | Alineamiento CDS/ventanas para clasificación de variantes |
| SPAdes | 4.0.0 | Ensamblaje de novo |
| RagTag | 2.1.0 | Scaffolding guiado por referencia |
| TGS-GapCloser | 1.2.1 | Relleno de gaps con lecturas Nanopore |
| NextPolish | 1.4.1 | Pulido con lecturas Illumina |
| QUAST | 5.0.2 | QC estructural de ensamblaje |
| BUSCO | 5.7.1 | Completitud de ensamblaje (`actinopterygii_odb10`) |
| BLAST | 2.14.1 | Búsqueda de secuencias local |
| Python | 3.9+ | Procesamiento de datos |
| pandas / matplotlib / scipy | — | Análisis y visualización |

### Entornos conda

```bash
# Conda del sistema — recorte, QC, gráficas
conda activate fastp_env        # fastp, multiqc, pandas, matplotlib

# Conda personal (miniconda3_crispresso)
conda activate crispresso2_env  # CRISPResso2, cas-offinder, pandas
conda activate liftoff_env      # Liftoff v1.5.1
conda activate crossmap_env     # CrossMap v0.7.3
conda activate primer3_env      # primer3_core v1.1.4 (para eprimer3)
conda activate nextpolish_env   # NextPolish v1.4.1
conda activate tgsgapcloser_env # TGS-GapCloser v1.2.1
```

## Datos de referencia

- **v1:** GCF_000633615.1 — https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000633615.1/
- **v2:** GCF_904066995.2 — https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_904066995.2/
- **Sitio de la sgRNA (v1):** `NC_024333.1:15922039-15922058`, hebra negativa
- **Locus bdnf (pseudogenoma v1):** `NC_024333.1:15923726-15938393` (−)

## Clúster

Scheduler SLURM en el clúster **hypatia** (Universidad de los Andes). Los scripts usan
`--partition=short` (2 días), `--partition=medium` (7 días) o `--partition=bigmem` según el
requerimiento — GATK HaplotypeCaller necesita medium, SPAdes co-assembly necesita bigmem. Ver
[docs/CLUSTER_ACCESS.md](docs/CLUSTER_ACCESS.md) para la cuenta compartida disponible para
colegas.

## Resultados clave

- ~1.4M SNPs y ~370K INDELs por muestra genoma completo, reflejando divergencia sustancial entre
  la población colombiana y la referencia Guanapo
- Ningún indel inducido por CRISPR detectado en ninguno de los 8 sitios off-target predichos por
  GATK; el único sitio con variantes es polimorfismo preexistente en el grupo Control
- 403 regiones hotspot de densidad elevada de variantes (FDR<0.05)
- Pseudogenoma colombiano: 99.5% de transferencia de anotación (Liftoff)
- Ensamblaje de novo colombiano: BUSCO 95.5% completo tras gap-filling + pulido
- Diseño de guías CRISPR completado para 8 genes candidatos (bdnf, agap3, grin1a/b, gria1a/b,
  gria2b, nlgn1) — ver [Guppy CRISPR Atlas](analysis/ko_guide_scan/report/guppy_crispr_atlas.html)
- Primers PCR diseñados y validados in-silico para bdnf (9 sitios on-/off-target)

Ver [docs/RESULTS.md](docs/RESULTS.md) para el mapa completo con rutas exactas.

## Citación

Si usas este pipeline, por favor cita las herramientas relevantes:

- **GATK:** Van der Auwera & O'Connor (2020). *Genomics in the Cloud*. O'Reilly.
- **CRISPResso2:** Clement et al. (2019). *Nature Biotechnology*.
- **CRISPOR:** Haeussler et al. (2016). *Genome Biology*.
- **Cas-OFFinder:** Bae et al. (2014). *Bioinformatics*.
- **BWA:** Li & Durbin (2009). *Bioinformatics*.
- **samtools:** Danecek et al. (2021). *GigaScience*.
- **SPAdes:** Bankevich et al. (2012). *Journal of Computational Biology*.
- **RagTag:** Alonge et al. (2022). *Genome Biology*.
- **Liftoff:** Shumate & Salzberg (2021). *Bioinformatics*.

## Autor

Diego Andrés Martínez
Ingeniero Biomédico y Biólogo
Universidad de los Andes, Colombia
da.martinez33@uniandes.edu.co | diegoandres3322@gmail.com

## Licencia

MIT License — ver archivo LICENSE.
