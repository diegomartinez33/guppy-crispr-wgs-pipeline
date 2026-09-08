# Tutorial: correr los pipelines parametrizados para un gen/guía propio

Esta guía es para colegas que quieran usar los dos pipelines pensados para reutilización directa
sobre **genes nuevos**, sin tocar código: **diseño de guías CRISPR** (KO/CRISPRi) y **diseño de
primers PCR**. Ambos ya aceptan el gen y la versión de genoma como parámetro de línea de
comandos. Para el detalle técnico de cada paso ver [PIPELINE.md](PIPELINE.md); para dónde están
los resultados ya generados, [RESULTS.md](RESULTS.md); para acceso al clúster y a los archivos
base necesarios, [CLUSTER_ACCESS.md](CLUSTER_ACCESS.md).

## 0. Requisitos previos (una sola vez por usuario)

1. Acceso a la cuenta compartida del clúster — ver [CLUSTER_ACCESS.md](CLUSTER_ACCESS.md).
2. Los genomas base (`reference/`) y el contenedor de CRISPOR
   (`codes/analysis/crispor_singularity/`) deben existir en tu copia del repositorio — son los
   únicos archivos pesados que se copian a la cuenta compartida, justamente para que este
   tutorial funcione sin tener que re-generarlos.
3. Solo para diseño de primers: instalar `primer3_core` una vez (no viene con el módulo EMBOSS):
   ```bash
   bash codes/analysis/setup_primer3.sh
   ```
   Esto crea el entorno conda `primer3_env`. No hace falta repetirlo salvo que se borre.

## 1. Elegir versión de genoma de referencia

Todo el pipeline soporta dos genomas: `v1` (Trinidad/Guanapo, hembra, 2014 — usado en todos los
resultados históricos de este proyecto) y `v2` (macho, PacBio+Hi-C, 2025 — el genoma RefSeq
actual, migración en progreso). Se selecciona con la variable de entorno `REF_VERSION`:

```bash
# Por defecto (sin especificar) usa v1 — rutas idénticas a como siempre ha sido el proyecto
sbatch codes/mapping/bwa_index.sh

# Para usar v2 explícitamente:
REF_VERSION=v2 sbatch codes/mapping/bwa_index.sh
```

`REF_VERSION` se propaga automáticamente a todos los scripts que hacen
`source codes/genome_versions.sh` (la mayoría del pipeline — ver la columna "Dual-genoma" en
[PIPELINE.md](PIPELINE.md) para saber cuáles). No es necesario editar ningún script para cambiar
de versión.

## 2. Diseño de guías para un gen nuevo

El escaneo de guías (`ko_guide_scan.py` para CRISPRko, `crispri_tss_scan.py` para CRISPRi) solo
necesita el símbolo del gen — busca el gen por su nombre (`gene=NOMBRE;`) en el GFF de la
referencia y del pseudogenoma, así que debe existir en ambas anotaciones.

```bash
module load minimap2
module load samtools/1.16.1   # orden importa: minimap2 primero, samtools después
module load singularity/3.7.1

# CRISPRko — candidatos de corte en el CDS
python3 codes/analysis/ko_guide_scan.py --gene mi_gen --population pseudogenome

# CRISPRi — candidatos de silenciamiento cerca del TSS
python3 codes/analysis/crispri_tss_scan.py --gene mi_gen --population pseudogenome
```

Parámetros:
- `--gene` (obligatorio): símbolo del gen tal como aparece en el GFF (p.ej. `bdnf`, `grin1a`).
- `--population`: `pseudogenome` (recomendado, preserva estructura exón/intrón — ver
  [PIPELINE.md §7](PIPELINE.md#7-genomas-poblacionales-colombianos)), `scaffolded` (ensamblaje
  de novo), o `pseudogenome_v2` (una vez exista, ver [RESULTS.md](RESULTS.md)).
- `--no-crispor`: omite la puntuación adicional de CRISPOR (el escaneo manual de PAMs y la
  clasificación de variantes poblacionales siempre corren, con o sin esta bandera).

Para correr varios genes seguidos, edita el arreglo al inicio de
`codes/analysis/run_ko_guide_scan.sh`:
```bash
GENES=(bdnf agap3 grin1a grin1b gria1a gria1b gria2b nlgn1 mi_gen_nuevo)
```
y ejecútalo con `bash codes/analysis/run_ko_guide_scan.sh`.

**Salida:** `analysis/ko_guide_scan/<gen>_<population>_guide_comparison.csv` (CRISPRko) y
`_crispri_candidates.csv` (CRISPRi), más los TSV crudos de CRISPOR si estaba disponible.
Para regenerar el reporte visual consolidado (como el Guppy CRISPR Atlas) después de añadir un
gen nuevo: `python3 codes/analysis/build_guide_report.py` (edita primero la lista `GENES` al
inicio del script).

**Si tu gen nuevo necesita un genoma que no está registrado en CRISPOR** (poco común — los 3
genomas del proyecto, `guppyRefTrinidad`/`guppyColPseudogenome`/`guppyRefMaleV2`, ya cubren
v1/pseudogenoma/v2), regístralo primero con `codes/analysis/crispor_add_genomes.sh` (o
`crispor_add_genome_v2.sh` como referencia del patrón).

## 3. Diseño de primers para un gen nuevo

Requiere un CSV de sitios en el mismo formato que produce `combine_offtargets.py`
(columnas: `chromosome,start,end,strand,offtarget_seq,mismatches,mit_score,cfd_score,locus,source,found_by_both`).
Si tu gen ya pasó por el descubrimiento de off-targets ([PIPELINE.md §5](PIPELINE.md#5-descubrimiento-y-análisis-de-off-targets)),
ese CSV ya existe en `crispresso${OUT_SUFFIX}/offtargets/combined/combined_offtargets.csv`.

```bash
sbatch --export=ALL,GENE=mi_gen,SITES_CSV=/ruta/a/mi_gen_offtargets.csv,REF_VERSION=v1 \
  codes/analysis/run_offtarget_primer_design.sh
```

O directamente en Python (útil para ajustar parámetros finos sin editar el wrapper):
```bash
module load emboss/6.6.0 minimap2 samtools/1.16.1
export EMBOSS_PRIMER3_CORE=<ruta_a>/envs/primer3_env/bin/primer3_core

python3 codes/analysis/design_offtarget_primers.py \
  --gene mi_gen --sites-csv mi_gen_offtargets.csv --ref-version v1 \
  --window-halfsize 500 --exclude-halfsize 75 \
  --product-min 200 --product-max 400 \
  --opt-tm 60 --min-tm 58 --max-tm 62 --min-gc 40 --max-gc 60 \
  --num-return 5 --specificity-mismatch-pct 10
```

Parámetros más útiles para ajustar:
- `--window-halfsize`: cuánta secuencia alrededor del sitio de corte se le da a `eprimer3` para
  buscar (default 500pb a cada lado).
- `--exclude-halfsize`: región alrededor del corte donde NINGÚN primer puede caer (default
  75pb) — evita diseñar un primer justo sobre el sitio donde se espera el indel.
- `--product-min`/`--product-max`: tamaño de amplicón deseado (default 200-400pb, buen rango
  para gel + secuenciación Illumina/Sanger).
- `--no-population-check`: omite el chequeo contra el pseudogenoma (más rápido, pero no detecta
  si un primer cae sobre una variante poblacional real).

**Salida:** `analysis/offtarget_primers/<gen>_<ref_version>_primers.csv` — una fila por par
candidato, con columnas de especificidad (`primersearch`) y variantes poblacionales en el sitio
de unión del primer. Revisa siempre la columna `design_status`: `OK` = candidato real,
`INPUT_ERROR_AMBIGUOUS_BASES` = el sitio tiene códigos IUPAC ambiguos en la referencia y no se
pudo diseñar nada (ver limitación en [PIPELINE.md §9](PIPELINE.md#9-diseño-de-primers-pcr)),
`NO_CANDIDATES_FOUND` = se buscó pero ningún candidato cumplió los parámetros (prueba relajando
`--min-tm`/`--max-tm`/`--min-gc`/`--max-gc`).

## 4. Preguntas frecuentes

**¿Puedo correr esto sobre una guía CRISPR que yo elija, en vez de un gen completo?** Los dos
pipelines trabajan a nivel de gen (buscan candidatos ellos mismos dentro del CDS/ventana TSS).
Si ya tienes una guía específica seleccionada y solo quieres sus off-targets/primers, usa
directamente `casoffinder.sh` + `crispor_offtarget_scan.sh` con tu secuencia de guía
(ver [PIPELINE.md §5](PIPELINE.md#5-descubrimiento-y-análisis-de-off-targets)) para generar tu
propio `combined_offtargets.csv`, y de ahí en adelante sigue el paso 3 de este tutorial.

**¿Qué pasa si mi gen no aparece en el GFF de la referencia?** `ko_guide_scan.py`/
`crispri_tss_scan.py` fallan con un error claro (`ERROR: gene 'X' not found in reference GFF`).
Verifica el símbolo exacto con `grep "gene=" reference/GCF_*_annotation.gff | grep -i mi_gen`.

**¿Necesito correr todo el pipeline de WGS (mapeo, GATK) para usar estos dos scripts?** No —
ambos solo necesitan los archivos de `reference/` (genomas + anotaciones + pseudogenoma) y,
para guías, el contenedor CRISPOR. No dependen de BAMs ni VCFs de ninguna muestra.
