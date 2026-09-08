# Mapa de resultados

Dónde está cada resultado ya generado, organizado por objetivo. ✅ = completo, 🔄 = en progreso,
⏳ = pendiente. Para el detalle de cómo se generó cada cosa ver [PIPELINE.md](PIPELINE.md).

## 1. Análisis de off-targets WGS (GATK + CRISPResso2) — ✅ completo (v1)

| Qué | Ruta |
|---|---|
| Resumen de variantes por muestra (conteos, Ti/Tv) | `codes/analysis/gatk_summary/gatk_variant_summary.csv` + `gatk_summary_barplot.png`, `gatk_titv_boxplot.png` |
| Genotipos en los 8 sitios off-target (solo OT4 tiene variantes, preexistentes en Control) | `codes/analysis/gatk_summary/offtarget_genotypes.csv` + `offtarget_genotype_heatmap.png` |
| % de edición on-target vs off-target por muestra/grupo | `codes/analysis/editing_comparison/editing_summary.csv` + `editing_heatmap.png`, `ontarget_barplot.png`, `offtarget_dotplot.png` |
| Lista final de los 8 off-targets (coordenadas, MIT/CFD, locus) | `crispresso/offtargets/combined/combined_offtargets.csv` + `combined_offtargets_igv.bed` |
| VCFs filtrados finales (SNP/INDEL, todo el genoma) | `gatk/trimmomatic/vcf_filtered/{snps,indels}_filtered.vcf.gz` |
| VCF restringido a los 8 loci off-target | `gatk/trimmomatic/vcf_offtargets/offtarget_variants.vcf.gz` |
| Reportes CRISPResso2 por muestra (on-target, 15 carpetas + `merged/`) | `crispresso/ontarget/trimmomatic/<SAMPLE>/CRISPResso_on_<SAMPLE>/` |
| Reportes CRISPRessoWGS por muestra (8 sitios off-target, 15 carpetas + `aggregate/`) | `crispresso/wgs/trimmomatic/<SAMPLE>/` |
| **Reporte visual consolidado** | [`analysis/reports/offtarget_wgs_report.html`](../analysis/reports/offtarget_wgs_report.html) |

**Hallazgo principal:** ningún indel inducido por CRISPR detectado en los 8 sitios off-target
por GATK; el único sitio con variantes (OT4) ya las tenía el grupo Control, es polimorfismo
poblacional preexistente, no un efecto de la edición.

## 2. Hotspots de variantes — ✅ completo, solo v1

| Qué | Ruta |
|---|---|
| Densidad de variantes por ventana (10kb/2kb) | `gatk/trimmomatic/hotspots/window_counts_annotated.csv`, `window_counts.tsv` |
| Regiones hotspot finales (403 fusionadas, FDR<0.05) | `gatk/trimmomatic/hotspots/hotspots.bed` |
| Anotación génica de hotspots (incl. ortólogos zebrafish, enriquecimiento gProfiler) | `hotspot_gene_summary.tsv`, `hotspot_gene_overlaps.tsv`, `hotspot_genes_zebrafish.txt`, `gProfiler_*.csv` |
| Manhattan plot genoma completo + densidad por cromosoma | `hotspot_manhattan_genome.png`, `hotspot_plots/` |
| 5 figuras resumen finales | `gatk/trimmomatic/hotspots/summary_plots/hotspot_summary_0{1..5}_*.png` |
| **Reporte visual consolidado** | [`analysis/reports/hotspots_report.html`](../analysis/reports/hotspots_report.html) |

**No corrido aún bajo v2** — `gatk/trimmomatic_v2/` no tiene subcarpeta `hotspots/` (ver
objetivo 7).

## 3. Pseudogenoma colombiano — ✅ completo (v1)

`reference/pseudogenome/` — genoma + todos los índices (`.fai`, `.dict`, BWA, minimap2) +
anotación Liftoff (99.5% transferencia, 26,264 genes, 0 huérfanos) + `.chain` para liftover
exacto. Detalle completo, método y limitaciones:
**[reference/pseudogenome/README.md](../reference/pseudogenome/README.md)**.

## 4. Ensamblaje de novo colombiano — ✅ completo (v1)

`reference/colombian_scaffolded_genome/` — genoma final (gap-filled + pulido) + índices +
anotación Liftoff transferida (bdnf: coverage=0.945, sequence_ID=0.923). QC comparativo en las 4
etapas (raw → polished → gapfilled → gapfilled+polished) en `assembly/qc_results/`. Detalle
completo, método y limitaciones:
**[reference/colombian_scaffolded_genome/README.md](../reference/colombian_scaffolded_genome/README.md)**.
Reporte visual comparativo de ambos genomas poblacionales (pseudogenoma + ensamblaje de novo):
[`analysis/reports/genome_resources_report.html`](../analysis/reports/genome_resources_report.html).

## 5. Diseño de guías CRISPR KO/CRISPRi (8 genes) — ✅ completo

`analysis/ko_guide_scan/` — un set de archivos por gen (bdnf, agap3, grin1a, grin1b, gria1a,
gria1b, gria2b, nlgn1): comparación de candidatos CRISPRko (`*_guide_comparison.csv`),
candidatos CRISPRi (`*_crispri_candidates.csv`), variantes de cuerpo completo del gen
(`*_gene_body_variants.csv`), y TSVs crudos de CRISPOR donde estuvo disponible.

**Reporte consolidado (Guppy CRISPR Atlas):**
[`analysis/ko_guide_scan/report/guppy_crispr_atlas.html`](../analysis/ko_guide_scan/report/guppy_crispr_atlas.html)
— también publicado como Artifact (enlace compartible bajo pedido).

**Limitación:** agap3, grin1a, gria1a no tienen puntajes de CRISPOR (códigos IUPAC ambiguos en
la referencia v1 hacen crashear `crispor.py`) — se repetirá contra v2 cuando esté disponible.
Detalle en [PIPELINE.md §8](PIPELINE.md#8-diseño-de-guías-crispr-ko--crispri-por-gen).

## 6. Diseño de primers PCR — 🔄 parcial (solo bdnf/v1)

| Qué | Ruta |
|---|---|
| Tabla final de primers (9 sitios: on-target + 8 off-target) | `analysis/offtarget_primers/bdnf_v1_primers.csv` |
| Salidas crudas de `eprimer3` por sitio | `analysis/offtarget_primers/raw/bdnf_v1/` |
| **Reporte visual** | [`analysis/reports/primer_design_report.html`](../analysis/reports/primer_design_report.html) |

**Pendiente:** los otros 7 genes candidatos (agap3, grin1a, grin1b, gria1a, gria1b, gria2b,
nlgn1) y la versión v2 — el script ya está listo y parametrizado (`--gene`/`--ref-version`), solo
falta correrlo (ver [TUTORIAL.md §3](TUTORIAL.md#3-diseño-de-primers-para-un-gen-nuevo)).

**Verificación de primers RT-qPCR existentes (2026-09-08)** — consulta puntual, no parte del
pipeline automático: 4 pares de primers ya usados en el laboratorio (bdnf + housekeeping
miosina/beta-actina/rpl13a) verificados contra v1, v2 y el pseudogenoma colombiano. Hallazgo
principal: **SNP colombiano real en `rpl_13a_F`** (confirmado también contra v2, que coincide con
el alelo original del primer) y **`miosina_guppy_F/R` sin sitio de unión identificable** en el
genoma de guppy (posible primer diseñado para otra especie). Detalle completo en la sección
"RT-qPCR" del [reporte de primers](../analysis/reports/primer_design_report.html) y en
`CLAUDE.md` ("RT-qPCR Primer Verification").

## 7. Migración al genoma de referencia v2 — 🔄 en progreso

| Etapa | Estado |
|---|---|
| Descubrimiento de off-targets (Cas-OFFinder + CRISPOR) para bdnf | ✅ hecho, cruzado (8/8 coinciden) — `crispresso_v2/offtargets/combined/combined_offtargets.csv` |
| Mapeo (BWA) + MarkDuplicates + HaplotypeCaller (15 muestras) | 🔄 en progreso — `gatk/trimmomatic_v2/{markdup,gvcf}/` existen; HaplotypeCaller corriendo |
| GenomicsDBImport / GenotypeGVCFs / VariantFiltration | ⏳ pendiente (dependencia de lo anterior) |
| Pseudogenoma v2 | ⏳ pendiente (necesita el VCF filtrado de arriba) |
| Ensamblaje de novo re-escalado contra v2 (Fase 2 de RagTag) | ⏳ pendiente (deliberadamente pospuesto hasta que termine lo anterior) |
| CRISPResso on-target/off-target bajo v2 | ⏳ pendiente — no existe `crispresso_v2/ontarget/` ni `wgs/` todavía |
| Hotspots bajo v2 | ⏳ pendiente |
| Diseño de guías/primers para bdnf ya soporta `--ref-version v2`/`--population pseudogenome_v2` | ✅ código listo, esperando el pseudogenoma v2 |

## 8. Archivos IGV — ✅ completo, solo v1

`igv_files/` — genoma pseudogenoma + anotación (bgzip+tabix) + BAMs fusionados por grupo (4) +
`features_of_interest.bed` (bdnf, sitio de guía, sitio de corte, 8 off-targets). Listo para
cargar directamente en IGV Desktop. Sin equivalente v2 todavía (consistente con el objetivo 7).

---

## Reportes visuales — resumen

Cada reporte existe en dos formas: publicado como Artifact (enlace compartible, privado por
defecto — compártelo desde el menú de la propia página) y como copia HTML autocontenida en el
repositorio (para enviar directamente a un colega sin necesidad de cuenta en claude.ai).

| Reporte | Objetivo | Artifact | Copia local |
|---|---|---|---|
| Guppy CRISPR Atlas | Diseño de guías KO/CRISPRi, 8 genes | (compartir bajo pedido) | [`analysis/ko_guide_scan/report/guppy_crispr_atlas.html`](../analysis/ko_guide_scan/report/guppy_crispr_atlas.html) |
| Off-Target WGS Report | Off-targets WGS (GATK + CRISPResso) | [enlace](https://claude.ai/code/artifact/3290263b-0f56-4d76-84fe-825d4d98110c) | [`analysis/reports/offtarget_wgs_report.html`](../analysis/reports/offtarget_wgs_report.html) |
| Variant Hotspots Report | Hotspots de variantes | [enlace](https://claude.ai/code/artifact/bd831a9e-276a-4aef-a8ef-d35be58f539c) | [`analysis/reports/hotspots_report.html`](../analysis/reports/hotspots_report.html) |
| Colombian Genome Resources | Pseudogenoma + ensamblaje de novo | [enlace](https://claude.ai/code/artifact/beb5bc85-ac67-4f3f-912f-ab4dc82a0d5d) | [`analysis/reports/genome_resources_report.html`](../analysis/reports/genome_resources_report.html) |
| Primer Design Report | Diseño de primers (bdnf/v1) | [enlace](https://claude.ai/code/artifact/eb740fdb-261b-41a5-b44b-e8530a82c215) | [`analysis/reports/primer_design_report.html`](../analysis/reports/primer_design_report.html) |
