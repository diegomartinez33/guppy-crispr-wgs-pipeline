# Acceso al clúster y archivos pesados

## Cuenta compartida

- **Clúster:** hypatia (Universidad de los Andes) — scheduler SLURM.
- **Usuario compartido:** `guppy-genome@hypatia.uniandes.edu.co`
- **Contraseña:** solicítala directamente a Diego Martínez (da.martinez33,
  diegoandres3322@gmail.com) — **no está en este repositorio ni se distribuye por GitHub**.

Esta cuenta existe específicamente para que los colegas puedan correr los dos pipelines
parametrizados para reutilización directa — **diseño de guías CRISPR** y **diseño de primers**
(ver [TUTORIAL.md](TUTORIAL.md)) — sin necesitar acceso a la cuenta personal del autor.

## Qué se copia a la cuenta compartida (y qué no)

El repositorio de código (`git clone`) se copia completo — es pequeño. La mayoría de los datos
generados por el pipeline **no** se versionan en git (`.gitignore` los excluye — ver comentario
en ese archivo, ~2.6TB en total hoy). De esos, solo una parte pequeña y bien definida se copia a
la cuenta compartida: exactamente lo que necesitan los dos pipelines pensados para uso de
colegas. El resto queda documentado aquí de forma informativa — si algún colega necesita
replicar una etapa completa del pipeline (QC/trimming, mapeo, GATK, ensamblaje de novo), lo hará
con sus propias muestras futuras siguiendo [PIPELINE.md](PIPELINE.md)/[TUTORIAL.md](TUTORIAL.md)
como guía, no reutilizando estos datos intermedios.

| Directorio | Tamaño aprox. | Qué contiene / lo genera | ¿Se copia? |
|---|---|---|---|
| `reference/` | 15G | Genomas v1+v2 + índices, pseudogenoma colombiano, genoma ensamblado de novo, BLAST dbs — la línea base de todo análisis posible; pseudogenoma y ensamblado ya llevan la información de variantes de la población colombiana | ✅ **Sí** |
| `codes/analysis/crispor_singularity/` | 6.5G | Contenedor Singularity de CRISPOR + genomas registrados — requerido por `ko_guide_scan.py`/`crispri_tss_scan.py` ([PIPELINE.md §8](PIPELINE.md#8-diseño-de-guías-crispr-ko--crispri-por-gen)) | ✅ **Sí** |
| `raw_fastq/` | 253G | FASTQ crudos de las 15 muestras ([PIPELINE.md §1](PIPELINE.md#1-qc-y-recorte-de-lecturas)) | ❌ No — no lo necesita ningún pipeline de colegas; para replicar QC/trimming, usar muestras propias |
| `trimmed_fastp/` | 251G | Lecturas recortadas con fastp (solo comparación, no usado río abajo) | ❌ No |
| `trimmed_trimmomatic/` | 241G | Lecturas recortadas con Trimmomatic (el recortador usado río abajo) | ❌ No — para remapear, usar muestras propias |
| `mapping/` | 784G | BAMs ordenados BWA-MEM, v1+v2, individuales+fusionados por grupo ([PIPELINE.md §2](PIPELINE.md#2-mapeo-al-genoma-de-referencia)) | ❌ No |
| `gatk/` | 610G | BAMs deduplicados, GVCFs, VCFs, workspace GenomicsDB (v1+v2) ([PIPELINE.md §3](PIPELINE.md#3-llamado-de-variantes--gatk)) | ❌ No — ninguno de los dos pipelines de colegas necesita BAMs/VCFs de muestra |
| `crispresso/` | 192M | Reportes CRISPResso2 (pequeño — mayormente texto/HTML) | ❌ No (ya está resumido en los reportes visuales, ver [RESULTS.md](RESULTS.md)) |
| `crispresso_v2/` | 58K | Migración v2, apenas iniciada | ❌ No |
| `assembly/` | 506G | Directorios de trabajo intermedios de SPAdes/RagTag/gap-filling/QC — las salidas FINALES de este proceso ya viven en `reference/colombian_scaffolded_genome/`, que sí se copia | ❌ No |
| `intermediate_files/` | 1.1G | Varios, no esencial | ❌ No |

**Total a copiar: ≈21.5G** (no los ~2.6TB completos). Ni el diseño de guías ni el de primers
necesitan nada más pesado — el CSV de sitios off-target (`combined_offtargets.csv`) ya es
pequeño y vive en git; `primer3_env`/el módulo EMBOSS son un entorno conda + módulo del clúster
que cada colega construye localmente con `codes/analysis/setup_primer3.sh`, no un archivo para
copiar.

## Cómo copiar (ejemplo, lo ejecuta cada usuario con sus propias credenciales)

```bash
# Desde tu propio usuario en hypatia, hacia la cuenta compartida:
rsync -avP --info=progress2 \
  reference/ \
  guppy-genome@hypatia.uniandes.edu.co:/ruta/destino/off-target_data/reference/

rsync -avP --info=progress2 \
  codes/analysis/crispor_singularity/ \
  guppy-genome@hypatia.uniandes.edu.co:/ruta/destino/off-target_data/codes/analysis/crispor_singularity/
```

Ajusta `/ruta/destino/` al home real de la cuenta compartida. Este agente no ejecuta esta copia
— requiere credenciales de la cuenta compartida que no se comparten en esta conversación, y es
una operación entre cuentas de usuario del clúster que corresponde iniciar y supervisar
directamente al dueño de los datos.
