# Colombian Guppy (*Poecilia reticulata*) Draft Scaffolded Genome

**Status:** Draft assembly, not reference-quality. Read the "Limitations" section
before using this genome for any downstream analysis.

**Author:** da.martinez33 | Universidad de los Andes (Colombia)
**Date assembled:** July 2026
**Contact:** diegoandres3322@gmail.com

---

## What this is

A population-specific *de novo* assembly of the Colombian *Poecilia reticulata*
(guppy) population, built to support CRISPR-Cas9 off-target analysis. It
captures structural variation (indels, novel sequence) between the Colombian
population and the standard reference genome that simple variant-calling
against the reference cannot recover.

This is one of **two** Colombian-population genome resources in this project —
they are complementary, not interchangeable:

| | This genome | `reference/pseudogenome/` |
|---|---|---|
| Method | *De novo* co-assembly + reference-guided scaffolding | Reference + majority-allele SNP/indel substitution (bcftools consensus) |
| Captures structural variants? | Yes (in principle) | No — single linear copy of the reference with point/short-indel edits only |
| Built from | 3 pooled Control individuals, de novo | 3 pooled Control individuals, variant calls only |

If you need a fast, low-risk drop-in replacement for the reference genome
(e.g., for read mapping or variant calling), the pseudogenome is usually the
safer choice. Use **this** genome if you specifically need structural
information not present in the reference, or you're validating something the
pseudogenome can't represent.

## How it was built

Pipeline: **SPAdes → RagTag → Liftoff**, scaffolded against the Trinidad/Guanapo
reference (GCF_000633615.1).

| Step | Tool | What it did |
|---|---|---|
| 1. De novo co-assembly | SPAdes v4.0.0 (`--isolate` mode) | Assembled short reads (150bp Illumina, 3 pooled Control individuals, ~47GB compressed) into contigs |
| 2. Reference-guided scaffolding | RagTag v2.1.0 | Ordered and oriented contigs into chromosome-scale scaffolds using the Trinidad reference as a guide |
| 3. Annotation transfer | Liftoff v1.5.1 | Transferred Trinidad reference gene models onto the new scaffold coordinates |

Input data: 3 Control replicate WGS samples (Control_MNP_I/II/III), Illumina
NovaSeq X, 2×150bp paired-end, trimmomatic-trimmed reads.

## Quality summary

Full reports: `assembly/qc_results/quast/` and `assembly/qc_results/busco/`
in the parent project. Plain-language summary below.

**Structural contiguity (QUAST):**
- N50 = **28.3 Mb** — genuinely chromosome-scale contiguity, achieved via
  reference-guided scaffolding (the raw short-read assembly alone had N50
  of only 669bp before scaffolding).
- Genome fraction = **82.8%** of the reference genome has a matching region
  in this assembly. The gap is expected — real divergence between the
  Colombian and Trinidad populations, a reference that is itself fragmented
  into 2,768 pieces, and short (<500bp) contigs filtered out before
  scaffolding.
- Duplication ratio = **1.078** — low/healthy; this assembly is only ~8%
  "larger" than the reference, meaning it mostly avoided assembling
  redundant duplicate copies of the same region despite being built from 3
  pooled individuals.

**Gene-content completeness (BUSCO, reference-independent):**
```
C:87.1% [S:86.2%, D:0.9%], F:6.2%, M:6.7%   (n=3,640 core ray-finned-fish genes)
```
- **87.1% complete** — a solid "good draft quality" score (reference-grade
  genomes typically score 95-99%; below ~80% is considered fragmented).
- Only 0.9% duplicated — consistent with the low QUAST duplication ratio;
  the assembly mostly represents genes once, as expected.
- 6.2% fragmented, 6.7% missing — some combination of real sequence
  divergence, coverage gaps, and short-read assembly limits in complex or
  repetitive regions.

## Limitations — read before using or sharing further

1. **This is a draft, not a reference-quality genome.** Chromosome-scale
   contiguity, but base-level accuracy is imperfect in places — of the
   3,169 "complete" BUSCO genes, 220 (~6.9%) contain internal stop codons,
   almost certainly small assembly errors (indels/frameshifts) rather than
   real biology. Don't draw conclusions from a single gene's sequence
   without checking it first.

2. **This represents a pooled composite of 3 individuals, not one fish's
   exact genome.** Some regions may reflect collapsed/averaged heterozygous
   variation across the 3 Control individuals rather than a single
   consistent haplotype. Not suitable for individual-level haplotype
   analysis.

3. **Coordinates are NOT interchangeable with the Trinidad reference
   (GCF_000633615.1).** Sequence names carry a `_RagTag` suffix and
   positions differ from the reference due to inserted gaps and
   rearrangements. Use a proper liftover/alignment, never assume 1:1
   coordinate correspondence.

4. **`Chr0_RagTag` (39Mb) is not a real chromosome.** It's an artificial
   concatenation of ~2,000 small contigs that RagTag could not confidently
   place, with no true biological order or continuity between them. Any
   analysis assuming real gene order, synteny, or structural continuity
   within Chr0 will produce nonsense — treat it as an "unplaced sequence
   bin," not a chromosome.

5. **~7% of the sequence is N (gap filler)**, inserted at every point RagTag
   joined two originally-separate fragments. This is not missing/deleted
   real sequence — it's a placeholder for unknown gap length.

6. **Annotation confidence varies by gene.** The GFF3 was transferred via
   Liftoff from the Trinidad reference. High confidence for genes
   individually checked (e.g., *bdnf*: 94.5% coverage, 92.3% sequence
   identity) — but genome-wide annotation quality was not verified
   gene-by-gene. Treat annotations in poorly-aligned or missing regions as
   lower-confidence.

7. **Built for a specific purpose**: this genome was assembled to support
   CRISPR-Cas9 off-target analysis for the Colombian guppy population. It
   may not be validated for other use cases (e.g., population genomics
   requiring precise per-base accuracy) without further work.

## Files in this directory

| File | Description |
|---|---|
| `colombian_scaffolded.fna` | The genome sequence (692MB, 692Mb total length) |
| `colombian_scaffolded.fna.{amb,ann,bwt,pac,sa}` | BWA index (pre-built, ready for read mapping) |
| `colombian_scaffolded.fna.fai` | samtools faidx index |
| `colombian_scaffolded.dict` | GATK/Picard sequence dictionary |
| `colombian_scaffolded.liftoff.gff3` | Gene annotations transferred from the Trinidad reference |
| `blast_db/colombian_scaffolded.n*` | Pre-built nucleotide BLAST database |

## Recommended citation / acknowledgment

If you use this genome, please note internally that it is a draft assembly
built via SPAdes v4.0.0 + RagTag v2.1.0 + Liftoff v1.5.1, scaffolded against
NCBI RefSeq assembly GCF_000633615.1 (Guppy_female_1.0_MT), and check with
da.martinez33 (diegoandres3322@gmail.com) before using it in any publication
or as a dependency for other analyses.

## Changelog

- **2026-07-15** — Initial draft assembly completed (SPAdes → RagTag →
  Liftoff). QUAST/BUSCO QC run. This README written.
