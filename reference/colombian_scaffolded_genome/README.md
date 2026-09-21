# Colombian Guppy (*Poecilia reticulata*) Draft Scaffolded Genome

**Status:** Draft assembly, not reference-quality. Read the "Limitations" section
before using this genome for any downstream analysis.

**Author:** da.martinez33 | Universidad de los Andes (Colombia)
**Date assembled:** July 2026
**Contact:** diegoandres3322@gmail.com

**Note (2026-09-20):** everything below describes the assembly scaffolded against the
v1 reference (GCF_000633615.1). A second assembly, scaffolded the same way against the
newer v2 reference (GCF_904066995.2), also exists at `reference/colombian_scaffolded_genome_v2/`
and is now fully complete, including QUAST QC. v2 wins on most quality metrics (BUSCO C:96.3%
vs v1's 95.5%; bdnf Liftoff coverage=0.980/sequence_ID=0.964 vs v1's 0.960/0.957; QUAST
misassemblies 13,768 vs v1's 23,914; NA50 159Kb vs v1's 114Kb) — the one exception is genome
fraction (v1: 92.2%, v2: 85.5%), which is a reference-composition effect rather than a real
assembly regression, explained in `analysis/reports/genome_resources_report.html`'s "v1 vs v2"
section. See [docs/RESULTS.md](../../docs/RESULTS.md#7-migration-to-the-v2-reference-genome--complete)
for the full comparison and CLAUDE.md §8 for build detail. This document is not duplicated per
version.

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

## Polishing experiment (2026-08) — did not improve the genome

Tried short-read polishing (NextPolish v1.4.1, 2 rounds, same Illumina
reads used to build the assembly) to see if it would fix the internal
stop codons found in some BUSCO genes. **Result: no meaningful
improvement, and slightly worse on several structural metrics.** The
polished genome was NOT adopted — `colombian_scaffolded.fna` (this
directory) remains the authoritative version.

| Metric | Original | Polished | Change |
|---|---|---|---|
| N50 | 28.31 Mb | 28.30 Mb | ≈ same |
| Genome fraction | 82.807% | 82.038% | **worse (-0.77pp)** |
| Duplication ratio | 1.078 | 1.087 | **worse** |
| # misassemblies | 7,589 | 8,043 | **worse (+6%)** |
| NA50 (aligned N50) | 198,159 | 193,546 | **worse** |
| Mismatches/100kbp | 567.95 | 557.72 | better (~1.8%) |
| Indels/100kbp | 137.54 | 136.60 | better (~0.7%) |
| BUSCO Complete | 87.1% (3169) | 87.1% (3173) | ≈ same |
| BUSCO genes w/ internal stop codons | 220 | 221 | ≈ same (no fix) |

Why it didn't help: NextPolish corrects disagreements between the
assembly and aligned reads — but the reads used for polishing were the
*same reads already used to build the assembly* with SPAdes. Where the
assembly already matched what those reads say, there was nothing to
correct, even where that sequence differs from the reference gene models
BUSCO checks against. The small increase in misassemblies suggests
polishing even introduced a few new small-scale artifacts. Takeaway: the
220 internal-stop-codon genes are more likely genuine population
divergence (or gene-prediction quirks) than assembly errors — and closing
that gap will need genuinely new information (long reads), not
reprocessing the same short reads. See the Known Issues section of
`CLAUDE.md` for the full technical debugging history (two real bugs found
along the way in NextPolish itself, both fixed) and the "Further
Improvement Options" roadmap for what's next.

## Gap-filling experiment (2026-08) — mixed result: big completeness gain, some structural cost

Unlike short-read polishing (above), which reprocessed the same Illumina
reads and found nothing new to correct, gap-filling adds genuinely new
information: ~2.4Gbp of existing Nanopore long reads (pooled telencephalon,
10 fish, same Colombian population, ~3.2x depth — too shallow for a full
hybrid reassembly, but enough to span some of the gaps RagTag left behind).
Tool: TGS-GapCloser v1.2.1, error-corrected with racon. Result: **56.3% of
gap regions filled** (179,359 of 318,572), a real, substantial change — but
QUAST shows a genuine trade-off, not a clean win.

| Metric | Original | Gap-filled | Change |
|---|---|---|---|
| N50 | 28.31 Mb | 29.43 Mb | better (+4%) |
| **Genome fraction** | 82.807% | **92.088%** | **much better (+9.3pp)** |
| Duplication ratio | 1.078 | 1.042 | better |
| N's per 100kbp | 7,145 | 1,564 | much better (-78%) |
| **# misassemblies** | 7,589 | **23,532** | **much worse (+210%)** |
| **NA50 (aligned N50)** | 198,159 | **115,382** | **worse (-42%)** |
| Mismatches/100kbp | 567.95 | 700.33 | worse |
| Indels/100kbp | 137.54 | 210.08 | worse |
| BUSCO Complete | 87.1% (3169) | **95.4% (3476)** | **much better** |
| BUSCO Missing | 6.7% (245) | **2.1% (73)** | **much better** |
| BUSCO Fragmented | 6.2% (226) | **2.5% (91)** | **much better** |

Why the trade-off: filling a gap replaces an `N` placeholder with real
sequence derived from Nanopore reads, corrected only with `racon` (long-read
consensus) — never checked against the high-precision Illumina reads. That
new sequence is where essentially all of the completeness gain (genome
fraction, BUSCO) comes from, but it's also where the new mismatches/indels
and most of the new misassemblies almost certainly live (a gap that was
previously invisible `N` can't be "misassembled" — once it's real sequence,
QUAST can now detect local disagreement there against the fragmented
reference). Net effect: substantially more complete, but with a
precision cost concentrated in the newly-filled regions specifically —
see "Targeted post-gap-fill polishing" below for the follow-up meant to
recover some of that precision without giving up the completeness gain.

**Update (2026-09-17): adopted.** Once the targeted post-gap-fill polishing below
completed and QC'd clean, the gap-filled+polished genome was promoted to the canonical
`colombian_scaffolded.fna` path (the pre-gap-fill scaffold was archived, not deleted, to
`pre_gapfill_archive/`). All downstream artifacts (Liftoff annotation, BWA/BLAST/samtools
indices) were rebuilt against it — see the Changelog entry below and CLAUDE.md §8 for the
"v1 de novo assembly promotion bug" writeup (this promotion was initially missed, then
found + fixed 2026-09-17).

## Targeted post-gap-fill polishing (2026-08/09) — complete, adopted

Rationale: whole-genome short-read polishing already failed once (see
"Polishing experiment" above) because it reprocessed reads the assembly
was already built from — there was nothing left to correct. The
newly-filled gap sequence is different: it was never built from or
compared against the Illumina reads at all (it came from Nanopore +
racon), so genuine disagreements between it and the Illumina reads should
exist, and NextPolish should have real work to do there this time.

Method: re-ran NextPolish (same tool/config as the original polishing
experiment, 2 rounds, same 3 Control samples' Illumina reads) on top of
the GAP-FILLED genome rather than the original. This is a whole-genome
NextPolish run, not literally restricted to the filled coordinates, but
the expected/intended effect is that it acts as a *targeted* correction in
practice: regions already consistent with the Illumina reads (i.e.
everything that was already in the pre-gap-fill genome) should see
little-to-no change, exactly as observed in the first polishing
experiment, while the newly-filled Nanopore-derived regions - which have
never been reconciled with Illumina data - are where real corrections
should land.

Result (2026-09-07 QC, job 716452): a clean, low-risk win on top of gap-filling — BUSCO
internal stop codons dropped 184→136 (-26%, the specific Nanopore-indel artifact this
polish targeted), Complete ticked up 95.4%→95.5%, at essentially zero structural cost
(genome fraction 92.088%→92.209%, misassemblies 23,532→23,914 — within noise, duplication
ratio unchanged at 1.042). Does not recover the structural precision lost during
gap-filling itself (that trade-off is inherent to TGS-GapCloser's long-read gap-filling and
was already the accepted cost — see "Gap-filling experiment" above). **Adopted as the
authoritative genome 2026-09-17** (see note near the top of this document).

## Files in this directory

| File | Description |
|---|---|
| `colombian_scaffolded.fna` | The final, adopted genome sequence (686MB) — gap-filled (TGS-GapCloser) + polished (NextPolish), promoted to this canonical path 2026-09-17 |
| `colombian_scaffolded.fna.{amb,ann,bwt,pac,sa}` | BWA index (pre-built, ready for read mapping) |
| `colombian_scaffolded.fna.fai` | samtools faidx index |
| `colombian_scaffolded.dict` | GATK/Picard sequence dictionary |
| `colombian_scaffolded.liftoff.gff3` | Gene annotations transferred from the Trinidad reference |
| `blast_db/colombian_scaffolded.n*` | Pre-built nucleotide BLAST database |
| `pre_gapfill_archive/` | The original pre-gap-fill RagTag scaffold + its own indices/annotation (superseded 2026-09-17, kept for reference rather than deleted) |

## Recommended citation / acknowledgment

If you use this genome, please note internally that it is a draft assembly
built via SPAdes v4.0.0 + RagTag v2.1.0 + Liftoff v1.5.1, scaffolded against
NCBI RefSeq assembly GCF_000633615.1 (Guppy_female_1.0_MT), and check with
da.martinez33 (diegoandres3322@gmail.com) before using it in any publication
or as a dependency for other analyses.

## Changelog

- **2026-07-15** — Initial draft assembly completed (SPAdes → RagTag →
  Liftoff). QUAST/BUSCO QC run. This README written.
- **2026-08-11** — Tested NextPolish short-read polishing (2 rounds).
  Result: no meaningful improvement (see "Polishing experiment" section
  above); genome not replaced. `colombian_scaffolded.fna` remains
  authoritative.
- **2026-08-31** — Tested TGS-GapCloser gap-filling with existing Nanopore
  data (56.3% of gaps filled). Result: mixed — genome fraction 82.8%→92.1%
  and BUSCO 87.1%→95.4% (both much better), but misassemblies 7,589→23,532
  and NA50 198K→115K (both worse), concentrated in the newly-filled
  regions (see "Gap-filling experiment" section above). Not yet adopted.
  Started a targeted post-gap-fill NextPolish run to recover precision in
  the newly-filled regions specifically (see "Targeted post-gap-fill
  polishing" section) — result pending.
- **2026-09-07** — Targeted post-gap-fill NextPolish QC completed (job
  716452): BUSCO internal stop codons 184→136 (-26%), Complete
  95.4%→95.5%, structural metrics essentially flat (genome fraction
  92.088%→92.209%, misassemblies 23,532→23,914, duplication ratio
  unchanged) — clean win, no new cost.
- **2026-09-17** — Promoted the gap-filled+polished genome to the
  canonical `colombian_scaffolded.fna` path (previously it sat only in
  `assembly/nextpolish_output_gapfilled/`, a bug not caught until this
  date — see CLAUDE.md §8 "v1 de novo assembly promotion bug"). The
  pre-gap-fill scaffold was archived to `pre_gapfill_archive/`, not
  deleted. Liftoff annotation and all indices rebuilt against the real
  final genome; bdnf Liftoff quality improved from
  coverage=0.945/sequence_ID=0.923 (stale scaffold) to
  coverage=0.960/sequence_ID=0.957 (actual final assembly). This genome
  is now genuinely authoritative, not just intended to be.
