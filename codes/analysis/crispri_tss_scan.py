#!/usr/bin/env python3
"""
CRISPRi guide scan: enumerate SpCas9 NGG candidates in a TSS-proximal
window (-50/+300 bp, per Gilbert et al. 2014 / Horlbeck et al. 2016 - the
window where dCas9-KRAB actually silences transcription; outside it,
efficacy drops sharply regardless of guide quality), classify each against
the Colombian pseudogenome exactly like ko_guide_scan.py does for the CDS
(IDENTICAL/PAM_BROKEN/SEED_VARIANT/DISTAL_VARIANT/NO_ALIGNMENT), and score
with CRISPOR where available.

Companion to ko_guide_scan.py (CRISPRko, CDS-only) - reuses its functions
directly rather than reimplementing gene lookup/alignment/classification.
Requested by the user 2026-09-06 after the "Guppy CRISPR Atlas" report only
covered CRISPRko; see that report's criteria section for why CRISPRi needs
a completely different genomic window (TSS-relative, not CDS-relative) and
cannot be derived from the CDS scan already done.

TSS definition: the mRNA feature's own start (not the CDS start, which is
downstream of any 5'UTR) - for a "-" strand gene this is the mRNA's END
coordinate (GFF always stores start<end regardless of strand). Same
representative transcript (chosen_tid) as ko_guide_scan.py's CDS scan, for
consistency between the two reports.

Window orientation: extracted in genomic (+ strand) coordinates, then
reverse-complemented for "-" strand genes so the returned sequence reads
5'->3' in the gene's own transcriptional direction, with position 0
corresponding to TSS-50 and the last position to TSS+300 - i.e.
tss_relative_position = -50 + local_0based_index for BOTH strands after
this correction.

Usage:
    python3 crispri_tss_scan.py --gene bdnf --population pseudogenome
"""
import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from ko_guide_scan import (  # noqa: E402
    GENOME_CHOICES, OUT_DIR, CRISPOR_SIF, CRISPOR_REF_GENOME_IDS, CRISPOR_GENOME_IDS,
    find_gene_features, faidx_seq, align_cs, parse_cs_variants, find_ngg_candidates,
    classify_candidate, revcomp, run_crispor,
)

WINDOW_UPSTREAM = 50    # bp upstream of TSS (Gilbert 2014 / Horlbeck 2016 CRISPRi window)
WINDOW_DOWNSTREAM = 300  # bp downstream of TSS


def tss_window(mrna_line, upstream, downstream):
    """Return (chrom, win_start, win_end, tss_genomic, strand) - genomic
    1-based inclusive coordinates for samtools faidx."""
    chrom, start, end, strand = mrna_line[0], int(mrna_line[3]), int(mrna_line[4]), mrna_line[6]
    if strand == "+":
        tss = start
        win_start, win_end = tss - upstream, tss + downstream
    else:
        tss = end
        win_start, win_end = tss - downstream, tss + upstream
    return chrom, max(1, win_start), win_end, tss, strand


def oriented_window_seq(fasta_path, chrom, win_start, win_end, strand):
    seq = faidx_seq(fasta_path, chrom, win_start, win_end)
    return revcomp(seq) if strand == "-" else seq


# ── TSS representative choice + isoform coverage (found 2026-09-14) ────────
#
# ko_guide_scan.py's CDS scan picks its representative isoform by longest
# CDS - reused here for consistency, but that criterion is meaningless for
# a TSS/promoter question. Investigated isoform TSS spread across the 8
# genes: most cluster within a few bp (irrelevant - one -50/+300bp window
# already covers them all), but bdnf (3,815bp spread, 6 distinct TSS) and
# agap3 (63,770bp spread) have GENUINE, independently RNA-seq-supported
# alternative transcription start sites - not weak/fragmentary Gnomon
# predictions: every bdnf isoform has "100% coverage...by RNAseq
# alignments" with 7-41 samples; agap3 similarly, 5-44 samples. Consistent
# with BDNF's well-documented multi-promoter architecture in other
# vertebrates (Aid et al. 2007, mouse/rat Bdnf gene structure).
#
# Sharpest finding: agap3's longest-CDS transcript (the one ko_guide_scan.py
# picks) uses the TIGHT TSS cluster (~21068190), but a DIFFERENT transcript
# 63kb further upstream (XM_008434265.2) has the STRONGEST RNA-seq support
# of all 7 isoforms (44 samples, vs 22 for the longest-CDS one) - the
# CDS-based choice was not even the best-evidenced promoter for this gene.
# So the CRISPRi representative is chosen independently here, by strongest
# RNA-seq support (not CDS length), and every other isoform's own TSS is
# checked against the resulting window so an uncovered alternative
# promoter is surfaced rather than silently ignored.


def rnaseq_sample_support(mrna_line):
    """Number of RNA-seq samples supporting all annotated introns, parsed
    from the mRNA's own NCBI Gnomon `model_evidence=` attribute - a direct
    measure of how well-supported this specific isoform's structure
    (including its TSS) is, independent of CDS length. Returns 0 if absent
    (e.g. a curated NM_/NR_ record with no model_evidence at all) rather
    than erroring - a safe, uninformative fallback, not a real 0-support
    claim."""
    m = re.search(r"including (\d+) samples? with support for all annotated introns", mrna_line[8])
    return int(m.group(1)) if m else 0


def choose_tss_representative(ref_mrnas, shared_tids):
    """Pick the shared transcript with the strongest independent RNA-seq
    support for ITS OWN TSS - not the longest CDS (ko_guide_scan.py's
    criterion, irrelevant to promoter usage). Ties fall back to longest
    mRNA span, then to Python's stable max() (first-seen) - deterministic,
    though ties this specific are not expected in practice."""
    return max(
        shared_tids,
        key=lambda t: (rnaseq_sample_support(ref_mrnas[t][0]), int(ref_mrnas[t][0][4]) - int(ref_mrnas[t][0][3])),
    )


def tss_isoform_coverage(ref_mrnas, ref_win_start, ref_win_end, chosen_tss):
    """For every annotated isoform of the gene (not just those shared with
    the population genome - the question is about the reference gene's
    real isoform diversity), does its own TSS fall inside the chosen
    window? Returns (covered_tids, uncovered_info) where uncovered_info is
    a list of (tid, tss, rnaseq_support, distance_from_chosen_tss) for any
    isoform whose TSS the window does NOT reach - i.e. a likely
    alternative promoter this design does not address at all."""
    covered, uncovered = [], []
    for tid, (mrna_line, _) in ref_mrnas.items():
        _, _, _, tss, _ = tss_window(mrna_line, WINDOW_UPSTREAM, WINDOW_DOWNSTREAM)
        if ref_win_start <= tss <= ref_win_end:
            covered.append(tid)
        else:
            uncovered.append((tid, tss, rnaseq_sample_support(mrna_line), abs(tss - chosen_tss)))
    uncovered.sort(key=lambda u: -u[2])  # strongest-supported alternative promoter first
    return covered, uncovered


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--population", choices=list(GENOME_CHOICES.keys()), default="pseudogenome")
    ap.add_argument("--no-crispor", action="store_true")
    args = ap.parse_args()

    gene = args.gene
    pop = GENOME_CHOICES[args.population]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tmp_prefix = f"/tmp/crispri_tss_scan_{gene}"

    print(f"=== {gene} CRISPRi (TSS -{WINDOW_UPSTREAM}/+{WINDOW_DOWNSTREAM}) vs {pop['label']} ===")

    ref_gene, ref_mrnas = find_gene_features(pop["ref_gff"], gene)
    pop_gene, pop_mrnas = find_gene_features(pop["gff"], gene)
    if ref_gene is None or pop_gene is None:
        sys.exit(f"ERROR: gene '{gene}' not found in one of the GFFs")

    shared_tids = set(ref_mrnas) & set(pop_mrnas)
    if not shared_tids:
        sys.exit(f"ERROR: no shared transcript ID for gene '{gene}'")
    chosen_tid = choose_tss_representative(ref_mrnas, shared_tids)
    ref_mrna, _ = ref_mrnas[chosen_tid]
    pop_mrna, _ = pop_mrnas[chosen_tid]
    chosen_support = rnaseq_sample_support(ref_mrna)
    print(f"Representative transcript: {chosen_tid} "
          f"(chosen by strongest RNA-seq support for its TSS: {chosen_support} samples)")

    ref_chrom, ref_win_start, ref_win_end, ref_tss, ref_strand = tss_window(
        ref_mrna, WINDOW_UPSTREAM, WINDOW_DOWNSTREAM
    )
    pop_chrom, pop_win_start, pop_win_end, pop_tss, pop_strand = tss_window(
        pop_mrna, WINDOW_UPSTREAM, WINDOW_DOWNSTREAM
    )
    print(f"Reference TSS: {ref_chrom}:{ref_tss} ({ref_strand}) - window {ref_win_start}-{ref_win_end}")
    print(f"{args.population} TSS: {pop_chrom}:{pop_tss} ({pop_strand}) - window {pop_win_start}-{pop_win_end}")

    # --- isoform TSS coverage: does this window also reach every OTHER
    # annotated isoform's own transcription start, or is there a genuinely
    # separate, uncovered alternative promoter? See note above
    # choose_tss_representative() for why this matters (bdnf, agap3). ---
    covered_tids, uncovered = tss_isoform_coverage(ref_mrnas, ref_win_start, ref_win_end, ref_tss)
    n_total_isoforms = len(ref_mrnas)
    print(f"TSS window covers {len(covered_tids)}/{n_total_isoforms} annotated isoform(s)' own TSS")
    if uncovered:
        print("NOT covered - likely alternative promoter(s), not addressed by any guide in this window:")
        for tid, tss_u, supp_u, dist in uncovered:
            print(f"  {tid}: TSS={tss_u} ({supp_u} RNA-seq samples), {dist}bp from the chosen TSS")

    ref_seq = oriented_window_seq(pop["ref_fasta"], ref_chrom, ref_win_start, ref_win_end, ref_strand)
    pop_seq = oriented_window_seq(pop["fasta"], pop_chrom, pop_win_start, pop_win_end, pop_strand)
    print(f"Reference window length: {len(ref_seq)}bp")
    print(f"{args.population} window length: {len(pop_seq)}bp")

    pos_c, cigar_c, cs_tag_c, flag_c = align_cs(ref_seq, pop_seq, tmp_prefix)
    if not cs_tag_c:
        sys.exit(f"ERROR: TSS windows for {gene} did not align (unmapped)")
    variants, pos_map = parse_cs_variants(cs_tag_c, pos_c)
    variants_by_pos = {v[0]: v for v in variants}
    print(f"Variants in TSS window: {len(variants)}")

    crispor_ref_scores, crispor_pop_scores = {}, {}
    if not args.no_crispor and CRISPOR_SIF.exists():
        print("\n=== Running CRISPOR (Singularity) ===")
        ref_genome_id = CRISPOR_REF_GENOME_IDS.get(pop["ref_fasta"])
        if ref_genome_id:
            crispor_ref_scores = run_crispor(
                ref_seq, ref_genome_id, tmp_prefix + "_crispor_ref", str(OUT_DIR / f"{gene}_reference_crispri")
            )
            print(f"CRISPOR reference guides scored: {len(crispor_ref_scores)}")
        pop_genome_id = CRISPOR_GENOME_IDS.get(args.population)
        if pop_genome_id:
            crispor_pop_scores = run_crispor(
                pop_seq, pop_genome_id, tmp_prefix + "_crispor_pop", str(OUT_DIR / f"{gene}_{args.population}_crispri")
            )
            print(f"CRISPOR {args.population} guides scored: {len(crispor_pop_scores)}")

    ref_candidates = find_ngg_candidates(ref_seq)
    rows = []
    for start0, strand, spacer, pam in ref_candidates:
        cls, qwin, note = classify_candidate(start0, spacer, pos_map, variants_by_pos, pop_seq)
        target = spacer + pam
        ref_c = crispor_ref_scores.get(target, {})
        pop_c = crispor_pop_scores.get(qwin, {}) if qwin else {}
        tss_rel_start = -WINDOW_UPSTREAM + start0
        tss_rel_end = tss_rel_start + 22
        rows.append({
            "gene": gene,
            "tss_rel_start": tss_rel_start,
            "tss_rel_end": tss_rel_end,
            "strand_in_window": strand,
            "ref_spacer": spacer,
            "ref_pam": pam,
            "classification": cls,
            "pop_window": qwin if qwin else "",
            "note": note,
            "ref_crispor_mitSpecScore": ref_c.get("mitSpecScore", ""),
            "ref_crispor_offtargetCount": ref_c.get("offtargetCount", ""),
            "ref_crispor_doench16_score": ref_c.get("Doench '16-Score", ""),
            "pop_crispor_mitSpecScore": pop_c.get("mitSpecScore", ""),
            "pop_crispor_offtargetCount": pop_c.get("offtargetCount", ""),
            # Constant per gene (the window is fixed once chosen, not
            # per-candidate) - kept per-row for easy filtering/joining,
            # same convention as ko_guide_scan.py's isoform_coverage_*.
            "tss_window_isoform_coverage_n": len(covered_tids),
            "tss_window_isoform_coverage_total": n_total_isoforms,
        })

    out_csv = OUT_DIR / f"{gene}_{args.population}_crispri_candidates.csv"
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else [])
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    counts = Counter(r["classification"] for r in rows)
    print(f"\n=== CRISPRi candidate summary ({WINDOW_UPSTREAM}bp upstream / {WINDOW_DOWNSTREAM}bp downstream of TSS) ===")
    print(f"Total NGG candidates in window: {len(rows)}")
    for cls in ["IDENTICAL", "PAM_BROKEN", "SEED_VARIANT", "DISTAL_VARIANT", "NO_ALIGNMENT"]:
        print(f"  {cls}: {counts.get(cls, 0)}")
    print(f"\nWritten: {out_csv}")


if __name__ == "__main__":
    main()
