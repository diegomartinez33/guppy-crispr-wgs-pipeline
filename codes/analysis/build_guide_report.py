#!/usr/bin/env python3
"""
Build a consolidated report of CRISPOR-predicted guides (not the manual
NGG scan) + population variant impact, for all genes already analyzed by
ko_guide_scan.py / crispri_tss_scan.py.

Design (changed 2026-09-08, per user request): the candidate guide list
and its metrics/off-targets now come directly from CRISPOR's own output
files (<gene>_<genome>_crispor_guides.tsv / _crispor_offs.tsv, written by
run_crispor() in ko_guide_scan.py and crispri_tss_scan.py) instead of this
script's earlier home-grown ranking of the manual NGG scan. The manual
scan's *_guide_comparison.csv / *_crispri_candidates.csv are still used,
but only as the source of population-safety classification
(IDENTICAL/PAM_BROKEN/SEED_VARIANT/DISTAL_VARIANT/EXON_JUNCTION_ARTIFACT)
and gene-structure position (CDS %/last-exon for KO, TSS offset for
CRISPRi) - CRISPOR's own guide enumeration doesn't know about exon
structure or the population genome at all, so this join is still required
for guide safety, not optional.

Fallback: agap3 (both contexts) and grin1a/gria1a (KO/CDS context only)
have NO CRISPOR guide files at all - the reference genome's CDS/TSS window
contains an IUPAC ambiguity code that crashes the CRISPOR container's
revComp()/Azimuth model (see CLAUDE.md). For those, falls back to the
manual-scan-only candidate list, clearly flagged in the output.
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from ko_guide_scan import GENOME_CHOICES, OUT_DIR, find_gene_features, exon_junction_boundaries  # noqa: E402

GENES = ["bdnf", "agap3", "grin1a", "grin1b", "gria1a", "gria1b", "gria2b", "nlgn1"]
POPULATION = "pseudogenome"
TOP_N = 5
TOP_OFFTARGETS = 3


def to_float(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def to_int(x):
    try:
        return int(x)
    except (ValueError, TypeError):
        return None


def last_exon_span(cds_lines, strand):
    lengths = [int(f[4]) - int(f[3]) + 1 for f in cds_lines]
    if strand == "-":
        lengths = list(reversed(lengths))
    total = sum(lengths)
    return total - lengths[-1] + 1, total


def load_tsv(path):
    if not path.exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def offtargets_by_guide_id(offs_path):
    by_guide = defaultdict(list)
    for row in load_tsv(offs_path):
        mm = to_int(row.get("mismatchCount"))
        cfd = to_float(row.get("cfdOfftargetScore"))
        if mm is None:
            continue
        row["_mm"] = mm
        row["_cfd"] = cfd if cfd is not None else -1.0
        by_guide[row["guideId"]].append(row)
    for gid in by_guide:
        by_guide[gid].sort(key=lambda r: (r["_mm"], -r["_cfd"]))
    return by_guide


def offtarget_summary(rows):
    return [
        {
            "chrom": r["chrom"],
            "start": to_int(r["start"]),
            "end": to_int(r["end"]),
            "strand": r["strand"],
            "mismatch_count": r["_mm"],
            "mit_score": to_float(r.get("mitOfftargetScore")),
            "cfd_score": r.get("_cfd") if r.get("_cfd", -1) >= 0 else None,
        }
        for r in rows[:TOP_OFFTARGETS]
    ]


def build_crispor_candidates(gene, ref_fasta_key, is_crispri, manual_rows, manual_key_fn,
                              position_fn, out_prefix):
    """Return (candidates, crispor_available). candidates come from CRISPOR's
    own guide TSV, joined with manual classification via target sequence."""
    guide_path = OUT_DIR / f"{out_prefix}_reference_crispor_guides.tsv" if not is_crispri else \
        OUT_DIR / f"{out_prefix}_reference_crispri_crispor_guides.tsv"
    offs_path = OUT_DIR / f"{out_prefix}_reference_crispor_offs.tsv" if not is_crispri else \
        OUT_DIR / f"{out_prefix}_reference_crispri_crispor_offs.tsv"
    pop_guide_path = OUT_DIR / f"{out_prefix}_pseudogenome_crispor_guides.tsv" if not is_crispri else \
        OUT_DIR / f"{out_prefix}_pseudogenome_crispri_crispor_guides.tsv"

    guides = load_tsv(guide_path)
    if not guides:
        return [], False

    offs_by_guide = offtargets_by_guide_id(offs_path)
    pop_guides = {g["targetSeq"]: g for g in load_tsv(pop_guide_path)}
    manual_by_seq = {manual_key_fn(r): r for r in manual_rows}

    candidates = []
    for g in guides:
        target = g["targetSeq"]
        manual = manual_by_seq.get(target)
        if manual is None:
            continue  # CRISPOR's own PAM scan found a candidate our manual scan didn't (boundary effect) - skip, we need the manual join for classification/position
        pop_g = pop_guides.get(target, {})
        cand = {
            "target_seq": target,
            "spacer": target[:20],
            "pam": target[20:23],
            "guide_id_crispor": g["guideId"],
            "classification": manual["classification"],
            "mit_spec_score": to_float(g.get("mitSpecScore")),
            "cfd_spec_score": to_float(g.get("cfdSpecScore")),
            "offtarget_count": to_int(g.get("offtargetCount")),
            "pop_offtarget_count": to_int(pop_g.get("offtargetCount")),
            "top_offtargets": offtarget_summary(offs_by_guide.get(g["guideId"], [])),
        }
        if not is_crispri:
            cand["doench16_score"] = to_float(g.get("Doench '16-Score"))
            cand["moreno_mateos_score"] = to_float(g.get("Moreno-Mateos-Score"))
        cand.update(position_fn(manual))
        candidates.append(cand)
    return candidates, True


def main():
    pop = GENOME_CHOICES[POPULATION]
    report = {}

    for gene in GENES:
        ref_gene, ref_mrnas = find_gene_features(pop["ref_gff"], gene)
        pop_gene, pop_mrnas = find_gene_features(pop["gff"], gene)
        shared_tids = set(ref_mrnas) & set(pop_mrnas)
        chosen_tid = max(shared_tids, key=lambda t: sum(int(f[4]) - int(f[3]) + 1 for f in ref_mrnas[t][1]))
        ref_strand = ref_gene[6]
        _, ref_cds = ref_mrnas[chosen_tid]
        last_start, last_end = last_exon_span(ref_cds, ref_strand)
        cds_len = last_end

        # ---- KO (CDS) ----
        ko_csv = OUT_DIR / f"{gene}_{POPULATION}_guide_comparison.csv"
        ko_rows = list(csv.DictReader(open(ko_csv)))
        ko_counts = {}
        for r in ko_rows:
            ko_counts[r["classification"]] = ko_counts.get(r["classification"], 0) + 1
        ko_variant_rows = [r for r in ko_rows if r["classification"] in ("PAM_BROKEN", "SEED_VARIANT", "DISTAL_VARIANT")]
        for r in ko_variant_rows:
            pos = int(r["cds_pos_1based"])
            r["_in_last_exon"] = last_start <= pos <= last_end

        def ko_key(r):
            return r["ref_spacer"] + r["ref_pam"]

        def ko_position(manual_row):
            pos = int(manual_row["cds_pos_1based"])
            return {
                "cds_pos_1based": pos,
                "frac_cds": round(pos / cds_len, 3),
                "strand": manual_row["strand"],
                "in_last_exon": last_start <= pos <= last_end,
            }

        ko_candidates, ko_crispor_available = build_crispor_candidates(
            gene, None, False, ko_rows, ko_key, ko_position, gene
        )
        if ko_crispor_available:
            ko_identical = [c for c in ko_candidates if c["classification"] == "IDENTICAL"]
            ko_identical.sort(key=lambda c: (
                c["in_last_exon"],
                -(c["mit_spec_score"] if c["mit_spec_score"] is not None else -1),
                c["offtarget_count"] if c["offtarget_count"] is not None else 9999,
                -(c["doench16_score"] if c["doench16_score"] is not None else -1),
            ))
            ko_top = ko_identical[:TOP_N]
        else:
            fallback = sorted(
                (r for r in ko_rows if r["classification"] == "IDENTICAL"),
                key=lambda r: (last_start <= int(r["cds_pos_1based"]) <= last_end, int(r["cds_pos_1based"]))
            )[:TOP_N]
            ko_top = [{
                "target_seq": r["ref_spacer"] + r["ref_pam"], "spacer": r["ref_spacer"], "pam": r["ref_pam"],
                "guide_id_crispor": None, "classification": "IDENTICAL",
                "mit_spec_score": None, "cfd_spec_score": None, "offtarget_count": None,
                "pop_offtarget_count": None, "top_offtargets": [], "doench16_score": None,
                "moreno_mateos_score": None, **ko_position(r),
            } for r in fallback]

        # ---- CRISPRi (TSS window) ----
        ci_csv = OUT_DIR / f"{gene}_{POPULATION}_crispri_candidates.csv"
        ci_rows = list(csv.DictReader(open(ci_csv)))
        ci_counts = {}
        for r in ci_rows:
            ci_counts[r["classification"]] = ci_counts.get(r["classification"], 0) + 1
        ci_variant_rows = [r for r in ci_rows if r["classification"] in ("PAM_BROKEN", "SEED_VARIANT", "DISTAL_VARIANT")]

        def ci_key(r):
            return r["ref_spacer"] + r["ref_pam"]

        def ci_position(manual_row):
            return {
                "tss_rel_start": int(manual_row["tss_rel_start"]),
                "tss_rel_end": int(manual_row["tss_rel_end"]),
                "strand_in_window": manual_row["strand_in_window"],
            }

        ci_candidates, ci_crispor_available = build_crispor_candidates(
            gene, None, True, ci_rows, ci_key, ci_position, gene
        )
        if ci_crispor_available:
            ci_identical = [c for c in ci_candidates if c["classification"] == "IDENTICAL"]
            ci_identical.sort(key=lambda c: (
                -(c["mit_spec_score"] if c["mit_spec_score"] is not None else -1),
                c["offtarget_count"] if c["offtarget_count"] is not None else 9999,
                abs(c["tss_rel_start"]),
            ))
            ci_top = ci_identical[:TOP_N]
        else:
            fallback = sorted(
                (r for r in ci_rows if r["classification"] == "IDENTICAL"),
                key=lambda r: abs(int(r["tss_rel_start"]))
            )[:TOP_N]
            ci_top = [{
                "target_seq": r["ref_spacer"] + r["ref_pam"], "spacer": r["ref_spacer"], "pam": r["ref_pam"],
                "guide_id_crispor": None, "classification": "IDENTICAL",
                "mit_spec_score": None, "cfd_spec_score": None, "offtarget_count": None,
                "pop_offtarget_count": None, "top_offtargets": [], **ci_position(r),
            } for r in fallback]

        report[gene] = {
            "transcript": chosen_tid,
            "cds_len": cds_len,
            "n_exons": len(ref_cds),
            "strand": ref_strand,
            "chrom": ref_gene[0],
            "gene_span": f"{ref_gene[3]}-{ref_gene[4]}",
            "classification_counts": ko_counts,
            "total_guides": len(ko_rows),
            "variant_affected_guides": ko_variant_rows,
            "top_ko_candidates": ko_top,
            "crispor_available": ko_crispor_available,
            "crispri": {
                "total_guides": len(ci_rows),
                "classification_counts": ci_counts,
                "variant_affected_guides": ci_variant_rows,
                "top_candidates": ci_top,
                "crispor_available": ci_crispor_available,
            },
        }
        print(f"{gene}: KO {len(ko_rows)} guides (crispor={ko_crispor_available}, top={len(ko_top)}), "
              f"CRISPRi {len(ci_rows)} guides (crispor={ci_crispor_available}, top={len(ci_top)})", file=sys.stderr)

    out_path = OUT_DIR / "report_data.json"
    with open(out_path, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"\nWritten: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
