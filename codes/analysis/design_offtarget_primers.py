#!/usr/bin/env python3
"""
Design PCR primer pairs around CRISPR on-target/off-target cut sites, for
wet-lab gel electrophoresis screening and targeted deep amplicon
sequencing (much higher depth than the WGS-based off-target analysis
already in this repo).

Input: a "combined_offtargets.csv"-shaped CSV (produced by
codes/CRISPResso/combine_offtargets.py) - columns chromosome, start, end,
strand, offtarget_seq, mismatches, mit_score, cfd_score, locus, source,
found_by_both. The row with mismatches==0 (source=CasOFFinder) is treated
as the on-target site; every other row is an off-target site, numbered
off_target_1..N in CSV row order.

Design is done against the REFERENCE genome ($REF, v1 female or v2 male
per --ref-version) - NOT the pseudogenome directly - matching this
project's established pattern (ko_guide_scan.py enumerates candidates on
the reference, then classifies against the pseudogenome). This keeps
primer coordinates anchored to the same genome used for the whole-genome
off-target discovery (Cas-OFFinder/CRISPOR both ran against the
reference), and works immediately for both v1 and v2 without needing the
v2 pseudogenome (not yet built as of this writing).

Reused verbatim from ko_guide_scan.py: faidx_seq(), revcomp(), align_cs(),
parse_cs_variants(), REF_FASTA_V1/V2 - no reimplementation.

Two external tools, both via `module load emboss/6.6.0`:
  - eprimer3: EMBOSS wrapper around the classic (boulder-IO) Primer3
    engine. IMPORTANT: needs the external primer3_core binary, which the
    module does NOT ship - see codes/analysis/setup_primer3.sh (one-time
    conda install) and the EMBOSS_PRIMER3_CORE env var this script expects
    to already be exported (done by run_offtarget_primer_design.sh).
  - primersearch: EMBOSS's in-silico PCR/specificity checker, no extra
    binary needed - used to confirm each candidate pair amplifies only
    the intended site genome-wide.

Reusable for future genes: --gene / --sites-csv / --ref-version are
sufficient - the script never needs a guide sequence, only the site
coordinates already produced by the existing off-target discovery
pipeline (same CSV shape for any gene).

Usage:
    python3 design_offtarget_primers.py --gene bdnf \
        --sites-csv crispresso/offtargets/combined/combined_offtargets.csv \
        --ref-version v1
"""
import argparse
import csv
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from ko_guide_scan import (  # noqa: E402
    REF_FASTA_V1, REF_FASTA_V2, faidx_seq, revcomp, align_cs, parse_cs_variants,
)

PROJECT_DIR = Path("/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data")
OUT_DIR = PROJECT_DIR / "analysis" / "offtarget_primers"

EMBOSS_PRIMER3_CORE = os.environ.get(
    "EMBOSS_PRIMER3_CORE",
    "/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso/envs/primer3_env/bin/primer3_core",
)

REF_BY_VERSION = {"v1": REF_FASTA_V1, "v2": REF_FASTA_V2}
PSEUDOGENOME_BY_VERSION = {
    "v1": PROJECT_DIR / "reference/pseudogenome/colombian_pseudogenome.fna",
    "v2": PROJECT_DIR / "reference/pseudogenome_v2/colombian_pseudogenome.fna",
}
# bcftools consensus's own .chain file (produced alongside the pseudogenome
# by make_pseudogenome.sh) - the authoritative ref->pseudogenome coordinate
# mapping. A naive "same nominal coordinates +/- small buffer" approach
# does NOT work here: the pseudogenome chromosome carrying bdnf is 5314bp
# longer than the reference's (net indel insertions), and the drift
# already accumulated by the bdnf locus (~15.9Mb in) was measured at
# ~2838bp - far past any reasonable fixed buffer. Confirmed via CrossMap
# directly on this cluster before wiring it in here.
CHAIN_BY_VERSION = {
    "v1": PROJECT_DIR / "reference/pseudogenome/colombian_pseudogenome.chain",
    "v2": PROJECT_DIR / "reference/pseudogenome_v2/colombian_pseudogenome.chain",
}
CROSSMAP_BIN = "/hpcfs/home/ing_civil/da.martinez33/miniconda3_crispresso/envs/crossmap_env/bin/CrossMap"


def liftover_region(chain_file, chrom, start1, end1, tmp_path):
    """1-based inclusive (start1, end1) -> 1-based inclusive pseudogenome
    coords via CrossMap bed (BED is 0-based half-open, hence the +/-1).
    Returns (chrom, start1, end1) or None if the region didn't lift over
    (e.g. falls in a gap/rearranged block - not expected here but handled
    defensively)."""
    bed_in = tmp_path.with_suffix(".liftin.bed")
    bed_out = tmp_path.with_suffix(".liftout.bed")
    with open(bed_in, "w") as fh:
        fh.write(f"{chrom}\t{start1 - 1}\t{end1}\tregion\n")
    cmd = [CROSSMAP_BIN, "bed", str(chain_file), str(bed_in), str(bed_out)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not bed_out.exists():
        return None
    with open(bed_out) as fh:
        line = fh.readline().strip()
    if not line or "Fail" in line:
        return None
    parts = line.split("\t")
    # CrossMap bed output (4-col input): chrom start end name -> chrom start end name
    # i.e. the lifted coordinates are the LAST 4 fields.
    out_chrom, out_start, out_end = parts[-4], int(parts[-3]), int(parts[-2])
    return out_chrom, out_start + 1, out_end


def fai_length(fasta_path, chrom):
    fai = Path(str(fasta_path) + ".fai")
    with open(fai) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == chrom:
                return int(parts[1])
    raise ValueError(f"chromosome '{chrom}' not found in {fai}")


def compute_cut_site(start, end, strand):
    """SpCas9 blunt cut, ~3bp upstream (PAM-proximal side) of the 20bp
    spacer's PAM-adjacent end. Not meant to be indel-breakpoint-precise -
    only used to center the primer-design window."""
    return end - 5 if strand == "+" else start + 5


def clamp_window(chrom_len, start, end):
    return max(1, start), min(chrom_len, end)


def count_ambiguous_bases(seq):
    """primer3_core (the classic boulder-IO engine eprimer3 shells out to)
    hard-rejects the ENTIRE input sequence - with exit code 0, no less,
    see run_eprimer3() - if it contains any IUPAC ambiguity code (only
    A/C/G/T/N tolerated, lowercase soft-masking is fine). The 2014
    short-read reference assembly genuinely has these scattered at
    unresolved-heterozygous-site positions (same root cause already
    documented for CRISPOR failures on agap3/grin1a/gria1a - see
    CLAUDE.md). Returns {code: count} for any of K/M/R/S/W/Y/B/D/H/V found
    (case-insensitive)."""
    counts = {}
    for c in seq.upper():
        if c not in "ACGTN":
            counts[c] = counts.get(c, 0) + 1
    return counts


def run_eprimer3(window_fa, out_file, excl_start, excl_end, args):
    env = os.environ.copy()
    env["EMBOSS_PRIMER3_CORE"] = EMBOSS_PRIMER3_CORE
    cmd = [
        "eprimer3", "-sequence", str(window_fa), "-outfile", str(out_file),
        "-task", "1", "-numreturn", str(args.num_return),
        "-excludedregion", f"{excl_start},{excl_end}",
        "-optsize", "20", "-minsize", "18", "-maxsize", "25",
        "-opttm", str(args.opt_tm), "-mintm", str(args.min_tm), "-maxtm", str(args.max_tm),
        "-maxdifftm", "3",
        "-ogcpercent", "50", "-mingc", str(args.min_gc), "-maxgc", str(args.max_gc),
        "-psizeopt", str((args.product_min + args.product_max) // 2),
        "-prange", f"{args.product_min}-{args.product_max}",
        "-auto",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    # eprimer3 can print "Error: ..." / "Died: ..." from the underlying
    # primer3_core and STILL exit 0 - confirmed on this cluster (an
    # IUPAC-ambiguity-code input silently produced "0 candidates" with no
    # error surfaced, indistinguishable from a legitimate empty search
    # result, until stdout/stderr text was checked directly). Don't trust
    # returncode alone.
    combined = (r.stdout or "") + (r.stderr or "")
    if r.returncode != 0 or "Error:" in combined or "Died:" in combined:
        print(f"WARNING: eprimer3 failed for {window_fa}: {combined.strip()[-500:]}", file=sys.stderr)
        return None  # distinct from [] (ran cleanly, found nothing)
    return parse_eprimer3(out_file)


# Matches one full numbered result entry (FORWARD then REVERSE, possibly
# separated by blank lines) - non-greedy + DOTALL so it doesn't bleed into
# the next numbered entry.
EPRIMER3_ENTRY_RE = re.compile(
    r"\d+\s+PRODUCT SIZE:\s*(\d+).*?"
    r"FORWARD PRIMER\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([ACGTN]+).*?"
    r"REVERSE PRIMER\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([ACGTN]+)",
    re.DOTALL,
)


def parse_eprimer3(path):
    text = Path(path).read_text()
    results = []
    for m in EPRIMER3_ENTRY_RE.finditer(text):
        (product_size, fwd_start, fwd_len, fwd_tm, fwd_gc, fwd_seq,
         rev_start, rev_len, rev_tm, rev_gc, rev_seq) = m.groups()
        results.append({
            "product_size": int(product_size),
            "fwd_local_start": int(fwd_start), "fwd_len": int(fwd_len),
            "fwd_tm": float(fwd_tm), "fwd_gc": float(fwd_gc), "fwd_seq": fwd_seq,
            "rev_local_start": int(rev_start), "rev_len": int(rev_len),
            "rev_tm": float(rev_tm), "rev_gc": float(rev_gc), "rev_seq": rev_seq,
        })
    return results


def run_primersearch(pairs, ref_fasta, mismatch_pct, tmp_dir):
    """pairs: list of (name, fwd_seq, rev_seq). Returns {name: [(amplimer_idx,
    length), ...]} - one entry per Amplimer found for that primer name."""
    pairs_file = tmp_dir / "pairs.txt"
    with open(pairs_file, "w") as fh:
        for name, fwd, rev in pairs:
            fh.write(f"{name}\t{fwd}\t{rev}\n")
    out_file = tmp_dir / "all.primersearch"
    cmd = [
        "primersearch", "-seqall", str(ref_fasta), "-infile", str(pairs_file),
        "-mismatchpercent", str(mismatch_pct), "-outfile", str(out_file), "-auto",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"ERROR: primersearch failed: {r.stderr.strip()[-800:]}")
    return parse_primersearch(out_file)


def parse_primersearch(path):
    text = Path(path).read_text()
    hits = {}
    for block in re.split(r"\n(?=Primer name )", text):
        m = re.match(r"Primer name (\S+)", block)
        if not m:
            continue
        name = m.group(1)
        amplimers = re.findall(r"Amplimer length:\s*(\d+)\s*bp", block)
        hits[name] = [int(x) for x in amplimers]
    return hits


def load_sites(csv_path):
    with open(csv_path) as fh:
        rows = list(csv.DictReader(fh))
    on_target_idx = next((i for i, r in enumerate(rows) if int(r["mismatches"]) == 0), None)
    if on_target_idx is None:
        sys.exit(f"ERROR: no mismatches==0 row (on-target) found in {csv_path}")
    sites = []
    off_i = 0
    for i, r in enumerate(rows):
        if i == on_target_idx:
            label, site_type = "on_target", "on_target"
        else:
            off_i += 1
            label, site_type = f"off_target_{off_i}", "off_target"
        sites.append({
            "label": label, "site_type": site_type,
            "chrom": r["chromosome"], "start": int(r["start"]), "end": int(r["end"]),
            "strand": r["strand"], "mismatches": r["mismatches"],
            "mit_score": r.get("mit_score", ""), "cfd_score": r.get("cfd_score", ""),
            "locus": r.get("locus", ""),
        })
    return sites


def check_population_variants(ref_fasta, pop_fasta, chain_file, chrom, ref_win_start, ref_win_end,
                               pop_buffer, tmp_prefix):
    """Locate the pseudogenome window via the reference->pseudogenome
    .chain file (bcftools consensus's own liftover mapping - NOT a fixed
    coordinate buffer, which is unsafe: confirmed on this cluster that the
    bdnf locus alone has ~2838bp of accumulated indel drift by ~15.9Mb into
    its chromosome, and the whole chromosome is 5314bp longer in the
    pseudogenome). Then align the two windows to find real variants.
    Returns (variants_by_pos, ok) - variants_by_pos keys are 1-based
    positions LOCAL to the reference window (1 = ref_win_start), matching
    primer local coordinates directly. ok=False if the pseudogenome/chain
    is missing for this ref_version, or the region didn't lift over/align -
    callers should treat that as "unknown", not "no variant"."""
    if not pop_fasta.exists() or not chain_file.exists():
        return {}, False
    lifted = liftover_region(chain_file, chrom, ref_win_start, ref_win_end, tmp_prefix)
    if lifted is None:
        print(f"WARNING: liftover failed for {chrom}:{ref_win_start}-{ref_win_end}", file=sys.stderr)
        return {}, False
    pop_chrom, pop_lift_start, pop_lift_end = lifted
    ref_seq = faidx_seq(ref_fasta, chrom, ref_win_start, ref_win_end)
    pop_win_start = max(1, pop_lift_start - pop_buffer)
    pop_win_end = pop_lift_end + pop_buffer
    pop_seq = faidx_seq(pop_fasta, pop_chrom, pop_win_start, pop_win_end)
    pos, cigar, cs_tag, flag = align_cs(ref_seq, pop_seq, str(tmp_prefix))
    if not cs_tag:
        print(f"WARNING: reference/pseudogenome window did not align for {chrom}:{ref_win_start}-{ref_win_end} "
              f"(lifted to {pop_chrom}:{pop_win_start}-{pop_win_end})", file=sys.stderr)
        return {}, False
    variants, _ = parse_cs_variants(cs_tag, pos)
    return {v[0]: v for v in variants}, True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--sites-csv", required=True, type=Path)
    ap.add_argument("--ref-version", choices=["v1", "v2"], default="v1")
    ap.add_argument("--window-halfsize", type=int, default=500)
    ap.add_argument("--exclude-halfsize", type=int, default=75)
    ap.add_argument("--product-min", type=int, default=200)
    ap.add_argument("--product-max", type=int, default=400)
    ap.add_argument("--opt-tm", type=float, default=60)
    ap.add_argument("--min-tm", type=float, default=58)
    ap.add_argument("--max-tm", type=float, default=62)
    ap.add_argument("--min-gc", type=float, default=40)
    ap.add_argument("--max-gc", type=float, default=60)
    ap.add_argument("--num-return", type=int, default=5)
    ap.add_argument("--specificity-mismatch-pct", type=float, default=10)
    ap.add_argument("--pop-window-buffer", type=int, default=100,
                     help="extra flanking bp on each side of the pseudogenome window vs the reference window")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--no-population-check", action="store_true")
    args = ap.parse_args()

    ref_fasta = REF_BY_VERSION[args.ref_version]
    pop_fasta = PSEUDOGENOME_BY_VERSION[args.ref_version]
    chain_file = CHAIN_BY_VERSION[args.ref_version]
    out_dir = args.out_dir
    raw_dir = out_dir / "raw" / f"{args.gene}_{args.ref_version}"
    raw_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = Path(f"/tmp/design_offtarget_primers_{args.gene}_{args.ref_version}")
    tmp_dir.mkdir(parents=True, exist_ok=True)

    sites = load_sites(args.sites_csv)
    print(f"=== {args.gene} ({args.ref_version}): {len(sites)} sites from {args.sites_csv} ===")

    site_candidates = {}   # label -> list of candidate dicts (with genomic coords, no specificity/variant info yet)
    all_pairs_for_search = []  # (unique_name, fwd_seq, rev_seq)

    site_design_status = {}  # label -> (status, note) for sites with zero usable candidates
    for site in sites:
        chrom, strand = site["chrom"], site["strand"]
        cut_pos = compute_cut_site(site["start"], site["end"], strand)
        chrom_len = fai_length(ref_fasta, chrom)
        win_start, win_end = clamp_window(chrom_len, cut_pos - args.window_halfsize, cut_pos + args.window_halfsize)

        window_seq = faidx_seq(ref_fasta, chrom, win_start, win_end)
        window_fa = tmp_dir / f"{site['label']}.fa"
        with open(window_fa, "w") as fh:
            fh.write(f">{args.gene}_{site['label']}\n{window_seq}\n")  # colon-free header, see docstring

        cut_local = cut_pos - win_start + 1
        excl_start = max(1, cut_local - args.exclude_halfsize)
        excl_end = min(len(window_seq), cut_local + args.exclude_halfsize)

        # Pre-flight check: primer3_core hard-rejects any IUPAC ambiguity
        # code in the WHOLE input (exit 0, silently 0 candidates otherwise -
        # see run_eprimer3()/count_ambiguous_bases()) - check before even
        # calling it, so the reason is reported precisely rather than
        # looking identical to "searched and found nothing".
        amb = count_ambiguous_bases(window_seq)
        if amb:
            detail = ", ".join(f"{k}x{v}" for k, v in sorted(amb.items()))
            note = f"reference window has {sum(amb.values())} IUPAC ambiguity code bases ({detail}) - primer3_core cannot process this input at all"
            print(f"  {site['label']} ({chrom}:{cut_pos}, {strand}): SKIPPED - {note}")
            site_design_status[site["label"]] = ("INPUT_ERROR_AMBIGUOUS_BASES", note)
            continue

        eprimer3_out = raw_dir / f"{site['label']}.eprimer3"
        candidates = run_eprimer3(window_fa, eprimer3_out, excl_start, excl_end, args)
        if candidates is None:
            note = "eprimer3/primer3_core reported an error - see raw output"
            print(f"  {site['label']} ({chrom}:{cut_pos}, {strand}): FAILED - {note}")
            site_design_status[site["label"]] = ("INPUT_ERROR_OTHER", note)
            continue
        print(f"  {site['label']} ({chrom}:{cut_pos}, {strand}): {len(candidates)} candidate pairs")
        if not candidates:
            site_design_status[site["label"]] = (
                "NO_CANDIDATES_FOUND",
                f"eprimer3 ran cleanly but found no primer pair meeting the constraints "
                f"(Tm {args.min_tm}-{args.max_tm}, GC {args.min_gc}-{args.max_gc}%, "
                f"product {args.product_min}-{args.product_max}bp) in this window",
            )

        for rank, c in enumerate(candidates, start=1):
            fwd_genomic_start = win_start + c["fwd_local_start"] - 1
            fwd_genomic_end = fwd_genomic_start + c["fwd_len"] - 1
            rev_genomic_start = win_start + c["rev_local_start"] - 1
            rev_genomic_end = rev_genomic_start + c["rev_len"] - 1
            pair_name = f"{site['label']}_r{rank}"
            row = {
                "gene": args.gene, "site_label": site["label"], "site_type": site["site_type"],
                "chromosome": chrom, "cut_pos": cut_pos, "mismatches_vs_guide": site["mismatches"],
                "mit_score": site["mit_score"], "cfd_score": site["cfd_score"], "locus": site["locus"],
                "window_start": win_start, "window_end": win_end,
                "primer_rank": rank, "product_size": c["product_size"],
                "forward_seq": c["fwd_seq"], "forward_tm": c["fwd_tm"], "forward_gc": c["fwd_gc"],
                "forward_genomic_start": fwd_genomic_start, "forward_genomic_end": fwd_genomic_end,
                "reverse_seq": c["rev_seq"], "reverse_tm": c["rev_tm"], "reverse_gc": c["rev_gc"],
                "reverse_genomic_start": rev_genomic_start, "reverse_genomic_end": rev_genomic_end,
                "design_status": "OK", "design_note": "",
                "_pair_name": pair_name, "_win_start": win_start, "_win_end": win_end, "_chrom": chrom,
            }
            site_candidates.setdefault(site["label"], []).append(row)
            all_pairs_for_search.append((pair_name, c["fwd_seq"], c["rev_seq"]))

    # ---- specificity check (one batched primersearch call for everything) ----
    print(f"\n=== Running primersearch ({len(all_pairs_for_search)} pairs vs {ref_fasta.name}) ===")
    hits_by_name = run_primersearch(all_pairs_for_search, ref_fasta, args.specificity_mismatch_pct, tmp_dir)

    # ---- population-variant check (one alignment per site, reused across its candidate ranks) ----
    variants_by_site = {}
    pop_check_ok = {}
    if not args.no_population_check:
        print(f"\n=== Checking primer footprints against pseudogenome ({pop_fasta}) ===")
        for site in sites:
            rows = site_candidates.get(site["label"], [])
            if not rows:
                continue
            win_start, win_end, chrom = rows[0]["_win_start"], rows[0]["_win_end"], rows[0]["_chrom"]
            variants, ok = check_population_variants(
                ref_fasta, pop_fasta, chain_file, chrom, win_start, win_end,
                args.pop_window_buffer, tmp_dir / f"popcheck_{site['label']}"
            )
            variants_by_site[site["label"]] = variants
            pop_check_ok[site["label"]] = ok
        n_ok = sum(pop_check_ok.values())
        print(f"  population check available for {n_ok}/{len(sites)} sites "
              f"({'pseudogenome not found: ' + str(pop_fasta) if n_ok == 0 else ''})")

    # ---- assemble final rows ----
    final_rows = []
    for site in sites:
        variants = variants_by_site.get(site["label"], {})
        checked = pop_check_ok.get(site["label"], False)
        rows_for_site = site_candidates.get(site["label"], [])
        if not rows_for_site and site["label"] in site_design_status:
            status, note = site_design_status[site["label"]]
            final_rows.append({
                "gene": args.gene, "site_label": site["label"], "site_type": site["site_type"],
                "chromosome": site["chrom"], "cut_pos": compute_cut_site(site["start"], site["end"], site["strand"]),
                "mismatches_vs_guide": site["mismatches"], "mit_score": site["mit_score"],
                "cfd_score": site["cfd_score"], "locus": site["locus"],
                "window_start": "", "window_end": "", "primer_rank": "", "product_size": "",
                "forward_seq": "", "forward_tm": "", "forward_gc": "",
                "forward_genomic_start": "", "forward_genomic_end": "",
                "reverse_seq": "", "reverse_tm": "", "reverse_gc": "",
                "reverse_genomic_start": "", "reverse_genomic_end": "",
                "design_status": status, "design_note": note,
                "specificity_amplimer_count": "", "specificity_status": "", "specificity_extra_hits": "",
                "primer_variant_flag": "", "primer_variant_note": "",
            })
            continue
        for row in rows_for_site:
            amplimers = hits_by_name.get(row["_pair_name"], [])
            n_amp = len(amplimers)
            if n_amp == 1:
                spec_status = "SPECIFIC"
            elif n_amp == 0:
                spec_status = "ERROR_NO_HIT"
            else:
                spec_status = "WARN_MULTI_HIT"
            row["specificity_amplimer_count"] = n_amp
            row["specificity_status"] = spec_status
            row["specificity_extra_hits"] = "; ".join(f"{l}bp" for l in amplimers[1:]) if n_amp > 1 else ""

            if not args.no_population_check:
                win_start = row["_win_start"]
                fwd_local = range(row["forward_genomic_start"] - win_start + 1, row["forward_genomic_end"] - win_start + 2)
                rev_local = range(row["reverse_genomic_start"] - win_start + 1, row["reverse_genomic_end"] - win_start + 2)
                fwd_hit = [p for p in fwd_local if p in variants]
                rev_hit = [p for p in rev_local if p in variants]
                if not checked:
                    row["primer_variant_flag"] = "UNKNOWN"
                    row["primer_variant_note"] = "pseudogenome not available for this ref-version yet"
                elif fwd_hit and rev_hit:
                    row["primer_variant_flag"] = "BOTH"
                elif fwd_hit:
                    row["primer_variant_flag"] = "FORWARD"
                elif rev_hit:
                    row["primer_variant_flag"] = "REVERSE"
                else:
                    row["primer_variant_flag"] = "NONE"
                notes = [f"{variants[p][1]} at ref-window pos {p} ({variants[p][2]}->{variants[p][3]})"
                         for p in (fwd_hit + rev_hit)]
                row["primer_variant_note"] = "; ".join(notes)
            else:
                row["primer_variant_flag"] = "SKIPPED"
                row["primer_variant_note"] = ""

            for k in ("_pair_name", "_win_start", "_win_end", "_chrom"):
                row.pop(k, None)
            final_rows.append(row)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"{args.gene}_{args.ref_version}_primers.csv"
    fieldnames = list(final_rows[0].keys()) if final_rows else []
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(final_rows)

    print(f"\n=== Summary ===")
    print(f"Total candidate pairs: {len(final_rows)}")
    for site in sites:
        n = len(site_candidates.get(site["label"], []))
        if site["label"] in site_design_status:
            status, note = site_design_status[site["label"]]
            print(f"  {site['label']}: {status} - {note}")
            continue
        n_spec = sum(1 for r in final_rows if r["site_label"] == site["label"] and r["specificity_status"] == "SPECIFIC")
        n_clean = sum(1 for r in final_rows if r["site_label"] == site["label"] and r["primer_variant_flag"] == "NONE")
        print(f"  {site['label']}: {n} candidates, {n_spec} specific, {n_clean} without population variants in primer footprint")
    print(f"\nWritten: {out_csv}")


if __name__ == "__main__":
    main()
