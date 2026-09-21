#!/usr/bin/env python3
"""
Verify existing RT-qPCR primer pairs against the guppy genome(s) and the
Colombian pseudogenome - formalizes the ad-hoc checklist worked out by
hand for bdnf/actb2/rpl13a/myosin (see CLAUDE.md, "RT-qPCR Primer
Verification"), so it can be re-run consistently for any future primer
set without repeating the same gaps that were found and fixed there:

  1. LOCATE by BLAST (blastn-short) against the whole reference genome,
     not a hand-picked gene list - a primer's true target gene is found
     from where it actually binds, not guessed from a gene-symbol regex
     (that regex missed guppy's un-symbolled LOC######## myosin heavy
     chain genes entirely - BLAST does not have this blind spot).
  2. ANNOTATE by spatial overlap with the GFF (gene + exon features at
     the hit position), regardless of gene symbol/naming - resolves the
     same LOC-gene blind spot from the annotation side too.
  3. REALIGN with an indel-aware pairwise aligner (EMBOSS `water`) once a
     candidate gene/transcript is found, never a fixed-length Hamming
     scan - a scan reports a single 1bp indel as "3 mismatches" and can
     make a near-perfect primer look unrelated to its real target.
  4. CROSS-CHECK every requested genome version (v1/v2) independently -
     a primer's binding site can be well-assembled in one version and
     fragmentary/absent in another (seen for v1's copy of a myosin
     paralog).
  5. POPULATION CHECK via exact liftover (.chain + CrossMap) of the
     matched genomic footprint to the Colombian pseudogenome - each ref
     version checked against its own pseudogenome (v1 against
     reference/pseudogenome/, v2 against reference/pseudogenome_v2/), for
     exon-only positions - a primer that matches genomic DNA is not
     automatically an mRNA primer; only the sequence within a real exon
     reaches the mRNA/cDNA RT-qPCR actually amplifies.

Reused verbatim (no reimplementation): REF_FASTA_V1/V2, REF_GFF_V1/V2,
faidx_seq(), revcomp() from ko_guide_scan.py; REF_BY_VERSION,
PSEUDOGENOME_BY_VERSION, CHAIN_BY_VERSION, CROSSMAP_BIN,
liftover_region() from design_offtarget_primers.py.

Usage:
    python3 verify_rtqpcr_primers.py \
        --primers-csv my_rtqpcr_primers.csv \
        --ref-versions v1,v2 \
        --population-check \
        --out-dir analysis/rtqpcr_verification

primers.csv columns: pair_name,forward_seq,reverse_seq
"""
import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from ko_guide_scan import REF_FASTA_V1, REF_FASTA_V2, REF_GFF_V1, REF_GFF_V2, faidx_seq, revcomp  # noqa: E402
from design_offtarget_primers import (  # noqa: E402
    REF_BY_VERSION, PSEUDOGENOME_BY_VERSION, CHAIN_BY_VERSION, CROSSMAP_BIN, liftover_region,
)

PROJECT_DIR = Path("/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data")
OUT_DIR = PROJECT_DIR / "analysis" / "rtqpcr_verification"
GFF_BY_VERSION = {"v1": REF_GFF_V1, "v2": REF_GFF_V2}
BLAST_DB_DIR = {
    "v1": PROJECT_DIR / "reference" / "blast_db_v1",
    "v2": PROJECT_DIR / "reference" / "blast_db_v2",
}

MIN_IDENTITY_PCT = 70.0     # BLAST hit quality floor - below this, treat as noise
MIN_LEN_FRAC = 0.6          # hit must cover at least this fraction of the primer
MAX_GENE_CANDIDATES = 30    # distinct genes realigned via water per primer per genome


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


# ── Step 0: BLAST db (built once, reused across runs/primers) ──────────────

def ensure_blast_db(ref_version):
    db_dir = BLAST_DB_DIR[ref_version]
    db_prefix = db_dir / Path(REF_BY_VERSION[ref_version]).stem
    if (db_prefix.with_suffix(".nhr")).exists() or Path(str(db_prefix) + ".nhr").exists():
        return db_prefix
    db_dir.mkdir(parents=True, exist_ok=True)
    print(f"  (building BLAST db for {ref_version} - one-time, ~10-30s)", file=sys.stderr)
    r = run(["makeblastdb", "-in", str(REF_BY_VERSION[ref_version]), "-dbtype", "nucl",
             "-out", str(db_prefix), "-title", ref_version])
    if r.returncode != 0:
        sys.exit(f"ERROR: makeblastdb failed for {ref_version}:\n{r.stderr}")
    return db_prefix


# ── Step 1: locate by BLAST ─────────────────────────────────────────────────

def blast_locate_all(primer_seq, db_prefix):
    """Returns ALL hits passing the length floor (by bitscore desc) - no
    hit-count cap here. Important for primers that cross a real exon-exon
    junction: their best RAW GENOMIC hit is only the larger of the two
    exonic sub-fragments (often well under full primer length), which can
    easily rank below dozens of short, coincidental, non-genic matches
    elsewhere in a large genome. If hits were capped before checking gene
    overlap, the true (but weak-looking) hit can be discarded before ever
    getting a chance - verified against a real case (rpl_13a_R): capping
    at the top 25 raw hits silently dropped its true rpl13a hit in favor
    of an unrelated gene's coincidental 16bp match. Capping instead
    happens per-GENE, after annotating (see verify_primer)."""
    r = run(["blastn", "-task", "blastn-short", "-db", str(db_prefix),
             "-word_size", "7", "-evalue", "1000", "-perc_identity", str(MIN_IDENTITY_PCT),
             "-outfmt", "6 sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand",
             "-max_target_seqs", "10"],
            input=f">q\n{primer_seq}\n")
    if r.returncode != 0 or not r.stdout.strip():
        return []
    min_len = len(primer_seq) * MIN_LEN_FRAC
    hits = []
    for line in r.stdout.strip().split("\n"):
        f = line.split("\t")
        sseqid, pident, length = f[0], float(f[1]), int(f[2])
        if length < min_len:
            continue
        sstart, send = int(f[7]), int(f[8])
        hits.append({
            "chrom": sseqid, "pident": pident, "length": length,
            "mismatch": int(f[3]), "gapopen": int(f[4]),
            "sstart": min(sstart, send), "send": max(sstart, send),
            "sstrand": f[11], "bitscore": float(f[10]),
        })
    hits.sort(key=lambda h: -h["bitscore"])
    return hits


# ── Step 2: annotate by spatial GFF overlap (symbol-agnostic) ──────────────

# A whole-genome GFF here is 1.3-1.6M lines / 350-470MB (v1/v2). Both
# genes_overlapping() and representative_transcript() used to re-open and
# linearly scan the ENTIRE file on every call - and they're called once per
# BLAST hit and once per candidate gene respectively (up to MAX_GENE_CANDIDATES
# per primer). For a primer with many scattered BLAST hits (e.g. no clean
# genomic match, like a junction-spanning primer) this meant hundreds of
# full-file re-reads per primer - the actual bottleneck (hours), not the
# `water` alignments. Fix: parse the GFF into an in-memory index ONCE per
# file (cached for the life of the process, i.e. once per ref version for
# the whole run), then serve every lookup from memory.
_GFF_INDEX_CACHE = {}


def load_gff_index(gff_path):
    gff_path = str(gff_path)
    if gff_path in _GFF_INDEX_CACHE:
        return _GFF_INDEX_CACHE[gff_path]
    print(f"  (indexing {gff_path} - one-time per run)", file=sys.stderr)
    genes_by_chrom = {}
    mrna_by_gene = {}   # (chrom, gene_symbol) -> [(transcript_id, strand), ...]
    exons_by_tx = {}    # transcript_id -> [(start, end), ...]
    with open(gff_path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, ftype, attrs = f[0], f[2], f[8]
            if ftype == "gene":
                gene_m = re.search(r"gene=([^;]+)", attrs)
                desc_m = re.search(r"description=([^;]+)", attrs)
                genes_by_chrom.setdefault(chrom, []).append({
                    "gene": gene_m.group(1) if gene_m else "?",
                    "description": (desc_m.group(1).replace("%2C", ",") if desc_m else ""),
                    "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                })
            elif ftype == "mRNA":
                gene_m = re.search(r"gene=([^;]+)", attrs)
                id_m = re.search(r"ID=([^;]+)", attrs)
                if gene_m and id_m:
                    mrna_by_gene.setdefault((chrom, gene_m.group(1)), []).append((id_m.group(1), f[6]))
            elif ftype == "exon":
                parent_m = re.search(r"Parent=([^;]+)", attrs)
                if parent_m:
                    exons_by_tx.setdefault(parent_m.group(1), []).append((int(f[3]), int(f[4])))
    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(key=lambda g: g["start"])
    idx = {"genes_by_chrom": genes_by_chrom, "mrna_by_gene": mrna_by_gene, "exons_by_tx": exons_by_tx}
    _GFF_INDEX_CACHE[gff_path] = idx
    return idx


def genes_overlapping(gff_path, chrom, start, end):
    """Any gene feature whose interval overlaps [start,end] on this chrom,
    regardless of symbol (catches LOC######## genes a symbol-regex would
    miss). Returns list of dicts."""
    idx = load_gff_index(gff_path)
    return [g for g in idx["genes_by_chrom"].get(chrom, [])
            if not (g["start"] > end or g["end"] < start)]


def representative_transcript(gff_path, gene_symbol, chrom):
    """Longest-by-summed-exon-length mRNA for this gene (matches the
    'longest CDS as representative isoform' convention already used in
    ko_guide_scan.py). Returns (transcript_id, exon_list[(start,end)], strand)
    or None."""
    idx = load_gff_index(gff_path)
    mrna_list = idx["mrna_by_gene"].get((chrom, gene_symbol), [])
    if not mrna_list:
        return None
    best_tid, best_len, best_strand = None, -1, None
    for tid, strand in mrna_list:
        exons = idx["exons_by_tx"].get(tid, [])
        total = sum(e - s + 1 for s, e in exons)
        if total > best_len:
            best_tid, best_len, best_strand = tid, total, strand
    if best_tid is None or best_len <= 0:
        return None
    return best_tid, sorted(idx["exons_by_tx"][best_tid]), best_strand


def build_spliced_mrna(fasta_path, chrom, exons, strand):
    segs = [faidx_seq(fasta_path, chrom, s, e) for s, e in exons]
    if strand == "-":
        segs = [revcomp(s) for s in reversed(segs)]
    return "".join(segs)


def exon_membership(exons, pos_start, pos_end):
    """True if [pos_start,pos_end] falls entirely within a single exon."""
    for s, e in exons:
        if s <= pos_start and pos_end <= e:
            return True, (s, e)
    return False, None


# ── Step 3: indel-aware realignment (EMBOSS water) ──────────────────────────

def water_align(seq_a, seq_b, tmp_prefix):
    fa_a, fa_b, out = f"{tmp_prefix}.a.fa", f"{tmp_prefix}.b.fa", f"{tmp_prefix}.water"
    Path(fa_a).write_text(f">a\n{seq_a}\n")
    Path(fa_b).write_text(f">b\n{seq_b}\n")
    r = run(["water", "-asequence", fa_a, "-bsequence", fa_b,
             "-gapopen", "10", "-gapextend", "0.5", "-outfile", out, "-auto"])
    if r.returncode != 0 or not Path(out).exists():
        return None
    text = Path(out).read_text()
    ident_m = re.search(r"# Identity:\s*(\d+)/(\d+)\s*\(([\d.]+)%\)", text)
    gaps_m = re.search(r"# Gaps:\s*(\d+)/(\d+)", text)
    if not ident_m:
        return None
    return {
        "identity_n": int(ident_m.group(1)), "identity_len": int(ident_m.group(2)),
        "identity_pct": float(ident_m.group(3)),
        "gaps": int(gaps_m.group(1)) if gaps_m else 0,
    }


# ── Per-primer verification ─────────────────────────────────────────────────

def verify_primer(name, seq, ref_version, tmp_prefix):
    out = {"primer": name, "ref_version": ref_version, "sequence": seq}
    db_prefix = ensure_blast_db(ref_version)
    hits = blast_locate_all(seq, db_prefix)
    if not hits:
        out["status"] = "NO_MATCH"
        return out

    gff = GFF_BY_VERSION[ref_version]

    # Collect one candidate gene per distinct BLAST hit region (smallest
    # overlapping gene wins per hit - see genes_overlapping()), deduplicated
    # by gene symbol, keeping the best (highest-bitscore) hit per gene.
    candidates = {}  # gene_symbol -> (hit, gene_dict)
    for hit in hits:
        genes = genes_overlapping(gff, hit["chrom"], hit["sstart"], hit["send"])
        if not genes:
            continue
        gene = min(genes, key=lambda g: g["end"] - g["start"])
        if gene["gene"] not in candidates:
            candidates[gene["gene"]] = (hit, gene)

    if not candidates:
        best_hit = hits[0]
        out["chrom"], out["hit_start"], out["hit_end"] = best_hit["chrom"], best_hit["sstart"], best_hit["send"]
        out["blast_pident"], out["blast_length"] = best_hit["pident"], best_hit["length"]
        out["status"] = "GENOMIC_MATCH_NO_GENE"
        out["gene"] = ""
        return out

    # Do NOT trust raw genomic BLAST ranking as the final answer: a primer
    # that spans a real exon-exon junction has NO good contiguous genomic
    # hit by definition (only the larger of its two exonic sub-fragments),
    # which can easily score lower than an unrelated, coincidentally-longer
    # genomic match elsewhere. The only reliable arbiter is testing each
    # candidate gene's own spliced mRNA and keeping whichever gives the
    # best (longest, then highest-identity) alignment - so every candidate
    # gene is tried, not just the top-bitscore hit's gene. Capped to the
    # top MAX_GENE_CANDIDATES distinct genes (by their best hit's
    # bitscore) to bound runtime against large gene families (e.g. the
    # myosin heavy chain tandem cluster can surface a dozen paralogs).
    ranked_candidates = sorted(candidates.items(), key=lambda kv: -kv[1][0]["bitscore"])
    best = None
    for gene_symbol, (hit, gene) in ranked_candidates[:MAX_GENE_CANDIDATES]:
        tx = representative_transcript(gff, gene_symbol, hit["chrom"])
        if tx is None:
            continue
        tx_id, exons, strand = tx
        mrna = build_spliced_mrna(REF_BY_VERSION[ref_version], hit["chrom"], exons, strand)
        is_exonic, _ = exon_membership(exons, hit["sstart"], hit["send"])

        if seq in mrna or revcomp(seq) in mrna:
            candidate_result = {
                "gene": gene_symbol, "gene_description": gene["description"],
                "transcript": tx_id, "chrom": hit["chrom"],
                "hit_start": hit["sstart"], "hit_end": hit["send"],
                "blast_pident": hit["pident"], "blast_length": hit["length"],
                "exon_membership": "single_exon" if is_exonic else "spans_junction",
                "mrna_match": "EXACT", "mrna_identity_pct": 100.0, "mrna_gaps": 0,
                "_score": 1.0,
            }
            best = candidate_result
            break  # can't beat an exact match
        # Try BOTH orientations and keep the better one - `water` (Smith-
        # Waterman) almost always returns SOME local alignment even between
        # unrelated sequences, so a plain `a or b` short-circuit on the
        # forward orientation practically never falls through to revcomp.
        # Verified case: miosina_guppy_R against its own true target
        # (LOC145552582) scored only 56.0% forward but 90.0% reverse-
        # complemented (the correct, already-documented result) - the old
        # `or` form would have kept the wrong 56.0% forward alignment.
        #
        # Ranking metric: identity_n / max(identity_len, len(seq)) - the
        # fraction of the PRIMER's own length that is correctly, non-
        # redundantly explained by this alignment. Neither identity_pct nor
        # identity_n alone is safe here - both were tried and both broke on
        # a real case:
        #   - identity_n alone: `water` can stitch a long, heavily gapped
        #     "identity block" (e.g. 19/39=48.7%, 19 gaps) that racks up
        #     more raw identical bases than a short, clean, real match
        #     (e.g. 18/20=90%, 0 gaps) - stk26 (unrelated gene) outranked
        #     the true myosin target this way.
        #   - identity_pct alone: a short coincidental fragment can be a
        #     "perfect" match over just PART of the primer (e.g. 12/12=
        #     100%, covering only 12 of the primer's 19bp) and outrank a
        #     real match covering nearly the whole primer with one indel
        #     (e.g. 19/20=95%) - sgsm2 (unrelated gene) outranked the true
        #     myosin target this way. The max(...) denominator makes this
        #     metric robust to both failure modes: it penalizes alignments
        #     shorter than the primer (denominator floors at len(seq)) AND
        #     alignments padded longer than the primer by excess gaps
        #     (denominator then uses identity_len instead).
        def _coverage_score(a):
            return a["identity_n"] / max(a["identity_len"], len(seq))

        aln_candidates = [a for a in (water_align(seq, mrna, tmp_prefix),
                                       water_align(revcomp(seq), mrna, tmp_prefix)) if a is not None]
        if not aln_candidates:
            continue
        aln = max(aln_candidates, key=_coverage_score)
        candidate_result = {
            "gene": gene_symbol, "gene_description": gene["description"],
            "transcript": tx_id, "chrom": hit["chrom"],
            "hit_start": hit["sstart"], "hit_end": hit["send"],
            "blast_pident": hit["pident"], "blast_length": hit["length"],
            "exon_membership": "single_exon" if is_exonic else "spans_junction",
            "mrna_match": "PARTIAL", "mrna_identity_pct": aln["identity_pct"],
            "mrna_gaps": aln["gaps"],
            "_score": _coverage_score(aln),
        }
        if best is None or candidate_result["_score"] > best["_score"]:
            best = candidate_result

    if best is None:
        top_hit = hits[0]
        out["chrom"], out["hit_start"], out["hit_end"] = top_hit["chrom"], top_hit["sstart"], top_hit["send"]
        out["status"] = "GENE_MATCH_NO_TRANSCRIPT"
        return out

    out.update(best)
    out["status"] = "OK" if best["mrna_match"] == "EXACT" else "WARN_PARTIAL_MATCH"
    return out


def population_check(row, tmp_prefix):
    """Runs against whichever genome version the row was verified under
    (v1 -> reference/pseudogenome/, v2 -> reference/pseudogenome_v2/, both
    built by make_pseudogenome.sh - see CLAUDE.md item #8 migration
    status). Runs for both EXACT and PARTIAL matches - for PARTIAL (e.g.
    the myosin primers, no single paralog is a perfect match),
    `hit_start`/`hit_end` are only the raw BLAST-matched core of the
    primer (shorter than the full primer - see `blast_length` vs the
    primer's own length), not its whole footprint, since the exact
    aligned span within the mRNA isn't recovered by `water_align()`. This
    still answers a real question - whether a Colombian-specific variant
    sits on top of the already-known reference-vs-primer differences in
    that core - just not over the full primer length; flagged explicitly
    in `population_note` so it isn't confused with the exact, whole-primer
    check used for EXACT matches."""
    ref_version = row.get("ref_version")
    if ref_version not in PSEUDOGENOME_BY_VERSION or row.get("mrna_match") not in ("EXACT", "PARTIAL"):
        row["population_status"] = "not_checked"
        return
    partial_coverage_note = ""
    if row["mrna_match"] == "PARTIAL":
        partial_coverage_note = (
            f"covers only the {row['blast_length']}bp exact BLAST core of the "
            f"{len(row['sequence'])}bp primer, not its full footprint; "
        )
    chain = CHAIN_BY_VERSION[ref_version]
    lifted = liftover_region(chain, row["chrom"], row["hit_start"], row["hit_end"], Path(tmp_prefix))
    if lifted is None:
        row["population_status"] = "liftover_failed"
        row["population_note"] = partial_coverage_note.rstrip("; ")
        return
    pg_chrom, pg_start, pg_end = lifted
    ref_seq = faidx_seq(REF_BY_VERSION[ref_version], row["chrom"], row["hit_start"], row["hit_end"])
    pg_seq = faidx_seq(PSEUDOGENOME_BY_VERSION[ref_version], pg_chrom, pg_start, pg_end)
    row["pseudogenome_region"] = f"{pg_chrom}:{pg_start}-{pg_end}"
    row["population_status"] = "IDENTICAL" if ref_seq == pg_seq else "VARIANT_FOUND"
    if ref_seq != pg_seq:
        diffs = [f"{i+1}:{a}>{b}" for i, (a, b) in enumerate(zip(ref_seq, pg_seq)) if a != b]
        row["population_note"] = partial_coverage_note + "; ".join(diffs)
    elif partial_coverage_note:
        row["population_note"] = partial_coverage_note.rstrip("; ")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--primers-csv", required=True, type=Path,
                     help="columns: pair_name,forward_seq,reverse_seq")
    ap.add_argument("--ref-versions", default="v1,v2", help="comma-separated, e.g. v1,v2")
    ap.add_argument("--population-check", action="store_true",
                     help="liftover exact single-exon v1 hits to the pseudogenome and compare")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()

    ref_versions = args.ref_versions.split(",")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = args.out_dir / "tmp"
    tmp_dir.mkdir(exist_ok=True)

    rows = []
    with open(args.primers_csv) as fh:
        for r in csv.DictReader(fh):
            for role, seq in (("F", r["forward_seq"]), ("R", r["reverse_seq"])):
                primer_name = f"{r['pair_name']}_{role}"
                print(f"=== {primer_name} ===")
                for rv in ref_versions:
                    tmp_prefix = tmp_dir / f"{primer_name}_{rv}"
                    row = verify_primer(primer_name, seq.strip().upper(), rv, tmp_prefix)
                    row["pair_name"] = r["pair_name"]
                    row["role"] = role
                    if args.population_check:
                        population_check(row, tmp_prefix)
                    rows.append(row)
                    ident = row.get("mrna_identity_pct", "")
                    pop = row.get("population_status", "")
                    print(f"  {rv}: {row['status']:28} gene={row.get('gene','-'):18} "
                          f"mRNA_identity={ident}%  population={pop}")

    fieldnames = ["pair_name", "primer", "role", "sequence", "ref_version", "status",
                  "chrom", "hit_start", "hit_end", "blast_pident", "blast_length",
                  "gene", "gene_description", "transcript", "exon_membership",
                  "mrna_match", "mrna_identity_pct", "mrna_gaps",
                  "population_status", "pseudogenome_region", "population_note"]
    out_csv = args.out_dir / "rtqpcr_primer_verification.csv"
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print(f"\nWritten: {out_csv}")


if __name__ == "__main__":
    main()
