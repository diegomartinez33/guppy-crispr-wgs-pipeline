#!/usr/bin/env python3
"""
ko_guide_scan.py — Compare CRISPR knockout guide candidates (SpCas9, NGG PAM)
between the NCBI Trinidad reference and the Colombian population genome
(pseudogenome by default) for a given gene.

Why: guides designed purely against the reference can be invalidated (PAM
destroyed) or weakened (seed-region variant) by real population-level
variants — exactly the kind of population-specific validation this
project's whole assembly/pseudogenome effort was built to support.

What it does, per gene:
  1. Look up the gene in both GFF3 annotations (by `gene=NAME;` — Liftoff
     preserves the reference's gene naming, confirmed 2026-08-31 for
     bdnf/agap3/grin1a/grin1b/gria1a/gria1b/gria2b).
  2. Extract the full gene span (both genomes) and align it with minimap2
     --cs to report all variants (SNPs/indels) across the gene body.
  3. Extract the CDS of the first-listed mRNA isoform (both genomes,
     strand-aware, exons concatenated in transcript order) - the target
     for knockout guide design.
  4. Align CDS-to-CDS (minimap2 --cs) for precise coordinate mapping
     between reference and population CDS positions.
  5. Enumerate all NGG-PAM candidate guides (20bp spacer + NGG, both
     strands) in the reference CDS. For each, classify using the CDS
     alignment:
       IDENTICAL      - same 23bp window in the population, guide unaffected
       PAM_BROKEN      - population variant destroys the NGG PAM
       SEED_VARIANT     - variant in the 10bp proximal to the PAM (most
                          critical for Cas9 binding/cutting)
       DISTAL_VARIANT   - variant in the distal 10bp of the spacer
       NO_ALIGNMENT     - reference window has no aligned counterpart
                          (falls in a gap/indel/unaligned region)
       EXON_JUNCTION_ARTIFACT - window straddles an exon-exon junction in
                          the spliced CDS used for PAM scanning - not a
                          real contiguous genomic sequence, so not an
                          actually editable Cas9 target regardless of what
                          the window's classification would otherwise be
                          (found 2026-09-06; corrected one gria2b false
                          positive in the 2026-08-31 results, see below).
  6. Separately report population-only candidates (NGG sites created by a
     population variant, absent from the reference) - not visible if you
     only ever design against the reference.
  7. (Optional, on by default when available) Run CRISPOR (via Singularity,
     2026-09 addition) against both the reference and population CDS, using
     custom genomes registered in the container (guppyRefTrinidad /
     guppyColPseudogenome). This ADDS real BWA-based off-target counts
     (mitSpecScore, offtargetCount) and efficiency scores (Doench '16 /
     Rule Set 2) to each guide already found by the manual PAM scan above -
     it does not replace step 5's variant classification, which is the only
     thing that tells you a guide is invalidated/weakened by a population
     variant. Guides are matched between the manual scan and CRISPOR by
     exact 23bp spacer+PAM sequence (both use the same NGG/20bp definition,
     reported 5'->3' on the guide's own strand).

No pysam/biopython dependency (avoids touching crispresso2_env, which
CRISPResso2 depends on) - uses samtools faidx (subprocess) for extraction
and minimap2 --cs (subprocess) for alignment, both already used elsewhere
in this project. CRISPOR runs inside a Singularity container (Docker itself
needs the docker group, not available on this cluster) - see
codes/analysis/crispor_add_genomes.sh for how the custom genomes were
registered.

Usage:
    python3 ko_guide_scan.py --gene bdnf --population pseudogenome
    python3 ko_guide_scan.py --gene grin1a --population scaffolded
    python3 ko_guide_scan.py --gene bdnf --no-crispor   # manual scan only
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

PROJECT_DIR = Path("/hpcfs/home/ing_civil/da.martinez33/UBC/off-target_data")

# v1 = GCF_000633615.1 (Trinidad/Guanapo, female, short-read, 2014) -
# suppressed by NCBI as of 2026, kept for continuity/comparison.
# v2 = GCF_904066995.2 (male, PacBio+Hi-C, 2025) - current NCBI RefSeq
# reference genome for the species. Each population entry below is paired
# with the reference it was actually built against (a v2 pseudogenome must
# be compared to the v2 reference, not v1) - see CLAUDE.md, "Migration to
# GCF_904066995.2 (v2)".
REF_FASTA_V1 = PROJECT_DIR / "reference/GCF_000633615.1_Guppy_female_1.0_MT_genomic.fna"
REF_GFF_V1 = PROJECT_DIR / "reference/GCF_000633615.1_annotation.gff"
REF_FASTA_V2 = PROJECT_DIR / "reference/GCF_904066995.2_P_reticulata-male-v2_genomic.fna"
REF_GFF_V2 = PROJECT_DIR / "reference/GCF_904066995.2_annotation.gff"

GENOME_CHOICES = {
    "pseudogenome": {
        "ref_fasta": REF_FASTA_V1,
        "ref_gff": REF_GFF_V1,
        "fasta": PROJECT_DIR / "reference/pseudogenome/colombian_pseudogenome.fna",
        "gff": PROJECT_DIR / "reference/pseudogenome/colombian_pseudogenome.liftoff.gff3",
        "label": "Colombian pseudogenome (bcftools consensus + Liftoff) vs v1 reference",
    },
    "scaffolded": {
        "ref_fasta": REF_FASTA_V1,
        "ref_gff": REF_GFF_V1,
        "fasta": PROJECT_DIR / "reference/colombian_scaffolded_genome/colombian_scaffolded.fna",
        "gff": PROJECT_DIR / "reference/colombian_scaffolded_genome/colombian_scaffolded.liftoff.gff3",
        "label": "Colombian de novo scaffolded genome (SPAdes+RagTag+Liftoff) vs v1 reference",
    },
    "pseudogenome_v2": {
        "ref_fasta": REF_FASTA_V2,
        "ref_gff": REF_GFF_V2,
        "fasta": PROJECT_DIR / "reference/pseudogenome_v2/colombian_pseudogenome.fna",
        "gff": PROJECT_DIR / "reference/pseudogenome_v2/colombian_pseudogenome.liftoff.gff3",
        "label": "Colombian pseudogenome (bcftools consensus + Liftoff) vs v2 reference (GCF_904066995.2)",
    },
}

OUT_DIR = PROJECT_DIR / "analysis" / "ko_guide_scan"

# CRISPOR (Singularity) - custom genomes registered via crispor_add_genomes.sh.
# guppyRefTrinidad/guppyRefMaleV2 are the two registered reference genomes;
# population-side scoring is skipped gracefully (with a printed note) for
# any population without a registered CRISPOR genome id (None below).
CRISPOR_SIF = PROJECT_DIR / "codes/analysis/crispor_singularity/crispor_v5.2c_amd64.sif"
CRISPOR_GENOMES_DIR = PROJECT_DIR / "codes/analysis/crispor_singularity/genomes"
CRISPOR_REF_GENOME_IDS = {
    REF_FASTA_V1: "guppyRefTrinidad",
    REF_FASTA_V2: "guppyRefMaleV2",
}
CRISPOR_GENOME_IDS = {
    "pseudogenome": "guppyColPseudogenome",
    "scaffolded": None,  # not registered as a CRISPOR genome (yet)
    "pseudogenome_v2": "guppyColPseudogenomeV2",
}

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"ERROR running {' '.join(cmd)}:\n{r.stderr}")
    return r.stdout


# ── GFF parsing ──────────────────────────────────────────────────────────


def transcript_id_of(mrna_fields):
    """Liftoff preserves the original RefSeq mRNA ID (e.g. 'rna-XM_008405157.2')
    across genomes, but NOT its listing order in the GFF - strip the
    'rna-' prefix so the same transcript can be matched by ID between the
    reference and population GFFs regardless of order."""
    m = re.search(r"ID=([^;]+)", mrna_fields[8])
    tid = m.group(1)
    return tid[4:] if tid.startswith("rna-") else tid


def find_gene_features(gff_path, gene_name):
    """Return (gene_line, {transcript_id: (mrna_line, [cds_lines_sorted])})."""
    gene_line = None
    mrna_lines = []
    with open(gff_path) as fh:
        lines = fh.readlines()

    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        feat_type, attrs = fields[2], fields[8]
        if feat_type == "gene" and f"gene={gene_name};" in attrs or (
            feat_type == "gene" and attrs.endswith(f"gene={gene_name}")
        ):
            gene_line = fields
        elif feat_type == "mRNA" and f"gene={gene_name};" in attrs:
            mrna_lines.append(fields)

    if gene_line is None:
        return None, {}

    mrna_by_tid = {}
    for mrna in mrna_lines:
        mrna_id_match = re.search(r"ID=([^;]+)", mrna[8])
        if not mrna_id_match:
            continue
        mrna_id = mrna_id_match.group(1)
        tid = transcript_id_of(mrna)
        cds_lines = []
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[2] == "CDS" and f"Parent={mrna_id}" in fields[8]:
                cds_lines.append(fields)
        if cds_lines:
            cds_lines.sort(key=lambda f: int(f[3]))
            mrna_by_tid[tid] = (mrna, cds_lines)

    return gene_line, mrna_by_tid


# ── Sequence extraction ──────────────────────────────────────────────────


def faidx_seq(fasta_path, chrom, start, end):
    """1-based inclusive coordinates, matching GFF convention."""
    region = f"{chrom}:{start}-{end}"
    out = run(["samtools", "faidx", str(fasta_path), region])
    seq_lines = out.splitlines()[1:]
    return "".join(seq_lines).upper()


def extract_cds(fasta_path, chrom, cds_lines, strand):
    """Concatenate CDS segments in transcript (5'->3') order."""
    segments = []
    for f in cds_lines:
        start, end = int(f[3]), int(f[4])
        segments.append(faidx_seq(fasta_path, chrom, start, end))
    if strand == "-":
        segments = [revcomp(s) for s in reversed(segments)]
    return "".join(segments)


# ── Alignment (minimap2 --cs) ────────────────────────────────────────────


def align_cs(ref_seq, qry_seq, tmp_prefix):
    ref_fa = f"{tmp_prefix}.ref.fa"
    qry_fa = f"{tmp_prefix}.qry.fa"
    with open(ref_fa, "w") as fh:
        fh.write(f">ref\n{ref_seq}\n")
    with open(qry_fa, "w") as fh:
        fh.write(f">qry\n{qry_seq}\n")
    out = run(["minimap2", "-a", "--cs", "-x", "asm20", ref_fa, qry_fa])
    cs_tag = None
    cigar = None
    pos = None
    flag = None
    for line in out.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        flag = int(f[1])
        if flag & 4:  # unmapped
            continue
        pos = int(f[3])  # 1-based leftmost REF position
        cigar = f[5]
        for tag in f[11:]:
            if tag.startswith("cs:Z:"):
                cs_tag = tag[5:]
        break
    return pos, cigar, cs_tag, flag


def parse_cs_variants(cs_tag, ref_start):
    """Walk a short-form `cs` tag, return list of (ref_pos_1based, kind,
    ref_allele, alt_allele) and a ref_pos -> qry_pos mapping array (for
    positions covered by the alignment)."""
    variants = []
    ref_pos = ref_start
    qry_pos = 1
    # position map: ref_pos (1-based) -> qry_pos (1-based), only for
    # positions actually present in both (identical or substituted)
    pos_map = {}
    tokens = re.findall(r":[0-9]+|\*[a-z]{2}|\+[a-z]+|-[a-z]+", cs_tag)
    for tok in tokens:
        op = tok[0]
        if op == ":":
            n = int(tok[1:])
            for i in range(n):
                pos_map[ref_pos + i] = qry_pos + i
            ref_pos += n
            qry_pos += n
        elif op == "*":
            ref_b, alt_b = tok[1], tok[2]
            variants.append((ref_pos, "SNP", ref_b.upper(), alt_b.upper()))
            pos_map[ref_pos] = qry_pos
            ref_pos += 1
            qry_pos += 1
        elif op == "+":
            ins = tok[1:]
            variants.append((ref_pos, "INS", "-", ins.upper()))
            qry_pos += len(ins)
        elif op == "-":
            dele = tok[1:]
            variants.append((ref_pos, "DEL", dele.upper(), "-"))
            ref_pos += len(dele)
    return variants, pos_map


# ── PAM enumeration ──────────────────────────────────────────────────────


def find_ngg_candidates(seq):
    """Return list of (spacer_start_0based, strand, spacer20, pam3)."""
    candidates = []
    n = len(seq)
    for i in range(n - 2):
        if seq[i + 1 : i + 3] == "GG" and i >= 20:
            spacer = seq[i - 20 : i]
            pam = seq[i : i + 3]
            candidates.append((i - 20, "+", spacer, pam))
    rc = revcomp(seq)
    n_rc = len(rc)
    for i in range(n_rc - 2):
        if rc[i + 1 : i + 3] == "GG" and i >= 20:
            spacer = rc[i - 20 : i]
            pam = rc[i : i + 3]
            # convert rc-coordinate back to + strand spacer_start (0-based)
            orig_start = n - (i + 3)
            candidates.append((orig_start, "-", spacer, pam))
    return candidates


def exon_junction_boundaries(cds_lines, strand):
    """find_ngg_candidates() scans the CONCATENATED CDS (exons spliced
    together, introns removed) - a 23bp candidate window whose position
    straddles two exons in that concatenated string does not correspond to
    any real, contiguous stretch of genomic DNA, so Cas9 cannot actually
    target it (it cuts genomic DNA, not spliced mRNA). Found via user
    question 2026-09-06 - audited against all 7 previously-reported genes:
    353/2699 candidates crossed a junction, all but one already classified
    IDENTICAL (harmless); one gria2b DISTAL_VARIANT call was a genuine
    false positive, corrected by this check.

    Returns the set of 1-based positions in the concatenated CDS that are
    the LAST base of an exon (transcript 5'->3' order) - a junction sits
    between position b and b+1."""
    lengths = [int(f[4]) - int(f[3]) + 1 for f in cds_lines]  # cds_lines sorted ascending by genomic start
    if strand == "-":
        lengths = list(reversed(lengths))  # extract_cds() reverses segment order for "-" strand
    boundaries = set()
    cum = 0
    for length in lengths[:-1]:
        cum += length
        boundaries.add(cum)
    return boundaries


def crosses_exon_junction(start0, boundaries):
    """start0: 0-based start of the 23bp window. Window covers 1-based
    positions [start0+1, start0+23]. Crosses a junction at boundary b (last
    base of an exon) if the window includes both b and b+1."""
    start1 = start0 + 1
    end1 = start0 + 23
    return any(start1 <= b <= end1 - 1 for b in boundaries)


def classify_candidate(ref_start0, ref_spacer, pos_map, variants_by_pos, qry_seq, junction_boundaries_set=frozenset()):
    """ref_start0: 0-based start of the 23bp window (spacer+PAM) in the
    reference CDS. Returns (classification, qry_window_or_None, note)."""
    if crosses_exon_junction(ref_start0, junction_boundaries_set):
        return "EXON_JUNCTION_ARTIFACT", None, "window spans an exon-exon junction in the spliced CDS - not a real contiguous genomic target"

    window_ref_positions = list(range(ref_start0 + 1, ref_start0 + 24))  # 1-based
    pam_positions = window_ref_positions[20:23]
    seed_positions = window_ref_positions[10:20]
    distal_positions = window_ref_positions[0:10]

    touched = [p for p in window_ref_positions if p in variants_by_pos]
    if not all(p in pos_map for p in window_ref_positions):
        return "NO_ALIGNMENT", None, "window not fully aligned (gap/indel nearby)"

    if any(p in pam_positions for p in touched):
        return "PAM_BROKEN", None, f"variant at ref pos {[p for p in touched if p in pam_positions]}"
    if any(p in seed_positions for p in touched):
        qstart = pos_map[window_ref_positions[0]]
        qend = pos_map[window_ref_positions[-1]]
        qwin = qry_seq[qstart - 1 : qend]
        return "SEED_VARIANT", qwin, f"variant at ref pos {[p for p in touched if p in seed_positions]}"
    if any(p in distal_positions for p in touched):
        qstart = pos_map[window_ref_positions[0]]
        qend = pos_map[window_ref_positions[-1]]
        qwin = qry_seq[qstart - 1 : qend]
        return "DISTAL_VARIANT", qwin, f"variant at ref pos {[p for p in touched if p in distal_positions]}"

    qstart = pos_map[window_ref_positions[0]]
    qend = pos_map[window_ref_positions[-1]]
    qwin = qry_seq[qstart - 1 : qend]
    return "IDENTICAL", qwin, ""


# ── CRISPOR (Singularity, optional) ─────────────────────────────────────


def run_crispor(seq, genome_id, tmp_fa_prefix, out_tsv_prefix):
    """Run crispor.py (SpCas9 NGG defaults) inside the Singularity container
    against a registered custom genome, for real BWA-based off-target
    counts/scores and efficiency scores - this only ADDS scores to guides
    already found by the manual PAM scan, it never decides classification.
    Returns {23bp targetSeq (guide's own strand): {column_name: value}}."""
    in_fa = f"{tmp_fa_prefix}.fa"
    guide_out = f"{out_tsv_prefix}_crispor_guides.tsv"
    offt_out = f"{out_tsv_prefix}_crispor_offs.tsv"
    with open(in_fa, "w") as fh:
        fh.write(f">cds\n{seq}\n")
    cmd = [
        "singularity", "exec",
        "-B", f"{CRISPOR_GENOMES_DIR}:/data/genomes",
        str(CRISPOR_SIF),
        "/data/www/crispor/crispor.py",
        "-g", "/data/genomes",
        genome_id, in_fa, guide_out,
        "-o", offt_out,
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"WARNING: CRISPOR failed for genome '{genome_id}': {r.stderr.strip()[-800:]}", file=sys.stderr)
        return {}
    if not Path(guide_out).exists():
        print(f"WARNING: CRISPOR produced no guide output for genome '{genome_id}'", file=sys.stderr)
        return {}
    scores = {}
    with open(guide_out) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                continue
            row = dict(zip(header, fields))
            scores[row.get("targetSeq", "")] = row
    return scores


# ── Main ──────────────────────────────────────────────────────────────────


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene", required=True, help="Gene symbol, e.g. bdnf, grin1a")
    ap.add_argument(
        "--population",
        choices=list(GENOME_CHOICES.keys()),
        default="pseudogenome",
        help="Which Colombian genome resource to compare against (default: pseudogenome)",
    )
    ap.add_argument(
        "--no-crispor",
        action="store_true",
        help="Skip CRISPOR (Singularity) scoring even if the container/genomes are "
        "available. The manual PAM scan + variant classification always runs "
        "regardless of this flag.",
    )
    args = ap.parse_args()

    gene = args.gene
    pop = GENOME_CHOICES[args.population]
    REF_FASTA = pop["ref_fasta"]
    REF_GFF = pop["ref_gff"]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tmp_prefix = f"/tmp/ko_guide_scan_{gene}"

    print(f"=== {gene} vs {pop['label']} ===")

    ref_gene, ref_mrnas = find_gene_features(REF_GFF, gene)
    pop_gene, pop_mrnas = find_gene_features(pop["gff"], gene)

    if ref_gene is None:
        sys.exit(f"ERROR: gene '{gene}' not found in reference GFF ({REF_GFF})")
    if pop_gene is None:
        sys.exit(f"ERROR: gene '{gene}' not found in {args.population} GFF ({pop['gff']})")
    if not ref_mrnas:
        sys.exit(f"ERROR: gene '{gene}' found but has no mRNA/CDS in reference GFF")
    if not pop_mrnas:
        sys.exit(f"ERROR: gene '{gene}' found but has no mRNA/CDS in {args.population} GFF")

    ref_chrom, ref_strand = ref_gene[0], ref_gene[6]
    ref_start, ref_end = int(ref_gene[3]), int(ref_gene[4])
    pop_chrom, pop_strand = pop_gene[0], pop_gene[6]
    pop_start, pop_end = int(pop_gene[3]), int(pop_gene[4])

    print(f"Reference: {ref_chrom}:{ref_start}-{ref_end} ({ref_strand})")
    print(f"{args.population}: {pop_chrom}:{pop_start}-{pop_end} ({pop_strand})")

    # --- full gene-body variant scan ---
    ref_gene_seq = faidx_seq(REF_FASTA, ref_chrom, ref_start, ref_end)
    pop_gene_seq = faidx_seq(pop["fasta"], pop_chrom, pop_start, pop_end)
    pos, cigar, cs_tag, flag = align_cs(ref_gene_seq, pop_gene_seq, tmp_prefix + "_genebody")
    gene_variants = []
    if cs_tag:
        gene_variants, _ = parse_cs_variants(cs_tag, pos)
    print(f"Gene-body variants (full span, incl. introns): {len(gene_variants)}")

    # --- CDS extraction: pick a transcript ID present in BOTH genomes,
    # preferring the one with the longest reference CDS (proxy for the
    # principal/canonical isoform - these GFFs don't mark one explicitly).
    # Liftoff preserves transcript IDs but NOT their listing order, so
    # matching by order alone (e.g. "first mRNA in each file") silently
    # compares two DIFFERENT isoforms - confirmed 2026-08-31 for bdnf
    # (ref's first-listed transcript != pseudogenome's first-listed one).
    shared_tids = set(ref_mrnas) & set(pop_mrnas)
    if not shared_tids:
        sys.exit(
            f"ERROR: no shared transcript ID between reference and {args.population} "
            f"for gene '{gene}' (ref has {len(ref_mrnas)}, pop has {len(pop_mrnas)})"
        )
    chosen_tid = max(shared_tids, key=lambda t: sum(int(f[4]) - int(f[3]) + 1 for f in ref_mrnas[t][1]))
    print(f"Representative transcript: {chosen_tid} ({len(shared_tids)} shared isoforms available)")
    ref_mrna, ref_cds = ref_mrnas[chosen_tid]
    pop_mrna, pop_cds = pop_mrnas[chosen_tid]
    ref_cds_seq = extract_cds(REF_FASTA, ref_chrom, ref_cds, ref_strand)
    pop_cds_seq = extract_cds(pop["fasta"], pop_chrom, pop_cds, pop_strand)
    print(f"Reference CDS length: {len(ref_cds_seq)}bp ({len(ref_cds)} exons)")
    print(f"{args.population} CDS length: {len(pop_cds_seq)}bp ({len(pop_cds)} exons)")

    # Candidates whose 23bp window straddles an exon-exon junction in the
    # concatenated CDS aren't real contiguous genomic DNA - see
    # exon_junction_boundaries() docstring.
    ref_junction_boundaries = exon_junction_boundaries(ref_cds, ref_strand)

    # --- CDS-to-CDS alignment for precise guide-window mapping ---
    pos_c, cigar_c, cs_tag_c, flag_c = align_cs(ref_cds_seq, pop_cds_seq, tmp_prefix + "_cds")
    if not cs_tag_c:
        sys.exit(f"ERROR: CDS sequences for {gene} did not align (unmapped) - cannot compare guides")
    cds_variants, pos_map = parse_cs_variants(cs_tag_c, pos_c)
    variants_by_pos = {v[0]: v for v in cds_variants}
    print(f"CDS variants: {len(cds_variants)}")

    # --- CRISPOR (Singularity) scoring - optional, additive only ---
    crispor_ref_scores, crispor_pop_scores = {}, {}
    if args.no_crispor:
        print("\nCRISPOR scoring skipped (--no-crispor)")
    elif not CRISPOR_SIF.exists():
        print(f"\nNOTE: CRISPOR image not found at {CRISPOR_SIF} - skipping CRISPOR scoring "
              "(manual PAM scan results above are unaffected)")
    else:
        print("\n=== Running CRISPOR (Singularity) ===")
        ref_genome_id = CRISPOR_REF_GENOME_IDS.get(REF_FASTA)
        if ref_genome_id:
            crispor_ref_scores = run_crispor(
                ref_cds_seq, ref_genome_id, tmp_prefix + "_crispor_ref", str(OUT_DIR / f"{gene}_reference")
            )
            print(f"CRISPOR reference guides scored: {len(crispor_ref_scores)}")
        else:
            print(f"NOTE: no CRISPOR genome registered for reference {REF_FASTA} - "
                  "skipping reference-side CRISPOR scoring")
        pop_genome_id = CRISPOR_GENOME_IDS.get(args.population)
        if pop_genome_id:
            crispor_pop_scores = run_crispor(
                pop_cds_seq, pop_genome_id, tmp_prefix + "_crispor_pop", str(OUT_DIR / f"{gene}_{args.population}")
            )
            print(f"CRISPOR {args.population} guides scored: {len(crispor_pop_scores)}")
        else:
            print(f"NOTE: no CRISPOR genome registered for '{args.population}' yet - "
                  "skipping population-side CRISPOR scoring")

    # --- enumerate + classify guide candidates ---
    ref_candidates = find_ngg_candidates(ref_cds_seq)
    rows = []
    for start0, strand, spacer, pam in ref_candidates:
        cls, qwin, note = classify_candidate(
            start0, spacer, pos_map, variants_by_pos, pop_cds_seq, ref_junction_boundaries
        )
        target = spacer + pam
        ref_c = crispor_ref_scores.get(target, {})
        pop_c = crispor_pop_scores.get(qwin, {}) if qwin else {}
        rows.append(
            {
                "gene": gene,
                "cds_pos_1based": start0 + 1,
                "strand": strand,
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
                "pop_crispor_doench16_score": pop_c.get("Doench '16-Score", ""),
            }
        )

    pop_junction_boundaries = exon_junction_boundaries(pop_cds, pop_strand)
    pop_only_candidates = find_ngg_candidates(pop_cds_seq)
    ref_windows = {(s0, st) for s0, st, _, _ in ref_candidates}
    # crude novel-site detection: pop candidate whose mapped ref position
    # has no equivalent ref candidate at all (only meaningful near variants)
    novel_rows = []
    inv_pos_map = {}
    for rp, qp in pos_map.items():
        inv_pos_map.setdefault(qp, rp)
    for start0, strand, spacer, pam in pop_only_candidates:
        if crosses_exon_junction(start0, pop_junction_boundaries):
            continue  # not a real contiguous genomic target, see exon_junction_boundaries()
        qpos = start0 + 1
        if qpos not in inv_pos_map:
            novel_rows.append(
                {
                    "gene": gene,
                    "cds_pos_1based_pop": qpos,
                    "strand": strand,
                    "pop_spacer": spacer,
                    "pop_pam": pam,
                }
            )

    # --- write outputs ---
    import csv

    out_csv = OUT_DIR / f"{gene}_{args.population}_guide_comparison.csv"
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else [])
        w.writeheader()
        w.writerows(rows)

    out_variants_csv = OUT_DIR / f"{gene}_{args.population}_gene_body_variants.csv"
    with open(out_variants_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ref_pos_1based_in_gene_span", "type", "ref_allele", "alt_allele"])
        w.writerows(gene_variants)

    # --- summary ---
    from collections import Counter

    counts = Counter(r["classification"] for r in rows)
    print("\n=== Guide candidate summary (reference CDS) ===")
    print(f"Total NGG candidates in reference CDS: {len(rows)}")
    for cls in ["IDENTICAL", "PAM_BROKEN", "SEED_VARIANT", "DISTAL_VARIANT", "NO_ALIGNMENT", "EXON_JUNCTION_ARTIFACT"]:
        print(f"  {cls}: {counts.get(cls, 0)}")
    print(f"Population-only novel PAM sites (not present in reference): {len(novel_rows)}")
    if crispor_ref_scores:
        matched = sum(1 for r in rows if r["ref_crispor_mitSpecScore"] != "")
        print(f"CRISPOR reference scores matched by exact sequence: {matched}/{len(rows)} guides")
    if crispor_pop_scores:
        matched_pop = sum(1 for r in rows if r["pop_crispor_mitSpecScore"] != "")
        print(f"CRISPOR {args.population} scores matched by exact sequence: {matched_pop}/{len(rows)} guides")
    print(f"\nWritten: {out_csv}")
    print(f"Written: {out_variants_csv}")


if __name__ == "__main__":
    main()
