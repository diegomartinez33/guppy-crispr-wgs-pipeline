#!/usr/bin/env python3
"""
Audit ko_guide_scan.py's already-reported guide candidates for the "exon
junction artifact" risk raised by the user: find_ngg_candidates() scans the
CONCATENATED CDS (exons spliced together, introns removed), so a 23bp
candidate window whose position straddles two exons in that concatenated
string does not correspond to any real, contiguous stretch of genomic DNA -
Cas9 cannot actually target it (it targets genomic DNA, not spliced mRNA).
Any such candidate in the existing *_guide_comparison.csv outputs would be
a biologically meaningless entry, however classified.

For each gene/population already analyzed, recomputes the exact same
transcript selection and CDS exon layout ko_guide_scan.py used (reusing its
own functions, not reimplementing the logic), determines exon-junction
boundary positions within the concatenated reference CDS, and checks every
already-reported candidate row against them.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from ko_guide_scan import GENOME_CHOICES, OUT_DIR, find_gene_features  # noqa: E402

GENES = ["bdnf", "agap3", "grin1a", "grin1b", "gria1a", "gria1b", "gria2b"]
POPULATION = "pseudogenome"


def junction_boundaries(cds_lines, strand):
    """Return the set of 1-based positions in the concatenated CDS that are
    the LAST base of an exon (i.e., position b -> a junction sits between
    b and b+1), in transcript (5'->3') order."""
    lengths = [int(f[4]) - int(f[3]) + 1 for f in cds_lines]  # cds_lines sorted ascending by genomic start
    if strand == "-":
        lengths = list(reversed(lengths))  # extract_cds() reverses segment order for "-" strand
    boundaries = set()
    cum = 0
    for length in lengths[:-1]:
        cum += length
        boundaries.add(cum)
    return boundaries, lengths


def crosses_junction(start0, boundaries):
    """start0: 0-based start of the 23bp window. Window covers 1-based
    positions [start0+1, start0+23]. Crosses a junction at boundary b (last
    base of an exon) if the window includes both b and b+1."""
    start1 = start0 + 1
    end1 = start0 + 23
    return any(start1 <= b <= end1 - 1 for b in boundaries)


def main():
    pop = GENOME_CHOICES[POPULATION]
    total_flagged = 0
    for gene in GENES:
        ref_gene, ref_mrnas = find_gene_features(pop["ref_gff"], gene)
        pop_gene, pop_mrnas = find_gene_features(pop["gff"], gene)
        if ref_gene is None or pop_gene is None:
            print(f"{gene}: SKIP (gene not found)")
            continue
        ref_strand = ref_gene[6]
        shared_tids = set(ref_mrnas) & set(pop_mrnas)
        if not shared_tids:
            print(f"{gene}: SKIP (no shared transcript)")
            continue
        chosen_tid = max(shared_tids, key=lambda t: sum(int(f[4]) - int(f[3]) + 1 for f in ref_mrnas[t][1]))
        _, ref_cds = ref_mrnas[chosen_tid]
        boundaries, lengths = junction_boundaries(ref_cds, ref_strand)

        csv_path = OUT_DIR / f"{gene}_{POPULATION}_guide_comparison.csv"
        if not csv_path.exists():
            print(f"{gene}: SKIP (no CSV at {csv_path})")
            continue

        flagged = []
        with open(csv_path) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                start0 = int(row["cds_pos_1based"]) - 1
                if crosses_junction(start0, boundaries):
                    flagged.append(row)

        print(f"{gene}: transcript={chosen_tid} exons={len(lengths)} exon_lengths={lengths} "
              f"junction_boundaries={sorted(boundaries)}")
        print(f"  Candidates flagged as exon-junction artifacts: {len(flagged)}")
        for row in flagged:
            print(f"    cds_pos={row['cds_pos_1based']} strand={row['strand']} "
                  f"classification={row['classification']} spacer={row['ref_spacer']} pam={row['ref_pam']}")
        total_flagged += len(flagged)

    print(f"\n=== TOTAL exon-junction-artifact candidates across all genes: {total_flagged} ===")


if __name__ == "__main__":
    main()
