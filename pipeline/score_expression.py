#!/usr/bin/env python3
"""Screen designs for expression liability and hard-filter failures.

Expression is an unscored gate: Adaptyv performs cell-free expression, and a design
that does not express is measured on none of the three judged axes. In Adaptyv's
EGFR round 1, 27% of designs never expressed; in the Bits to Binders CD20 campaign,
42% were not recovered.

Two independent competitions found that amino-acid composition predicts expression
far better than any structure-confidence metric:

                                    EGFR (expression)   CD20 (recovery)
    Glu content / K+E helix fraction      0.77              0.91
    Lys content                           0.73                -
    alpha-helix propensity                0.67                -
    ipTM                                  0.58           not predictive
    pLDDT                                 0.55           not predictive

Bits to Binders' retrospective filter lifted recovery from 56.8% to 81.2%, and
adding a no-cysteine rule raised the functional rate from 5.9% to 10.6%.

Nothing else in this pipeline looks at expression at all, and it is cheap, orthogonal
to every other signal here, and computable from sequence alone.

SCOPE NOTE. Optimisation-track designs differ from wild-type at only 10 of ~227
residues, so their global composition barely moves and these scores will cluster
tightly -- the useful output there is the hard-filter check, not the ranking. The
composition signal earns its keep on de novo designs, where the sequence is new.

Usage:
    score_expression.py --sequences <csv with a sequence column> [--column heavy_sequence]
    score_expression.py --cdrh3 <csv with a cdrh3 column>   # splices into trastuzumab VH
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# trastuzumab variable domains; CDRH3 occupies heavy-chain positions 99-108 (1-based)
TRASTUZUMAB_VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRY"
                  "ADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
TRASTUZUMAB_VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPS"
                  "RFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
CDRH3_START, CDRH3_END = 99, 108  # 1-based inclusive

# Chou-Fasman helix propensity
HELIX_PROPENSITY = {
    "A": 1.42, "C": 0.70, "D": 1.01, "E": 1.51, "F": 1.13, "G": 0.57, "H": 1.00,
    "I": 1.08, "K": 1.16, "L": 1.21, "M": 1.45, "N": 0.67, "P": 0.57, "Q": 1.11,
    "R": 0.98, "S": 0.77, "T": 0.83, "V": 1.06, "W": 1.08, "Y": 0.69,
}
CHARGE = {"K": 1.0, "R": 1.0, "H": 0.1, "D": -1.0, "E": -1.0}

# thresholds adapted from the Bits to Binders retrospective filter
MAX_REPEAT_RUN = 4        # longest run of one residue
MAX_GLU_RUN = 3           # consecutive glutamates


def longest_run(seq: str, residues: str | None = None) -> int:
    best = run = 0
    previous = ""
    for residue in seq:
        if residue == previous and (residues is None or residue in residues):
            run += 1
        else:
            run = 1 if (residues is None or residue in residues) else 0
        previous = residue
        best = max(best, run)
    return best


def assess(sequence: str) -> dict:
    n = len(sequence) or 1
    counts = {aa: sequence.count(aa) for aa in set(sequence)}
    glu = counts.get("E", 0)
    lys = counts.get("K", 0)
    cys = counts.get("C", 0)

    helix = sum(HELIX_PROPENSITY.get(aa, 1.0) for aa in sequence) / n
    charge = sum(CHARGE.get(aa, 0.0) for aa in sequence)

    # hard filters, which are disqualifying rather than merely unfavourable
    hard: list[str] = []
    # antibody variable domains carry one structural disulfide, so an even count is
    # expected; an odd count leaves an unpaired cysteine, which abguard rejects
    if cys % 2 == 1:
        hard.append(f"odd-cysteine-count({cys})")
    if len(sequence) > 250:
        hard.append(f"too-long({len(sequence)})")
    if re.search(r"N[^P][ST]", sequence):
        hard.append("n-glycosylation-sequon")

    warn: list[str] = []
    if longest_run(sequence) > MAX_REPEAT_RUN:
        warn.append(f"repeat-run({longest_run(sequence)})")
    if longest_run(sequence, "E") > MAX_GLU_RUN:
        warn.append(f"glu-run({longest_run(sequence, 'E')})")

    return {
        "length": len(sequence),
        "glu_frac": round(glu / n, 4),
        "lys_frac": round(lys / n, 4),
        "ke_frac": round((glu + lys) / n, 4),
        "cys_count": cys,
        "helix_propensity": round(helix, 3),
        "net_charge": round(charge, 1),
        "hard_fail": ";".join(hard),
        "warnings": ";".join(warn),
        "expressible": int(not hard),
    }


def splice_cdrh3(cdrh3: str) -> str:
    if len(cdrh3) != CDRH3_END - CDRH3_START + 1:
        raise ValueError(f"expected {CDRH3_END - CDRH3_START + 1} residues, got {len(cdrh3)}")
    return TRASTUZUMAB_VH[:CDRH3_START - 1] + cdrh3 + TRASTUZUMAB_VH[CDRH3_END:]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--sequences", type=Path, help="CSV of full sequences")
    source.add_argument("--cdrh3", type=Path, help="CSV of CDRH3 10-mers to splice in")
    parser.add_argument("--column", default=None, help="sequence column name")
    parser.add_argument("--paired", action="store_true",
                        help="append trastuzumab's light chain before scoring")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    path = args.sequences or args.cdrh3
    rows = list(csv.DictReader(path.open()))
    if not rows:
        print(f"error: no rows in {path}", file=sys.stderr)
        return 1

    column = args.column or ("cdrh3" if args.cdrh3 else "heavy_sequence")
    if column not in rows[0]:
        print(f"error: column '{column}' not in {path} (have: {', '.join(rows[0])})",
              file=sys.stderr)
        return 1

    # the wild-type construct, as the baseline every design is compared against
    wt_seq = TRASTUZUMAB_VH + (TRASTUZUMAB_VL if args.paired else "")
    wt = assess(wt_seq)
    print(f"wild-type trastuzumab {'VH+VL' if args.paired else 'VH'} ({wt['length']} aa)")
    print(f"  Glu {wt['glu_frac']:.3f}  Lys {wt['lys_frac']:.3f}  K+E {wt['ke_frac']:.3f}  "
          f"Cys {wt['cys_count']}  helix {wt['helix_propensity']}  charge {wt['net_charge']:+.1f}")
    print(f"  hard fails: {wt['hard_fail'] or 'none'}   warnings: {wt['warnings'] or 'none'}\n")

    failed = 0
    for row in rows:
        value = row[column].strip().upper()
        sequence = splice_cdrh3(value) if args.cdrh3 else value
        if args.paired and args.cdrh3:
            sequence += TRASTUZUMAB_VL
        row.update(assess(sequence))
        failed += not row["expressible"]

    print(f"scored {len(rows):,} designs")
    print(f"  hard-filter failures: {failed:,} ({100 * failed / len(rows):.1f}%)")
    reasons: dict[str, int] = {}
    for row in rows:
        for reason in filter(None, row["hard_fail"].split(";")):
            key = reason.split("(")[0]
            reasons[key] = reasons.get(key, 0) + 1
    for reason, count in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print(f"    {reason}: {count:,}")

    warned = sum(1 for row in rows if row["warnings"])
    print(f"  composition warnings: {warned:,}")

    if args.out:
        with args.out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote {len(rows):,} rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
