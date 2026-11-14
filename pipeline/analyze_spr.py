#!/usr/bin/env python3
"""What the organisers' SPR data says about beating trastuzumab.

The competition scores affinity against measured trastuzumab. This script asks the
provided SPR table two questions that decide whether that is realistic:

  1. How noisy is the assay? Wild-type appears several times, so its own replicate
     spread is the floor below which an "improvement" is not measurable.
  2. How does affinity track the number of CDRH3 substitutions? The novelty gate
     forces >= 3 substitutions, so the measured penalty at >= 3 is the cost of
     being allowed to submit at all.

Data: trastuzumab_her2_spr_all_768.csv from s3://apheris-oslo-hackathon/data/
her2-trastuzumab/, from Chinery et al. (doi:10.1101/2024.03.26.586756). Cite them.

Usage:
    analyze_spr.py [--spr hackathon-data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv]
                   [--designs WGPPYFYAMD WASSSEYYFD]
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path

# Trastuzumab's wild-type CDRH3, heavy chain B99-108.
WILD_TYPE = "WGGDGFYAMD"

# abguard's substitution floors, by the track it assigns (see HANDOFF.md).
NOVELTY_FLOOR_OPTIMISATION = 3
NOVELTY_FLOOR_DE_NOVO = 2

DEFAULT_SPR = Path("hackathon-data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv")


def substitutions(seq: str, reference: str = WILD_TYPE) -> int:
    """Edit distance to the reference.

    Most CDRH3s here are 10 residues like wild-type, where this is just the
    Hamming distance — but the SPR table also holds 9- and 16-residue loops
    (the README's "always 10 residues" describes the FACS table, not this one),
    so indels have to count. Levenshtein, inlined to keep this dependency-free.
    """
    if len(seq) == len(reference):
        return sum(a != b for a, b in zip(seq, reference))

    prev = list(range(len(reference) + 1))
    for i, a in enumerate(seq, 1):
        cur = [i]
        for j, b in enumerate(reference, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a != b)))
        prev = cur
    return prev[-1]


def load_measured(path: Path) -> list[tuple[float, str]]:
    """Return [(kd_nM, cdrh3)] for rows with a quantified KD, tightest first.

    Only 190 of the 768 rows have one; the rest are no_reportable_kd,
    inconclusive or not_reported, and silently dropping them is the point —
    a missing KD is not a weak binder, it is an unmeasured one.
    """
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    measured = [
        (float(r["kd_nM"]), r["cdrh3"])
        for r in rows
        if r.get("kd_nM") not in (None, "", "NA", "nan")
    ]
    measured.sort()
    print(f"{len(rows)} rows, {len(measured)} with a quantified KD", file=sys.stderr)
    return measured


def report_assay_floor(measured: list[tuple[float, str]]) -> float:
    """Wild-type's replicate spread. Returns the median wild-type KD."""
    wt = [kd for kd, seq in measured if seq == WILD_TYPE]
    if not wt:
        sys.exit("wild-type not found in the SPR table — cannot calibrate")

    lo, hi, med = min(wt), max(wt), statistics.median(wt)
    print("## Assay floor — wild-type measured against itself\n")
    print(f"{len(wt)} wild-type replicates: {lo:.4f} - {hi:.4f} nM, median {med:.3f} nM")
    print(f"**{hi / lo:.1f}x spread on an identical sequence.**")
    print("An affinity gain smaller than that is not distinguishable from noise.\n")
    return med


def report_substitution_cost(measured: list[tuple[float, str]], wt_median: float) -> None:
    print("## Affinity vs. number of CDRH3 substitutions\n")
    print("| subs | n | median KD (nM) | best KD (nM) | fold worse than WT |")
    print("| ---: | ---: | ---: | ---: | ---: |")
    for d in range(len(WILD_TYPE) + 1):
        group = [kd for kd, seq in measured if substitutions(seq) == d]
        if not group:
            continue
        med = statistics.median(group)
        flag = " ←" if d == 0 else ""
        print(
            f"| {d} | {len(group)} | {med:.3f} | {min(group):.4f} "
            f"| {med / wt_median:.0f}x{flag} |"
        )
    print()


def report_novelty_tension(measured: list[tuple[float, str]], wt_median: float) -> None:
    """The central conflict: novelty demands substitutions, affinity punishes them."""
    eligible = [(kd, seq) for kd, seq in measured if substitutions(seq) >= NOVELTY_FLOOR_OPTIMISATION]
    better = [(kd, seq) for kd, seq in eligible if kd < wt_median]

    print("## The novelty/affinity tension, measured\n")
    print(
        f"Of {len(measured)} variants with a quantified KD, {len(eligible)} clear the "
        f">= {NOVELTY_FLOOR_OPTIMISATION}-substitution novelty floor."
    )
    print(
        f"**Of those, {len(better)} beat wild-type's median KD "
        f"({len(better) / len(eligible) * 100:.1f}%).**\n"
    )
    if better:
        print("| KD (nM) | CDRH3 | subs |")
        print("| ---: | --- | ---: |")
        for kd, seq in sorted(better):
            print(f"| {kd:.4f} | `{seq}` | {substitutions(seq)} |")
        print()


def report_designs(measured: list[tuple[float, str]], designs: list[str], wt_median: float) -> None:
    """Place our designs in the measured distribution by substitution count."""
    if not designs:
        return
    print("## Our designs, placed in the measured distribution\n")
    print("No design here has been measured — this is the prior for its substitution count.\n")
    print("| design CDRH3 | subs | median KD at that count | best ever seen | in table? |")
    print("| --- | ---: | ---: | ---: | --- |")
    for seq in designs:
        d = substitutions(seq)
        group = [kd for kd, s in measured if substitutions(s) == d]
        seen = "**yes — not novel**" if any(s == seq for _, s in measured) else "no"
        med = f"{statistics.median(group):.1f}" if group else "-"
        best = f"{min(group):.4f}" if group else "-"
        print(f"| `{seq}` | {d} | {med} | {best} | {seen} |")
    print(f"\nWild-type median for comparison: {wt_median:.3f} nM.\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spr", type=Path, default=DEFAULT_SPR, help=f"default: {DEFAULT_SPR}")
    ap.add_argument(
        "--designs",
        nargs="*",
        default=["WGPPYFYAMD", "WASSSEYYFD"],
        help="design CDRH3 windows to place in the distribution (default: the two refolded validation designs)",
    )
    args = ap.parse_args()

    if not args.spr.exists():
        sys.exit(f"{args.spr} not found — sync hackathon-data/ first (see TASKS.md)")

    measured = load_measured(args.spr)
    print(f"# Trastuzumab SPR — what it takes to beat wild-type\n")
    print(f"Source: `{args.spr}`, Chinery et al. doi:10.1101/2024.03.26.586756\n")
    wt_median = report_assay_floor(measured)
    report_substitution_cost(measured, wt_median)
    report_novelty_tension(measured, wt_median)
    report_designs(measured, args.designs, wt_median)


if __name__ == "__main__":
    main()
