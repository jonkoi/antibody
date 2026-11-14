#!/usr/bin/env python3
"""Calibrate an ipTM cut-off against measured non-binders (Lane B).

The question: if we rank designs on interface confidence, what fraction of known
NON-binders would clear the threshold we pick?

Data: Smorodina et al. (zenodo 18390239) folded a panel of nanobodies against a
panel of antigens with AlphaFold3, Boltz-2 and Chai-1. In each
`*_matrix_clean.csv`, rows are nanobody PDB ids and columns antigen PDB ids, so
a cell is one predicted complex. A cell whose row id equals its column id is the
nanobody with its own antigen -- a COGNATE (real) pair. Every other cell is a
SHUFFLED pair: a nanobody paired with an antigen it is not known to bind.

Shuffled pairs are the negatives here. That is the authors' own design, but it is
an assumption, not an SPR measurement: a shuffled pair is presumed non-binding
because nobody reported it binding. The SPR tables in
sigma2-ab-ag-champloo/spr-supplementary/ measure only 29 pairs and their README
states plainly that "no reportable fit" must NOT be read as confirmed
non-binding, so they cannot carry a false-positive rate.

No third-party dependencies, by design: stdlib only.

Usage:
    python3 pipeline/calibrate_confidence.py
    python3 pipeline/calibrate_confidence.py --zip <path to iptm_confidence_scores.zip>
"""

from __future__ import annotations

import argparse
import csv
import io
import zipfile
from pathlib import Path

DEFAULT_ZIP = Path("hackathon-data/zenodo-18390239-ab-ag-complexes/iptm_confidence_scores.zip")

# ipTM values worth reporting a false-positive rate for: 0.7 is the threshold
# quoted in our own README, 0.82-0.90 is the range our n=100 shortlist occupies.
REPORT_AT = [0.5, 0.6, 0.7, 0.75, 0.8, 0.82, 0.85, 0.9, 0.95]
TARGET_FPRS = [0.20, 0.10, 0.05, 0.01]


def load_matrix(raw: bytes) -> tuple[list[tuple[str, str, float]], dict]:
    """Return [(nanobody, antigen, iptm)] plus a shape/parse report."""
    reader = csv.reader(io.StringIO(raw.decode()))
    header = next(reader)
    antigens = header[1:]
    cells: list[tuple[str, str, float]] = []
    rows = 0
    blank = 0
    unparsable = 0
    for row in reader:
        if not row or not row[0].strip():
            continue
        rows += 1
        nanobody = row[0].strip()
        for antigen, value in zip(antigens, row[1:]):
            text = value.strip()
            if not text or text.upper() in {"NA", "NAN", "NONE"}:
                blank += 1
                continue
            try:
                cells.append((nanobody, antigen.strip(), float(text)))
            except ValueError:
                unparsable += 1
    report = {
        "rows": rows,
        "cols": len(antigens),
        "cells": len(cells),
        "blank": blank,
        "unparsable": unparsable,
        "nanobodies": len({c[0] for c in cells}),
        "antigens": len({c[1] for c in cells}),
    }
    return cells, report


def split_labels(cells) -> tuple[list[float], list[float]]:
    """Cognate (row id == column id) vs shuffled."""
    cognate = [v for nb, ag, v in cells if nb == ag]
    shuffled = [v for nb, ag, v in cells if nb != ag]
    return cognate, shuffled


def auc(pos: list[float], neg: list[float]) -> float:
    """Mann-Whitney U / ROC AUC, ties counted as half. O(n log n)."""
    if not pos or not neg:
        return float("nan")
    merged = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    # average ranks over ties
    rank_sum_pos = 0.0
    i = 0
    while i < len(merged):
        j = i
        while j + 1 < len(merged) and merged[j + 1][0] == merged[i][0]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0  # ranks are 1-based
        rank_sum_pos += avg_rank * sum(1 for k in range(i, j + 1) if merged[k][1] == 1)
        i = j + 1
    n_pos, n_neg = len(pos), len(neg)
    return (rank_sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def auc_ci(a: float, n_pos: int, n_neg: int) -> tuple[float, float, float]:
    """Hanley-McNeil standard error and 95% CI for an AUC.

    Matters here because there are only ~90-150 cognate pairs, so an AUC gap of
    a few hundredths is not a real difference between models or features.
    """
    if not (0 < a < 1) or n_pos < 2 or n_neg < 2:
        return float("nan"), float("nan"), float("nan")
    q1 = a / (2 - a)
    q2 = 2 * a * a / (1 + a)
    var = (a * (1 - a) + (n_pos - 1) * (q1 - a * a) + (n_neg - 1) * (q2 - a * a)) / (n_pos * n_neg)
    se = var**0.5
    return se, a - 1.96 * se, a + 1.96 * se


def rate_at_or_above(values: list[float], threshold: float) -> float:
    return sum(1 for v in values if v >= threshold) / len(values) if values else float("nan")


def threshold_for_fpr(neg: list[float], target: float) -> float:
    """Lowest ipTM whose false-positive rate is <= target."""
    ordered = sorted(neg, reverse=True)
    allowed = int(len(ordered) * target)  # how many negatives we tolerate
    if allowed >= len(ordered):
        return float("-inf")
    # the (allowed+1)-th highest negative must fall below the cut
    return ordered[allowed] + 1e-12 if allowed < len(ordered) else float("inf")


def quantile(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    s = sorted(values)
    pos = q * (len(s) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (pos - lo)


def describe(name: str, values: list[float]) -> str:
    if not values:
        return f"  {name:<9} (none)"
    return (
        f"  {name:<9} n={len(values):<6} min={min(values):.3f} "
        f"p25={quantile(values, 0.25):.3f} median={quantile(values, 0.5):.3f} "
        f"p75={quantile(values, 0.75):.3f} max={max(values):.3f}"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zip", type=Path, default=DEFAULT_ZIP, help=f"default: {DEFAULT_ZIP}")
    args = ap.parse_args()

    if not args.zip.exists():
        raise SystemExit(
            f"not found: {args.zip}\n"
            "sync it first:  aws s3 sync s3://apheris-oslo-hackathon/data/ hackathon-data/"
        )

    print(f"source: {args.zip}")
    print(
        "labels: cognate = nanobody PDB id equals antigen PDB id (the matrix diagonal);\n"
        "        shuffled = every other cell, presumed non-binding by the study's design.\n"
    )

    with zipfile.ZipFile(args.zip) as z:
        names = sorted(n for n in z.namelist() if n.endswith(".csv"))
        per_model = {}
        for name in names:
            model = Path(name).name.replace("_matrix_clean.csv", "")
            cells, report = load_matrix(z.read(name))
            cognate, shuffled = split_labels(cells)
            per_model[model] = (cognate, shuffled)

            print(f"=== {model} ===")
            print(
                f"  matrix {report['rows']}x{report['cols']} -> {report['cells']} scored cells "
                f"({report['blank']} blank, {report['unparsable']} unparsable); "
                f"{report['nanobodies']} nanobodies x {report['antigens']} antigens"
            )
            print(describe("cognate", cognate))
            print(describe("shuffled", shuffled))
            a = auc(cognate, shuffled)
            se, lo, hi = auc_ci(a, len(cognate), len(shuffled))
            print(
                f"  ROC AUC (cognate vs shuffled) = {a:.3f}  "
                f"95% CI [{lo:.3f}, {hi:.3f}]  (Hanley-McNeil, SE {se:.3f})"
            )
            print()

    print("=" * 78)
    print("FALSE-POSITIVE RATE OF AN ipTM CUT-OFF")
    print("  share of known non-binding (shuffled) pairs scoring >= the cut-off,")
    print("  with the share of real pairs retained in brackets")
    print("=" * 78)
    header = "  ipTM  " + "".join(f"{m:>22}" for m in per_model)
    print(header)
    for t in REPORT_AT:
        row = f"  {t:<6.2f}"
        for model, (cognate, shuffled) in per_model.items():
            fpr = rate_at_or_above(shuffled, t)
            tpr = rate_at_or_above(cognate, t)
            row += f"{f'{fpr:6.1%} fp ({tpr:4.0%} kept)':>22}"
        print(row)

    print()
    print("=" * 78)
    print("CUT-OFF NEEDED FOR A GIVEN FALSE-POSITIVE RATE")
    print("=" * 78)
    print("  target FPR" + "".join(f"{m:>22}" for m in per_model))
    for target in TARGET_FPRS:
        row = f"  {target:<10.0%}"
        for model, (cognate, shuffled) in per_model.items():
            t = threshold_for_fpr(shuffled, target)
            tpr = rate_at_or_above(cognate, t)
            row += f"{f'ipTM>={t:.3f} ({tpr:4.0%} kept)':>22}"
        print(row)

    print()
    print("What this does and does not license:")
    print("  - These are unbound VHH-antigen predictions. Our refold hands the model the")
    print("    HER2 template with framework and CDRH1/H2 frozen, so our ipTM values are")
    print("    not on this scale and a threshold read off this table CANNOT be applied to")
    print("    our designs directly. What transfers is the ordering of models by AUC and")
    print("    the shape of the trade-off.")
    print("  - Negatives are shuffled pairs, not SPR-confirmed non-binders. The measured")
    print("    SPR set (spr-supplementary/, 29 rows) is too small and its authors state")
    print("    'no reportable fit' is not confirmed non-binding.")


if __name__ == "__main__":
    main()
