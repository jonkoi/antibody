#!/usr/bin/env python3
"""Rank CDRH3 candidates by predicted HER2 binding, trained on the FACS screen.

The structural pipeline cannot fold 51,417 candidates, and the confidence scores
it would rank them with are poorly calibrated anyway (see calibrate_iptm.py).
This provides the cheap sequence-level filter that cuts the legal set down to a
foldable shortlist, trained on the only ranking signal in the challenge that is
grounded in measured binding: the 524,346-variant FACS screen.

MODEL
A position-specific log-odds score (naive Bayes over one-hot positions):

    score(s) = sum_i [ log P(s_i | high) - log P(s_i | not high) ]

with Laplace smoothing. 200 parameters against ~524k labelled examples, so it is
heavily over-determined and needs no regularisation beyond the prior. It assumes
positions contribute independently, which is false for a binding loop -- epistasis
is real. It is therefore a *filter*, not a final ranking: it narrows 51,417 to a
shortlist that a structural method can afford to evaluate properly.

Labels follow the convention in the published splits: high == binder, medium and
low == non-binder.

HONEST EVALUATION
Reported on a held-out split, and additionally stratified by edit distance from
wild-type, because performance in the distance-3 shell is what actually matters --
that is where every legal candidate lives. A model that scores well overall but
poorly at distance 3 is useless to us.

Usage:
    rank_candidates.py --facs <facs.csv.zip> --candidates <legal.csv> [--out ranked.csv]
    rank_candidates.py --facs <facs.csv.zip> --evaluate-only
"""
from __future__ import annotations

import argparse
import csv
import io
import math
import random
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
PSEUDOCOUNT = 1.0
POSITIVE_CLASS = "high"


def load_facs(path: Path) -> list[tuple[str, int, int]]:
    """(cdrh3, label, edit_distance) from the zipped FACS table."""
    rows: list[tuple[str, int, int]] = []
    with zipfile.ZipFile(path) as archive:
        with archive.open(archive.namelist()[0]) as handle:
            for row in csv.DictReader(io.TextIOWrapper(handle, "utf-8")):
                seq = row["cdrh3"].strip().upper()
                label = int(row["binding_class"].strip().lower() == POSITIVE_CLASS)
                try:
                    distance = int(row["edit_distance"])
                except (KeyError, ValueError):
                    distance = -1
                rows.append((seq, label, distance))
    return rows


def fit(rows: list[tuple[str, int, int]], length: int) -> list[dict[str, float]]:
    """Per-position log-odds of each residue, high versus rest."""
    high = [defaultdict(float) for _ in range(length)]
    rest = [defaultdict(float) for _ in range(length)]
    n_high = n_rest = 0

    for seq, label, _ in rows:
        if len(seq) != length:
            continue
        target = high if label else rest
        for i, residue in enumerate(seq):
            target[i][residue] += 1
        if label:
            n_high += 1
        else:
            n_rest += 1

    if not n_high or not n_rest:
        raise ValueError("need both classes present to fit")

    denom_high = n_high + PSEUDOCOUNT * len(ALPHABET)
    denom_rest = n_rest + PSEUDOCOUNT * len(ALPHABET)

    weights: list[dict[str, float]] = []
    for i in range(length):
        column = {}
        for residue in ALPHABET:
            p_high = (high[i][residue] + PSEUDOCOUNT) / denom_high
            p_rest = (rest[i][residue] + PSEUDOCOUNT) / denom_rest
            column[residue] = math.log(p_high) - math.log(p_rest)
        weights.append(column)
    return weights


def score(seq: str, weights: list[dict[str, float]]) -> float:
    return sum(weights[i].get(residue, 0.0) for i, residue in enumerate(seq))


def auc(positive: list[float], negative: list[float]) -> float:
    """Tie-corrected Mann-Whitney AUC."""
    if not positive or not negative:
        return float("nan")
    merged = sorted([(v, 1) for v in positive] + [(v, 0) for v in negative])
    ranks: list[float] = [0.0] * len(merged)
    i = 0
    while i < len(merged):
        j = i
        while j + 1 < len(merged) and merged[j + 1][0] == merged[i][0]:
            j += 1
        mid = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = mid
        i = j + 1
    rank_sum = sum(r for r, (_, lab) in zip(ranks, merged) if lab == 1)
    n_pos, n_neg = len(positive), len(negative)
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def evaluate(rows: list[tuple[str, int, int]], length: int, seed: int = 0) -> None:
    rng = random.Random(seed)
    shuffled = rows[:]
    rng.shuffle(shuffled)
    cut = int(0.8 * len(shuffled))
    train, test = shuffled[:cut], shuffled[cut:]

    weights = fit(train, length)
    scored = [(score(s, weights), lab, dist) for s, lab, dist in test if len(s) == length]

    overall = auc([v for v, lab, _ in scored if lab], [v for v, lab, _ in scored if not lab])
    print(f"held-out AUC (n={len(scored):,}): {overall:.3f}")

    print("\nby edit distance from wild-type:")
    print(f"  {'dist':>4}  {'n':>8}  {'%high':>6}  {'AUC':>6}")
    buckets: dict[int, list[tuple[float, int]]] = defaultdict(list)
    for value, lab, dist in scored:
        buckets[dist].append((value, lab))
    for dist in sorted(b for b in buckets if b >= 0):
        values = buckets[dist]
        pos = [v for v, lab in values if lab]
        neg = [v for v, lab in values if not lab]
        rate = 100 * len(pos) / len(values)
        a = auc(pos, neg)
        flag = "   <- the legal shell" if dist == 3 else ""
        print(f"  {dist:>4}  {len(values):>8,}  {rate:>5.1f}%  {a:>6.3f}{flag}")

    print("\nmost discriminative positions (spread of log-odds across residues):")
    full = fit(rows, length)
    spread = sorted(
        ((max(col.values()) - min(col.values()), i) for i, col in enumerate(full)),
        reverse=True)
    for value, i in spread:
        best = max(full[i], key=full[i].get)
        worst = min(full[i], key=full[i].get)
        print(f"  position {i + 1:>2}: spread {value:5.2f}   favours {best}, disfavours {worst}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--facs", type=Path, required=True)
    parser.add_argument("--candidates", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--length", type=int, default=10)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--evaluate-only", action="store_true")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    rows = load_facs(args.facs)
    n_pos = sum(lab for _, lab, _ in rows)
    print(f"training rows: {len(rows):,}  ({n_pos:,} high = {100 * n_pos / len(rows):.1f}%)\n")

    evaluate(rows, args.length, args.seed)

    if args.evaluate_only:
        return 0
    if not args.candidates:
        print("\nno --candidates given; nothing to rank", file=sys.stderr)
        return 0

    weights = fit(rows, args.length)
    with args.candidates.open() as handle:
        candidates = list(csv.DictReader(handle))
    for row in candidates:
        row["score"] = f"{score(row['cdrh3'], weights):.4f}"
    candidates.sort(key=lambda r: float(r["score"]), reverse=True)

    print(f"\n\nranked {len(candidates):,} candidates")
    print(f"\ntop {args.top}:")
    header = f"  {'rank':>4}  {'cdrh3':<12}{'score':>8}  {'substitutions':<22}{'DG':>3}{'P/G':>4}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for n, row in enumerate(candidates[:args.top], 1):
        print(f"  {n:>4}  {row['cdrh3']:<12}{float(row['score']):>8.2f}  "
              f"{row.get('substitutions', ''):<22}"
              f"{'Y' if row.get('removes_dg') == '1' else '-':>3}"
              f"{'Y' if row.get('pro_for_gly') == '1' else '-':>4}")

    both = [r for r in candidates
            if r.get("removes_dg") == "1" and r.get("pro_for_gly") == "1"
            and r["cdrh3"][0] == "W"]
    if both:
        print(f"\ntop {args.top} among the {len(both):,} that also remove DG, place Pro at a Gly,")
        print("and retain the conserved Trp:")
        for n, row in enumerate(both[:args.top], 1):
            print(f"  {n:>4}  {row['cdrh3']:<12}{float(row['score']):>8.2f}  "
                  f"{row.get('substitutions', '')}")

    if args.out:
        with args.out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(candidates[0]))
            writer.writeheader()
            writer.writerows(candidates)
        print(f"\nwrote {len(candidates):,} ranked rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
