#!/usr/bin/env python3
"""Score legal candidates against the substitutions found in measured winners.

Of 182 non-wild-type CDRH3 variants with a measured KD, only seven beat
trastuzumab's median affinity of 0.146 nM. Those seven are the most direct
evidence in the entire challenge about which positions tolerate change and which
substitutions help, and nothing else in our pipeline uses them.

Two signals come out of them:

  POSITIONS THEY CHANGE     3, 4, 6, 8, 9   -- tolerant, and sometimes improving
  POSITIONS THEY NEVER TOUCH  1, 2, 5, 7, 10 -- untouched in all seven winners

The second is as informative as the first. Position 1 is the conserved tryptophan
at the base of the loop; positions 2, 5, 7 and 10 are left alone by every variant
that improved on trastuzumab, which is evidence they are load-bearing.

This matters because it contradicts a heuristic we had been applying. Wild-type is
glycine-rich (positions 2, 3, 5) and replacing a glycine with proline rigidifies the
loop, which should help thermostability. Our developability shortlist put proline at
position 5 in 1,412 of 1,927 candidates -- but position 5 is never changed by any
measured winner. **Position 3 is the glycine the winners do change (G3S, twice), so
it is the defensible rigidification site; position 5 is not.**

The score is deliberately simple and interpretable: reward substitutions at
winner-tolerant positions, reward the specific residues winners chose, and penalise
any change at a position all seven left alone.

Usage:
    winner_prior.py --spr <spr.csv> --candidates <legal.csv> [--out scored.csv]
"""
from __future__ import annotations

import argparse
import collections
import csv
from pathlib import Path

WILD_TYPE = "WGGDGFYAMD"

# weights are ordinal, chosen to express a clear ordering rather than fitted:
# an exact winner substitution > any change at a tolerant position > no change,
# and a change at a never-touched position is a real penalty.
W_EXACT_MATCH = 3.0       # same residue a winner used at that position
W_TOLERANT_POSITION = 1.0  # position a winner changed, different residue
W_CONSERVED_PENALTY = -2.0  # position no winner ever changed


def load_winners(spr_csv: Path, wild_type: str) -> tuple[list[tuple[str, float]], float]:
    rows = list(csv.DictReader(spr_csv.open()))
    measured = [
        (r["cdrh3"].strip().upper(), float(r["kd_nM"]))
        for r in rows
        if r["kd_nM"].strip() and len(r["cdrh3"].strip()) == len(wild_type)
    ]
    wt_values = sorted(v for s, v in measured if s == wild_type)
    if not wt_values:
        raise ValueError("wild-type not found in the SPR table; cannot set a baseline")
    baseline = wt_values[len(wt_values) // 2]
    winners = sorted(((s, v) for s, v in measured if s != wild_type and v < baseline),
                     key=lambda x: x[1])
    return winners, baseline


def build_prior(winners: list[tuple[str, float]], wild_type: str):
    favoured: dict[int, collections.Counter] = collections.defaultdict(collections.Counter)
    for seq, _ in winners:
        for i, (a, b) in enumerate(zip(seq, wild_type)):
            if a != b:
                favoured[i][a] += 1
    tolerant = set(favoured)
    conserved = set(range(len(wild_type))) - tolerant
    return favoured, tolerant, conserved


def score(seq: str, wild_type: str, favoured, tolerant, conserved) -> tuple[float, list[str]]:
    total = 0.0
    notes: list[str] = []
    for i, (a, b) in enumerate(zip(seq, wild_type)):
        if a == b:
            continue
        if i in conserved:
            total += W_CONSERVED_PENALTY
            notes.append(f"!{b}{i + 1}{a}")       # changed a position winners never touch
        elif a in favoured[i]:
            total += W_EXACT_MATCH * favoured[i][a]
            notes.append(f"*{b}{i + 1}{a}")       # exactly a winner's substitution
        else:
            total += W_TOLERANT_POSITION
            notes.append(f"{b}{i + 1}{a}")
    return total, notes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--spr", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--wild-type", default=WILD_TYPE)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--top", type=int, default=25)
    args = parser.parse_args()

    wt = args.wild_type
    winners, baseline = load_winners(args.spr, wt)
    favoured, tolerant, conserved = build_prior(winners, wt)

    print(f"trastuzumab median KD {baseline} nM; {len(winners)} variants beat it\n")
    for seq, value in winners:
        subs = " ".join(f"{wt[i]}{i+1}{seq[i]}" for i in range(len(wt)) if seq[i] != wt[i])
        print(f"  {seq}  {value:>8.4f} nM   {subs}")

    print(f"\ntolerant positions (changed by a winner): "
          f"{', '.join(str(i + 1) for i in sorted(tolerant))}")
    for i in sorted(tolerant):
        print(f"    pos {i + 1} ({wt[i]}): " +
              ", ".join(f"{a}x{c}" for a, c in favoured[i].most_common()))
    print(f"conserved positions (untouched by all winners): "
          f"{', '.join(str(i + 1) for i in sorted(conserved))}")

    rows = list(csv.DictReader(args.candidates.open()))
    for row in rows:
        value, notes = score(row["cdrh3"], wt, favoured, tolerant, conserved)
        row["winner_score"] = f"{value:.1f}"
        row["winner_notes"] = " ".join(notes)
        row["touches_conserved"] = int(any(n.startswith("!") for n in notes))
    rows.sort(key=lambda r: float(r["winner_score"]), reverse=True)

    clean = [r for r in rows if r["touches_conserved"] == 0]
    print(f"\n{len(rows):,} candidates scored; {len(clean):,} change only tolerant positions")

    print(f"\ntop {args.top} (none touching a conserved position):")
    print(f"  {'cdrh3':<12}{'winner':>7}{'seqmodel':>10}  notes   (* = exact winner substitution)")
    for row in clean[:args.top]:
        seqm = row.get("score", "")
        seqm = f"{float(seqm):.2f}" if seqm else "-"
        print(f"  {row['cdrh3']:<12}{float(row['winner_score']):>7.1f}{seqm:>10}  {row['winner_notes']}")

    if args.out:
        with args.out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote {len(rows):,} rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
