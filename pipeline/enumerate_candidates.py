#!/usr/bin/env python3
"""Enumerate CDRH3 variants that satisfy the hackathon's novelty requirement.

abguard requires a submitted CDRH3 to sit at Hamming distance >= 3 from
trastuzumab's own loop and from every CDRH3 in the provided datasets. Separately,
the FACS data shows the strong-binder rate falling with substitution count
(90.0% at 3, 74.6% at 4, 59.4% at 5), so distance 3 is simultaneously the
smallest legal change and the one that best preserves binding. This enumerates
the complete legal set at that distance, exactly.

ALGORITHM
Brute force is 823,080 x 525,063 comparisons. Three observations collapse it:

1. A candidate C at distance 3 from wild-type differs from WT only at its three
   mutated positions P, so for any reference R the distance splits as

       d(C, R) = A + (3 - m)
         A = positions OUTSIDE P where R differs from WT  (independent of C)
         m = positions INSIDE  P where R agrees with C

   C is illegal exactly when m >= A + 1, and any reference with A >= 3 is
   irrelevant for that P.

2. A depends on R only through which positions R differs from WT -- a 10-bit
   mask. There are at most 1024 distinct masks, so A is computed 1024 times per
   position set, not 525,063 times.

3. Each surviving reference generates its forbidden triples directly by choosing
   which positions match it, rather than by scanning all 6,859 candidates.

Exact, no approximation, pure standard library.

Usage:
    enumerate_candidates.py --facs <path.zip> --spr <path.csv> [--out survivors.csv]
"""
from __future__ import annotations

import argparse
import csv
import io
import itertools
import sys
import time
import zipfile
from collections import defaultdict
from pathlib import Path

WILD_TYPE = "WGGDGFYAMD"
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
MIN_DISTANCE = 3


def load_references(facs_zip: Path, spr_csv: Path, length: int) -> tuple[set[str], dict[int, int]]:
    """Unique CDRH3s of the given length from both provided tables.

    References of a different length are excluded by construction: Hamming
    distance is undefined across lengths, and any edit path between sequences of
    different length contains at least one insertion or deletion, so a 9- or
    16-residue reference cannot sit within 3 substitutions of a 10-residue
    candidate. The counts excluded are returned so the report can state them.
    """
    refs: set[str] = set()
    skipped: dict[int, int] = defaultdict(int)

    def take(seq: str) -> None:
        seq = seq.strip().upper()
        if len(seq) == length:
            refs.add(seq)
        else:
            skipped[len(seq)] += 1

    with zipfile.ZipFile(facs_zip) as archive:
        with archive.open(archive.namelist()[0]) as handle:
            for row in csv.DictReader(io.TextIOWrapper(handle, "utf-8")):
                take(row["cdrh3"])

    with spr_csv.open() as handle:
        for row in csv.DictReader(handle):
            take(row["cdrh3"])

    return refs, dict(skipped)


def enumerate_legal(wild_type: str, refs: set[str], distance: int,
                    allowed: list[int] | None = None) -> list[dict]:
    """`allowed` restricts which positions may be substituted (0-based).

    Used to confine the search to the positions that measured winners actually
    change, which both shrinks the space and keeps every candidate inside the
    region the data supports.
    """
    n = len(wild_type)
    alts = {i: [c for c in ALPHABET if c != wild_type[i]] for i in range(n)}
    sites = sorted(allowed) if allowed is not None else list(range(n))

    # group references by which positions differ from wild-type
    by_mask: dict[int, list[str]] = defaultdict(list)
    for ref in refs:
        mask = 0
        for i, (a, b) in enumerate(zip(ref, wild_type)):
            if a != b:
                mask |= 1 << i
        by_mask[mask].append(ref)
    by_mask[0].append(wild_type)  # wild-type is itself a reference

    survivors: list[dict] = []
    gly = {i for i, c in enumerate(wild_type) if c == "G"}
    wt_has_dg = "DG" in wild_type

    for positions in itertools.combinations(sites, distance):
        p_mask = sum(1 << p for p in positions)

        # A per mask, keeping only masks that can still produce a violation
        live: list[tuple[int, list[str]]] = []
        for mask, group in by_mask.items():
            a = bin(mask & ~p_mask).count("1")
            if a < distance:
                live.append((a, group))

        forbidden: set[tuple[str, ...]] = set()
        for a, group in live:
            needed = a + 1          # matches inside P required for a violation
            if needed > distance:
                continue
            for ref in group:
                r = [ref[p] for p in positions]
                # choose which positions match the reference; the rest range freely
                for size in range(needed, distance + 1):
                    for matched in itertools.combinations(range(distance), size):
                        if any(r[i] == wild_type[positions[i]] for i in matched):
                            continue  # candidate must differ from WT there, so cannot match
                        free = [i for i in range(distance) if i not in matched]
                        choices = [alts[positions[i]] for i in free]
                        for combo in itertools.product(*choices) if free else [()]:
                            triple = [None] * distance
                            for i in matched:
                                triple[i] = r[i]
                            for i, residue in zip(free, combo):
                                triple[i] = residue
                            forbidden.add(tuple(triple))

        for combo in itertools.product(*(alts[p] for p in positions)):
            if combo in forbidden:
                continue
            seq = list(wild_type)
            for pos, residue in zip(positions, combo):
                seq[pos] = residue
            variant = "".join(seq)
            survivors.append({
                "cdrh3": variant,
                "substitutions": " ".join(
                    f"{wild_type[p]}{p + 1}{variant[p]}" for p in positions),
                "removes_dg": int(wt_has_dg and "DG" not in variant),
                "pro_for_gly": int(any(variant[p] == "P" and p in gly for p in positions)),
            })

    return survivors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--facs", type=Path, required=True)
    parser.add_argument("--spr", type=Path, required=True)
    parser.add_argument("--wild-type", default=WILD_TYPE)
    parser.add_argument("--distance", type=int, default=MIN_DISTANCE)
    parser.add_argument("--positions", type=str,
                        help="restrict substitutions to these 1-based positions, "
                             "e.g. '3,4,6,8,9' for the winner-tolerant set")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    wt = args.wild_type
    allowed = None
    if args.positions:
        allowed = [int(p) - 1 for p in args.positions.split(",")]
        if min(allowed) < 0 or max(allowed) >= len(wt):
            print("error: --positions out of range", file=sys.stderr)
            return 1
    started = time.time()
    refs, skipped = load_references(args.facs, args.spr, len(wt))
    print(f"references at length {len(wt)}: {len(refs):,}   "
          f"(excluded by length: {skipped or 'none'})")

    sites = allowed if allowed is not None else list(range(len(wt)))
    space = len(list(itertools.combinations(sites, args.distance))) \
        * (len(ALPHABET) - 1) ** args.distance
    print(f"candidate space at distance {args.distance}: {space:,}")

    survivors = enumerate_legal(wt, refs, args.distance, allowed)
    elapsed = time.time() - started

    dg = sum(s["removes_dg"] for s in survivors)
    pro = sum(s["pro_for_gly"] for s in survivors)
    both = [s for s in survivors if s["removes_dg"] and s["pro_for_gly"]]

    pct = 100 * len(survivors) / space if space else 0
    print(f"\nclear the novelty gate: {len(survivors):,}  ({pct:.1f}% of the space)")
    print(f"  also remove DG:       {dg:,}")
    print(f"  also Pro-for-Gly:     {pro:,}")
    print(f"  both:                 {len(both):,}")
    print(f"elapsed: {elapsed:.1f}s")

    if both:
        print("\nexamples satisfying both developability criteria:")
        for row in both[:10]:
            print(f"  {row['cdrh3']}   {row['substitutions']}")

    if args.out and survivors:
        with args.out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(survivors[0]))
            writer.writeheader()
            writer.writerows(survivors)
        print(f"\nwrote {len(survivors):,} rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
