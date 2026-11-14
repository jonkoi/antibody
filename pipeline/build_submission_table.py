#!/usr/bin/env python3
"""Build a valid `abguard` submission table from a list of CDRH3 designs.

`abguard team-validate` / `novelty-report` / `liability-report` all read the same
table, which must carry exactly these columns:

    sequence_id, kind, heavy_sequence, light_sequence, design_description

Rows are partitioned on `kind` (`vhh` / `scfv` / `paired`) and each partition is
validated against its own dataframely schema. Our HER2 designs are `paired`
(heavy + light), the light chain being wild-type trastuzumab in every design.

DESIGN WINDOW vs IMGT CDRH3 -- read this before trusting a distance
------------------------------------------------------------------
Our design window is the 10 residues at heavy-chain positions 99-108:

    ...TAVYYCSR[WGGDGFYAMD]YWGQGTLVTVSS
                ^^^^^^^^^^   <- what the design model is allowed to change

abguard numbers with IMGT, and its CDRH3 is the 13-residue span 105-117, i.e.
heavy-chain positions 97-109 of this construct:

    ...TAVYYC[SRWGGDGFYAMDY]WGQGTLVTVSS
               ^^^^^^^^^^^^^  <- what abguard extracts, scores and compares

So abguard sees three residues we never designed: S97 and R98 on the N-side and
Y109 on the C-side. They are invariant in every provided table row and in every
design we emit, so Hamming distance over the 13-mer equals Hamming distance over
the 10-mer *for the trastuzumab library* -- but they do count toward liability
motifs and toward comparisons against non-trastuzumab reference entries.

Usage
-----
    # from explicit 10-mers
    build_submission_table.py -o out.csv WGPPYFYAMD WGSLTFYAMD

    # from a CSV column (default column name: cdrh3)
    build_submission_table.py -o out.csv --from-csv data/legal_3mutants_ranked.csv \
                              --column cdrh3 --sample 50 --seed 0

    # from a file of one 10-mer per line, or stdin
    build_submission_table.py -o out.csv --from-lines probes.txt

Each row's `design_description` records the substitutions against wild-type, so
abguard's output can be read back against the design rationale.
"""
from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

# Wild-type trastuzumab, as given in the challenge data.
WT_HEAVY = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
WT_LIGHT = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)

# Design window: heavy-chain positions 99-108 (1-based, inclusive) == the
# 10 residues abguard sees at IMGT CDRH3 offsets 2..11 of its 13-mer.
WINDOW_START = 99
WINDOW_LEN = 10
# Zero-based slice bounds into WT_HEAVY.
_LO = WINDOW_START - 1
_HI = _LO + WINDOW_LEN

WT_CDRH3_10 = WT_HEAVY[_LO:_HI]
# What abguard actually extracts and compares (IMGT 105-117).
WT_CDRH3_IMGT13 = WT_HEAVY[_LO - 2 : _HI + 1]

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

assert WT_CDRH3_10 == "WGGDGFYAMD", WT_CDRH3_10
assert WT_CDRH3_IMGT13 == "SRWGGDGFYAMDY", WT_CDRH3_IMGT13
assert len(WT_HEAVY) == 120 and len(WT_LIGHT) == 107

COLUMNS = [
    "sequence_id",
    "kind",
    "heavy_sequence",
    "light_sequence",
    "design_description",
]


def graft(cdrh3: str) -> str:
    """Substitute `cdrh3` into the wild-type heavy chain design window."""
    cdrh3 = cdrh3.strip().upper()
    if len(cdrh3) != WINDOW_LEN:
        raise ValueError(
            f"{cdrh3!r} is {len(cdrh3)} residues; the design window is {WINDOW_LEN}. "
            "Pass the 10-mer, not abguard's 13-residue IMGT CDRH3."
        )
    bad = sorted(set(cdrh3) - AMINO_ACIDS)
    if bad:
        raise ValueError(f"{cdrh3!r} has non-amino-acid characters: {bad}")
    return WT_HEAVY[:_LO] + cdrh3 + WT_HEAVY[_HI:]


def substitutions(cdrh3: str) -> list[str]:
    """Substitutions against wild-type, in heavy-chain numbering (e.g. 'G101P')."""
    cdrh3 = cdrh3.strip().upper()
    return [
        f"{wt}{WINDOW_START + i}{mut}"
        for i, (wt, mut) in enumerate(zip(WT_CDRH3_10, cdrh3))
        if wt != mut
    ]


def describe(cdrh3: str, note: str | None = None) -> str:
    subs = substitutions(cdrh3)
    n = len(subs)
    parts = [
        f"CDRH3 {cdrh3} grafted at heavy 99-108",
        f"{n} substitution{'' if n == 1 else 's'} vs trastuzumab"
        + (f" ({' '.join(subs)})" if subs else " (wild-type)"),
    ]
    if note:
        parts.append(note)
    return "; ".join(parts)


def build_rows(
    designs: list[tuple[str, str, str | None]],
) -> list[dict[str, str]]:
    """designs: (sequence_id, cdrh3_10mer, note) -> submission rows."""
    rows = []
    seen: set[str] = set()
    for seq_id, cdrh3, note in designs:
        if seq_id in seen:
            raise ValueError(f"duplicate sequence_id {seq_id!r}; abguard requires unique IDs")
        seen.add(seq_id)
        rows.append(
            {
                "sequence_id": seq_id,
                "kind": "paired",
                "heavy_sequence": graft(cdrh3),
                "light_sequence": WT_LIGHT,
                "design_description": describe(cdrh3, note),
            }
        )
    return rows


def write_table(rows: list[dict[str, str]], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    delim = "\t" if out.suffix.lower() in {".tsv", ".tab"} else ","
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter=delim)
        writer.writeheader()
        writer.writerows(rows)


def _load_csv_column(path: Path, column: str) -> list[str]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if column not in (reader.fieldnames or []):
            raise SystemExit(
                f"{path} has no column {column!r}; available: {reader.fieldnames}"
            )
        return [row[column] for row in reader if row.get(column)]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cdrh3", nargs="*", help="CDRH3 10-mers given directly")
    ap.add_argument("--from-csv", type=Path, help="read 10-mers from a CSV column")
    ap.add_argument("--column", default="cdrh3", help="column name for --from-csv (default: cdrh3)")
    ap.add_argument("--from-lines", type=Path, help="read 10-mers one per line ('-' for stdin)")
    ap.add_argument("--sample", type=int, help="take a random sample of this many")
    ap.add_argument("--seed", type=int, default=0, help="sampling seed (default: 0)")
    ap.add_argument("--prefix", default="design", help="sequence_id prefix (default: design)")
    ap.add_argument("--include-wt", action="store_true", help="prepend a wild-type trastuzumab row")
    ap.add_argument("-o", "--output", type=Path, required=True, help="output table (.csv / .tsv)")
    args = ap.parse_args(argv)

    seqs: list[str] = [s.strip().upper() for s in args.cdrh3]
    if args.from_csv:
        seqs += [s.strip().upper() for s in _load_csv_column(args.from_csv, args.column)]
    if args.from_lines:
        fh = sys.stdin if str(args.from_lines) == "-" else args.from_lines.open()
        with fh:
            seqs += [ln.strip().upper() for ln in fh if ln.strip() and not ln.startswith("#")]
    if not seqs:
        raise SystemExit("no CDRH3s given; pass them as arguments or use --from-csv/--from-lines")

    if args.sample is not None and args.sample < len(seqs):
        seqs = random.Random(args.seed).sample(seqs, args.sample)

    designs: list[tuple[str, str, str | None]] = []
    if args.include_wt:
        designs.append(("wt_trastuzumab", WT_CDRH3_10, "wild-type control"))
    width = max(3, len(str(len(seqs))))
    for i, s in enumerate(seqs, start=1):
        designs.append((f"{args.prefix}_{i:0{width}d}_{s}", s, None))

    rows = build_rows(designs)
    write_table(rows, args.output)
    print(f"wrote {len(rows)} rows to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
