#!/usr/bin/env python3
"""Build an abguard submission table from a batch of RFD4-Proteina designs.

Each design only varies a 10-residue CDRH3 window, so a full submission row is
the wild-type Fv with that window spliced in. Emits the five columns
`team-validate` wants, plus a wild-type control row that MUST be rejected --
if the control passes, the novelty gate is not actually on.

Validate the result on the VM with:
    abguard team-validate <out>.csv --reference-db ~/team15-claude/her2_novelty.duckdb

Usage:
    make_submission.py --designs runs/<run>/design/generation/inference_samples \
                       --reference inputs/12-her2/her2_fv_1n8z.cif \
                       --out runs/<run>/submission_all.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

WILD_TYPE_WINDOW = "WGGDGFYAMD"


def parse_cif_chains(path: Path) -> dict[str, dict[int, str]]:
    """Return {chain: {seq_id: one_letter}} from an mmCIF's CA atoms."""
    chains: dict[str, dict[int, str]] = {}
    with path.open() as fh:
        header: list[str] = []
        in_loop = False
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                header.append(s.split(".", 1)[1])
                in_loop = True
                continue
            if in_loop:
                if s.startswith(("#", "loop_", "_")) or not s:
                    in_loop = False
                    continue
                f = s.split()
                if len(f) < len(header):
                    continue
                row = dict(zip(header, f))
                if row.get("label_atom_id") != "CA":
                    continue
                chain = row.get("label_asym_id") or row.get("auth_asym_id")
                try:
                    seq_id = int(row.get("label_seq_id") or row.get("auth_seq_id"))
                except (TypeError, ValueError):
                    continue
                aa = AA3TO1.get((row.get("label_comp_id") or "").upper())
                if chain and aa:
                    chains.setdefault(chain, {})[seq_id] = aa
    return chains


def chain_sequence(residues: dict[int, str]) -> str:
    return "".join(residues[k] for k in sorted(residues))


def window_of(residues: dict[int, str], start: int, end: int) -> str:
    return "".join(residues.get(i, "-") for i in range(start, end + 1))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--designs", type=Path, required=True, help="folder of <id>/<id>_generated.cif")
    ap.add_argument("--reference", type=Path, required=True, help="wild-type Fv mmCIF (1N8Z-derived)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--chain", default="B", help="heavy chain id (default B)")
    ap.add_argument("--light-chain", default="A", help="light chain id (default A)")
    ap.add_argument("--window", nargs=2, type=int, default=(99, 108), metavar=("START", "END"))
    ap.add_argument("--description", default="RFD4-Proteina CDRH3 redesign vs HER2 domain IV, B99-108 free")
    args = ap.parse_args()

    start, end = args.window
    ref = parse_cif_chains(args.reference)
    if args.chain not in ref or args.light_chain not in ref:
        sys.exit(f"reference lacks chain {args.chain} or {args.light_chain}")

    wt_heavy = chain_sequence(ref[args.chain])
    wt_light = chain_sequence(ref[args.light_chain])
    wt_window = window_of(ref[args.chain], start, end)

    if wt_window not in wt_heavy:
        sys.exit(f"wild-type window {wt_window} not found in heavy chain")
    if wt_heavy.count(wt_window) != 1:
        sys.exit(f"window {wt_window} is ambiguous in the heavy chain")
    if wt_window != WILD_TYPE_WINDOW:
        print(f"warning: window {wt_window} != expected {WILD_TYPE_WINDOW}", file=sys.stderr)

    rows = []
    skipped = []
    for cif in sorted(args.designs.glob("*/*_generated.cif")):
        design_id = cif.parent.name
        chains = parse_cif_chains(cif)
        if args.chain not in chains:
            skipped.append((design_id, "no heavy chain"))
            continue
        win = window_of(chains[args.chain], start, end)
        if "-" in win or len(win) != end - start + 1:
            skipped.append((design_id, f"incomplete window {win!r}"))
            continue
        subs = sum(a != b for a, b in zip(win, wt_window))
        rows.append({
            "sequence_id": design_id,
            "kind": "paired",
            "heavy_sequence": wt_heavy.replace(wt_window, win),
            "light_sequence": wt_light,
            "design_description": f"{args.description}; window {win}; {subs} subs from wild-type",
        })

    # Control: wild-type itself. abguard MUST reject this row. If it passes, the
    # novelty gate is off and every other verdict in the file is meaningless.
    rows.append({
        "sequence_id": "control_wildtype_trastuzumab",
        "kind": "paired",
        "heavy_sequence": wt_heavy,
        "light_sequence": wt_light,
        "design_description": "wild-type trastuzumab control, MUST be rejected",
    })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {len(rows)} rows ({len(rows) - 1} designs + 1 control) to {args.out}")
    print(f"wild-type window {args.chain}{start}-{end}: {wt_window}")
    for design_id, why in skipped:
        print(f"  skipped {design_id}: {why}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
