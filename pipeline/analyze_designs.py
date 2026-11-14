#!/usr/bin/env python3
"""Join RFD4-Proteina designs to their refold scores and the wild-type reference.

For each design this reports:
  * the designed-window sequence and its substitutions against wild-type
  * window identity, and whether it clears abguard's novelty floor
  * the refold confidence metrics, if a scored run is supplied

Usage:
    analyze_designs.py --designs runs/<run>/designs \
                       --reference inputs/12-her2/her2_fv_1n8z.cif \
                       [--scored runs/<run>/scored/refold_metrics.csv] \
                       [--chain B] [--window 99 108]
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# abguard rejects a design whose CDRH3 is too close to trastuzumab's or to any
# reference CDRH3. Reported by `abguard --help`:
#   trastuzumab variant (>=70% identity to a known antibody) -> >= 3 substitutions
#   de novo (<70% identity)                                  -> >= 2 substitutions
NOVELTY_FLOOR_OPTIMISATION = 3
NOVELTY_FLOOR_DE_NOVO = 2

# WARNING. The distance this script reports is to WILD-TYPE ONLY, and that is a
# necessary-but-not-sufficient condition. The real gate measures distance to the
# nearest of 525,063 reference CDRH3s, and the 524k-variant FACS table samples
# the near-wild-type neighbourhood densely, so a design can sit >= 3 from
# trastuzumab and still be 2 from some measured variant.
#
# This is not hypothetical: her2_cdrh3_2 (WGPPYFYAMD) is 3 substitutions from
# wild-type, which this script called novel, and abguard REJECTED it -- it is 2
# mutations from SRWGMPHFYAMDY in tz_her2_facs_3classes_524346.
#
# Treat the column below as a cheap pre-filter only. The verdict comes from:
#   abguard team-validate <submission>.csv --reference-db her2_novelty.duckdb
# See TASKS.md, "The novelty gate is much harder than arithmetic suggests".


def parse_cif_chains(path: Path) -> dict[str, dict[int, str]]:
    """Return {chain_id: {seq_id: one_letter}} from an mmCIF's CA atoms.

    Reads the _atom_site loop by column name rather than position, so it
    tolerates the differing column orders that Boltz-2, RFD4 and RCSB emit.
    """
    lines = path.read_text().splitlines()
    header: list[str] = []
    rows: list[list[str]] = []

    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        j = i + 1
        block: list[str] = []
        while j < len(lines) and lines[j].strip().startswith("_"):
            block.append(lines[j].strip())
            j += 1
        if any(h.startswith("_atom_site.") for h in block):
            header = block
            while j < len(lines):
                stripped = lines[j].strip()
                if not stripped or stripped.startswith(("loop_", "_", "#")):
                    break
                rows.append(lines[j].split())
                j += 1
            break
        i = j

    if not header:
        raise ValueError(f"no _atom_site loop found in {path}")

    idx = {name: n for n, name in enumerate(header)}

    def field(row: list[str], *names: str) -> str | None:
        for name in names:
            if name in idx and idx[name] < len(row):
                return row[idx[name]]
        return None

    chains: dict[str, dict[int, str]] = {}
    for row in rows:
        if field(row, "_atom_site.label_atom_id") != "CA":
            continue
        chain = field(row, "_atom_site.auth_asym_id", "_atom_site.label_asym_id")
        seq_id = field(row, "_atom_site.label_seq_id", "_atom_site.auth_seq_id")
        comp = field(row, "_atom_site.label_comp_id")
        if chain is None or seq_id is None or not seq_id.lstrip("-").isdigit():
            continue
        chains.setdefault(chain, {})[int(seq_id)] = AA3TO1.get(comp or "", "X")
    return chains


def window_sequence(residues: dict[int, str], start: int, end: int) -> str:
    return "".join(residues.get(pos, "-") for pos in range(start, end + 1))


def substitutions(wild_type: str, variant: str) -> list[tuple[int, str, str]]:
    """Positional substitutions, 1-based within the window."""
    return [
        (n + 1, wt, mut)
        for n, (wt, mut) in enumerate(zip(wild_type, variant))
        if wt != mut
    ]


def load_scores(path: Path) -> dict[str, dict[str, str]]:
    with path.open() as handle:
        return {row["query_id"]: row for row in csv.DictReader(handle)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--designs", type=Path, required=True,
                        help="directory of <design_id>/<design_id>_generated.cif")
    parser.add_argument("--reference", type=Path, required=True,
                        help="wild-type mmCIF the design was conditioned on")
    parser.add_argument("--scored", type=Path,
                        help="refold_metrics.csv from `module-run refold-score`")
    parser.add_argument("--chain", default="B", help="designed chain (default: B)")
    parser.add_argument("--window", type=int, nargs=2, default=(99, 108),
                        metavar=("START", "END"),
                        help="designed residue range, inclusive (default: 99 108)")
    parser.add_argument("--json", action="store_true", help="emit JSON instead of a table")
    args = parser.parse_args()

    start, end = args.window

    reference_chains = parse_cif_chains(args.reference)
    if args.chain not in reference_chains:
        print(f"error: chain {args.chain} not in {args.reference} "
              f"(found: {', '.join(sorted(reference_chains))})", file=sys.stderr)
        return 1
    wt_window = window_sequence(reference_chains[args.chain], start, end)

    scores = load_scores(args.scored) if args.scored else {}

    records = []
    for cif in sorted(args.designs.glob("*/*_generated.cif")):
        design_id = cif.parent.name
        chains = parse_cif_chains(cif)
        if args.chain not in chains:
            continue
        variant = window_sequence(chains[args.chain], start, end)
        subs = substitutions(wt_window, variant)
        identity = (len(wt_window) - len(subs)) / len(wt_window) if wt_window else 0.0
        row = scores.get(design_id, {})

        def number(key: str) -> float | None:
            raw = row.get(key, "")
            try:
                return float(raw)
            except (TypeError, ValueError):
                return None

        records.append({
            "design_id": design_id,
            "window_sequence": variant,
            "n_substitutions": len(subs),
            "substitutions": [f"{wt}{start + pos - 1}{mut}" for pos, wt, mut in subs],
            "window_identity": round(identity, 3),
            "clears_optimisation_floor": len(subs) >= NOVELTY_FLOOR_OPTIMISATION,
            "clears_de_novo_floor": len(subs) >= NOVELTY_FLOOR_DE_NOVO,
            "refolded": bool(row),
            "min_ipae": number("min_ipae"),
            "iptm": number("iptm"),
            "binder_plddt": number("binder_plddt"),
            "backbone_rmsd": number("backbone_rmsd"),
        })

    if not records:
        print(f"error: no *_generated.cif under {args.designs}", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps({"wild_type_window": wt_window,
                          "window": [start, end],
                          "chain": args.chain,
                          "designs": records}, indent=2))
        return 0

    print(f"wild-type {args.chain}{start}-{end}: {wt_window}\n")
    header = f"{'design':<16}{'window':<14}{'subs':>5}{'ident':>8}{'vs_wt':>7}" \
             f"{'min_ipae':>10}{'iptm':>8}{'pLDDT':>8}"
    print(header)
    print("-" * len(header))
    # refolded designs first, then by fewest substitutions (closest to wild-type)
    for rec in sorted(records, key=lambda r: (not r["refolded"], r["n_substitutions"])):
        fmt = lambda v, p=3: f"{v:.{p}f}" if v is not None else "-"
        print(f"{rec['design_id']:<16}{rec['window_sequence']:<14}"
              f"{rec['n_substitutions']:>5}{rec['window_identity']:>8.0%}"
              f"{('yes' if rec['clears_optimisation_floor'] else 'NO'):>7}"
              f"{fmt(rec['min_ipae']):>10}{fmt(rec['iptm']):>8}"
              f"{fmt(rec['binder_plddt'], 1):>8}")

    unrefolded = sum(1 for r in records if not r["refolded"])
    if unrefolded:
        print(f"\n{unrefolded} design(s) not refolded — metrics shown as '-'.")
    print(
        "\n'vs_wt' = clears the >=3-substitution floor AGAINST WILD-TYPE ONLY."
        "\n  This is a pre-filter, NOT the novelty verdict. The real gate also"
        "\n  compares against 525,063 reference CDRH3s, and designs that pass"
        "\n  here have been rejected there. Confirm every submission with:"
        "\n    abguard team-validate <submission>.csv --reference-db her2_novelty.duckdb"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
