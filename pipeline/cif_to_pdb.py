#!/usr/bin/env python3
"""Convert an mmCIF to PDB, for tools that will not read mmCIF (e.g. StaB-ddG).

Dependency-free. Keeps chain ids and author residue numbers EXACTLY as the mmCIF
gives them, because StaB-ddG's request names residues by chain and number and its
schema warns that the template does not establish those identifiers -- so a
silent renumbering here would produce a confident answer about the wrong residue.

Refuses to write a file it cannot represent faithfully: PDB has fixed-width
columns, so a chain id longer than one character or a residue number outside
-999..9999 is an error rather than a truncation.

Usage:
    cif_to_pdb.py in.cif out.pdb [--chains A B C D]
    cif_to_pdb.py in.cif out.pdb --verify B 99 TRP
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# PDB-standard element guess from a PDBx atom id, used only when the mmCIF omits
# type_symbol. Two-letter elements first so CL/BR/ZN are not read as carbon.
TWO_LETTER = ("CL", "BR", "ZN", "MG", "MN", "FE", "CA", "NA", "SE")


def read_atom_site(path: Path) -> list[dict[str, str]]:
    """Parse the _atom_site loop. Handles quoted values and the '.'/'?' nulls."""
    header: list[str] = []
    rows: list[dict[str, str]] = []
    in_loop = False
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                header.append(s.split(".", 1)[1])
                in_loop = True
                continue
            if in_loop:
                if not s or s.startswith(("#", "loop_", "_", "data_")):
                    in_loop = False
                    continue
                fields = split_cif_line(s)
                if len(fields) != len(header):
                    # Never skip silently: a dropped row is a missing residue or
                    # a missing chain, and downstream tools will not notice.
                    sys.exit(
                        f"{path}: malformed _atom_site row -- {len(fields)} fields, "
                        f"expected {len(header)}:\n  {s[:160]}"
                    )
                rows.append(dict(zip(header, fields)))
    if not rows:
        sys.exit(f"{path}: no _atom_site rows found")
    return rows


def split_cif_line(line: str) -> list[str]:
    """Split on whitespace, respecting single and double quotes.

    An empty quoted value ('' -- which this file uses for ig_region on the HER2
    chains) MUST yield an empty token, not nothing. Dropping it shortens the row
    by one field, the row then looks malformed against the header, and the whole
    chain disappears -- silently producing a 'complex' with no antigen in it.
    """
    out: list[str] = []
    cur = ""
    quote = ""
    had_quote = False
    for ch in line:
        if quote:
            if ch == quote:
                quote = ""
            else:
                cur += ch
        elif ch in "'\"":
            quote = ch
            had_quote = True
        elif ch.isspace():
            if cur or had_quote:
                out.append(cur)
                cur = ""
                had_quote = False
        else:
            cur += ch
    if cur or had_quote:
        out.append(cur)
    return out


def null(v: str | None) -> str:
    return "" if v in (None, ".", "?") else v


def guess_element(atom_id: str) -> str:
    a = atom_id.strip().upper().lstrip("0123456789")
    for e in TWO_LETTER:
        if a.startswith(e):
            return e
    return a[:1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cif", type=Path)
    ap.add_argument("pdb", type=Path)
    ap.add_argument("--chains", nargs="*", help="keep only these chains, in this order")
    ap.add_argument("--verify", nargs=3, action="append", metavar=("CHAIN", "RESNUM", "RESNAME"),
                    help="assert this chain/number/residue exists; repeatable")
    args = ap.parse_args()

    rows = read_atom_site(args.cif)

    # author numbering is what external tools mean by "residue 99"
    def chain_of(r):
        return null(r.get("auth_asym_id")) or null(r.get("label_asym_id"))

    def resnum_of(r):
        return null(r.get("auth_seq_id")) or null(r.get("label_seq_id"))

    keep = args.chains
    out_lines: list[str] = []
    serial = 0
    seen: dict[tuple[str, str], str] = {}
    skipped_alt = 0

    order = keep if keep else sorted({chain_of(r) for r in rows})
    for chain in order:
        chain_rows = [r for r in rows if chain_of(r) == chain]
        if not chain_rows:
            sys.exit(f"chain {chain!r} not present in {args.cif}")
        if len(chain) != 1:
            sys.exit(f"chain id {chain!r} does not fit PDB's one-character column")
        for r in chain_rows:
            alt = null(r.get("label_alt_id"))
            if alt and alt not in ("A", "1"):
                skipped_alt += 1
                continue
            resnum = resnum_of(r)
            try:
                rn = int(resnum)
            except ValueError:
                sys.exit(f"non-integer residue number {resnum!r} in chain {chain}")
            if not -999 <= rn <= 9999:
                sys.exit(f"residue number {rn} in chain {chain} does not fit PDB columns")

            resname = null(r.get("label_comp_id")) or null(r.get("auth_comp_id"))
            atom = null(r.get("label_atom_id")) or null(r.get("auth_atom_id"))
            icode = null(r.get("pdbx_PDB_ins_code"))
            elem = (null(r.get("type_symbol")) or guess_element(atom)).upper()
            seen[(chain, resnum)] = resname

            serial += 1
            # PDB atom-name convention: one-character elements start in column 14
            name = atom if len(atom) >= 4 else f" {atom:<3}"
            try:
                x, y, z = (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))
            except (KeyError, ValueError):
                sys.exit("missing or non-numeric coordinates")
            occ = null(r.get("occupancy")) or "1.00"
            bfac = null(r.get("B_iso_or_equiv")) or "0.00"
            out_lines.append(
                f"ATOM  {serial:>5} {name:<4}{'':1}{resname:>3} {chain}{rn:>4}{icode or ' ':1}   "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{float(occ):>6.2f}{float(bfac):>6.2f}          {elem:>2}"
            )
        out_lines.append(f"TER   {serial + 1:>5}      {resname:>3} {chain}{rn:>4}")
        serial += 1
    out_lines.append("END")

    for chain, resnum, resname in args.verify or []:
        got = seen.get((chain, resnum))
        if got is None:
            sys.exit(f"VERIFY FAILED: {chain}{resnum} is not in the output")
        if got.upper() != resname.upper():
            sys.exit(f"VERIFY FAILED: {chain}{resnum} is {got}, expected {resname}")
        print(f"verified {chain}{resnum} = {got}")

    args.pdb.write_text("\n".join(out_lines) + "\n")
    print(f"wrote {args.pdb}: {serial} records, chains {','.join(order)}")
    if skipped_alt:
        print(f"  skipped {skipped_alt} alternate-location atoms")
    print(f"  residues per chain: " + ", ".join(
        f"{c}={sum(1 for (ch, _) in seen if ch == c)}" for c in order))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
