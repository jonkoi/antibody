#!/usr/bin/env python3
"""Build StaB-ddG request payloads: a per-position scan, and our designs.

Two jobs, because they answer different questions and carry different weight:

  scan   -- every position of B99-108 mutated to Ala in the WILD-TYPE background,
            plus W99 to each residue our accepted designs actually use there.
            These are genuine point mutations, which is what StaB-ddG models, so
            this is the on-label job and the one that can settle whether W99 is
            load-bearing.

  designs -- our 14 shortlisted designs, each as ONE variant carrying all 5-9 of
            its substitutions. The API accepts multi-site variants, so the model
            scores them jointly rather than us summing per-mutation values -- but
            whether StaB-ddG's training covers multi-mutants is UNVERIFIED, so
            treat the output as ordinal and say so. See STRATEGY.md.

Partner A is the antibody (chains A light, B heavy); partner B is HER2 (chains C
and D, split at the unresolved 581-590 loop). Mutation sites are named by the
ORIGINAL PDB chain and author residue number, and every wild_type below is
checked against the structure before anything is written.

Usage:
    stab_ddg_jobs.py --pdb runs/stab-ddg/wildtype_complex.pdb \
                     --reference inputs/12-her2/her2_fv_1n8z.cif \
                     --shortlist runs/team15-claude-n100-seed1/shortlist.csv \
                     --outdir runs/stab-ddg
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import tarfile
from pathlib import Path

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

HEAVY_CHAIN = "B"
WINDOW = range(99, 109)
PARTNER_A = ["A", "B"]      # antibody: light + heavy
PARTNER_B = ["C", "D"]      # HER2 domain IV, split at the unresolved loop


def wild_type_window(reference: Path) -> dict[int, str]:
    """{residue_number: one_letter} for B99-108, read from the structure."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from cif_to_pdb import read_atom_site, null  # noqa: E402

    out: dict[int, str] = {}
    for r in read_atom_site(reference):
        if null(r.get("label_asym_id")) != HEAVY_CHAIN:
            continue
        try:
            n = int(null(r.get("auth_seq_id")))
        except (TypeError, ValueError):
            continue
        if n in WINDOW:
            aa = AA3TO1.get(null(r.get("label_comp_id")).upper())
            if aa:
                out[n] = aa
    missing = [n for n in WINDOW if n not in out]
    if missing:
        sys.exit(f"reference lacks B{missing} -- cannot build mutations")
    return out


def mutation(resnum: int, wt: str, mut: str) -> dict:
    return {
        "chain": HEAVY_CHAIN,
        "residue_number": resnum,
        "insertion_code": "",
        "wild_type": wt,
        "mutant": mut,
    }


def build_scan(wt: dict[int, str], designs: list[str]) -> list[dict]:
    """Alanine scan of the window, plus W99 to each residue our designs use."""
    variants = []
    for n in WINDOW:
        w = wt[n]
        if w == "A":
            continue  # Ala-to-Ala is not a mutation
        variants.append({"id": f"{w}{n}A", "mutations": [mutation(n, w, "A")]})

    first = wt[min(WINDOW)]           # W99 in wild-type trastuzumab
    seen = {"A"}
    for seq in designs:
        mut = seq[0]
        if mut != first and mut not in seen:
            seen.add(mut)
            variants.append({"id": f"{first}99{mut}", "mutations": [mutation(99, first, mut)]})
    return variants


def build_designs(wt: dict[int, str], rows: list[dict]) -> list[dict]:
    """One multi-site variant per shortlisted design."""
    variants = []
    for r in rows:
        seq = r["_window"]
        muts = [
            mutation(n, wt[n], seq[i])
            for i, n in enumerate(WINDOW)
            if seq[i] != wt[n]
        ]
        if not muts:
            continue
        variants.append({"id": r["sequence_id"], "mutations": muts})
    return variants


def write_job(outdir: Path, name: str, pdb: Path, variants: list[dict]) -> None:
    tar_path = outdir / f"{name}.tar"
    with tarfile.open(tar_path, "w") as tf:
        tf.add(pdb, arcname=pdb.name)      # archive root, no leading ./
    params = {
        "pdb_file": pdb.name,
        "partner_a": PARTNER_A,
        "partner_b": PARTNER_B,
        "variants": variants,
    }
    (outdir / f"{name}.params.json").write_text(json.dumps(params, indent=1))
    print(f"{name}: {len(variants)} variants -> {tar_path}, {name}.params.json")
    for v in variants[:3]:
        sites = ", ".join(f"{m['wild_type']}{m['residue_number']}{m['mutant']}" for m in v["mutations"])
        print(f"    {v['id']}: {sites}")
    if len(variants) > 3:
        print(f"    ... {len(variants) - 3} more")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", type=Path, required=True)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--shortlist", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    if not args.pdb.exists():
        sys.exit(f"{args.pdb} not found -- run cif_to_pdb.py first")
    args.outdir.mkdir(parents=True, exist_ok=True)

    wt = wild_type_window(args.reference)
    wt_seq = "".join(wt[n] for n in WINDOW)
    print(f"wild-type {HEAVY_CHAIN}{min(WINDOW)}-{max(WINDOW)}: {wt_seq}")

    rows = []
    with args.shortlist.open(newline="") as fh:
        for r in csv.DictReader(fh):
            desc = r["design_description"]
            if "CDRH3 " not in desc:
                sys.exit(f"cannot read window from: {desc[:80]}")
            r["_window"] = desc.split("CDRH3 ", 1)[1].split(";")[0].strip()
            if len(r["_window"]) != len(wt_seq):
                sys.exit(f"{r['sequence_id']}: window {r['_window']!r} is not {len(wt_seq)} residues")
            rows.append(r)
    print(f"{len(rows)} shortlisted designs")

    write_job(args.outdir, "scan", args.pdb, build_scan(wt, [r["_window"] for r in rows]))
    write_job(args.outdir, "designs", args.pdb, build_designs(wt, rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
