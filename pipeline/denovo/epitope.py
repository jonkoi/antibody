#!/usr/bin/env python3
"""Measure which HER2 residues a designed VHH contacts, in a BoltzGen complex.

Reference epitope = HER2 residues within 5 A of the trastuzumab Fab in 1N8Z.
Usage: epitope.py <structure.cif|pdb> [target_chains] [binder_chain]
"""
import math
import sys

# Measured trastuzumab epitope, HER2 mature numbering (from 1N8Z, <5 A)
TRAST_EPI = [557, 558, 560, 561, 569, 570, 571, 572, 573, 579,
             591, 593, 596, 598, 600, 602, 603, 604, 605]
CORE4 = [572, 573, 591, 593]  # the pipeline's hotspot atoms' residues

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M"}


def read_atoms(path):
    """-> list of (chain, resi, resn, atom, x, y, z)"""
    txt = open(path).read()
    out = []
    if "_atom_site." in txt:
        lines = txt.splitlines()
        hdr, start = [], None
        for i, l in enumerate(lines):
            if l.startswith("_atom_site."):
                hdr.append(l.strip().split(".", 1)[1])
                start = i + 1
        idx = {h: k for k, h in enumerate(hdr)}
        ch = "auth_asym_id" if "auth_asym_id" in idx else "label_asym_id"
        ri = "auth_seq_id" if "auth_seq_id" in idx else "label_seq_id"
        for l in lines[start:]:
            if not (l.startswith("ATOM") or l.startswith("HETATM")):
                continue
            p = l.split()
            try:
                out.append((p[idx[ch]], int(p[idx[ri]]), p[idx["label_comp_id"]],
                            p[idx["label_atom_id"]], float(p[idx["Cartn_x"]]),
                            float(p[idx["Cartn_y"]]), float(p[idx["Cartn_z"]])))
            except (ValueError, KeyError, IndexError):
                continue
    else:
        for l in txt.splitlines():
            if l.startswith(("ATOM", "HETATM")):
                out.append((l[21], int(l[22:26]), l[17:20].strip(), l[12:16].strip(),
                            float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return out


def chain_seq(atoms, chain):
    seen, seq = {}, []
    for a in atoms:
        if a[0] == chain and a[1] not in seen:
            seen[a[1]] = 1
            seq.append((a[1], AA3.get(a[2], "X")))
    seq.sort()
    return "".join(s for _, s in seq)


def contacts(atoms, tgt_chains, bnd_chain, cutoff=5.0):
    tg = [a for a in atoms if a[0] in tgt_chains]
    bd = [a for a in atoms if a[0] == bnd_chain]
    best = {}
    for a in tg:
        k = (a[0], a[1], a[2])
        for b in bd:
            d = math.dist(a[4:7], b[4:7])
            if d < best.get(k, 99):
                best[k] = d
    return {k: v for k, v in best.items() if v < cutoff}


def to_her2(chain, resi, mode):
    """Map a chain/resi in the design complex back to HER2 mature numbering."""
    if mode == "cont":      # single chain, pos 1 == HER2 489
        return resi + 488
    if chain == "C" or mode == "C":
        return resi + 488
    return resi + 590       # chain D, pos 1 == HER2 591


def score(path, tgt_chains, bnd_chain, mode):
    atoms = read_atoms(path)
    con = contacts(atoms, tgt_chains, bnd_chain)
    hit = set()
    for (ch, ri, rn), d in con.items():
        hit.add(to_her2(ch, ri, mode if mode == "cont" else ch))
    ep = set(TRAST_EPI)
    inter = hit & ep
    return {
        "binder_seq": chain_seq(atoms, bnd_chain),
        "n_contact": len(hit),
        "epitope_hits": sorted(inter),
        "epi_recall": round(len(inter) / len(ep), 3),          # of trast epitope covered
        "epi_precision": round(len(inter) / max(len(hit), 1), 3),  # of contacts on epitope
        "core4": sorted(hit & set(CORE4)),
        "off_epitope": sorted(hit - ep),
    }


if __name__ == "__main__":
    import json
    p = sys.argv[1]
    tc = list(sys.argv[2]) if len(sys.argv) > 2 else ["A"]
    bc = sys.argv[3] if len(sys.argv) > 3 else "H"
    mode = sys.argv[4] if len(sys.argv) > 4 else "cont"
    print(json.dumps(score(p, tc, bc, mode), indent=1))
