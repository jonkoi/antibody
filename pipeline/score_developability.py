#!/usr/bin/env python3
"""Developability and stability-risk scoring for CDRH3 designs.

No thermostability predictor exists in the hackathon toolchain: Nebius `stab-ddg`
computes *binding* ddG, FoldX and CamSol are licensed and not available, and no
dataset in the challenge carries a melting temperature. So this does not predict Tm.

What it does instead is score the two things that are mechanistically tied to
thermostability and developability, and that are decidable from sequence:

  1. Chemical liabilities — motifs that degrade on storage (deamidation,
     isomerization, oxidation, free cysteine, N-glycosylation sequons). These are
     the LAP Antibody Liability Reference classes that abguard also reports.
  2. Biophysical character — net charge, hydrophobicity, aromatic load and
     backbone rigidity, which drive aggregation and unfolded-state entropy.

Everything is reported as a delta against wild-type, because the competition scores
improvement over trastuzumab rather than absolute values.

Motifs are evaluated across the junction with flanking framework residues, not just
inside the designed window: a liability can be created at the boundary by a
substitution at either end.

Usage:
    score_developability.py --designs runs/<run>/designs \
                            --reference inputs/12-her2/her2_fv_1n8z.cif
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_designs import parse_cif_chains, window_sequence  # noqa: E402

# --- liability motifs ------------------------------------------------------
# (regex, label, penalty, mechanism). Penalties are ordinal risk weights, not
# energies: they rank designs against each other and mean nothing in isolation.
LIABILITIES: list[tuple[str, str, float, str]] = [
    (r"N[GS]",      "deamidation-high",   3.0, "Asn-Gly/Ser: fast deamidation to iso-Asp"),
    (r"N[TNH]",     "deamidation-low",    1.0, "Asn-Thr/Asn/His: slower deamidation"),
    (r"DG",         "isomerization-high", 3.0, "Asp-Gly: fast succinimide isomerization"),
    (r"D[SDP]",     "isomerization-low",  1.5, "Asp-Ser/Asp/Pro: slower isomerization"),
    (r"M",          "oxidation-met",      1.5, "Met: oxidation to sulfoxide"),
    (r"W",          "oxidation-trp",      1.0, "Trp: photo-oxidation"),
    (r"C",          "free-cysteine",     10.0, "Unpaired Cys: hard filter failure"),
    (r"N[^P][ST]",  "n-glycosylation",    4.0, "N-X-S/T sequon: heterogeneous glycan"),
]

# Kyte-Doolittle hydropathy
KD_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
CHARGE = {"K": 1.0, "R": 1.0, "H": 0.1, "D": -1.0, "E": -1.0}
AROMATIC = set("FWY")


def find_liabilities(sequence: str, offset: int) -> list[dict]:
    """Liability motifs in `sequence`, positions reported with `offset` added."""
    found = []
    for pattern, label, penalty, mechanism in LIABILITIES:
        for match in re.finditer(f"(?={pattern})", sequence):
            start = match.start()
            motif = re.match(pattern, sequence[start:])
            if motif is None:
                continue
            found.append({
                "motif": motif.group(0),
                "label": label,
                "position": start + offset,
                "penalty": penalty,
                "mechanism": mechanism,
            })
    return found


def biophysics(window: str) -> dict:
    n = len(window) or 1
    net_charge = sum(CHARGE.get(aa, 0.0) for aa in window)
    return {
        "net_charge": round(net_charge, 2),
        "mean_hydropathy": round(sum(KD_HYDROPATHY.get(aa, 0.0) for aa in window) / n, 3),
        "aromatic_fraction": round(sum(aa in AROMATIC for aa in window) / n, 3),
        "proline_count": window.count("P"),
        "glycine_count": window.count("G"),
    }


def assess(window: str, flank_left: str, flank_right: str) -> dict:
    context = flank_left + window + flank_right
    # only report liabilities that touch the designed window
    window_start, window_end = len(flank_left), len(flank_left) + len(window)
    hits = [
        lia for lia in find_liabilities(context, offset=0)
        if window_start <= lia["position"] + len(lia["motif"]) - 1 and lia["position"] < window_end
    ]
    for lia in hits:
        lia["in_flank"] = not (window_start <= lia["position"] < window_end)
        lia["position"] -= window_start  # 0-based within window, may be negative
    return {
        "liabilities": hits,
        "liability_penalty": round(sum(lia["penalty"] for lia in hits), 2),
        **biophysics(window),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--designs", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--chain", default="B")
    parser.add_argument("--window", type=int, nargs=2, default=(99, 108),
                        metavar=("START", "END"))
    parser.add_argument("--flank", type=int, default=2,
                        help="framework residues of context each side (default: 2)")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    start, end = args.window
    reference = parse_cif_chains(args.reference)[args.chain]
    wt_window = window_sequence(reference, start, end)
    flank_left = window_sequence(reference, start - args.flank, start - 1)
    flank_right = window_sequence(reference, end + 1, end + args.flank)

    wt = assess(wt_window, flank_left, flank_right)

    records = []
    for cif in sorted(args.designs.glob("*/*_generated.cif")):
        chains = parse_cif_chains(cif)
        if args.chain not in chains:
            continue
        window = window_sequence(chains[args.chain], start, end)
        scored = assess(window, flank_left, flank_right)
        removed = {l["label"] for l in wt["liabilities"]} - {l["label"] for l in scored["liabilities"]}
        introduced = {l["label"] for l in scored["liabilities"]} - {l["label"] for l in wt["liabilities"]}
        records.append({
            "design_id": cif.parent.name,
            "window": window,
            **scored,
            "liability_delta": round(scored["liability_penalty"] - wt["liability_penalty"], 2),
            "removed_liabilities": sorted(removed),
            "introduced_liabilities": sorted(introduced),
            "rigidity_delta": (scored["proline_count"] - wt["proline_count"])
                              - (scored["glycine_count"] - wt["glycine_count"]),
        })

    if not records:
        print(f"error: no *_generated.cif under {args.designs}", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps({"wild_type": {"window": wt_window, **wt}, "designs": records}, indent=2))
        return 0

    print(f"wild-type {args.chain}{start}-{end}: {wt_window}"
          f"   (context {flank_left}|{wt_window}|{flank_right})")
    print(f"wild-type liabilities: "
          f"{', '.join(f'{l['motif']}@{l['position']} {l['label']}' for l in wt['liabilities']) or 'none'}"
          f"  → penalty {wt['liability_penalty']}")
    print(f"wild-type biophysics: charge {wt['net_charge']:+.1f}, "
          f"hydropathy {wt['mean_hydropathy']}, aromatic {wt['aromatic_fraction']}, "
          f"{wt['proline_count']}P/{wt['glycine_count']}G\n")

    header = (f"{'design':<16}{'window':<13}{'liab':>6}{'Δliab':>7}"
              f"{'charge':>8}{'hydro':>8}{'rigid':>7}  notes")
    print(header)
    print("-" * (len(header) + 20))
    for rec in sorted(records, key=lambda r: (r["liability_delta"], -r["rigidity_delta"])):
        notes = []
        if rec["removed_liabilities"]:
            notes.append("removes " + "+".join(rec["removed_liabilities"]))
        if rec["introduced_liabilities"]:
            notes.append("ADDS " + "+".join(rec["introduced_liabilities"]))
        print(f"{rec['design_id']:<16}{rec['window']:<13}"
              f"{rec['liability_penalty']:>6.1f}{rec['liability_delta']:>+7.1f}"
              f"{rec['net_charge']:>+8.1f}{rec['mean_hydropathy']:>8.2f}"
              f"{rec['rigidity_delta']:>+7d}  {'; '.join(notes)}")

    print("\nΔliab < 0 is an improvement on wild-type. rigidity = ΔPro − ΔGly;")
    print("positive suggests a more rigid loop (lower unfolded-state entropy).")
    print("Ordinal risk weights, not energies: they rank designs, they do not predict Tm.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
