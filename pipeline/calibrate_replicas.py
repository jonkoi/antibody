#!/usr/bin/env python3
"""Does replica CONSISTENCY discriminate binders better than a single ipTM? (Lane B, part 2)

pipeline/calibrate_confidence.py showed that one ipTM per pair barely separates
real nanobody-antigen pairs from shuffled ones. This script asks whether the
spread across repeated predictions of the SAME pair carries signal that a single
score throws away.

Data: chai1_all_per_replica_confidence_metrics.zip (Smorodina et al.,
doi:10.11582/2026.owjnnp5q). One .npz per replica, named
    system_{nb_index}_{ag_index}_{NB_PDB}_{AG_PDB}_trunk_{t}_scores.model_idx_{m}.npz
holding scalars: aggregate_score, ptm, iptm, per_chain_ptm, per_chain_pair_iptm,
has_inter_chain_clashes, chain_chain_clashes.

Labels, as in part 1: NB_PDB == AG_PDB is a cognate (real) pair; anything else is
a shuffled pair, presumed non-binding by the study's design rather than measured.

Features compared, per pair, over its replicas:
    single    one replica's ipTM (what we get from a single fold)
    mean      mean ipTM
    max       best-of-N ipTM
    min       worst-of-N ipTM
    std       standard deviation of ipTM  (consistency itself)
    mean_m_sd mean - std                  (confident AND consistent)
    agg_mean  mean aggregate_score
    clash     fraction of replicas with an inter-chain clash

numpy is not installed on the dev box, so .npy scalars are parsed from the
format header directly. Stdlib only.

Usage:
    python3 pipeline/calibrate_replicas.py
    python3 pipeline/calibrate_replicas.py --limit 20000   # quick subsample
"""

from __future__ import annotations

import argparse
import ast
import io
import re
import struct
import zipfile
from collections import defaultdict
from pathlib import Path

DEFAULT_ZIP = Path("hackathon-data/sigma2-ab-ag-champloo/chai1_all_per_replica_confidence_metrics.zip")
NAME_RE = re.compile(r"^(system_\d+_\d+_([0-9A-Za-z]+)_([0-9A-Za-z]+))_trunk_(\d+)_scores\.model_idx_(\d+)\.npz$")

WANT = ("iptm", "aggregate_score", "has_inter_chain_clashes")

_DTYPES = {  # numpy descr -> (struct code, size)
    "<f4": ("<f", 4), "<f8": ("<d", 8), "<i4": ("<i", 4), "<i8": ("<q", 8),
    "|b1": ("<b", 1), "<u1": ("<B", 1), "|u1": ("<B", 1),
}


def read_npy_scalar(raw: bytes) -> float | None:
    """First element of a .npy array, without numpy."""
    if raw[:6] != b"\x93NUMPY":
        return None
    major = raw[6]
    if major == 1:
        (hlen,) = struct.unpack("<H", raw[8:10])
        start = 10
    else:
        (hlen,) = struct.unpack("<I", raw[8:12])
        start = 12
    header = ast.literal_eval(raw[start : start + hlen].decode("latin1").strip())
    descr = header["descr"]
    if descr not in _DTYPES:
        return None
    code, size = _DTYPES[descr]
    body = raw[start + hlen :]
    if len(body) < size:
        return None
    return float(struct.unpack(code, body[:size])[0])


def mean(xs: list[float]) -> float:
    return sum(xs) / len(xs)


def stdev(xs: list[float]) -> float:
    if len(xs) < 2:
        return 0.0
    m = mean(xs)
    return (sum((x - m) ** 2 for x in xs) / (len(xs) - 1)) ** 0.5


def auc(pos: list[float], neg: list[float]) -> float:
    """ROC AUC via rank sum, ties half-weighted."""
    if not pos or not neg:
        return float("nan")
    merged = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    rank_sum = 0.0
    i = 0
    while i < len(merged):
        j = i
        while j + 1 < len(merged) and merged[j + 1][0] == merged[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        rank_sum += avg * sum(1 for k in range(i, j + 1) if merged[k][1] == 1)
        i = j + 1
    n_p, n_n = len(pos), len(neg)
    return (rank_sum - n_p * (n_p + 1) / 2.0) / (n_p * n_n)


def auc_ci(a: float, n_pos: int, n_neg: int) -> tuple[float, float, float]:
    """Hanley-McNeil SE and 95% CI. With only ~150 cognate pairs, an AUC gap of a
    few hundredths is not a real difference between features."""
    if not (0 < a < 1) or n_pos < 2 or n_neg < 2:
        return float("nan"), float("nan"), float("nan")
    q1 = a / (2 - a)
    q2 = 2 * a * a / (1 + a)
    var = (a * (1 - a) + (n_pos - 1) * (q1 - a * a) + (n_neg - 1) * (q2 - a * a)) / (n_pos * n_neg)
    se = var**0.5
    return se, a - 1.96 * se, a + 1.96 * se


def quantile(values: list[float], q: float) -> float:
    s = sorted(values)
    pos = q * (len(s) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (pos - lo)


def fpr_at(neg: list[float], t: float) -> float:
    return sum(1 for v in neg if v >= t) / len(neg)


def tpr_at(pos: list[float], t: float) -> float:
    return sum(1 for v in pos if v >= t) / len(pos)


def threshold_for_fpr(neg: list[float], target: float) -> float:
    ordered = sorted(neg, reverse=True)
    allowed = int(len(ordered) * target)
    if allowed >= len(ordered):
        return float("-inf")
    return ordered[allowed] + 1e-12


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    ap.add_argument("--limit", type=int, default=0, help="read only the first N replica files")
    ap.add_argument("--min-replicas", type=int, default=10, help="pairs with fewer replicas are dropped")
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help="write the per-pair feature table here, so re-analysis does not repeat the ~25 min parse",
    )
    args = ap.parse_args()

    if not args.zip.exists():
        raise SystemExit(f"not found: {args.zip}")

    per_system: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    cognate_flag: dict[str, bool] = {}
    skipped = 0
    read = 0

    with zipfile.ZipFile(args.zip) as z:
        names = z.namelist()
        if args.limit:
            names = names[: args.limit]
        print(f"source: {args.zip}")
        print(f"replica files: {len(names)}")
        for name in names:
            m = NAME_RE.match(name)
            if not m:
                skipped += 1
                continue
            system, nb_pdb, ag_pdb = m.group(1), m.group(2), m.group(3)
            cognate_flag[system] = nb_pdb.upper() == ag_pdb.upper()
            try:
                inner = zipfile.ZipFile(io.BytesIO(z.read(name)))
            except Exception:
                skipped += 1
                continue
            for key in WANT:
                try:
                    val = read_npy_scalar(inner.read(f"{key}.npy"))
                except KeyError:
                    continue
                if val is not None:
                    per_system[system][key].append(val)
            read += 1

    print(f"parsed {read} files, skipped {skipped}")

    features: dict[str, tuple[list[float], list[float]]] = defaultdict(lambda: ([], []))
    table: list[dict] = []
    n_pairs = 0
    for system, vals in per_system.items():
        iptm = vals.get("iptm", [])
        if len(iptm) < args.min_replicas:
            continue
        n_pairs += 1
        agg = vals.get("aggregate_score", [])
        clash = vals.get("has_inter_chain_clashes", [])
        sd = stdev(iptm)
        row = {
            "single": iptm[0],
            "mean": mean(iptm),
            "max": max(iptm),
            "min": min(iptm),
            "std": sd,
            "mean_m_sd": mean(iptm) - sd,
            "agg_mean": mean(agg) if agg else float("nan"),
            "clash": (sum(clash) / len(clash)) if clash else float("nan"),
        }
        is_cog = cognate_flag[system]
        table.append({"system": system, "cognate": int(is_cog), "n_replicas": len(iptm), **row})
        for key, value in row.items():
            if value != value:  # NaN
                continue
            features[key][0 if is_cog else 1].append(value)

    n_cog = sum(1 for s in per_system if cognate_flag.get(s) and len(per_system[s].get("iptm", [])) >= args.min_replicas)
    print(f"pairs with >= {args.min_replicas} replicas: {n_pairs} ({n_cog} cognate, {n_pairs - n_cog} shuffled)")
    reps = [len(v.get("iptm", [])) for v in per_system.values()]
    if reps:
        print(f"replicas per pair: min={min(reps)} median={quantile([float(r) for r in reps], 0.5):.0f} max={max(reps)}")
    print()

    print("=" * 86)
    print("DISCRIMINATION BY FEATURE  (cognate vs shuffled, Chai-1)")
    print("=" * 86)
    print(f"  {'feature':<10}{'AUC':>7}{'95% CI':>18}  {'cog med':>9}{'shuf med':>10}   note")
    notes = {
        "single": "one fold, what our pipeline currently gets",
        "mean": "average over replicas",
        "max": "best-of-N, the optimistic read",
        "min": "worst-of-N, the pessimistic read",
        "std": "consistency alone (higher = more scattered)",
        "mean_m_sd": "confident AND consistent",
        "agg_mean": "Chai-1's own aggregate score",
        "clash": "share of replicas with an inter-chain clash",
    }
    order = sorted(features, key=lambda k: -abs(auc(*features[k]) - 0.5))
    for key in order:
        pos, neg = features[key]
        if not pos or not neg:
            continue
        a = auc(pos, neg)
        _se, lo, hi = auc_ci(a, len(pos), len(neg))
        flag = "" if (lo > 0.5 or hi < 0.5) else "  <- includes 0.5"
        print(
            f"  {key:<10}{a:>7.3f}{f'[{lo:.3f}, {hi:.3f}]':>18}  "
            f"{quantile(pos, 0.5):>9.3f}{quantile(neg, 0.5):>10.3f}   {notes.get(key, '')}{flag}"
        )
    print()
    print("  AUC 0.5 = coin flip. Below 0.5 means the feature runs the other way")
    print("  (for std and clash, lower is expected to be better, so < 0.5 is the useful direction).")
    print()

    best = max(features, key=lambda k: abs(auc(*features[k]) - 0.5))
    pos, neg = features[best]
    single_pos, single_neg = features["single"]
    print("=" * 86)
    print(f"BEST FEATURE '{best}' vs SINGLE ipTM, at matched false-positive rates")
    print("=" * 86)
    print(f"  {'target FPR':<12}{'single ipTM cut (kept)':>30}{f'{best} cut (kept)':>30}")
    for target in (0.20, 0.10, 0.05, 0.01):
        t_single = threshold_for_fpr(single_neg, target)
        t_best = threshold_for_fpr(neg, target)
        print(
            f"  {target:<12.0%}"
            f"{f'>={t_single:.3f} ({tpr_at(single_pos, t_single):.0%})':>30}"
            f"{f'>={t_best:.3f} ({tpr_at(pos, t_best):.0%})':>30}"
        )
    print()
    print("Caveats:")
    print("  - Negatives are shuffled pairs, not SPR-confirmed non-binders.")
    print("  - Chai-1 only: the Boltz-2 per-replica archive is 1.7 TB and the AlphaFold3")
    print("    shuffled archive 163 GB, so neither is in the bucket. Whether the same")
    print("    ordering holds for Boltz-2, which our pipeline folds with, is UNVERIFIED.")
    print("  - Cognate pairs are few, so AUCs carry real uncertainty; the 95% CI above is")
    print("    the honest read. A CI spanning 0.5 means the feature is not distinguishable")
    print("    from a coin flip on this data.")

    if args.out:
        import csv as _csv

        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", newline="") as fh:
            writer = _csv.DictWriter(fh, fieldnames=list(table[0].keys()))
            writer.writeheader()
            writer.writerows(sorted(table, key=lambda r: (-r["cognate"], r["system"])))
        print(f"\nper-pair features -> {args.out} ({len(table)} rows)")


if __name__ == "__main__":
    main()
