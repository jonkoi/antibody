#!/usr/bin/env python3
"""Shared data loading, encoding and split construction for the CDRH3 binding models.

Dependencies: numpy (see ``.venv-ml``; create with
``uv venv .venv-ml --python 3.12 && uv pip install --python .venv-ml/bin/python
numpy pandas scikit-learn scipy xgboost lightgbm torch``).

Two datasets:

FACS  -- 524,346 CDRH3 variants with a coarse 3-class binder call
         (high / medium / low).  Convention from the published splits:
         high == binder, medium+low == non-binder.  This is a *binder/non-binder*
         call, not an affinity.
SPR   -- 768 rows, 187 of which carry a real KD on a 10-residue CDRH3.
         This is the target that actually matters but it is tiny.

SPLITS
Random splits leak: the 10-mers form dense single-mutant neighbourhoods, so a
sequence's Hamming-1 neighbour very often sits on the other side of a random
cut.  ``cluster_labels`` builds single-linkage clusters at Hamming distance 1
(connected components of the Hamming-1 graph) so that whole neighbourhoods move
together.  ``kfold_indices`` can split either way; everything downstream reports
both.

The headline number is the AUC restricted to edit distance 3, because every
legal candidate lives in that shell.  It is computed from pooled out-of-fold
predictions so that all 1,490 distance-3 rows contribute, not just the ~300 that
land in one test split.
"""
from __future__ import annotations

import csv
import io
import os
import zipfile
from pathlib import Path

import numpy as np

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(ALPHABET)}
WT_CDRH3 = "WGGDGFYAMD"
LENGTH = 10
VH_PREFIX = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
             "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSR")
VH_SUFFIX = "YWGQGTLVTVSS"
VL_SEQUENCE = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
               "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")

REPO = Path(__file__).resolve().parent.parent
FACS_ZIP = REPO / "data/her2-trastuzumab/tz_her2_facs_3classes_524346.csv.zip"
SPR_CSV = REPO / "data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv"
CANDIDATES_CSV = REPO / "data/legal_3mutants_ranked.csv"
CACHE = REPO / "data/cache"

# Trastuzumab's own 8 SPR replicates span 0.0285-0.178 nM.  The competition
# threshold is their median.  Anything above the replicate floor is not
# distinguishable from wild-type by this assay.
WT_KD_MEDIAN_NM = 0.146
WT_KD_REPLICATES_NM = (0.152, 0.123, 0.146, 0.165, 0.178, 0.0378713, 0.02851802, 0.07916459)

CLASS_ORDINAL = {"low": 0, "medium": 1, "high": 2}


def vh_sequence(cdrh3: str) -> str:
    return VH_PREFIX + cdrh3 + VH_SUFFIX


def encode_int(seqs) -> np.ndarray:
    """(n, 10) int8 of alphabet indices."""
    out = np.empty((len(seqs), LENGTH), dtype=np.int8)
    for i, s in enumerate(seqs):
        for j, c in enumerate(s):
            out[i, j] = AA_INDEX[c]
    return out


def onehot(codes: np.ndarray) -> np.ndarray:
    """(n, 10) int codes -> (n, 200) float32 one-hot, position-major."""
    n = codes.shape[0]
    out = np.zeros((n, LENGTH * 20), dtype=np.float32)
    cols = codes.astype(np.int32) + (np.arange(LENGTH, dtype=np.int32) * 20)[None, :]
    out[np.arange(n)[:, None], cols] = 1.0
    return out


# ---------------------------------------------------------------- FACS

def load_facs(use_cache: bool = True):
    """-> dict with codes (n,10) int8, y (n,) int8 binder, ordinal (n,) int8, dist (n,) int8."""
    CACHE.mkdir(parents=True, exist_ok=True)
    cache = CACHE / "facs.npz"
    if use_cache and cache.exists():
        d = np.load(cache)
        return {k: d[k] for k in d.files}

    seqs, y, ordinal, dist = [], [], [], []
    with zipfile.ZipFile(FACS_ZIP) as archive:
        name = next(n for n in archive.namelist() if not n.startswith("__MACOSX"))
        with archive.open(name) as handle:
            for row in csv.DictReader(io.TextIOWrapper(handle, "utf-8")):
                s = row["cdrh3"].strip().upper()
                if len(s) != LENGTH:
                    continue
                cls = row["binding_class"].strip().lower()
                seqs.append(s)
                y.append(1 if cls == "high" else 0)
                ordinal.append(CLASS_ORDINAL[cls])
                dist.append(int(row["edit_distance"]))
    out = {
        "codes": encode_int(seqs),
        "y": np.asarray(y, dtype=np.int8),
        "ordinal": np.asarray(ordinal, dtype=np.int8),
        "dist": np.asarray(dist, dtype=np.int8),
    }
    np.savez_compressed(cache, **out)
    return out


def facs_sequences() -> list[str]:
    d = load_facs()
    return decode(d["codes"])


def decode(codes: np.ndarray) -> list[str]:
    arr = np.asarray(codes)
    table = np.array(list(ALPHABET))
    return ["".join(row) for row in table[arr.astype(np.int64)]]


# ---------------------------------------------------------------- SPR

def load_spr():
    """-> dict: seqs (list[str]), codes, log10kd, dist, is_wt.

    Only 10-residue CDRH3s with a reportable KD.  The 16-residue rows are a
    different antibody and are dropped; the 9-residue rows cannot be aligned to
    the 10-position encoding and are dropped too (3 of them carry a KD).
    Duplicate sequences are collapsed to the geometric mean of their KDs, which
    also folds the 8 wild-type replicates into one row.
    """
    by_seq: dict[str, list[float]] = {}
    with SPR_CSV.open() as handle:
        for row in csv.DictReader(handle):
            s = row["cdrh3"].strip().upper()
            kd = row["kd_nM"].strip()
            if len(s) != LENGTH or not kd:
                continue
            if any(c not in AA_INDEX for c in s):
                continue
            by_seq.setdefault(s, []).append(float(kd))

    seqs = sorted(by_seq)
    log10kd = np.array([np.mean(np.log10(by_seq[s])) for s in seqs])
    n_rep = np.array([len(by_seq[s]) for s in seqs])
    dist = np.array([sum(a != b for a, b in zip(s, WT_CDRH3)) for s in seqs], dtype=np.int8)
    return {
        "seqs": seqs,
        "codes": encode_int(seqs),
        "log10kd": log10kd,
        "n_rep": n_rep,
        "dist": dist,
        "is_wt": dist == 0,
    }


def spr_replicate_noise() -> float:
    """Std-dev of log10(KD) across the 8 wild-type replicates: the assay floor."""
    return float(np.std(np.log10(WT_KD_REPLICATES_NM), ddof=1))


# ---------------------------------------------------------------- candidates

def load_candidates():
    with CANDIDATES_CSV.open() as handle:
        rows = list(csv.DictReader(handle))
    seqs = [r["cdrh3"] for r in rows]
    return rows, seqs, encode_int(seqs)


# ---------------------------------------------------------------- splits

def cluster_labels(codes: np.ndarray, use_cache: bool = True) -> np.ndarray:
    """Single-linkage clusters at Hamming distance 1 over 10-mers.

    Connected components of the graph joining any two sequences that differ at
    exactly one position.  Found without an O(n^2) scan: for each position,
    blank that position out and group by the remaining 9 residues -- every pair
    inside such a group is a Hamming-1 neighbour.  Union-find over those groups.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    cache = CACHE / f"clusters_{codes.shape[0]}.npy"
    if use_cache and cache.exists():
        return np.load(cache)

    n = codes.shape[0]
    parent = np.arange(n, dtype=np.int64)

    def find(x: int) -> int:
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    codes64 = codes.astype(np.int64)
    for pos in range(LENGTH):
        keep = [p for p in range(LENGTH) if p != pos]
        key = np.zeros(n, dtype=np.int64)
        for p in keep:
            key = key * 20 + codes64[:, p]
        order = np.argsort(key, kind="stable")
        sorted_key = key[order]
        boundaries = np.flatnonzero(np.diff(sorted_key)) + 1
        for grp in np.split(order, boundaries):
            if grp.size < 2:
                continue
            root = find(int(grp[0]))
            for idx in grp[1:]:
                other = find(int(idx))
                if other != root:
                    parent[other] = root

    labels = np.array([find(i) for i in range(n)], dtype=np.int64)
    _, labels = np.unique(labels, return_inverse=True)
    np.save(cache, labels)
    return labels


def kfold_indices(n: int, k: int = 5, groups: np.ndarray | None = None, seed: int = 0):
    """Yield (train_idx, test_idx).  Grouped if ``groups`` is given.

    Groups are assigned to folds greedily by descending size so the folds stay
    balanced even when one cluster swallows most of the data.
    """
    rng = np.random.default_rng(seed)
    if groups is None:
        order = rng.permutation(n)
        for fold in range(k):
            test = order[fold::k]
            train = np.setdiff1d(order, test, assume_unique=True)
            yield train, test
        return

    uniq, counts = np.unique(groups, return_counts=True)
    order = np.argsort(-counts)
    load = np.zeros(k, dtype=np.int64)
    assign = {}
    for i in order:
        f = int(np.argmin(load))
        assign[uniq[i]] = f
        load[f] += counts[i]
    fold_of = np.array([assign[g] for g in groups])
    all_idx = np.arange(n)
    for fold in range(k):
        test = all_idx[fold_of == fold]
        train = all_idx[fold_of != fold]
        yield train, test


# ---------------------------------------------------------------- metrics

def auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Tie-corrected Mann-Whitney AUC."""
    labels = np.asarray(labels).astype(bool)
    n_pos = int(labels.sum())
    n_neg = labels.size - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="stable")
    ranked = scores[order]
    ranks = np.empty(scores.size, dtype=np.float64)
    i = 0
    while i < ranked.size:
        j = i
        while j + 1 < ranked.size and ranked[j + 1] == ranked[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return float((ranks[labels].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def auc_ci(scores: np.ndarray, labels: np.ndarray, n_boot: int = 2000, seed: int = 0):
    """(auc, lo, hi) with a 95% stratified-bootstrap interval.

    The replicate count is scaled down for large n: each replicate costs an
    O(n log n) sort, and once n is in the hundreds of thousands the interval is
    narrower than the third decimal anyway.  The strata we actually care about
    (distance 3, n=1,490) keep the full replicate count.
    """
    point = auc(scores, labels)
    labels = np.asarray(labels).astype(bool)
    pos = np.flatnonzero(labels)
    neg = np.flatnonzero(~labels)
    if pos.size == 0 or neg.size == 0:
        return point, float("nan"), float("nan")
    if scores.size > 20_000:
        n_boot = max(100, min(n_boot, 20_000_000 // scores.size))
    rng = np.random.default_rng(seed)
    vals = np.empty(n_boot)
    for b in range(n_boot):
        p = rng.choice(pos, pos.size, replace=True)
        q = rng.choice(neg, neg.size, replace=True)
        idx = np.concatenate([p, q])
        vals[b] = auc(scores[idx], labels[idx])
    return point, float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def average_precision(scores: np.ndarray, labels: np.ndarray) -> float:
    labels = np.asarray(labels).astype(int)
    order = np.argsort(-scores, kind="stable")
    lab = labels[order]
    tp = np.cumsum(lab)
    precision = tp / np.arange(1, lab.size + 1)
    denom = lab.sum()
    return float((precision * lab).sum() / denom) if denom else float("nan")


def ece(probs: np.ndarray, labels: np.ndarray, n_bins: int = 10):
    """Expected calibration error plus the per-bin table."""
    labels = np.asarray(labels).astype(float)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    idx = np.clip(np.digitize(probs, edges[1:-1]), 0, n_bins - 1)
    total = 0.0
    table = []
    for b in range(n_bins):
        m = idx == b
        if not m.any():
            continue
        conf = float(probs[m].mean())
        frac = float(labels[m].mean())
        total += m.sum() / probs.size * abs(conf - frac)
        table.append((edges[b], edges[b + 1], int(m.sum()), conf, frac))
    return float(total), table


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    def rank(x):
        order = np.argsort(x, kind="stable")
        r = np.empty(x.size, dtype=np.float64)
        sx = x[order]
        i = 0
        while i < sx.size:
            j = i
            while j + 1 < sx.size and sx[j + 1] == sx[i]:
                j += 1
            r[order[i:j + 1]] = (i + j) / 2.0 + 1.0
            i = j + 1
        return r
    ra, rb = rank(np.asarray(a, float)), rank(np.asarray(b, float))
    ra -= ra.mean()
    rb -= rb.mean()
    denom = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / denom) if denom else float("nan")


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "16")
    facs = load_facs()
    print(f"FACS: {facs['codes'].shape[0]:,} rows, "
          f"{facs['y'].mean():.3%} binder")
    print("  by edit distance:")
    for d in range(1, 11):
        m = facs["dist"] == d
        print(f"    {d:>2}: n={m.sum():>7,}  binder={facs['y'][m].mean():.1%}")
    labels = cluster_labels(facs["codes"])
    sizes = np.bincount(labels)
    print(f"\nHamming-1 single-linkage clusters: {sizes.size:,}; "
          f"largest {sizes.max():,} ({sizes.max()/labels.size:.1%} of data); "
          f"singletons {(sizes == 1).sum():,}")

    spr = load_spr()
    print(f"\nSPR: {len(spr['seqs'])} unique 10-mers with KD "
          f"({int(spr['n_rep'].sum())} measurements)")
    print(f"  log10(KD) range {spr['log10kd'].min():.2f} .. {spr['log10kd'].max():.2f}")
    print(f"  wild-type replicate noise: sd(log10 KD) = {spr_replicate_noise():.3f} "
          f"({10 ** spr_replicate_noise():.1f}x)")
    print("  by edit distance:", np.bincount(spr["dist"], minlength=11).tolist())
