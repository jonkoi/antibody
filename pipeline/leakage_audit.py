#!/usr/bin/env python3
"""How optimistic is the random-split AUC, and by how much?

Dependencies: numpy, scipy, scikit-learn.

Usage:
    .venv-ml/bin/python pipeline/leakage_audit.py [--sample 20000]

THE PROBLEM
Every one of the 51,417 legal candidates is at least 3 substitutions from all
525,063 library CDRH3s -- that is what the novelty filter enforces.  So at
prediction time the model never sees a sequence with a close relative in its
training set.

A random CV split does not reproduce that.  This script measures what fraction
of held-out rows have a Hamming-<= 2 neighbour in the training fold, shell by
shell, and then -- for the shells where enough rows survive the filter --
reports the AUC twice: on all held-out rows, and on only the rows with no close
training relative.  The difference is the optimism that the random-split number
carries.

It cannot be measured directly at distance 3, where *every* row has a close
training neighbour, so the distance-3 AUC has to be discounted by the gap seen
at distances 6-10 instead.  Stating that honestly is the point of this script.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "16")

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import features as F
from binding_data import auc, auc_ci, kfold_indices, load_facs
from neighbourhood import min_dist_ge3_mask


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", type=int, default=20000,
                   help="max rows per distance shell to audit (the filter is O(16k) per row)")
    ap.add_argument("--C", type=float, default=3.0)
    args = ap.parse_args()

    facs = load_facs()
    codes, y, dist = facs["codes"], facs["y"], facs["dist"]

    from sklearn.linear_model import LogisticRegression
    tr, te = next(kfold_indices(y.size, 5, seed=0))
    x = F.mutation_sparse_fast(codes, 2)
    print(f"fitting the pairwise model on {tr.size:,} rows ...", flush=True)
    model = LogisticRegression(C=args.C, max_iter=300).fit(x[tr], y[tr])
    score = model.decision_function(x[te])
    # te is a random permutation slice (binding_data.kfold_indices), not sorted,
    # so np.searchsorted(te, ...) silently returns garbage positions; build the
    # real inverse map from original row index -> position in te/score instead.
    te_pos = np.full(y.size, -1, dtype=np.int64)
    te_pos[te] = np.arange(te.size)

    rng = np.random.default_rng(0)
    print(f"\nheld-out rows with NO training sequence within Hamming distance 2")
    print(f"('novel shell' = the condition the candidates are actually in)\n")
    print(f"  {'dist':>4} {'held out':>9} {'audited':>8} {'novel':>8} {'novel %':>8} "
          f"{'AUC all':>9} {'AUC novel':>11} {'gap':>7}")
    for d in range(1, 11):
        idx = te[dist[te] == d]
        if idx.size < 10:
            continue
        audit = idx if idx.size <= args.sample else rng.choice(idx, args.sample, replace=False)
        keep = min_dist_ge3_mask(codes[audit], codes[tr])
        n_novel = int(keep.sum())
        sub = te_pos[audit]
        assert (sub >= 0).all()
        a_all = auc(score[sub], y[audit])
        if n_novel >= 50 and y[audit[keep]].min() != y[audit[keep]].max():
            a_nov, lo, hi = auc_ci(score[sub[keep]], y[audit[keep]], n_boot=2000, seed=d)
            gap = f"{a_all - a_nov:+.4f}"
            nov = f"{a_nov:.4f}"
        else:
            nov, gap = "   n/a", "    n/a"
        print(f"  {d:>4} {idx.size:>9,} {audit.size:>8,} {n_novel:>8,} "
              f"{100*keep.mean():>7.2f}% {a_all:>9.4f} {nov:>11} {gap:>7}", flush=True)

    print("\nRead this as: wherever the novel shell is big enough to measure, the")
    print("gap is the amount by which a random split flatters the model.  The")
    print("distance-3 row has no measurable novel shell, so its AUC must be")
    print("discounted by roughly the gap seen in the shells that do.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
