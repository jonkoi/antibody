#!/usr/bin/env python3
"""Predict log10(KD) for a CDRH3 variant from the 180-sequence SPR panel.

Dependencies: numpy, scipy, scikit-learn, xgboost (see .venv-ml).
Optional: ESM-2 embeddings produced by ``esm_embed.py`` (``--esm <dir>``).

Usage:
    .venv-ml/bin/python pipeline/train_affinity_model.py
    .venv-ml/bin/python pipeline/train_affinity_model.py --esm data/esm

WHY THIS IS THE MODEL THAT MATTERS, AND WHY IT IS HARD
The competition threshold is trastuzumab's median KD, 0.146 nM.  FACS high /
medium / low is a binder call at a single display concentration, not an
affinity, so no amount of FACS data answers the question directly.  Only the
SPR panel does, and it has 180 unique 10-residue CDRH3s with a fitted KD.

Three facts bound what can honestly be claimed, all measured here:

1. Trastuzumab's own eight SPR replicates span 0.029-0.178 nM, i.e.
   sd(log10 KD) = 0.31, a factor of two.  A model cannot be meaningfully better
   than that, and six of the seven variants that "beat 0.146 nM" sit inside
   wild-type's own replicate spread.
2. Only 27 of the 180 SPR sequences appear in the 524k FACS screen, and on that
   overlap the FACS class ranks KD at Spearman +0.19.  FACS "low" variants
   include a 0.45 nM binder; FACS "high" includes a 174 nM one.  So the
   FACS->SPR transfer that the pretraining story assumes is weak, and this
   script measures it rather than assuming it.
3. Edit distance from wild-type alone ranks log10(KD) at Spearman +0.42.  Any
   model must be compared against that, not against zero.

EVALUATION
Repeated stratified K-fold (default 10 folds x 20 repeats) plus leave-one-out.
Hyperparameters are chosen by an inner CV inside each training fold, so the
reported numbers are not tuned on the test rows.  Reported: Spearman, MAE and
RMSE in log10 units, with the spread across repeats.  Plus the only operational
question that matters -- if we take the model's top k, how many true sub-nM
binders do we get, against the base rate.
"""
from __future__ import annotations

import argparse
import json
import sys
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "16")

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import features as F
from binding_data import (REPO, load_spr, spearman, spr_replicate_noise,
                          WT_KD_MEDIAN_NM)

RESULTS_DIR = REPO / "runs/affinity-models"


# ------------------------------------------------------------------ features

def build_blocks(spr, esm_dir: Path | None):
    """Named feature blocks over the SPR sequences."""
    import facs_scorer

    codes = spr["codes"]
    blocks = {
        "dist": spr["dist"].astype(np.float64).reshape(-1, 1),
        "mut": F.mutation_sparse_fast(codes, 1).toarray().astype(np.float64),
        "phys": F.physchem(codes).astype(np.float64),
        "onehot": F.onehot(codes).astype(np.float64),
    }
    facs = facs_scorer.facs_features(codes, holdout_spr=True, cache_key="spr")
    blocks["facs"] = facs
    for i, name in enumerate(facs_scorer.FEATURE_NAMES):
        blocks[f"facs_{name}"] = facs[:, i:i + 1]

    if esm_dir is not None:
        path = esm_dir / "out_spr.npz"
        if path.exists():
            d = np.load(path)
            order = {s: i for i, s in enumerate(d["cdrh3"].tolist())}
            idx = np.array([order[s] for s in spr["seqs"]])
            blocks["esm_cdr"] = d["cdr"][idx].astype(np.float64)
            blocks["esm_vh"] = d["vh"][idx].astype(np.float64)
            blocks["esm_mm"] = d["mm"][idx].astype(np.float64).reshape(-1, 1)
            blocks["esm_res"] = d["res"][idx].astype(np.float64).reshape(len(idx), -1)
        else:
            print(f"  (no ESM file at {path}; skipping ESM blocks)", file=sys.stderr)
    return blocks


def assemble(blocks, names):
    missing = [n for n in names if n not in blocks]
    if missing:
        return None
    return np.concatenate([blocks[n] for n in names], axis=1)


# ------------------------------------------------------------------ estimators

def make_estimator(kind: str):
    """(factory, inner-CV grid).  Every factory returns a fresh sklearn pipeline."""
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import Ridge
    from sklearn.kernel_ridge import KernelRidge
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.decomposition import PCA

    if kind == "mean":
        from sklearn.dummy import DummyRegressor
        return (lambda **kw: DummyRegressor(strategy="mean")), [{}]
    if kind == "ridge":
        return ((lambda alpha: make_pipeline(StandardScaler(), Ridge(alpha=alpha))),
                [{"alpha": a} for a in (0.3, 1, 3, 10, 30, 100, 300, 1000, 3000)])
    if kind == "krr":
        return ((lambda alpha, gamma: make_pipeline(
                    StandardScaler(), KernelRidge(kernel="rbf", alpha=alpha, gamma=gamma))),
                [{"alpha": a, "gamma": g}
                 for a in (0.03, 0.1, 0.3, 1, 3, 10)
                 for g in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2)])
    if kind == "pca_ridge":
        return ((lambda n, alpha: make_pipeline(
                    StandardScaler(), PCA(n_components=n, random_state=0),
                    StandardScaler(), Ridge(alpha=alpha))),
                [{"n": n, "alpha": a} for n in (5, 10, 20, 40)
                 for a in (0.3, 1, 3, 10, 30, 100)])
    if kind == "rf":
        return ((lambda n_estimators, min_samples_leaf: RandomForestRegressor(
                    n_estimators=n_estimators, min_samples_leaf=min_samples_leaf,
                    random_state=0, n_jobs=16)),
                [{"n_estimators": 400, "min_samples_leaf": m} for m in (1, 2, 4, 8)])
    if kind == "gbm":
        from xgboost import XGBRegressor
        return ((lambda max_depth, n_estimators, learning_rate: XGBRegressor(
                    max_depth=max_depth, n_estimators=n_estimators,
                    learning_rate=learning_rate, subsample=0.8,
                    colsample_bytree=0.8, reg_lambda=2.0, n_jobs=16,
                    tree_method="hist")),
                [{"max_depth": d, "n_estimators": n, "learning_rate": 0.05}
                 for d in (2, 3, 4) for n in (200, 500)])
    raise ValueError(kind)


def _inner_select(factory, grid, x, y, seed=0, folds=5):
    """Pick the grid point with the lowest inner-CV MSE.  No test rows involved."""
    if len(grid) == 1:
        return grid[0]
    from sklearn.model_selection import KFold
    best, best_mse = grid[0], np.inf
    kf = list(KFold(folds, shuffle=True, random_state=seed).split(x))
    for params in grid:
        err = 0.0
        for tr, te in kf:
            m = factory(**params).fit(x[tr], y[tr])
            err += float(((m.predict(x[te]) - y[te]) ** 2).sum())
        if err < best_mse:
            best, best_mse = params, err
    return best


def repeated_cv(kind, x, y, folds=10, repeats=20, seed=0):
    """Out-of-fold predictions for each repeat; hyperparameters chosen inside."""
    from sklearn.model_selection import StratifiedKFold
    factory, grid = make_estimator(kind)
    strata = np.digitize(y, np.quantile(y, [0.2, 0.4, 0.6, 0.8]))
    preds = np.empty((repeats, y.size))
    for r in range(repeats):
        skf = StratifiedKFold(folds, shuffle=True, random_state=seed + r)
        for tr, te in skf.split(x, strata):
            params = _inner_select(factory, grid, x[tr], y[tr], seed=seed + r)
            m = factory(**params).fit(x[tr], y[tr])
            preds[r, te] = m.predict(x[te])
    return preds


def loo(kind, x, y, seed=0):
    factory, grid = make_estimator(kind)
    out = np.empty(y.size)
    for i in range(y.size):
        tr = np.delete(np.arange(y.size), i)
        params = _inner_select(factory, grid, x[tr], y[tr], seed=seed)
        out[i] = factory(**params).fit(x[tr], y[tr]).predict(x[i:i + 1])[0]
    return out


def metrics(pred, y):
    return {
        "spearman": spearman(pred, y),
        "mae": float(np.mean(np.abs(pred - y))),
        "rmse": float(np.sqrt(np.mean((pred - y) ** 2))),
    }


def enrichment(pred, y, k=10, threshold=0.0):
    """Of the k lowest predicted log10(KD), how many are truly below threshold."""
    hit = y < threshold
    top = np.argsort(pred)[:k]
    base = float(hit.mean())
    got = float(hit[top].mean())
    return {"k": k, "threshold_log10": threshold, "base_rate": base,
            "precision_at_k": got,
            "enrichment": got / base if base > 0 else float("nan")}


# ------------------------------------------------------------------ main

CONFIGS = [
    ("mean", ["dist"], "predict the mean (null)"),
    ("ridge", ["dist"], "edit distance only"),
    ("ridge", ["facs"], "FACS model scores only (4 feat)"),
    ("ridge", ["dist", "facs"], "distance + FACS scores"),
    ("ridge", ["mut"], "single-mutation coding (190)"),
    ("ridge", ["dist", "mut"], "distance + single-mutation coding"),
    ("ridge", ["dist", "facs", "mut"], "distance + FACS + mutation coding"),
    ("ridge", ["dist", "phys"], "distance + physicochemical"),
    ("ridge", ["dist", "facs", "phys", "mut"], "distance + FACS + phys + mutation"),
    ("krr", ["dist", "facs", "mut"], "kernel ridge (RBF): dist + FACS + mutation"),
    ("rf", ["dist", "facs", "phys", "mut"], "random forest, all tabular"),
    ("gbm", ["dist", "facs", "phys", "mut"], "gradient boosting, all tabular"),
    ("ridge", ["esm_mm"], "ESM-2 masked-marginal only"),
    ("ridge", ["dist", "esm_mm"], "distance + ESM-2 masked marginal"),
    ("pca_ridge", ["esm_cdr"], "ESM-2 CDRH3 embedding, PCA+ridge"),
    ("pca_ridge", ["esm_res"], "ESM-2 per-residue embedding, PCA+ridge"),
    ("ridge", ["dist", "facs", "esm_mm"], "distance + FACS + ESM mm"),
    ("pca_ridge", ["dist", "facs", "esm_cdr"], "distance + FACS + ESM embedding"),
    ("krr", ["dist", "facs", "mut", "esm_mm"], "kernel ridge, + ESM mm"),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--esm", type=Path, default=REPO / "data/esm")
    ap.add_argument("--folds", type=int, default=10)
    ap.add_argument("--repeats", type=int, default=20)
    ap.add_argument("--loo", action="store_true", help="also run leave-one-out")
    ap.add_argument("--tag", default="main")
    args = ap.parse_args()

    spr = load_spr()
    y = spr["log10kd"]
    noise = spr_replicate_noise()
    thresh = float(np.log10(WT_KD_MEDIAN_NM))
    sub_nm = float(np.log10(1.0))

    print(f"SPR panel: {y.size} unique 10-residue CDRH3s with a fitted KD")
    print(f"  log10(KD): mean {y.mean():+.2f}  sd {y.std(ddof=1):.2f}  "
          f"range {y.min():+.2f} .. {y.max():+.2f}")
    print(f"  wild-type replicate noise: sd(log10 KD) = {noise:.3f} "
          f"-- the floor no model can beat")
    print(f"  below 0.146 nM (the competition bar): "
          f"{int((y < thresh).sum())}/{y.size} = {(y < thresh).mean():.1%}")
    print(f"  below 1 nM: {int((y < sub_nm).sum())}/{y.size} = {(y < sub_nm).mean():.1%}")
    print(f"  Spearman(edit distance, log10 KD) = {spearman(spr['dist'].astype(float), y):+.3f}")

    blocks = build_blocks(spr, args.esm if args.esm.exists() else None)
    print(f"\nfeature blocks: " + ", ".join(f"{k}({v.shape[1]})" for k, v in blocks.items()))

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    print(f"\n{'='*104}")
    print(f"{'model':<48}{'rho':>14}{'MAE':>13}{'RMSE':>7}{'P@10<1nM':>10}{'enr':>6}")
    print(f"{'':48}{'(mean+-sd over repeats)':>14}")
    print("=" * 104)
    for kind, names, label in CONFIGS:
        x = assemble(blocks, names)
        if x is None:
            continue
        preds = repeated_cv(kind, x, y, args.folds, args.repeats)
        per = [metrics(p, y) for p in preds]
        rho = np.array([m["spearman"] for m in per])
        mae = np.array([m["mae"] for m in per])
        rmse = np.array([m["rmse"] for m in per])
        mean_pred = preds.mean(0)
        enr = enrichment(mean_pred, y, k=10, threshold=sub_nm)
        row = dict(kind=kind, blocks=names, label=label, n_features=int(x.shape[1]),
                   rho=float(rho.mean()), rho_sd=float(rho.std()),
                   mae=float(mae.mean()), mae_sd=float(mae.std()),
                   rmse=float(rmse.mean()), **enr)
        if args.loo:
            row["loo"] = metrics(loo(kind, x, y), y)
        rows.append(row)
        print(f"{label[:47]:<48}{rho.mean():>6.3f}+-{rho.std():.3f}"
              f"{mae.mean():>8.3f}+-{mae.std():.3f}{rmse.mean():>7.3f}"
              f"{enr['precision_at_k']:>10.2f}{enr['enrichment']:>6.1f}")
        np.save(RESULTS_DIR / f"pred_{kind}_{'+'.join(names)}.npy", mean_pred)

    with (RESULTS_DIR / f"summary_{args.tag}.json").open("w") as fh:
        json.dump({"noise_floor_log10": noise, "n": int(y.size), "rows": rows}, fh, indent=2)

    best = max(rows, key=lambda r: r["rho"])
    print(f"\nbest by Spearman: {best['label']}  rho={best['rho']:.3f}+-{best['rho_sd']:.3f}, "
          f"MAE={best['mae']:.3f} log10 units")
    print(f"assay replicate floor is MAE ~{noise * 0.8:.2f} log10 units "
          f"(0.8 x sd for a half-normal); null model MAE = {rows[0]['mae']:.3f}")
    print(f"\nwrote {RESULTS_DIR / f'summary_{args.tag}.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
