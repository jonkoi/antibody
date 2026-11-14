#!/usr/bin/env python3
"""Model comparison for ranking CDRH3 variants by HER2 binding (FACS target).

Dependencies (create the venv once):
    uv venv .venv-ml --python 3.12
    uv pip install --python .venv-ml/bin/python numpy pandas scikit-learn scipy xgboost
    uv pip install --python .venv-ml/bin/python --index-url \
        https://download.pytorch.org/whl/cpu torch

Usage:
    .venv-ml/bin/python pipeline/train_binding_model.py --models all
    .venv-ml/bin/python pipeline/train_binding_model.py --models logodds,pair_lr --folds 5

WHAT IS BEING MEASURED
The headline is the AUC restricted to edit distance 3, because every legal
candidate is a 3-mutant.  It is computed from pooled out-of-fold predictions so
that all 1,490 distance-3 rows contribute (a single 20% split leaves only ~300).
A 95% stratified-bootstrap interval is reported with it, because 1,490 rows at
90% positive leaves only ~149 negatives and the interval is wide.

Three evaluation regimes, in increasing strictness:

  cv        5-fold random CV, pooled out-of-fold.  Directly comparable to the
            0.717 that the incumbent log-odds ranker reports at distance 3.
  shell     train on edit distance >= 5, test on distance <= 4.  Removes the
            strongest form of leakage (a distance-3 test row having a
            Hamming-1 training neighbour that is also a distance-3 row).
  spr       zero-shot transfer: Spearman of the FACS-trained score against the
            180 real SPR log10(KD) values.  The FACS score is a binder call;
            this asks whether it carries affinity information at all.  Sequences
            with an SPR measurement are dropped from FACS training first.

A note on what could not be done.  Deployment conditions cannot be reproduced
inside FACS: every legal candidate is >= 3 substitutions from all 525,063
library sequences, but in FACS *every* row at edit distance <= 4 has a
Hamming-<= 2 neighbour somewhere in a random 80% training fold (measured, not
assumed -- see neighbourhood.py).  The "novel shell" test set is therefore empty
at distance 3.  The optimism this leaves in the distance-3 CV number is
quantified at distance 6-10, where the novel shell is large enough to measure,
by ``--novel-gap``.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "16")

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))

import features as F
from binding_data import (LENGTH, REPO, auc, auc_ci, average_precision, ece,
                          kfold_indices, load_facs, load_spr, spearman)

RESULTS_DIR = REPO / "runs/binding-models"


# ------------------------------------------------------------------ models
# Every model exposes fit(codes, y, ordinal) and score(codes) -> higher = more
# likely to be a binder.  ``proba`` is True when score() returns a calibrated
# probability, which is what the calibration table needs.

class LogOdds:
    """The incumbent: position-independent naive-Bayes log-odds.  200 params."""
    name = "logodds"
    label = "log-odds (incumbent floor)"
    proba = False

    def __init__(self, pseudocount: float = 1.0):
        self.pc = pseudocount

    def fit(self, codes, y, ordinal=None):
        y = y.astype(bool)
        w = np.empty((LENGTH, 20))
        for p in range(LENGTH):
            hi = np.bincount(codes[y, p], minlength=20) + self.pc
            lo = np.bincount(codes[~y, p], minlength=20) + self.pc
            w[p] = np.log(hi / hi.sum()) - np.log(lo / lo.sum())
        self.w = w
        return self

    def score(self, codes):
        return self.w[np.arange(LENGTH)[None, :], codes.astype(np.int64)].sum(1)


class SingleLR:
    """L2 logistic regression on the 10x19 WT-centred single-mutation coding."""
    name = "single_lr"
    label = "logistic, singles only (190 feat)"
    proba = True
    order = 1

    # C swept over {0.3, 1, 3, 10} on pooled-OOF distance-3 AUC: 0.878, 0.887,
    # 0.890, 0.891 for the order-2 model.  Flat above 1, so C=3 is taken rather
    # than the arg-max, to avoid tuning on the reported metric.
    def __init__(self, C: float = 3.0, max_iter: int = 300):
        self.C, self.max_iter = C, max_iter

    def _x(self, codes):
        return F.mutation_sparse_fast(codes, self.order)

    def fit(self, codes, y, ordinal=None):
        from sklearn.linear_model import LogisticRegression
        self.m = LogisticRegression(C=self.C, max_iter=self.max_iter,
                                    solver="lbfgs", n_jobs=-1)
        self.m.fit(self._x(codes), y)
        return self

    def score(self, codes):
        return self.m.predict_proba(self._x(codes))[:, 1]


class PairLR(SingleLR):
    """Singles + all 45 position pairs, WT-centred.  16,435 features."""
    name = "pair_lr"
    label = "logistic, singles+pairs (16,435 feat)"
    order = 2


class TripleLR(SingleLR):
    """Singles + pairs + all 120 triples.  839,515 features."""
    name = "triple_lr"
    label = "logistic, +triples (839,515 feat)"
    order = 3


class OneHotPairLR(SingleLR):
    """Plain one-hot singles + pair co-occurrence, not WT-centred.

    Same span as PairLR; included to show that the WT-centred *prior* is what
    helps, not extra capacity.
    """
    name = "onehot_pair_lr"
    label = "logistic, one-hot singles+pairs (18,200 feat)"
    order = 2

    def _x(self, codes):
        codes = codes.astype(np.int64)
        n = codes.shape[0]
        rows, cols = [], []
        r = np.arange(n)
        for p in range(LENGTH):
            rows.append(r)
            cols.append(p * 20 + codes[:, p])
        off = LENGTH * 20
        for blk, (p, q) in enumerate(F.PAIRS):
            rows.append(r)
            cols.append(off + blk * 400 + codes[:, p] * 20 + codes[:, q])
        row = np.concatenate(rows)
        col = np.concatenate(cols)
        return sparse.coo_matrix((np.ones(row.size, np.float32), (row, col)),
                                 shape=(n, off + len(F.PAIRS) * 400)).tocsr()


class PairRidgeOrdinal:
    """Ridge regression on the ordinal label (low=0, medium=1, high=2).

    Uses the extra resolution that the 3-class FACS call carries and that the
    binary high-vs-rest collapse throws away.  Closed form, so it is also the
    cheapest way to fit 16,435 features on 524k rows.
    """
    name = "pair_ridge_ord"
    label = "ridge on 3-class ordinal, singles+pairs"
    proba = False
    order = 2

    def __init__(self, alpha: float = 30.0):
        self.alpha = alpha

    def fit(self, codes, y, ordinal=None):
        from sklearn.linear_model import Ridge
        target = ordinal.astype(np.float64) if ordinal is not None else y.astype(np.float64)
        self.m = Ridge(alpha=self.alpha, solver="sparse_cg", max_iter=2000)
        self.m.fit(F.mutation_sparse_fast(codes, self.order), target)
        return self

    def score(self, codes):
        return self.m.predict(F.mutation_sparse_fast(codes, self.order))


class XGB:
    """Gradient boosting on one-hot + physicochemical features."""
    name = "xgb"
    label = "XGBoost (one-hot + physchem)"
    proba = True

    def __init__(self, n_estimators=800, max_depth=8, lr=0.08, physchem=True):
        self.kw = dict(n_estimators=n_estimators, max_depth=max_depth,
                       learning_rate=lr, subsample=0.8, colsample_bytree=0.8,
                       min_child_weight=5, reg_lambda=1.0, tree_method="hist",
                       n_jobs=16, eval_metric="logloss")
        self.physchem = physchem

    def _x(self, codes):
        parts = [F.onehot(codes)]
        if self.physchem:
            parts.append(F.physchem(codes))
        return np.concatenate(parts, axis=1)

    def fit(self, codes, y, ordinal=None):
        from xgboost import XGBClassifier
        self.m = XGBClassifier(**self.kw)
        self.m.fit(self._x(codes), y)
        return self

    def score(self, codes):
        return self.m.predict_proba(self._x(codes))[:, 1]


class XGBOrdinal(XGB):
    """Same trees, regressing the ordinal label instead of the binary one."""
    name = "xgb_ord"
    label = "XGBoost regressor on 3-class ordinal"
    proba = False

    def fit(self, codes, y, ordinal=None):
        from xgboost import XGBRegressor
        kw = dict(self.kw)
        kw["eval_metric"] = "rmse"
        self.m = XGBRegressor(**kw)
        self.m.fit(self._x(codes), ordinal.astype(np.float32))
        return self

    def score(self, codes):
        return self.m.predict(self._x(codes))


class CNN:
    """The published companion-repo architecture, reimplemented in torch.

    Conv1D(400, kernel 5, same, relu) -> Dropout(0.2) -> MaxPool1D(2, stride 1)
    -> Flatten -> Dense(300, relu) -> Dense(1, sigmoid), Adam lr=7.5e-5, BCE,
    input (10, 20) one-hot.
    """
    name = "cnn"
    label = "CNN, published architecture"
    proba = True

    def __init__(self, epochs=20, batch=256, lr=7.5e-5, device=None, seed=0,
                 ordinal_head=False, hidden=300, filters=400, val_frac=0.05):
        self.epochs, self.batch, self.lr, self.seed = epochs, batch, lr, seed
        self.device = device
        self.ordinal_head = ordinal_head
        self.hidden, self.filters, self.val_frac = hidden, filters, val_frac

    def _build(self):
        import torch
        from torch import nn

        class Net(nn.Module):
            def __init__(self, filters, hidden):
                super().__init__()
                self.conv = nn.Conv1d(20, filters, kernel_size=5, padding="same")
                self.drop = nn.Dropout(0.2)
                self.pool = nn.MaxPool1d(2, stride=1)
                self.fc1 = nn.Linear(filters * (LENGTH - 1), hidden)
                self.fc2 = nn.Linear(hidden, 1)

            def forward(self, x):                       # x: (n, 10, 20)
                h = torch.relu(self.conv(x.transpose(1, 2)))
                h = self.pool(self.drop(h)).flatten(1)
                return self.fc2(torch.relu(self.fc1(h))).squeeze(-1)

        return Net(self.filters, self.hidden)

    def fit(self, codes, y, ordinal=None):
        import torch
        torch.manual_seed(self.seed)
        torch.set_num_threads(16)
        dev = torch.device(self.device or ("cuda" if torch.cuda.is_available() else "cpu"))
        self._dev = dev
        net = self._build().to(dev)
        opt = torch.optim.Adam(net.parameters(), lr=self.lr)

        target = (ordinal.astype(np.float32) / 2.0 if self.ordinal_head
                  else y.astype(np.float32))
        x = torch.from_numpy(F.onehot_3d(codes))
        t = torch.from_numpy(target)

        rng = np.random.default_rng(self.seed)
        perm = rng.permutation(len(t))
        n_val = int(self.val_frac * len(t))
        val, tr = perm[:n_val], perm[n_val:]
        xv, tv = x[val].to(dev), t[val].to(dev)
        loss_fn = torch.nn.BCEWithLogitsLoss()

        best, best_state, patience = np.inf, None, 0
        for epoch in range(self.epochs):
            net.train()
            order = torch.from_numpy(rng.permutation(len(tr)))
            idx = torch.from_numpy(tr)[order]
            for s in range(0, idx.numel(), self.batch):
                b = idx[s:s + self.batch]
                xb, tb = x[b].to(dev), t[b].to(dev)
                opt.zero_grad()
                loss_fn(net(xb), tb).backward()
                opt.step()
            net.eval()
            with torch.no_grad():
                vl = float(loss_fn(net(xv), tv))
            if vl < best - 1e-5:
                best, patience = vl, 0
                best_state = {k: v.detach().clone() for k, v in net.state_dict().items()}
            else:
                patience += 1
                if patience >= 4:
                    break
        if best_state is not None:
            net.load_state_dict(best_state)
        net.eval()
        self.net = net
        self.epochs_used = epoch + 1
        return self

    def score(self, codes):
        import torch
        out = []
        x = torch.from_numpy(F.onehot_3d(codes))
        with torch.no_grad():
            for s in range(0, x.shape[0], 8192):
                out.append(torch.sigmoid(self.net(x[s:s + 8192].to(self._dev))).cpu().numpy())
        return np.concatenate(out)


class CNNOrdinal(CNN):
    name = "cnn_ord"
    label = "CNN, 3-class ordinal head"

    def __init__(self, **kw):
        kw["ordinal_head"] = True
        super().__init__(**kw)


class MLP(CNN):
    """Dense net straight on the 200-feature one-hot, no convolution.

    Tests whether the convolution in the published architecture earns its place:
    a 10-residue loop has no translational symmetry to exploit.
    """
    name = "mlp"
    label = "MLP on one-hot (512-256)"

    def __init__(self, epochs=25, batch=512, lr=1e-3, **kw):
        super().__init__(epochs=epochs, batch=batch, lr=lr, **kw)

    def _build(self):
        import torch
        from torch import nn

        class Net(nn.Module):
            def __init__(self):
                super().__init__()
                self.body = nn.Sequential(
                    nn.Linear(LENGTH * 20, 512), nn.ReLU(), nn.Dropout(0.2),
                    nn.Linear(512, 256), nn.ReLU(), nn.Dropout(0.2),
                    nn.Linear(256, 1))

            def forward(self, x):
                return self.body(x.flatten(1)).squeeze(-1)

        return Net()


MODELS = {m.name: m for m in [LogOdds, SingleLR, PairLR, TripleLR, OneHotPairLR,
                              PairRidgeOrdinal, XGB, XGBOrdinal, CNN, CNNOrdinal, MLP]}


# ------------------------------------------------------------------ evaluation

def report(scores, y, dist, label, proba=False, boot=2000, quiet=False):
    """Print the stratified table and return the row for the summary."""
    overall, lo, hi = auc_ci(scores, y, n_boot=boot)
    row = {"model": label, "auc_overall": overall, "n": int(y.size)}
    if not quiet:
        print(f"\n  {label}")
        print(f"    overall AUC {overall:.4f} [{lo:.4f}, {hi:.4f}]   n={y.size:,}")
        print(f"    {'dist':>4} {'n':>8} {'%bind':>6} {'AUC':>7} "
              f"{'95% CI':>18} {'AP':>6}")
    for d in range(1, 11):
        m = dist == d
        if m.sum() < 10 or y[m].min() == y[m].max():
            continue
        a, l, h = auc_ci(scores[m], y[m], n_boot=boot, seed=d)
        ap = average_precision(scores[m], y[m])
        if d == 3:
            row.update(auc_d3=a, auc_d3_lo=l, auc_d3_hi=h, ap_d3=ap, n_d3=int(m.sum()))
        if d == 4:
            row.update(auc_d4=a)
        if not quiet:
            flag = "  <- the legal shell" if d == 3 else ""
            print(f"    {d:>4} {m.sum():>8,} {100*y[m].mean():>5.1f}% {a:>7.4f} "
                  f"  [{l:.3f}, {h:.3f}]  {ap:>6.3f}{flag}")
    if proba and not quiet:
        e, table = ece(scores, y)
        row["ece"] = e
        print(f"    calibration: ECE {e:.4f}")
        print(f"      {'bin':>12} {'n':>8} {'pred':>6} {'actual':>7}")
        for a0, a1, n, conf, frac in table:
            print(f"      {a0:.1f}-{a1:.1f}    {n:>8,} {conf:>6.3f} {frac:>7.3f}")
    elif proba:
        row["ece"] = ece(scores, y)[0]
    return row


def run_cv(model_cls, facs, folds=5, seed=0, kw=None, verbose=True):
    """Pooled out-of-fold scores from k-fold random CV."""
    codes, y, ordinal, dist = facs["codes"], facs["y"], facs["ordinal"], facs["dist"]
    oof = np.full(y.size, np.nan)
    t0 = time.time()
    for i, (tr, te) in enumerate(kfold_indices(y.size, folds, seed=seed)):
        m = model_cls(**(kw or {})).fit(codes[tr], y[tr], ordinal[tr])
        oof[te] = m.score(codes[te])
        if verbose:
            print(f"      fold {i+1}/{folds} done ({time.time()-t0:.0f}s)", flush=True)
    assert not np.isnan(oof).any()
    return oof, time.time() - t0


def run_shell(model_cls, facs, kw=None, train_min_dist=5):
    """Train on edit distance >= train_min_dist, test on everything closer."""
    codes, y, ordinal, dist = facs["codes"], facs["y"], facs["ordinal"], facs["dist"]
    tr = np.flatnonzero(dist >= train_min_dist)
    te = np.flatnonzero(dist < train_min_dist)
    m = model_cls(**(kw or {})).fit(codes[tr], y[tr], ordinal[tr])
    return m.score(codes[te]), te


def run_spr_transfer(model_cls, facs, spr, kw=None):
    """Fit on FACS with the SPR sequences removed, then correlate with log10(KD)."""
    from neighbourhood import pack
    codes, y, ordinal = facs["codes"], facs["y"], facs["ordinal"]
    spr_keys = set(pack(spr["codes"]).tolist())
    keep = np.array([k not in spr_keys for k in pack(codes).tolist()])
    m = model_cls(**(kw or {})).fit(codes[keep], y[keep], ordinal[keep])
    s = m.score(spr["codes"])
    # higher score should mean tighter binding, i.e. lower log10(KD)
    rho = spearman(s, -spr["log10kd"])
    nonwt = ~spr["is_wt"]
    rho_nonwt = spearman(s[nonwt], -spr["log10kd"][nonwt])
    return rho, rho_nonwt, int(keep.sum()), s


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", default="logodds,single_lr,pair_lr",
                    help="comma list or 'all'; choices: " + ",".join(MODELS))
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--regimes", default="cv,shell,spr")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--tag", default="main")
    ap.add_argument("--device", default=None)
    args = ap.parse_args()

    names = list(MODELS) if args.models == "all" else args.models.split(",")
    bad = [n for n in names if n not in MODELS]
    if bad:
        print(f"unknown models: {bad}", file=sys.stderr)
        return 2
    regimes = set(args.regimes.split(","))

    facs = load_facs()
    spr = load_spr()
    print(f"FACS {facs['y'].size:,} rows, {facs['y'].mean():.1%} binder; "
          f"SPR {len(spr['seqs'])} sequences with KD")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    summary = []
    for name in names:
        cls = MODELS[name]
        kw = {}
        if args.device and name in ("cnn", "cnn_ord", "mlp"):
            kw["device"] = args.device
        print(f"\n{'='*78}\n{cls.label}  [{name}]\n{'='*78}", flush=True)
        row = {"name": name, "label": cls.label}

        if "cv" in regimes:
            oof, secs = run_cv(cls, facs, args.folds, args.seed, kw)
            np.save(RESULTS_DIR / f"oof_{name}.npy", oof)
            r = report(oof, facs["y"], facs["dist"], f"{args.folds}-fold CV (pooled OOF)",
                       proba=cls.proba, boot=args.boot)
            row.update({f"cv_{k}": v for k, v in r.items() if k != "model"})
            row["cv_seconds"] = round(secs, 1)

        if "shell" in regimes:
            s, te = run_shell(cls, facs, kw)
            r = report(s, facs["y"][te], facs["dist"][te],
                       "shell extrapolation (train dist>=5, test dist<=4)",
                       proba=cls.proba, boot=args.boot)
            row.update({f"shell_{k}": v for k, v in r.items() if k != "model"})

        if "spr" in regimes:
            rho, rho_nw, n_tr, s = run_spr_transfer(cls, facs, spr, kw)
            np.save(RESULTS_DIR / f"spr_score_{name}.npy", s)
            print(f"\n    SPR transfer (trained on {n_tr:,} FACS rows, SPR sequences held out)")
            print(f"      Spearman(score, -log10 KD): all {len(spr['seqs'])} = {rho:+.3f}; "
                  f"non-WT = {rho_nw:+.3f}")
            row.update(spr_rho=rho, spr_rho_nonwt=rho_nw)

        summary.append(row)
        with (RESULTS_DIR / f"summary_{args.tag}.json").open("w") as fh:
            json.dump(summary, fh, indent=2)

    print(f"\n\n{'='*100}\nSUMMARY -- distance-3 AUC is the number that counts\n{'='*100}")
    nan = float("nan")
    print(f"{'model':<42}{'CV d3':>7} {'95% CI':<16}{'CV all':>8}"
          f"{'shell d3':>9}{'SPRrho':>8}{'ECE':>7}{'sec':>7}")
    for r in summary:
        ci = f"[{r.get('cv_auc_d3_lo', nan):.3f},{r.get('cv_auc_d3_hi', nan):.3f}]"
        print(f"{r['label'][:41]:<42}{r.get('cv_auc_d3', nan):>7.4f} {ci:<16}"
              f"{r.get('cv_auc_overall', nan):>8.4f}"
              f"{r.get('shell_auc_d3', nan):>9.4f}"
              f"{r.get('spr_rho_nonwt', nan):>8.3f}"
              f"{r.get('cv_ece', nan):>7.3f}"
              f"{r.get('cv_seconds', nan):>7.0f}")
    print(f"\nwrote {RESULTS_DIR / f'summary_{args.tag}.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
