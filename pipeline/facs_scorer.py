#!/usr/bin/env python3
"""Fit the FACS binder models on the full screen and score arbitrary CDRH3s.

Dependencies: numpy, scipy, scikit-learn, xgboost.

Exists so that the SPR affinity regression and the final candidate scoring both
use *the same* FACS-derived features, fitted once and cached, and so that the
sequences carrying an SPR measurement can be excluded from the FACS fit in one
place rather than in two.

``facs_features(codes, holdout_spr=True)`` returns a (n, k) array of FACS-model
scores.  With ``holdout_spr=True`` the 27 CDRH3s that appear in both the FACS
screen and the SPR panel are dropped from the training set, so the returned
scores are honest features for an SPR model.
"""
from __future__ import annotations

import numpy as np

import features as F
from binding_data import REPO, load_facs, load_spr
from neighbourhood import pack

CACHE = REPO / "data/cache"
FEATURE_NAMES = ("logodds", "pair_lr_logit", "pair_ridge_ord", "single_lr_logit")


def _fit_models(codes, y, ordinal):
    from sklearn.linear_model import LogisticRegression, Ridge
    from train_binding_model import LogOdds

    out = {}
    out["logodds"] = LogOdds().fit(codes, y)
    x1 = F.mutation_sparse_fast(codes, 1)
    out["single_lr_logit"] = LogisticRegression(C=3.0, max_iter=300).fit(x1, y)
    x2 = F.mutation_sparse_fast(codes, 2)
    out["pair_lr_logit"] = LogisticRegression(C=3.0, max_iter=300).fit(x2, y)
    out["pair_ridge_ord"] = Ridge(alpha=30.0, solver="sparse_cg",
                                  max_iter=2000).fit(x2, ordinal.astype(np.float64))
    return out


def _apply(models, codes):
    cols = []
    for name in FEATURE_NAMES:
        m = models[name]
        if name == "logodds":
            cols.append(m.score(codes))
        elif name == "single_lr_logit":
            cols.append(m.decision_function(F.mutation_sparse_fast(codes, 1)))
        elif name == "pair_lr_logit":
            cols.append(m.decision_function(F.mutation_sparse_fast(codes, 2)))
        else:
            cols.append(m.predict(F.mutation_sparse_fast(codes, 2)))
    return np.stack(cols, axis=1)


_CACHED: dict[bool, object] = {}


def get_models(holdout_spr: bool = True):
    if holdout_spr in _CACHED:
        return _CACHED[holdout_spr]
    facs = load_facs()
    codes, y, ordinal = facs["codes"], facs["y"], facs["ordinal"]
    if holdout_spr:
        spr_keys = set(pack(load_spr()["codes"]).tolist())
        keep = np.fromiter((k not in spr_keys for k in pack(codes).tolist()),
                           bool, codes.shape[0])
        codes, y, ordinal = codes[keep], y[keep], ordinal[keep]
    models = _fit_models(codes, y, ordinal)
    _CACHED[holdout_spr] = models
    return models


def facs_features(codes: np.ndarray, holdout_spr: bool = True,
                  cache_key: str | None = None) -> np.ndarray:
    if cache_key:
        path = CACHE / f"facs_feat_{cache_key}_{int(holdout_spr)}.npy"
        if path.exists():
            return np.load(path)
    out = _apply(get_models(holdout_spr), codes)
    if cache_key:
        CACHE.mkdir(parents=True, exist_ok=True)
        np.save(CACHE / f"facs_feat_{cache_key}_{int(holdout_spr)}.npy", out)
    return out
