#!/usr/bin/env python3
"""Score all 51,417 legal 3-mutant CDRH3 candidates and write the ranked table.

Dependencies: numpy, scipy, scikit-learn (see .venv-ml).

Usage:
    .venv-ml/bin/python pipeline/score_candidates_ml.py \
        --out data/legal_3mutants_scored_ml.csv

TWO STAGES, BECAUSE THE TWO DATASETS ANSWER DIFFERENT QUESTIONS

Stage 1 -- will it still bind at all?
    P(FACS "high") from the WT-centred pairwise epistatic logistic model fitted
    on all 524,346 FACS variants.  Pooled out-of-fold AUC in the distance-3
    shell is 0.889 [0.855, 0.920] against 0.684 for the incumbent
    position-independent log-odds model.  This stage is well validated.

Stage 2 -- how tight, in nM?
    Ridge regression of log10(KD) on edit distance plus the WT-centred
    single-mutation coding, fitted on the 180-sequence SPR panel.  This stage is
    weak and is reported with an interval, not a point estimate.  The interval
    is a split-conformal band taken from the repeated-CV residuals, so its
    stated coverage is empirical rather than assumed.

Combined -- P(beats trastuzumab)
    p_bind * Phi((log10(0.146) - mu) / sigma), with sigma the conformal residual
    scale.  It comes out small for every candidate, which is the honest answer:
    only 1 of the 19 measured distance-3 variants in the SPR panel beats 0.146 nM
    and the model cannot resolve a 2-fold affinity difference when its own
    residual spread is a factor of ten.

    The FACS score is deliberately NOT used to rank affinity.  On the 27
    sequences present in both datasets it ranks KD at Spearman +0.19, and within
    the distance-3 shell (n=19) the correlation is negative.  Measured, in
    train_affinity_model.py.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))

import facs_scorer
import features as F
from binding_data import (REPO, WT_CDRH3, WT_KD_MEDIAN_NM, load_candidates,
                          load_spr, spr_replicate_noise)

AFFINITY_BLOCKS = ("dist", "mut")
AFFINITY_ALPHA_GRID = (1, 3, 10, 30, 100, 300, 1000)


def affinity_blocks(codes, dist):
    return np.concatenate([
        dist.astype(np.float64).reshape(-1, 1),
        F.mutation_sparse_fast(codes, 1).toarray().astype(np.float64),
    ], axis=1)


def fit_affinity(spr):
    """Ridge on [distance, single-mutation coding]; alpha by 10-fold CV.

    Returns (model, sigma) where sigma is the residual sd from repeated
    out-of-fold prediction -- the honest predictive spread, not the in-sample one.
    """
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import KFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    x = affinity_blocks(spr["codes"], spr["dist"])
    y = spr["log10kd"]

    def build(alpha):
        return make_pipeline(StandardScaler(), Ridge(alpha=alpha))

    best, best_err = None, np.inf
    for alpha in AFFINITY_ALPHA_GRID:
        err = 0.0
        for rep in range(10):
            for tr, te in KFold(10, shuffle=True, random_state=rep).split(x):
                err += float(((build(alpha).fit(x[tr], y[tr]).predict(x[te]) - y[te]) ** 2).sum())
        if err < best_err:
            best, best_err = alpha, err

    resid = []
    for rep in range(10):
        oof = np.empty(y.size)
        for tr, te in KFold(10, shuffle=True, random_state=100 + rep).split(x):
            oof[te] = build(best).fit(x[tr], y[tr]).predict(x[te])
        resid.append(oof - y)
    resid = np.concatenate(resid)
    sigma = float(np.std(resid, ddof=1))
    q90 = float(np.quantile(np.abs(resid), 0.90))
    return build(best).fit(x, y), best, sigma, q90


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=REPO / "data/legal_3mutants_scored_ml.csv")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    rows, seqs, codes = load_candidates()
    dist = np.full(len(seqs), 3, dtype=np.int64)
    bad = [s for s in seqs if sum(a != b for a, b in zip(s, WT_CDRH3)) != 3]
    assert not bad, f"{len(bad)} candidates are not 3-mutants"
    print(f"{len(seqs):,} candidates, all at edit distance 3")

    print("fitting the FACS epistatic model on all 524,346 variants ...", flush=True)
    models = facs_scorer.get_models(holdout_spr=False)
    x2 = F.mutation_sparse_fast(codes, 2)
    p_bind = models["pair_lr_logit"].predict_proba(x2)[:, 1]
    logit = models["pair_lr_logit"].decision_function(x2)
    logodds = models["logodds"].score(codes)
    ord_score = models["pair_ridge_ord"].predict(x2)

    print("fitting the SPR affinity model on 180 measured KDs ...", flush=True)
    spr = load_spr()
    aff, alpha, sigma, q90 = fit_affinity(spr)
    mu = aff.predict(affinity_blocks(codes, dist))
    print(f"  ridge alpha={alpha}; out-of-fold residual sd = {sigma:.3f} log10 units "
          f"(90% of |residuals| within {q90:.3f}, i.e. a factor of {10**q90:.0f})")
    print(f"  wild-type replicate noise for comparison: {spr_replicate_noise():.3f} log10")

    bar = float(np.log10(WT_KD_MEDIAN_NM))
    p_beat_given_bind = stats.norm.cdf((bar - mu) / sigma)
    p_beat = p_bind * p_beat_given_bind

    print(f"\npredicted log10(KD) over the candidates: "
          f"{mu.min():+.2f} .. {mu.max():+.2f} (median {np.median(mu):+.2f}); "
          f"the bar is {bar:+.2f}")
    print(f"P(beats 0.146 nM): max {p_beat.max():.4f}, median {np.median(p_beat):.4f}")
    print(f"candidates with P(FACS binder) > 0.9: {(p_bind > 0.9).sum():,} "
          f"({(p_bind > 0.9).mean():.1%})")

    order = np.argsort(-p_beat)
    fields = ["rank", "cdrh3", "substitutions", "removes_dg", "pro_for_gly",
              "p_beat_wt", "facs_p_bind", "pred_log10_kd", "pred_kd_nM",
              "pred_kd_nM_lo90", "pred_kd_nM_hi90",
              "facs_pair_logit", "facs_ordinal", "logodds_incumbent"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for rank, i in enumerate(order, 1):
            w.writerow({
                "rank": rank, "cdrh3": seqs[i],
                "substitutions": rows[i].get("substitutions", ""),
                "removes_dg": rows[i].get("removes_dg", ""),
                "pro_for_gly": rows[i].get("pro_for_gly", ""),
                "p_beat_wt": f"{p_beat[i]:.6f}",
                "facs_p_bind": f"{p_bind[i]:.5f}",
                "pred_log10_kd": f"{mu[i]:.4f}",
                "pred_kd_nM": f"{10 ** mu[i]:.4g}",
                "pred_kd_nM_lo90": f"{10 ** (mu[i] - q90):.4g}",
                "pred_kd_nM_hi90": f"{10 ** (mu[i] + q90):.4g}",
                "facs_pair_logit": f"{logit[i]:.4f}",
                "facs_ordinal": f"{ord_score[i]:.4f}",
                "logodds_incumbent": f"{logodds[i]:.4f}",
            })
    print(f"\nwrote {args.out}")

    print(f"\ntop {args.top} by P(beats trastuzumab):")
    print(f"  {'#':>3}  {'cdrh3':<11}{'subs':<22}{'P(beat)':>9}{'P(bind)':>9}"
          f"{'pred KD nM':>11}  {'90% band nM':>20} DG P/G")
    for rank, i in enumerate(order[:args.top], 1):
        print(f"  {rank:>3}  {seqs[i]:<11}{rows[i].get('substitutions',''):<22}"
              f"{p_beat[i]:>9.4f}{p_bind[i]:>9.4f}{10**mu[i]:>11.3g}"
              f"   {10**(mu[i]-q90):>8.3g} - {10**(mu[i]+q90):<8.3g}"
              f"  {'Y' if rows[i].get('removes_dg')=='1' else '-'}"
              f"  {'Y' if rows[i].get('pro_for_gly')=='1' else '-'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
