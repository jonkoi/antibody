#!/usr/bin/env python3
"""Calibrate ipTM as a binder/non-binder discriminator.

Reads the ipTM matrices from Smorodina et al. (Zenodo 18390239,
iptm_confidence_scores.zip): one CSV per folding model, rows = nanobody PDB ID,
columns = antigen PDB ID. Cells on the label-matched diagonal are cognate
(real) pairs; every other cell is a shuffled non-cognate pair.

For each model it reports:
  1. cognate vs non-cognate ipTM distributions (n, mean, median, std, quartiles)
  2. ROC AUC (rank-based Mann-Whitney U, tie-corrected) + Hanley-McNeil SE
  3. confusion / precision / recall at ipTM > 0.5 and > 0.7
  4. top-N precision vs the chance baseline (the "rank and pick" question)

Pure standard library: no numpy/pandas/sklearn/scipy required.

Usage:
    python3 -I calibrate_iptm.py <dir-with-*_matrix_clean.csv> [--csv out.csv]
"""

import csv
import math
import os
import sys

MODELS = [
    ("AF3", "af3_matrix_clean.csv"),
    ("Boltz-2", "boltz2_matrix_clean.csv"),
    ("Chai-1", "chai1_matrix_clean.csv"),
]

THRESHOLDS = (0.5, 0.7)
TOP_NS = (10, 20, 50, 91, 100, 200, 500)

# ipTM values our own HER2 designs scored, to be placed in each distribution.
QUERY_SCORES = (0.90, 0.946, 0.965)


# ---------------------------------------------------------------- loading


def load_matrix(path):
    """Return (col_labels, {row_label: {col_label: float|None}}, n_bad)."""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    header = rows[0]
    col_labels = [c.strip() for c in header[1:]]
    matrix = {}
    n_bad = 0
    for raw in rows[1:]:
        if not raw or not raw[0].strip():
            continue
        rlabel = raw[0].strip()
        cells = {}
        for j, clabel in enumerate(col_labels):
            raw_val = raw[j + 1].strip() if j + 1 < len(raw) else ""
            try:
                cells[clabel] = float(raw_val)
            except ValueError:
                cells[clabel] = None
                n_bad += 1
        matrix[rlabel] = cells
    return col_labels, matrix, n_bad


def split_cognate(col_labels, matrix):
    """Split cells into (cognate, non_cognate) by LABEL match, not position."""
    row_labels = list(matrix)
    colset = set(col_labels)
    rowset = set(row_labels)
    shared = [r for r in row_labels if r in colset]

    cognate, non_cognate = [], []
    for r in row_labels:
        for c in col_labels:
            v = matrix[r][c]
            if v is None:
                continue
            (cognate if r == c else non_cognate).append(v)

    diag_in_order = all(
        r == c for r, c in zip(row_labels, col_labels)
    ) and len(row_labels) == len(col_labels)

    meta = {
        "n_rows": len(row_labels),
        "n_cols": len(col_labels),
        "n_shared_labels": len(shared),
        "rows_not_in_cols": sorted(rowset - colset),
        "cols_not_in_rows": sorted(colset - rowset),
        "positional_diagonal_matches_labels": diag_in_order,
    }
    return cognate, non_cognate, meta


# ---------------------------------------------------------------- stats


def quantile(sorted_vals, q):
    """Linear-interpolation quantile (same convention as numpy default)."""
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def describe(vals):
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return {}
    mean = sum(s) / n
    var = sum((v - mean) ** 2 for v in s) / (n - 1) if n > 1 else 0.0
    return {
        "n": n,
        "mean": mean,
        "median": quantile(s, 0.50),
        "std": math.sqrt(var),
        "min": s[0],
        "q25": quantile(s, 0.25),
        "q75": quantile(s, 0.75),
        "q90": quantile(s, 0.90),
        "max": s[-1],
    }


def auc_mannwhitney(pos, neg):
    """Tie-corrected ROC AUC via mid-ranks. Equals P(pos>neg)+0.5*P(tie)."""
    n1, n0 = len(pos), len(neg)
    if n1 == 0 or n0 == 0:
        return float("nan"), float("nan")

    combined = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    ranks = [0.0] * len(combined)
    i = 0
    while i < len(combined):
        j = i
        while j + 1 < len(combined) and combined[j + 1][0] == combined[i][0]:
            j += 1
        mid = (i + j) / 2.0 + 1.0  # 1-based mid-rank
        for k in range(i, j + 1):
            ranks[k] = mid
        i = j + 1

    r1 = sum(rk for rk, (_, lab) in zip(ranks, combined) if lab == 1)
    u1 = r1 - n1 * (n1 + 1) / 2.0
    auc = u1 / (n1 * n0)

    # Hanley & McNeil (1982) standard error.
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc * auc / (1.0 + auc)
    se_sq = (
        auc * (1 - auc) + (n1 - 1) * (q1 - auc**2) + (n0 - 1) * (q2 - auc**2)
    ) / (n1 * n0)
    return auc, math.sqrt(max(se_sq, 0.0))


def threshold_stats(pos, neg, thr):
    tp = sum(1 for v in pos if v > thr)
    fn = len(pos) - tp
    fp = sum(1 for v in neg if v > thr)
    tn = len(neg) - fp
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision == precision and recall == recall and (precision + recall) > 0
        else float("nan")
    )
    return {
        "thr": thr,
        "tp": tp,
        "fn": fn,
        "fp": fp,
        "tn": tn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "fpr": fp / (fp + tn) if (fp + tn) else float("nan"),
    }


def top_n_precision(pos, neg, ns):
    """Rank every pair by ipTM descending; what fraction of top N is cognate?

    Ties straddling the cut are handled by expectation (the fraction of cognate
    among the tied block is credited proportionally), which is what you would
    get on average from an arbitrary tie-break.
    """
    labelled = sorted(
        [(v, 1) for v in pos] + [(v, 0) for v in neg], key=lambda t: -t[0]
    )
    total = len(labelled)
    prevalence = len(pos) / total if total else float("nan")

    out = []
    for n in ns:
        if n > total:
            continue
        # Expected number of positives in the top n with proportional ties.
        taken = 0
        hits = 0.0
        i = 0
        while taken < n and i < total:
            j = i
            while j + 1 < total and labelled[j + 1][0] == labelled[i][0]:
                j += 1
            block = labelled[i : j + 1]
            block_pos = sum(lab for _, lab in block)
            room = n - taken
            if len(block) <= room:
                hits += block_pos
                taken += len(block)
            else:
                hits += block_pos * room / len(block)
                taken += room
            i = j + 1
        prec = hits / n
        out.append(
            {
                "n": n,
                "expected_cognate": hits,
                "precision": prec,
                "chance": prevalence,
                "enrichment": prec / prevalence if prevalence else float("nan"),
            }
        )
    return prevalence, out


# ---------------------------------------------------------------- reporting


def fmt(x, nd=3):
    if x != x:
        return "  n/a"
    return f"{x:.{nd}f}"


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    data_dir = argv[1]
    csv_out = None
    if "--csv" in argv:
        csv_out = argv[argv.index("--csv") + 1]

    results = []
    for name, fname in MODELS:
        path = os.path.join(data_dir, fname)
        if not os.path.exists(path):
            print(f"!! missing {path}", file=sys.stderr)
            continue
        col_labels, matrix, n_bad = load_matrix(path)
        cognate, non_cognate, meta = split_cognate(col_labels, matrix)
        results.append(
            {
                "model": name,
                "file": fname,
                "n_bad_cells": n_bad,
                "meta": meta,
                "cognate": cognate,
                "non_cognate": non_cognate,
            }
        )

    line = "=" * 78

    # --- 0. sanity -------------------------------------------------------
    print(line)
    print("0. MATRIX SANITY CHECK (diagonal identified by LABEL match)")
    print(line)
    for r in results:
        m = r["meta"]
        print(
            f"{r['model']:<9} {m['n_rows']}x{m['n_cols']} matrix | "
            f"shared labels: {m['n_shared_labels']} | "
            f"non-numeric/empty cells: {r['n_bad_cells']}"
        )
        print(
            f"{'':<9} cognate (diagonal) n={len(r['cognate'])}, "
            f"non-cognate n={len(r['non_cognate'])} | "
            f"row order == col order: {m['positional_diagonal_matches_labels']}"
        )
        if m["rows_not_in_cols"]:
            print(f"{'':<9} rows absent from columns: {m['rows_not_in_cols']}")
        if m["cols_not_in_rows"]:
            print(f"{'':<9} columns absent from rows: {m['cols_not_in_rows']}")

    # --- 1. distributions ------------------------------------------------
    print()
    print(line)
    print("1. ipTM DISTRIBUTIONS: cognate (real) vs non-cognate (shuffled)")
    print(line)
    hdr = (
        f"{'model':<9} {'class':<12} {'n':>6} {'mean':>7} {'median':>7} "
        f"{'std':>7} {'min':>6} {'q25':>6} {'q75':>6} {'q90':>6} {'max':>6}"
    )
    print(hdr)
    print("-" * len(hdr))
    for r in results:
        for cls, vals in (("cognate", r["cognate"]), ("non-cognate", r["non_cognate"])):
            d = describe(vals)
            print(
                f"{r['model']:<9} {cls:<12} {d['n']:>6} {fmt(d['mean']):>7} "
                f"{fmt(d['median']):>7} {fmt(d['std']):>7} {fmt(d['min'],2):>6} "
                f"{fmt(d['q25'],2):>6} {fmt(d['q75'],2):>6} {fmt(d['q90'],2):>6} "
                f"{fmt(d['max'],2):>6}"
            )
        print()

    # --- 2. AUC ----------------------------------------------------------
    print(line)
    print("2. ROC AUC (ipTM as the only score; 0.5 = coin flip)")
    print(line)
    print(f"{'model':<9} {'AUC':>7} {'SE':>7} {'95% CI':>18} {'delta-median':>13}")
    print("-" * 58)
    for r in results:
        auc, se = auc_mannwhitney(r["cognate"], r["non_cognate"])
        lo, hi = auc - 1.96 * se, auc + 1.96 * se
        dmed = describe(r["cognate"])["median"] - describe(r["non_cognate"])["median"]
        r["auc"], r["auc_se"] = auc, se
        print(
            f"{r['model']:<9} {fmt(auc):>7} {fmt(se):>7} "
            f"{'[' + fmt(lo) + ', ' + fmt(hi) + ']':>18} {fmt(dmed):>13}"
        )

    # --- 3. thresholds ---------------------------------------------------
    print()
    print(line)
    print("3. THRESHOLD ANALYSIS (positive call = ipTM > threshold)")
    print(line)
    hdr = (
        f"{'model':<9} {'thr':>5} {'TP':>5} {'FN':>5} {'FP':>6} {'TN':>6} "
        f"{'precision':>10} {'recall':>8} {'FPR':>7} {'F1':>6}"
    )
    print(hdr)
    print("-" * len(hdr))
    for r in results:
        for thr in THRESHOLDS:
            t = threshold_stats(r["cognate"], r["non_cognate"], thr)
            print(
                f"{r['model']:<9} {thr:>5.1f} {t['tp']:>5} {t['fn']:>5} "
                f"{t['fp']:>6} {t['tn']:>6} {fmt(t['precision']):>10} "
                f"{fmt(t['recall']):>8} {fmt(t['fpr']):>7} {fmt(t['f1'],2):>6}"
            )
        print()

    # --- 4. top-N --------------------------------------------------------
    print(line)
    print("4. RANK-BY-ipTM AND TAKE TOP N: fraction that are real binders")
    print(line)
    for r in results:
        prev, rows = top_n_precision(r["cognate"], r["non_cognate"], TOP_NS)
        print(f"{r['model']}  (chance / prevalence = {prev*100:.2f}% cognate)")
        print(
            f"  {'top N':>6} {'exp. real binders':>18} {'precision':>10} "
            f"{'vs chance':>10}"
        )
        for row in rows:
            print(
                f"  {row['n']:>6} {row['expected_cognate']:>18.1f} "
                f"{row['precision']*100:>9.1f}% {row['enrichment']:>9.1f}x"
            )
        print()

    # --- 5. what is a given score worth? ---------------------------------
    print(line)
    print("5. CALIBRATION OF OUR OWN DESIGN SCORES")
    print(line)
    print(
        "For a query ipTM q: what fraction of SHUFFLED (known non-binding) pairs\n"
        "score at least q, and what fraction of real cognate pairs do?\n"
        "A score only means something if non-cognate pairs rarely reach it."
    )
    print()
    hdr = (
        f"{'model':<9} {'query':>6} {'%non-cognate >= q':>19} "
        f"{'%cognate >= q':>15} {'LR+':>7}"
    )
    print(hdr)
    print("-" * len(hdr))
    for r in results:
        for q in QUERY_SCORES:
            nc = sum(1 for v in r["non_cognate"] if v >= q) / len(r["non_cognate"])
            co = sum(1 for v in r["cognate"] if v >= q) / len(r["cognate"])
            lr = co / nc if nc > 0 else float("inf")
            lr_s = "inf" if lr == float("inf") else f"{lr:.1f}"
            print(
                f"{r['model']:<9} {q:>6.3f} {nc*100:>18.2f}% "
                f"{co*100:>14.1f}% {lr_s:>7}"
            )
        print()

    # --- optional machine-readable ---------------------------------------
    if csv_out:
        with open(csv_out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(
                [
                    "model",
                    "n_cognate",
                    "n_non_cognate",
                    "cognate_mean",
                    "cognate_median",
                    "noncognate_mean",
                    "noncognate_median",
                    "auc",
                    "auc_se",
                ]
                + [
                    f"{k}_{thr}"
                    for thr in THRESHOLDS
                    for k in ("tp", "fn", "fp", "precision", "recall")
                ]
            )
            for r in results:
                dc, dn = describe(r["cognate"]), describe(r["non_cognate"])
                row = [
                    r["model"],
                    dc["n"],
                    dn["n"],
                    dc["mean"],
                    dc["median"],
                    dn["mean"],
                    dn["median"],
                    r["auc"],
                    r["auc_se"],
                ]
                for thr in THRESHOLDS:
                    t = threshold_stats(r["cognate"], r["non_cognate"], thr)
                    row += [t["tp"], t["fn"], t["fp"], t["precision"], t["recall"]]
                w.writerow(row)
        print(f"wrote {csv_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
