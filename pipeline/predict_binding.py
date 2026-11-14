#!/usr/bin/env python3
"""Rank CDRH3 designs with a model trained on MEASURED binding (Lane D).

Why: findings 7 and 8 (findings/lane-b-confidence.md) ruled out folding
confidence as a way to rank this pool, and min_ipae is not calibrated against
anything measured. This gives a ranking signal that is grounded in experiment.

Data: Mehta et al. (zenodo 10831512) FACS-sorted trastuzumab CDRH3 variants,
524,346 sequences labelled 1 = high binder (32.8%). The authors' own
clonotype-aware splits are used as shipped (367,044 / 78,651 / 78,651): similar
CDRH3s are kept inside one split, so held-out performance is not inflated by
near-duplicates of training rows.

Model: logistic regression on a 10-residue CDRH3, stdlib only.
    unary     position x amino acid                (10 x 20  =   200 weights)
    pairwise  position pair x amino acid pair      (45 x 400 = 18,000 weights)
Only 10 unary + 45 pairwise features are active per sequence, so an epoch over
367k rows is cheap without numpy. AdaGrad, L2, epoch chosen on validation AUC.
The test split is scored once, at the end.

Baselines reported alongside, because a complicated model that cannot beat a
position-weight matrix is not worth trusting:
    majority  predict the majority class
    pwm       per-position log-odds, positions treated as independent

Usage:
    python3 pipeline/predict_binding.py
    python3 pipeline/predict_binding.py --epochs 8 --score runs/team15-claude-n100-seed1/shortlist.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import re
from collections import defaultdict
from pathlib import Path

DEFAULT_DATA = Path("hackathon-data/zenodo-10831512-trastuzumab-ngs/extracted")
WILD_TYPE = "WGGDGFYAMD"  # heavy B99-108, the designed window
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AA)}
LEN = 10
PAIRS = [(i, j) for i in range(LEN) for j in range(i + 1, LEN)]
PAIR_IDX = {p: k for k, p in enumerate(PAIRS)}
N_UNARY = LEN * len(AA)
N_WEIGHTS = N_UNARY + len(PAIRS) * len(AA) * len(AA)


def featurise(seq: str) -> list[int] | None:
    """Active feature indices for one CDRH3, or None if unusable."""
    if len(seq) != LEN or any(c not in AA_IDX for c in seq):
        return None
    idx = [p * len(AA) + AA_IDX[c] for p, c in enumerate(seq)]
    for (i, j), k in PAIR_IDX.items():
        idx.append(N_UNARY + k * 400 + AA_IDX[seq[i]] * len(AA) + AA_IDX[seq[j]])
    return idx


def load_split(path: Path) -> tuple[list[str], list[int], list[str], int]:
    seqs: list[str] = []
    labels: list[int] = []
    clusters: list[str] = []
    dropped = 0
    with path.open() as fh:
        for row in csv.DictReader(fh):
            seq = row["seq"].strip().upper()
            if featurise(seq) is None:
                dropped += 1
                continue
            seqs.append(seq)
            labels.append(int(row["label"]))
            clusters.append(row.get("cluster_centre", ""))
    return seqs, labels, clusters, dropped


def auc(scores: list[float], labels: list[int]) -> float:
    pairs = sorted(zip(scores, labels))
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if not n_pos or not n_neg:
        return float("nan")
    rank_sum = 0.0
    i = 0
    while i < len(pairs):
        j = i
        while j + 1 < len(pairs) and pairs[j + 1][0] == pairs[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        rank_sum += avg * sum(1 for k in range(i, j + 1) if pairs[k][1] == 1)
        i = j + 1
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def auc_ci(a: float, n_pos: int, n_neg: int) -> tuple[float, float]:
    """Hanley-McNeil 95% CI."""
    if not (0 < a < 1) or n_pos < 2 or n_neg < 2:
        return float("nan"), float("nan")
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    se = ((a * (1 - a) + (n_pos - 1) * (q1 - a * a) + (n_neg - 1) * (q2 - a * a)) / (n_pos * n_neg)) ** 0.5
    return a - 1.96 * se, a + 1.96 * se


def train_pwm(seqs: list[str], labels: list[int]) -> list[float]:
    """Per-position log-odds of high vs low binding, Laplace-smoothed."""
    pos = [[1.0] * len(AA) for _ in range(LEN)]
    neg = [[1.0] * len(AA) for _ in range(LEN)]
    for seq, y in zip(seqs, labels):
        table = pos if y == 1 else neg
        for p, c in enumerate(seq):
            table[p][AA_IDX[c]] += 1.0
    weights = [0.0] * N_UNARY
    for p in range(LEN):
        tot_p, tot_n = sum(pos[p]), sum(neg[p])
        for a in range(len(AA)):
            weights[p * len(AA) + a] = math.log((pos[p][a] / tot_p) / (neg[p][a] / tot_n))
    return weights


def score_pwm(seq: str, weights: list[float]) -> float:
    return sum(weights[p * len(AA) + AA_IDX[c]] for p, c in enumerate(seq))


def sigmoid(z: float) -> float:
    if z >= 0:
        return 1.0 / (1.0 + math.exp(-z)) if z < 60 else 1.0
    return math.exp(z) / (1.0 + math.exp(z)) if z > -60 else 0.0


def train_logreg(
    feats: list[list[int]],
    labels: list[int],
    val_feats: list[list[int]],
    val_labels: list[int],
    epochs: int,
    lr: float,
    l2: float,
    seed: int,
) -> tuple[list[float], float, int]:
    """AdaGrad logistic regression. Returns (best weights, best val AUC, best epoch)."""
    w = [0.0] * N_WEIGHTS
    g2 = [1e-8] * N_WEIGHTS
    bias = 0.0
    bias_g2 = 1e-8
    order = list(range(len(feats)))
    rng = random.Random(seed)
    best = (list(w), -1.0, 0, 0.0)

    for epoch in range(1, epochs + 1):
        rng.shuffle(order)
        loss = 0.0
        for n in order:
            idx = feats[n]
            z = bias + sum(w[i] for i in idx)
            p = sigmoid(z)
            y = labels[n]
            loss -= math.log(max(p, 1e-12)) if y else math.log(max(1 - p, 1e-12))
            g = p - y
            bias_g2 += g * g
            bias -= lr * g / math.sqrt(bias_g2)
            for i in idx:
                gi = g + l2 * w[i]
                g2[i] += gi * gi
                w[i] -= lr * gi / math.sqrt(g2[i])
        val_scores = [bias + sum(w[i] for i in idx) for idx in val_feats]
        a = auc(val_scores, val_labels)
        print(f"  epoch {epoch:>2}  train logloss {loss / len(feats):.4f}   val AUC {a:.4f}")
        if a > best[1]:
            best = (list(w), a, epoch, bias)

    print(f"  best epoch {best[2]} (val AUC {best[1]:.4f})")
    return best[0] + [best[3]], best[1], best[2]


def nearest_in_train(seq: str, train_seqs: set[str], by_len: list[str]) -> tuple[int, str]:
    """Minimum Hamming distance to any training CDRH3, and one witness."""
    if seq in train_seqs:
        return 0, seq
    best, witness = LEN + 1, ""
    for other in by_len:
        d = 0
        for a, b in zip(seq, other):
            if a != b:
                d += 1
                if d >= best:
                    break
        if d < best:
            best, witness = d, other
            if best == 1:
                break
    return best, witness


def extract_cdrh3(heavy: str) -> str | None:
    """The designed 10-residue window: IMGT CDRH3 is SR + window + Y."""
    m = re.search(r"YYCSR(.{10})Y", heavy)
    return m.group(1) if m else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--epochs", type=int, default=6)
    ap.add_argument("--lr", type=float, default=0.3)
    ap.add_argument("--l2", type=float, default=1e-6)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--score", type=Path, default=Path("runs/team15-claude-n100-seed1/shortlist.csv"))
    ap.add_argument("--unary-only", action="store_true", help="ablation: drop the pairwise features")
    args = ap.parse_args()

    for name in ("train.csv", "val.csv", "test.csv"):
        if not (args.data_dir / name).exists():
            raise SystemExit(
                f"missing {args.data_dir / name}\n"
                "extract the authors' splits from zenodo-10831512-trastuzumab-ngs/data.zip:\n"
                "  data/affinity_data/her2/her2_aff_large/processed/remove_overlap/"
                "clonotype_split/0.7_0.15_0.15/{train,val,test}.csv"
            )

    print("=" * 78)
    print("DATA  (authors' clonotype-aware splits, used as shipped)")
    print("=" * 78)
    splits = {}
    for name in ("train", "val", "test"):
        seqs, labels, clusters, dropped = load_split(args.data_dir / f"{name}.csv")
        splits[name] = (seqs, labels, clusters)
        n_pos = sum(labels)
        print(
            f"  {name:<6} n={len(seqs):>7}  positives={n_pos:>6} ({n_pos / len(seqs):.1%})  "
            f"clusters={len({c for c in clusters if c}):>6}  dropped={dropped}"
        )
    tr_seqs, tr_labels, _ = splits["train"]
    va_seqs, va_labels, _ = splits["val"]
    te_seqs, te_labels, te_clusters = splits["test"]

    overlap = len(set(te_seqs) & set(tr_seqs))
    print(f"  exact CDRH3 overlap train n test: {overlap}")
    print()
    print("  where the measured mass sits, by substitutions from trastuzumab's window:")
    buckets: dict[int, list[int]] = defaultdict(list)
    for s, y in zip(tr_seqs, tr_labels):
        buckets[sum(1 for a, b in zip(s, WILD_TYPE) if a != b)].append(y)
    print(f"    {'subs':>4}{'train rows':>12}{'% of train':>12}{'high-binder rate':>18}")
    for d in sorted(buckets):
        ys = buckets[d]
        print(f"    {d:>4}{len(ys):>12}{len(ys) / len(tr_seqs):>11.1%}{sum(ys) / len(ys):>17.1%}")
    print("    (our accepted designs carry 5-9 substitutions -- check they are not off the end)")
    print()

    print("=" * 78)
    print("BASELINES")
    print("=" * 78)
    maj = 1 if sum(tr_labels) * 2 > len(tr_labels) else 0
    print(f"  majority class = {maj}; accuracy on test = {sum(1 for y in te_labels if y == maj) / len(te_labels):.1%}, AUC = 0.500")
    # Is the task just "count the mutations"? If distance from wild-type alone
    # scores well, a model that looks clever may be riding that one number --
    # and our designs all sit at 5-9 substitutions, off the end of where most of
    # the measured mass lies, so we would be extrapolating on exactly that axis.
    dist_scores = [-sum(1 for a, b in zip(s, WILD_TYPE) if a != b) for s in te_seqs]
    a_dist = auc(dist_scores, te_labels)
    lo, hi = auc_ci(a_dist, sum(te_labels), len(te_labels) - sum(te_labels))
    print(f"  edit distance from wild-type alone  test AUC {a_dist:.4f}  95% CI [{lo:.4f}, {hi:.4f}]")

    pwm = train_pwm(tr_seqs, tr_labels)
    pwm_test = [score_pwm(s, pwm) for s in te_seqs]
    a_pwm = auc(pwm_test, te_labels)
    lo, hi = auc_ci(a_pwm, sum(te_labels), len(te_labels) - sum(te_labels))
    print(f"  pwm (independent positions)  test AUC {a_pwm:.4f}  95% CI [{lo:.4f}, {hi:.4f}]")
    print()

    print("=" * 78)
    print("LOGISTIC REGRESSION" + ("  (unary only -- ablation)" if args.unary_only else "  (unary + pairwise)"))
    print("=" * 78)
    def feats_of(seqs: list[str]) -> list[list[int]]:
        out = []
        for s in seqs:
            f = featurise(s)
            out.append(f[:LEN] if args.unary_only else f)
        return out

    tr_feats = feats_of(tr_seqs)
    va_feats = feats_of(va_seqs)
    w_and_bias, best_val, _ = train_logreg(
        tr_feats, tr_labels, va_feats, va_labels, args.epochs, args.lr, args.l2, args.seed
    )
    w, bias = w_and_bias[:-1], w_and_bias[-1]

    te_feats = feats_of(te_seqs)
    te_scores = [bias + sum(w[i] for i in idx) for idx in te_feats]
    a_lr = auc(te_scores, te_labels)
    n_pos = sum(te_labels)
    lo, hi = auc_ci(a_lr, n_pos, len(te_labels) - n_pos)
    print()
    print(f"  TEST AUC {a_lr:.4f}  95% CI [{lo:.4f}, {hi:.4f}]   (scored once, after epoch choice on val)")

    preds = [1 if sigmoid(s) >= 0.5 else 0 for s in te_scores]
    tp = sum(1 for p, y in zip(preds, te_labels) if p == 1 and y == 1)
    fp = sum(1 for p, y in zip(preds, te_labels) if p == 1 and y == 0)
    fn = sum(1 for p, y in zip(preds, te_labels) if p == 0 and y == 1)
    acc = sum(1 for p, y in zip(preds, te_labels) if p == y) / len(te_labels)
    prec = tp / (tp + fp) if tp + fp else float("nan")
    rec = tp / (tp + fn) if tp + fn else float("nan")
    print(f"  at threshold 0.5: accuracy {acc:.1%}  precision {prec:.1%}  recall {rec:.1%}")

    # Per-cluster macro AUC: the test split holds few clusters, so the row count
    # overstates the effective sample size.
    by_cluster: dict[str, tuple[list[float], list[int]]] = defaultdict(lambda: ([], []))
    for s, y, c in zip(te_scores, te_labels, te_clusters):
        by_cluster[c][0].append(s)
        by_cluster[c][1].append(y)
    usable = [(sc, ys) for sc, ys in by_cluster.values() if 0 < sum(ys) < len(ys)]
    if usable:
        macro = sum(auc(sc, ys) for sc, ys in usable) / len(usable)
        print(
            f"  macro AUC over {len(usable)} test clusters with both classes: {macro:.4f} "
            f"(of {len(by_cluster)} clusters)"
        )
    print()

    if not args.score or not args.score.exists():
        print(f"(no design table at {args.score}; skipping design scoring)")
        return

    print("=" * 78)
    print(f"SCORING DESIGNS  {args.score}")
    print("=" * 78)
    train_set = set(tr_seqs)
    train_list = list(train_set)
    all_measured = train_set | set(va_seqs) | set(te_seqs)

    rows = []
    with args.score.open() as fh:
        for row in csv.DictReader(fh):
            heavy = row.get("heavy_sequence", "")
            cdrh3 = extract_cdrh3(heavy)
            if not cdrh3:
                continue
            rows.append((row.get("sequence_id", "?"), cdrh3))
    rows.append(("control_wildtype", WILD_TYPE))

    scored = []
    for name, cdrh3 in rows:
        f = featurise(cdrh3)
        if f is None:
            continue
        f = f[:LEN] if args.unary_only else f
        p = sigmoid(bias + sum(w[i] for i in f))
        d, witness = nearest_in_train(cdrh3, train_set, train_list)
        scored.append((p, name, cdrh3, d, witness, cdrh3 in all_measured))

    scored.sort(reverse=True)
    print(f"  {'design':<22}{'CDRH3':<12}{'P(high)':>9}{'nn dist':>9}  nearest measured CDRH3")
    for p, name, cdrh3, d, witness, measured in scored:
        flag = "  <- IN training data" if measured else ""
        print(f"  {name:<22}{cdrh3:<12}{p:>9.3f}{d:>9}  {witness}{flag}")

    print()
    print("Reading this honestly:")
    print(f"  - Test AUC {a_lr:.4f}. The pwm baseline gets {a_pwm:.4f}; if the gap is small,")
    print("    the interactions are not earning their keep and the simpler model should be used.")
    print("  - 'nn dist' is the Hamming distance to the nearest TRAINING CDRH3. A design at")
    print("    distance 1-2 sits inside measured territory, so its prediction is interpolation;")
    print("    at larger distances the model is extrapolating and the number is weaker evidence.")
    print("  - The label is a FACS binding class, not an affinity. It cannot rank designs by KD,")
    print("    and the competition scores affinity by SPR.")
    print("  - Trastuzumab's own window scores as a control: a model that cannot rate wild-type")
    print("    highly is not measuring what we want.")


if __name__ == "__main__":
    main()
