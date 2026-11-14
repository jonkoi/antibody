#!/usr/bin/env python3
"""Feature encodings for 10-residue CDRH3 variants.

Dependencies: numpy, scipy.

TWO CODINGS

``onehot``  -- the textbook (10, 20) indicator, 200 features.  What the
              published CNN takes.

``mutation``-- WT-centred.  A feature fires only where the sequence differs
              from trastuzumab:

                order 1: (position, substituted residue)   10 x 19 = 190
                order 2: (pos p mut, pos q mut) for p<q     45 x 361 = 16,245
                order 3: triples                           120 x 6,859 = 823,080

              The span of orders 1+2 is the same as one-hot singles + pairs, but
              the *regularisation* is not: L2 on the mutation coding shrinks
              towards "behaves like wild-type", which is the correct prior for a
              deep mutational scan.  A sequence at edit distance d activates
              exactly d + C(d,2) + C(d,3) features, so the matrix is very
              sparse (<= 55 non-zeros per row at order 2).

              This is also the natural basis for the problem: every legal
              candidate is a 3-mutant, so its score decomposes into 3 single
              effects, 3 pairwise couplings and 1 three-way coupling.

``physchem`` -- 5 Georgiev/VHSE-style scales per position (50 features) plus
              whole-loop aggregates, for the tree models.
"""
from __future__ import annotations

import itertools

import numpy as np
from scipy import sparse

from binding_data import ALPHABET, AA_INDEX, LENGTH, WT_CDRH3

WT_CODES = np.array([AA_INDEX[c] for c in WT_CDRH3], dtype=np.int64)
PAIRS = list(itertools.combinations(range(LENGTH), 2))
TRIPLES = list(itertools.combinations(range(LENGTH), 3))

# Five orthogonal physicochemical scales (Georgiev PC1-PC5, z-scored).
_SCALES = {
    "A": (0.57, 3.37, -3.66, 2.34, -1.07), "C": (2.66, -1.52, -3.29, -3.77, 2.96),
    "D": (-2.46, -0.66, -0.57, 0.14, 0.75), "E": (-3.08, 3.45, 0.05, 0.62, -0.49),
    "F": (3.12, 0.68, 2.40, -0.35, -0.88), "G": (0.15, -3.49, -2.97, 2.06, 0.70),
    "H": (-0.39, 1.00, -0.89, -2.25, 0.99), "I": (3.10, 0.37, 0.26, 1.04, -0.05),
    "K": (-3.89, 1.47, 1.95, 1.17, 0.53), "L": (2.72, 1.88, 1.92, 5.33, 0.08),
    "M": (1.89, 3.88, -1.57, -3.58, -2.55), "N": (-2.02, -1.92, 0.04, -0.65, 1.61),
    "P": (-0.58, -4.33, -0.02, -0.21, -8.31), "Q": (-2.54, 1.82, -0.82, -1.85, 0.09),
    "R": (-2.80, 0.31, 2.84, 0.25, 0.20), "S": (-1.10, -2.05, -2.19, 1.36, 1.78),
    "T": (-0.65, -1.60, -1.39, 0.63, 1.35), "V": (2.64, 0.03, -0.67, 2.34, 0.64),
    "W": (1.89, -0.09, 4.21, -2.77, 0.72), "Y": (0.79, -2.62, 4.11, -0.63, 1.47),
}
SCALE_MAT = np.array([_SCALES[a] for a in ALPHABET], dtype=np.float32)
_HYDRO = {"A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2,
          "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5,
          "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3}
_CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}
_VOL = {"A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9, "G": 60.1,
        "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7, "M": 162.9, "N": 114.1,
        "P": 112.7, "Q": 143.8, "R": 173.4, "S": 89.0, "T": 116.1, "V": 140.0,
        "W": 227.8, "Y": 193.6}
HYDRO = np.array([_HYDRO[a] for a in ALPHABET], dtype=np.float32)
CHARGE = np.array([_CHARGE.get(a, 0.0) for a in ALPHABET], dtype=np.float32)
VOLUME = np.array([_VOL[a] for a in ALPHABET], dtype=np.float32)


def onehot(codes: np.ndarray) -> np.ndarray:
    n = codes.shape[0]
    out = np.zeros((n, LENGTH * 20), dtype=np.float32)
    cols = codes.astype(np.int64) + (np.arange(LENGTH) * 20)[None, :]
    out[np.arange(n)[:, None], cols] = 1.0
    return out


def onehot_3d(codes: np.ndarray) -> np.ndarray:
    """(n, 10, 20) float32 -- the CNN's input shape."""
    return onehot(codes).reshape(-1, LENGTH, 20)


def mutation_dim(order: int) -> int:
    d = LENGTH * 19
    if order >= 2:
        d += len(PAIRS) * 19 * 19
    if order >= 3:
        d += len(TRIPLES) * 19 * 19 * 19
    return d


def mutation_sparse(codes: np.ndarray, order: int = 2) -> sparse.csr_matrix:
    """WT-centred interaction features as a CSR matrix.

    Single-mutation column for (p, a) is p*19 + rank(a) where rank skips the
    wild-type residue at p.  Pair and triple blocks index the Cartesian product
    of those within-position ranks.
    """
    codes = codes.astype(np.int64)
    n = codes.shape[0]
    # within-position rank of the substituted residue, 0..18; -1 if wild-type
    shift = (codes > WT_CODES[None, :]).astype(np.int64)
    is_mut = codes != WT_CODES[None, :]
    rank = np.where(is_mut, codes - shift, -1)

    off_pair = LENGTH * 19
    off_tri = off_pair + len(PAIRS) * 361

    indptr = np.zeros(n + 1, dtype=np.int64)
    cols_all: list[np.ndarray] = []
    rows_nnz = np.zeros(n, dtype=np.int64)
    mut_positions = [np.flatnonzero(is_mut[i]) for i in range(n)]

    for i in range(n):
        mp = mut_positions[i]
        r = rank[i, mp]
        cols = [mp * 19 + r]
        if order >= 2 and mp.size >= 2:
            pi, qi = np.triu_indices(mp.size, k=1)
            p, q = mp[pi], mp[qi]
            pair_block = np.array([PAIR_INDEX[(int(a), int(b))] for a, b in zip(p, q)])
            cols.append(off_pair + pair_block * 361 + r[pi] * 19 + r[qi])
        if order >= 3 and mp.size >= 3:
            combos = np.array(list(itertools.combinations(range(mp.size), 3)))
            a, b, c = combos[:, 0], combos[:, 1], combos[:, 2]
            tri_block = np.array([TRIPLE_INDEX[(int(mp[x]), int(mp[y]), int(mp[z]))]
                                  for x, y, z in combos])
            cols.append(off_tri + tri_block * 6859 + r[a] * 361 + r[b] * 19 + r[c])
        merged = np.concatenate(cols) if cols else np.empty(0, dtype=np.int64)
        cols_all.append(merged)
        rows_nnz[i] = merged.size

    indptr[1:] = np.cumsum(rows_nnz)
    indices = np.concatenate(cols_all) if cols_all else np.empty(0, dtype=np.int64)
    data = np.ones(indices.size, dtype=np.float32)
    return sparse.csr_matrix((data, indices, indptr),
                             shape=(n, mutation_dim(order)))


PAIR_INDEX = {pq: i for i, pq in enumerate(PAIRS)}
TRIPLE_INDEX = {pqr: i for i, pqr in enumerate(TRIPLES)}


def mutation_sparse_fast(codes: np.ndarray, order: int = 2) -> sparse.csr_matrix:
    """Vectorised equivalent of ``mutation_sparse`` (identical output).

    Builds the COO triplets with numpy instead of a Python loop over rows, which
    matters at 524k rows.  Each (row, position-pair) slot is emitted only when
    both positions are mutated, so the non-zero pattern is data dependent; the
    trick is to emit every candidate slot and then mask.
    """
    codes = codes.astype(np.int64)
    n = codes.shape[0]
    shift = (codes > WT_CODES[None, :]).astype(np.int64)
    rank = codes - shift
    is_mut = codes != WT_CODES[None, :]

    off_pair = LENGTH * 19
    off_tri = off_pair + len(PAIRS) * 361

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []

    r_idx, p_idx = np.nonzero(is_mut)
    rows.append(r_idx)
    cols.append(p_idx * 19 + rank[r_idx, p_idx])

    if order >= 2:
        for blk, (p, q) in enumerate(PAIRS):
            m = is_mut[:, p] & is_mut[:, q]
            if not m.any():
                continue
            ri = np.flatnonzero(m)
            rows.append(ri)
            cols.append(off_pair + blk * 361 + rank[ri, p] * 19 + rank[ri, q])

    if order >= 3:
        for blk, (p, q, s) in enumerate(TRIPLES):
            m = is_mut[:, p] & is_mut[:, q] & is_mut[:, s]
            if not m.any():
                continue
            ri = np.flatnonzero(m)
            rows.append(ri)
            cols.append(off_tri + blk * 6859
                        + rank[ri, p] * 361 + rank[ri, q] * 19 + rank[ri, s])

    row = np.concatenate(rows)
    col = np.concatenate(cols)
    data = np.ones(row.size, dtype=np.float32)
    mat = sparse.coo_matrix((data, (row, col)), shape=(n, mutation_dim(order)))
    return mat.tocsr()


def onehot_interaction_dim(order: int) -> int:
    d = LENGTH * 20
    if order >= 2:
        d += len(PAIRS) * 400
    if order >= 3:
        d += len(TRIPLES) * 8000
    return d


def onehot_interaction_sparse(codes: np.ndarray, order: int = 2) -> sparse.csr_matrix:
    """Plain one-hot singles plus residue-pair (and triple) co-occurrence.

    Not WT-centred: the column for positions (p, q) is indexed by the actual
    residue pair, so every row fires exactly one column per position pair
    regardless of how far the sequence is from trastuzumab.

    This spans the same space as ``mutation_sparse_fast`` at the same order, but
    it generalises better in practice -- L2 here shrinks towards "no preference
    between residues", whereas on the WT-centred coding it shrinks towards
    "behaves exactly like trastuzumab", which turns out to be the wrong prior
    once a sequence carries several substitutions.  Measured, not assumed: see
    the model comparison in train_binding_model.py, where this coding wins on
    distance-3 AUC and on SPR transfer.
    """
    codes = codes.astype(np.int64)
    n = codes.shape[0]
    r = np.arange(n)
    rows, cols = [], []
    for p in range(LENGTH):
        rows.append(r)
        cols.append(p * 20 + codes[:, p])
    off = LENGTH * 20
    if order >= 2:
        for blk, (p, q) in enumerate(PAIRS):
            rows.append(r)
            cols.append(off + blk * 400 + codes[:, p] * 20 + codes[:, q])
        off += len(PAIRS) * 400
    if order >= 3:
        for blk, (p, q, s) in enumerate(TRIPLES):
            rows.append(r)
            cols.append(off + blk * 8000
                        + codes[:, p] * 400 + codes[:, q] * 20 + codes[:, s])
    row = np.concatenate(rows)
    col = np.concatenate(cols)
    return sparse.coo_matrix((np.ones(row.size, np.float32), (row, col)),
                             shape=(n, onehot_interaction_dim(order))).tocsr()


def physchem(codes: np.ndarray) -> np.ndarray:
    """Per-position scales plus whole-loop aggregates and WT-contrast terms."""
    codes = codes.astype(np.int64)
    n = codes.shape[0]
    pcs = SCALE_MAT[codes].reshape(n, LENGTH * 5)                 # 50
    hyd = HYDRO[codes]                                            # (n,10)
    chg = CHARGE[codes]
    vol = VOLUME[codes]
    is_mut = (codes != WT_CODES[None, :]).astype(np.float32)
    aa_counts = np.zeros((n, 20), dtype=np.float32)
    for a in range(20):
        aa_counts[:, a] = (codes == a).sum(axis=1)
    agg = np.stack([
        hyd.sum(1), hyd.mean(1), chg.sum(1), np.abs(chg).sum(1),
        vol.sum(1), vol.std(1), is_mut.sum(1),
        (VOLUME[codes] - VOLUME[WT_CODES][None, :]).sum(1),
        (HYDRO[codes] - HYDRO[WT_CODES][None, :]).sum(1),
        (CHARGE[codes] - CHARGE[WT_CODES][None, :]).sum(1),
    ], axis=1).astype(np.float32)
    return np.concatenate([pcs, hyd, chg, vol, is_mut, aa_counts, agg], axis=1)


def _self_check() -> None:
    rng = np.random.default_rng(1)
    codes = rng.integers(0, 20, size=(200, LENGTH))
    codes[0] = WT_CODES
    for order in (1, 2, 3):
        a = mutation_sparse(codes[:60], order)
        b = mutation_sparse_fast(codes[:60], order)
        assert (a != b).nnz == 0, f"order {order} mismatch"
    x = mutation_sparse_fast(codes, 2)
    d = (codes != WT_CODES[None, :]).sum(1)
    expect = d + d * (d - 1) // 2
    assert np.array_equal(np.diff(x.indptr), expect), "order-2 nnz wrong"
    x3 = mutation_sparse_fast(codes, 3)
    expect3 = expect + d * (d - 1) * (d - 2) // 6
    assert np.array_equal(np.diff(x3.indptr), expect3), "order-3 nnz wrong"
    assert x.getrow(0).nnz == 0, "wild-type row should be all zeros"
    assert onehot_3d(codes).sum() == codes.shape[0] * LENGTH
    for order, per_row in ((1, 10), (2, 55), (3, 175)):
        oi = onehot_interaction_sparse(codes, order)
        assert oi.shape == (codes.shape[0], onehot_interaction_dim(order))
        assert np.all(np.diff(oi.indptr) == per_row), f"one-hot order {order} nnz"
    print(f"features self-check passed\n"
          f"  WT-centred : order-2 {mutation_dim(2):,}  order-3 {mutation_dim(3):,}\n"
          f"  one-hot    : order-2 {onehot_interaction_dim(2):,}  "
          f"order-3 {onehot_interaction_dim(3):,}\n"
          f"  physchem   : {physchem(codes[:2]).shape[1]}")


if __name__ == "__main__":
    _self_check()
