#!/usr/bin/env python3
"""Hamming-neighbourhood tools for leakage-matched evaluation.

Dependencies: numpy.

WHY THIS EXISTS
The 51,417 legal candidates clear a novelty filter: each is at least 3
substitutions from trastuzumab *and* from every one of the 525,063 CDRH3s in
the provided library.  So at prediction time the model never sees a sequence
within Hamming distance 2 of anything it was trained on.

A random CV split does not reproduce that.  Most held-out 10-mers have a
Hamming-1 or Hamming-2 neighbour in the training fold, and a model with 18,200
pairwise features can exploit it.  The AUC that comes out is therefore an
optimistic estimate of what we will get on the candidates.

``min_dist_ge3_mask`` marks the test rows that are at least 3 substitutions from
every training row.  Evaluating on those -- the "novel shell" -- matches
deployment exactly.  It is the number we trust.

Single linkage at Hamming distance 1 was tried first and is unusable: the
library is dense enough that the Hamming-1 graph percolates into one component
holding 84.5% of the data, so grouped splits degenerate.
"""
from __future__ import annotations

import numpy as np

LENGTH = 10
BASE = 20
_POW = (BASE ** np.arange(LENGTH - 1, -1, -1)).astype(np.int64)


def pack(codes: np.ndarray) -> np.ndarray:
    """(n, 10) int codes -> (n,) int64 base-20 key."""
    return codes.astype(np.int64) @ _POW


def _neighbours_le2(keys: np.ndarray) -> np.ndarray:
    """All distance-1 and distance-2 variants of each key.

    (n,) -> (n, 10*19 + 45*19*19) == (n, 16435).  Changing position p from
    residue r to r' shifts the key by (r' - r) * 20**(9-p); the two-position
    case is the sum of two such shifts.  Working in key space avoids
    materialising sequences.
    """
    n = keys.size
    digits = np.empty((n, LENGTH), dtype=np.int64)
    rest = keys.copy()
    for p in range(LENGTH - 1, -1, -1):
        digits[:, p] = rest % BASE
        rest //= BASE

    # per-position deltas to every other residue: (n, 10, 19)
    alt = np.arange(BASE, dtype=np.int64)[None, None, :]              # (1,1,20)
    delta = (alt - digits[:, :, None]) * _POW[None, :, None]          # (n,10,20)
    mask = alt != digits[:, :, None]
    d1 = delta[np.broadcast_to(mask, delta.shape)].reshape(n, LENGTH, BASE - 1)

    out = [keys[:, None] + d1.reshape(n, -1)]
    for p in range(LENGTH):
        for q in range(p + 1, LENGTH):
            pair = (keys[:, None, None]
                    + d1[:, p, :].reshape(n, BASE - 1, 1)
                    + d1[:, q, :].reshape(n, 1, BASE - 1))
            out.append(pair.reshape(n, -1))
    return np.concatenate(out, axis=1)


def min_dist_ge3_mask(test_codes: np.ndarray, train_codes: np.ndarray,
                      chunk: int = 256) -> np.ndarray:
    """True where the test row is >= 3 substitutions from every training row.

    Exhaustive and exact: enumerates the 16,435 sequences within distance 2 of
    each test row and looks every one up in the (sorted) training key set.
    Also excludes exact matches.
    """
    train_keys = np.unique(pack(train_codes))
    test_keys = pack(test_codes)
    keep = np.ones(test_keys.size, dtype=bool)
    for start in range(0, test_keys.size, chunk):
        block = test_keys[start:start + chunk]
        hits = np.isin(block, train_keys)
        nb = _neighbours_le2(block)
        hits |= np.isin(nb, train_keys).any(axis=1)
        keep[start:start + chunk] = ~hits
    return keep


def deployment_split(test_codes: np.ndarray, all_codes: np.ndarray,
                     chunk: int = 256) -> np.ndarray:
    """Training mask that reproduces the candidates' novelty condition.

    Returns a boolean mask over ``all_codes``, True where that row is at least
    3 substitutions from *every* row of ``test_codes``.  Training on the masked
    rows and testing on ``test_codes`` is the only protocol in this repository
    that matches deployment: the 51,417 legal candidates are >= 3 substitutions
    from all 525,063 library sequences, so the model will never have a close
    relative of a candidate in its training set.

    Built in the cheap direction: enumerate the distance-<= 2 neighbourhood of
    the (small) test set once, then test the (large) corpus for membership.
    """
    test_keys = pack(test_codes)
    blocks = [np.unique(test_keys)]
    for start in range(0, test_keys.size, chunk):
        blocks.append(np.unique(_neighbours_le2(test_keys[start:start + chunk]).ravel()))
    forbidden = np.unique(np.concatenate(blocks))
    return ~np.isin(pack(all_codes), forbidden)


def _self_check() -> None:
    rng = np.random.default_rng(0)
    codes = rng.integers(0, 20, size=(5, LENGTH))
    nb = _neighbours_le2(pack(codes))
    assert nb.shape == (5, 10 * 19 + 45 * 19 * 19), nb.shape
    # decode a sample of neighbours and confirm the distance is 1 or 2
    for i in range(5):
        for k in rng.integers(0, nb.shape[1], size=50):
            key = int(nb[i, k])
            digits = []
            for _ in range(LENGTH):
                digits.append(key % BASE)
                key //= BASE
            digits = np.array(digits[::-1])
            d = int((digits != codes[i]).sum())
            assert d in (1, 2), d
    assert len(np.unique(nb[0])) == nb.shape[1], "neighbour keys should be distinct"
    print("neighbourhood self-check passed")


if __name__ == "__main__":
    _self_check()
