#!/usr/bin/env python3
"""ESM-2 embeddings and zero-shot likelihood scores for trastuzumab VH variants.

Runs on the hackathon GPU box (H100).  Dependencies: torch (cu124), transformers.
Reads a text file of 10-residue CDRH3 sequences, one per line, and writes an
.npz holding, per sequence:

    res     (n, 10, 1280) float16 -- per-residue representation at the ten CDRH3
                                     positions of the full 120-residue VH
    cdr     (n, 1280)      float16 -- mean over those ten positions
    vh      (n, 1280)      float16 -- mean over the whole VH
    mm      (n,)           float32 -- masked-marginal pseudo-log-likelihood of
                                     the CDRH3, summed over its ten positions

Only the CDRH3 varies; the 108 framework residues and the whole light chain are
constant, so a whole-VH mean pool is dominated by constant signal.  The
per-residue slice at the ten varying positions is the informative part; the
pooled vectors are kept because they are small and the SPR regression has only
180 training points, where 12,800 features is far too many.

The masked-marginal score is the standard ESM zero-shot fitness proxy: mask each
CDRH3 position in turn and read off the log-probability the model assigns to the
residue actually there.  It measures "does this look like a real antibody loop",
which is a different question from "does it bind HER2", and is included as an
orthogonal feature rather than as a ranker.

Usage:
    python esm_embed.py --seqs cdrh3.txt --out emb.npz [--model facebook/esm2_t33_650M_UR50D]
                        [--batch 64] [--no-res] [--no-mm]
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import torch

VH_PREFIX = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
             "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSR")
VH_SUFFIX = "YWGQGTLVTVSS"
CDR_START = len(VH_PREFIX)          # 98
CDR_LEN = 10


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seqs", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--model", default="facebook/esm2_t33_650M_UR50D")
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--no-res", action="store_true", help="skip per-residue output")
    ap.add_argument("--no-mm", action="store_true", help="skip masked marginals")
    args = ap.parse_args()

    from transformers import AutoTokenizer, EsmForMaskedLM

    cdrs = [l.strip().upper() for l in args.seqs.read_text().split() if l.strip()]
    assert all(len(c) == CDR_LEN for c in cdrs), "expected 10-residue CDRH3s"
    seqs = [VH_PREFIX + c + VH_SUFFIX for c in cdrs]
    print(f"{len(seqs):,} sequences, length {len(seqs[0])}", flush=True)

    tok = AutoTokenizer.from_pretrained(args.model)
    model = EsmForMaskedLM.from_pretrained(args.model, torch_dtype=torch.float16)
    model = model.eval().cuda()
    dim = model.config.hidden_size
    mask_id = tok.mask_token_id

    # token index of VH residue i is i + 1 (the <cls> token sits at 0)
    cdr_tok = np.arange(CDR_START, CDR_START + CDR_LEN) + 1

    n = len(seqs)
    res = np.zeros((n, CDR_LEN, dim), np.float16) if not args.no_res else None
    cdr_mean = np.zeros((n, dim), np.float16)
    vh_mean = np.zeros((n, dim), np.float16)
    mm = np.zeros(n, np.float32)

    t0 = time.time()
    for start in range(0, n, args.batch):
        batch = seqs[start:start + args.batch]
        enc = tok(batch, return_tensors="pt", padding=True)
        ids = enc["input_ids"].cuda()
        attn = enc["attention_mask"].cuda()
        with torch.no_grad():
            out = model(input_ids=ids, attention_mask=attn, output_hidden_states=True)
            h = out.hidden_states[-1]                       # (b, L, dim)
            sl = h[:, cdr_tok, :]
            if res is not None:
                res[start:start + len(batch)] = sl.float().cpu().numpy().astype(np.float16)
            cdr_mean[start:start + len(batch)] = sl.mean(1).float().cpu().numpy().astype(np.float16)
            real = attn.unsqueeze(-1).to(h.dtype)
            vh_mean[start:start + len(batch)] = (
                (h * real).sum(1) / real.sum(1)).float().cpu().numpy().astype(np.float16)

            if not args.no_mm:
                total = torch.zeros(len(batch), device="cuda")
                for k, ti in enumerate(cdr_tok):
                    masked = ids.clone()
                    true = masked[:, ti].clone()
                    masked[:, ti] = mask_id
                    logits = model(input_ids=masked, attention_mask=attn).logits[:, ti, :]
                    lp = torch.log_softmax(logits.float(), dim=-1)
                    total += lp.gather(1, true.unsqueeze(1)).squeeze(1)
                mm[start:start + len(batch)] = total.cpu().numpy()

        if (start // args.batch) % 50 == 0:
            done = start + len(batch)
            rate = done / max(time.time() - t0, 1e-9)
            print(f"  {done:,}/{n:,}  {rate:.0f} seq/s  "
                  f"eta {(n - done) / max(rate, 1e-9) / 60:.1f} min", flush=True)

    payload = {"cdrh3": np.array(cdrs), "cdr": cdr_mean, "vh": vh_mean, "mm": mm}
    if res is not None:
        payload["res"] = res
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.out, **payload)
    print(f"wrote {args.out} in {time.time() - t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
