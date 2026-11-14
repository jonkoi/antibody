#!/usr/bin/env python3
"""Score trastuzumab CDRH3 variants against HER2 with JapanFold's AF2-IG model.

AF2-IG (AlphaFold2 initial-guess) re-predicts a designed complex from the
coordinates you give it and reports ipTM / interface pAE for the target-binder
interface.  It is deterministic (no MSA, no sampling, no seed), which is why we
use it for ranking instead of a replica-averaged fold.

THE INPUT FORMAT THAT WORKS (discovered empirically, 2026-11-14)
----------------------------------------------------------------
`af2ig` takes exactly ONE target chain and ONE binder chain, and
`target.chain` must be a SINGLE CHARACTER -- "A,C,D" and ["A","C","D"] are both
refused with
    400 "target.chain must be a single-character chain id"
Any chain in the pasted structure that is neither the target chain nor the
binder chain is silently DROPPED (the job's `residues` field proves it: pasting
all four chains of 1N8Z with chain: C gave residues=212 = 92 + 120).

So a two-sided complex has to be merged into one target chain.  The framing
that works, benchmarked below, is

    target chain 'T' = antibody light chain (A, numbered 1-107)
                     + HER2 489-580 (C, original numbering)
                     + HER2 591-607 (D, original numbering)
    binder chain 'B' = heavy chain, 120 aa, with the designed CDRH3 at 99-108

The large gaps in the residue numbering (107 -> 489, 580 -> 591) are what tell
AlphaFold these are separate segments, so no spurious peptide bonds are formed.

Measured on a panel of SPR-labelled CDRH3 variants (interface pAE, lower is
better; wild type = WGGDGFYAMD):

    target = HER2 489-580 only          ipTM 0.147  ipae 23.8   useless
    target = HER2 489-580 + 591-607     ipTM 0.251  ipae 21.8   weak ordering
    target = light chain + both HER2    ipTM 0.863  ipae 12.1   clean separation

Dropping the light chain destroys the signal: trastuzumab's epitope spans both
HER2 fragments, and AF2 needs VL present to dock VH confidently.

Only the backbone matters.  N/CA/C/O-only and all-heavy-atom structures give
the same answer to within 0.001 ipTM, so we paste backbone only and halve the
request body.

API limits that shape the batching
----------------------------------
  * request body <= 8 MB, and gzipped bodies are rejected  -> ~60 targets/job
    at 114 kB of PDB text per target
  * max_residues 1024 for af2ig (we use 336)
  * one account is held on two *waits* at a time, so we poll one job per call

Usage
-----
    export JAPANFOLD_API_KEY=...            # never hardcoded here
    python3 pipeline/score_af2ig.py panel                 # validation panel
    python3 pipeline/score_af2ig.py spr   --limit 100      # SPR-labelled set
    python3 pipeline/score_af2ig.py rank  --limit 3000     # ranked candidates
    python3 pipeline/score_af2ig.py pae   --top 200        # keep PAE matrices

Every mode appends to data/af2ig_scores.csv and skips design ids already in it,
so an interrupted run is resumed by re-running the same command.
"""

import argparse
import csv
import hashlib
import json
import os
import sys
import time
import uuid
from collections import OrderedDict

import requests

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIF = os.path.join(REPO, "inputs", "12-her2", "her2_fv_1n8z.cif")
SPR = os.path.join(REPO, "data", "her2-trastuzumab", "trastuzumab_her2_spr_all_768.csv")
RANKED = os.path.join(REPO, "data", "legal_3mutants_ranked.csv")
PROD = os.path.join(REPO, "runs", "team15-claude-n100-seed1", "designs_summary.csv")
OUT = os.path.join(REPO, "data", "af2ig_scores.csv")
PAE_DIR = os.path.join(REPO, "data", "af2ig_pae")

WT_CDRH3 = "WGGDGFYAMD"
CDRH3_START = 98          # 0-based index of heavy-chain position 99
CDRH3_END = 108

FIELDS = ["design_id", "cdrh3", "source", "label", "kd_nM", "prior_score",
          "iptm", "interface_pae", "ipae", "plddt", "ptm",
          "chain_ptm_target", "chain_ptm_binder", "n_residues", "binder_residues",
          "status", "job_id"]

AA3to1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
          "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
          "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
          "TYR": "Y", "VAL": "V"}
BACKBONE = ("N", "CA", "C", "O")


# --------------------------------------------------------------------------- #
# structure handling
# --------------------------------------------------------------------------- #
def parse_cif(path=CIF):
    """Minimal mmCIF atom_site reader -> {chain: {(seq_id, comp): [atoms]}}."""
    lines = open(path).read().splitlines()
    cols, start = [], None
    for i, line in enumerate(lines):
        if line.startswith("_atom_site."):
            cols.append(line.strip()[len("_atom_site."):])
            start = i
    if start is None:
        raise SystemExit("no _atom_site loop in %s" % path)
    idx = {c: k for k, c in enumerate(cols)}
    chains = OrderedDict()
    for line in lines[start + 1:]:
        if line.startswith("#") or not line.strip():
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        f = line.split()
        if f[idx["type_symbol"]] == "H":
            continue
        name = f[idx["label_atom_id"]]
        if name not in BACKBONE:
            continue
        key = (int(f[idx["label_seq_id"]]), f[idx["label_comp_id"]])
        xyz = tuple(float(f[idx["Cartn_" + a]]) for a in "xyz")
        chains.setdefault(f[idx["label_asym_id"]], OrderedDict()).setdefault(key, []).append(
            (name, xyz, f[idx["type_symbol"]]))
    return chains


def chain_seq(chains, c):
    return "".join(AA3to1.get(k[1], "X") for k in chains[c])


# The one framing that works: light chain + both HER2 fragments merged into a
# single target chain 'T', heavy chain kept as binder chain 'B'.
TARGET_GROUPS = [("T", [("A", 1), ("C", 489), ("D", 591)]),
                 ("B", [("B", 1)])]
TARGET_CHAIN = "T"
BINDER_CHAIN = "B"


def build_pdb(chains):
    out, ser = [], 1
    for newc, parts in TARGET_GROUPS:
        for src, start in parts:
            for i, (key, atoms) in enumerate(chains[src].items()):
                resi, comp = start + i, key[1]
                for name, (x, y, z), el in atoms:
                    out.append("ATOM  %5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
                               % (ser, name, comp, newc, resi, x, y, z, el.rjust(2)))
                    ser += 1
        out.append("TER")
    return "\n".join(out) + "\nEND\n"


def target_yaml(pdb, binder_seq):
    body = "".join("      " + l + "\n" for l in pdb.splitlines())
    return ("target:\n  structure: |\n%s  chain: %s\nbinder:\n  sequence: %s\n  chain: %s\n"
            % (body, TARGET_CHAIN, binder_seq, BINDER_CHAIN))


def heavy_with(cdrh3, wt_heavy):
    if len(cdrh3) != len(WT_CDRH3):
        raise ValueError("CDRH3 must be %d aa, got %r" % (len(WT_CDRH3), cdrh3))
    return wt_heavy[:CDRH3_START] + cdrh3 + wt_heavy[CDRH3_END:]


# --------------------------------------------------------------------------- #
# API client
# --------------------------------------------------------------------------- #
class JapanFold(object):
    MAX_BODY = 8 * 1024 * 1024

    def __init__(self):
        key = os.environ.get("JAPANFOLD_API_KEY")
        if not key:
            raise SystemExit("JAPANFOLD_API_KEY is not set; export it and re-run "
                             "(get one at https://japanfold.aiand.com/account/token)")
        self.base = os.environ.get("JAPANFOLD_BASE_URL",
                                   "https://api.japanfold.aiand.com").rstrip("/")
        self.s = requests.Session()
        self.s.headers.update({"Authorization": "Bearer " + key,
                               "X-JapanFold-Client": "skill"})

    def credits(self):
        return self.s.get(self.base + "/v1/credits", timeout=60).json()

    def submit(self, payload, idem=None):
        body = json.dumps(payload)
        if len(body) > self.MAX_BODY:
            raise ValueError("body %.1f MB exceeds the API's 8 MB limit" % (len(body) / 1e6))
        h = {"Content-Type": "application/json",
             "Idempotency-Key": idem or str(uuid.uuid4())}
        r = self.s.post(self.base + "/v1/predictions", data=body, headers=h, timeout=900)
        if r.status_code not in (200, 202):
            raise RuntimeError("submit %s: %s" % (r.status_code, r.text[:500]))
        return r.json()

    def job(self, jid, wait=0):
        p = {"wait": wait} if wait else None
        return self.s.get(self.base + "/v1/jobs/" + jid, params=p,
                          timeout=(wait or 0) + 90).json()

    def results(self, jid):
        return self.s.get(self.base + "/v1/jobs/" + jid + "/results", timeout=300).json()

    def usage(self, jid):
        return self.s.get(self.base + "/v1/jobs/" + jid + "/usage", timeout=60).json()

    def artifact(self, jid, path, dest):
        r = self.s.get("%s/v1/jobs/%s/artifacts/%s" % (self.base, jid, path),
                       timeout=600, allow_redirects=True)
        r.raise_for_status()
        with open(dest, "wb") as fh:
            fh.write(r.content)
        return dest


# --------------------------------------------------------------------------- #
# candidate sources
# --------------------------------------------------------------------------- #
def load_spr(limit=None, seed=1):
    """SPR-labelled 10mer CDRH3s on the wild-type framework (ground truth)."""
    import random
    wt_heavy = chain_seq(parse_cif(), "B")
    rows = [r for r in csv.DictReader(open(SPR))
            if len(r["cdrh3"]) == 10
            and r["vh_sequence"][:CDRH3_START] == wt_heavy[:CDRH3_START]
            and r["vh_sequence"][CDRH3_END:] == wt_heavy[CDRH3_END:]]
    seen, binders, nonb = set(), [], []
    for r in rows:
        if r["cdrh3"] in seen:
            continue
        seen.add(r["cdrh3"])
        if r["binding_result"] == "quantified_binding" and r["kd_nM"]:
            binders.append(r)
        elif r["binding_result"] == "no_reportable_kd":
            nonb.append(r)
    random.Random(seed).shuffle(binders)
    random.Random(seed).shuffle(nonb)
    half = (limit or 200) // 2
    out = []
    for r in binders[:half]:
        out.append(dict(design_id="spr_B_" + r["cdrh3"], cdrh3=r["cdrh3"],
                        source="spr", label="binder", kd_nM=r["kd_nM"], prior_score=""))
    for r in nonb[:half]:
        out.append(dict(design_id="spr_N_" + r["cdrh3"], cdrh3=r["cdrh3"],
                        source="spr", label="nonbinder", kd_nM="", prior_score=""))
    return out


def load_panel():
    """Wild type plus the sanity-check controls."""
    spec = [("wt", WT_CDRH3, "wild_type", "0.0285"),
            ("pos_WGSQGFYAYD", "WGSQGFYAYD", "binder", "0.0234"),
            ("pos_WGSDGWYAMD", "WGSDGWYAMD", "binder", "0.0608"),
            ("pos_WGGQGFYAMD", "WGGQGFYAMD", "binder", "0.09"),
            ("neg_WGGDGFEAYD", "WGGDGFEAYD", "nonbinder", ""),
            ("neg_WGGDGDRAID", "WGGDGDRAID", "nonbinder", ""),
            ("neg_WGGDDFNIMD", "WGGDDFNIMD", "nonbinder", ""),
            ("neg_WGGDGGWAED", "WGGDGGWAED", "nonbinder", ""),
            ("scram_reverse", "DMAYFGDGGW", "scramble", ""),
            ("scram_charged", "EKEKDRDKRE", "scramble", "")]
    return [dict(design_id=i, cdrh3=c, source="panel", label=l, kd_nM=k, prior_score="")
            for i, c, l, k in spec]


def load_ranked(limit=None, offset=0, spread=0, path=None):
    """Candidates from the prior-ranked CSV (already sorted by `score`, best first).

    `spread` takes that many of the `limit` from evenly spaced positions across
    the whole file instead of the head, so the prior-vs-AF2IG rank correlation
    is not measured on a truncated range.
    """
    all_rows = list(csv.DictReader(open(path or RANKED)))
    def mk(i, r):
        return dict(design_id="cand_%06d_%s" % (i, r["cdrh3"]), cdrh3=r["cdrh3"],
                    source="ranked", label="", kd_nM="", prior_score=r["score"])
    n_head = (limit - spread) if limit else len(all_rows)
    out = [mk(i, all_rows[i]) for i in range(offset, min(offset + n_head, len(all_rows)))]
    if spread:
        step = max(1, (len(all_rows) - offset) // spread)
        for i in range(offset, len(all_rows), step):
            out.append(mk(i, all_rows[i]))
    return out


def load_production():
    out = []
    if not os.path.exists(PROD):
        return out
    for r in csv.DictReader(open(PROD)):
        out.append(dict(design_id="prod_" + r["design_id"], cdrh3=r["cdrh3"],
                        source="production", label="", kd_nM="", prior_score=""))
    return out


# --------------------------------------------------------------------------- #
# scoring driver
# --------------------------------------------------------------------------- #
def already_scored():
    done = {}
    if os.path.exists(OUT):
        for r in csv.DictReader(open(OUT)):
            if r.get("status") == "ok":
                done[r["design_id"]] = r
    return done


def append_rows(rows):
    """Append under an advisory lock so several workers can share the CSV."""
    import fcntl
    new = not os.path.exists(OUT)
    with open(OUT, "a", newline="") as fh:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX)
        try:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            if new:
                w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in FIELDS})
            fh.flush()
        finally:
            fcntl.flock(fh.fileno(), fcntl.LOCK_UN)


def batches(items, pdb, per_job):
    """Split into jobs small enough for the 8 MB body limit."""
    per_target = len(target_yaml(pdb, "A" * 120)) + 120
    cap = max(1, int((JapanFold.MAX_BODY - 50000) / per_target))
    n = min(per_job, cap)
    for i in range(0, len(items), n):
        yield items[i:i + n]


def score(items, jf, pdb, wt_heavy, per_job=50, in_flight=3, write_pae=False,
          name_prefix="af2ig", quiet=False):
    """Submit batches, collect rows, append to the CSV as each job lands."""
    todo = list(batches(items, pdb, per_job))
    pending, done_rows, spent = [], [], 0.0
    bi = 0
    while todo or pending:
        while todo and len(pending) < in_flight:
            chunk = todo.pop(0)
            bi += 1
            targets = [{"name": d["design_id"],
                        "content": target_yaml(pdb, heavy_with(d["cdrh3"], wt_heavy))}
                       for d in chunk]
            params = {"write_pae": True} if write_pae else {}
            payload = {"model": "af2ig", "name": "%s_b%03d" % (name_prefix, bi),
                       "targets": targets}
            if params:
                payload["params"] = params
            # Idempotency key = hash of what is being submitted, so a retry after a
            # timeout returns the first job instead of billing twice, while two
            # workers that legitimately submit different batches never collide.
            idem = hashlib.sha256(
                ("af2ig|" + "|".join(sorted(d["design_id"] + ":" + d["cdrh3"] for d in chunk))
                 + "|pae=%d" % bool(write_pae)).encode()).hexdigest()
            job = jf.submit(payload, idem=idem)
            pending.append((job["id"], chunk))
            if not quiet:
                print("[submit] %s %d targets (%d batches left)"
                      % (job["id"][:8], len(chunk), len(todo)), flush=True)
        if not pending:
            break
        jid, chunk = pending.pop(0)
        job = jf.job(jid, wait=90)
        while job.get("status") in ("queued", "running"):
            job = jf.job(jid, wait=90)
        by_id = {d["design_id"]: d for d in chunk}
        rows = []
        if job.get("status") == "succeeded":
            res = jf.results(jid)
            for row in res.get("rows", []):
                d = by_id.get(row["id"], {})
                cp = row.get("chains_ptm") or {}
                rows.append(dict(d, iptm=row.get("iptm"), interface_pae=row.get("interface_pae"),
                                 ipae=row.get("ipae"), plddt=row.get("plddt"), ptm=row.get("ptm"),
                                 chain_ptm_target=cp.get("0"), chain_ptm_binder=cp.get("1"),
                                 n_residues=row.get("n_residues"),
                                 binder_residues=row.get("binder_residues"),
                                 status=row.get("status", "ok"), job_id=jid))
            if write_pae:
                os.makedirs(PAE_DIR, exist_ok=True)
                for a in res.get("artifacts", []):
                    if a.get("type") == "pae":
                        jf.artifact(jid, a["path"], os.path.join(PAE_DIR, a["path"]))
            if job.get("partial"):
                print("[PARTIAL] %s delivered %s of %s"
                      % (jid[:8], job["partial"].get("delivered"),
                         job["partial"].get("requested")), flush=True)
        else:
            print("[FAILED] %s %s: %s" % (jid[:8], job.get("status"),
                                          str(job.get("error"))[:300]), flush=True)
            for d in chunk:
                rows.append(dict(d, status="job_" + str(job.get("status")), job_id=jid))
        try:
            spent += float(jf.usage(jid).get("amount") or 0)
        except Exception:
            pass
        append_rows(rows)
        done_rows.extend(rows)
        ok = sum(1 for r in rows if r.get("status") == "ok")
        if not quiet:
            print("[done]   %s %d/%d ok   cumulative $%.4f   total scored %d"
                  % (jid[:8], ok, len(chunk), spent, len(done_rows)), flush=True)
    return done_rows, spent


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["panel", "spr", "rank", "production", "pae", "credits"])
    ap.add_argument("--limit", type=int, default=None, help="candidates to score")
    ap.add_argument("--offset", type=int, default=0, help="skip this many ranked rows")
    ap.add_argument("--spread", type=int, default=0,
                    help="rank mode: take this many of --limit evenly across the whole file")
    ap.add_argument("--per-job", type=int, default=50, help="targets per job (<=8 MB body)")
    ap.add_argument("--in-flight", type=int, default=3, help="jobs submitted before waiting")
    ap.add_argument("--top", type=int, default=200, help="pae mode: how many top designs")
    ap.add_argument("--write-pae", action="store_true")
    ap.add_argument("--candidates", default=None,
                    help="rank mode: read candidates from this CSV instead of "
                         "data/legal_3mutants_ranked.csv (needs a `cdrh3` column)")
    ap.add_argument("--force", action="store_true", help="re-score designs already in the CSV")
    a = ap.parse_args()

    jf = JapanFold()
    if a.mode == "credits":
        c = jf.credits()
        print("balance $%.2f  rate %s" % (c["balance"], c["rate"]["description"]))
        return

    chains = parse_cif()
    wt_heavy = chain_seq(chains, "B")
    assert wt_heavy[CDRH3_START:CDRH3_END] == WT_CDRH3, "CDRH3 is not at 99-108"
    pdb = build_pdb(chains)

    if a.mode == "panel":
        items = load_panel()
    elif a.mode == "spr":
        items = load_spr(a.limit or 100)
    elif a.mode == "production":
        items = load_production() + load_panel()[:1]
    elif a.mode == "rank":
        # wild type first (the baseline every design is compared against), then the
        # candidates; the production batch has its own mode so parallel rank workers
        # do not all re-submit it.
        items = load_panel()[:1] + load_ranked(a.limit, a.offset, a.spread, a.candidates)
    elif a.mode == "pae":
        done = already_scored()
        scored = [r for r in done.values() if r.get("interface_pae")]
        scored.sort(key=lambda r: float(r["interface_pae"]))
        items = [dict(design_id=r["design_id"] + "_pae", cdrh3=r["cdrh3"], source="pae",
                      label=r.get("label", ""), kd_nM="", prior_score=r.get("prior_score", ""))
                 for r in scored[:a.top]]
        a.write_pae = True

    if not a.force:
        done = already_scored()
        seen, items2 = set(), []
        for d in items:
            if d["design_id"] in done or d["design_id"] in seen:
                continue
            seen.add(d["design_id"])
            items2.append(d)
        skipped = len(items) - len(items2)
        items = items2
        if skipped:
            print("skipping %d already-scored designs" % skipped)

    if not items:
        print("nothing to do")
        return
    before = jf.credits()["balance"]
    print("%d designs to score, %d per job, balance $%.2f" % (len(items), a.per_job, before))
    t0 = time.time()
    rows, spent = score(items, jf, pdb, wt_heavy, per_job=a.per_job,
                        in_flight=a.in_flight, write_pae=a.write_pae,
                        name_prefix="af2ig_" + a.mode)
    after = jf.credits()["balance"]
    ok = sum(1 for r in rows if r.get("status") == "ok")
    print("\n%d/%d scored in %.1f min, billed $%.4f (balance $%.2f -> $%.2f)"
          % (ok, len(items), (time.time() - t0) / 60.0, spent, before, after))
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
