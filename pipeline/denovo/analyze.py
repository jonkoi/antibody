#!/usr/bin/env python3
"""Join BoltzGen metrics with structural target-fold fidelity and epitope contacts."""
import glob, json, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from epitope import read_atoms, chain_seq

SP = os.path.dirname(os.path.abspath(__file__))
REF = "/home/duy.khoi.tran@esmartsystems.com/antibody/inputs/12-her2/her2_fv_1n8z.cif"
EPI = [557, 558, 560, 561, 569, 570, 571, 572, 573, 579,
       591, 593, 596, 598, 600, 602, 603, 604, 605]
CORE4 = [572, 573, 591, 593]

ref = read_atoms(REF)
R = {a[1]: a[4:7] for a in ref if a[3] == "CA" and a[0] in ("C", "D")}


def kabsch(P, Q):
    P = np.array(P, float); Q = np.array(Q, float)
    P = P - P.mean(0); Q = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Rm = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt((((P @ Rm.T) - Q) ** 2).sum(1).mean()))


def analyse(path, arm):
    atoms = read_atoms(path)
    chains = {}
    for a in atoms:
        chains.setdefault(a[0], []).append(a)
    lens = {c: len(chain_seq(atoms, c)) for c in chains}
    if arm == "cont":
        tgt = {"A": 488}
        bnd = "B"
    else:
        tgt = {"A": 488, "B": 590}
        bnd = "C"
    if bnd not in chains:
        return None
    # target fold fidelity
    D = {}
    for c, off in tgt.items():
        for a in chains.get(c, []):
            if a[3] == "CA":
                D[a[1] + off] = a[4:7]
    com = sorted(set(R) & set(D))
    rmsd = kabsch([D[x] for x in com], [R[x] for x in com]) if len(com) > 20 else None
    ce = [x for x in com if x in EPI]
    erms = kabsch([D[x] for x in ce], [R[x] for x in ce]) if len(ce) > 4 else None
    # contacts
    bd = [a[4:7] for a in chains[bnd]]
    B = np.array(bd)
    hit = set()
    for c, off in tgt.items():
        for a in chains.get(c, []):
            d = np.sqrt(((B - np.array(a[4:7])) ** 2).sum(1)).min()
            if d < 5.0:
                hit.add(a[1] + off)
    ep = set(EPI)
    inter = hit & ep
    return {
        "seq": chain_seq(atoms, bnd),
        "tgt_rmsd": None if rmsd is None else round(rmsd, 2),
        "epi_rmsd": None if erms is None else round(erms, 2),
        "n_contact": len(hit),
        "n_epi_hit": len(inter),
        "epi_recall": round(len(inter) / len(ep), 3),
        "epi_precision": round(len(inter) / max(len(hit), 1), 3),
        "core4": len(hit & set(CORE4)),
        "epi_hits": sorted(inter),
    }


if __name__ == "__main__":
    rowsf = sys.argv[1] if len(sys.argv) > 1 else SP + "/rows_w1.json"
    rows = json.load(open(rowsf))
    byjob = {}
    for d in rows:
        byjob.setdefault(d["_job"], {})[os.path.basename(d["structure"])] = d
    out = []
    for job, m in byjob.items():
        arm = "split" if job.startswith("split") else "cont"
        for fn, d in m.items():
            p = os.path.join(SP, "designs", job, fn)
            if not os.path.exists(p):
                continue
            a = analyse(p, arm)
            if a is None:
                continue
            rec = {"id": "%s_r%s" % (job, d["final_rank"]), "job": job, "arm": arm,
                   "rank": int(d["final_rank"]),
                   "iptm": round(float(d["design_to_target_iptm"]), 3),
                   "minpae": round(float(d["min_design_to_target_pae"]), 2),
                   "dsasa": round(float(d["delta_sasa_refolded"])),
                   "hbonds": int(d["plip_hbonds_refolded"]),
                   "liab": int(d["liability_num_violations"]),
                   "ptm": round(float(d["design_ptm"]), 3),
                   "len": len(a["seq"])}
            rec.update(a)
            out.append(rec)
    json.dump(out, open(SP + "/analysis.json", "w"), indent=1)
    print("analysed", len(out))
    ok = [r for r in out if r["tgt_rmsd"] is not None and r["tgt_rmsd"] < 3.0]
    print("target folded correctly (<3A):", len(ok), "/", len(out))
    import statistics as st
    for label, sel in (("all", out), ("tgt_ok", ok)):
        if not sel: continue
        print(label, "n=%d" % len(sel),
              "iptm med %.3f max %.3f" % (st.median(r["iptm"] for r in sel), max(r["iptm"] for r in sel)),
              "| epi_recall med %.2f max %.2f" % (st.median(r["epi_recall"] for r in sel), max(r["epi_recall"] for r in sel)),
              "| core4>=3: %d" % sum(1 for r in sel if r["core4"] >= 3))
    print()
    print("TOP by (core4, epi_recall, iptm):")
    for r in sorted(out, key=lambda r: (-r["core4"], -r["epi_recall"], -r["iptm"]))[:15]:
        print("  %-22s arm=%-5s iptm=%.3f pae=%5.2f rmsd=%5s epiR=%.2f epiP=%.2f core4=%d dSASA=%4d hb=%d liab=%d len=%d"
              % (r["id"], r["arm"], r["iptm"], r["minpae"], r["tgt_rmsd"], r["epi_recall"],
                 r["epi_precision"], r["core4"], r["dsasa"], r["hbonds"], r["liab"], r["len"]))
