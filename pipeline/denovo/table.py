#!/usr/bin/env python3
"""Join everything into one ranked candidate table."""
import csv, json, os, sys, statistics as st

SP = os.path.dirname(os.path.abspath(__file__))
norm = lambda s: s.replace("-", "_").replace("_SS", "")


def load():
    an = {norm(r["id"]): r for r in json.load(open(SP + "/analysis.json"))}
    tm = {norm(r["ID"]): float(r["NanoMelt Tm (C)"])
          for r in csv.DictReader(open(SP + "/w1_tm.csv"))}
    tms = {norm(r["ID"]): float(r["NanoMelt Tm (C)"])
           for r in csv.DictReader(open(SP + "/w1_ss_tm.csv"))}
    nv = {norm(r["sequence_id"]): r for r in json.load(open(SP + "/w1_novelty.json"))}
    nvs = {norm(r["sequence_id"]): r for r in json.load(open(SP + "/w1_ss_novelty.json"))}
    ss = {norm(k): v for k, v in json.load(open(SP + "/ss_map.json")).items()}
    valid = {norm(r["sequence_id"]) for r in json.load(open(SP + "/w1_ss_valid.json"))}
    af = {}
    for f in ("af2ig_orig.json",):
        p = os.path.join(SP, f)
        if os.path.exists(p):
            af.update({norm(k): v for k, v in json.load(open(p)).items()})
    out = []
    for k, a in an.items():
        r = dict(a)
        r["key"] = k
        r["tm_orig"] = tm.get(k)
        r["tm_ss"] = tms.get(k)
        r["ident"] = (nv.get(k) or {}).get("identity_to_trastuzumab")
        r["ident_ss"] = (nvs.get(k) or {}).get("identity_to_trastuzumab")
        r["class_ss"] = (nvs.get(k) or {}).get("novelty_class")
        r["cdrh3"] = (nv.get(k) or {}).get("cdrh3")
        r["cdrh3_dist"] = (nv.get(k) or {}).get("cdrh3_min_distance")
        r["abguard_pass"] = k in valid
        r["seq_ss"] = (ss.get(k) or {}).get("ss")
        a2 = af.get(k, {})
        r["af_iptm"] = a2.get("iptm"); r["af_ipae"] = a2.get("interface_pae")
        out.append(r)
    return out


if __name__ == "__main__":
    rows = load()
    # de-dup identical sequences, keep best epitope recall
    byseq = {}
    for r in rows:
        byseq.setdefault(r["seq"], []).append(r)
    uniq = [max(v, key=lambda r: (r["epi_recall"], r["iptm"])) for v in byseq.values()]
    print("unique designs:", len(uniq), "of", len(rows), "rows")

    def ok(r):
        return (r["abguard_pass"] and r["seq_ss"] and r["tm_ss"]
                and r["ident_ss"] is not None and r["ident_ss"] < 0.70)
    g = [r for r in uniq if ok(r)]
    print("pass de-novo gate WITH restored disulfide (abguard ok + ident_ss<0.70 + Tm):", len(g))

    def sc(r):
        return (r["core4"], r["epi_recall"], -(r["tgt_rmsd"] or 99) < -3.0, r["tm_ss"])
    g.sort(key=sc, reverse=True)
    hdr = ("rank id arm bgIPTM afIPTM afIPAE tgtRMSD epiRMSD epiR epiP core4 dSASA hb "
           "len Tm_SS identSS cdrh3 d3")
    print()
    print("%-3s %-19s %-5s %6s %6s %6s %7s %7s %5s %5s %5s %5s %3s %4s %6s %7s %-22s %3s" % tuple(hdr.split()))
    for i, r in enumerate(g[:25], 1):
        print("%-3d %-19s %-5s %6.3f %6.3f %6.1f %7.2f %7.2f %5.2f %5.2f %5d %5d %3d %4d %6.1f %7.3f %-22s %3s"
              % (i, r["key"], r["arm"], r["iptm"], r["af_iptm"] or 0, r["af_ipae"] or 0,
                 r["tgt_rmsd"] or 0, r["epi_rmsd"] or 0, r["epi_recall"], r["epi_precision"],
                 r["core4"], r["dsasa"], r["hbonds"], r["len"], r["tm_ss"], r["ident_ss"],
                 r["cdrh3"] or "?", r["cdrh3_dist"]))
    json.dump(g, open(SP + "/candidates.json", "w"), indent=1)
    print()
    print("Tm_SS of gated set: median %.1f  max %.1f  n>67.2: %d  n>70: %d"
          % (st.median(r["tm_ss"] for r in g), max(r["tm_ss"] for r in g),
             sum(1 for r in g if r["tm_ss"] > 67.2), sum(1 for r in g if r["tm_ss"] > 70)))
    print("target refolded <3A in gated set: %d/%d"
          % (sum(1 for r in g if r["tgt_rmsd"] and r["tgt_rmsd"] < 3.0), len(g)))
