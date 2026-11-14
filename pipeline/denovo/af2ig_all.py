#!/usr/bin/env python3
"""Submit AF2-IG scoring for every BoltzGen design, in batches, as compact PDB."""
import glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jf import post
from epitope import read_atoms, chain_seq

SP = os.path.dirname(os.path.abspath(__file__))
ELEM = {"C": "C", "N": "N", "O": "O", "S": "S", "H": "H"}


def to_pdb(atoms, keep):
    out = []
    i = 0
    for a in atoms:
        if a[0] not in keep:
            continue
        i += 1
        nm = a[3]
        an = (" %-3s" % nm) if len(nm) < 4 else nm
        out.append("ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
                   % (i, an, a[2], a[0], a[1], a[4], a[5], a[6], ELEM.get(nm[0], "C")))
    out.append("END")
    return "\n".join(out)


def build(path, tgt, bnd, use_seq=None):
    atoms = read_atoms(path)
    pdb = to_pdb(atoms, set(tgt) | {bnd})
    seq = use_seq or chain_seq(atoms, bnd)
    return ("target:\n  structure: |\n" +
            "".join("    " + l + "\n" for l in pdb.splitlines()) +
            "  chain: %s\n" % tgt[0] +
            "binder:\n  sequence: %s\n  chain: %s\n" % (seq, bnd)), seq


if __name__ == "__main__":
    rows = json.load(open(SP + "/rows_w1.json"))
    ssmap = json.load(open(SP + "/ss_map.json"))
    variant = sys.argv[1] if len(sys.argv) > 1 else "orig"   # orig | ss
    targets = []
    for d in rows:
        job = d["_job"]
        arm = "split" if job.startswith("split") else "cont"
        tgt, bnd = (["A", "B"], "C") if arm == "split" else (["A"], "B")
        p = os.path.join(SP, "designs", job, os.path.basename(d["structure"]))
        if not os.path.exists(p):
            continue
        did = "%s_r%s" % (job, d["final_rank"])
        seq = None
        if variant == "ss":
            e = ssmap.get(did)
            if not e or e["ss"].count("C") != 2:
                continue
            seq = e["ss"]
        y, s = build(p, tgt, bnd, seq)
        targets.append({"name": did.replace("-", "_") + ("_SS" if variant == "ss" else ""),
                        "content": y})
    # de-duplicate by name
    seen, uniq = set(), []
    for t in targets:
        if t["name"] in seen:
            continue
        seen.add(t["name"]); uniq.append(t)
    print("targets", len(uniq), "variant", variant)
    jobs = []
    B = 16
    for k in range(0, len(uniq), B):
        batch = uniq[k:k + B]
        r = post("/v1/predictions", {"model": "af2ig", "name": "dn-af2ig-%s-%d" % (variant, k),
                                     "params": {"output_format": "pdb"}, "targets": batch},
                 timeout=300)
        if "id" in r:
            print("OK batch", k, r["id"], len(batch))
            jobs.append(r["id"])
        else:
            print("FAIL batch", k, json.dumps(r)[:400])
    json.dump(jobs, open(SP + "/af2ig_jobs_%s.json" % variant, "w"))
