#!/usr/bin/env python3
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jf import post

SP = os.path.dirname(os.path.abspath(__file__))

# HER2 mature numbering. Continuous domain IV construct 489-607 (119 aa).
DOM4 = ("CHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTCFGPEADQ"
        "CVACAHYKDPPFCVARCPSGVKPDLSYMPIWKFPDEEGACQPCPIN")
assert len(DOM4) == 119
C_SEQ = DOM4[:92]      # 489-580 (chain C in reference cif)
D_SEQ = DOM4[102:]     # 591-607 (chain D)
assert len(C_SEQ) == 92 and len(D_SEQ) == 17

# Measured trastuzumab epitope (<5 A to Fab in 1N8Z), as 1-based positions of DOM4
EPI_FULL = [557, 558, 560, 561, 569, 570, 571, 572, 573, 579,
            591, 593, 596, 598, 600, 602, 603, 604, 605]
EPI_CORE = [571, 572, 573, 591, 593, 596]
to_dom4 = lambda rs: ",".join(str(r - 488) for r in rs)
to_c = lambda rs: ",".join(str(r - 488) for r in rs if r <= 580)
to_d = lambda rs: ",".join(str(r - 590) for r in rs if r >= 591)

JOBS = []


def spec_cont(epi):
    return (f"entities:\n"
            f"  - protein:\n      id: A\n      sequence: {DOM4}\n"
            f"      binding_types:\n        binding: {to_dom4(epi)}\n"
            f"  - protein:\n      id: H\n      sequence: 115..130\n")


def spec_split(epi):
    return (f"entities:\n"
            f"  - protein:\n      id: C\n      sequence: {C_SEQ}\n"
            f"      binding_types:\n        binding: {to_c(epi)}\n"
            f"  - protein:\n      id: D\n      sequence: {D_SEQ}\n"
            f"      binding_types:\n        binding: {to_d(epi)}\n"
            f"  - protein:\n      id: H\n      sequence: 115..130\n")


PLAN = [
    ("cont-full-s1", spec_cont(EPI_FULL), 1),
    ("cont-full-s2", spec_cont(EPI_FULL), 2),
    ("cont-full-s3", spec_cont(EPI_FULL), 3),
    ("cont-core-s4", spec_cont(EPI_CORE), 4),
    ("cont-core-s5", spec_cont(EPI_CORE), 5),
    ("split-full-s6", spec_split(EPI_FULL), 6),
    ("split-full-s7", spec_split(EPI_FULL), 7),
    ("split-core-s8", spec_split(EPI_CORE), 8),
]

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
B = int(sys.argv[2]) if len(sys.argv) > 2 else 20

out = []
for name, spec, seed in PLAN:
    r = post("/v1/designs", {"model": "boltzgen", "protocol": "nanobody-anything",
                             "name": "dn-" + name, "spec": spec,
                             "params": {"num_designs": N, "budget": B, "seed": seed}})
    if "id" in r:
        print("OK  ", name, r["id"], r["status"], "residues", r.get("residues"))
        out.append({"name": name, "id": r["id"], "seed": seed, "spec": spec})
    else:
        print("FAIL", name, json.dumps(r)[:600])

with open(os.path.join(SP, "bg_jobs.json"), "w") as f:
    json.dump(out, f, indent=1)
print(len(out), "jobs submitted")
