#!/usr/bin/env python3
"""Build the labelled benchmark panel used to bake off every model.

Ground truth, strongest first:
  * FACS high/low labels (524k variants, Chinery et al.) -- sampled stratified
    by edit distance, because the submission regime is distance 3-4 and a model
    that only separates distance 2 from distance 9 is useless to us.
  * SPR KD for 182 non-WT variants (continuous label).
Writes runs/model-zoo/panel.csv.  Deterministic (seed 0).
"""
import os, sys, zipfile
import numpy as np, pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FACS = os.path.join(REPO, "data/her2-trastuzumab/tz_her2_facs_3classes_524346.csv.zip")
SPR = os.path.join(REPO, "data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv")
OUT = os.path.join(REPO, "runs/model-zoo/panel.csv")
WT = "WGGDGFYAMD"
N_PER = {3: 30, 4: 30, 5: 20, 6: 20}   # per class per distance


def main():
    rng = np.random.default_rng(0)
    z = zipfile.ZipFile(FACS)
    f = pd.read_csv(z.open("tz_her2_facs_3classes_524346.csv"))
    rows = []
    for d, n in N_PER.items():
        for cls in ("high", "low"):
            sub = f[(f.edit_distance == d) & (f.binding_class == cls)]
            k = min(n, len(sub))
            take = sub.iloc[rng.choice(len(sub), k, replace=False)]
            for _, r in take.iterrows():
                rows.append(dict(id=f"facs_d{d}_{cls}_{r.cdrh3}", cdrh3=r.cdrh3, vh=r.vh_sequence,
                                 vl=r.vl_sequence, source="facs", facs=cls, dist=d, kd_nM=np.nan))
    s = pd.read_csv(SPR)
    q = s[s.binding_result == "quantified_binding"].drop_duplicates("cdrh3")
    for _, r in q.iterrows():
        d = sum(a != b for a, b in zip(r.cdrh3, WT)) + abs(len(r.cdrh3) - len(WT))
        rows.append(dict(id=f"spr_{r.cdrh3}", cdrh3=r.cdrh3, vh=r.vh_sequence, vl=r.vl_sequence,
                         source="spr", facs="", dist=d, kd_nM=r.kd_nM))
    p = pd.DataFrame(rows)
    # keep one row per distinct (vh): SPR rows duplicating a FACS row keep the SPR label
    p = p.sort_values("source", ascending=False).drop_duplicates("vh").sort_values(["source", "dist", "id"])
    p["y_binder"] = np.where(p.facs == "high", 1, np.where(p.facs == "low", 0, np.nan))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    p.to_csv(OUT, index=False)
    print(p.groupby(["source", "dist", "facs"]).size().to_string())
    print("total", len(p), "->", OUT)


if __name__ == "__main__":
    main()
