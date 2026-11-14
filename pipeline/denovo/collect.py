#!/usr/bin/env python3
"""Fetch results + structures for finished BoltzGen jobs."""
import json, os, sys, urllib.request, urllib.error
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jf import get, BASE, KEY

SP = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(SP, "designs")
os.makedirs(OUT, exist_ok=True)

jobs = json.load(open(os.path.join(SP, "bg_jobs.json")))
allrows = []
for j in jobs:
    st = get("/v1/jobs/" + j["id"])
    if st.get("status") != "succeeded":
        print("SKIP", j["name"], st.get("status"), (st.get("error") or "")[:200])
        continue
    res = get("/v1/jobs/%s/results" % j["id"])
    if allrows == [] :
        print("RESULT KEYS:", list(res.keys()))
        d0 = (res.get("designs") or [{}])[0]
        print("DESIGN ROW KEYS:", list(d0.keys()))
        print("SAMPLE:", json.dumps(d0)[:1500])
    ds = res.get("designs") or []
    arts = {a.get("name") or a.get("path") or "": a for a in (res.get("artifacts") or [])}
    print(j["name"], "designs:", len(ds), "artifacts:", len(arts))
    jd = os.path.join(OUT, j["name"])
    os.makedirs(jd, exist_ok=True)
    for a in (res.get("artifacts") or []):
        nm = a.get("name") or os.path.basename(a.get("url", ""))
        if not nm or not nm.endswith((".cif", ".pdb")):
            continue
        p = os.path.join(jd, nm)
        if os.path.exists(p):
            continue
        url = a["url"]
        if url.startswith("/"):
            url = BASE + url
        r = urllib.request.Request(url, headers={"Authorization": "Bearer " + KEY})
        try:
            with urllib.request.urlopen(r, timeout=120) as f, open(p, "wb") as g:
                g.write(f.read())
        except urllib.error.HTTPError as e:
            print("  dl fail", nm, e.code)
    for d in ds:
        d["_job"] = j["name"]
        d["_jobid"] = j["id"]
        d["_seed"] = j["seed"]
        allrows.append(d)

json.dump(allrows, open(os.path.join(SP, "designs_raw.json"), "w"), indent=1)
print("total design rows:", len(allrows))
