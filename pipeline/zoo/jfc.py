"""Tiny JapanFold REST client (stdlib only, key from env, never printed)."""
import json, os, time, urllib.request, urllib.error

BASE = os.environ.get("JAPANFOLD_BASE_URL", "https://api.japanfold.aiand.com").rstrip("/")
KEY = os.environ.get("JAPANFOLD_API_KEY", "")

DOM4 = ("CHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQNGSVTCFGPEADQ"
        "CVACAHYKDPPFCVARCPSGVKPDLSYMPIWKFPDEEGACQPCPIN")   # HER2 489-607, 119 aa
EPI_FULL = [557, 558, 560, 561, 569, 570, 571, 572, 573, 579,
            591, 593, 596, 598, 600, 602, 603, 604, 605]      # <5 A to Fab in 1N8Z
EPI_IDX = [r - 489 for r in EPI_FULL]                          # 0-based index into DOM4


def req(method, path, body=None, timeout=300, tries=4):
    data = json.dumps(body).encode() if body is not None else None
    last = None
    for k in range(tries):
        r = urllib.request.Request(BASE + path, data=data, method=method, headers={
            "Authorization": "Bearer " + KEY, "Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(r, timeout=timeout) as f:
                raw = f.read()
                return json.loads(raw) if raw else {}
        except urllib.error.HTTPError as e:
            txt = e.read().decode()[:1500]
            if e.code in (429, 500, 502, 503, 504) and k < tries - 1:
                time.sleep(3 * (k + 1)); last = (e.code, txt); continue
            return {"_http_error": e.code, "body": txt}
        except Exception as e:                      # network blip
            last = str(e); time.sleep(3 * (k + 1))
    return {"_http_error": -1, "body": str(last)}


def fetch(path, timeout=300):
    r = urllib.request.Request(BASE + path, headers={"Authorization": "Bearer " + KEY})
    with urllib.request.urlopen(r, timeout=timeout) as f:
        return f.read()


def complex_fasta(vh, vl, her2=DOM4):
    return f">H\n{vh}\n>L\n{vl}\n>T\n{her2}\n"


def submit(model, targets, name, params=None):
    body = {"model": model, "name": name, "targets": targets}
    if params: body["params"] = params
    return req("POST", "/v1/predictions", body)


def job(jid, wait=0):
    return req("GET", f"/v1/jobs/{jid}" + (f"?wait={wait}" if wait else ""), timeout=wait + 60)


def results(jid):
    return req("GET", f"/v1/jobs/{jid}/results")
