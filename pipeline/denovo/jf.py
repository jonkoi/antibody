#!/usr/bin/env python3
"""Minimal JapanFold client (stdlib only)."""
import json, os, sys, urllib.request, urllib.error

BASE = os.environ.get("JAPANFOLD_BASE_URL", "https://api.japanfold.aiand.com")
KEY = os.environ["JAPANFOLD_API_KEY"]


def req(method, path, body=None, timeout=180):
    data = json.dumps(body).encode() if body is not None else None
    r = urllib.request.Request(BASE + path, data=data, method=method,
                               headers={"Authorization": "Bearer " + KEY,
                                        "Content-Type": "application/json",
                                        "X-JapanFold-Client": "skill"})
    try:
        with urllib.request.urlopen(r, timeout=timeout) as f:
            return json.loads(f.read())
    except urllib.error.HTTPError as e:
        return {"_http_error": e.code, "body": e.read().decode()[:3000]}


def post(path, body, timeout=180):
    return req("POST", path, body, timeout)


def get(path, timeout=180):
    return req("GET", path, None, timeout)


if __name__ == "__main__":
    print(json.dumps(get(sys.argv[1]), indent=1)[:6000])
