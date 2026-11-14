#!/usr/bin/env python3
"""Dependency-free reconnaissance of the Zenodo 10831512 trastuzumab-NGS pretrained models.

Designed to run with `python3 -I` on a bare interpreter (no h5py / numpy / tensorflow).

- HDF5: walks the superblock + object headers just enough to list root-level
  attributes and group names, and scrapes the Keras `model_config` / `keras_version`
  JSON that Keras 2.x writes as a root attribute.
- Pickle: uses stdlib `pickletools` to disassemble without executing anything.

Usage:
    python3 -I inspect_pretrained.py h5   <file.h5>
    python3 -I inspect_pretrained.py pkl  <file.pkl>
    python3 -I inspect_pretrained.py csv  <file.csv>
"""
import json
import re
import sys


# ----------------------------------------------------------------- HDF5 ------
def h5_report(path):
    with open(path, "rb") as fh:
        blob = fh.read()

    print(f"== {path}  ({len(blob):,} bytes)")
    print(f"   magic: {blob[:8]!r}  -> {'HDF5' if blob[:8] == b'\\x89HDF\\r\\n\\x1a\\n' else 'NOT HDF5'}")

    # Keras 2.x writes these as root-group attributes, stored as plain UTF-8.
    for key in ("keras_version", "backend"):
        m = re.search(key.encode() + rb"[\x00-\x20]{0,32}?([\x20-\x7e]{1,24})", blob)
        if m:
            print(f"   {key}: {m.group(1).decode('ascii', 'replace')}")

    # model_config: a JSON document. Find the first balanced {...} after the key.
    idx = blob.find(b"model_config")
    if idx == -1:
        print("   model_config: ABSENT -> weights-only file (save_weights, not save_model)")
    else:
        start = blob.find(b"{", idx)
        depth, end = 0, None
        for i in range(start, min(start + 4_000_000, len(blob))):
            c = blob[i]
            if c == 0x7B:
                depth += 1
            elif c == 0x7D:
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break
        if end:
            cfg = json.loads(blob[start:end].decode("utf-8", "replace"))
            print("   model_config: present")
            summarize_keras_config(cfg)

    # HDF5 stores link/group names as readable strings in object headers.
    names = sorted(set(re.findall(rb"[A-Za-z_][A-Za-z0-9_/:.]{2,48}", blob[:200_000])))
    interesting = [
        n.decode() for n in names
        if any(t in n.lower() for t in (b"dense", b"conv", b"lstm", b"gru", b"embed", b"pool",
                                        b"drop", b"flatten", b"batch", b"norm", b"kernel",
                                        b"bias", b"layer", b"attention", b"input"))
    ]
    if interesting:
        print(f"   layer/weight names found in header: {interesting}")

    # layer_names attribute (weights-only files still carry this)
    m = re.search(rb"layer_names(.{0,4000}?)backend", blob, re.S)
    if m:
        toks = re.findall(rb"[a-z_]+_?\d*", m.group(1))
        print(f"   layer_names blob tokens: {[t.decode() for t in toks][:40]}")


def summarize_keras_config(cfg):
    cls = cfg.get("class_name")
    print(f"     model class: {cls}")
    layers = cfg.get("config", {}).get("layers", [])
    print(f"     {len(layers)} layers:")
    for L in layers:
        c = L.get("config", {})
        bits = []
        for k in ("batch_input_shape", "units", "filters", "kernel_size", "activation",
                  "rate", "pool_size", "input_dim", "output_dim"):
            if k in c and c[k] is not None:
                bits.append(f"{k}={c[k]}")
        print(f"       - {L.get('class_name'):<16} {c.get('name',''):<20} {', '.join(bits)}")


# --------------------------------------------------------------- pickle ------
def pkl_report(path, limit=400_000):
    import io
    import pickletools

    with open(path, "rb") as fh:
        head = fh.read(limit)
    print(f"== {path}")
    sink = io.StringIO()
    try:
        pickletools.dis(io.BytesIO(head), out=sink)
    except Exception as exc:  # truncated stream is expected
        sink.write(f"\n[truncated: {type(exc).__name__}: {exc}]\n")
    text = sink.getvalue()

    globals_ = sorted(set(re.findall(r"GLOBAL\s+'(\S+ \S+)'", text)))
    stacks = sorted(set(re.findall(r"SHORT_BINUNICODE '([a-zA-Z_][\w.]*)'[\s\S]{0,80}?STACK_GLOBAL", text)))
    print("   GLOBAL / STACK_GLOBAL symbols (classes the pickle will instantiate):")
    for g in globals_:
        print(f"     {g}")

    # Module paths referenced anywhere in the raw bytes - cheapest estimator fingerprint.
    mods = sorted(set(re.findall(rb"(?:sklearn|flaml|lightgbm|xgboost|catboost|numpy|scipy|pandas)"
                                 rb"[\w.]{0,60}", head)))
    print("   module paths in first chunk:")
    for m in mods:
        print(f"     {m.decode()}")

    for key in ("best_estimator", "best_config", "best_loss", "_estimator_type",
                "classes_", "feature_name", "n_features"):
        if key.encode() in head:
            print(f"   contains key: {key}")


# ------------------------------------------------------------------ csv ------
def csv_report(path, nrows=3):
    with open(path, "r", errors="replace") as fh:
        header = fh.readline().rstrip("\n")
        cols = header.split(",")
        print(f"== {path}")
        print(f"   {len(cols)} columns: {cols[:40]}{' ...' if len(cols) > 40 else ''}")
        for i in range(nrows):
            line = fh.readline()
            if not line:
                break
            vals = line.rstrip("\n").split(",")
            print(f"   row{i}: " + ", ".join(
                f"{c}={v}" for c, v in list(zip(cols, vals))[:12]))
        n = 2 + sum(1 for _ in fh)
    print(f"   total lines (incl header): {n:,}")


if __name__ == "__main__":
    kind, target = sys.argv[1], sys.argv[2]
    {"h5": h5_report, "pkl": pkl_report, "csv": csv_report}[kind](target)
