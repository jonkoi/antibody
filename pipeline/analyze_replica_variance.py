#!/usr/bin/env python3
"""Run-to-run (replica) variance analysis for structure-prediction scoring.

Context: our candidate ranking uses a SINGLE folding run per design. This script
quantifies how much repeat predictions of the *same* input disagree, using the
Smorodina et al. "Champloo" benchmark data (106 cognate / 11,342 non-cognate
nanobody-antigen pairs, each folded 50x by AlphaFold3, Boltz-2 and Chai-1).

Two independent read-outs:

  mode `rmsd`  - between-replica C-alpha RMSD (Angstrom), streamed on stdin from
                 rmsd_between_replicas.zip (one CSV per folding model).
                 Columns: structure_id_1,rmsd,structure_id_2,name,type,split_label
  mode `iptm`  - per-replica Chai-1 ipTM read straight out of
                 chai1_all_per_replica_confidence_metrics.zip (one .npz per
                 replica). Gives sd(ipTM) across the 50 replicas of one input,
                 i.e. the noise floor of the score we actually rank on.

No third-party deps: stdlib only (the analysis VM has no numpy/pandas), and the
.npy scalars inside each .npz are decoded by hand.
"""

import ast
import math
import struct
import sys
import zipfile
from array import array
from collections import defaultdict

# ---------------------------------------------------------------- statistics


def quantile(sorted_seq, q):
    """Linear-interpolated quantile of an already-sorted sequence."""
    if not sorted_seq:
        return float("nan")
    if len(sorted_seq) == 1:
        return float(sorted_seq[0])
    pos = q * (len(sorted_seq) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(sorted_seq) - 1)
    frac = pos - lo
    return float(sorted_seq[lo]) * (1 - frac) + float(sorted_seq[hi]) * frac


def describe(values, thresholds=()):
    """n/mean/sd/quartiles/max (+ tail fractions) for an iterable of floats."""
    v = sorted(values)
    n = len(v)
    if n == 0:
        return None
    mean = math.fsum(v) / n
    var = math.fsum((x - mean) ** 2 for x in v) / (n - 1) if n > 1 else 0.0
    out = {
        "n": n,
        "mean": mean,
        "sd": math.sqrt(var),
        "min": float(v[0]),
        "p25": quantile(v, 0.25),
        "median": quantile(v, 0.50),
        "p75": quantile(v, 0.75),
        "p90": quantile(v, 0.90),
        "p95": quantile(v, 0.95),
        "p99": quantile(v, 0.99),
        "max": float(v[-1]),
    }
    for t in thresholds:
        out["frac>%g" % t] = sum(1 for x in v if x > t) / n
    return out


def fmt(label, d, prec=3):
    if d is None:
        return "%-26s (no data)" % label
    keys = ["n", "mean", "sd", "min", "p25", "median", "p75", "p90", "p95", "p99", "max"]
    parts = []
    for k in keys:
        parts.append("%s=%s" % (k, d[k] if k == "n" else round(d[k], prec)))
    for k in sorted(d):
        if k.startswith("frac>"):
            parts.append("%s=%.4f" % (k, d[k]))
    return "%-26s %s" % (label, "  ".join(parts))


# ---------------------------------------------------------------- mode: rmsd


def system_of(name):
    """'system_100_17_8SOS_6HER_trunk_4' -> 'system_100_17_8SOS_6HER'."""
    i = name.find("_trunk_")
    return name[:i] if i >= 0 else name


def run_rmsd(model, dump=None):
    """Stream the pairwise-RMSD CSV on stdin; per-model + per-pair summaries."""
    all_vals = {"real": array("f"), "shuffled": array("f")}
    # per input system: n, sum, sumsq, max, #>2A, #>5A
    per_sys = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0, 0])
    sys_type = {}
    bad = 0

    inp = sys.stdin
    header = inp.readline().rstrip("\n").split(",")
    ix_rmsd = header.index("rmsd")
    ix_name = header.index("name")
    ix_type = header.index("type")

    for line in inp:
        f = line.rstrip("\n").split(",")
        try:
            r = float(f[ix_rmsd])
        except (ValueError, IndexError):
            bad += 1
            continue
        if r != r:  # NaN
            bad += 1
            continue
        t = f[ix_type]
        syst = system_of(f[ix_name])
        if t in all_vals:
            all_vals[t].append(r)
        else:
            all_vals.setdefault(t, array("f")).append(r)
        sys_type[syst] = t
        a = per_sys[syst]
        a[0] += 1
        a[1] += r
        a[2] += r * r
        if r > a[3]:
            a[3] = r
        if r > 2.0:
            a[4] += 1
        if r > 5.0:
            a[5] += 1

    print("=" * 100)
    print("MODEL: %s   (between-replica RMSD, Angstrom)   skipped/NaN rows: %d" % (model, bad))
    print("=" * 100)
    thr = (1.0, 2.0, 5.0, 10.0)
    for t in sorted(all_vals):
        print(fmt("pairwise RMSD [%s]" % t, describe(all_vals[t], thr)))

    # per-input-system view: is instability concentrated in a few hard inputs?
    for t in sorted(set(sys_type.values())):
        means, maxes, fr2 = [], [], []
        for s, a in per_sys.items():
            if sys_type[s] != t or a[0] == 0:
                continue
            means.append(a[1] / a[0])
            maxes.append(a[3])
            fr2.append(a[4] / a[0])
        print(fmt("per-system MEAN rmsd [%s]" % t, describe(means, (1.0, 2.0, 5.0))))
        print(fmt("per-system MAX  rmsd [%s]" % t, describe(maxes, (2.0, 5.0))))
        print(fmt("per-system frac>2A [%s]" % t, describe(fr2)))
        n_sys = len(means)
        if n_sys:
            tot_over2 = sum(per_sys[s][4] for s in per_sys if sys_type[s] == t)
            ms = sorted(means, reverse=True)
            k = max(1, n_sys // 10)
            print(
                "  [%s] systems=%d ; share of all >2A comparisons living in the "
                "worst-decile systems: see per-system mean spread (top-10%% mean "
                "rmsd >= %.2fA) ; total >2A comparisons=%d"
                % (t, n_sys, ms[k - 1], tot_over2)
            )

    # fraction of all unstable comparisons contributed by the worst systems
    for t in sorted(set(sys_type.values())):
        ss = sorted(((per_sys[s][4], s) for s in per_sys if sys_type[s] == t),
                    reverse=True)
        tot = sum(c for c, _ in ss)
        if not tot:
            continue
        for frac in (0.05, 0.10, 0.25, 0.50):
            k = max(1, int(len(ss) * frac))
            print("  [%s] worst %3d%% of systems hold %5.1f%% of all >2A comparisons"
                  % (t, frac * 100, 100.0 * sum(c for c, _ in ss[:k]) / tot))

    if dump:
        with open(dump, "w") as fh:
            fh.write("system,type,n,mean_rmsd,max_rmsd,frac_gt2\n")
            for s, a in per_sys.items():
                if a[0]:
                    fh.write("%s,%s,%d,%.4f,%.4f,%.4f\n"
                             % (s, sys_type[s], a[0], a[1] / a[0], a[3], a[4] / a[0]))
        print("  wrote per-system table -> %s" % dump)
    sys.stdout.flush()


# ---------------------------------------------------------------- mode: iptm


def read_npy_scalar(blob):
    """Decode a tiny .npy holding one float (no numpy available)."""
    if blob[:6] != b"\x93NUMPY":
        raise ValueError("not a .npy")
    major = blob[6]
    if major == 1:
        hlen = struct.unpack("<H", blob[8:10])[0]
        start = 10
    else:
        hlen = struct.unpack("<I", blob[8:12])[0]
        start = 12
    head = ast.literal_eval(blob[start:start + hlen].decode("latin1").strip())
    dt = head["descr"]
    data = blob[start + hlen:]
    code = {"<f4": "<f", "<f8": "<d", "|b1": "<?", "<i8": "<q", "<i4": "<i"}[dt]
    return struct.unpack(code, data[:struct.calcsize(code)])[0]


def parse_chai_name(member):
    """system_100_17_8SOS_6HER_trunk_4_scores.model_idx_9.npz -> (system, kind)."""
    base = member.split("/")[-1]
    stem = base.split("_trunk_")[0]
    parts = stem.split("_")  # system, idxA, idxB, pdbA, pdbB
    kind = "real" if len(parts) >= 5 and parts[1] == parts[2] else "shuffled"
    return stem, kind


def run_iptm(zip_path, limit_systems=None, dump=None):
    """sd(ipTM) across the replicas of each input, from the Chai-1 npz bundle."""
    zf = zipfile.ZipFile(zip_path)
    names = [n for n in zf.namelist() if n.endswith(".npz")]
    per_sys = defaultdict(list)
    kinds = {}
    for n in names:
        s, k = parse_chai_name(n)
        kinds[s] = k
        per_sys[s].append(n)

    systems = sorted(per_sys)
    if limit_systems:
        real = [s for s in systems if kinds[s] == "real"][:limit_systems]
        shuf = [s for s in systems if kinds[s] == "shuffled"][:limit_systems]
        systems = real + shuf

    print("=" * 100)
    print("CHAI-1 per-replica ipTM  (systems=%d of %d, npz files=%d)"
          % (len(systems), len(per_sys), len(names)))
    print("=" * 100)

    rows = []  # (system, kind, n, mean, sd, min, max, range)
    for s in systems:
        vals = []
        for member in per_sys[s]:
            try:
                with zipfile.ZipFile(zf.open(member)) as inner:
                    vals.append(float(read_npy_scalar(inner.read("iptm.npy"))))
            except Exception:
                continue
        if len(vals) < 2:
            continue
        d = describe(vals)
        rows.append((s, kinds[s], d["n"], d["mean"], d["sd"], d["min"], d["max"],
                     d["max"] - d["min"]))

    for kind in ("real", "shuffled"):
        sub = [r for r in rows if r[1] == kind]
        if not sub:
            continue
        print("-- %s pairs: %d systems, replicas/system median=%g" %
              (kind, len(sub), quantile(sorted(r[2] for r in sub), 0.5)))
        print(fmt("  sd(ipTM) within input", describe([r[4] for r in sub],
                                                     (0.01, 0.019, 0.02, 0.05)), 4))
        print(fmt("  range(ipTM) within in.", describe([r[7] for r in sub],
                                                       (0.019, 0.05, 0.1)), 4))
        print(fmt("  mean ipTM per input", describe([r[3] for r in sub]), 4))
        # noise floor conditional on a high mean ipTM (our designs score ~0.95)
        hi = [r for r in sub if r[3] >= 0.85]
        if hi:
            print(fmt("  sd(ipTM) | mean>=0.85", describe([r[4] for r in hi],
                                                          (0.01, 0.019, 0.02)), 4))
        hi9 = [r for r in sub if r[3] >= 0.90]
        if hi9:
            print(fmt("  sd(ipTM) | mean>=0.90", describe([r[4] for r in hi9],
                                                          (0.01, 0.019, 0.02)), 4))
        # how often does one replica mis-order two inputs whose true means differ?
        print(resolve_report([(r[3], r[4], r[2]) for r in sub], kind))
        # is the noise level the same at every score level? (heteroscedasticity)
        print("  sd(ipTM) as a function of the input's mean ipTM [%s]:" % kind)
        edges = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.01]
        for lo, hi_e in zip(edges, edges[1:]):
            b = [r[4] for r in sub if lo <= r[3] < hi_e]
            if b:
                bs = sorted(b)
                print("    mean ipTM %.2f-%.2f : n=%4d  mean sd=%.4f  median sd=%.4f"
                      "  p90 sd=%.4f  max sd=%.4f"
                      % (lo, hi_e, len(b), math.fsum(b) / len(b), quantile(bs, 0.5),
                         quantile(bs, 0.9), bs[-1]))

    if dump:
        with open(dump, "w") as fh:
            fh.write("system,type,n_replicas,mean_iptm,sd_iptm,min_iptm,max_iptm,range_iptm\n")
            for r in rows:
                fh.write("%s,%s,%d,%.5f,%.5f,%.5f,%.5f,%.5f\n" % r)
        print("wrote per-system ipTM table -> %s" % dump)
    sys.stdout.flush()


def tm_d0(length):
    if length <= 15:
        return 1.0
    return max(1.0, 1.24 * (length - 15) ** (1.0 / 3.0) - 1.8)


def af3_replica_metrics(doc):
    """Per-replica interface confidence from one AF3 confidences.json.

    The released AF3 files carry pae / plddt / contact_probs but NOT iptm, so we
    rebuild an ipTM-shaped interface score from the PAE matrix (the ipSAE-style
    approximation: mean of the TM kernel 1/(1+(e/d0)^2) over inter-chain token
    pairs, taken asymmetrically per chain and maximised, as ipTM does). Absolute
    values are not AF3's ipTM; the run-to-run SPREAD is what we use.
    """
    chains = doc["token_chain_ids"]
    pae = doc["pae"]
    uniq = sorted(set(chains))
    if len(uniq) != 2:
        return None
    a = [i for i, c in enumerate(chains) if c == uniq[0]]
    b = [i for i, c in enumerate(chains) if c == uniq[1]]
    d0 = tm_d0(min(len(a), len(b)))
    out = {}
    best = -1.0
    inter = []
    for rows, cols in ((a, b), (b, a)):
        tot = 0.0
        cnt = 0
        for i in rows:
            ri = pae[i]
            for j in cols:
                e = ri[j]
                tot += 1.0 / (1.0 + (e / d0) ** 2)
                inter.append(e)
                cnt += 1
        if cnt:
            best = max(best, tot / cnt)
    out["iface_tm"] = best
    out["iface_pae"] = math.fsum(inter) / len(inter) if inter else float("nan")
    inter.sort()
    out["iface_pae_p05"] = quantile(inter, 0.05) if inter else float("nan")
    plddt = doc.get("atom_plddts") or []
    out["plddt"] = math.fsum(plddt) / len(plddt) if plddt else float("nan")
    return out


def run_af3pae(tar_path, dump=None, max_systems=None):
    """Stream the AF3 per-replica tar.zst; replica spread of interface metrics."""
    import json
    import subprocess
    import tarfile

    proc = subprocess.Popen(["zstd", "-dc", tar_path], stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL)
    tf = tarfile.open(fileobj=proc.stdout, mode="r|")
    per_sys = defaultdict(lambda: defaultdict(list))
    seen = 0
    for member in tf:
        if not member.isfile() or not member.name.endswith(".json"):
            continue
        parts = member.name.split("/")
        syst = parts[-2]
        if max_systems and syst not in per_sys and len(per_sys) >= max_systems:
            continue
        fh = tf.extractfile(member)
        if fh is None:
            continue
        try:
            doc = json.loads(fh.read())
        except Exception:
            continue
        m = af3_replica_metrics(doc)
        if not m:
            continue
        for k, v in m.items():
            per_sys[syst][k].append(v)
        seen += 1
        if seen % 500 == 0:
            print("  ... %d replicas, %d systems" % (seen, len(per_sys)), flush=True)
    tf.close()
    proc.stdout.close()
    proc.wait()

    print("=" * 100)
    print("AF3 per-replica interface confidence (PAE-derived); %d replicas, %d systems"
          % (seen, len(per_sys)))
    print("=" * 100)
    rows = []
    for s, d in per_sys.items():
        if len(d.get("iface_tm", [])) < 2:
            continue
        t = describe(d["iface_tm"])
        p = describe(d["iface_pae"])
        l = describe(d["plddt"])
        rows.append((s, t["n"], t["mean"], t["sd"], t["max"] - t["min"],
                     p["mean"], p["sd"], l["mean"], l["sd"]))
    for key, idx, lo_hi in (("iface_tm sd", 3, None), ("iface_tm range", 4, None),
                            ("iface_pae mean", 5, None), ("iface_pae sd", 6, None),
                            ("mean pLDDT", 7, None), ("pLDDT sd", 8, None)):
        print(fmt(key, describe([r[idx] for r in rows]), 4))
    print("\n  replica spread of the interface score, split by the system's own level:")
    edges = [0.0, 0.3, 0.5, 0.7, 0.8, 0.9, 1.01]
    for lo, hi in zip(edges, edges[1:]):
        b = [r for r in rows if lo <= r[2] < hi]
        if not b:
            continue
        sds = sorted(r[3] for r in b)
        cv = [r[3] / r[2] for r in b if r[2] > 0]
        print("    iface_tm %.1f-%.1f : n=%3d  mean sd=%.4f  median sd=%.4f  "
              "p90 sd=%.4f  median CV=%.3f"
              % (lo, hi, len(b), math.fsum(sds) / len(sds), quantile(sds, 0.5),
                 quantile(sds, 0.9), quantile(sorted(cv), 0.5) if cv else float("nan")))
    if dump:
        with open(dump, "w") as fh:
            fh.write("system,n,iface_tm_mean,iface_tm_sd,iface_tm_range,"
                     "iface_pae_mean,iface_pae_sd,plddt_mean,plddt_sd\n")
            for r in rows:
                fh.write("%s,%d,%.5f,%.5f,%.5f,%.4f,%.4f,%.3f,%.3f\n" % r)
        print("wrote -> %s" % dump)
    sys.stdout.flush()


def load_iptm_matrix(zip_path, member):
    """{(pdb1,pdb2): iptm} from iptm_confidence_scores/*_matrix_clean.csv."""
    z = zipfile.ZipFile(zip_path)
    rows = z.read(member).decode().splitlines()
    hdr = rows[0].split(",")[1:]
    out = {}
    for row in rows[1:]:
        f = row.split(",")
        for h, v in zip(hdr, f[1:]):
            try:
                out[(f[0].upper(), h.upper())] = float(v)
            except ValueError:
                pass
    return out


def run_join(model, rmsd_csv, iptm_zip):
    """Does a confident single-run score imply a reproducible structure?"""
    mat = load_iptm_matrix(iptm_zip, "iptm_confidence_scores/%s_matrix_clean.csv" % model)
    recs = []
    with open(rmsd_csv) as fh:
        fh.readline()
        for line in fh:
            s, t, n, mean_r, max_r, f2 = line.rstrip("\n").split(",")
            p = s.split("_")
            if len(p) < 5:
                continue
            key = (p[3].upper(), p[4].upper())
            if key in mat:
                recs.append((mat[key], float(mean_r), float(max_r), float(f2), t))
    print("=" * 100)
    print("%s: aggregated ipTM vs between-replica structural spread (%d systems matched)"
          % (model, len(recs)))
    print("=" * 100)
    edges = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.01]
    for t in ("real", "shuffled"):
        sub = [r for r in recs if r[4] == t]
        if not sub:
            continue
        print("-- %s pairs (n=%d)" % (t, len(sub)))
        for lo, hi in zip(edges, edges[1:]):
            b = [r for r in sub if lo <= r[0] < hi]
            if not b:
                continue
            mr = sorted(x[1] for x in b)
            print("   ipTM %.1f-%.1f : n=%5d  median mean-replica-RMSD=%6.2fA  "
                  "p90=%6.2fA  frac of systems w/ mean RMSD<2A=%.3f"
                  % (lo, hi, len(b), quantile(mr, 0.5), quantile(mr, 0.9),
                     sum(1 for x in mr if x < 2.0) / len(mr)))
    sys.stdout.flush()


def norm_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def misrank(delta, sd, k=1):
    """P(a k-replica-mean comparison ranks two designs the wrong way round)."""
    se = sd * math.sqrt(2.0 / k)
    return norm_cdf(-delta / se)


def resolve_report(stats, kind):
    """Minimum resolvable difference and replica count, from the pooled noise sd."""
    sds = [sd for _, sd, _ in stats if sd == sd]
    if not sds:
        return ""
    pooled = math.sqrt(math.fsum(s * s for s in sds) / len(sds))
    med = quantile(sorted(sds), 0.5)
    hi = [sd for m, sd, _ in stats if m >= 0.90]
    pooled_hi = math.sqrt(math.fsum(s * s for s in hi) / len(hi)) if hi else float("nan")
    lines = ["  [%s] pooled within-input sd(ipTM)=%.4f (median %.4f, high-scorers %.4f)"
             % (kind, pooled, med, pooled_hi)]
    for use in (pooled, pooled_hi):
        if use != use:
            continue
        # two-sided 95% CI on a difference of two independent k-replica means
        lsd = 1.96 * math.sqrt(2.0) * use
        lines.append("    sd=%.4f -> 1 replica resolves only gaps > %.3f; "
                     "gap 0.019 needs k >= %d replicas per design"
                     % (use, lsd, max(1, math.ceil(2 * (1.96 * use / 0.019) ** 2))))
    return "\n".join(lines)


# ---------------------------------------------------------------- entry point

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "rmsd":
        run_rmsd(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
    elif mode == "join":
        run_join(sys.argv[2], sys.argv[3], sys.argv[4])
    elif mode == "af3pae":
        run_af3pae(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None,
                   int(sys.argv[4]) if len(sys.argv) > 4 else None)
    elif mode == "iptm":
        lim = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != "-" else None
        run_iptm(sys.argv[2], lim, sys.argv[4] if len(sys.argv) > 4 else None)
    else:
        raise SystemExit("usage: analyze_replica_variance.py "
                         "{rmsd MODEL [dump.csv] | iptm ZIP [N|-] [dump.csv]}")
