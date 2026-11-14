# Execution — technical

14 Nov 2026. The commands, paths and gotchas. The argument and the decisions are
in **[STRATEGY.md](STRATEGY.md)**; the checklist is [TASKS.md](TASKS.md).

---

## The workflow

### What we built vs. what we should run

The pipeline folds first and discovers novelty failures afterwards. That wastes
the GPU on designs that cannot be submitted. **Move the novelty gate before the
fold.**

```
  CURRENT  (wasteful)
                                                      ┌──────────────┐
   RFD4 ──▶ refold-prepare ──▶ fold top 20 ──▶ score ─▶│ abguard      │──▶ 1 of 3
   design      (rank)          (GPU, slow)             │ team-validate│    survive
                                                       └──────────────┘
                                             GPU already spent ───────┘

  CORRECTED  (novelty first)
                 ┌──────────────┐
   RFD4 ──▶ ─────│ abguard      │──▶ survivors ──▶ refold-prepare ──▶ fold ──▶ score
   design   100  │ team-validate│     (~?)           (rank)          (GPU)
                 │ --reference-db│
                 └──────────────┘
                   CPU, seconds          GPU spent only on submittable designs
```

### The parallel plan

Only the three analysis lanes are genuinely independent. The GPU is a single
serialized resource and no orchestration changes that.

```
 TIER 0 ── humans, minutes ── unblocks everything
   ├─ confirm deadline ─────────────┐
   ├─ re-share corrupt deck ────────┤
   └─ CHASE NEBIUS CREDENTIALS ─────┼──▶ unblocks stab-ddg (thermostability)
                                    │
 ══════════════════════════════════ │ ═══════════════════════════════════════
                                    │
 GPU LANE (serial — one H100)       │   CPU LANES (parallel — independent data)
   ┌────────────────────────┐       │    ┌──────────────────────────────────┐
   │ batch n=100 RUNNING    │       │    │ A  feasible-novel region         │
   │   ↓                    │       │    │    FACS 524k + reference DB      │
   │ novelty-filter (CPU)   │       │    │    ★ pitch centrepiece           │
   │   ↓                    │       │    ├──────────────────────────────────┤
   │ fold survivors         │       │    │ B  calibrate cut-off             │
   │   ↓                    │       │    │    ~8,100 shuffled pairs         │
   │ rank on min_ipae       │       │    ├──────────────────────────────────┤
   └────────────┬───────────┘       │    │ C  re-derive findings 1-5        │
                │                   │    │    adversarial, raw data only    │
                │                   │    └──────────────┬───────────────────┘
                │                   │                   │
 PhD LANE (humans, parallel)        │                   │
   ┌────────────────────────┐       │                   │
   │ epitope hotspots  ◀────┼───────┘                   │
   │ loop drift?            │                           │
   │ double proline?        │                           │
   └────────────┬───────────┘                           │
                │                                       │
                └───────────────┬───────────────────────┘
                                ▼
                    ┌───────────────────────┐
                    │  VALIDATION GATE      │   re-run every script,
                    │  (me, not the agents) │   check numbers vs raw data
                    └───────────┬───────────┘
                                ▼
                    ┌───────────────────────┐
                    │  PITCH  (single author)│  ← written last, holds everything
                    │  + submission table    │
                    └───────────────────────┘
```

**Why not a dev/test/writer hierarchy.** Our deliverables are analyses, not
features: a script can run perfectly and still produce a meaningless number, so
functional tests verify the wrong thing. Lane C replaces the test role — an
independent re-derivation from raw data, *not shown our conclusions*, because an
agent handed a hypothesis tends to return it confirmed. We have already caught
two wrong claims in HANDOFF and a real bug in `abguard`, so this is not
theoretical.

**The validation gate.** Reports come back equally confident whether right or
wrong, so each lane must (1) leave a reproducible, dependency-free script in
`pipeline/`, matching `analyze_spr.py`'s style, (2) print the raw numbers it
claims with row counts, and (3) state what it could not verify. Then the numbers
get re-run and spot-checked against the data. That is the only validation that
means anything here.

---

## Access

```bash
ssh -i ~/.ssh/team-15 hackathon@89.169.96.144        # user: hackathon, team-15
```

One H100 80GB, `docker`, `aws` (credentials work on the VM only), and `abguard`
at `~/.local/bin/abguard`.

- **The home directory is shared across the team** — same user, same key, one GPU.
  Teammates' folders: `igor/`, `foundry-assistant/`. Keep our work in
  `~/team15-claude/` and `~/rfd4/runs/team15-claude-*`.
- **Check `nvidia-smi` before launching anything on the GPU.**
- `ssh host "cmd"` does not source `.bashrc`, so `abguard`, `uv` and `claude` look
  missing. Use `ssh host 'bash -lc "cmd"'`.
- Long jobs need `setsid nohup ... &` or they die when the SSH session closes.

---

## Data

Local at `hackathon-data/` (gitignored), and on the VM at `~/hackathon-data/`.

| Path | Contents |
| --- | --- |
| `her2-trastuzumab/tz_her2_facs_3classes_524346.csv.zip` | 524,346 variants; `high` 172,149 / `medium` 190,011 / `low` 162,186; `edit_distance` 1–10 |
| `her2-trastuzumab/trastuzumab_her2_spr_all_768.csv` | 768 rows, **190 with `kd_nM`**; wild-type appears 8× |
| `sigma2-ab-ag-champloo/` | per-sample ipTM/pLDDT/PAE + `spr-supplementary/` (**29 + 12 SPR rows, NOT 454, and not usable as negatives**) |
| `zenodo-18390239-ab-ag-complexes/` | scores for 1.8M predicted nanobody–antigen complexes |
| `zenodo-10831512-trastuzumab-ngs/` | ready-made train/val/test splits of the FACS data |

Gotchas:

- **`unzip` is not installed** on the dev box — use Python's `zipfile`.
- The SPR affinity column is **`kd_nM`**, not `KD`.
- **CDRH3s are not all 10 residues** in the SPR table (20 are 9, 8 are 16), so
  substitution counts need edit distance, not Hamming. The "always 10 residues"
  claim applies to the FACS table only.
- Wild-type **is** in the SPR table (8 replicates) though the README says
  trastuzumab is absent — that claim is about the FACS table.
- **Cite Chinery et al., doi:10.1101/2024.03.26.586756** for both tables.

---

## abguard

```bash
cd ~/team15-claude
abguard team-validate abguard_submission.csv --reference-db her2_novelty.duckdb
```

Submission table columns: `sequence_id, kind, heavy_sequence, light_sequence,
design_description`. `kind` is `vhh` / `scfv` / `paired`; ours is `paired`.

**The reference DB is already built** at `~/team15-claude/her2_novelty.duckdb` —
525,063 distinct CDRH3s, exactly the number the brief quotes, in 37s:

```bash
abguard reference-ingest her2_novelty.duckdb \
  ~/hackathon-data/her2-trastuzumab/tz_her2_facs_3classes_524346.csv.zip \
  ~/hackathon-data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv
abguard reference-build her2_novelty.duckdb
```

Gotchas — two of these contradict HANDOFF:

- **`--reference-db` is required for the novelty gate.** Without it, prior art is
  *not* checked: the tool warns `--no-novelty: prior art was not checked` and
  returns `is_novel_enough = null`. HANDOFF's claim that a 610,109-sequence
  reference is used by default does not hold in practice.
- **`abguard novelty-report` is buggy** — it ignores `--reference-db` and emits
  the `--no-novelty` warning even when the flag is absent, comparing against
  trastuzumab alone. **Use `team-validate` for any verdict.**
- **The IMGT CDRH3 is not our design window.** abguard scores `SRWGPPYFYAMDY`
  (13 residues) while the model designs `WGPPYFYAMD` (10, B99–108) — IMGT adds
  `SR` before and `Y` after. Those three are frozen at wild-type so counts are
  unaffected, but `identity_to_trastuzumab` is computed over the whole paired
  sequence (96.9–98.7%), not the window. Everything lands in the **optimisation
  track**, needing ≥3 mutations; nothing built this way can reach de novo.
- Use `--drop` to discard non-compliant rows instead of failing the whole table.
- `abguard dashboard` launches a marimo liability UI — needs an SSH port-forward.

### Current verdict

```
ERROR  1 of 3 row(s) pass the schema
  - [novelty] control_wildtype_trastuzumab: 0 mutation(s) from trastuzumab; needs >= 3
  - [novelty] her2_cdrh3_2: CDRH3 SRWGPPYFYAMDY is 2 mutation(s)
      from SRWGMPHFYAMDY in tz_her2_facs_3classes_524346; needs >= 3
```

`her2_cdrh3_3` passes. The wild-type control fails as designed, which is how we
know the gate works.

---

## The running batch

Launched ~11:55 CET: 100 designs, seed 1, weights `her2-novel-2026-11-09`, top 20
to Boltz-2.

```bash
cd ~/rfd4 && OUT=runs/team15-claude-n100-seed1 N=100 TOPN=20 SEED=1 \
  setsid nohup bash ~/team15-claude/run_design.sh all \
  > ~/team15-claude/batch-n100.log 2>&1 &
```

Progress: `tail -f ~/team15-claude/batch-n100.log`.
Output: `~/rfd4/runs/team15-claude-n100-seed1/`.

Then:

```bash
python3 pipeline/analyze_designs.py \
  --designs runs/<run>/designs \
  --reference inputs/12-her2/her2_fv_1n8z.cif \
  --scored runs/<run>/scored/refold_metrics.csv
```

Ranking notes: the `pass` column is **always empty** for the HER2 example — the
scorer cannot measure absolute movement. Rank on `min_ipae`, which showed real
spread (0.310 vs 0.331), **not** on ipTM or pLDDT, which are saturated.

---

## Batch results — done

```bash
# 1. all 100 designs -> submission table (+ wild-type control that must fail)
python3 pipeline/make_submission.py \
  --designs runs/team15-claude-n100-seed1/design/generation/inference_samples \
  --reference inputs/12-her2/her2_fv_1n8z.cif \
  --out runs/team15-claude-n100-seed1/submission_all.csv

# 2. novelty gate, on the VM
abguard team-validate submission_all.csv --reference-db her2_novelty.duckdb \
  --drop -o accepted.csv          # 71 of 100 accepted

# 3. refold the designs WE chose, not refold-prepare's ranking
IDS="her2_cdrh3_60 her2_cdrh3_94 ..." bash ~/team15-claude/refold_subset.sh
```

Artefacts: `submission_all.csv` (101 rows), `accepted.csv` (71),
`subset-scored/refold_metrics.csv` (15 refolded), **`shortlist.csv` (14 unique,
W99-retaining first — this is the submission pool).**

[pipeline/refold_subset.sh](pipeline/refold_subset.sh) exists because
`refold-prepare --top-n` ranks on interface H-bonds and hotspot obedience, which
picked almost only high-substitution designs and wasted 30% of the fold on
novelty-rejected ones. It runs prepare over all 100, filters `holo/queries.json`
to chosen ids, then folds and scores.

### Numbers

| subs | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| designs | 2 | 6 | 5 | 3 | 17 | 35 | 28 | 4 |
| accepted | 0 | 0 | 0 | 1 | 14 | 26 | 26 | 4 |
| rate | 0% | 0% | 0% | 33% | 82% | 74% | 93% | 100% |

- `her2_cdrh3_50` and `_69` were **0 mutations** from a FACS variant — the model
  re-invented already-measured sequences.
- W99 retention: 8% of accepted vs 79% of rejected; SPR median 18.8 nM with W99
  vs 147 nM without.
- `_85` and `_93` are byte-identical (`STAAVGYAMD`); deduplicated in the shortlist.
- Candidate `min_ipae` spans 0.332–0.395 and ipTM 0.82–0.90 — more spread than the
  n=4 validation run, so `min_ipae` does discriminate at this scale.

## The remaining analysis

### A — Feasible-novel region — **ANSWERED**, see Batch results above.
Empty below 5 substitutions; viable at 5–6. No further work needed.

### B — Calibrate the cut-off — **DONE**, see [findings/lane-b-confidence.md](findings/lane-b-confidence.md)

```bash
python3 pipeline/calibrate_confidence.py   # reproduced; numbers below verified
python3 pipeline/calibrate_replicas.py
```

Negatives are the **shuffled pairs** in
`zenodo-18390239-ab-ag-complexes/iptm_confidence_scores.zip` (~8,100 per model),
presumed non-binding by the study's design. ROC AUC: **AF3 0.754, Boltz-2 0.645,
Chai-1 0.568.** Boltz-2 — the model we fold with — at ipTM ≥ 0.70 gives **94.3%
false positives** (AF3: 3.0%); a 5% FPR needs ipTM ≥ 0.951 and keeps 18% of real
pairs. Our candidates at ipTM 0.82–0.90 sit inside the shuffled interquartile
range. Repeat folding does not help: all replica aggregations land at AUC
0.510–0.543.

**Do not transfer the threshold.** Our refold is template-guided with the
framework and CDRH1/H2 frozen, so our ipTM is not on this scale. What transfers
is the model ordering and the shape of the trade-off.

**The "454 SPR-measured non-binders, 21 at ipTM ≥ 0.7" claim was false** and is
corrected throughout. `spr-supplementary/` holds 29 and 12 rows; "454" appears
nowhere in the data; and its README states "N/A denotes no reportable affinity
fit, NOT confirmed non-binding", so those rows cannot serve as negatives. Only 3
rows carry an AF3 ipTM at all.

### C — Adversarial re-derivation

Independently re-derive the five findings in STRATEGY.md from the raw tables.
**Must not be shown our conclusions** — give it the data and the question only.
Report confirmed / refuted / could-not-verify.

---

## Reproducing the findings

```bash
python3 pipeline/analyze_spr.py          # noise floor, substitution cost, novelty tension
python3 pipeline/analyze_designs.py ...  # window diff + metric join (pre-filter only)
```

`analyze_designs.py`'s `vs_wt` column is a **pre-filter, not a verdict** — it
compares to wild-type alone, and it passed `her2_cdrh3_2`, which abguard rejected.

---

## Sources and secrets

Full detail in [CONNECT.md](CONNECT.md). AMASS is configured and its key verified;
it needs a restart to load:

```bash
cd ~/antibody && set -a && . ./.env && set +a && claude
```

**Nebius is live.** `credentials.private.json` was recovered from
`participant-12.zip` and installed at the repo root (mode 600, gitignored by
`*.private.json`). Tooling vendored to `vendor/nebius/` (`forge.py`, `SKILL.md`,
`RECIPES.md`).

```bash
python3 vendor/nebius/forge.py models            # 34 of 59 available on our key
python3 vendor/nebius/forge.py schema stab-ddg --full
```

**`stab-ddg` is StaB-ddG — *binding* ddG, not thermostability.** Single operation
`predict-binding-ddg`; takes a complex PDB, two partner chain sets, and explicit
point mutations:

```json
{"pdb_file": "complex.pdb", "partner_a": ["A"], "partner_b": ["B"],
 "variants": [{"id": "R27P", "mutations": [{"chain": "A", "residue_number": 27,
   "insertion_code": "", "wild_type": "R", "mutant": "P"}]}]}
```

Gotchas from `SKILL.md`:

- **Needs a `.pdb` in a tar**, and our structures are mmCIF — conversion required,
  with chain ids and residue numbers verified against the structure (the template
  does not establish them).
- **~12 min first call**, so send long jobs early and come back to the run id.
- Runs are recorded in `nebius-runs/` in the **current folder**, which is how
  `resume` finds the key — always work from the same directory.
- `422` names the offending field: fix that one field, retry once, then stop.
- **Do not use `openfold2`** — it burns 55 min and fails. Use `openfold3`/`boltz2`.
- `--no-wait` returns the shell; `fetch` verifies checksums.

Secrets: `.env` and `credentials.private.json`, both gitignored. Treat the AMASS
key and the VM SSH key as **exposed** — both were pasted into a chat transcript.
