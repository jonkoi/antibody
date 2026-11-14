# Tasks — Team 15

Working checklist as of **14 Nov 2026, ~11:40 CET**. Derived from
[HANDOFF.md](HANDOFF.md), plus what we learned inspecting this machine.

Strategy and reasoning live in HANDOFF.md — this file is only *what to do next*.

---

## Where we are running

This box (`temp-cpu-server`, europe-west1-b) is **not the hackathon VM**:

| | this box | hackathon VM |
| --- | --- | --- |
| GPU | none (`nvidia-smi` absent) | 1× H100 80GB |
| Docker | not installed | all Foundry images local |
| `abguard` | not installed | in a uv venv, symlinked to `~/.local/bin/abguard` |
| AWS credentials | none | unknown — untested |
| repo | `~/antibody` (this repo) | run artefacts already rescued from it |

**Everything in the pipeline runs on the VM, not here.** This box is for editing,
analysis and the pitch.

Installed here so far: AWS CLI v2.40.3 at `~/.local/bin/aws` (not on `PATH` by
default). `unzip` is absent — extract zips with
`python3 -c "import zipfile;zipfile.ZipFile('x.zip').extractall('.')"`.

---

## Blocker 0 — RESOLVED

The VM is reachable and the data is in hand.

```bash
ssh -i ~/.ssh/team-15 hackathon@89.169.96.144     # user: hackathon, team-15
```

Confirmed on the VM: 1× H100 80GB, `docker`, `aws` (with working credentials for
`s3://apheris-oslo-hackathon/`), and `abguard` at `~/.local/bin/abguard`. The
organisers' data was already synced there at `~/hackathon-data` (4.7 GB), and
`~/rfd4` holds `Data/`, `Settings/`, `Scripts/`, `Tutorial.html` and the numbered
docs `00-start.md` … `18-vm.md`.

**The home directory is shared across the whole team** — same user, same key, one
GPU. Teammates' folders already there: `igor/`, `foundry-assistant/`. Work in your
own folder (e.g. `~/<yourname>/`) so nobody overwrites anybody.

Pulled down to this box by `rsync` over the same SSH key (no AWS credentials
needed here):

- `rfd4-docs/` — the full tutorial, 3 MB, **done**
- `hackathon-data/` — the 4.7 GB data set

Both are gitignored (`.gitignore:26-27`).

**Remember:** `ssh host "cmd"` does not source `.bashrc`, so `claude`, `uv` and
`abguard` look missing. Use `ssh host 'bash -lc "cmd"'`.

### What is in the data

From `hackathon-data/README.md`. Both HER2 tables come from
[Chinery et al.](https://www.biorxiv.org/content/10.1101/2024.03.26.586756v2) —
**cite them when we use the data.**

| Folder | What it holds | Size |
| --- | --- | --- |
| `her2-trastuzumab/` | the two trastuzumab tables — FACS classes + SPR affinities | 6 MB |
| `zenodo-18390239-ab-ag-complexes/` | scores for 1.8M predicted nanobody–antigen complexes | 0.6 GB |
| `sigma2-ab-ag-champloo/` | per-sample ipTM/pLDDT/PAE for those predictions, plus SPR for a few pairs | 1.8 GB |
| `zenodo-10831512-trastuzumab-ngs/` | ready-made train/val/test splits of the FACS data | 2.5 GB |

- `tz_her2_facs_3classes_524346.csv.zip` — **524,346** variants, each a different
  10-residue CDRH3, labelled `high` (172,149) / `medium` (190,011) / `low`
  (162,186) HER2 binding, with `edit_distance` from wild-type `WGGDGFYAMD`.
  Classes say *whether* a variant binds, not how tightly. Trastuzumab itself is
  not in the table.
- `trastuzumab_her2_spr_all_768.csv` — 768 variants measured by SPR; **only 190
  rows have a KD.** Lower KD = tighter binding. This is the affinity ground truth.
- These two tables are the challenge's novelty reference: a design must differ
  from all **525,063** unique CDRH3s in them.
- Note `unzip` is absent on this box — read the zip with Python's `zipfile`, or
  pandas directly.

---

## NEW FINDING — the novelty gate is much harder than arithmetic suggests

**`her2_cdrh3_2` is rejected.** Run on the VM:

```bash
cd ~/team15-claude
abguard team-validate abguard_submission.csv --reference-db her2_novelty.duckdb
```

```
ERROR  submission is not compliant (1 of 3 row(s) pass the schema)
  - [novelty] control_wildtype_trastuzumab: CDRH3 SRWGGDGFYAMDY is 0 mutation(s)
      from SRWGGDGFYAMDY in trastuzumab; a novel design needs at least 3
  - [novelty] her2_cdrh3_2: CDRH3 SRWGPPYFYAMDY is 2 mutation(s)
      from SRWGMPHFYAMDY in tz_her2_facs_3classes_524346; a novel design needs at least 3
```

Design_2 is 3 substitutions from trastuzumab — our arithmetic said it cleared the
floor — but it is only **2 mutations from a FACS variant** we had never compared
against. Only `her2_cdrh3_3` survives. The wild-type control failed as designed,
which confirms the gate is working.

### Why this is structural, not bad luck

The FACS table holds 524,346 variants densely covering the near-wild-type
neighbourhood (1,490 at edit distance 3, 6,392 at 4). **Any conservative design
is likely within 2 mutations of something already measured.** Combine with the
SPR finding below and the trap is explicit:

- affinity needs **few** substitutions (median KD degrades ~1000× by 7 subs)
- novelty needs distance from **525,063 known** CDRH3s, which saturate exactly
  the few-substitution region

**These pull in opposite directions, and that is the real story of this
challenge.** It is a much stronger pitch point than "we made a tight binder".

### What follows — do these first

- [x] **Build the real novelty reference.** Done, on the VM at
      `~/team15-claude/her2_novelty.duckdb`: **525,063 distinct CDRH3s**, exactly
      the number the brief quotes, built in 37s from the two HER2 tables:
      ```bash
      abguard reference-ingest her2_novelty.duckdb \
        ~/hackathon-data/her2-trastuzumab/tz_her2_facs_3classes_524346.csv.zip \
        ~/hackathon-data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv
      abguard reference-build her2_novelty.duckdb
      ```
- [ ] **Novelty-filter BEFORE spending fold compute.** Currently we fold the top
      20, then discover rejections. Reverse it: generate → abguard → fold only
      survivors. This is the single biggest efficiency win available.
- [ ] **Re-check every design we intend to submit** with `team-validate
      --reference-db`. Never trust substitution arithmetic again.
- [ ] Use `--drop` to discard non-compliant rows instead of failing the whole
      table, once we have a real batch.

### Two corrections to HANDOFF

- **"abguard ships its own reference of 610,109 CDRH3s, used by default" is
  wrong in practice.** Without `--reference-db`, prior art is *not* checked — the
  tool warns `--no-novelty: prior art was not checked` and returns
  `is_novel_enough = null`. You must pass `--reference-db` explicitly. Our
  reference has 525,063 entries, not 610,109.
- **`abguard novelty-report` appears broken**: it ignores `--reference-db` and
  emits the `--no-novelty` warning even when that flag is absent, comparing only
  against trastuzumab. **Use `team-validate` for any real verdict.**

### The IMGT window question is answered

HANDOFF asked whether abguard's CDRH3 is the same residues as our design window.
**It is not.** abguard scores `SRWGPPYFYAMDY` (13 residues) while the model
designs `WGPPYFYAMD` (10, B99–108) — IMGT adds `SR` before and `Y` after. Those
three residues are frozen at wild-type, so substitution *counts* are unaffected,
but note `identity_to_trastuzumab` is computed over the whole paired sequence
(98.7% for design_2), not the window. Everything lands in the **optimisation
track** requiring ≥3 mutations; nothing we design this way will ever reach the
de novo track, which needs a different scaffold entirely.

---

## NEW FINDING — the SPR data says affinity is close to unwinnable

Run `python3 pipeline/analyze_spr.py` (new, dependency-free) to reproduce. This
came out of the organisers' own `trastuzumab_her2_spr_all_768.csv` and it should
change what we submit and what we pitch.

**1. The assay's own floor is 6.2×.** Wild-type trastuzumab appears 8 times in
the table, measured at 0.0285 – 0.178 nM — a 6.2-fold spread on an *identical
sequence*. Any affinity "improvement" smaller than that is not distinguishable
from noise. Note the README's claim that trastuzumab is absent applies to the
FACS table; it is very much in the SPR table, and those replicates are a gift.

**2. Affinity collapses monotonically with substitution count.** Median KD by
number of CDRH3 substitutions from wild-type:

| subs | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| median KD (nM) | 0.135 | 4.5 | 20.5 | 15.0 | 60.7 | 92.7 | 93.0 | 146.0 | 248.5 | 190.7 | 203.0 |
| fold worse | — | 33× | 152× | 112× | 451× | 689× | 691× | 1085× | 1848× | 1418× | 1509× |

**3. The novelty gate and affinity are in direct conflict, quantified.** 153 of
the 190 measured variants clear abguard's ≥3-substitution floor. **Only 2 of
those 153 (1.3%) beat wild-type's median KD** — and the better one, `WGSQGFYAYD`
at 0.0234 nM, falls inside wild-type's own replicate spread, so it is arguably
not an improvement at all. The empirical hit rate for "novel *and* tighter than
trastuzumab" is approximately zero.

**4. This condemns our current top-ranked design.** `refold-prepare` ranked
`her2_cdrh3_3` (`WASSSEYYFD`, 7 subs) first, and it had the best ipTM (0.965) and
min_ipae (0.310). But the measured median at 7 substitutions is **146 nM, ~1085×
worse than wild-type**, and the best 7-sub variant ever measured is 2.9 nM — still
20× worse than wild-type. `her2_cdrh3_2` (`WGPPYFYAMD`, 3 subs) sits where the
only successes live. **The confidence metrics and the empirical data disagree,
and the empirical data is measurements.**

### What follows

- [ ] **Re-rank on substitution count, not just confidence.** Prefer 3–4 subs.
      Every design we submit on the affinity axis should be as conservative as
      the novelty gate permits.
- [ ] **Do not claim an affinity win smaller than 6.2×.** Judges who know this
      data will catch it. Framing a result as "within assay noise of wild-type,
      with thermostability improved" is both honest and stronger.
- [ ] **This is the pitch's spine.** It is a quantitative, data-grounded argument
      for why we went after thermostability rather than affinity — made from the
      organisers' own data, not from a hunch. Judgment is 40%.
- [ ] **Mine the 50 variants at edit distance 1–2 and 1,490 at 3** in the FACS
      table (`edit_distance` column) — the sparsely sampled, conservative region
      is exactly where the novelty gate forces us and where the SPR winners are.
- [ ] Cite **Chinery et al., doi:10.1101/2024.03.26.586756** whenever we use
      either table.

---

## The four sources — connection status on this box

HANDOFF records all four as "connected and verified working" — but that was from
the VM/original workstation. **On this box, three of the four are not usable yet.**

| # | Source | What it gives us | Auth | Status here |
| --- | --- | --- | --- | --- |
| 1 | **AMASS** | literature + trials evidence (the NB46 thermostability argument) | `.env` → `AMASS_API_KEY` | **configured + key verified**, needs restart |
| 2 | **JapanFold** | folding | `.env` → `JAPANFOLD_API_KEY` | key **missing**, MCP not configured |
| 3 | **Nebius Scientific AI** | `stab-ddg` thermostability, `protein-scoring` | `credentials.private.json` | file **missing** |
| 4 | **Apheris Foundry** | RFD4-Proteina, Boltz-2 etc. — the whole pipeline | none (local Docker) | **working, VM only** |

Plus a fifth that HANDOFF does not list as a platform: **AWS/S3** for the data.
Credentials exist on the VM and work; this box has none, but does not need them
— we pull over SSH instead.

Full setup instructions, endpoints and re-issue points: **[CONNECT.md](CONNECT.md)**.

- [x] **AMASS MCP configured.** [.mcp.json](.mcp.json) points at
      `https://mcp.amass.tech/mcp` and reads the key from the environment, so no
      secret is committed. Key verified live.
- [ ] **Restart Claude Code to load it** — MCP servers are read at startup:
      ```bash
      cd ~/antibody && set -a && . ./.env && set +a && claude
      ```
      Until then the `search_amass_*` tools do not exist, and
      [summery.md](summery.md) (the AMASS evidence-scout skill) refuses to run
      without them. That skill is what builds the pitch's evidence section.
- [ ] **Get the JapanFold key** from `japanfold.aiand.com/account/token`, add to
      `.env`, add a block to `.mcp.json`. Lowest priority — five folding models
      already run on the VM.
- [ ] **Recover `credentials.private.json`** from `participant-12.zip`. Not on
      this box and **not on the VM either** — it is on a teammate's machine.
      This one matters: Nebius `stab-ddg` is how we score thermostability (task
      5), the strategic differentiator. 34/59 models were available.
- [x] **Reach the VM** — done, see Blocker 0 above.

---

## Critical path

Ordered. Items 1–3 are cheap and unblock judgement; item 4 is long-running and
should be started before pitch work, not after.

### 1. Verify the deadline

- [ ] **Confirm the 18:00 submission deadline with the organisers.** HANDOFF
      records it as assumed, not confirmed. Everything below is scheduled
      against it, so a wrong assumption here invalidates the plan.

### 2. Answer the epitope question — blocks scaling

- [ ] **Decide whether the four hotspot atoms are the right epitope definition.**
      The run aims at HER2 Lys593 NZ, Leu573 CD1, Val591 CG2, Pro572 CD (from
      1N8Z), set in `C_HOT` in [inputs/12-her2/her2_cdrh3.yaml](inputs/12-her2/her2_cdrh3.yaml).
      That set is the *entire* operational definition of "target the trastuzumab
      epitope" and every design inherits it. Needs a structural biologist.
      Write-up: <https://claude.ai/artifact/K1Kd6wCkfJ2v9yUYExp6M1>
- [ ] If the set is wrong, fix `C_HOT` **before** the production batch.

### 3. Run `abguard` on the designs — DONE, and it rejected one

Results and consequences are in the novelty finding at the top of this file.
Artefacts: submission table at
[runs/her2-validation-seed1/abguard_submission.csv](runs/her2-validation-seed1/abguard_submission.csv)
(includes a wild-type control, which correctly failed), reference DB on the VM at
`~/team15-claude/her2_novelty.duckdb`.

Still open from this task:

- [ ] **Run `liability-report`** — we have not looked at LAP sequence
      liabilities at all yet, and they are part of what `team-validate` enforces
      (unpaired cysteines, ≤250 aa, IMGT-numberable). Design_3 passed the schema,
      so nothing is currently blocking, but we have never seen the liability
      summary.
- [ ] **Try `abguard dashboard`** — it launches a marimo liability dashboard in
      the browser. This is the closest thing to a ready-made UI for the data, and
      it is already installed on the VM. Would need an SSH port-forward.

The original questions, now answered:

- [x] Input table built (`kind = paired`; heavy + light from
      `refold_metrics.csv`), with a wild-type control row added.
- [x] Ran on both refolded designs. **design_2 rejected, design_3 passed.**
- [x] **Track: optimisation**, needing ≥3 CDRH3 mutations.
      `identity_to_trastuzumab` is 96.9–98.7%, computed over the whole paired
      sequence, so nothing we produce this way can reach the de novo track.
- [x] **IMGT window differs** — `SR` + our 10 residues + `Y`. Counts unaffected
      because those three are frozen, but do not assume the windows match.
- [x] **A reference DB *is* required** — the opposite of what HANDOFF says. Built
      it; see the novelty finding above.

### 4. Production batch — RUNNING

- [x] **Launched** 14 Nov ~11:55 CET: 100 designs, seed 1, weights
      `her2-novel-2026-11-09`, top 20 to be refolded with Boltz-2. GPU confirmed
      at 98%. Command used:
      ```bash
      cd ~/rfd4 && OUT=runs/team15-claude-n100-seed1 N=100 TOPN=20 SEED=1 \
        setsid nohup bash ~/team15-claude/run_design.sh all \
        > ~/team15-claude/batch-n100.log 2>&1 &
      ```
      Progress: `tail -f ~/team15-claude/batch-n100.log`.
      Output: `~/rfd4/runs/team15-claude-n100-seed1/`.
      `setsid` matters — without it the run dies when the SSH session closes.
- [ ] **When it finishes, novelty-filter the 100 designs before trusting any
      ranking** (see the novelty finding above — we now expect a high rejection
      rate among the conservative designs that are worth submitting).
- [ ] Optionally re-run with a second seed for more candidates in the surviving
      region.

**Shared GPU.** It was idle when we launched, but one other session was logged
in. Check `nvidia-smi` before launching anything else, and keep work inside
`~/team15-claude/` and `~/rfd4/runs/team15-claude-*` so nobody is overwritten.
- [ ] Optionally sweep a second folding model over the same designs for
      independent evidence:
      ```bash
      FOLDER=protenix-v2 bash pipeline/run_design.sh predict score
      ```
- [ ] Analyse with:
      ```bash
      python3 pipeline/analyze_designs.py \
        --designs runs/<run>/designs \
        --reference inputs/12-her2/her2_fv_1n8z.cif \
        --scored runs/<run>/scored/refold_metrics.csv
      ```
      (Verified working against the validation run.)

### 5. Score thermostability — the differentiator

HANDOFF's argument: S-tier needs affinity *and* epitope preservation *and*
thermostability at once; most teams will chase affinity; no thermostability
predictor ships in the Apheris stack, so this is where we can separate.

- [ ] Run Nebius `stab-ddg` on the batch survivors. Available and unused.
- [ ] Consider `protein-scoring` too.
- [ ] Note the format lever: VHH scaffolds are intrinsically more thermostable
      **and** express far better in Adaptyv's cell-free system — expression is an
      unscored gate that zeroes you out entirely.

### 6. Calibrate the confidence cut-off

- [ ] **Do not rank on ipTM or pLDDT.** Validation scored ipTM 0.95–0.97 and
      pLDDT 96–99 because the refold was handed HER2 as a template with
      framework and other CDRs frozen and only 10 residues moving. These will
      not discriminate across a large pool.
- [x] **DONE** — and the premise was wrong. There is no 454-non-binder set:
      `spr-supplementary/` holds 29 and 12 rows, and its README says "N/A denotes
      no reportable affinity fit, NOT confirmed non-binding", so they cannot be
      negatives. Calibrated instead against ~8,100 **shuffled** pairs in
      `zenodo-18390239-ab-ag-complexes/`. Boltz-2 ROC AUC **0.645**; at
      ipTM ≥ 0.70 it gives **94.3% false positives**. See
      [findings/lane-b-confidence.md](findings/lane-b-confidence.md).
- [ ] Rank on `min_ipae` primarily (0.310 vs 0.331 in validation — real spread,
      but n=2 is not a distribution). The `pass` column is **always empty** for
      the HER2 example; the scorer cannot measure absolute movement.

### 7. Assemble the submission

- [ ] **Filter by novelty first, then rank.** Known binders fail the novelty
      filter, so filtering after ranking wastes the slots.
- [ ] Assemble up to 20 ranked designs, split across whatever tracks abguard
      assigns.
- [ ] Pitch placement converts 4 design slots into up to 16, multiplying both
      experimental pools — so the pitch is not separate from the science prize.

---

## Pitch — three of five prize pools

Decided today, and **Judgment alone is 40%.**

- [ ] Build the pitch. Lead with the thermostability argument: it is the
      bottleneck, therefore the differentiator.
- [ ] Use the saturated-confidence finding as evidence of judgement — knowing
      *why* ipTM 0.96 is not a good result is exactly what "Judgment" rewards.
- [ ] Bring the AMASS supporting evidence: nanobody NB46 competes with
      trastuzumab for HER2 *and* retains affinity to 90 °C, so two of the three
      axes are demonstrably satisfiable in VHH format simultaneously.
- [ ] Pools: optimisation (€800/500/200), de novo (€800/500/200), best pitch
      ($500), AMASS best-use, pitch placement.

---

## Not started

- [ ] **No VHH / de novo designs generated.** Only trastuzumab CDRH3 redesign so
      far, despite the de novo pool being worth the same as optimisation. RFD4's
      HER2 weights only redesign CDRH3 — to touch other CDRs or framework, use
      the `antibody-2026-10-30` weights.
- [ ] **Open question 2:** did the designed loop stay on the epitope, or drift? A
      folding model will happily return a confident complex with the paratope on
      the wrong surface patch, and ipTM cannot tell. Refolded structures are in
      `runs/her2-validation-seed1/refold_boltz2/`. A browser structure viewer
      (NGL/Mol*) over those files would make this askable of a structural
      biologist.
- [ ] **Open question 3:** is design_2's double proline (`GDG` → `PPY`) clever or
      a liability? Rigidifying the loop could help affinity and thermostability,
      or distort the paratope.

---

## Housekeeping

- [ ] **`summery.md` is misnamed and probably should not be committed.** It is a
      513-line copy of the AMASS `amass-biomedical-evidence-scout` skill
      definition, not a summary. `.gitignore` lists it under "third-party
      material, kept locally but not redistributed", but commit `83b8058` added
      it anyway, so it is tracked despite the ignore rule. Decide: rename, or
      `git rm --cached summery.md`.
- [ ] **The key pasted into chat is now in that transcript.** `~/.ssh/team-15`
      (ed25519, `oslo-team-15`,
      `SHA256:WUnkJfmqhsQxhA/0RhKD4QfShiGXAsrYDhBw+wW9zT0`) is mode 600 and
      gitignored, but treat it as exposed. The VM dies tonight, so the practical
      risk window is hours — still, do not paste the Nebius, AMASS or JapanFold
      secrets the same way.
- [ ] **The VM shuts down end of Saturday 14 November.** Everything needed is
      already in this repo except the Docker images, which cannot be moved. Pull
      any new run artefacts down **before** shutdown.
