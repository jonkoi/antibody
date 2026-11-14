# Agent work log — log1

Running log of agent sessions on this repo, written so the next agent inherits the
experience rather than rediscovering it. Newest session at the bottom.

**Companion log: [`log0.md`](log0.md)** covers the data analysis and the `pipeline/` code.
**Read both.** This log covers platform reconnaissance, external resources, and the VM
rescue. The two logs were written by two workers running concurrently in this same working
tree — see Session 2 below, which explains how that was discovered.

**Conventions.** Each session gets a `## Session N` heading with date, what was asked,
what was found, what was *not* worth doing, and dead ends. Dead ends are as valuable as
findings — record them.

---

## Session 1 — 14 Nov 2026

### Asked
1. Does the work so far use the HER2-finetuned RFD4-Proteina model from the challenge deck?
2. Go online across all four platforms (AMASS, JapanFold, Nebius Scientific AI, Apheris
   Foundry) and find data or resources that help.
3. Start maintaining this log.

### Q1 answered: yes, and it is the only design model used
Provenance is in the run artifacts, not just the scripts — `ckpt_name =
rfd4-proteina-her2-novel-2026-11-09-ema.ckpt` on all four rows of
`runs/her2-validation-seed1/design_metrics/per_example_metrics.csv`.

The deck's single "HER2 Finetuned" box is really **two** checkpoints. Organisers' own
benchmark table (`ai-antibodies-hackathon-apheris.pdf` p.4–5):

| Checkpoint | Known binders | Pass novelty | Interface err | Used? |
| --- | --- | --- | --- | --- |
| `her2-novel-2026-11-09` | 4% | 86% | 6.8 Å | **yes** |
| `her2-2026-11-07` (conservative/balanced/exploratory) | 21% / 21% / 0.3% | 0% / 6% / 74% | 9.4 / 9.1 / 9.1 Å | no |
| `antibody-2026-10-30` (not finetuned) | 0% | 100% | 8.2 Å | no |
| `release-2026-10-19` (general) | — | — | — | no |

`her2-novel` is the right default — best interface error *and* high novelty pass rate.
Do not spend time on `her2-2026-11-07` except possibly its exploratory setting.
`antibody-2026-10-30` is the one genuinely unused capability that matters: it is required
for the de novo track and for touching anything other than CDRH3.

### Deck deltas (challenge deck is NEWER than the repo's strategy section)
`2026-11_Oslo_hackathon_Apheris_challenge.pdf` was created 14 Nov 10:20 CET, after the
09:50 validation run. It revises numbers HANDOFF.md reasons from:
- Prizes now **$1,500 pot per track, $720/$360/$270**. HANDOFF.md:122 still cites
  €800/500/200 and five pools. Deck shows two tracks splitting evenly.
- Pitch placement concretely = 12/8/4 extra designs for 1st/2nd/3rd + 2 from each of 6
  judges. So 4 slots → at most 16. HANDOFF's claim checks out.
- Slide says optimisation needs "≥2 mutations away from every sequence in the provided
  data". **The slide is imprecise.** Older PDF says "more than 2", i.e. ≥3, which matches
  README and `abguard --help`. Trust abguard, not the slide.

### Repo hygiene finding
`summery.md` is misnamed — it is a verbatim copy of the `amass-biomedical-evidence-scout`
skill, not a project summary. Identical to `amass-biomedical-evidence-scout.skill` and
`.claude/skills/my-skill/amass-biomedical-evidence-scout/SKILL.md`. README.md and
HANDOFF.md are the real handoff documents.

### Platform survey — all four touched

| Platform | Reached? | How |
| --- | --- | --- |
| AMASS | yes | MCP, live literature returns |
| JapanFold | yes | MCP, `list_models` + `get_credits` |
| Nebius Scientific AI | yes | `forge.py models` / `forge.py schema` |
| Apheris Foundry | **no** | local Docker on the hackathon VM only — see below |

**Apheris Foundry is not reachable from this machine.** There is no Docker here, the VM's
IP was deliberately never recorded in the repo (HANDOFF.md says so), and `~/.ssh/known_hosts`
is hashed so the IP cannot be recovered from it. The key is present and valid
(`~/.ssh/team-15`, ED25519, comment `oslo-team-15`). To use Foundry again someone must
supply the IP from the organisers' onboarding mail. The VM shuts down end of today
regardless, so **treat Foundry as a closing window, and anything still needed from it as
urgent**. Everything else below runs without it.

---

## THE HEADLINE FINDING: the thermostability axis is no longer blocked

README.md says "**A thermostability predictor cannot be trained on provided data**" and
"no thermostability predictor ships in the Apheris stack". Both are still literally true.
But the unstated conclusion — that the axis therefore can't be scored — is **wrong**.
There is public experimental Tm data and there are free, downloadable, state-of-the-art
nanobody Tm predictors. HANDOFF.md calls thermostability "the bottleneck and therefore the
differentiator"; this is the lever.

### ⚠️ The one caveat that governs everything here
**Every one of these tools is nanobody/VHH-only.** None applies to trastuzumab's Fab/IgG
format. So this helps the **de novo VHH track only** — and it makes that track
strategically much stronger, which is the same direction HANDOFF.md's strategy paragraph
already argued on format grounds. For the optimisation track the thermostability axis
remains unscored; `score_developability.py`'s mechanistic risk reporting is still the
honest answer there.

### Resource 1 — NanoMelt (recommended; use this first)
State of the art, pip-installable, trained models bundled, no training needed.

- **Install:** `pip install nanomelt` (PyPI, currently **1.5.0** — note the paper is 1.4.0)
- **Run:** `nanomelt predict -i <sequence|fasta> -o out.csv [-align] [-ncpu N]`
- **Needs first:** `conda install -c conda-forge openmm pdbfixer biopython` and
  `conda install -c bioconda anarci`. Python ≥3.8, tested to 3.12.
- **Webserver, no install at all:** <https://www-cohsoftware.ch.cam.ac.uk/>
- **Repo (datasets + notebooks live here, NOT in the PyPI sdist):**
  `https://gitlab.developers.cam.ac.uk/ch/sormanni/nanomelt.git` → `/datasets`, `/notebooks`
- **Paper:** mAbs 2025, <https://doi.org/10.1080/19420862.2024.2442750>, PMID 39772905
- **Architecture:** Ridge stack over GPR(one-hot), GPR(VHSE), SVR(ESM-1b), SVR(ESM-2 t30).
  Trained on 640 apparent-Tm measurements (literature + 129 new in-house).
- **🔑 Embedding note that matters:** the bundled model file is `SVR_esm2_t30.joblib` —
  **ESM-2 t30 is the 150M model, not 15B.** It generates its own embeddings via `fair-esm`.
  So there is no large-model dependency and no need for Nebius embeddings at all.
- **LICENSE: CC BY-NC-SA 4.0 — non-commercial.** Fine for a hackathon entry, but it is a
  real restriction and should be **disclosed in the pitch** rather than quietly used.
- Verified contents by downloading the sdist and listing it (did not execute it). Bundled:
  `GPR_onehot.joblib`, `GPR_vhse.joblib`, `SVR_esm1b.joblib`, `SVR_esm2_t30.joblib`,
  `ridge_ensemb.joblib`. Also ships `modules/alignment/liabilities.py`, which overlaps
  `pipeline/score_developability.py` — worth diffing before writing more liability code.

### Resource 2 — NbThermo (the underlying experimental dataset)
- Repo: <https://github.com/Valdes-Tresanco-MS/NbThermo> · Site:
  <https://valdes-tresanco-ms.github.io/NbThermo> (Angular SPA — **WebFetch gets an empty
  shell, needs a real browser or the repo**; don't waste another fetch on it)
- Paper: PMID 37042467, <https://doi.org/10.1093/database/baad021>
- 564 nanobodies from 65 publications. **Tm range 39–98 °C, median 67.2 °C.** Aho numbering
  via ANARCI. Also carries **KD in nM and expression yield in mg/l** — so affinity and
  expression, not just Tm.
- **Finding that cuts against our README:** comparing consensus frameworks of Tm ≤60 °C vs
  Tm ≥75 °C nanobodies found only 2–4 differing positions, and those residues appear in
  both groups. Their conclusion: "the Fw structure alone does not define Nb thermostability,
  implying that **the CDRs make an essential contribution**", with the CDR3–framework
  interface specifically implicated (Kinoshita et al.). README.md asserts "CDRH3 edits
  rarely improve Tm" — NbThermo's data says CDR3 is in fact a major Tm determinant. That
  reframes CDR3 redesign from "can't help Tm" to "is one of the main handles on Tm".
- Also: CDR3 length correlates only **poorly** with Tm (so don't use length as a proxy),
  and Tm vs refolding ratio correlates moderately (r ≈ 0.63).

### Resource 3 — TEMPRO (alternative; useful mainly for its dataset)
- Repo: <https://github.com/Jerome-Alvarez/TEMPRO> (the guessed `sdean-lab/TEMPRO` 404s —
  do not retry it). Paper PMID 39154093, DOI 10.1038/s41598-024-70101-6.
- `sdab_data.xlsx` = the processed 567-sequence set with Tm. Repo also has the embedding
  generator, pre-trained models, and the AlphaFold2 PDBs.
- DNN on ESM-2 embeddings. **MAE 4.03 °C / RMSE 5.66 °C / R² 71.4 %** with ESM-2 **15B**;
  **ESM-2 3B gives MAE 4.24 °C**, i.e. nearly identical. Nebius has `esm2-3b`, so TEMPRO is
  reproducible on our stack if we ever want to retrain. NanoMelt is still the better
  ready-made choice.
- Negative result worth keeping: AlphaFold2 per-region pLDDT correlates with Tm at only
  ρ ≈ 0.24 (best region, FR3). **Do not expect pLDDT to predict Tm.** Physicochemical
  features and NetSurfP3 also underperformed embeddings.

### Resource 4 — NBsTem (newest, strongest reported numbers)
- PMID 41564239, ACS Appl. Mater. Interfaces, Jan 2026. Dual model: `NBsTem_Tm` on 514
  experimental Tm + `NBsTem_Q` on an MD-derived indicator over 704 structures.
- **Pearson 0.83 on an external test**, reported to beat three competitors. Web server
  ("NBsTem"). Proposes a concrete screening gate: **Tm > 65 °C and Qclass IV**.
- Handles nanobodies with missing residues. Not yet verified as downloadable — only the
  paper was read.

### Other Tm work found, not yet chased
`NbBayesLM` (Bayesian PLM, PMC13272114); sparse-autoencoder interpretability of fine-tuned
PLMs for nanobody Tm (Bioinformatics 2026, repo `matsunagalab/paper_nanobody-thermostability-sae`,
reuses TEMPRO's split and seed).

### Engineering moves from the literature (design levers, not just scoring)
- **CDR grafting onto a thermostable framework** (PMID 40613834): grafting Nb20's CDRs onto
  the A4.2m framework gave **+10 °C Tm** and +55 % refolding; further mutations reached
  68 °C Tm and >82 % refolding. A directly usable de-novo strategy — pick a high-Tm
  framework, graft the designed CDRs.
- **Disulfides** (PMID 29784954, Kunz et al., 195 citations): an optional second disulfide
  suppresses aggregation on heat denaturation; up to three disulfides reported giving up to
  **+19 °C**. Counterintuitive bonus: *lengthening* a CDR loop contributed to stability.
  ⚠️ But abguard rejects **unpaired** cysteines, and nanobodies with two disulfides lose
  function in the cytosol — so any added cysteine must be a designed, paired bond.
- **NB46** (PMID 42911746, J Med Chem, Oct 2026, JuFo 3): retains HER2 affinity to **90 °C**
  *and* competitively blocks trastuzumab. Already cited in HANDOFF.md; this is the
  existence proof that high Tm + trastuzumab-epitope binding co-occur in VHH format.
- **Trastuzumab-epitope VHH with sub-nanomolar affinity** (PMID 37632647) — sequence is
  published. Useful as a positive control, but **it is a known binder, so it will fail the
  novelty gate**. Use it to calibrate, never submit it.

---

## JapanFold: the capability gaps in the Apheris stack, and they are cheap

**Balance: $1,100** (= ~4,230 processor-hours at $0.26/h). This is effectively unlimited for
our purposes and is the most underused asset we have.

What JapanFold has that Apheris Foundry does not:

1. **Boltz-2 with real affinity prediction.** Foundry's Boltz-2 is refold/confidence only;
   JapanFold's `boltz2` lists `affinity` in its caps. This is a direct read on the affinity
   axis. ~$0.0045 per target (p90).
2. **AF2-IG (`af2ig`) — the discriminator we actually need.** Described as "the filter a
   binder-design pipeline uses to tell a real design from a plausible one". Takes a design
   you already have (target structure + binder sequence) and re-predicts the interface,
   reporting ipTM and interface pAE. Deterministic: no MSA, no sampling, same input → same
   answer. **$0.00035–0.00081 per target.** Scoring 1,000 designs costs under a dollar.
   - 🔑 **With `write_pae: true` it returns the full PAE matrix and the docs say this is
     "the one Adaptyv's ipSAE pipeline reads."** Adaptyv is the lab that will test the
     submissions. Ranking on the metric the testing lab itself computes is the single
     highest-leverage thing on this list, and it directly fixes the saturated-confidence
     problem HANDOFF.md flags (iptm 0.95–0.97 across the board, no discrimination).
3. **OpenDDE · Antibody-Antigen (`opendde-abag`)** — an antibody-antigen-specific
   checkpoint on the Protenix-v2 stack, 13,921 jobs run. A genuinely independent second
   opinion on our complexes, specialised for exactly our problem shape.
4. **`use_tfg_guidance`** — "contact pairs or **epitope residues the antibody must reach**",
   enforced as rigid-body restraints during diffusion. This speaks directly to **open
   question 2** (did the designed loop stay on the epitope or drift?): we can fold with the
   epitope enforced and compare against folding it free.
5. **BoltzGen `nanobody-anything` protocol** — ranked, sequenced VHH binders with
   confidence metrics, for the de novo track. Also `antibody-anything`, `protein-redesign`.
6. **ESMC (300M/600M/6B) and SaProt (650M/1.3B) embeddings.** SaProt's card says
   "especially for stability/function prediction" — a candidate feature set if we retrain
   a Tm model rather than using NanoMelt off the shelf.

Hard limits to design around: `max_residues` 1024 (AF2-IG, RFD4, BoltzGen), 1920 for
Boltz-2; `max_binder_length` 200; `max_rfd4_designs` 8 per job; `max_runtime_design_s` 2700.

## Nebius: 34 of 59 models usable on our key
Relevant and unused:
- **`stab-ddg`** — schema read this session. `operation: predict-binding-ddg`, takes
  `pdb_file`, `partner_a`/`partner_b` chain lists, and a `variants` array of
  `{chain, residue_number, insertion_code, wild_type, mutant}`. For us: partner_a = A(light)
  + B(heavy), partner_b = C + D (HER2). **Confirms README's point — this is *binding* ΔΔG,
  i.e. the affinity axis, not folding stability.** Source pinned to
  `github.com/LDeng0205/StaB-ddG`.
- **`antifold`** — antibody-specific CDR inverse folding. An independent sequence proposer
  for a fixed RFD4 backbone; a cheap way to diversify CDRH3 at the novelty floor.
- **`nanobodybuilder2`** / **`abodybuilder3`** — fast VHH / paired-Fv structures.
- **`esm2-3b`**, **`esmc-*`**, **`saprot-*`**, **`igbert`** (antibody-specific embeddings).
- **`bioemu2`** — conformational ensembles. Possible route to a *physical* loop-rigidity
  readout for design_2's double proline (open question 3), rather than a score.
- MD: `gromacs`, `amber`, `namd`, `openmm`, plus `mdprecheck`.

⚠️ `openfold2` is flagged UNRELIABLE as of 12 Nov (endpoint fails to start) — don't build on it.

---

## Recommended next actions, revised in light of the above
Ordered by value per unit effort. Items 1–3 are new; the rest refine HANDOFF.md's list.

1. **Score every design with JapanFold AF2-IG using `write_pae: true`.** Costs cents, is
   deterministic, and yields the PAE matrix Adaptyv's own ipSAE pipeline consumes. Use it
   as the primary ranking signal instead of saturated Foundry iptm/pLDDT. Calibrate the
   cut-off on the 454 SPR-measured non-binders, which is HANDOFF action 6 — AF2-IG makes it
   affordable.
2. **If pursuing the de novo VHH track, score candidates with NanoMelt.** Either the
   webserver (zero install) or `pip install nanomelt`. This turns thermostability from an
   unscored axis into a ranked one, for the one track where it applies. Disclose the
   non-commercial license.
3. **Re-examine the "CDRH3 edits rarely improve Tm" premise.** NbThermo's framework/CDR
   analysis argues the opposite. This is a pitch-level argument, not just an implementation
   detail — Judgment is 40 % of the pitch score.
4. Verify the deadline (unchanged, still unconfirmed).
5. Answer open question 1 (epitope hotspot definition) — still blocks scaling.
6. Run `abguard` on the two validation designs — still the cheapest high-value test.
7. Production batch, if the VM is still alive. **This is the only item that needs Foundry,
   so it is the only one with a hard deadline today.**

## Dead ends and traps — do not repeat these
- `github.com/Sormannilab/NanoMelt`, `github.com/Sormannilab`, `github.com/sdean-lab/TEMPRO`,
  `github.com/sdean-lab` — all permanent 404. Guessed names. The real ones are in the log above.
- The NbThermo site is a client-rendered Angular app; WebFetch returns a 9-byte `" FrontEnd"`
  shell. Use the GitHub repo instead.
- The PyPI "source" convenience URL
  (`files.pythonhosted.org/packages/source/n/nanomelt/...`) silently returns **0 bytes**.
  Get the real hashed URL from `https://pypi.org/pypi/nanomelt/json` → `urls[].url`.
- `Read`ing `summery.md` returns the amass skill, not a summary — see repo hygiene above.
- `mcp__amass__get_amass_biomedcore_record` takes `{type, value}` (type one of
  `amassId|pmid|doi`), **not** `record_id`, and `includeFulltext` is camelCase.
- A `cd`-prefixed Bash call can get its sandbox verdict invalidated if the shell moved since;
  the fix is to re-issue the identical call, which succeeds.

---

## Session 2 — 14 Nov 2026, ~12:45–13:10 CET

### Asked
The user supplied the VM connection details that Session 1 said were missing, then asked
what was going on. So: get onto Apheris Foundry, the fourth platform, before the VM dies.

### VM access — resolved, record this
```
IP   89.169.96.144     user  hackathon     key  ~/.ssh/team-15   (ED25519, "oslo-team-15")
ssh -i ~/.ssh/team-15 hackathon@89.169.96.144
```
Host `computeinstance-e00y4btfkg7hzqhhcb`. **H100 80GB completely idle** (0 MiB used, no
compute processes), 354 GB free of 484 GB. The ssh banner warns about non-post-quantum key
exchange — cosmetic, ignore it.

Mechanics that cost time:
- **`rsync`/`scp` do not expand a remote `~` inside quotes.** `hackathon@IP:'~/dir'` fails
  with "No such file or directory" and rsync exits 23 / broken pipe. Use absolute
  `/home/hackathon/...`. This is the single thing that wasted the most time this session.
- HANDOFF's gotcha is right: wrap remote commands as `ssh host 'bash -lc "cmd"'` or
  `abguard`, `uv` etc. appear missing.
- `rsync --info=progress2` emits ~30 KB of carriage-return spam into the tool result. Use
  `--info=stats2` and `| tail`.

### 🔴 THE FINDING: there is a second worker in this same working tree
At session start `HEAD` was `779ebb4`. Partway through, `HEAD` had moved to **`214499a`
"Enumerate the legal design space; calibrate interface confidence"**, and files I never
created had appeared (`data/legal_3mutants.csv` tracked; `log0.md`,
`pipeline/rank_candidates.py`, `pipeline/analyze_replica_variance.py` untracked;
`pipeline/run_design.sh` modified).

I initially mistook `data/legal_3mutants.csv` for something my own rsync had pulled down.
It was not. **If you see unexplained files or a moved HEAD, check `git log` before
assuming your own tool calls did it.** Consequences for anyone working here:
- **Do not `git checkout`/`reset`/`clean`, and do not force-add broadly.** Another worker's
  uncommitted edits are live in this tree.
- Read `log0.md` before acting; it is the other worker's log and it is substantive.
- Co-ordinate through the logs, which is presumably why the user asked for `log1.md` by name.

### Foundry work: the production batch had already been run, and was unrescued
HANDOFF.md says "**No production batch**" and lists launching one as next-action 4. That is
**out of date**. On the VM:

| path | state |
| --- | --- |
| `~/rfd4/runs/her2` | the 4-design validation run (09:51), already in the repo |
| `~/rfd4/runs/team15-claude-n100-seed1` | **100 designs, complete, 128 MB (11:30)** |

The batch has `design/generation/inference_metrics/per_example_metrics.csv`,
`prepared/selection.csv`, `scored-boltz2/refold_metrics.csv`, Boltz-2 refolds for 100
designs, **plus a second `subset/` pass** with its own prepared/selection and
scored-boltz2/refold_metrics. Verified 100 directories under
`design/generation/inference_samples`.

**Rescued into the repo** (all verified present):

| destination | size | what |
| --- | --- | --- |
| `runs/team15-claude-n100-seed1/` | 128 M | the production batch, complete |
| `analysis/analysis-iptm/` | — | `calibrate_iptm.py`, `out.txt`, `summary.csv`, and the three 91×91 ipTM matrices |
| `analysis/analysis-variance/` | 4.5 M | `analyze_replica_variance.py`, `band.py`, `mat.py`, `peek.py`, `top.py`, per-system RMSD/ipTM CSVs, probe logs |
| `analysis/enumerate_candidates.py` | 7 K | the legal-set enumerator |
| `analysis/her2-sanity.log`, `analysis/data-sync.log` | — | run logs |
| `data/her2-trastuzumab/` | 6 M | `tz_her2_facs_3classes_524346.csv.zip`, `trastuzumab_her2_spr_all_768.csv` (gitignored) |
| `data/spr-supplementary/` | — | `spr_table_s1_pairwise.csv`, `spr_crossreactivity_14target.csv`, 2 READMEs |
| `data/README.upstream.md` | 10 K | organisers' dataset README (gitignored) |

`.gitignore` already anticipated `data/her2-trastuzumab/` and `data/README.upstream.md`, so
the layout I chose matches what the other worker had planned.

### Still on the VM, deliberately not rescued
Re-downloadable from public sources — record the identifiers, not the bytes:
- `hackathon-data/zenodo-10831512-trastuzumab-ngs/data.zip` — **2.5 GB**, Zenodo 10831512.
  log0 establishes its 42 `.h5` files are not usable; don't re-pull it for those. It does
  hold `data/proteinmpnn_output/design_3000H3s_temp0pt3.fa` (3,000 ProteinMPNN VH designs),
  which is the one genuinely interesting thing in it.
- `hackathon-data/sigma2-ab-ag-champloo/` — **1.3 GB + 487 MB**, Sigma2 archive
  `10.11582/2026.owjnnp5q`. The ipTM matrices we actually needed are already extracted into
  `analysis/analysis-iptm/data/`.
- `hackathon-data/zenodo-18390239-ab-ag-complexes/` — Zenodo 18390239, 569 MB + 16 MB.
- `~/igor/` (4.7 GB) — an S3 mirror of `s3://apheris-oslo-hackathon/rfd4-proteina/docs/`,
  same content as `~/rfd4/`. Nothing original.
- `~/analysis-models/` (82 MB) — `data/affinity_data/{her2,influenza,example}`. **Not
  rescued and not yet inspected. If anything is still wanted off the VM, look here first.**

All of `s3://apheris-oslo-hackathon/` is reachable from the VM with `aws s3 ls` (the VM has
`.aws` credentials). Those credentials are VM-local and die with it.

### `abguard` is present and working
`/home/hackathon/.local/bin/abguard` — runs, prints usage. **Still never run on a design**,
which is now next-action 4 in log0 *and* item 6 in log1's Session 1 list. It is the cheapest
unresolved question in the project and it needs the VM, so it is on today's clock.

### Convergence worth noting between the two logs
log0's finding 3 and log1's Session 1 recommendation were reached independently and agree:

- log0 measured, on the organisers' own 91×91 benchmark, that **Boltz-2 ipTM has ROC AUC
  0.645, passes 94.3% of known non-binders at a 0.7 threshold, and its top-10 ranked
  contains zero real binders**, while **AlphaFold3 reaches AUC 0.754 and 51.8% precision in
  its top 10**.
- log1 independently recommended JapanFold's **AF2-IG**, an AlphaFold2 initial-guess filter,
  on the grounds that it is deterministic, costs ~$0.0005/design against a $1,100 balance,
  and emits the PAE matrix Adaptyv's own ipSAE pipeline reads.

Same conclusion from two directions: **the AF2/AF3 family discriminates binders and Boltz-2
does not.** The Foundry pipeline ranks on Boltz-2 `min_ipae`/`iptm`, so its ranking of the
100-design batch should be treated as unreliable until re-scored. That makes re-scoring
`runs/team15-claude-n100-seed1` with AF2-IG the highest-value compute action available, and
it does **not** need the VM — JapanFold is hosted, so it survives the shutdown.

### Revised urgent list — what actually dies with the VM today
Only two things genuinely require Foundry:
1. **Run `abguard`** on the top candidates from `data/legal_3mutants_ranked.csv` and on the
   batch designs. Resolves the novelty rule, the 10-mer/11-mer/IMGT registration question
   (log0, "Worth knowing"), and which track each design lands in.
2. **Anything else wanted out of `~/analysis-models/`** or `s3://apheris-oslo-hackathon/`.

Everything else — AF2-IG re-scoring, NanoMelt, Nebius `stab-ddg`/`antifold`, all AMASS work —
is hosted and unaffected by the shutdown. The 100-design batch is already safe in the repo.

---

## Session 3 — 14 Nov 2026, ~13:10–13:45 CET

### Asked
Inherit log1 and find state-of-the-art deep-learning methods for antibody design against a
HER2-like antigen.

### Where the output went
**[`docs/sota-antibody-design.md`](docs/sota-antibody-design.md)** — the survey itself:
calibration campaigns, method tables for de novo and optimisation, the ranking/filtering
literature, a stack-availability map, and the open questions. This entry records what it
changes and what it corrects. Read the survey before acting on any of the below.

### 🔴 THE HEADLINE FINDING: the Adaptyv competition was won on the framework, not the CDR
Adaptyv Bio — **the lab that will test our designs** — ran a crowdsourced EGFR binder
competition and published the post-mortem (bioRxiv 10.1101/2025.04.17.648362). It is the
closest published analogue to this challenge that exists.

**The R2 winner (Cradle, KD 1.21 nM, 8.2× better than the cetuximab control) took cetuximab,
reformatted it as an scFv, and mutated only framework residues — keeping the CDRs.** The best
fully de novo design in the same round was 82 nM. Only ~10 of ~58 validated binders across
both rounds came from de novo/hallucination methods, and the hallucinated share of true
binders *fell* from 40% to 12% between rounds.

Independently, `log1` Session 1's thermostability search had already turned up the other half
of this: **Ng et al., mAbs 2025 (PMID 40653772) designed and experimentally validated a
stabilising *framework* mutation on trastuzumab itself — VH `S85N+R87T` — which raises Tm by
DSF and retains both HER2 binding and ADCC.** Rosetta `relax` + `pmut_scan` on PDB 7mn8.

**Why this is actionable for us:** abguard scores novelty on CDRH3 only (HANDOFF.md:163), so
framework edits are orthogonal to the novelty gate. The play is therefore: spend the mandatory
≥3 CDRH3 substitutions as cheaply as possible — `log0`'s distance-3 shell, the 90%-binding
tier — and take the *differentiation* from a framework mutation with published experimental Tm
and function data. That is an A-tier design by construction (epitope preserved, Tm improved,
affinity retained), which is exactly the strategy `log0`'s next-action 7 argues for, and it
closes the axis `README.md` calls unscorable for Fab format.

**Three honest blockers, all unresolved:**
1. **Whether the rules permit framework edits at all is unconfirmed.** The CDRH3-only window
   in our pipeline is a *checkpoint* limitation (the HER2 RFD4 weights only redesign CDRH3 —
   HANDOFF.md:165), not a stated rule. Nobody has checked the submission rules on this.
   **This is the single question whose answer unlocks or kills the finding.**
   *Update, same session:* the other worker's `data/abguard/novelty_AB.csv` (written 13:15)
   **empirically settles the abguard half of this.** Its schema is
   `cdrh3, cdrh3_confidence, identity_to_trastuzumab, novelty_class,
   required_cdrh3_mutations, cdrh3_min_distance, …` — novelty is scored on **CDRH3 distance
   alone**, and abguard takes a full `heavy_sequence`, so a framework-mutated heavy chain is
   accepted input. Track assignment keys on `identity_to_trastuzumab`, which VH `S85N+R87T`
   moves from 1.0 to ≈0.991 (2 substitutions in ~227 variable-domain residues) — nowhere near
   the 0.70 optimisation/de-novo boundary. **So framework edits are confirmed orthogonal to
   both the novelty gate and the track assignment.** What remains is purely a *rules* question
   for the organisers, not a tooling question. Nobody has yet run abguard on a
   framework-mutated probe; that is a one-row test.
2. **`S85N+R87T` introduces a novel N-linked glycosylation sequon.** The paper's own authors
   cannot exclude that ectopic glycans drive their result. Worse, our own
   `pipeline/score_developability.py:50` penalises `N[^P][ST]` at weight **4.0** — the
   heaviest liability it scores. The one validated stabilising mutation for our exact molecule
   trips our own worst flag. Disclose this; do not hide it.
3. The actual ΔTm is in a DSF figure panel, not the text. Magnitude unknown.

### 🔴 This challenges log1 Session 1's recommendation #1 (rank on AF2-IG ipSAE)
Session 1 called ranking on AF2-IG + ipSAE "the single highest-leverage thing on this list".
The new evidence is mixed rather than confirming:

- **For:** ipSAE is genuinely what Adaptyv consumes — JapanFold's `write_pae` help text says
  the PAE matrix is "the one Adaptyv's ipSAE pipeline reads", and Adaptyv selected **600 of
  ~1,196 Nipah designs by best Boltz-2 ipSAE score** (+400 expert-curated, +200 community vote).
- **Against:** Bits to Binders (bioRxiv 10.64898/2026.03.03.709355, 12,000 designs) measured
  ipSAE directly and found it — with LIS, SAP, Rosetta InterfaceAnalyzer, pDockQ and
  ProteinMPNN likelihood — **"lacked substantial predictive power."** Adaptyv themselves noted
  ipSAE **reproducibility varies across hardware and seeds**, and added expert-selected designs
  to mitigate it. And in the EGFR competition the eventual winner was ranked **322nd–623rd** by
  the organisers' own in-silico metrics.
- **Resolution, and it is a narrowing not a retraction:** AF2-IG remains the right *primary*
  signal for our task — deterministic, $0.00084/target, and `log0` measured the AF2/AF3 family
  as the only one that discriminates antibody binders (AUC 0.754 vs Boltz-2's 0.645) on a
  benchmark with real non-binders, which is the most task-relevant evidence anyone has. But
  rank on an **ensemble**, keep **diversity** in the submission set, and never let a metric cut
  a design the measured FACS/SPR data likes. Bits to Binders is a CAR-T proliferation assay on
  80-aa miniproteins, so its negative result transfers only partly — but it transfers.

### 🔑 The one genuinely new, cheap, orthogonal filter: amino-acid composition
Two independent competitions, different targets, different assays, same answer —
**Glu/Lys composition beats every structure-confidence metric at predicting whether a design
survives to be measured at all:**

| Signal | Adaptyv EGFR (AUC, expression) | Bits to Binders (AUC, recovery) |
| --- | --- | --- |
| Glu content / K+E α-helix fraction | **0.77** | **0.91** |
| Lys content | **0.73** | — |
| DNA entropy / AA repeats / EE repeats | — | **0.90 / 0.88 / 0.82** |
| **ipTM** | **0.58** | "did not predict" |
| **pLDDT** | **0.55** | "did not predict" |
| ESM-2 PLL | **0.44** | 0.85, but attributed to entropy/repeat artefacts |

Bits to Binders' retrospective filter lifted recovery 56.8% → 81.2%; adding "no cysteines"
lifted functional rate 5.9% → 10.6%. This is pure-stdlib computable, costs nothing, and is
orthogonal to everything in `pipeline/` today. Their one-line summary: **"design for
expression, not just binding."**

⚠️ **This cuts against Session 1's disulfide thermostability lever.** Session 1 recorded
engineered second disulfides giving up to +19 °C. Bits to Binders found cysteine-containing
designs *depleted*, and says the result "validates the cautionary choice many model developers
have made in explicitly avoiding cysteines." Both can hold — a designed, paired bond is not a
stray Cys, and abguard rejects unpaired Cys anyway — but the prior on "add a disulfide" should
drop.

### 🔴 Antibody ΔΔG prediction is not fit for ranking — drop `stab-ddg` as a signal
Three 2025–26 papers, independently:
- **Graphinity** (Nat Comp Sci 2025, PMID 40629090): reaches Pearson 0.87 but "appears to be
  overtraining on the few hundred experimental data points available and performance is not
  robust to train-test cut-offs"; concludes **orders of magnitude more data** are needed.
- **AbAgym** (mAbs 2025, PMID 41277441): 324k mutations, 68 DMS datasets — force fields
  "modest", ML models **"only marginally better than random."**
- **MochiBind** (iScience 2026, PMID 42564542): existing predictors "lack generalization …
  deteriorating for antibodies targeting antigens absent from training data or datasets lacking
  non-binders" — our exact situation.

`log0` already noted `stab-ddg` is *binding* ΔΔG, not folding. This is the stronger objection:
the whole class is unreliable. **The cautionary HER2 case is published too:** NAZ-mAb variants
(PMID 42909185, Oct 2026, Kyte/Greiff, **Oslo**) ranked by FoldX interaction energy on a
p95HER2 antibody came back **2–6× worse** by Biacore (326 nM / 94.5 nM vs 53.3 nM parent)
despite better ELISA signal. Same city, same antigen family, one month old. Excellent slide.

### Prior art we did not know about — matters for the pitch
- **DualGPT-AB** (Nat Computational Science, Apr 2026, PMID 41986730) is the closest published
  analogue to our optimisation track: dual-stage conditional GPT + RL for multi-property CDRH3
  generation, **against HER2**, where **8 of 100 sampled designs bound HER2 well and wet-lab
  validation showed enhanced tumoricidal activity vs Herceptin.** No code findable. Cite it as
  the benchmark and the proof the problem is tractable; we cannot run it.
- **AB-Gen** (PMID 37364719) did HER2 CDRH3 library design with GPT+RL in 2023; code on Zenodo.
- **The de novo affinity ceiling is real and quantified.** Best first-round de novo
  antibody/VHH affinities with no experimental maturation: RFantibody tens–hundreds nM
  (single-digit nM *only* after OrthoRep), Chai-2 82 nM, BoltzGen nanobodies 13–14.4 nM,
  Germinal "nanomolar". **Trastuzumab is 0.146 nM.** No published de novo antibody campaign has
  come within two orders of magnitude in one round. The de novo track cannot be won on affinity
  against a trastuzumab benchmark — only on novelty, epitope and Tm. This independently
  supports `log0`'s "bank A-tier, don't gamble on affinity" conclusion.
- **Published negative result that saves us effort:** antibody-specific language models
  (AntiBERTy, AbLang2, AntiBERTa2) carry **no framework stability signal** — they just recover
  the germline/WT residue with high probability. Generic ESM-2 and structure-based Rosetta were
  both more informative. Do not reach for an AbLM for stability.
- **ESM-2 log-probabilities did not correlate with KD at all** in the EGFR competition; **ESM3
  and ESM-C did** correlate strongly with binding category. If we use a PLM likelihood it should
  be ESM-C. ⚠️ But JapanFold exposes ESM-C as *embeddings* only (mean/max/cls pooled), not a
  pseudo-log-likelihood — the LM head is not exposed, so this is not a drop-in.
- **pDockQ/mpDockQ and interface area** were the only metrics relevant across all molecule
  types in the EGFR ANOVA (F ≈ 9.9–19.4), beating ipTM/iPAE. Both are computable from the PAE
  matrix AF2-IG already returns.

### Corrections to earlier claims in this repo
- **`log1` Session 1: "max_residues 1024 (AF2-IG, RFD4, BoltzGen)" is wrong for RFD4.**
  `list_models` gives RFD4 `max_residues` **512**. 1024 is the global `reference_residues` and
  holds for AF2-IG and BoltzGen; RFdiffusion3 and PXDesign allow 1536, Boltz-2 1920.
- **`log1` Session 1 implied AF2-IG is the model Adaptyv's ipSAE pipeline reads.** ipSAE is
  computed from *any* model's PAE matrix; Adaptyv ranked the Nipah competition on **Boltz-2**
  ipSAE. AF2/AF3-family PAE is still the better-founded input per `log0`'s benchmark, but the
  "Adaptyv uses AF2-IG" reading is not supported.
- **PDB 7mn8** (cryo-EM trastuzumab Fab–HER2, 3.45 Å; + 8q6j to fill HER2's missing residues)
  is an **independent second structure of the interface we define from 1N8Z**. `log0`
  next-action 3 wants a structural biologist to check our four hotspot atoms — superposing
  7mn8 on 1N8Z answers much of that computationally, today, for free.

### Capability notes from `list_models` (free call, starts nothing)
- BoltzGen protocols include **`nanobody-anything`** (VHH), **`antibody-anything`**, and
  **`protein-redesign`** = "re-design residues of an existing binder" — the last is directly
  relevant to the optimisation track and nobody has tried it. ~$0.04/job p90.
- RFD4's card states upstream uses **`high/front_loaded`** for antibody/nanobody CDRs.
  `pipeline/run_design.sh:24` already sets this. No action — recording so nobody "fixes" it.
- **OpenDDE · Antibody-Antigen** has 13,671 jobs behind it and is an antibody-specific
  second opinion at $0.011/target.

### Dead ends and traps — do not repeat these
- **Chai-2 is the best-reported de novo antibody method and is unusable.** 15.5% hit rate over
  52 targets (20.0% VHH, 13.7% scFv) at ≤20 designs each, but the paper states **no** API, no
  weights, no licence. It also never verified epitope specificity ("competitive binding assays
  and laboratory 3D structure determination should be performed"), never mentions HER2, and
  uses trastuzumab only as a polyreactivity control. Treat 15.5% as the field ceiling, not an
  option. Don't spend time looking for access.
- **The "BoltzGen 3% / RFdiffusion 22% / Mosaic 89%" Nipah per-method hit rates are
  unverified and should not be repeated.** The competition pages give design counts only
  (BoltzGen 121, BindCraft 30, RFdiffusion 19, Germinal 14) with no binder counts. The real
  numbers are in the collection download, not on any page.
- **`adaptyvbio.com/blog/po104/` is about EGFR Round 2, not Nipah**, despite surfacing in
  Nipah searches. Its 14% / BindCraft 9% figures are EGFR figures.
- **No public code for DualGPT-AB** despite it being the closest prior art. Checked PyPI-style
  searches, GitHub-by-author, and the publisher page; only secondary coverage.
- **Germinal is open (`SantiagoMille/germinal`) but licence-encumbered** — depends on PyRosetta
  and IgLM, both restricted. `OpenGerminal` is an Apache-2.0 reimplementation that swaps them
  out. Check which you need before building on it.
- PubMed returns **HTTP 203 + a cookie interstitial** to automated fetches; use the AMASS MCP
  or bioRxiv directly. bioRxiv's `/content/<doi>v1.full` HTML works where the PDF is awkward.

### Revised next actions, merged with log0's list
Only the ones this session changes; `log0`'s 1–7 otherwise stand.

1. **Confirm with the organisers whether the rules permit framework mutations.** The
   cheapest high-value question in the project right now — it unlocks or kills the headline
   finding above. Framework edits are already confirmed orthogonal to abguard; check the
   rules and run `abguard` on a trastuzumab variant carrying VH `S85N+R87T` plus a legal
   distance-3 CDRH3 to see what track it lands in and whether it passes.
2. **Add the composition/expression filter** (Glu/Lys fraction, K+E helix fraction, repeat and
   entropy checks, cysteine count) to `pipeline/score_developability.py`. Stdlib, free,
   orthogonal, and twice-replicated against wet-lab outcomes.
3. **Superpose 7mn8 on 1N8Z** to validate the four-hotspot epitope definition — partially
   answers `log0` next-action 3 without needing a structural biologist.
4. **Drop `stab-ddg` from any ranking plan.** Keep it as a reported number if you like, never
   as a selection criterion.
5. **If the de novo track is pursued, use BoltzGen `nanobody-anything` on JapanFold** — the
   only open, wet-lab-validated, already-accessible de novo VHH method on our stack. Expect
   10–20 nM, not sub-nM; design the pitch accordingly.
6. **Consider BoltzGen `protein-redesign`** for the optimisation track. Untried, cheap, and
   the only generative route we have that is built for "improve this existing binder".

---

## Session 4 — 14 Nov 2026, ~14:05–14:15 CET

### Asked
Inherit this log and proceed. ~4h to the (unconfirmed) 18:00 deadline.

### State found on arrival: a second, fast-moving worker has nearly finished the submission
This session found the tree far ahead of where Session 3 left it, with no corresponding log
entries yet — `log0.md`'s last commit (`c23b4e2`, 13:31) predates all of it. In the ~40
minutes between that commit and this session starting, uncommitted work produced: an
expression/ML scoring stack (`pipeline/{binding_data,features,train_binding_model,
train_affinity_model,score_candidates_ml,leakage_audit,esm_embed,facs_scorer}.py`), a full
de novo VHH campaign (`pipeline/denovo/`, `runs/denovo-vhh/`, 21 JapanFold jobs, $2.61
spent), AF2-IG scoring of the production batch (`data/af2ig_scores.csv`, still being
written as this session started), and — appearing literally while this session was
reading the tree (mtime 14:09, after Session 4 began) — **`runs/final-shortlist/
submission_table_18.csv`, an 18-design optimisation-track submission table with wild-type
control, ML rank, AF2-IG ipTM/iPAE, and liability flags per row.** All of the above is
**still untracked/uncommitted** as of this writing. **Do not touch `pipeline/`, `data/`, or
`runs/` paths outside this log while that worker may still be writing — read-only is safe,
editing is not.**

### What that work resolved, cross-checked against Sessions 1–3's open questions

**1. The de novo VHH route is dead, and it is a real negative result, not a dead end from
lack of trying.** `runs/denovo-vhh/README.md` (213 lines, uncommitted, explicitly flagged
"nothing here is committed" — **rescue this before it can be lost**) ran BoltzGen
`nanobody-anything` for real: 21 JapanFold jobs, 840 designs requested, 160 analysed, 88
unique sequences, $2.61 of $1,100 spent. Verdict: **BoltzGen's epitope specification does
not steer the design.** Under Boltz-2 unguided co-folding — the only method tested that
passes its own positive control (reproduces the real trastuzumab–HER2 epitope 19/19
residues) — all 6 top candidates dock to a different patch of HER2 domain IV (~495–563),
contacting 0–3 of 19 epitope residues and 0/4 hotspot atoms, while scoring a confident
0.76–0.92 ipTM there. BoltzGen's own internal ipTM is actively misleading: corr(target
distortion, reported confidence) = **+0.53**, i.e. worse geometry scores higher. AF2-IG
could not discriminate on this target either (crystal interface scores indistinguishable
from median design — confounded by VH-alone not being a binder without VL).
- This reverses Session 1's and Session 3 action 5's recommendation to pursue BoltzGen
  `nanobody-anything` for the de novo track. **Do not spend remaining JapanFold budget on
  this route.** If retried later, the README's own fix suggestion is geometric epitope
  constraint (RFdiffusion3/PXDesign hotspot restraints, not sequence annotation) — a
  bigger job than remains before 18:00.
- **Novelty is not the blocker.** 44/85 unique designs classify `de_novo` under abguard
  (`identity_to_trastuzumab` < 0.70), 75/85 pass full validation. The README also ran a
  check this log should have run itself: of 120 real, published NbThermo nanobodies, **0
  pass a whole-chain "<70% identity to any known antibody" reading** (min identity to
  nearest neighbour 0.730, median 0.983) — so that reading of the 70% rule cannot be what
  is enforced; abguard's identity-to-trastuzumab-only definition is the only operational
  one, and it is satisfiable. Worth confirming with organisers anyway, but it is not
  currently blocking anything.
- **Thermostability is not the blocker either — and this is the direct payoff of Session
  1's resource list.** NanoMelt was actually installed and validated (not just scoped):
  MAE 1.59 °C, Pearson r 0.916 against 40 NbThermo nanobodies (caveat: likely inside
  NanoMelt's own training set, so this validates the install, not generalisation — quote
  the paper's held-out numbers in the pitch, not this number). **The CC BY-NC-SA 4.0
  licence must still be disclosed**, per Session 1.
- **Confirms Session 1's disulfide lever, with a number attached.** All 160 raw BoltzGen
  designs came back cysteine-free (it deletes the canonical VHH C22–C92 framework
  disulfide). Restoring it costs 2 framework substitutions, no CDR contact, and gains
  **median +7.0 °C** (range +3.5 to +12.1, 85/85 improved): 60.6 → 69.7 °C median, max
  76.1 °C. This is a clean, citable, reproducible-on-our-stack result even though the
  overall route is dead — worth a slide regardless ("we found and fixed a thermostability
  bug BoltzGen introduces by default").

**2. The epitope hotspot question (log0 next-action 3 / open question 1, and this log's
Session 3 action 3) is resolved — and resolved more directly than the 7mn8 superposition
this session set out to do.** The same README measured actual 5 Å contacts in 1N8Z
directly: HER2 chain C residues `557 558 560 561 569 570 571 572 573 579`, chain D
residues `591 593 596 598 600 602 603 604 605`. The pipeline's four hotspot atoms resolve
to **Pro572 CD, Phe573 CD1, Ile591 CG2, Lys593 NZ** — all four fall inside that measured
contact set. (The task brief's residue *names* for two of them, "Leu573"/"Val591", are
wrong; the atom names and positions are right.) **The existing epitope definition is
sound; no further structural check is needed.** This session had planned to cross-check
via 7mn8 (confirmed reachable: `files.rcsb.org/download/7MN8.pdb` returns 200, but no
Biopython on this box, only numpy/scipy/sklearn in `.venv-ml/`) but that is now redundant
work — skipped.

**3. The optimisation track has a finished submission artefact.** `runs/final-shortlist/
submission_table_18.csv` (18 designs + wild-type control) and
`optimisation_shortlist_18.csv` (ranking detail: `p_beat_wt`, `pred_kd_nM` with 90% interval,
`removes_dg`, `pro_for_gly`, `net_charge`, `af2ig_iptm`/`af2ig_ipae`, `ml_rank_of_51417`) —
consistent with log0's accumulated findings (distance-3 preference, DG-hotspot removal,
winner-tolerant positions) and this log's AF2-IG recommendation, now actually computed. No
de novo designs are in it, consistent with finding 1 above. Not yet inspected row-by-row
for the charge/expression tension log0 flagged in finding 0c (net charge +4.2 to +6.2) —
worth one pass before submission if there is time, but do not block on it.

### What remains, revised given ~4h to deadline and the above
Only items that are still open and that this session did not already see done:
1. **Rescue `runs/denovo-vhh/README.md` and its data into a commit before it can be lost.**
   It explicitly says it is uncommitted; it is the only record of a genuinely valuable
   negative result (de novo route dead) and a positive one (disulfide fix, NanoMelt
   validated). Coordinate with whoever is still writing to `runs/` before committing out
   from under them.
2. **Confirm the 18:00 deadline** — still unconfirmed by anyone, four sessions in.
3. **Disclose the NanoMelt CC BY-NC-SA 4.0 licence in the pitch** if any de novo material
   (even just the disulfide finding) is used as a slide.
4. The framework-mutation rules question (this log's Session 3, action 1) remains
   unresolved and now moot if the de novo route is abandoned and the submission is
   optimisation-track-only — the `S85N+R87T` framework-mutation finding only pays off if
   framework edits are confirmed legal, and nobody has asked.
5. Given the submission table already exists, the remaining highest-value use of time is
   likely **pitch construction** (Judgment = 40%), not more scoring — echoing log0's own
   priority #2.

### Dead ends / non-findings this session
- Planned to superpose 7mn8 on 1N8Z to validate the epitope hotspots; turned out already
  answered directly and more rigorously (measured 5 Å contacts in 1N8Z itself) by the
  concurrent de novo VHH work. Did not need to run it.
- No Biopython on this box (`.venv-ml/` has numpy/scipy/sklearn only) — if the 7mn8
  superposition is ever wanted, budget time to install or hand-roll Kabsch on CA atoms.
