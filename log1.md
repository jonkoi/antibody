# Agent work log — log1

Running log of agent sessions on this repo, written so the next agent inherits the
experience rather than rediscovering it. Newest session at the bottom.

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
