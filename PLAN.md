# Who does what — final run to 18:00

Deadline **18:00 Oslo time (CET) — confirmed.** Revised 13:35 CET with **4.4
hours left.** All times below are CET.

Discussion: [STRATEGY.md](STRATEGY.md) · Commands: [EXECUTION.md](EXECUTION.md) ·
Checklist: [TASKS.md](TASKS.md)

---

## Timeline, working backwards

| By (CET) | What must be true |
| --- | --- |
| **16:30** | **All analysis stops.** Whatever is not finished does not go in the deck. |
| 16:30 | PhDs have returned the design pick (see findings/phd-brief.md) |
| **17:00** | Deck assembly starts — **from a pre-drafted skeleton**, not from scratch |
| 17:40 | One person owns the 5-minute version; read it aloud once |
| **17:50** | Submitted |

**One hour for the deck is tight.** Mitigation: the deck outline and every number
are frozen in this file before 16:30, so 17:00 is assembly and rehearsal only. If
anything is still being computed at 16:30 it does not go in — no exceptions, or
we will be editing slides at 17:55.

The submission axis is already safe: 14/14 valid against the 525,063-sequence
gate. Everything from here is design *selection* and deck, not generation.

---

## Master (one person + this session)

The master does not do analysis. The master owns correctness and the story.

1. **Validate every worker result before it enters the deck.** Re-run the script,
   check the number against raw data. This has already caught one false claim
   (the phantom "454 non-binders") that had spread into five files.
2. **Own the submission table.** One file, `team-validate` clean against the
   525,063-sequence reference, wild-type control still failing.
3. **Own the deck narrative.** One author. The findings come from workers; the
   argument does not get written by committee.
4. **Keep a single source of truth.** When a worker corrects something, fix it
   everywhere — the 454 error survived because it lived in five places.
5. **Decide what is *not* claimed.** See "What we must not claim" below.

---

## Three kinds of worker — and which is which

This was underspecified in the first version of this plan. "Worker" means three
different things here, and the distinction decides what can actually be run in
parallel:

| Kind | Good for | Limits |
| --- | --- | --- |
| **Spawned subagent** (master spawns it) | bounded, well-specified jobs; **anything needing a mind that has not seen our conclusions** | no MCP tools (so no AMASS); cannot make judgement calls; returns only a report, so its reasoning is not inspectable |
| **Peer Claude session** (own user, own permissions) | long-running independent lanes; work needing its own approvals | will not act on a peer's authorisation for outward-facing calls — correctly |
| **Human** (us, the PhDs) | judgement, biology, anything irreversible | scarce; do not spend on anything a script can do |

**A spawned subagent is the only clean way to run W3.** It inherits none of this
session's context, so it cannot confirm our hypotheses back to us — which is the
entire point of an adversarial re-derivation. Spawn it with the *data and the
questions*, never the answers, and explicitly forbid it from reading
STRATEGY.md / EXECUTION.md / PLAN.md / TASKS.md / findings/.

**Verify what subagents return.** A subagent's report arrives in the same
confident register whether it is right or wrong, and its working is not visible.
Re-run its script before its numbers enter the deck — the same gate that caught
the phantom "454 non-binders".

## Workers — parallel, independent, each leaves a script

Rule for all: stdlib only, script in `pipeline/`, print raw numbers **and row
counts**, and state what you could not verify. A result without a reproducible
script does not go in the deck.

### W1 — FACS binding model *(in progress, temp-user-dc)*
Measured-data ranking signal, the only one we have. **Done:** test AUC 0.9842,
macro-AUC over test clusters 0.9306, zero train/test CDRH3 overlap, wild-type
control P(high)=1.000. **Holding the per-design ranking** because the mid-table
moved between training runs — correct call. Deliver either a seed-stable ranking
or an explicit "only the top 3 are stable, the rest unranked".

### W2 — StaB-ddG per-position scan of B99–108 — **DONE, and it moved the answer**

Alanine scan of the whole window in the wild-type background, plus W99→Ser and
W99→Tyr. Ran in 118 s, not the 12 min advertised. Binding ddG, kcal/mol:

| mutation | ddG | ± MC error |
| --- | ---: | ---: |
| **Y105A** | **0.866** | 0.114 |
| F104A | 0.462 | 0.104 |
| W99Y | 0.312 | 0.061 |
| **W99A** | **0.269** | 0.140 |
| G101A | 0.216 | 0.035 |
| **W99S** | **0.207** | 0.148 |
| D108A | 0.153 | 0.104 |
| G103A / G100A | 0.042 / 0.028 | |
| M107A / D102A | −0.087 / −0.140 | |

**Three things follow.**

1. **Tyr105, not Trp99, is the load-bearing residue** — 3× the effect and well
   outside W99A's error. W99→Ser, the substitution 12 of our 14 designs make,
   is 0.207 ± 0.148: indistinguishable from zero.
2. **All 14 designs preserve Y105. Only 1 of 14 preserves F104.** So the residue
   that matters most is already protected everywhere, and the real exposure is a
   position nobody had looked at.
3. **Every single-position effect is ≤0.87 kcal/mol — at most ~4.3× in Kd, i.e.
   inside the 6.2× assay noise floor.** No single-residue choice in this window
   can produce a *measurable* affinity change. That is an independent
   confirmation, from structure rather than from the SPR table, that
   CDRH3-tip redesign is the wrong lever for affinity.

Second job (`designs-ddg`, the 14 as multi-site variants) submitted. The API
accepts multi-site variants so the model scores them jointly rather than us
summing — better than feared — but whether StaB-ddG's training covers
multi-mutants is **unverified**, so that output stays ordinal and flagged.

### W3 — Adversarial re-derivation (Lane C) *(RUNNING — spawned subagent)*
Re-derives the five headline findings from the raw tables with no sight of our
conclusions. Both Claude sessions are disqualified (we have both read them), so
this is a **freshly spawned subagent** — it inherits no context, which is the
whole point. Given the data and the questions only, and forbidden from reading
our analysis docs. Deliver: confirmed / refuted / cannot-be-determined, with its
own numbers and sample sizes. **Re-run its scripts before trusting it.**

### W4 — AMASS evidence pull
Its own prize pool, currently untouched. Needs a restart to load
(`cd ~/antibody && set -a && . ./.env && set +a && claude`), then the evidence
scout skill. Target the thermostability case: the NB46 precedent (competes with
trastuzumab for HER2, holds affinity to 90 °C) and VHH thermostability in
general. Deliver cited records — this doubles as deck content and as the pool-4
submission.

### W5 — Liability report *(RUNNING — spawned subagent)*
Never run, and part of what `team-validate` enforces. `abguard liability-report`
plus `--all-liabilities` on the 14 shortlisted designs, and the citation for the
LAP liability reference. Bounded and CPU-only, so a good subagent job.
`abguard dashboard` (marimo, needs an SSH port-forward) is optional extra.

### Not doing
De novo track (needs different weights, a fresh batch and new novelty checks on a
GPU that dies tonight). Repeat folding (proven worthless: every replica
aggregation lands at AUC 0.510–0.543). Chasing an affinity record (1.3% odds).

---

## The two PhDs — decide, don't compute

Everything here is a judgement call that no model in our stack can make. Each
needs a yes/no we can act on, not a discussion.

### Q0 — Is W99 load-bearing? **Now three-way evidence, and it has shifted.**

**ANSWERED — W99 is real but weak. Neither of the two positions we held was
right.** Independent stratified analysis (findings/lane-c-independent.md):

- Raw W99 effect **8.11×** (p=2e-08) → **2.43× adjusted for mutation count**
  (p=0.015). It **survives** stratification (van Elteren p=5.4e-04), so it is not
  purely a confound — but half the raw advantage was mutation load, only 2 of 7
  strata are individually significant, and **two strata point the other way**.
- **2.43× sits barely above the ~2.5× best-case assay resolution.** So W99 is a
  weak-to-moderate prior, not a hard rule. `FGSRGFYAYA` (0.471 nM, 5 mutations)
  proves it is not required.
- **StaB-ddG independently agrees**: W99→Ser costs 0.207 ± 0.148 kcal/mol
  (≈1.4×). Two unrelated methods, same conclusion: small, maybe real, not
  decisive.

**Practical answer: do not let W99 drive the pick.** It is worth ~2.4×, below
what the assay can show, and keeping it costs a `TrpOx` developability liability.
Select on the convergent evidence instead (see findings/shortlist-record.md).

The earlier two-line framing, kept for the record:

- **StaB-ddG (W2):** W99→Ser costs 0.207 ± 0.148 kcal/mol — zero within error.
  **Tyr105 is the real hotspot at 0.866**, and all 14 designs already keep it.
- **The FACS model (W1):** both stable-high designs (`_19`, `_85`) *lost* W99,
  and `_60` — a W99-retainer — scores 0.133, near the bottom.
- **The confound is real:** W99-keeping SPR variants carry median **4**
  substitutions versus **6** for W99-losing ones, so much of the 18.8-vs-147 nM
  gap is mutation load, not W99.

**But it does not fully survive stratification, and that is the live question.**
Within each substitution count the W99-keepers are still usually tighter, and at
**6 substitutions — our target band — it is 13.1 nM versus 147.0 nM, an 11×
difference** (n=10 vs 19). So: is that residual an artefact of small n, or a real
effect that StaB-ddG underestimates because it scores a static structure?

**What we need from you:** a yes/no on whether to treat W99 as decisive. If no
(which is where the evidence now points), the submission leads with `_19` and
`_85` on the measured-data tiering, and `_60` drops — note it also destroys F104,
the second hotspot, which may be exactly why the FACS model scores it 0.133.

The original framing, for context:
The conserved tryptophan at the base of CDRH3 is kept by **8%** of our
novelty-accepted designs but **79%** of the rejected ones, and SPR median
affinity is **18.8 nM with it versus 147 nM without**. The novelty gate
systematically strips out the residue the measured data says you need.

Only **2 of our 14** candidates keep it: `her2_cdrh3_60` (`WSYAYDYAMD`, 5 subs)
and `her2_cdrh3_94` (`WTAAVDYAHD`, 6 subs). They rank 11th and 13th of 15 on
`min_ipae`.

**The decision:** is that SPR correlation mechanistic, or is it confounded by how
the Chinery library was built? If mechanistic → `_60` and `_94` are the
submission and the other twelve are padding. If confounded → rank on `min_ipae`
and submit `_71`/`_88`. W2's ddG scan will give an independent line, but the call
is yours.

### Q1 — Are the four hotspot atoms the right epitope definition?
**Read the corrected identities below — two of four were mislabelled in every
document until 13:00 today**, so any earlier discussion of this question was
based on wrong residue names.

Verified directly against `inputs/12-her2/her2_fv_1n8z.cif`:

| selector in `her2_cdrh3.yaml` | residue actually there | previously called |
| --- | --- | --- |
| `C/*/572/CD` | **PRO 572** | Pro572 ✓ |
| `C/*/573/CD1` | **PHE 573** | "Leu573" ✗ |
| `D/*/591/CG2` | **ILE 591** | "Val591" ✗ |
| `D/*/593/NZ` | **LYS 593** | Lys593 ✓ |

The selections still resolved to real atoms — PHE has a CD1, ILE has a CG2 — so
the model was aimed at genuine coordinates and the geometry is sound. Only the
*description* was wrong. But that description is what this question asks you to
validate, and "is Leu573 the right hotspot?" is unanswerable when the residue is
Phe573.

**That set is our entire operational definition of "targets the trastuzumab
epitope"**, and all 100 designs inherit it. If it is wrong, everything downstream
is aimed at the wrong patch. Blocks any future batch.

### Q2 — Did the loop stay on the epitope, or drift?
A folding model returns a confident complex with the paratope on the wrong
surface patch, and ipTM cannot tell — we now know Boltz-2's ipTM admits 94.3% of
*shuffled* pairs at the 0.70 threshold, so confidence is no guard at all here.
**Structures are committed**: `runs/team15-claude-n100-seed1/shortlist-structures/`,
28 files, `*_designed.cif` and `*_refolded.cif` for each candidate. Needs eyes.
*(A browser viewer can be built if that helps — ask.)*

### Q3 — Is the motif collapse real or an artifact?
The accepted designs converged hard: **11 of 15 begin `ST`**, many are
`ST..V[DG]YAMD`, and two were byte-identical. Is `ST`-initiated CDRH3 plausible
for HER2 binding, or is this the model exploiting the novelty gate by abandoning
the W99 anchor? Related: is rigidifying the loop with prolines (`GDG`→`PPY`)
clever or a liability? If the motif is wrong, **all 14 candidates fail together**
— our real diversity is much lower than the count suggests.

### Q4 — Sanity-check the thermostability argument
We cannot score it (see below), so we argue it: VHH scaffolds are intrinsically
more thermostable and express far better in Adaptyv's cell-free system, where
expression is an unscored gate that zeroes you out entirely. Is that defensible
as stated?

---

## The presentation

Judgment is **40%**. The findings *are* the pitch — lead with what we learned
about the problem, not with what we built.

**1. The trap, in one figure.** Two curves running opposite ways:

| subs from WT | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| novelty acceptance (our 100) | **0%** | **0%** | **0%** | 33% | 82% | 74% | 93% | 100% |
| measured high-binder rate | 98.2% | 90.9% | 72.7% | 58.5% | 42.2% | 31.7% | 28.5% | 24.5% |

Caveat to state if asked: the SPR arm of this is **survivorship-censored** — only
3.8–5.6% of 9–10-mutation variants yield a KD at all, versus 100% at 0–1, so the
real affinity cliff is *steeper* than any table of medians shows. The trustworthy
span is the first four mutations: **0.095 → 34 nM, ~360×**.

The novelty gate forces designs into exactly the regime where measured binding is
least likely. **Zero of our 13 conservative designs survived it.** Two designs
were 0 mutations from an already-measured variant — the model re-invented known
sequences.

**2. Why the usual metrics cannot save you.** Rank the folding models by ROC AUC
against ~8,190 non-pairs: **AF3 0.754, Boltz-2 0.645, Chai-1 0.568** (lower CI
bound 0.506 — barely a coin flip). Boltz-2 is what we fold with.

Then the sharper point: **a fixed ipTM threshold is meaningless across models.**
At ipTM ≥ 0.70, Boltz-2 admits **94.3%** of known non-pairs while Chai-1 admits
**0.2%** — yet Chai-1 is the *worse* model. The threshold measures the model's
scale, not its quality: Boltz-2's non-pair median is already 0.884, Chai-1's
positives top out at 0.814. Our candidates sit at 0.82–0.90, inside the range of
random pairs. Folding repeatedly does not help either (AUC 0.510–0.543).

Do **not** present the 94.3%/0.2% comparison as a model ranking — it inverts the
real one. Present AUC, then use the threshold figures to show why thresholds
travel badly. Validation: the AF3 diagonal reproduces three cognate ipTM values
printed in the supplementary PDF (8W90 0.12, 8AV2 0.15, 9ETJ 0.46) to the digit.

**3. What we did instead.** Novelty-gate first, then fold only survivors — the
default ranking wasted 30% of our GPU on designs the gate rejected. Then rank on
measured data: a FACS-trained model, test AUC 0.9842, macro-AUC 0.9306 across
clusters, zero train/test overlap. Our designs are Hamming 3–4 from any measured
sequence yet sit in a regime covering 95% of training rows — novel, but densely
covered. Base rate there is ~28–32%, so a high score is real lift.

**4. The honest part — and this is the Judgment 40%.**
- The assay's own floor is **~7×** (wild-type measured 8 times, 0.0285–0.178 nM;
  log10 SD 0.307), so we claim no affinity win smaller than that. The replicates
  are also **bimodal** — split 3.44× by printed decimal precision, suggesting two
  merged runs — so even 6.2× was generous.
- **Nothing in the measured data beats trastuzumab by a resolvable margin.** Two
  of 150 variants at ≥3 mutations are nominally tighter; the best beats
  wild-type's own best replicate by 1.22× against a ~7× floor. Zero real wins.
- **Nothing in either stack predicts thermostability.** One of the three scoring
  axes is unaddressed. We say so.
- W99: the gate strips the anchor residue; we flag it as unresolved rather than
  pretending the ranking settles it.
- We corrected **five** of our own inherited claims mid-project: a phantom "454
  measured non-binders" that had spread through five files; abguard's novelty
  gate being on by default (it is not); `stab-ddg` being a thermostability model
  (it predicts binding ddG); the IMGT CDRH3 matching our design window (it does
  not); and **two of the four epitope hotspot residues being misnamed** — "Leu573"
  is Phe573 and "Val591" is Ile591. Each was found by checking a claim against
  the raw data rather than inheriting it. That is the method, and it is worth a
  slide of its own.

**5. What we would do next.** Novelty-aware design (constrain generation against
the 525,063 reference *during* sampling, not after), VHH format for
thermostability and expression, and an assay-aware target: match wild-type within
noise while clearing novelty, rather than chase a record.

### What we must not claim
- No affinity win below 6.2×.
- No thermostability number.
- No ipTM threshold without naming the model.
- No per-design ranking that is not seed-stable.
- Nothing from the prize-pool map as fact until the PDF is readable.
