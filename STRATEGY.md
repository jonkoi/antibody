# What to do to win — for discussion

14 Nov 2026. For the whole team, including the two bio PhDs.

This is the **discussion** document: the argument, the decisions and who does
what. The commands, paths and scripts are in **[EXECUTION.md](EXECUTION.md)**.

Also: [TASKS.md](TASKS.md) (checklist) · [HANDOFF.md](HANDOFF.md) (state)

> **Deadline 18:00 today — confirmed.** If that is Oslo local time (CET) and this
> machine is on UTC, check the hour you actually have left.
>
> **The organisers' PDF slide deck**
> (`2026-11_Oslo_hackathon_Apheris_challenge.pdf`) is still unreadable — 302,382
> UTF-8 replacement bytes, 0 of 159 streams decompress. HANDOFF.md is fine; this
> is only about that PDF. So the prize pools below come from our own notes rather
> than the organisers' text, and are the one thing here not independently
> verified.

---

## The decision

Three of the five prize pools are decided **today by the pitch**, Judgment is 40%
of it, and pitch placement multiplies our design slots from 4 to 16. The two
experimental pools are decided weeks later — and we have measured that they are a
**1.3% lottery**.

**So the pitch is the primary deliverable, our findings are what make it strong,
and thermostability is the axis we argue rather than chase.**

| Pool | Decided | Realistic? |
| --- | --- | --- |
| Best pitch ($500) | **today** | **high** |
| AMASS best use | **today** | **high** |
| Pitch placement | **today** | **high, and multiplies the two below** |
| Optimisation (€800/500/200) | weeks later, by SPR | low |
| De novo (€800/500/200) | weeks later, by SPR | very low — we have nothing built |

Designs are scored on **affinity**, **epitope preservation** and
**thermostability**. S-tier improves all three at once.

---

## What we know that others probably don't

All measured from the organisers' own data. This is the product, as much as the
designs are.

**1. The assay cannot see small wins.** Wild-type trastuzumab appears eight times
in the SPR table, measured between 0.0285 and 0.178 nM — a **6.2-fold spread on
an identical sequence**. Any affinity improvement smaller than that is
indistinguishable from noise.

**2. Affinity collapses as you mutate.** Median binding gets ~1000× weaker by
seven substitutions, and the decline is monotonic. Conservative designs bind;
adventurous ones don't.

**3. "Novel *and* better than trastuzumab" is a 1.3% event.** Of 153 measured
variants that clear the novelty floor, two beat wild-type's median — and the
better of those sits *inside* wild-type's own replicate spread, so it may not be
an improvement at all.

**4. The novelty filter and affinity are in direct conflict.** This is the
important one. Our design `her2_cdrh3_2` is three substitutions from
trastuzumab — which we thought cleared the gate — and it was **rejected**,
because it is only two mutations from a variant in the FACS table. That table's
524,346 variants densely cover the near-wild-type region, which is exactly where
binding survives. **The rules push us away from the only place designs work.**

**5. The trap is now measured directly, on our own 100 designs.** We ran the real
novelty gate over the whole batch. **71 of 100 passed — but acceptance depends
entirely on how adventurous the design is:**

| Substitutions from wild-type | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **Novelty acceptance** | **0%** | **0%** | **0%** | 33% | 82% | 74% | 93% | 100% |
| Median measured affinity (nM) | 20.5 | 15.0 | 60.7 | 92.7 | 93.0 | 146 | 249 | 191 |

**Not one of our 13 conservative designs (≤4 substitutions) survived. Zero of
thirteen.** Acceptance climbs monotonically as designs get more adventurous;
measured affinity collapses monotonically over the same range. **The two curves
run in exactly opposite directions**, which is the trap stated as a number.

**But the accepted designs are not extrapolations** — an earlier version of this
document implied the viable region was simply empty, which overstated it. The
measured high-binder rate by substitution count (from the FACS training data) is
98.2%, 90.9%, 72.7%, 58.5%, 42.2%, 31.7%, 28.5%, 24.5% at 2–9 substitutions, and
**our accepted designs sit at 5–9 substitutions, a regime holding 95.2% of all
training rows.** They are Hamming 3–4 from any individual measured sequence yet
squarely inside the measured distribution: novel, but densely covered. The base
rate there is ~28–32%, so a high model score is genuine lift, not noise.

Two of the rejected designs were **zero mutations** from a variant already in the
FACS table — the model re-invented sequences that have already been measured.

**6. The pipeline's own ranking makes this worse.** `refold-prepare` ranks on
interface H-bonds and hotspot obedience. In our batch it selected almost
exclusively high-substitution designs, and **30% of what it picked was later
rejected by the novelty gate** — GPU time spent on designs that cannot be
submitted, in the substitution range where binding does not survive anyway.

**7. Boltz-2's confidence score barely distinguishes real pairs from random
ones.** Measured on ~8,100 shuffled antibody–antigen pairs, ROC AUC by model:
**AF3 0.754, Boltz-2 0.645, Chai-1 0.568.** Boltz-2 is the model we fold with.
The false-positive rate of an ipTM cut-off is damning:

| ipTM cut-off | Boltz-2 false positives | AF3 false positives |
| --- | ---: | ---: |
| ≥ 0.70 | **94.3%** | 3.0% |
| ≥ 0.82 | 76.6% | 1.1% |
| ≥ 0.90 | 39.8% | 0.4% |

Our candidates sit at ipTM 0.82–0.90 — **inside the interquartile range of
shuffled, non-binding pairs.** Getting Boltz-2 to a 5% false-positive rate needs
ipTM ≥ 0.951, which discards 82% of genuine pairs. This is the hard evidence that
ranking on `min_ipae` rather than ipTM was right, and that **any confidence
threshold quoted without naming the folding model is meaningless.**

Two caveats to keep: the negatives are *shuffled* pairs, presumed non-binding by
the study's design, not SPR-measured ones. And our refold is handed the HER2
template with framework and CDRH1/H2 frozen, so our ipTM is not on this scale —
what transfers is the ordering of models and the shape of the trade-off, not a
number to paste into a filter.

**Correction:** earlier versions of this document, README.md and HANDOFF.md all
claimed "454 SPR-measured non-binders, 21 scoring ipTM ≥ 0.7". **That set does
not exist.** `spr-supplementary/` holds 29 and 12 rows, "454" appears nowhere in
the data, and its README explicitly says "N/A denotes no reportable affinity fit,
NOT confirmed non-binding" — those rows must not be used as negatives. The claim
was inherited and repeated without being checked.

**8. Folding a design repeatedly does not help rank it.** Across 14,721 pairs ×
50 replicas, every aggregation of ipTM (max, mean, std, mean−std) lands at AUC
0.510–0.543 with confidence intervals spanning 0.5. Best-of-50 is
indistinguishable from a single fold. **Do not spend GPU on repeat folds hoping
averaged confidence will discriminate.** (Chai-1 only; the Boltz-2 per-replica
archive is 1.7 TB and not in the bucket.)

---

---

## Where that leaves us: aim at six substitutions

The two curves above cross in a narrow band, and that band is the whole strategy.

**Six substitutions is the sweet spot.** Novelty acceptance is 82%, and the best
6-substitution variant ever measured binds at **0.185 nM** — against wild-type's
median of 0.135 nM. That difference is **well inside the assay's own 6.2× noise
floor**, so a good 6-substitution design is, as far as this assay can tell,
*as tight as trastuzumab* while clearly clearing the novelty gate.

Below five substitutions nothing survives novelty. Above seven, measured affinity
is gone. We have **15 novelty-accepted designs at 5–6 substitutions**, and they
are the submission pool. Only two had been folded, because the pipeline's ranking
ignored the rest — the other thirteen are folding now.

That is the submission thesis: *not* "we beat trastuzumab", but **"we found the
only region where novel and functional overlap, and we have designs in it."**

### But the gate also strips out the conserved tryptophan — read this before picking

W99, the tryptophan anchoring the base of CDRH3, is where the three signals
collide:

| | keeps W99 | loses W99 |
| --- | ---: | ---: |
| Our designs **accepted** by novelty (n=71) | **8%** | 92% |
| Our designs **rejected** by novelty (n=29) | **79%** | 21% |
| Median measured affinity, SPR | **18.8 nM** | 147 nM |

**The novelty gate selects against the one residue the measured data says you
need.** Keeping W99 costs you novelty; losing it costs you ~8× in median
affinity. This is the same trap as before, now at single-residue resolution, and
it is the sharpest thing we have found.

**Consequence: our confidence ranking and our biology disagree, and we should
trust the biology.** Only **2 of the 15** candidates keep W99 —
`her2_cdrh3_60` (`WSYAYDYAMD`, 5 subs) and `her2_cdrh3_94` (`WTAAVDYAHD`, 6
subs). They rank 11th and 13th of 15 on `min_ipae`, yet they are the only ones
that satisfy novelty *and* retain the anchor residue. The shortlist is ordered
W99-first for exactly this reason, not by confidence score.

**For the PhDs — this is question zero.** Is W99 genuinely load-bearing here, or
is the SPR correlation confounded by the library's design? If it is load-bearing,
`_60` and `_94` are our submission and the other twelve are padding. If it is not,
we should rank on `min_ipae` and submit `_71` and `_88`. Nothing we can compute
settles it.

**One honest limitation:** the accepted designs collapsed onto a narrow motif —
11 of 15 begin `ST`, and two were byte-identical (deduplicated, so 14 remain). If
that motif is wrong, our candidates fail together. Real diversity is lower than
the count suggests.

---

## What we have, and what's missing

**Have:** a validated design pipeline · a 100-design batch running now · the real
novelty reference, all 525,063 sequences · a working submission validator · all
the data · the AMASS evidence engine · **two bio PhDs** · one H100 **that dies
tonight**.

**Missing:** authoritative rules (the deck is still corrupt) · a confirmed
deadline · any de novo design · epitope validation · a calibrated quality
threshold · **and any way at all to predict thermostability — see below.**

### Correction: we cannot score thermostability, and `stab-ddg` is not what we thought

Nebius credentials are now installed and working (34 of 59 models). But
**`stab-ddg` is StaB-ddG, which predicts *binding* ddG for point mutations across
two partner chains — not melting temperature or folding stability.** Its only
operation is `predict-binding-ddg`.

HANDOFF, and the earlier version of this document, treated it as our
thermostability scorer. That was wrong. **No model in either stack predicts
thermostability.** Two consequences:

- **Thermostability becomes an axis we *argue*, not one we score.** The case rests
  on format (VHH scaffolds are intrinsically more thermostable and express far
  better in Adaptyv's cell-free system, where expression is an unscored gate that
  zeroes you out) and on the NB46 precedent. In the pitch we should say plainly
  that we could not measure it — admitting the gap is stronger than implying a
  number we do not have.
- **`stab-ddg` is promoted to something more useful: ranking our shortlist.** It
  scores mutations across the antibody–HER2 interface, which is the affinity axis,
  and it is **independent of ipTM and `min_ipae`** — the metrics we have proven
  are saturated. It is the only discriminating signal available to us.

**And it can settle the W99 question — this is its right use.** W99→X in the
wild-type background *is* a point mutation, exactly what the model was built for.
A per-position ddG scan across B99–108 in the wild-type background would turn our
W99 evidence from correlational (8% vs 79% retention; SPR median 18.8 vs 147 nM)
into an independent mechanistic line. **Prioritise this**; ~12 min queue, so send
it early.

**But using it to rank the 14 is off-label, and we must say so.** StaB-ddG scores
**point** mutations. Our designs carry 5–9 simultaneous substitutions in a
10-residue loop, and summing per-mutation ddG assumes additivity — least safe in
CDRH3, where the substituted residues contact each other and reshape the backbone
together. Treat any such ranking as ordinal at best, check whether StaB-ddG's own
validation covers multi-mutants, and state the limitation rather than hide it.

**Say plainly in the pitch that one of the three scoring axes is unaddressed.**
Nothing in either stack predicts thermostability, and letting StaB-ddG stand in
for it would be misrepresentation.

---

## Who does what

**Organisers / team chat:** deadline **confirmed 18:00**. Nebius credentials
installed. Still outstanding: **ask whether the submission wants the structures**
(we hold 111MB that dies with the VM — see TASKS.md), and re-share the PDF deck
through a channel that does not re-encode text, or just confirm the prize pools
and scoring axes in a message.

**The two PhDs — where a model cannot substitute.** Full write-up with sequences
and figures: <https://claude.ai/artifact/K1Kd6wCkfJ2v9yUYExp6M1>

1. **Are the four hotspot atoms the right epitope definition?** We aimed the model
   at HER2 Lys593 NZ, Phe573 CD1, Ile591 CG2 and Pro572 CD, taken from 1N8Z. That
   set is our *entire* operational definition of "targets the trastuzumab
   epitope", and every design inherits it. **This blocks any larger batch** — the
   most valuable question in the project.
2. **Did the designed loop stay on the epitope, or drift off it?** A folding model
   will happily return a confident complex with the paratope on the wrong surface
   patch, and the confidence scores cannot tell. This needs eyes on the
   structures. *(We can build a browser viewer for them — just ask.)*
3. **Is the double proline substitution clever or a liability?** Rigidifying the
   loop could help both affinity and thermostability, or distort the paratope.
   Treat it as a transferable principle — that particular design is already
   novelty-rejected.
4. **Sanity-check the thermostability argument**, especially the VHH format case
   and the NB46 precedent (competes with trastuzumab for HER2, holds affinity to
   90 °C), which shows two of the three axes are satisfiable at once.

**Analysis.** The novelty-rate and feasible-region questions are **answered
above** — that was the pitch centrepiece and it came out sharper than hoped. Still
**also answered** — the confidence calibration is done (findings 7 and 8 above),
and it both vindicated ranking on `min_ipae` and killed the idea of repeat folds.
What remains is running `stab-ddg` to rank the 14 candidates on binding ddG and
to settle the W99 question. Details in [EXECUTION.md](EXECUTION.md).

---

## What we are deliberately not doing

- **Not chasing an affinity record.** We measured the odds: 1.3%.
- **Not attempting the de novo track.** It needs different weights, a fresh batch
  and new novelty checks, on a GPU that disappears tonight. We will argue the VHH
  case in the pitch instead of building one.
- **Not claiming an affinity win smaller than 6.2×.** Judges who know this data
  will check. "Within assay noise of wild-type, with thermostability improved" is
  both honest and a stronger story.

---

## The pitch

Spine, in order:

1. The novelty/affinity trap, with the numbers.
2. Why the confidence metrics lie here — calibrated against the organisers' own
   measured non-binders.
3. The 6.2× noise floor, and therefore what an honest claim looks like.
4. Thermostability as the real differentiator, with the NB46 precedent.
5. What we would do next, and why.

Point 5 is the Judgment 40%. Knowing *why* a 0.96 confidence score is not a
result is exactly what that rewards.
