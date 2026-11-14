# Pitch deck — skeleton with frozen numbers

**Assembly starts 17:00 CET. Submit by 17:50.** Every number here is verified and
final; 17:00 is layout and rehearsal, not writing. If a number is not in this
file by 16:30 it does not go on a slide.

Judging: **Judgment is 40%.** So the deck leads with what we learned about the
problem, not with what we built. Ten slides, ~5 minutes.

Sources: [findings/](findings/) · full evidence in [STRATEGY.md](STRATEGY.md)

---

## Slide 1 — Title

**We set out to beat trastuzumab. We found out why almost nobody can.**

Team 15 · Trastuzumab CDRH3 redesign vs HER2 domain IV · Oslo AI × Antibodies

*Speaker note: do not open with the pipeline. Open with the finding.*

---

## Slide 2 — The task, in one line

Redesign the tip of trastuzumab's CDRH3 — **10 residues, heavy chain 99–108** —
to bind HER2 domain IV better, while differing from **525,063 known CDRH3
sequences** by ≥3 mutations.

Everything else frozen: light chain, heavy framework, CDRH1/H2, HER2 backbone,
and the relative pose.

---

## Slide 3 — THE FINDING ★ (the one slide that matters)

**The novelty filter and binding affinity pull in opposite directions.**

| substitutions from wild-type | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **passed the novelty gate** (our 100 designs) | **0%** | **0%** | **0%** | 33% | 82% | 74% | 93% | 100% |
| **measured high-binder rate** (524,346 variants) | 98.2% | 90.9% | 72.7% | 58.5% | 42.2% | 31.7% | 28.5% | 24.5% |

Two monotonic curves, running opposite ways, from **two independent sources** —
our novelty gate and the organisers' measured library.

- **Zero of our 13 conservative designs survived.**
- **Two designs were rejected for being 0 mutations from an already-measured
  sequence** — the model re-invented known variants.
- The gate forces every submission into the regime where measured binding is
  least likely.

*Speaker note: this is the whole pitch. If you only land one slide, land this.*

---

## Slide 4 — And the target is already near-optimal

From the organisers' SPR data (n=187 with a fitted K_D):

- Wild-type trastuzumab: **0.095–0.135 nM**
- Affinity falls **~1.69× per mutation** (ρ=+0.483)
- Of 150 variants at ≥3 mutations, **two are nominally tighter. Zero beat
  wild-type by a margin the assay can resolve.**

**The assay's own floor is ~7×** — wild-type measured 8 times spans 6.2×
(0.0285–0.178 nM), and the replicates are bimodal, splitting 3.44× by printed
precision, which suggests two merged runs.

**So "beat trastuzumab on affinity" is not a realistic target. We re-aimed at
matching it within noise while clearing novelty.**

*Caveat if asked: the SPR table is survivorship-censored — only 3.8–5.6% of
9–10-mutation variants yield a K_D at all, versus 100% at 0–1. The real cliff is
steeper than the medians show. Trustworthy span: first four mutations,
0.095 → 34 nM, ~360×.*

---

## Slide 5 — Why the standard metrics cannot rescue you

**Rank the folding models by ROC AUC** against ~8,190 mismatched pairs:

| model | ROC AUC | ipTM ≥ 0.70 admits | median, non-pairs |
| --- | ---: | ---: | ---: |
| AF3 | **0.754** | 3.0% | 0.180 |
| **Boltz-2** (ours) | **0.645** | **94.3%** | 0.884 |
| Chai-1 | **0.568** | 0.2% | 0.191 |

**Chai-1 has the best false-positive rate and the worst AUC.** A fixed ipTM
threshold measures the model's *scale*, not its quality — Boltz-2's non-pair
median is already 0.884; Chai-1's true positives top out at 0.814.

**Our designs sit at ipTM 0.82–0.90 — inside the range of random pairings.**

Repeat folding doesn't help either: across 14,721 pairs × 50 replicas, every
aggregation of ipTM lands at **AUC 0.510–0.543**.

> **A confidence threshold quoted without naming the folding model is
> meaningless.**

*Do NOT present the 94.3% vs 0.2% as a model ranking — it inverts the real one.
Lead with AUC, then use the thresholds to make the point about transferability.*

---

## Slide 6 — What we did differently

**1. Novelty-gate first, then fold.** The default pipeline folds the top 20, then
discovers rejections — it wasted **30% of our GPU** on designs that cannot be
submitted, all in the substitution range where binding does not survive.

**2. Rank on measured data, not confidence.** A logistic model (position one-hot
+ pairwise interactions) trained on 524,346 measured variants, on the authors'
clonotype-aware splits:

| model | test AUC |
| --- | ---: |
| edit distance from WT alone | 0.753 |
| PWM, independent positions | 0.921 |
| **unary + pairwise interactions** | **0.984** |

Macro-AUC within test clusters **0.931** — the honest number, since 662 of 706
clusters are single-class. Zero exact train/test CDRH3 overlap. Wild-type scores
P(high) = 1.000 as a control.

**3. Check developability.** All 14 candidates carry **fewer liabilities (5–7)
than approved trastuzumab itself (8)**, and introduce no new high-severity
liability. *(LAP reference, PLOS Comput Biol 2024;20(3):e1011881.)*

---

## Slide 7 — The result

**14 designs, all `team-validate` clean against the 525,063-sequence gate.**
Class `novel`, 97.4–97.8% identity to trastuzumab, minimum distance 3–4
mutations.

> **FILL AT 16:30 — the PhDs' pick goes here.**
> Lead candidate: `_61` (`SLGSADYAFD`) — best predicted binding ΔΔG of the 14
> (0.42 kcal/mol, ~2× of wild-type, inside the assay noise floor), high tier on
> the measured-data model (0.840), and tied-cleanest on developability (5
> liabilities, none in CDRH3).

*If the PhDs choose differently, use their pick and say why — a biologist's call
over a model's is itself the Judgment point.*

---

## Slide 8 — Where we are honest ★ (this is the 40%)

**Five claims we inherited and disproved by checking them against raw data:**

| claim | reality |
| --- | --- |
| "454 SPR-measured non-binders, 21 at ipTM ≥ 0.7" | **that set does not exist** — the file holds 29 rows, and its README forbids using them as negatives |
| abguard checks novelty by default | it does **not** — without `--reference-db` prior art is never checked |
| `stab-ddg` scores thermostability | it scores **binding** ΔΔG for point mutations |
| the IMGT CDRH3 is our design window | it is `SR` + our 10 residues + `Y` |
| epitope hotspots are Leu573 / Val591 | they are **PHE 573 / ILE 591** — our own annotation was wrong; the organisers' spec, which drove the run, was correct |

**Three things we will not claim:**

1. **No thermostability number.** Verified across **all 59** Nebius models and
   the whole Apheris stack: `stab-ddg` is the only stability-adjacent one, and it
   scores *binding* ΔΔG. One of the three scoring axes is **entirely
   unaddressed**, and a sequence-likelihood proxy dressed up as a Tm would be
   worse than saying so.

   **And the VHH argument is about a scaffold we did not build.** Our 14
   submissions are paired Fv designs sharing **220 of 227 residues (96.9%)** with
   each other and with wild-type trastuzumab — the entire light chain, heavy
   framework, CDRH1 and CDRH2 are identical, and only 7 positions vary. So large
   thermostability differences *among* our candidates are unlikely, though we
   cannot assert invariance: CDR and CDR–CDR interactions can carry much of the
   stability *(PMID 24143255)*. Nanobody Tm predictors exist (TEMPRO, NanoMelt,
   NBsTem) but **all are nanobody-specific and none transfers to a
   CDRH3-redesigned IgG scaffold**, so quoting one as our number would be wrong.
   Thermostability is the axis we could not score — and the scaffold change is
   what we would do about it, not something we did.
2. **No affinity win below ~7×**, the assay's own floor.
3. **No per-design ranking we cannot reproduce.** Our measured-data model has
   test AUC 0.984 but its per-design order moves between training seeds
   (Spearman 0.88–0.91), so we report **tiers, not a ranking**, and one design
   (`_77`) gets **no number at all** — it moved 0.229 → 0.790 on nothing but a
   reshuffled training order.

**And our two best signals disagree.** Spearman(ΔΔG, FACS) = **+0.522**: the
designs the measured-data model likes best are the ones the energetics model says
bind worst. Likely because **only 37 of 733 SPR sequences appear in the FACS
training set**, so that model has never really been tested against affinity. We
hedge across both rather than pick a winner.

---

## Slide 9 — What we would do next

1. **Novelty-aware generation.** Constrain sampling against the 525,063
   reference *during* design, not after. Our 30% rejection rate is pure waste,
   and the rejections are concentrated exactly where binding survives.
2. **Target Y105, not W99.** An alanine scan of the window says **Y105 is the
   anchor (0.87 kcal/mol)**, three times the conserved W99 (0.27). All 14 of our
   designs preserve it; only 1 of 14 preserves F104 (0.46), the second. *Nobody
   was looking at either.*
3. **Change format, not the loop.** Every single-position ΔΔG in this window is
   ≤0.87 kcal/mol — **at most ~4.3× in K_D, inside the assay floor.** No
   single-residue choice here can produce a measurable affinity change. CDRH3-tip
   redesign is the wrong lever. VHH scaffolds are intrinsically more thermostable
   and express far better in cell-free systems, where **expression is an unscored
   gate that zeroes you out entirely**.
4. **Calibrate per model.** Never transfer an ipTM threshold between folding
   models.

**The NB46 precedent — stated carefully.** Both parts confirmed in print: NB46
competitively blocks trastuzumab on HER2 *and* retains HER2 affinity after
exposure to **90 °C**. *(Zhao et al., J Med Chem, 2026-10-22; PMID 42911746,
doi:10.1021/acs.jmedchem.6c02664.)*

**Say "retains affinity after 90 °C exposure", never "Tm = 90 °C".** That figure
comes from a radiolabeling-robustness series, not a melting-temperature
measurement. **No Tm value for NB46 exists in the literature.** A judge who reads
the paper will catch the difference.

VHH thermostability generally: apparent Tm **60–80 °C**, some to 85 °C, versus
**47–69 °C** for scFvs *(Dumoulin et al., Protein Sci 2002; PMID 11847273)*.

**Three caveats to state before a judge finds them:**
- **Reversibility is largely a myth** — in a ~70-nanobody survey the large
  majority aggregate irreversibly on heat denaturation *(PMID 29784954)*.
- **Stability does not buy affinity** — explicitly uncorrelated across
  nanobodies from 18–100 °C *(PMID 39593207)*. Do not imply stability is a proxy
  for potency.
- **Stability is not predictable from framework sequence** — NbThermo's 564
  curated Tm values show no framework pattern; the variable loops dominate
  *(PMID 37042467)*.

**Clinical reality check:** all three anti-HER2 nanobodies that reached patients
(2Rs15d, NM-02/RAD201, ¹⁷⁷Lu-RAD202) are **radio-conjugates**. No naked
anti-HER2 VHH therapeutic is approved or in registrational trials. Pitch the
scaffold as a **theranostic / stability play**, not a Herceptin replacement.

---

## Slide 10 — Close

**We did not beat trastuzumab. We measured why the competition is hard, and we
know which of our numbers to trust.**

- The novelty gate and affinity are **structurally opposed** — 0% acceptance
  below 5 substitutions, and measured binding collapses above it.
- The standard confidence metric is **inside the noise of random pairings** for
  the model we fold with.
- Our 14 submissions are **novelty-verified** and **cleaner than the approved
  drug** on developability.
- One scoring axis we **could not address at all**, and we say so.

---

## Pre-flight, 17:40

- [ ] Every number traceable to a file in `findings/`
- [ ] Slide 5 presents **AUC first**, thresholds second
- [ ] Slide 7 shows the PhDs' actual pick
- [ ] No thermostability number anywhere
- [ ] No affinity claim below ~7×
- [ ] `_77` carries no score
- [ ] Chinery et al. (doi:10.1101/2024.03.26.586756) cited for both HER2 tables
- [ ] One person has read it aloud once, start to finish
