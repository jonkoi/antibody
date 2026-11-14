# Decision brief for our structural biologists

**Oslo AI × Antibodies hackathon, Team 15. Submission closes 18:00 CET today.
We need your answers by 16:30.**

We redesigned trastuzumab's CDRH3 against HER2 domain IV. The computation is
done: 100 designs generated, novelty-filtered, folded, and scored by four
independent methods. **What we cannot do is decide.** Five questions below. Each
is a judgement call where our methods either disagree or fall below their own
resolution limit.

Nothing here needs you to read code. Section 6 lists what we plan to run next, so
you can redirect us before we spend the remaining time.

---

## 1. First, what we would get wrong without you

Three errors we have already caught by checking inherited claims against raw
data. We mention them so you calibrate how much to trust the rest:

- **Two of our four epitope hotspot residues were misnamed** in every document
  until an hour ago. `C/573/CD1` is **PHE 573**, not "Leu573"; `D/591/CG2` is
  **ILE 591**, not "Val591". The selections hit real atoms, so the designs were
  aimed at genuine coordinates — but the labels were wrong.
- A set of "454 SPR-measured non-binders" that we were calibrating against **does
  not exist** (the file has 29 rows, and its README forbids using them as
  negatives).
- `stab-ddg` is a **binding**-ddG predictor, not a thermostability predictor, so
  one of the competition's three scoring axes is unaddressed by anything we have.

---

## 2. The constraint that shapes everything

Designs are filtered by `abguard` against **525,063 known CDRH3 sequences**: a
submission must differ from every one of them by ≥3 mutations. Our 100 designs
hit this wall in a very specific way.

| substitutions from wild-type | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **passed the novelty gate** | **0%** | **0%** | **0%** | 33% | 82% | 74% | 93% | 100% |
| measured high-binder rate (FACS) | 98.2% | 90.9% | 72.7% | 58.5% | 42.2% | 31.7% | 28.5% | 24.5% |

**Not one conservative design survived.** Two were rejected for being *zero*
mutations from a sequence already in the measured library. The filter forces us
into the 5–9 substitution regime, where measured binding is least likely.

Measured affinity context, from the organisers' SPR data (n=187 with a KD):
wild-type trastuzumab is **0.095–0.135 nM**; affinity falls ~**1.69× per
mutation**; and of 150 variants at ≥3 mutations, **zero beat wild-type by a
margin the assay can resolve** (its own 8 replicates span **6.2×**, so the honest
detection floor is ~7×).

**We are therefore not trying to beat trastuzumab. We are trying to match it
within noise while clearing novelty.**

---

## 3. The four signals, and where they conflict

Each design was scored by four independent methods. Full table in section 5.

| method | what it measures | status |
| --- | --- | --- |
| **StaB-ddG** | binding ΔΔG vs wild-type, from structure | all 14 predicted *worse* than WT, 2.0–19.0× |
| **FACS model** | P(high binder), trained on 524k measured variants | test AUC 0.984; per-design ranking **not** seed-stable |
| **abguard liabilities** | developability (LAP reference) | all 14 **cleaner than trastuzumab itself** |
| **Boltz-2 ipTM / min_ipae** | fold confidence | **unusable for ranking** — see below |

**Two problems you should know about before reading the numbers.**

**(a) The two measured-data signals anti-correlate.** Spearman(ΔΔG, FACS) =
**+0.522** — the designs the FACS model likes best are the ones StaB-ddG says
bind worst. The most likely reason: **only 37 of 733 SPR sequences appear in the
FACS training set at all**, so the FACS model has essentially never been tested
against measured affinity, while StaB-ddG is a static-structure method applied to
5–9 simultaneous mutations, which is outside its validated regime. *We do not
know which to trust.*

**(b) Fold confidence is worthless here.** Against ~8,190 mismatched
antibody–antigen pairs, ROC AUC is AF3 0.754, **Boltz-2 0.645**, Chai-1 0.568.
Boltz-2 is what we folded with, and its *non-pair* median ipTM is 0.884 — our
designs sit at 0.82–0.90, inside the range of random pairs. We have excluded it
from selection.

---

## 4. The five decisions

### Q1 — Are these the right four epitope atoms? *(blocks any future work)*

We aimed the design model at four HER2 atoms from PDB **1N8Z**, and this set is
our *entire* operational definition of "targets the trastuzumab epitope":

| selector | residue (verified in the structure) |
| --- | --- |
| `C/572/CD` | **PRO 572** |
| `C/573/CD1` | **PHE 573** |
| `D/591/CG2` | **ILE 591** |
| `D/593/NZ` | **LYS 593** |

HER2 is split into chain C (489–580) and chain D (591–607) because the 581–590
loop is unresolved in 1N8Z. **Is this the trastuzumab epitope, or have we aimed
at the wrong patch?** Everything downstream inherits this.

### Q2 — Did the designed loops stay on the epitope, or drift?

A folding model will return a confident complex with the paratope on the wrong
surface — and we have just shown its confidence score cannot detect that. **This
needs eyes on structures.** 28 files, designed and refolded for each of the 14
candidates:

`runs/team15-claude-n100-seed1/shortlist-structures/her2_cdrh3_<N>_{designed,refolded}.cif`

### Q3 — Is the motif collapse real, or an artefact?

The surviving designs converged hard: **11 of 15 begin `ST`**, most are
`ST··V[DG]YAMD`, and two were byte-identical. Is an `ST`-initiated CDRH3 plausible
for HER2 binding, or is the model exploiting the novelty filter by abandoning the
conserved tryptophan at the base of the loop? **If the motif is wrong, all 14
candidates fail together** — our real diversity is far lower than 14.

### Q4 — W99 and Y105: which residue actually matters?

Wild-type CDRH3 is `WGGDGFYAMD`. Two findings that we would like checked against
your intuition:

- **W99 (the conserved Trp) is weaker than expected.** SPR shows W99-keeping
  variants binding 8.1× tighter — but **adjusted for mutation load that drops to
  2.43×** (p=0.015; it survives stratification at p=5×10⁻⁴, so it is real, but two
  of seven strata reverse). StaB-ddG independently puts W99→Ser at **0.21 ± 0.15
  kcal/mol ≈ 1.4×**. Both land near or below the assay's ~2.5× best-case
  resolution. And `FGSRGFYAYA` (0.47 nM, 5 mutations) binds well *without* it.
- **Y105 looks like the real hotspot.** Alanine scan of the whole window in the
  wild-type background: **Y105A = 0.87 kcal/mol**, three times W99A (0.27) and the
  largest effect measured. F104A is second at 0.46. **All 14 designs preserve
  Y105; only 1 of 14 preserves F104.**
- Caveat we cannot resolve: *every* single-position ΔΔG in this window is
  ≤0.87 kcal/mol (≤4.3× in Kd), i.e. inside the assay noise floor. Either no
  residue here is a dominant hotspot, or the method is systematically
  under-calling. **Which is it?**

### Q5 — Which designs do we submit? *(yours to decide)*

We have up to 4 slots (16 if the pitch places). We are deliberately **not**
choosing, because the two measured-data signals disagree and the choice is
biological.

---

## 5. The 14 candidates

All 14 are `team-validate` clean against the 525,063-sequence gate, class
`novel`, 97.4–97.8% identity to trastuzumab, minimum distance 3–4 mutations.
"ΔΔG" is predicted binding penalty vs wild-type (lower is better); "FACS" is
P(high binder); "liab" counts developability liabilities (**wild-type trastuzumab
scores 8**).

| design | CDRH3 window | subs | ΔΔG (fold worse) | FACS | liab | CDR3 Trp | notes |
| --- | --- | ---: | ---: | ---: | ---: | :-: | --- |
| **_61** | SLGSADYAFD | 6 | **0.42 (2.0×)** | 0.840 | **5** | — | **best on every signal** |
| _88 | SLSLGNYAFD | 6 | 0.64 (3.0×) | 0.069 | **5** | — | cleanest, FACS dislikes |
| _60 | WSYAYDYAMD | 5 | 0.72 (3.4×) | 0.133 | 7 | yes | keeps W99, breaks F104 |
| _4 | YWSLHYYAMD | 6 | 1.24 (8.2×) | 0.011 | 7 | yes | bottom on FACS |
| _6 | STSAYGYAMD | 6 | 1.27 (8.5×) | 0.108 | 7 | — | |
| _94 | WTAAVDYAHD | 6 | 1.31 (9.1×) | 0.863 | 6 | yes | keeps W99, high FACS |
| _71 | STSTYFYAFD | 6 | 1.35 (9.7×) | 0.572 | 6 | — | **only design keeping F104** |
| _75 | STFSVDYAMD | 6 | 1.38 (10.3×) | 0.092 | 6 | — | |
| _77 | STAAVDYAMD | 6 | 1.48 (12.1×) | *unstable* | 6 | — | **no usable score** |
| _15 | STSAVDYAMD | 6 | 1.50 (12.5×) | 0.062 | 7 | — | |
| _7 | STASVDYAMD | 6 | 1.53 (13.1×) | 0.819 | 6 | — | |
| _85 | STAAVGYAMD | 6 | 1.70 (17.7×) | **0.973** | 6 | — | FACS top-2, ΔΔG worst |
| _19 | STASVGYAMD | 6 | 1.72 (18.2×) | **0.998** | 6 | — | FACS top-1, ΔΔG worst |
| _96 | STSSVGYAMD | 6 | 1.75 (19.0×) | 0.904 | 7 | — | |

**The structure of the choice:** `_19`/`_85` are the best on measured binding data
and the worst on energetics. `_61`/`_88` are the reverse. **`_61` is the only
design in the top tier of both**, and is tied-cleanest on developability.

Only **three** designs (`_61`, `_88`, `_60`) are predicted to land within the
assay's ~7× noise floor of wild-type. Everything the FACS model ranks top is
predicted 18–19× worse.

---

## 6. What we plan to run next — redirect us if this is wrong

Time is the binding constraint: **analysis stops 16:30**, deck at 17:00, submit by
17:50. In parallel, unless you say otherwise:

| we will | why | abandon if |
| --- | --- | --- |
| **Not generate more designs** | the GPU dies tonight and novelty-filtering is the bottleneck, not generation | you think the `ST` motif (Q3) is wrong — then we need a fresh batch urgently |
| **Not attempt the de novo track** | needs different weights, a new batch and fresh novelty checks | — |
| **Not chase an affinity record** | zero of 150 measured variants beat WT resolvably | — |
| Pull literature evidence on VHH thermostability | thermostability is the one axis nothing predicts; we must argue it, not score it | you know a better precedent than nanobody NB46 |
| Build the pitch around the novelty/affinity trade-off | it is the strongest thing we found, and judging is 40% "judgment" | you disagree that it is real |

**What we will not claim**, and please hold us to this: no affinity improvement
smaller than ~7×; no thermostability number at all; no ipTM threshold without
naming the folding model; no per-design ranking that is not reproducible across
training seeds.

---

## 7. The fastest version, if you have ten minutes

1. **Q1** — are PRO572 / PHE573 / ILE591 / LYS593 the trastuzumab epitope? *(yes/no)*
2. **Q4** — is Y105 the real anchor rather than W99? *(yes/no)*
3. **Q3** — is `ST··V[DG]YAMD` a plausible HER2-binding CDRH3? *(yes/no)*
4. **Q5** — pick up to 4 from the table. If you want one: **`_61`**. If you want a
   hedge across the disagreeing signals: **`_61` + `_19` + `_88` + `_94`**.

Q2 (loop drift) is the one worth more than ten minutes if you have it, because no
metric we have can substitute for looking.
