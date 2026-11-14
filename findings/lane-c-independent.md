# Lane C — independent re-derivation (adversarial)

Run by a freshly spawned agent with **no access to our conclusions** and no sight
of STRATEGY.md, EXECUTION.md, PLAN.md, TASKS.md or findings/. Given the raw data
and five questions only. Stdlib-only scripts written from scratch.

**Verdict: our directional claims all hold. Three of our numbers were wrong or
overstated, and one of our deck slides was actively misleading.**

---

## 1. CONFIRMED, and our numbers were too generous

| our claim | independent finding |
| --- | --- |
| assay noise floor **6.2×** | **~7×** is the honest floor for single measurements |
| affinity falls with mutation count | confirmed, ρ=+0.483 (n=187), **+1.69× per mutation** |
| "2 of 153 beat wild-type (1.3%)" | **2 of 150 nominally — and ZERO by a resolvable margin** |
| ipTM barely separates real pairs | confirmed, and **externally validated** |

**The noise floor is worse than we said, and the wild-type replicates are
bimodal.** The 8 values split by printed decimal precision: five at 3 dp
(geometric mean 0.152 nM) versus three at 7–8 dp (gmean 0.044 nM) — **3.44×
apart**, which looks like two merged runs. Log10 SD 0.307 overall, 0.141
within-batch. So a two-variant difference needs **~7×** to be believable, or
~2.5× at best if both sat in the same batch. There is no run/date column, so the
batch reading **cannot be confirmed from this data** — but it is enough that
**any sub-7× affinity claim is unsafe.** Use ~7×, not 6.2×.

**Nothing beats trastuzumab by a margin this assay can resolve.** Of 150 variants
at ≥3 mutations with a KD, two measure nominally tighter: `WGSQGFYAYD` (3 mut,
0.0234 nM) and `WGGSGPYLFD` (4 mut, 0.095 nM). The first beats wild-type's best
replicate by **1.22× against a ~7× floor**; the second lands exactly on the
wild-type geometric mean. **Zero reproducible wins.** Only 10 of 150 even reach
KD < 1 nM. Our "1.3%" should become "nominally 1.3%, actually zero".

**The substitution table is survivorship-censored, so the cliff is steeper than
we showed.** Quantification rate falls from 100% at 0–1 mutations to 44% at 3,
26% at 6, and **3.8–5.6% at 9–10**. The apparent flattening at 8–10 mutations is
survivor bias — only the rare winners get a KD at all. The trustworthy part is
the first four mutations: **0.095 → 34 nM, ~360×**. Lead with that.

---

## 2. THE DECK SLIDE WE HAD WRONG

We were going to present, as a model comparison:

> Boltz-2 admits **94.3%** of known non-binders at ipTM ≥ 0.70, versus AF3's 3.0%.

**That ranking is backwards as a quality claim.** The three models live on
completely different ipTM scales, so a fixed threshold measures the scale, not
the model:

| model | ROC AUC (95% CI) | FPR @ 0.70 | TPR @ 0.70 | median, true pairs | median, non-pairs |
| --- | --- | ---: | ---: | ---: | ---: |
| **AF3** | **0.754** (0.696–0.813) | 3.0% | 47.3% | 0.500 | 0.180 |
| Boltz-2 | 0.645 (0.584–0.707) | **94.3%** | 98.9% | 0.911 | 0.884 |
| **Chai-1** | **0.568** (0.506–0.630) | **0.2%** | 7.9% | 0.208 | 0.191 |

**Chai-1 has the best false-positive rate and the worst AUC** — its positives top
out at 0.814, so a 0.70 cut-off rejects almost everything (TPR 7.9%). Boltz-2's
*non-pair* median is already 0.884, so 0.70 is below the noise.

**The correct claim, and it is a better one:** rank models by **AUC — AF3 0.754,
Boltz-2 0.645, Chai-1 0.568 (lower CI bound 0.506, barely above a coin flip)** —
and present the 94.3% as proof that **a fixed ipTM threshold is meaningless
across folding models**, which is the actual lesson. At a 5% FPR operating point
you recover 49.5% of real pairs with AF3, 17.6% with Boltz-2, 16.9% with Chai-1.
With a 1.1% base rate, AF3 at 0.70 gives 14.9% precision — six false positives
per true one.

Scale-free cross-check, same ordering: the cognate antigen ranks #1 of 91 for AF3
in 35.2% of cases, Boltz-2 5.5%, Chai-1 10.1%.

**External validation of the setup** (stronger than anything we had): the AF3
diagonal reproduces all three cognate ipTM values printed in the supplementary
PDF — 8W90 = 0.12, 8AV2 = 0.15, 9ETJ = 0.46 — to the stated digit, and 11 of the
12 non-cognate SPR pairs land in 0.75–0.90 exactly as the PDF claims. The matrix
orientation and the diagonal-is-cognate convention are therefore confirmed by
data outside the matrix itself.

**Also: the matrix is 91×91 — 91 true pairs, ~8,190 non-pairs — not the "106 real
/ 11,342 shuffled" the README advertises.** About 16 systems are missing from the
`_clean` files. Chai-1 is missing 240 cells including 2 diagonals (5E7F, 8E0E),
AF3 is missing 53. Dropping NAs silently changes the denominator per model.

---

## 3. W99: REAL BUT WEAK — this corrects us in both directions

We had swung from "W99 is decisive" to "W99 is not the discriminator". Neither is
right.

| | raw | adjusted for mutation count |
| --- | ---: | ---: |
| W99-keeping vs losing, gmean KD | **8.11×** | **2.43×** |
| significance | p = 2.0e-08 | p = 0.015 |

**The effect survives stratification** (van Elteren z = −3.46, **p = 5.4e-04**),
but roughly half the raw advantage was mutation-load confounding: AUC falls from
0.742 to 0.688. Only 2 of 7 strata are individually significant, and **two strata
point the other way**. `FGSRGFYAYA` (0.471 nM, 5 mutations) proves W99 is not
required.

**So W99 is worth ~2.4–3×, which sits barely above the ~2.5× best-case assay
resolution. It is a weak-to-moderate prior, not a hard rule.**

This now *agrees* with StaB-ddG, which put W99→Ser at 0.207 ± 0.148 kcal/mol
(≈1.4×). Two independent methods both say: small, possibly real, not decisive.
That is the honest answer to Q0, and it is neither of the two positions we held.

---

## 4. A caveat on the FACS model we had missed

**The FACS and SPR sets barely overlap: only 37 of 733 SPR 10-mers (27 of the 180
with a KD) appear in the 524,346-row FACS table.** So a FACS-trained classifier
is *essentially untested against measured affinities* — it predicts a binding
*class* from cell-surface display, validated only against more of the same assay.

This is the most plausible explanation for the ddG/FACS anti-correlation
(Spearman +0.522): the two signals are not measuring the same quantity, and only
one of them has ever been checked against SPR. **Neither should be presented as
the authoritative ranking.**

---

## 5. Other traps worth knowing

- **Three CDRH3 lengths in the SPR file**: 740 ten-mers, **20 nine-mers**, 8
  sixteen-mers. `zip()` silently truncates, so Hamming distance is wrong for the
  9-mers; 3 of them carry a KD and would land in a bogus mutation bin.
- **A different antibody hides in the file**: 8 identical rows with CDRH3
  `ARDRIVDDYYYYGMDV` and the only other light chain, no KDs. Filter on
  `vl_sequence`, not CDRH3 length.
- **"No reportable fit" is not a non-binder** — and the same applies to the 460
  `no_reportable_kd` and 104 `inconclusive` rows, and to the FACS `low` class.
  None are usable true negatives. (We already corrected this.)
- **Some AF3 "non-pairs" may be genuine pairs** — e.g. 8E0E×8DAM at 0.93/0.91
  looks like one antigen deposited twice. But the 246 high off-diagonal cells
  span 71 of 91 nanobodies, so leakage explains only a handful; the rest is real
  non-specific confidence.
- **FACS data is unusually clean**: all 524,346 CDRH3s distinct, all exactly 10
  residues, `edit_distance` a perfect Hamming match in every row, wild-type
  absent. Safe to trust.
