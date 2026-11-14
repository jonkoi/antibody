# Lane D (W1) — a ranking signal grounded in measured binding

Done 14 Nov 2026. Script: [../pipeline/predict_binding.py](../pipeline/predict_binding.py),
stdlib only. Re-run it to reproduce every number here.

**Why:** findings 7 and 8 ([lane-b-confidence.md](lane-b-confidence.md)) ruled out
folding confidence as a way to rank this pool, and `min_ipae` is calibrated
against nothing measured. This is the only ranking signal we have that rests on
experiment.

**Answer in one line:** the model predicts measured FACS binding class well and
reproducibly (test AUC 0.9842, macro-AUC 0.9306 across clusters), but it does
**not** yield a seed-stable per-design ranking. Two designs separate robustly;
the rest resolve only into tiers.

---

## Data

Mehta et al. ([zenodo 10831512](https://doi.org/10.5281/zenodo.10831512), CC BY
4.0) FACS-sorted trastuzumab CDRH3 variants. The authors' clonotype-aware splits
are used **as shipped** — similar CDRH3s stay inside one split, so held-out
performance is not inflated by near-duplicates.

| split | n | positives | clusters |
| --- | ---: | ---: | ---: |
| train | 367,044 | 120,505 (32.8%) | 13,829 |
| val | 78,651 | 25,822 (32.8%) | 251 |
| test | 78,651 | 25,822 (32.8%) | 706 |

Exact CDRH3 overlap between train and test: **0**.

Label is 1 = high binder. It is a **FACS binding class, not an affinity**, so it
cannot rank designs by KD — and the competition scores affinity by SPR.

## Model

Logistic regression on the 10-residue window, AdaGrad, epoch chosen on
validation, test split scored once at the end.

- unary: position × amino acid (200 weights)
- pairwise: position pair × amino acid pair (18,000 weights)

Only 55 features are active per sequence, so an epoch over 367k rows takes ~40 s
without numpy.

## Performance, against baselines

| model | test AUC (95% CI) |
| --- | --- |
| majority class | 0.500 |
| edit distance from wild-type alone | 0.7534 [0.7496, 0.7573] |
| position-weight matrix (independent positions) | 0.9209 [0.9186, 0.9233] |
| logistic regression, unary only | 0.9437 [0.9417, 0.9457] |
| **logistic regression, unary + pairwise** | **0.9842 [0.9831, 0.9853]** |
| — macro-AUC over test clusters with both classes | **0.9306** (44 of 706) |

**Quote 0.9306, not 0.9842, for ranking designs.** 662 of 706 test clusters are
single-class, so the row count overstates the effective sample size and the
overall AUC is mostly *between-family* discrimination — the easy problem. The
macro figure is within-family discrimination, which is what ranking our designs
actually requires.

Two shortcuts ruled out:

- **Not just counting mutations.** Distance from wild-type alone reaches 0.7534,
  far short of 0.9842.
- **The interactions earn their keep.** +0.04 over unary-only, +0.06 over the PWM.
  In the unary-only ablation the shortlist's scores collapse to a flat 0.13–0.65.

Control: trastuzumab's own window `WGGDGFYAMD` scores **P(high) = 1.000** in every
run.

## Our designs are novel but densely covered

Training-set coverage by substitutions from trastuzumab's window:

| subs | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| train rows | 226 | 596 | 2,982 | 17,258 | 59,392 | 115,589 | 107,389 | 49,399 | 14,213 |
| % of train | 0.1% | 0.2% | 0.8% | 4.7% | 16.2% | 31.5% | 29.3% | 13.5% | 3.9% |
| high-binder rate | 98.2% | 90.9% | 72.7% | 58.5% | 42.2% | 31.7% | 28.5% | 24.5% | 21.8% |

Our accepted designs carry 5–9 substitutions — **95.2% of all training rows**.
Every shortlisted design sits at Hamming distance 3–4 from its nearest measured
CDRH3, none closer. So they are far from any individual measured sequence while
squarely inside the measured distribution: novel sequences in a densely covered
regime, which is the good case and not the extrapolation initially feared.

### The novelty trap, from two independent sources

Measured high-binder rate falls monotonically as substitutions rise, while
abguard acceptance rises monotonically over the same axis (acceptance from the
n=100 batch):

| subs from WT | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| novelty acceptance | 0% | 0% | 0% | 33% | 82% | 74% | 93% | 100% |
| measured high-binder rate | 98.2% | 90.9% | 72.7% | 58.5% | 42.2% | 31.7% | 28.5% | 24.5% |

Two curves running in opposite directions, from the training data and the novelty
gate respectively. The gate forces designs into the regime where measured binding
is **least** likely — base rate ~28–32% where our accepted designs live. A design
scoring 0.95+ there is real lift over that base rate. It also explains the 0%
acceptance at ≤4 substitutions as the same trade-off seen from the other side.

## Seed stability — why there is no per-design ranking

Three independent training runs (seeds 1–3), identical in every other respect.
Performance is reproducible: test AUC **0.9842 / 0.9844 / 0.9849**, wild-type
control 1.000 every time. The per-design ordering is not: Spearman between seed
rankings is only **0.881, 0.908, 0.881**.

| design | CDRH3 | seed 1 | seed 2 | seed 3 | mean | spread |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| her2_cdrh3_19 | STASVGYAMD | 0.998 | 0.999 | 0.997 | 0.998 | 0.002 |
| her2_cdrh3_85 | STAAVGYAMD | 0.958 | 0.968 | 0.992 | 0.973 | 0.034 |
| her2_cdrh3_96 | STSSVGYAMD | 0.923 | 0.942 | 0.847 | 0.904 | 0.095 |
| her2_cdrh3_94 | WTAAVDYAHD | 0.755 | 0.884 | 0.951 | 0.863 | 0.196 |
| her2_cdrh3_61 | SLGSADYAFD | 0.652 | 0.956 | 0.911 | 0.840 | 0.304 |
| her2_cdrh3_7 | STASVDYAMD | 0.617 | 0.890 | 0.951 | 0.819 | 0.334 |
| her2_cdrh3_71 | STSTYFYAFD | 0.514 | 0.586 | 0.617 | 0.572 | 0.103 |
| her2_cdrh3_77 | STAAVDYAMD | 0.229 | 0.234 | 0.790 | 0.418 | 0.561 |
| her2_cdrh3_60 | WSYAYDYAMD | 0.022 | 0.073 | 0.303 | 0.133 | 0.281 |
| her2_cdrh3_6 | STSAYGYAMD | 0.034 | 0.042 | 0.247 | 0.108 | 0.213 |
| her2_cdrh3_75 | STFSVDYAMD | 0.009 | 0.108 | 0.158 | 0.092 | 0.149 |
| her2_cdrh3_88 | SLSLGNYAFD | 0.017 | 0.182 | 0.008 | 0.069 | 0.174 |
| her2_cdrh3_15 | STSAVDYAMD | 0.024 | 0.039 | 0.123 | 0.062 | 0.099 |
| her2_cdrh3_4 | YWSLHYYAMD | 0.012 | 0.012 | 0.010 | 0.011 | 0.002 |

`her2_cdrh3_77` moves 0.229 → 0.790 under nothing but a reshuffled training
order. A number that unstable is not evidence about a design.

### What may be reported

Tiers, not a ranking:

- **Stable high — `her2_cdrh3_19`, `her2_cdrh3_85`.** Means 0.998 and 0.973 with
  spreads of 0.002 and 0.034. Top two in all three runs. These separate robustly.
- **High but unordered — `_96`, `_94`, `_61`, `_7`.** Means 0.82–0.90, spreads up
  to 0.33. All plausibly good; their order among themselves is not reproducible.
- **Mid — `_71`.** Mean 0.572.
- **Unusable — `_77`.** Spread 0.561. Report no number for it.
- **Stable low — `_4`** (0.011, spread 0.002), with `_15`, `_88`, `_75`, `_6`,
  `_60` all below 0.14 on average.

So: a model that generalises well, a measured-data basis for preferring two
designs, and no licence for a 14-way ordering.

## What this does not establish

- **FACS class is not affinity.** It cannot rank by KD, and the competition
  measures affinity by SPR.
- **No thermostability.** Nothing here speaks to it; one of the three scoring axes
  stays unaddressed.
- **Within-family discrimination rests on 44 clusters.** The honest macro-AUC of
  0.9306 comes from the only test clusters that contain both classes.
- **The designs are extrapolations at the sequence level** even though their
  substitution-count regime is dense: none is closer than Hamming 3 to a measured
  CDRH3, so each prediction is a generalisation, not a lookup.

## Reproducing

```bash
# the authors' splits, extracted from zenodo-10831512-trastuzumab-ngs/data.zip
# path: data/affinity_data/her2/her2_aff_large/processed/remove_overlap/
#       clonotype_split/0.7_0.15_0.15/{train,val,test}.csv
python3 pipeline/predict_binding.py                    # ~4 min, 6 epochs
python3 pipeline/predict_binding.py --unary-only       # ablation
python3 pipeline/predict_binding.py --seed 2           # stability
```

Cite Mehta et al., [zenodo 10831512](https://doi.org/10.5281/zenodo.10831512),
and Chinery et al., [doi:10.1101/2024.03.26.586756](https://doi.org/10.1101/2024.03.26.586756)
for the FACS and SPR tables.
