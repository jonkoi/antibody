# Lane B — how far can we trust a folding model's confidence?

Done 14 Nov 2026. Scripts: [../pipeline/calibrate_confidence.py](../pipeline/calibrate_confidence.py),
[../pipeline/calibrate_replicas.py](../pipeline/calibrate_replicas.py). Both stdlib
only; re-run them to reproduce every number here.

**Answer in one line:** Boltz-2's ipTM — the metric our pipeline produces — barely
separates real antibody–antigen pairs from deliberately mismatched ones
(AUC 0.645), and at the ipTM our shortlist occupies, **77% of known non-binders
would pass**. Confidence cannot select our designs. Ranking on `min_ipae` instead
of ipTM was right.

---

## Correction first: the "454 non-binders" do not exist

Both [../README.md](../README.md) and [../EXECUTION.md](../EXECUTION.md) state
that the organisers ship *454 VHH–antigen pairs measured by SPR as true
non-binders, 21 of which still score ipTM ≥ 0.7*. Nothing in the bucket supports
this, and the number should not reach the pitch.

`hackathon-data/sigma2-ab-ag-champloo/spr-supplementary/` holds **29** rows, not
454 — 6 cognate and 23 non-cognate pairs, of which 4 have a fitted KD, 12 gave no
reportable fit and 12 were never measured. The string "454" appears nowhere in
`hackathon-data/README.md`.

More importantly those rows *cannot* carry a false-positive rate. The authors
state: **"N/A denotes no reportable affinity fit, NOT confirmed non-binding"**,
and the data README repeats it — *"do not use those rows as true negatives"*. Only
3 of the 29 rows carry an AF3 ipTM at all; the rest appear solely inside a figure
in the supplementary PDF.

What the bucket *does* have is `zenodo-18390239-ab-ag-complexes/iptm_confidence_scores.zip`:
ipTM matrices in which the nanobody panel is folded against the whole antigen
panel. Cells where the nanobody's PDB id equals the antigen's are **cognate**
(real) pairs; all others are **shuffled** pairs — presumed non-binding by the
study's design. That gives ~8,100–8,200 negatives per folding model instead of 29.

Call them shuffled pairs, never "SPR-measured non-binders".

---

## Part 1 — a single ipTM, three folding models

91 × 91 matrix per model (Smorodina et al., zenodo 18390239, CC BY 4.0).

| model | cognate n | shuffled n | ROC AUC (95% CI) | cognate median | shuffled median |
| --- | ---: | ---: | --- | ---: | ---: |
| AlphaFold3 | 91 | 8,137 | **0.754** [0.696, 0.813] | 0.500 | 0.180 |
| **Boltz-2** (ours) | 91 | 8,190 | **0.645** [0.584, 0.707] | 0.911 | 0.884 |
| Chai-1 | 89 | 7,952 | **0.568** [0.506, 0.630] | 0.208 | 0.191 |

CIs are Hanley–McNeil; with ~90 cognate pairs a gap of a few hundredths is not a
real difference. AF3 beating Boltz-2 is only marginally resolved (the CIs almost
touch), so treat the *ordering* as suggestive and the *level* as the finding.

False-positive rate of an ipTM cut-off — the share of shuffled pairs that pass,
with the share of real pairs kept in brackets:

| cut-off | AlphaFold3 | **Boltz-2** | Chai-1 |
| --- | --- | --- | --- |
| ipTM ≥ 0.70 | 3.0% (47% kept) | **94.3%** (99% kept) | 0.2% (8% kept) |
| ipTM ≥ 0.82 | 1.1% (38% kept) | **76.6%** (87% kept) | 0.0% (0% kept) |
| ipTM ≥ 0.90 | 0.4% (14% kept) | **39.8%** (62% kept) | 0.0% (0% kept) |
| ipTM ≥ 0.95 | 0.0% (0% kept) | 5.6% (22% kept) | 0.0% (0% kept) |

Three conclusions.

**1. Our shortlist's ipTM carries no information.** Boltz-2's shuffled pairs have
median ipTM 0.884 and an interquartile range of 0.826–0.920. The n=100 shortlist
sits at 0.82–0.90 — *inside* that range. Deliberately mismatched pairs score what
our designs score.

**2. A defensible Boltz-2 cut-off costs almost all the real pairs.** 5% false
positives needs ipTM ≥ 0.951, retaining only 18% of genuine binders.

**3. A threshold without a model name is meaningless.** ipTM ≥ 0.82 is a 1.1%
false-positive rate on AF3 and 76.6% on Boltz-2. The 0.7 figure the "21 of 454"
claim was being used to justify is the worst case: 94% of known non-binders clear
it under Boltz-2.

---

## Part 2 — does repeating the fold help? No.

If one score is noisy, averaging 50 should help. It does not.

`chai1_all_per_replica_confidence_metrics.zip`: 735,172 per-replica `.npz` files,
14,721 pairs (153 cognate, 14,568 shuffled), 50 replicas each. Per pair we
computed single / mean / max / min / std / mean−std of ipTM, mean
`aggregate_score`, and the inter-chain clash rate.

| feature | AUC (95% CI) | |
| --- | --- | --- |
| max (best-of-50) | 0.543 [0.496, 0.590] | spans 0.5 |
| mean | 0.542 [0.495, 0.589] | spans 0.5 |
| single fold | 0.539 [0.492, 0.586] | spans 0.5 |
| std (consistency) | 0.510 [0.464, 0.556] | spans 0.5 |
| clash rate | 0.500 | no signal |

Every confidence interval includes 0.5. Best-of-50 is indistinguishable from a
single fold, and both are indistinguishable from a coin flip on this data.

**Actionable:** do not spend GPU folding each design repeatedly in the hope that
averaged or best-of-N confidence will rank them. On this evidence it buys nothing.

---

## What this does not establish

- **Negatives are shuffled, not measured.** A shuffled pair is presumed
  non-binding because nobody reported it binding. The rates above are relative to
  that assumption.
- **Part 2 is Chai-1 only** — and Chai-1 is the weakest of the three. The Boltz-2
  per-replica archive is 1.7 TB and the AF3 shuffled one 163 GB, so neither is in
  the bucket. Whether replica statistics are equally useless for **Boltz-2**,
  which we actually fold with, is **unverified**.
- **Absolute thresholds do not transfer to our pipeline.** Our refold is handed
  the HER2 template with the framework and CDRH1/CDRH2 frozen, and only ten
  residues ever move, so our ipTM is not on the same scale as an unbound
  VHH–antigen prediction. What transfers is the model ordering and the shape of
  the trade-off — not a number to paste into our filter.
- **`min_ipae` itself is not calibrated here.** These matrices carry ipTM only.
  Calibrating a PAE-derived interface score (the paper uses ipSAE) needs the PAE
  matrices in the NIRD per-replica archives; the Chai-1 `.npz` files hold scalars
  only, no PAE.

## Reproducing

```bash
python3 pipeline/calibrate_confidence.py                 # seconds
python3 pipeline/calibrate_replicas.py                   # ~25 min, 735k files
python3 pipeline/calibrate_replicas.py --limit 20000     # quick check
python3 pipeline/calibrate_replicas.py --out runs/analysis/chai1_replica_features.csv
```

The `--out` table is per-pair features, so any re-analysis skips the long parse.

Cite Smorodina et al., [doi:10.64898/2026.03.02.709004](https://doi.org/10.64898/2026.03.02.709004)
(dataset [zenodo 18390239](https://doi.org/10.5281/zenodo.18390239),
[doi:10.11582/2026.owjnnp5q](https://doi.org/10.11582/2026.owjnnp5q)); AlphaFold3
outputs are non-commercial research only.
