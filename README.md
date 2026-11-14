# Trastuzumab CDRH3 redesign — Team 15

Oslo Antibodies × AI Hackathon, 14 November 2026.

Pipeline for redesigning trastuzumab's CDRH3 against HER2 domain IV, plus the
validation run that established it works end to end.

The competition scores a design on three axes against the measured trastuzumab
baseline — **affinity** (SPR), **epitope preservation** (competitive binding) and
**thermostability**. S-tier means improving on all three.

**Picking this up mid-stream?** Start with [HANDOFF.md](HANDOFF.md) for state and next
actions, then [log0.md](log0.md) (data analysis, pipeline code) and [log1.md](log1.md)
(platform survey, external thermostability resources). Dead ends are recorded there —
read them before repeating work.
open questions, next actions and gotchas.

---

## Pipeline

Four stages. Stage 3 is the load-bearing one: the folding model never sees what
the design model intended, so agreement between them is independent evidence
rather than the design model marking its own homework.

| Stage | Tool | What it does |
| --- | --- | --- |
| 1. design | RFD4-Proteina (`module-run design`) | Proposes CDRH3 backbone + sequence, everything else frozen |
| 2. prepare | `module-run refold-prepare` | Selects top-N, writes folding inputs incl. HER2 template |
| 3. fold | Boltz-2 / Protenix / OpenFold3 / Mol-1 | Independently re-folds each design with HER2 |
| 4. score | `module-run refold-score` | Compares refold to design → confidence metrics |

```bash
# on the hackathon VM, from ~/rfd4 (holds Data/ and Settings/)
N=100 TOPN=20 SEED=1 bash pipeline/run_design.sh all

# sweep a second folding model over the same designs
FOLDER=protenix-v2 bash pipeline/run_design.sh predict score
```

`run_design.sh` generalises the organisers' 4-design example into a parameterised
batch driver and adds a fold-model sweep. Weights, design count, seed, sequence
randomness (`TSRK`) and folding model are all environment knobs.

### Analysis

```bash
# sequence diff + refold metrics
python3 pipeline/analyze_designs.py \
  --designs runs/her2-validation-seed1/designs \
  --reference inputs/12-her2/her2_fv_1n8z.cif \
  --scored runs/her2-validation-seed1/scored/refold_metrics.csv

# liabilities, biophysics and loop rigidity vs wild-type
python3 pipeline/score_developability.py \
  --designs runs/her2-validation-seed1/designs \
  --reference inputs/12-her2/her2_fv_1n8z.cif
```

Both parse mmCIF directly with no external dependencies. The first diffs each
design's window against wild-type and joins to refold metrics; the second scores
chemical liabilities and biophysical character.

---

## The findings that drive selection

### 1. Binding collapses with mutation count, and the novelty floor sits on the knee

From the 524,346-variant FACS table, high-binder rate against CDRH3 edit distance
from trastuzumab:

| mutations | 1 | 2 | **3** | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| high binders | 96.0% | 97.5% | **90.0%** | 74.6% | 59.4% | 43.9% | 33.6% | 26.8% | 18.9% | 15.8% |

One and two mutations bind best but **fail the novelty gate**. Three is the tightest
design that passes, and it retains 90% high-binder rate. Every further mutation costs
roughly 15 percentage points.

The validation run produced designs at 3, 7, 7 and 9 mutations — three of four in the
19–34% regime when they could have been at 90%. **Nothing in the generation step
constrains mutation count, and it should.**

The constraint is harder than "mutate exactly 3 positions": abguard requires ≥3
substitutions from trastuzumab *and* from every one of the 525,063 provided CDRH3s,
so a candidate that is a one-residue extension of a known 2-mutant is rejected.
Selecting 3-mutants is a constrained search over the provided set, not sampling.

### 2. Wild-type trastuzumab CDRH3 is itself a developability liability

`WGGDGFYAMD` carries an **Asp-Gly isomerization hotspot** (the worst chemical
liability class), a Met oxidation site, and **three glycines in ten residues** — an
unusually flexible loop.

This makes the thermostability axis more winnable than it first appears, provided
designs are selected for it. Almost any redesign of positions 102–103 removes the DG,
and any substitution away from Gly rigidifies the loop.

In the validation run, **design_2 (`WGPPYFYAMD`) replaced two glycines with two
prolines** — the textbook loop-thermostabilisation move, reducing unfolded-state
backbone entropy. It is simultaneously at the novelty floor (3 mutations, 90%
high-binder regime) and liability-reduced. That combination is the selection target.

**Caveat, stated plainly:** removing a DG hotspot improves *chemical* stability, not
melting temperature. Only the Gly→Pro rigidification bears directly on the measured
thermostability axis. And a rigidified loop can equally lock a non-binding
conformation — that needs structural review, not a score.

### 3. Only 6.2% of the minimal-change space is legal

`pipeline/enumerate_candidates.py` enumerates the complete distance-3 set exactly
and filters it against every provided CDRH3. Pure standard library, 5.6 seconds.

| | count |
| --- | --- |
| candidates at distance 3 | 823,080 |
| **clear the novelty gate** | **51,417 (6.2%)** |
| ...also remove the DG hotspot | 29,124 |
| ...also place Pro at a wild-type Gly | 3,842 |
| ...both | 2,652 |
| ...both, retaining the conserved Trp | **1,927** |

**93.8% of minimal-change designs are illegal** — the 525,042 provided sequences
densely blanket trastuzumab's neighbourhood. A generative model sampling this space
blind wastes most of its output, because nothing tells it where the forbidden
regions are. Enumerating first and folding only legal candidates is strictly better
use of GPU time.

The position profile of the 1,927 shortlist converges on the same levers the
developability analysis identified independently: position 5 is substituted in all
of them, and the proline lands at position 5 (1,412) or position 3 (565) — two of
the three wild-type glycines.

Brute force would be 823,080 x 525,042 comparisons. Three observations collapse it:
the distance to a reference splits into a part independent of the candidate plus a
part confined to the three mutated positions; that first part depends on the
reference only through a 10-bit mask, so it is computed 1,024 times per position set
rather than 525,042; and each surviving reference generates its forbidden
combinations directly instead of being compared against all 6,859 candidates.

### 4. Do not rank on Boltz-2 interface confidence

`pipeline/calibrate_iptm.py` scores the provided 91x91 ipTM matrices, where the
diagonal holds real nanobody-antigen pairs and everything off it is a shuffled
non-binder (base rate 1.1%).

| model | ROC AUC | passes at ipTM > 0.7 | top-10 ranked |
| --- | --- | --- | --- |
| AlphaFold3 | **0.754** | 236 non-binders; rejects 48 of 91 real | **51.8% real (47x)** |
| Boltz-2 | 0.645 | **7,723 non-binders — 94.3% of them** | **0% real** |
| Chai-1 | 0.568 | 17 | finds only 7 of 89 real |

**Boltz-2's top 10 contained no real binders at all**, and its median confidence for
pairs that demonstrably do not bind is 0.884. Our validation run used Boltz-2, and
its 0.946 / 0.965 scores are matched by 7.5% of known non-binders — a likelihood
ratio of 3.4. AF3's maximum across the whole benchmark is 0.93, so those numbers are
not even on its scale.

AF3 is a poor classifier but a useful *ranker*: 31-fold enrichment in its top 100.
Apheris ships **OpenFold3-OpenBind**, an AF3 reproduction, and Protenix v2. Refold
with those rather than Boltz-2, use confidence for relative ordering within one
target, and never as a pass/fail gate.

**Caveat:** that benchmark is nanobody-antigen pairs folded without templates, while
we fold a CDRH3-redesigned Fv against HER2 *with* a template, and OpenFold3-OpenBind
reproduces AF3 rather than being it. The transfer is plausible, not demonstrated.

### Why there is no learned stability model here

Nebius `stab-ddg` computes *binding* ΔΔG, not folding stability. FoldX and CamSol are
licensed and unavailable. No dataset in the challenge carries a melting temperature —
FACS has binding classes, SPR has `kd_nM` only. **A thermostability predictor cannot
be trained on provided data**, so `score_developability.py` reports mechanistic risk
with ordinal weights and explicitly does not predict Tm.

---

## Design constraints

Only **heavy chain positions 99–108** — a 10-residue window at the tip of CDRH3 —
were free to change. Held fixed: the entire light chain, the heavy framework and
CDRH1/CDRH2, HER2 domain IV backbone, and the relative pose of antibody and antigen.

The model was aimed at four HER2 contact atoms taken from PDB **1N8Z**:
Lys593 NZ, Leu573 CD1, Val591 CG2, Pro572 CD. That hotspot set is the entire
operational definition of "target the trastuzumab epitope" in this run.

HER2 is split into chains C (489–580) and D (591–607) because input validation
rejects a numbering gap inside one chain — the 581–590 loop is unresolved in 1N8Z.

---

## Validation run

`runs/her2-validation-seed1/` — 4 designs, seed 1, `her2-novel-2026-11-09`
weights, top 2 refolded with Boltz-2. ~4 min wall on one H100.

Wild-type window (B99–108): `WGGDGFYAMD`

| design | window | subs | identity | min_ipae ↓ | iptm ↑ | binder_pLDDT ↑ |
| --- | --- | --- | --- | --- | --- | --- |
| her2_cdrh3_2 | `WGPPYFYAMD` | 3 | 70% | 0.331 | 0.946 | 95.6 |
| her2_cdrh3_3 | `WASSSEYYFD` | 7 | 30% | 0.310 | 0.965 | 98.5 |
| her2_cdrh3_1 | `DSPAAGYAFD` | 7 | 30% | not refolded | | |
| her2_cdrh3_0 | `DLSLSDFSFD` | 9 | 10% | not refolded | | |

Two observations from four designs — both need a real batch to confirm:

**The conserved Trp tracks selection.** Designs 0 and 1 both lost W99; designs 2
and 3 kept it, and those are exactly the two `refold-prepare` ranked into the top
2. W99 anchors the base of CDRH3.

**design_2 sits on the novelty knife-edge.** Three substitutions is exactly
abguard's floor for the optimisation track, and 70% window identity is exactly the
boundary the task description words as "under 70%". Designs this conservative need
checking with abguard rather than by arithmetic — and the IMGT CDRH3 that abguard
scores is not necessarily the same 10-residue window the model was given.

### Confidence metrics are saturated here

iptm 0.95–0.97 and pLDDT 96–99 are not evidence of excellent designs. The refold
was handed HER2 as a template, framework and both other CDRs were frozen at
wild-type, and only ten residues ever moved — the folding model is largely
reproducing a complex it was given most of the answer to.

Consequences for selection:

1. **iptm and pLDDT will not discriminate** across a large pool. Ranking several
   hundred designs that all score 0.94–0.97 is ranking on noise.
2. **min_ipae has some spread** (0.310 vs 0.331) but two points is not a
   distribution.
3. **Confidence is not binding.** The organisers ship 454 VHH–antigen pairs
   measured by SPR as true non-binders, 21 of which still score ipTM ≥ 0.7.
   Cut-offs should be calibrated against that set rather than taken from raw
   confidence.

---

## Novelty gate

Every submission is filtered by `abguard`, which assigns the track itself —
teams do not choose it. From `abguard --help`:

- identity to nearest known antibody **≥ 70%** → optimisation track → CDRH3 must
  differ from trastuzumab *and from every reference CDRH3* by **≥ 3** substitutions
- identity **< 70%** → de novo track → **≥ 2** substitutions from trastuzumab

abguard ships its own reference of **610,109** CDRH3s and uses it by default, so
no reference needs building. Also enforced: VHH/scFv/Fab only, ≤ 250 aa, no
unpaired cysteines, IMGT-numberable at ≥ 0.1 confidence.

---

## Layout

```
pipeline/
  run_design.sh            batch driver (design → prepare → fold → score)
  run-her2.upstream.sh     organisers' 4-design example, kept for reference
  analyze_designs.py       sequence diff + metric join
  score_developability.py  liabilities, biophysics, loop rigidity vs wild-type
inputs/12-her2/            design spec and 1N8Z-derived input structure
runs/her2-validation-seed1/
  designs/                 4 generated structures
  design_metrics/          per-design generation metrics
  prepared/                refold inputs + HER2 template
  refold_boltz2/           Boltz-2 outputs (structures, scores, PAE)
  scored/                  refold_metrics.csv
```

## Environment

One H100 80GB on Nebius. Every Foundry model is a local Docker image — there is no
Foundry Hub on the VM, so `foundry-toolkit` is not used and the images are invoked
with `docker run` directly.

Inputs under `inputs/` and `run-her2.upstream.sh` come from the hackathon bucket
(`s3://apheris-oslo-hackathon/`) and are the organisers' material, included for
reproducibility.
