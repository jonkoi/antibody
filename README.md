# Trastuzumab CDRH3 redesign — Team 15

Oslo Antibodies × AI Hackathon, 14 November 2026.

Pipeline for redesigning trastuzumab's CDRH3 against HER2 domain IV, plus the
validation run that established it works end to end.

The competition scores a design on three axes against the measured trastuzumab
baseline — **affinity** (SPR), **epitope preservation** (competitive binding) and
**thermostability**. S-tier means improving on all three.

**Picking this up mid-stream?** Start with [HANDOFF.md](HANDOFF.md) — current state,
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
python3 pipeline/analyze_designs.py \
  --designs runs/her2-validation-seed1/designs \
  --reference inputs/12-her2/her2_fv_1n8z.cif \
  --scored runs/her2-validation-seed1/scored/refold_metrics.csv
```

Parses mmCIF directly (no external dependencies), diffs each design's designed
window against wild-type extracted from the input structure, and joins the
substitution counts to the refold metrics.

---

## Design constraints

Only **heavy chain positions 99–108** — a 10-residue window at the tip of CDRH3 —
were free to change. Held fixed: the entire light chain, the heavy framework and
CDRH1/CDRH2, HER2 domain IV backbone, and the relative pose of antibody and antigen.

The model was aimed at four HER2 contact atoms taken from PDB **1N8Z**:
Lys593 NZ, Phe573 CD1, Ile591 CG2, Pro572 CD. That hotspot set is the entire
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
3. **Confidence is not binding.** Calibrated against ~8,100 shuffled
   VHH–antigen pairs, Boltz-2's ipTM has ROC AUC **0.645**, and an ipTM ≥ 0.70
   cut-off admits **94.3%** of them. (An earlier version of this README claimed
   "454 SPR-measured non-binders, 21 at ipTM ≥ 0.7" — that set does not exist;
   see findings/lane-b-confidence.md.) Cut-offs should be calibrated rather than taken from raw
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
