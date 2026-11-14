# Trastuzumab CDRH3 redesign — Team 15

Oslo Antibodies × AI Hackathon, 14 November 2026.

**Submission: [SUBMISSION.csv](SUBMISSION.csv)** — 4 designs, `abguard
team-validate` clean (**4/4 valid**) against all 525,063 reference CDRH3s.
**Pitch: [DECK.md](DECK.md)** · **Evidence: [findings/](findings/)**

---

## What we set out to do, and what actually happened

Redesign the 10-residue tip of trastuzumab's CDRH3 (heavy chain 99–108) to bind
HER2 domain IV better, while differing from every one of 525,063 known CDRH3s by
≥3 mutations. Judged on **affinity**, **epitope preservation** and
**thermostability**.

We built a generative pipeline, ran 100 designs, found that the approach **does
not work**, then found out why — and what does.

### 1. The novelty filter and binding affinity pull in opposite directions

| substitutions from wild-type | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| passed the novelty gate (our 100) | **0%** | **0%** | **0%** | 33% | 82% | 74% | 93% | 100% |
| measured high-binder rate (524,346 variants) | 98.2% | 90.9% | 72.7% | 58.5% | 42.2% | 31.7% | 28.5% | 24.5% |

Not one of our 13 conservative designs survived. Two were rejected for being
*zero* mutations from an already-measured sequence — the model re-invented known
variants.

### 2. But we over-read that, and the correction is the real result

Five designs produced by a coworker, **enumerated rather than generated**, sit at
**3 substitutions and all five pass the same gate**. So 0% acceptance was a
property of **our generator's output distribution**, not of the design space:
RFD4-Proteina scatters mutations across the binding interface, which is exactly
where the measured library samples densest.

**Generate-then-filter cannot find the conservative region. Constraint-aware
enumeration can.**

### 3. Only 4 of the 10 designable residues touch HER2

Computed from 1N8Z: **W99, G103, F104, Y105**. The other six contact nothing, so
≥3 mutations can be placed entirely on them — 2,091,995 such sequences exist at
3–4 substitutions. Our generative designs preserved at most **2 of 4** contacts;
zero of 71 preserved three or more. The enumerated set preserves **3 of 4**.

### 4. Alanine scanning is blind to the interface's tightest contact

Our ΔΔG scan rated G103A at 0.042 kcal/mol — essentially zero. But G103's
**backbone carbonyl** makes the **2.42 Å hydrogen bond to Lys593, the tightest
contact in the whole interface.** Alanine retains the backbone carbonyl. We would
have gone on mutating it.

### 5. Nothing beats trastuzumab by a margin the assay can resolve

Wild-type is measured 8 times spanning **6.2×**, and those replicates are bimodal
(3.44× apart by printed precision, suggesting two merged runs), so the honest
detection floor is **~7×**. Of 150 measured variants at ≥3 mutations, two are
nominally tighter and **zero** beat it resolvably.

### 6. Boltz-2's ipTM cannot rank these designs

ROC AUC against ~8,190 mismatched pairs: AF3 0.754, **Boltz-2 0.645**, Chai-1
0.568. Boltz-2's *non-pair* median ipTM is 0.884, and our designs sat at
0.82–0.90 — inside the range of random pairings. A fixed threshold measures the
model's scale, not its quality: Chai-1 has the best false-positive rate and the
worst AUC.

### 7. Thermostability is unaddressed

Verified across **all 59** Nebius models and the Apheris stack: nothing predicts
it. `stab-ddg` is StaB-ddG, which scores *binding* ΔΔG. One of the three judged
axes has no number from us, and we say so rather than substitute a proxy.

---

## The submission

| # | design | CDRH3 | subs | contacts kept | ΔΔG (fold vs WT) | liab | role |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| 1 | cw_3 | `WGGMGPYWMD` | 3 | 3/4 | **0.17 ± 0.17 (1.3×)** | 8 | best energetics |
| 2 | cw_2 | `WGGKGPYWMD` | 3 | 3/4 | 0.42 (2.0×) | 8 | best predicted KD |
| 3 | her2_cdrh3_71 | `STSTYFYAFD` | 6 | 2/4 | 1.35 (9.7×) | 6 | **only design retaining Phe104** |
| 4 | her2_cdrh3_61 | `SLGSADYAFD` | 6 | 1/4 | 0.42 (2.0×) | **5** | independent pattern, cleanest |

ΔΔG is StaB-ddG binding penalty — a model neither design approach trained. cw_2
and cw_3 land inside the ~7× noise floor, and cw_3 is statistically
**indistinguishable from wild-type**. All four pass the hard expression filters
(zero unpaired cysteines, zero N-glycosylation sequons); approved trastuzumab
itself scores 8 liabilities.

**Refold check.** All five enumerated designs were folded with Boltz-2 against
the HER2 template alongside a wild-type control. CDRH3 stays **2.2–3.0 Å** from
the four epitope hotspots (wild-type control: 2.63 Å), so there is no loop drift,
and backbone phi at residue 104 lands at **−59° to −78°** — the proline-compatible
range — so `F104P` does not clash with the loop.

**Known risk, stated plainly.** cw_2 and cw_3 share `F104P` and `A106W`. F104 is
peripheral to HER2 (nearest atom **4.27 Å**) but packs against the light chain
(Tyr49/Leu46/Tyr55), and `A106W` gives each design a second CDR3 tryptophan — a
documented oxidation and aggregation hotspot (PMID 30608647, 31678251). Slots 3
and 4 therefore hold designs that fail *independently*, slot 3 being the only
design we have that retains Phe104.

---

## Reproducing it

```bash
# novelty reference: 525,063 CDRH3s from the two HER2 tables
abguard reference-ingest her2_novelty.duckdb \
  hackathon-data/her2-trastuzumab/tz_her2_facs_3classes_524346.csv.zip \
  hackathon-data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv
abguard reference-build her2_novelty.duckdb

# the gate. --reference-db is REQUIRED: without it prior art is never checked,
# and `abguard novelty-report` ignores the flag entirely.
abguard team-validate SUBMISSION.csv --reference-db her2_novelty.duckdb

# analyses, stdlib only, no numpy/pandas
python3 pipeline/analyze_spr.py           # assay floor, cost per substitution
python3 pipeline/calibrate_confidence.py  # ipTM ROC AUC per folding model
python3 pipeline/predict_binding.py       # FACS binding model, test AUC 0.984
python3 pipeline/cif_to_pdb.py --help     # mmCIF → PDB for StaB-ddG
python3 pipeline/make_submission.py --help
```

Pipeline: `design` (RFD4-Proteina) → **novelty gate** → `fold` (Boltz-2) →
`score`. **Run the gate before folding** — the default `refold-prepare` ranking
spent 30% of our GPU on designs the gate later rejected, all in the substitution
range where binding does not survive.

---

## Layout

```
SUBMISSION.csv          the 4 designs, team-validate clean
DECK.md                 the pitch, every number frozen and sourced
findings/               the evidence, one file per lane
pipeline/               tooling, stdlib only
inputs/12-her2/         organisers' design spec and 1N8Z-derived structure
inputs/12-her2-hotspot/ variant spec pinning F104 and Y105
runs/                   run artefacts and result tables
STRATEGY.md   the argument       EXECUTION.md  the commands
PLAN.md       who does what      CONNECT.md    credentials and the VM
HANDOFF.md log0.md log1.md       earlier state and recorded dead ends
```

---

## Corrections we made to our own work

Each was found by checking a claim against raw data, and several came from
methods or people outside the work. They are listed because the habit mattered
more than any single design:

- a "454 SPR-measured non-binder" set we were calibrating against **does not
  exist** — 29 rows, and its README forbids using them as negatives
- `abguard` does **not** check novelty without `--reference-db`
- `stab-ddg` scores **binding** ΔΔG, not thermostability
- the IMGT CDRH3 is `SR` + our design window + `Y`, not the window itself
- two of four epitope hotspots were misnamed in our prose from the first commit:
  **PHE 573** not "Leu573", **ILE 591** not "Val591". The organisers' spec names
  no residues and was always correct, so the designs were aimed at real
  coordinates — only our description was wrong
- **"the gate makes conservative designs impossible" — disproved** by the
  enumerated set
- our ΔΔG scan is **blind to G103's backbone-mediated contact**
- we attributed a 3.42 Å contact to F104; that is the **Phe573↔Tyr105** stack,
  and F104's nearest HER2 atom is 4.27 Å
- an earlier epitope table reported per-residue distances that had been written
  rather than computed; the computed values are W99 4.05 Å, F104 4.27 Å,
  Y105 3.32 Å

---

Both HER2 tables are from **Chinery et al., doi:10.1101/2024.03.26.586756** —
cite them when using this data. Liabilities per the **LAP Antibody Liability
Reference, PLOS Comput Biol 2024;20(3):e1011881**. Structure: **PDB 1N8Z**
(Cho et al., Nature 2003, doi:10.1038/nature01392).
