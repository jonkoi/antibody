# The 14 shortlisted designs — authoritative record

All numbers verified against `runs/team15-claude-n100-seed1/liability_report.json`
(produced by `abguard` on the VM) and the StaB-ddG and FACS-model outputs. Sorted
by the FACS model's mean score.

| design | IMGT CDRH3 | novel | min dist | liabilities | high-sev | CDR3 Trp | FACS mean |
| --- | --- | :-: | ---: | ---: | ---: | :-: | ---: |
| her2_cdrh3_19 | SRSTASVGYAMDY | ✓ | 4 | 6 | 2 | — | **0.998** |
| her2_cdrh3_85 | SRSTAAVGYAMDY | ✓ | 4 | 6 | 2 | — | **0.973** |
| her2_cdrh3_96 | SRSTSSVGYAMDY | ✓ | 3 | 7 | 2 | — | 0.904 |
| her2_cdrh3_94 | SRWTAAVDYAHDY | ✓ | 3 | 6 | 2 | **yes** | 0.863 |
| her2_cdrh3_61 | SRSLGSADYAFDY | ✓ | 3 | **5** | 2 | — | 0.840 |
| her2_cdrh3_7 | SRSTASVDYAMDY | ✓ | 3 | 6 | 2 | — | 0.819 |
| her2_cdrh3_71 | SRSTSTYFYAFDY | ✓ | 3 | 6 | 2 | — | 0.572 |
| her2_cdrh3_60 | SRWSYAYDYAMDY | ✓ | 3 | 7 | 2 | **yes** | 0.133 |
| her2_cdrh3_6 | SRSTSAYGYAMDY | ✓ | 3 | 7 | 2 | — | 0.108 |
| her2_cdrh3_75 | SRSTFSVDYAMDY | ✓ | 4 | 6 | 2 | — | 0.092 |
| her2_cdrh3_88 | SRSLSLGNYAFDY | ✓ | 3 | **5** | 2 | — | 0.069 |
| her2_cdrh3_15 | SRSTSAVDYAMDY | ✓ | 3 | 7 | 2 | — | 0.062 |
| her2_cdrh3_4 | SRYWSLHYYAMDY | ✓ | 3 | 7 | 2 | **yes** | 0.011 |
| her2_cdrh3_77 | SRSTAAVDYAMDY | ✓ | 3 | 6 | 2 | — | *unusable* |

**`team-validate` verdict: `shortlist: 14/14 design(s) valid`**, with the
525,063-CDRH3 novelty gate active via `--reference-db`. `--no-novelty` was not
used. All 14 are class `novel`, identity to trastuzumab 97.4–97.8%, minimum
distance to any known CDRH3 **3 or 4 mutations** — so every one clears the
≥3 floor, several with nothing to spare.

`_77` is listed but carries **no score**: it moved 0.229 → 0.790 between training
runs on nothing but a reshuffled training order. Printing a number for it would
imply a precision that does not exist.

---

## Developability: our designs are cleaner than the approved drug

87 liability occurrences across 14 designs, **0 hard errors**, and the headline:

**Wild-type trastuzumab scores 8 liabilities. Every one of our designs scores
5–7.** No design introduces a *new* high-severity liability — all 2 high-severity
hits per design are the `Isom` (`DT`, CDRH1) and `DeAmdH` (`NG`, CDRH2) motifs
inherited from trastuzumab's own framework, identical in all 14 and in wild-type.

Reference: **LAP Antibody Liability Reference, PLOS Comput Biol 2024;20(3):e1011881.**

Liability classes present — none in the designed window except as noted:

| tag | severity | occurrences | designs |
| --- | --- | ---: | ---: |
| TrpOx | medium | 17 | 14 |
| Isom | **high** | 14 | 14 |
| DeAmdH | **high** | 14 | 14 |
| DeAmdM | medium | 14 | 14 |
| DeAmdL | low | 14 | 14 |
| MetOx | medium | 10 | 10 |
| FragM | medium | 4 | 4 |

Absent entirely: FragH, mCys, xCys, Ngly, Hydro, IntBind.

Note a tool quirk: `TrpOx` fires on the conserved FR4 tryptophan of `WGQG` at
IMGT 118 in all 14 and in wild-type, although the reference scopes TrpOx to CDRs.
Ignore that one; it is not a design property.

---

## The W99 question is now answered four ways, and they agree

Keeping the conserved tryptophan looked like the right call on the SPR
correlation alone. Four independent lines now say it is not the discriminator:

1. **Binding energetics (StaB-ddG).** W99→Ser costs **0.207 ± 0.148 kcal/mol** —
   zero within error. **Tyr105 is the real hotspot at 0.866**, and *all 14
   designs already preserve it*. Only 1 of 14 preserves F104 (0.462), the second.
2. **Measured-binding model (FACS, test AUC 0.9842).** Both stable-high designs
   lost W99; `_60`, a W99-retainer, scores 0.133 — near the bottom. `_60` also
   destroys F104, which may be precisely why.
3. **The SPR correlation is substantially confounded.** W99-keeping variants
   carry median **4** substitutions versus **6** for W99-losing ones. *But it does
   not fully vanish under stratification — at 6 substitutions it is still 13.1 nM
   vs 147.0 nM (n=10 vs 19), an 11× gap.* This residual is the one live thread.
4. **Developability.** The three designs with a tryptophan in CDRH3 — `_60`,
   `_94`, `_4` — are the only ones carrying a genuine `TrpOx` liability. A CDR3
   Trp is an oxidation and aggregation hotspot, unlike the FR4 Trp everyone has.
   **Retaining W99 carries a developability cost and buys no measurable affinity.**

**Recommended reading:** W99 is not decisive. Lead with `_19` and `_85`.

Every single-position effect in the window is ≤0.87 kcal/mol — at most ~4.3× in
Kd, **inside the 6.2× assay noise floor**. No single-residue choice here can
produce a measurable affinity change, which is independent structural
confirmation that CDRH3-tip redesign is the wrong lever for affinity.

---

## The two measured-data signals DISAGREE — and that decides the pick

StaB-ddG scored all 14 designs as multi-site variants (the API scores them
jointly; we did not sum anything). Every design is predicted to bind *worse* than
wild-type, from 2.0x to 19.0x:

| design | ddG kcal/mol | ±SE | fold worse than WT | FACS mean |
| --- | ---: | ---: | ---: | ---: |
| **her2_cdrh3_61** | **0.42** | 0.23 | **2.0x** | 0.840 |
| her2_cdrh3_88 | 0.64 | 0.23 | 3.0x | 0.069 |
| her2_cdrh3_60 | 0.72 | 0.21 | 3.4x | 0.133 |
| her2_cdrh3_4 | 1.24 | 0.15 | 8.2x | 0.011 |
| her2_cdrh3_6 | 1.27 | 0.15 | 8.5x | 0.108 |
| her2_cdrh3_94 | 1.31 | 0.23 | 9.1x | 0.863 |
| her2_cdrh3_71 | 1.35 | 0.19 | 9.7x | 0.572 |
| her2_cdrh3_75 | 1.38 | 0.25 | 10.3x | 0.092 |
| her2_cdrh3_77 | 1.48 | 0.23 | 12.1x | *unusable* |
| her2_cdrh3_15 | 1.50 | 0.22 | 12.5x | 0.062 |
| her2_cdrh3_7 | 1.53 | 0.24 | 13.1x | 0.819 |
| her2_cdrh3_85 | 1.70 | 0.17 | 17.7x | 0.973 |
| her2_cdrh3_19 | 1.72 | 0.18 | 18.2x | 0.998 |
| her2_cdrh3_96 | 1.75 | 0.18 | 19.0x | 0.904 |

**Spearman(ddG, FACS) = +0.522 (n=13).** Higher ddG means worse binding, so a
*positive* correlation means **the designs the FACS model likes best are the ones
StaB-ddG says bind worst.** Our two independent measured-data-grounded signals
actively disagree.

Three readings, and we cannot currently distinguish them:
1. The FACS model is motif-matching the `ST..V[DG]YAMD` family that dominates the
   accepted pool, and that motif is energetically expensive.
2. StaB-ddG's multi-mutant extrapolation is unreliable — **unverified** for
   5-9 simultaneous substitutions, which is exactly our regime.
3. They measure different things: FACS is cell-surface binding (avidity,
   expression included), ddG is interface binding energy on a static structure.

**The 6.2x assay noise floor equals ddG 1.08 kcal/mol.** Only three designs are
predicted to land within it: **`_61` (2.0x), `_88` (3.0x), `_60` (3.4x)**.
Everything the FACS model ranked top is predicted 18-19x worse — far outside.

### `_61` is the pick, because it is the only design strong on every signal

| signal | `her2_cdrh3_61` (`SLGSADYAFD`) |
| --- | --- |
| StaB-ddG | **best of 14** (0.42, within the noise floor) |
| FACS model | 0.840 — high tier |
| Liabilities | **5, tied cleanest**, none in CDRH3 |
| CDR3 Trp | none |
| Novelty | `novel`, distance 3, 14/14 valid |

No other design is good on all four. `_19`/`_85` win on FACS alone and are the
worst on energetics; `_88`/`_60` are good on energetics and bottom-tier on FACS.

**Recommended submission: `_61` first.** Then choose between `_19`/`_85` (trust
the measured FACS model) and `_88` (trust the energetics) — and say in the pitch
that we are deliberately hedging across two signals that disagree, rather than
pretending one is authoritative.

---

## If you want a risk-based cut

- **Submit first:** `_19`, `_85` — stable-high on measured data, 6 liabilities, no
  CDR3 Trp, and both at minimum distance **4** rather than 3.
- **Cleanest on developability:** `_61` and `_88` — 5 liabilities each, zero
  CDRH3-specific. `_61` also scores 0.840. **`_61` is the best combined pick.**
- **Cut on risk if trimming:** `_60`, `_94`, `_4` (CDR3 Trp) — and `_60` and `_4`
  are also bottom-tier on the FACS model.
- **Do not submit with a number:** `_77`.
