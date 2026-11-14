# De novo VHH track — HER2 domain IV / trastuzumab epitope

Run date: 14 Nov 2026. Scripts: `pipeline/denovo/`. Nothing here is committed.

## Headline

**This route does not work, and we know why.** BoltzGen `nanobody-anything` ignores
the epitope specification. Under the one scoring method that passes its own positive
control — Boltz-2 unguided co-folding with MSA, which reproduces the trastuzumab
crystal epitope **19/19 residues exactly** — all 6 top candidates dock to a
*different* patch of HER2 domain IV (residues ~495–563), contacting **0 of 4
hotspots and 0–3 of the 19 epitope residues**. The epitope axis fails, and with it
the whole rationale for the track.

Novelty is **not** the blocker (44/85 pass abguard's de-novo gate) and
thermostability is **not** the blocker (median 69.7 °C, max 76.1 °C after one fix).
Epitope targeting is.

### The decisive experiment

| | HER2 target fold | epitope recall | 4 hotspots | Boltz-2 ipTM |
|---|---|---|---|---|
| **trastuzumab Fv (control)** | 0.77 Å | **19/19 = 1.00** | **4/4** | 0.975 |
| split_full_s6_r19 | 1.20 Å | 3/19 = 0.16 | 0/4 | 0.916 |
| split_core_s8_r13 | 1.26 Å | 1/19 = 0.05 | 0/4 | 0.890 |
| cont_core_s5_r14 | 0.99 Å | 1/19 = 0.05 | 0/4 | 0.883 |
| split_full_s6_r18 | 1.12 Å | 0/19 | 0/4 | 0.768 |
| cont_full_s1_r15 | 1.03 Å | 0/19 | 0/4 | 0.759 |
| split_full_s6_r4 | 1.08 Å | 0/19 | 0/4 | 0.475 |

The target folds to ~1 Å in every case, so this is not a bad-target artefact. The
designs bind HER2 confidently (ipTM 0.76–0.92, well over the 0.5 pass line) — at the
wrong site. Every one lands on the same N-terminal patch of domain IV
(≈495–563) regardless of which epitope residues were specified, which is the
signature of `binding_types: binding` acting as weak conditioning rather than a
constraint.

The epitope contacts seen in BoltzGen's *own* refolded complexes (recall up to 0.79)
do not survive an independent fold. They were artefacts of the design-stage geometry.

### The novelty scare, and why it is a false alarm

Measured against 519 real nanobodies (NbThermo), all 85 unique designs are
**74.8–87.0% identical** (median 82.1%), none under 70%. BoltzGen designs CDRs onto
one fixed consensus VHH framework — every design starts
`EVQLVESGGGLVQAGGSLRLSCAASGFTFSSY…` — and the framework is ~68% of the sequence,
so whole-chain identity is pinned above 70% whatever the CDRs do.

**But a whole-chain reading of the rule is unsatisfiable by anything.** For 120
randomly chosen real, published, independently derived NbThermo nanobodies, the
maximum identity to any *other* real nanobody is **0.730–1.000 (median 0.983,
minimum 0.730)**. **0 of 120 real nanobodies would pass** a whole-chain
"<70% identity to any known antibody" test. Our designs at 0.82 are *less*
nanobody-like than a typical real nanobody is to its nearest neighbour.

So the rule cannot mean whole-chain identity against a nanobody database. The only
operational definition we have is `abguard`, the organisers' own tool: identity to
**trastuzumab** plus a CDRH3 Hamming check against 610k reference CDRH3s. Under
that gate **44/85 designs classify `de_novo` (identity < 0.70) and 75/85 pass full
validation.**

Metric calibration: my Biopython global-alignment identity agrees with abguard's
`identity_to_trastuzumab` to within 0.03 (mine ~0.015 higher), so both figures sit
on the same scale as abguard's 0.70 threshold.

**Residual risk, unverified:** we do not know what Adaptyv's database check actually
compares — whole chain, CDRs only, or a window. Ask the organisers. If it is
anything close to whole-chain against a nanobody DB, no VHH entry from any team can
pass, so it almost certainly is not.

## The epitope, measured not assumed

HER2 residues within 5 Å of the trastuzumab Fab in `inputs/12-her2/her2_fv_1n8z.cif`
(mature numbering):

    chain C: 557 558 560 561 569 570 571 572 573 579
    chain D: 591 593 596 598 600 602 603 604 605

The pipeline's four hotspot atoms resolve to **Pro572, Phe573, Ile591, Lys593** —
not "Leu573 / Val591" as the task description says. The atom names are right
(Phe has CD1, Ile has CG2); the residue names in the task text are wrong. All four
are genuine epitope contacts, so the existing hotspot set is sound.

## What was run

| route | jobs | designs |
|---|---|---|
| BoltzGen `nanobody-anything`, HER2 489-607 continuous, full/core epitope | 11 | 220 |
| same, split chains C(489-580)+D(591-607) as in the reference cif | 5 | 100 |
| HER2 domain III+IV (320-607) larger context | 3 | 60 |
| lysozyme positive-control target | 2 | 40 |

840 designs requested across 21 jobs; wave 1 (8 jobs) fully analysed = 160 reported
designs, **88 unique sequences**. Total JapanFold spend **$2.61** of $1,100.

`seed` does not vary BoltzGen output on this API — jobs with seeds 1/2/3 and 17/18
returned byte-identical metrics. Vary the *spec*, not the seed.

**Unresolved / still running when this was written:** the two lysozyme
positive-control design jobs (would have told us whether BoltzGen's low internal
ipTM is HER2-specific or general) and the three HER2 domain III+IV larger-context
jobs (the one remaining chance that a bigger target context fixes the epitope
targeting). Both are now secondary: the Boltz-2 test shows the designs bind HER2
confidently at ipTM 0.76–0.92, just at the wrong site, so BoltzGen's internal ipTM
was misleading in *both* directions and the larger-context arm would need the same
Boltz-2 epitope check before it meant anything. Job IDs are in `bg_jobs_w3.json`
if anyone wants to pick them up.

## Three negative results worth keeping

1. **0/320 designs clear ipTM > 0.5** on BoltzGen's own refold
   (median 0.17, max 0.41; min design-to-target PAE median 15 Å).
2. **Interface confidence is an artefact of target distortion.** Only 42/160 designs
   refolded HER2 within 3 Å of the crystal structure. corr(target CA-RMSD, ipTM) =
   **+0.53** — the *more* the target is deformed, the *higher* the reported interface
   confidence. Designs with correct target geometry score worse (mean ipTM 0.142)
   than badly distorted ones (0.206).
3. **AF2-IG cannot score this target.** The real 1N8Z crystal interface
   (HER2 + trastuzumab VH) scores ipTM 0.147 / interface PAE 23.8 — indistinguishable
   from the median design (0.134 / 25.2). Confounded (VH alone is not a binder
   without VL), but it means AF2-IG as configured gives no discrimination here.
   The one design above 0.5 (`split_full_s6_r11`, 0.602) has epitope recall 0.00 —
   it binds the wrong face.

4. **The epitope specification does not steer the design** — see the headline table.
   This is the finding that closes the route.

## The misleading positive: OpenDDE-abAg with epitope guidance

Antibody-antigen checkpoint, `use_tfg_guidance: true`, epitope pinned to
C571/572/573 + D591/593/596. Read the HER2-C ↔ antibody-heavy **chain-pair** ipTM,
not the global ipTM (which the H↔L pair inflates for the control).

| | global ipTM | HER2-C ↔ heavy pair ipTM |
|---|---|---|
| **trastuzumab Fv (control)** | 0.655 | **0.277** |
| split_full_s6_r18 | 0.496 | **0.481** |
| split_full_s6_r19 | 0.416 | 0.272 |
| cont_full_s1_r15 | 0.401 | 0.257 |
| split_core_s8_r13 | 0.419 | 0.227 |

Read naively this says the designs dock to the trastuzumab epitope as confidently as
trastuzumab itself. **It does not.** The guidance *told* the model where to put the
VHH, so this measures only that the epitope is geometrically compatible with a VHH
paratope — which it is. Boltz-2, asked the same question without guidance, puts
every design somewhere else. Treat TFG-guided ipTM as a compatibility check, never
as evidence of epitope preservation.

## Thermostability: the disulfide finding

Every one of the 160 designs came back **cysteine-free** — BoltzGen removed the
canonical VHH C22–C92 framework disulfide. Restoring it (2 framework substitutions,
no CDR contact) gains a NanoMelt **median +7.0 °C** (range +3.5 to +12.1,
85/85 improved): median 60.6 → 69.7 °C, max 76.1 °C, 20 designs above 70 °C.
Cost: +1.6% identity to trastuzumab, which pushes some designs over the 0.70
abguard line.

NanoMelt install validated against 40 NbThermo nanobodies with measured Tm:
MAE 1.59 °C, Pearson r 0.916. **Caveat: NbThermo is literature-derived and these
40 are almost certainly inside NanoMelt's own 640-measurement training set, so this
validates the installation, not generalisation.** The paper's held-out numbers are
what to quote. Licence **CC BY-NC-SA 4.0 (non-commercial) — must be disclosed.**

## Gates

| check | result |
|---|---|
| abguard schema + liabilities + CDRH3 novelty | 75/85 pass (with disulfide) |
| unpaired cysteines | none (0 Cys before, canonical pair after) |
| ≤250 aa | all pass (115–130 aa) |
| abguard `identity_to_trastuzumab` < 0.70 (the operative gate) | **44/85 pass** |
| whole-chain identity to any known nanobody < 0.70 | 0/85 — but 0/120 *real* nanobodies pass either, so this is not the rule |

## Best candidate, if a slot is spent anyway

`cont_full_s1_r15` + restored disulfide (125 aa, NanoMelt Tm 69.3 °C):

    EVQLVESGGGLVQAGGSLRLSCAASGFTFSSYAMGWVRQAPGKGRELVAAISSDGSSTYYADSVKGRFTISRDNAKNTLYLQMNSLKPEDTAVYYCARLRYSSSYGADLSDFDYWGQGTQVTVSS

Passes every abguard gate: schema, liabilities, ≤250 aa, one canonical disulfide,
`identity_to_trastuzumab` 0.688 → class `de_novo`, CDRH3 `ARLRYSSSYGADLSDFDY`
13 mutations from the nearest of 610k known CDRH3s. In BoltzGen's own refold it
reproduces the target epitope at 0.94 Å and contacts 15/19 epitope residues with
956 Å² buried and 5 H-bonds — the best in the set on those measures.

**But Boltz-2 puts it on HER2 495–537, nowhere near the epitope.** So it is a
plausible VHH that probably binds HER2 domain IV somewhere, is likely to express,
and is likely ~69 °C. Expected outcome: **B-tier (thermostability only), possibly
zero if it does not compete with trastuzumab.** It is not an A-tier play.

Alternate if you want the highest epitope odds: `split_full_s6_r19` — the only
design Boltz-2 places adjacent to the epitope (contacts 560/561/569, grazing the
rim), Tm 73.4 °C, `identity_to_trastuzumab` 0.625 (the widest novelty margin in the
set), Boltz-2 ipTM 0.916.

## What to ask / what to do next

- **The 70% rule**: ask the organisers what it is measured against. 0 of 120 real
  published nanobodies pass a whole-chain reading, so a whole-chain reading cannot
  be the rule — but confirm it rather than assume.
- `seed` is inert on this API — grow a pool by varying the spec (epitope subset,
  binder length, target construct), not the seed.
- **Use Boltz-2 unguided co-folding as the epitope filter** for any future binder
  work on this target. It is the only method here that passed its own positive
  control (19/19 on trastuzumab), it costs ~$0.005 a design, and it would have
  rejected this whole batch in one job before any of the other scoring was done.
  Do not use BoltzGen's internal ipTM, AF2-IG, or TFG-guided OpenDDE for this.
- If the route is retried: the fix is to constrain the epitope *geometrically*, not
  by sequence annotation — e.g. RFdiffusion3 / PXDesign with real hotspot atoms on
  the 1N8Z structure, then graft the resulting paratope onto a VHH scaffold. That is
  a bigger job than one afternoon.
- Untested Tm lever if ever needed: graft designed CDRs onto a high-Tm VHH framework
  (PMID 40613834 reports +10 °C from exactly this move).
