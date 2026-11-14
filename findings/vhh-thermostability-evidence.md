# VHH thermostability as an engineering strategy for HER2 antibodies

**Evidence briefing for the Oslo AI × Antibodies pitch, Team 15. Compiled 2026-11-14
from Amass BiomedCore / TrialCore. Every claim carries its amassId plus PMID/DOI/NCT.**

---

## Bottom line

1. **The scaffold claim is real but narrower than the folklore.** Camelid VHHs do reach
   apparent melting temperatures (60–80 °C, some to 85 °C) that have "never been reported
   for any functional conventional antibody fragment," while scFvs typically sit at 47–69 °C.
   But the *reversibility* story is oversold: in a 70-nanobody survey, the large majority
   aggregate irreversibly on heat denaturation.
2. **NB46 is a strong, specific precedent — and both of your claims about it are confirmed
   in print.** It competitively blocks trastuzumab on HER2, *and* it retains HER2 affinity
   at temperatures up to 90 °C. Source: a Journal of Medicinal Chemistry paper dated
   2026-10-22, i.e. three weeks old.
3. **Every HER2 nanobody in the clinic is a radio-conjugate, not a conventional
   therapeutic.** Three programmes have reached patients (2Rs15d, NM-02/RAD201,
   177Lu-RAD202). No naked anti-HER2 VHH therapeutic is approved or in registrational
   trials. Pitch the scaffold as a *theranostic / stability* play, not as a Herceptin
   replacement.

---

## 1. Are VHH scaffolds intrinsically more thermostable than scFv/Fab?

### The supporting evidence

The foundational quantitative comparison is Dumoulin et al. (Protein Sci 2002), who
characterised six camelid VHH binders across four orthogonal denaturation methods
(guanidinium, urea, temperature, pressure):

- apparent Tm **60–80 °C**, conformational stability ΔG(H₂O) **30–60 kJ/mol**
- chemical unfolding fully reversible, two-state (N ↔ U)
- explicit claim: *"Such high thermodynamic stability has never been reported for any
  functional conventional antibody fragment, even when engineered antigen binders are
  considered."*
- `AMBC_IrgR8teRYDnY9w0cVEfJDgVLvO8` · PMID 11847273 · doi:10.1110/ps.34602 · 626 citations

Functional retention at extreme temperature is independently replicated:

- Llama VHHs bind antigen after a **90 °C** heat shock, where murine mAbs do not. CD/ITC/NMR
  showed most of VHH-R2 is *unfolded* at 80 °C, but adding antigen to the denatured protein
  at 80 °C regenerated the native complex spectrum — **antigen-induced refolding**.
  `AMBC_EcS1EOBwbsTstCUPnC41ZyBaMFC` · PMID 15778955 · doi:10.1002/prot.20378
- Semisynthetic llama sdAb library: binders to vaccinia, lysozyme, cholera toxin, ricin and
  SEB, with "enhanced thermal stability in comparison to conventional IgG **and scFv**
  antibodies"; sdAbs bind antigens at >90 °C and refold after denaturation.
  `AMBC_W8zsa2vharbpP3NewDcRBd1HCIg` · PMID 17165813 · doi:10.1021/ac0610053
- Anti-ochratoxin nanobodies retained antigen binding after **95 °C / 5 min** or
  **90 °C / 75 min**, benchmarked head-to-head against mAb 6H8.
  `AMBC_3RsLuuxayHJx8oKJ779ILPQx4Yn` · PMID 28107910 · doi:10.1016/j.talanta.2016.11.039
- sdAb A3 (anti-SEB): Tm **85 °C** with sub-nanomolar affinity (0.14 nM); CDR-swap
  experiments localise much of the stability to CDR2 and CDR–CDR interactions, not the
  framework. `AMBC_U8akxb8L4OV4CqWTMc9PWCYHkqR` · PMID 24143255 · doi:10.1371/journal.pone.0077678

### The scFv/Fab baseline you are comparing against

- Anti-carbohydrate Se155-4 scFv: Tm **60 °C**, raised to 69 °C only by engineering an
  interchain disulfide. `AMBC_23gmFFXjaELRkWsp0sPiNsPxif2` · PMID 8543036 · doi:10.1016/0014-5793(95)01325-3
- Anti-NP scFvs by DSC: apparent Tm **64–66 °C** (Tyr95-type) and **47–48 °C** (Gly95-type).
  `AMBC_6ofUIwA9OQIhu8ryTO41T6kF2CM` · PMID 28249223 · doi:10.1016/j.molimm.2017.02.010
- Fab vs scFv dissection: mutual stabilisation across VH/VL and CH1/CL interfaces, with the
  interchain disulfide essential; Fabs gain **kinetic** stabilisation (very slow unfolding)
  rather than higher equilibrium stability. Relevant because it means "Fab is more stable"
  and "Fab has a higher Tm" are different statements.
  `AMBC_UbuiiiqWGx5CJL3fTGRstMrcIM8` · PMID 15769469 · doi:10.1016/j.jmb.2005.01.053
- Orthogonal 2025 profiling (nanoDSF/DLS/SEC/SAXS/CD/EM) of IgG vs scFv vs tandem bi-scFv:
  engineered fragments show "increased aggregation propensity and reduced conformational
  stability" versus full-length IgG.
  `AMBC_KgotbB8JuXngENs9CbgBVz2BJBL` · PMID 41079721 · doi:10.3389/fphar.2025.1667210
- Review-level synthesis of the VHH-vs-scFv comparison (higher solubility, stability,
  smaller size, lower production cost), 206 citations.
  `AMBC_YTjXMostTwzQfnjFQx2Dg6mYKzJ` · PMID 34863296 · doi:10.1186/s40364-021-00332-6

### Three caveats you should state before a judge finds them

- **Reversibility is mostly a myth.** Kunz et al. characterised ~70 nanobodies and found
  irreversible aggregation on heat denaturation for the *large majority*, with
  consequences for stability and immunogenicity. An optional second disulfide bond
  suppresses aggregation; a longer CDR loop helps further. Nanobodies with two disulfides,
  however, tend to lose function in the cytosol.
  `AMBC_5B9K3jhVhCjmovlNTqCmSdETHHG` · PMID 29784954 · doi:10.1038/s41598-018-26338-z · 195 citations
- **Thermostability does not buy you affinity.** Across three nanobodies analysed from
  18–100 °C, thermodynamic stability showed *no correlation* with target binding affinity.
  Do not let the pitch imply stability is a proxy for potency.
  `AMBC_OrGlPaR0YBpRriW4bt5Fbajx7jt` · PMID 39593207 · doi:10.1002/1873-3468.15064
- **Stability is not predictable from framework sequence.** NbThermo (564 curated Tm values)
  reports no apparent framework sequence pattern separating low- from high-Tm nanobodies;
  the highly variable loops dominate. Llama and camel Tm distributions are similar.
  `AMBC_OR8D1qsp2aJpC6kT6NBNMU9klfw` · PMID 37042467 · doi:10.1093/database/baad021

### Directly useful to our pipeline

If we want a *computed* stability axis (recall `stab-ddg` is a binding-ddG predictor and
leaves thermostability unaddressed), three trained predictors exist:

| tool | training data | reported accuracy | record |
| --- | --- | --- | --- |
| TEMPRO | 567 seqs (incl. NbThermo) | MAE 4.03 °C, RMSE 5.66 °C | `AMBC_8yFWQQblsRUsYJF6QYZV4BY3104` · PMID 39154093 |
| NanoMelt | 640 Tm (129 new) | ensemble over seq embeddings, state-of-the-art | `AMBC_4p5hpBu388L12EXJjuWsRFDy0rL` · PMID 39772905 |
| NBsTem | 514 experimental Tm + 704 MD structures | Pearson 0.83 external test | `AMBC_HhzUY6m5QUypsvb4H36IYwhlitK` · PMID 41564239 |

Caveat: all three are nanobody-specific. None transfers to a CDRH3-redesigned IgG scaffold
like ours, so quoting them as *our* stability number would be wrong.

---

## 2. NB46 — the precedent, in depth

**Both parts of the claim are confirmed.** NB46 appears in exactly two publications, both
in J Med Chem, both from Shaoli Song's group (Shanghai).

### 2a. The thermostability + competition paper — the one that matters

> "NB46 retained HER2 affinity at temperatures up to **90 °C**, demonstrating exceptional
> thermostability. … **Blocking studies confirmed that NB46 and trastuzumab competitively
> bind to HER2.**"

- `AMBC_QrqakQNKusdXCwo4XxWBsdgLhLU` · **PMID 42911746** · **doi:10.1021/acs.jmedchem.6c02664**
- *Evaluation of 68Ga Labeling Strategies and Thermostability of HER2-Targeted Nanobody NB46*
- Zhao L, Bai L, Zhu J, Sun L, Liu X, Lu X, Zhang J, Xu X, Song S
- J Med Chem, **2026-10-22**. JuFo 3. 0 citations (too new).

What the paper actually did, and the limits of it:

- Compared **random vs site-specific ⁶⁸Ga labeling**. Both conjugates kept HER2 affinity
  comparable to unmodified NB46; no significant difference in cellular uptake, PK,
  micro-PET/CT or biodistribution in HER2+ and HER2− tumour-bearing mice.
- **>95% radiochemical yields at all tested temperatures.** The 90 °C figure comes from this
  radiolabeling-robustness series — it is *functional retention after thermal exposure
  during conjugation*, *not* a measured Tm. There is no Tm value for NB46 in the record.
- Elevated labeling temperature did not compromise tumour targeting but gave a mild increase
  in renal uptake, most notably the random conjugate at 90 °C.
- The trastuzumab competition is established by **blocking studies**, not by structure or
  epitope mapping. The paper does not state which HER2 domain NB46 binds.
- Amass holds abstract only (`hasFulltext: false`), so exact KD values and the thermal
  protocol are not in the record. **If the pitch needs a number beyond "90 °C", someone has
  to pull the PDF.**

### 2b. The parent radiotheranostics paper

- `AMBC_zvxHn9OjB8s2xX37uG5EkaDVMc` · **PMID 41963252** · **doi:10.1021/acs.jmedchem.6c00837**
- *Development and Preclinical Evaluation of a Novel 68Ga/177Lu-Labeled Nanobody for
  Radiotheranostics of HER2-Positive Tumors* · J Med Chem, 2026-04-10. 2 citations.
- NB46 engineered with a **C-terminal Gly-Ser-Cys (GSC) tripeptide** for site-specific NOTA
  conjugation, radiolabeled with ⁶⁸Ga and ¹⁷⁷Lu.
- In HER2+ xenografts: rapid blood clearance, high tumour accumulation, low non-target
  uptake; [¹⁷⁷Lu]Lu-NOTA-NB46 showed prolonged tumour retention.
- **Single-dose [¹⁷⁷Lu]Lu-NOTA-NB46 inhibited tumour growth and extended survival
  dose-dependently, with good tolerability. Kidneys were the dose-limiting organ.**
- Status: **preclinical, mouse only.** The paper's own conclusion is that it "supports
  further translation into clinical investigation" — i.e. NB46 has not entered a trial.
  No TrialCore cross-links on either record (`referencesTrialCore: []`).

### 2c. Why NB46 is the right precedent for *our* design problem

We redesigned trastuzumab's CDRH3 against HER2 domain IV. The relevant question is whether a
single-domain scaffold can hold the trastuzumab epitope at all — and the literature says the
answer depends entirely on the clone:

| VHH / sdAb | relation to trastuzumab epitope | affinity | record |
| --- | --- | --- | --- |
| **NB46** | **competes** (blocking) | not stated in record | PMID 42911746 |
| Nikkhoi VHH | **binds the same epitope** as trastuzumab | sub-nanomolar | `AMBC_B8fWxL8oMPmGiT52sD3PuMX7Ote` · PMID 37632647 |
| 2Rs15d (clinical) | does **not** compete with trastuzumab or pertuzumab | low nM | `AMBC_ZOTqq4bXeYRnjgVlalNVWThvwql` · PMID 21478264 |
| NM-02 / RAD201 (clinical) | **non-overlapping** epitope (X-ray) | — | `AMBC_U1ejblWQqpAcxfe3HF7dy4EjrZ4` · PMID 37459326 |
| Hussack panel (6 VHHs) | all three epitopes **distinct** from trastuzumab | 1–51 nM | `AMBC_CbPmbuFS7r9HzqdPaiTwMZJPbRl` · PMID 30518413 |

So: the trastuzumab-competitive VHH epitope is **precedented but unusual** — most isolated
anti-HER2 VHHs land elsewhere on the ectodomain. NB46 and the Nikkhoi clone are the two
cases that put a single domain on trastuzumab's site, and NB46 is the only one of those with
a published 90 °C result. That pairing is the pitch.

---

## 3. Clinical and preclinical HER2 nanobody programmes

### In patients

**2Rs15d (Brussels / VUB — the most mature)**
- Preclinical selection from a 38-nanobody panel: `AMBC_ZOTqq4bXeYRnjgVlalNVWThvwql` · PMID 21478264 · doi:10.1096/fj.10-180331
- ⁶⁸Ga-NOTA-2Rs15d tracer development: `AMBC_EIWU76ylBLlELqQSQezLFRRVyrQ` · PMID 23487015 · doi:10.2967/jnumed.112.111021
- **Phase I**, 20 women, primary/metastatic breast carcinoma: no adverse reactions, effective
  dose 0.043 mSv/MBq (~4.6 mSv/patient), 90 min imaging advised.
  `AMBC_GZSCaQebQwK6sPV4h7KeXjlVcLP` · PMID 26449837 · doi:10.2967/jnumed.115.162024 · **435 citations**
- **Phase II**, 20 patients: repeatability coefficient **21.8%** (comparable to ¹⁸F-FDG); safe;
  uptake even in HER2-low/negative patients; revealed interlesional heterogeneity that FDG
  missed. `AMBC_M1PUUI51bwAnijRYt0ORTecyQyp` · PMID 38302159 · doi:10.2967/jnumed.123.266254

**NM-02 → RAD201 (Shanghai General Hospital / Nanomab)**
- Phase I, 10 breast cancer patients, ⁹⁹ᵐTc-NM-02 SPECT/CT: no drug-related adverse
  reactions, effective dose 6.56×10⁻³ mSv/MBq, SUVmax 1.5 cutoff for HER2 positivity.
  `AMBC_YS5Nb9FBRwms9kn9enWqj9Fx3iI` · PMID 34328338 · **NCT04040686** (`AMTC_PFacFYy3JdTrbt3i5L3ASxyIdmU`, completed, n=10)
- 30-patient imaging trial: uptake correlated with HER2 expression in newly diagnosed
  patients (not in treated patients); no non-specific uptake in inflamed tissue; revealed
  intra-/intertumoral heterogeneity in 9/30.
  `AMBC_AzomIxTaqb1z2Kk0enwqYA1uIV4` · PMID 38459598 · **NCT04674722** (`AMTC_CPokJ8NVSw2MoNDv6tZhNe2T8Ms`, EARLY_PHASE1, completed 2024-10-30, n=43, also ¹⁸⁸Re therapy arm)
- Renamed RAD201; German series of 6 heavily pretreated patients, kidney dose cut ~50% with
  Gelofusine. `AMBC_BfmJWluZ3KxLUeqbPcdBHYhxHIX` · PMID 36482076 · doi:10.1007/s00259-022-06066-3
- Humanisation + X-ray structure + developability profiling — the cleanest published
  "VHH-as-drug-candidate" dossier for HER2, and it reports favourable thermal and physical
  drug-like properties. `AMBC_U1ejblWQqpAcxfe3HF7dy4EjrZ4` · PMID 37459326 · doi:10.1371/journal.pone.0288259

**¹⁷⁷Lu-RAD202 — HEAT trial (Radiopharm Theranostics). The only HER2 sdAb *therapeutic* in humans.**
- **NCT06824155** (`AMTC_3jVE9rXfVW42lAvHI0a21KCI41p`): first-in-human Phase 0/1, open-label,
  **RECRUITING**, started 2025-02-12, n=30 estimated, Australia. Also registered as
  ACTRN12625000191493 (`AMTC_RiV2SkFM9lJToxeEPjuTRo1FQHz`, n=24 actual).
- AACR 2026 late-breaker: as of Nov 2025, **3 participants treated at 30 mCi**; lesion
  dosimetry up to **3.6 Gy/lesion**; no serious treatment-related AEs (two Grade 1:
  dysgeusia, pleural effusion); trial now enrolling at 75 mCi.
  `AMBC_NYZwEijJNAiiyR6m82WlxtGS1pj` · doi:10.1158/1538-7445.am2026-ct046 (no PMID — conference abstract)

### Preclinical, worth one line each

- **NB46** ⁶⁸Ga/¹⁷⁷Lu radiotheranostic pair — see §2b. PMID 41963252.
- Click-chemistry "click-to-label" HER2VHH platform: one nanobody, four isotopes
  (¹⁸F/⁶⁸Ga/⁹⁹ᵐTc/¹⁷⁷Lu) via TCO–tetrazine, avoiding harsh labeling conditions; ¹⁷⁷Lu arm
  inhibited growth dose-dependently. `AMBC_4ZPqWvtCQp6oR3E9ZlfifuHHMOR` · PMID 42419170 · doi:10.1016/j.ejmech.2026.119104
- Biparatopic nanobody binder A9F5-H2F5-Fc targeting HER2 ECD I+II: outperformed
  trastuzumab+pertuzumab on internalisation/degradation and synergised with trastuzumab in
  resistant cells. `AMBC_TCt4e6wWzZGZwA4TXHEEdmH5u3S` · PMID 41221272 · doi:10.3389/fimmu.2025.1711448
- Four affinity-matured anti-HER2 nanobodies, three epitopes; A9B5-Fc (domain I) beat
  pertuzumab as a trastuzumab partner under resistance.
  `AMBC_Lf1tLFJbBEbZZt2jYxtvARPgyXg` · PMID 37954614 · doi:10.3389/fimmu.2023.1292839
- 2Rs15d nanobody–drug conjugate (MMAE, cathepsin-B-sensitive), site-selective.
  `AMBC_8DK9vqkM1RIpy5QUrchKoxMAos4` · PMID 39860260 · doi:10.3390/molecules30020391
- Nanobody-based CAR-T is an established format across >10 targets, with anti-BCMA
  nanobody-CAR-T clinically comparable to scFv-CAR-T.
  `AMBC_TZKfaNVNE7NlP3OMOXjwEpsktgG` · PMID 35468841 · doi:10.1186/s40364-022-00371-7

### What is *not* there

No approved anti-HER2 nanobody therapeutic, and no Phase 2/3 HER2 nanobody therapeutic
trial in TrialCore. DrugCore has no harmonised anti-HER2 VHH entry (searching it returns
unrelated protein drugs) — so "22K drugs and none of them is this" is itself a finding if
we want to argue whitespace. The HER2-targeted single-domain field in the clinic is
**imaging-first, radiotherapy-second**.

---

## 4. What this licenses us to say at 17:00 — and what it does not

**Safe to claim:**
- VHH scaffolds reach apparent Tm 60–80 °C where comparable scFvs sit at 47–69 °C, a gap
  documented since 2002 and never matched by a conventional antibody fragment. (PMID 11847273)
- A HER2-targeting nanobody, NB46, retains HER2 affinity after exposure to 90 °C and
  competitively blocks trastuzumab. (PMID 42911746, published 2026-10-22)
- Three anti-HER2 single-domain programmes have reached patients; one, ¹⁷⁷Lu-RAD202, is an
  actively recruiting first-in-human therapeutic. (NCT06824155)

**Do not claim:**
- *"Nanobodies refold reversibly."* Most do not. (PMID 29784954)
- *"NB46 has a Tm of 90 °C."* The record shows functional retention up to 90 °C in a
  radiolabeling series, not a measured Tm.
- *"NB46 binds domain IV"* or *"NB46 shares trastuzumab's epitope."* It *competes*; the
  epitope is not mapped in the record. Most anti-HER2 VHHs bind elsewhere.
- *"A thermostable scaffold will be a better binder."* Stability and affinity were explicitly
  uncorrelated. (PMID 39593207)
- Any transfer of nanobody Tm predictors to our IgG-scaffold designs.

**Honest framing of the gap in our own work:** nothing here addresses the thermostability of
*our* CDRH3-redesigned trastuzumab variants. The VHH literature is an argument for a
*different scaffold choice*, not evidence about the designs we generated. If the pitch uses
this section, it should be framed as "the stability axis we could not score, and the scaffold
the field uses to solve it."

---

## 5. Provenance

All records retrieved from Amass BiomedCore and TrialCore on 2026-11-14. Searches that
returned the key evidence: `NB46`; `nanobody retains antigen binding after heating to 90
degrees reversible thermal unfolding refolding`; `nanobody VHH thermostability melting
temperature compared with scFv and Fab fragments`; `anti-HER2 nanobody single domain
antibody epitope competition with trastuzumab affinity`; `99mTc-NM-02 ... first-in-human`;
`68Ga anti-HER2 nanobody 2Rs15d clinical PET imaging patients`; TrialCore `NCT06824155`.

Two records are conference abstracts with no PMID (AACR 2026 CT046; ASCO 2015 e11600) and
are cited by DOI. Neither NB46 paper has fulltext in Amass; both are abstract-level only.
