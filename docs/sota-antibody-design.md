# SOTA deep-learning antibody design for HER2-like antigens — survey, Nov 2026

Compiled 14 Nov 2026. Companion to [`../log1.md`](../log1.md) Session 3, which records how
this was gathered and what it changes. Reference material: the log holds the decisions.

**Scope.** "HER2-like antigen" here means a large, well-characterised, structurally solved
cell-surface protein antigen with a known therapeutic epitope and an existing
sub-nanomolar clinical antibody. That framing matters: it is what makes the *optimisation*
problem different from the *de novo* problem, and most method papers only address one.

**Selection rule.** Methods are included only where there is wet-lab outcome data, or where
the method is a scorer whose predictive power has been *measured against* wet-lab data.
In-silico-only benchmarks are excluded — there are hundreds and they do not discriminate.
Where a number comes from a preprint, vendor blog or secondary coverage rather than a
peer-reviewed primary source, it is marked.

---

## Part 1 — The calibration set: four campaigns with measured outcomes

These are the only honest priors on "what hit rate should we expect". Note how much of the
variance is assay definition rather than method quality.

| Campaign | Designs tested | Outcome | Hit definition | Source quality |
| --- | --- | --- | --- | --- |
| **Adaptyv EGFR R1** (2024) | 201 | 146 expressed (73%), **5 binders (2.5%)** | BLI | preprint, peer-review pending |
| **Adaptyv EGFR R2** (2024) | 400 | 378 expressed (95%), **53 binders (13%)** | BLI, 3 replicates | same |
| **Adaptyv Nipah** (2026) | ~1,196 | 1,028 expressed (86%), **111 binders (9.3%)**, 28 ≤ single-digit nM, best 0.37 nM | BLI proxy KD, 1:1 fit | competition site, internally inconsistent (111 vs 116) |
| **Bits to Binders / CD20** (2026) | 11,796 synthesised | 6,811 recovered (58%), **707 CD20-specific (5.9%)**; team rates **0.6%–38.4%** | pooled CAR-T proliferation, FDR 5% | preprint |
| **Chai-2** (2025) | ≤20 per target × 52 targets | **15.5%** (20.0% VHH, 13.7% scFv); ≥1 binder for 26/52 targets | 1-point BLI >300% background | preprint, no availability |

**🔑 Adaptyv Bio runs the first three and will test our designs.** The EGFR competition is
the closest published analogue to this challenge: crowdsourced, scored in-silico, a fixed
number tested, same lab, same assay.

### What the Adaptyv EGFR competition winner actually did
**R2 winner (Cradle, KD 1.21 nM, 8.2× better than the cetuximab control at 9.94 nM) started
from cetuximab, reformatted it as an scFv, and mutated only framework residues — keeping
the CDRs intact.** Best fully de novo binder in the same round: **82 nM** (BindCraft).
Only ~10 of ~58 validated binders across both rounds came from de novo or hallucination
methods, and the hallucinated share of true binders *fell* from 40% (R1) to 12% (R2).

Kruskal-Wallis across design strategies found no significant difference in median KD
(p = 0.374), so this is "the best designs came from optimising an existing lead", not
"optimisation beats de novo on average". The authors' own summary: *"current design and
filtering strategies optimize for binding probability but not affinity."* That is the same
gap `log0` found in the FACS data, reached from the other direction.

### The affinity ceiling of de novo antibody design
Collecting the best first-round, no-experimental-maturation affinities actually reported
for **antibody or VHH formats**:

| Method | Best reported de novo antibody/VHH affinity |
| --- | --- |
| RFantibody (Nature 2025) | tens–hundreds nM initially; single-digit nM **only after OrthoRep affinity maturation** |
| Chai-2 (preprint) | 82 nM (VHH, CCL2); 3.3 µM (scFv) |
| BoltzGen (blog/secondary) | 13 nM (PD-L1), 14.4 nM (IL-7Rα) nanobodies |
| Germinal (Nat Biotech 2026) | "nanomolar", 43–101 designs per antigen |

**Trastuzumab is 0.146 nM (median of 8 SPR replicates in our own data).** No published de
novo antibody or VHH campaign has produced a binder within two orders of magnitude of that
in a single round. The de novo track therefore cannot be won on affinity against a
trastuzumab benchmark; it can only be won on novelty, epitope and thermostability.
*Caveat: the Nipah competition's best overall design was 0.37 nM, but that entry is
classed "Other", not antibody/VHH, and the KD is a standardised 1:1 proxy fit.*

---

## Part 2 — Methods, by what they are for

### A. De novo epitope-targeted antibody / VHH design

| Method | Venue | Format | Wet-lab evidence | Available? | On our stack? |
| --- | --- | --- | --- | --- | --- |
| **RFantibody** (RFdiffusion fine-tune) | **Nature, Nov 2025** (PMID 41193805) | VHH, scFv, full mAb | 4 epitopes, cryo-EM confirmed poses, all 6 CDRs atomically accurate in one scFv | open | RFdiffusion3 / RFD4 on JapanFold + Nebius |
| **BoltzGen** | preprint 2025 | nanobody, protein, peptide | nM for 6/9 novel targets at ≤15 designs/target (lab blog); Nipah competition 121 designs entered, per-method rate unpublished | **open** | ✅ **JapanFold `nanobody-anything`, `antibody-anything`, `protein-redesign`** |
| **Germinal** | **Nat Biotechnology, Jun 2026** (PMID 42337361) | nanobody, scFv | functional binders on 4/4 targets, 43–101 designs each; 4–22% success (preprint) | **open**, `SantiagoMille/germinal`; PyRosetta + IgLM licence encumbered → `OpenGerminal` Apache-2.0 | self-host only |
| **Chai-2** | preprint 2025 | VHH, scFv, miniprotein | 15.5% over 52 targets, ≤20 designs each | **no availability stated anywhere in the paper** | ✗ |
| **tFold System** | **Nat Communications, Dec 2025** (PMID 41365923) | mAb | nanomolar by SPR after phage display, 4 antigens, epitope-specific | open | self-host |
| **NanoDesigner** | J Cheminformatics 2025 | nanobody | ~2× de novo success rate via EM-style iterative docking↔CDR refinement | open | self-host |
| **RFD4-Proteina** (HER2 fine-tune) | organisers' checkpoint | CDRH3 only | organisers' benchmark only | hackathon | ✅ already in use |
| **PXDesign** | preprint 2025 | backbone only, no sequence | — | open | ✅ JapanFold |

**Chai-2 is the best-reported method and the least usable.** It names no API, no weights and
no licence; it also never verified epitope specificity (the paper says competitive binding
and structure determination *"should be performed"*), does not mention HER2, and uses
trastuzumab only as a polyreactivity control. Treat its 15.5% as the field's upper bound,
not as something we can call.

**BoltzGen is the usable one**, and `nanobody-anything` plus `protein-redesign` are already
sitting on our JapanFold key at ~$0.04/job (p90).

### B. Optimising an existing antibody — the optimisation track

| Method | Venue | Antigen | Result | Code |
| --- | --- | --- | --- | --- |
| **DualGPT-AB** | **Nat Computational Science, Apr 2026** (PMID 41986730) | **HER2** | **8 of 100 sampled designs bind HER2 well; wet-lab confirmed enhanced tumoricidal activity vs Herceptin** | none found |
| **AB-Gen** | Genomics Proteomics Bioinformatics 2023 (PMID 37364719) | **HER2** | GPT + RL over CDRH3; 509 of N generated passed all property filters; in-silico only | Zenodo 10.5281/zenodo.7657016 |
| **Cradle (framework-only optimisation)** | Adaptyv EGFR R2 | EGFR | **1.21 nM, competition winner, 8.2× over control** | proprietary |
| **Ng et al. framework stabilisation** | **mAbs 2025** (PMID 40653772) | **HER2 / trastuzumab** | see below | Rosetta `pmut_scan` |
| **MicroMutate pipeline** | Brief Bioinform 2025 (PMID 40889116) | H7N9 HA, DR5 | 4.62× and 2.07× affinity from single points, starting sub-nM | — |
| **NAZ-mAb variants** | Curr Res Struct Biol, Oct 2026 (PMID 42909185) | **p95HER2** | **cautionary: FoldX-ranked variants were all WORSE** (326 nM, 94.5 nM vs 53.3 nM parent) | — |

**DualGPT-AB is the closest published prior art to our optimisation track** — multi-property
conditional CDRH3 generation against HER2, wet-lab validated against Herceptin itself,
published 7 months ago in a Nature journal. No code was findable. It is pitch-relevant
(we are not doing something speculative) and a benchmark we should cite, not a tool.

**The NAZ-mAb paper is the cautionary twin**, and it is from Oslo (Kyte/Greiff). Fixed-backbone
FoldX interaction-energy ranking on a HER2-family antibody produced two variants that
*looked* better by ELISA and were **2–6× worse by Biacore**. The authors frame it explicitly
as "a cautionary benchmark … highlighting the limitations of fixed-backbone
interaction-energy ranking as a predictor of binding affinity and yield." Published a month
ago, same city, same antigen family.

### C. Thermostability

**Fab / paired-Fv format — no longer unscored.** `README.md` says no folding-stability
predictor exists for Fab format in our stack, and `log1` Session 1 found the nanobody Tm
tools are VHH-only. Both remain true. But Ng et al. 2025 did it on **trastuzumab itself**:

- **Method:** Rosetta `relax` on PDB **7mn8** (cryo-EM trastuzumab Fab–HER2, 3.45 Å) to
  generate n = 200 backbones, then `pmut_scan` for ΔΔG_folding at every VH/VL position × 19
  amino acids, then a second scan over double mutants.
- **Validated output: VH `S85N+R87T` raises Tm (DSF) and retains both HER2 binding and
  ADCC.** Two other designs failed informatively: `R50S+R59N` (near the CDRs) lost binding
  and ADCC; `A40R+K43D` had production yield too low to test.
- **⚠️ `S85N+R87T` introduces a novel N-linked glycosylation sequon.** The authors cannot
  exclude that ectopic glycans drive their result. Our own `pipeline/score_developability.py`
  penalises an `N[^P][ST]` sequon at weight 4.0 — the heaviest liability it scores. So the
  one validated stabilising mutation for our exact molecule trips our own worst liability
  flag. That tension is real and must be stated in any pitch that uses it.
- **Published negative result worth having: antibody-specific language models are useless
  for framework stability.** AntiBERTy, AbLang2 and AntiBERTa2 all simply recover the
  wild-type/germline residue in the framework with high probability, carrying no stability
  signal; generic ESM-2 (t30 150M) was more varied, and structure-based Rosetta more varied
  still. Do not reach for an AbLM here.
- Also: 7mn8 (+ 8q6j to fill HER2's missing residues) is an **independent second structure
  of the same interface** we define from 1N8Z. `log0` next-action 3 asks for a structural
  biologist to check our four hotspot atoms; superposing 7mn8 onto 1N8Z answers a large part
  of that computationally.

**VHH format:** NanoMelt / NbThermo / TEMPRO / NBsTem — already surveyed in `log1` Session 1.

### D. Ranking and filtering — where this project actually stands to lose

This is the part of the literature that matters most to us, because `log0` already showed
our ranking signal is near-worthless and the whole pipeline's value rests on the reranking.

**🔴 Replicated across two independent competitions: what predicts *expression* is amino-acid
composition, not structure confidence.**

| Signal | Adaptyv EGFR (ROC AUC, expression) | Bits to Binders (ROC AUC, recovery) |
| --- | --- | --- |
| Glu content / K+E α-helix fraction | **0.77** (`composition_GLU_binder`) | **0.91** (K+E helix residue fraction) |
| Lys content | **0.73** | — (folded into the above) |
| DNA sequence entropy | — | **0.90** |
| DNA AA-dinucleotide repeats | — | **0.88** |
| EE (glutamate) repeats | — | **0.82** |
| α-helix propensity | 0.67 | — |
| **ipTM** | **0.58** | "did not predict experimental outcomes" |
| **pLDDT** | **0.55** | "did not predict experimental outcomes" |
| **ESM-2 PLL** | **0.44** | 0.85, but authors attribute it to correlation with entropy/repeat artefacts |

Two different labs, two different targets, two different assays, same answer: **Glu/Lys
composition beats every structure-confidence metric at predicting whether a design will
even survive to be measured.** Bits to Binders' retrospective filter (GC ≥60%, DNA entropy
<1.9, ≥45 AA repeats, ≥8 E repeats, ≥30% residues in a K+E helix) lifted recovery from
56.8% → 81.2%; adding a "no cysteines" filter raised functional rate from 5.9% → 10.6%.

**What predicts *binding* — weakly, and not what we assumed:**
- Adaptyv EGFR, binding ROC AUC: `esm_pll` highest; Foldseek similarity to known PDB
  structures `fident_pdb` 0.68 / `alntmscore_pdb` 0.66; pLDDT 0.66; **ipTM 0.64**. The
  authors warn this may just reflect optimised known binders being overrepresented.
- ANOVA, binders vs non-binders: **pDockQ/mpDockQ and interface area were the only metrics
  relevant across all molecule types** (F ≈ 9.9–19.4). ipTM/iPAE/ESM-logprob had high F
  only for peptides.
- **ESM-2 log-probabilities did not correlate with KD at all; ESM3 and ESM-C did correlate
  strongly with binding category.** If we use a PLM likelihood, it should be ESM-C, not ESM-2.
  ✅ **Directly runnable:** Nebius `esmc-300m` / `esmc-600m` produce "embeddings **or variant
  scores**", i.e. the likelihood head is exposed there. JapanFold's ESM-C is embeddings-only,
  so use the Nebius lane for this, not JapanFold.
- Adaptyv's own conclusion on their ranking metrics (AF2 iPAE, ipTM): *"valuable preliminary
  screens but insufficiently predictive of true experimental binding probability or
  affinity … one could potentially filter out a large amount of viable binders."*

**🔴 ipSAE — the metric `log1` Session 1 recommended as the primary ranking signal — has
mixed support.**
- For it: it is genuinely what Adaptyv consumes. JapanFold's `write_pae` help text says the
  PAE matrix is *"the one Adaptyv's ipSAE pipeline reads"*, and Adaptyv selected **600 of
  ~1,196 Nipah designs by best Boltz-2 ipSAE score** (the rest: 400 expert-curated, 200
  community vote).
- Against it: Bits to Binders measured ipSAE directly and found it, along with LIS, SAP,
  Rosetta InterfaceAnalyzer, pDockQ and ProteinMPNN likelihood, **"lacked substantial
  predictive power."** Adaptyv themselves noted ipSAE **reproducibility varied across
  hardware and seeds**, and added expert-selected designs specifically to mitigate that.
  Per-arm outcomes (ipSAE-selected vs expert-selected) have not been published.
- Note the model detail: Adaptyv ranked Nipah on **Boltz-2** ipSAE, whereas `log0` measured
  Boltz-2 ipTM at ROC AUC 0.645 and AF3 at 0.754 on an antibody–antigen benchmark with real
  non-binders. ipSAE from an AF2/AF3-family PAE is the better-founded choice, and `log0`'s
  benchmark is the most task-relevant evidence we have.

**🔴 Antibody ΔΔG prediction is not fit for ranking. Two independent 2025–26 papers:**
- **Graphinity** (Nat Comp Sci 2025, PMID 40629090) reaches Pearson 0.87 but *"appears to be
  overtraining on the few hundred experimental data points available and performance is not
  robust to train-test cut-offs"*; concludes **orders of magnitude more data are needed**.
- **AbAgym** (mAbs Nov 2025, PMID 41277441) curated 324k mutations from 68 DMS datasets and
  benchmarked ΔΔG methods: force fields achieved *"modest performance"*, ML models performed
  *"only marginally better than random."*
- **MochiBind** (iScience 2026, PMID 42564542) adds that existing affinity predictors
  *"lack generalization, with performance deteriorating for antibodies targeting antigens
  absent from training data or datasets lacking non-binders"* — and that a sequence-only
  model matches structure-based ones including Boltz-2 and Graphinity.
- ⇒ **Nebius `stab-ddg` should not be used as a ranking signal.** `log0` already noted it
  measures binding rather than folding ΔΔG; this is the stronger objection.

---

## Part 3 — What is runnable on our stack today

| Capability | Best available method | Where | Cost (p90) |
| --- | --- | --- | --- |
| De novo VHH, ranked + sequenced | **BoltzGen `nanobody-anything`** | JapanFold | ~$0.04/job |
| De novo antibody | BoltzGen `antibody-anything` | JapanFold | ~$0.04/job |
| **Redesign residues of an existing binder** | BoltzGen `protein-redesign` | JapanFold | ~$0.04/job |
| CDRH3 redesign on HER2 | RFD4-Proteina `her2-novel` | Foundry (VM) | — |
| All-atom antibody/nanobody CDR design | RFD4 `rfd4-binder`, recipe `high/front_loaded` | JapanFold | ~$0.026/job |
| Deterministic design rescoring + PAE matrix | **AF2-IG `write_pae: true`** | JapanFold | **$0.00084/target** |
| Antibody-antigen co-fold, second opinion | **OpenDDE · Antibody-Antigen** (13,671 jobs run) | JapanFold | $0.011/target |
| Epitope-restrained folding | `use_tfg_guidance` on OpenDDE | JapanFold | — |
| AF3-family ranking (`log0`: AUC 0.754) | `openfold3`, `protenix-v2` — **⚠️ real AlphaFold3 is NOT on our key** | Nebius | — |
| Antibody inverse folding | `antifold` | Nebius | — |
| Fab ΔΔG_folding scan | **nothing — not on any platform we have** (see note below) | — | — |
| VHH Tm | NanoMelt (CC BY-NC-SA) | local pip | free |
| PLM embeddings | ESM-C 300M/600M/6B, SaProt 650M/1.3B | JapanFold | ~$0.006/job |

**⚠️ Provenance of this table.** The JapanFold rows were verified this session against a
live `list_models` call. The Nebius rows were verified against `forge.py models` on our key
(37 models). Everything else is inherited from `log1` Session 1 or from skill descriptions.
Three rows were wrong on the first pass, because they were taken from skill descriptions
rather than checked: real **AlphaFold3 is not on our Nebius key** (only `openfold3` and
`protenix-v2`), there is **no Rosetta/PyRosetta/Flex ddG/FoldX** on it either, and **ESM-C
variant scores are available** on Nebius after all. All three are fixed above.

**🔴 The Rosetta ΔΔG_folding scan that Ng et al. used is NOT runnable on our stack.**
`forge.py models` shows no Rosetta, PyRosetta, Flex ddG or FoldX on our Nebius key; the only
ddG model is `stab-ddg`, which is *binding* ddG. `log0`'s original statement — "there is no
folding-stability predictor in the Apheris or Nebius stacks for Fab format" — stands, and an
earlier draft of this document contradicted it in error. **Consequence for the headline
finding in Part 4:** we can *use* Ng et al.'s published, experimentally validated
`VH S85N+R87T`, but we cannot *extend* their scan to search for better framework mutations
ourselves. The framework lever is limited to what the literature already validated.

**Correction to `log1` Session 1:** RFD4's `max_residues` is **512**, not 1024. The 1024
figure is the global `reference_residues` and applies to AF2-IG and BoltzGen. RFdiffusion3
and PXDesign allow 1536; Boltz-2 1920.

---

## Part 4 — The three things this survey changes

1. **The highest-value unused lever is the framework, not the CDR.** Adaptyv's own
   competition was won by keeping CDRs and mutating framework; the one experimentally
   validated Tm improvement for trastuzumab (`VH S85N+R87T`) is a framework mutation that
   preserves HER2 binding and ADCC. Our novelty gate is scored on CDRH3, so framework edits
   are orthogonal to it — spend the mandatory ≥3 CDRH3 substitutions as cheaply as possible
   (distance 3 = `log0`'s 90% binding tier) and take the differentiation from the framework.
   Blockers: the HER2 RFD4 checkpoint cannot touch framework (so these are hand-specified
   literature mutations, not generated); the sequon liability above; and whether the rules
   permit framework edits at all is **unconfirmed**.
2. **Add a composition/expression filter before ranking.** Glu/Lys composition beat every
   structure-confidence metric at predicting expression in two independent competitions
   (AUC 0.77 and 0.91 vs 0.55–0.58 for ipTM/pLDDT). This is cheap, stdlib-computable, and
   orthogonal to everything in `pipeline/` today. "Design for expression, not just binding."
3. **Lower the confidence placed on any single ranking metric, ipSAE included.** AF2-IG with
   `write_pae` is still the right primary signal for our task — it is deterministic, costs
   $0.00084/target, and `log0` measured the AF2/AF3 family as the only family that
   discriminates antibody binders. But Bits to Binders measured ipSAE as non-predictive,
   Adaptyv found it hardware/seed-sensitive, and the EGFR competition winner was ranked
   322nd–623rd by the organisers' own metrics. Rank on an ensemble, keep diversity in the
   submission set, and never let a metric cut designs the measured FACS/SPR data likes.

## Open questions this survey did not resolve
- Per-method hit rates in the Nipah competition exist only in the downloadable collection
  (`proteinbase.com/api/proteins/download?collectionId=019be357-…&slug=nipah-binder-competition-results`),
  not on any page. The secondary claims of "BoltzGen 3%, RFdiffusion 22%, Mosaic 89%" are
  **unverified and should not be repeated** — the pages show design counts only
  (BoltzGen 121, BindCraft 30, RFdiffusion 19, Germinal 14), with no binder counts.
- Whether Adaptyv's ipSAE-selected arm outperformed the expert-selected arm: unpublished.
- The actual ΔTm of `VH S85N+R87T` is in a figure panel (DSF melting curves), not the text.
- DualGPT-AB code availability — not found; check the paper's Code Availability section.
