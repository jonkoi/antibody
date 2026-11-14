# Combined process and reasoning — presentation source of truth

Snapshot: remote `ab8ade5` + local `bfb5dfe`, reconciled on 14 November 2026. Read this before older findings, whose later sections sometimes reverse earlier conclusions.

## What is actually in the list

`team15_01.csv` and `SUBMISSION.csv` now contain the same **546 screen-passing designs**: 369 optimisation, 177 de novo. The remote list's first **24 sequences and order are preserved**, and all remaining unique, validated designs are appended as backups. These are a ranked inventory, not 546 equally promising antibodies. All 24 original priorities are optimisation designs.

The union was 653 unique candidates: our previous 639 plus 14 new shortlisted entries; ten of the colleague's 24 already matched full sequences in our list. Fresh live abguard rejected 107. Exact format/heavy/light identity is the deduplication key; framework and disulfide variants remain distinct. The broad 51,417-sequence enumeration pool is **not** exported wholesale. Only explicitly shortlisted/evaluated proposals are included. Two controls are separate.

Reconciliation: [priority crosswalk](../runs/reconciliation/priority24_crosswalk.csv), [full evidence CSV](../runs/reconciliation/ranked_with_evidence.csv), [summary](../runs/reconciliation/sync_summary.json). The remote 24-row CSV and previous local rankings are preserved in `runs/reconciliation/`.

## How the work developed

| Workstream | Process and reason | What it established | Important limit / source |
|---|---|---|---|
| Original generation | HER2-finetuned RFD4; 4-design pilot then 100-design production; later hotspot-constrained seed-7 batch | Generated novel sequence/backbone proposals; conservative proposals in the first production batch often failed novelty | Generator confidence is not measured binding. `runs/team15-claude-n100-seed1/`, seed-7 sequence/metrics records in `runs/all-ranked/source/` |
| Exact sequence search | Enumerate three changes in the 10-residue window; remove candidates too close to provided CDRH3s; apply winner-position and liability heuristics | 51,417 sequences in the minimal-change legal search set; conservative novel designs exist | The search varied all ten positions. The five F104P/A106W leaders retain three of four crystal-contact residue identities, not four. `pipeline/enumerate_candidates.py`, `log0.md` |
| Binding-class models | Both approaches learn unary/pairwise sequence effects from measured FACS classes | Colleague's clonotype-aware model: test AUC 0.9842; within mixed-class cluster macro-AUC 0.9306 across 44 clusters. Seed-sensitive per-design order motivates tiers/diversity | FACS high is an assay class, not KD. Our sklearn/CV and colleague's split model use different evaluations; their AUCs are not directly comparable. `pipeline/predict_binding.py`, `findings/lane-d-binding-model.md`, `runs/ranked/facs_scores.log` |
| Affinity model | Separate small SPR regression, with broad uncertainty, combines with binding-class predictions | Gave a candidate-priority hypothesis around F104P/A106W | Historical scores targeted 0.146 nM; actual eight-replicate median is 0.1345 nM. Output labels now state the historical threshold. No new candidate has measured SPR affinity. `pipeline/score_candidates_ml.py`, `runs/reconciliation/verified_facts.json` |
| Novelty and liabilities | Real abguard with explicit provided-HER2 reference DB, plus LAP reports | Fresh union: 546/653 accepted; track assigned by the tool | Eligibility is not efficacy. A local sequon heuristic was stricter than abguard; sequons are reported, not universally disqualifying. `runs/all-ranked/{accepted,novelty,liabilities}.json` |
| Structural contact analysis | Crystal 1N8Z; four targeted HER2 hotspots; compare refolds against controls | Candidate-specific poses can disagree with confidence scores; target-contact checks catch off-site predictions | Retained sequence identity, four-hotspot coverage, 19-residue recall and a specific O–NZ distance are different measurements. None is a competitive-binding assay. `findings/interface-hbond-check.md`, `runs/final-shortlist/` |
| Binding energetics | StaB-ddG single-site scan and joint multi-site predictions on crystal complex | Extra computational evidence favouring some different candidates; raw output includes model Monte Carlo SE | Predicts binding ΔΔG, not folding/Tm. SE is not experimental uncertainty. Multi-mutant generalisation unverified; comparisons with an SPR noise estimate do not validate the predictor. `runs/stab-ddg/` |
| Framework stability | Transfer published trastuzumab VH S85N/R87T onto shortlisted loops | K/M variants pass screening and retain strict epitope contacts in 3/3 diffusion samples each | No Tm measurement; samples share a seed/MSA. NLT sequon may matter; a distant attachment site does not guarantee absence of glycan effects. `runs/framework-hail-mary/` |
| Earlier de novo VHHs | BoltzGen nanobody protocol, later restore canonical framework cysteines; NanoMelt scoring | Tested earlier leaders docked elsewhere despite confidence. Restoration increased predicted Tm in the model | Predicted Tm gains are not measured thermal improvements. Known off-site candidates remain low-ranked long shots at the user's request. `runs/denovo-vhh/README.md` |
| One-shot de novo retry | Antibody-specific RFD4, complete VHH framework, all 3 CDRs redesigned with explicit geometric hotspots; 12 designs | All 12 pass; 7 classify de novo. Candidate 1 has partial epitope overlap in 4/6 samples over two seeds | 0/6 strict Fab-footprint passes; no KD. NanoMelt 69.68°C is not a proven gain vs trastuzumab. CDRH3 NRS sequon may alter binding. `runs/denovo-hail-mary/decision.json` |
| Literature | AMASS and primary papers informed hotspots, oxidation risk, framework transfer, VHH thermal predictors | Helped choose interventions and qualify claims | Crystal geometry was computed from coordinates; AMASS did not itself measure it. No claim of exhaustive absence of prior art. `docs/sota-antibody-design.md`, literature-backed findings |

## Why the leading priorities differ from our earlier ranking

The earlier local broad order put reproducible K/M epitope contacts first. The colleague added an atom-pair geometry check, StaB-ddG outputs, a different FACS model and a diversity policy, elevating `SLGSADYAFD` from local rank 174 to team rank 2. This adds evidence; it does not demonstrate an affinity improvement.

| Team rank | CDRH3 / variant | Main reason | Counter-evidence to retain |
|---:|---|---|---|
| 1 | WGGRGPYWMD | Low predicted binding ΔΔG ~0.20; reported template O103–NZ593 3.48 Å; liability trade-off | Only single-run contact evidence in the verified local broad check; untracked replica summary not counted |
| 2 | SLGSADYAFD | Different generative motif; predicted ΔΔG ~0.42; reported O–NZ 2.51 Å; five liabilities | Six CDRH3 changes; a short model distance is not proof of stronger binding |
| 3 | WGHDGPYWMD | Strongest historical sequence-affinity score among these | Latest matched experiment only 1/3 strict contact samples on original and framework constructs |
| 4 | WGPDGPYWMD | Diversity/loop-rigidification hypothesis; reported O–NZ 3.20 Å | No measured Tm; proline-containing output having compatible backbone angles does not validate binding |
| 5–8 | M, framework-M, K, framework-K | Alternative energetics/contact/stability hypotheses; framework variants preserve relevant controls | M/K parent template O–NZ distances are long, while independent whole-epitope contacts remain reproducible. These criteria are not equivalent. Framework ΔΔG shown in older notes was the parent's prediction |
| 9–24 | Diverse FACS-ranked backups | Model-only reserve; minimum pairwise 10-mer distance 3; 14/16 retain all four contact residue identities | Sequence retention does not prove structural contacts; no new SPR/Tm data |
| Remaining | Full validated reserve, including de novo long shots | Completeness and optional slots; original heuristic order retained | Screening alone does not indicate quality; known off-site and untested candidates are labelled |

The first-24 order is a preserved team decision, not a newly calibrated mathematical ranking. Track-specific de novo rank remains separate. Neither a single H-bond-like distance nor low model ΔΔG should veto or prove a candidate in isolation.

## Presentation facts checked again from raw data

- Full FACS table: 524,346 rows. At exactly three substitutions: **1,341/1,490 = 90.0% high class**. The 90.9% number in older slides refers to the training split, not the full table.
- SPR: 190 quantified rows overall; 187 with 10-residue CDRH3s. WT eight values span 0.02851802–0.178 nM; **median 0.1345**, geometric mean 0.09539 nM. Historical probability columns explicitly use the old **0.146 nM** threshold.
- A ~7.10-fold contrast follows from a rough 95% single-measurement comparison using pooled log10 SD 0.307. This is a small, potentially mixed-batch estimate, not a universal assay detection limit or proof that smaller improvements cannot exist.
- Secondary 16 priorities: **14/16** preserve four contact-residue identities; minimum pairwise Hamming distance **3**.
- Fold-confidence benchmark: shuffled negatives are presumed mismatches, not SPR-confirmed non-binders. Report AF3/Boltz-2/Chai-1 AUC 0.754/0.645/0.568 with uncertainty; do not transfer model thresholds. The 50-replica no-improvement result was Chai-1, not a demonstration about every model.
- **No newly designed antibody has experimental KD, competition or Tm measurements in this repository.** All model-derived ΔΔG, structures and NanoMelt values must be labelled predictions.

Recompute scalar facts with `python3 runs/reconciliation/audit_facts.py`. Raw prediction CSVs, saved structures, validation reports and historical list snapshots remain linked in the evidence CSVs.
