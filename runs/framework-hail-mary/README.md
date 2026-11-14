# Framework stabilisation hail Mary

Prepared 14 November 2026. This is an exploratory panel, not the final submission.

Hypothesis: combine the experimentally studied trastuzumab VH S85N+R87T framework pair with three novel CDRH3s already shortlisted for affinity and epitope retention. Framework and CDR changes may interact; additivity is unproven.

Source: Ng et al., mAbs 2025, DOI https://doi.org/10.1080/19420862.2025.2532117 (PMID 40653772). The paper reports increased DSF melting temperature and retained HER2 cell binding/function for the framework pair on its original antibody background. It does not validate these new combinations or establish improved SPR affinity. Authors explicitly cannot exclude an effect of ectopic glycans.

## Construct audit

The local 120-residue VH has S85 and R87 at one-based sequence positions; the six WT identities in the paper's three named pairs (A40/K43, R50/R59, S85/R87) all match this sequence indexing. Our edits are explicitly sequence positions 85 and 87, not an automatic interpretation of IMGT position labels. The new local motif is NLT at 85–87 (WT SLR). No light-chain edits; lengths remain VH 120 + VL 107. `sequence_audit.json` lists every edit and sequon.

Panel: three original-framework comparators, the three corresponding framework combinations, and framework-only WT-CDRH3 control. `framework_candidates_3.csv` contains only the three experimental combinations.

## Live validator result

`validation.log`, `accepted.json`, and `liabilities.json` come from the VM's real abguard, with `/home/hackathon/team15-claude/her2_novelty.duckdb` explicitly selected.

- All three combinations and all three matched originals pass: 6/7 rows.
- Framework-only WT-CDRH3 control correctly fails novelty (0 substitutions).
- Each combination adds an N-glycosylation liability in FR3. It is reported, not rejected by this validator. The local expression heuristic's sequon rejection must not be described as an organiser hard rule.
- Tool acceptance does not confirm organisers' policy on framework edits or literature-derived combinations.

## Independent structural check

`check_epitope.py` submits an unguided Boltz-2 co-fold with MSA: the seven panel rows plus native WT control, 3 diffusion samples each, seed 42. Light/heavy/HER2-C/HER2-D chains A/B/C/D; target fragments copied by sequence from the existing WT control structure (92 + 17 residues). No epitope guidance or template is supplied. Request and job ID are saved. Structural contacts cannot predict melting temperature or the contribution of glycans; these inputs do not contain the new glycan.

## Interpretation and remaining questions

Before experimental-slot allocation, confirm framework edits are permitted and identify the actual expression system, glycosylation capability, and assayed format. The published stability gain may not transfer to an unglycosylated variable fragment. Do not transfer the old CDR-only affinity probabilities as calibrated probabilities for the framework combinations. Do not call any design S-tier before experimental measurements.

If checks remain favourable, the Lys-CDRH3 framework combination is the first balanced exploratory choice; the His-CDRH3 combination represents the strongest original model affinity hypothesis. Preserve original-framework comparators. Do not replace the main shortlist wholesale.

## Completed structural result

Job `92eb351f2567f7fccf0e61ebd7e5f27e` succeeded: 8 targets × 3 diffusion samples. Total charge **$0.021805** (`usage.json`). Raw structures, API result rows, request parameters and contact measurements are saved locally. `python3 runs/framework-hail-mary/check_epitope.py analyse` reproduces the contact tables without network access.

Success means ≥85% of the 19 crystal epitope residues contacted within 5 Å and all four core hotspot residues contacted, counting either antibody chain. These are residue contacts, not specific hotspot-atom distance tests. Samples share seed/MSA and are not three independent-seed experiments.

| CDRH3/background | Original framework | S85N+R87T framework |
|---|---:|---:|
| WGGKGPYWMD | 3/3 (minimum recall 1.000) | 3/3 (minimum recall 1.000) |
| WGGMGPYWMD | 3/3 (minimum recall 1.000) | 3/3 (minimum recall 0.947) |
| WGHDGPYWMD | 1/3 (minimum recall 0.263) | 1/3 (minimum recall 0.316) |
| WT CDRH3 control | 3/3 (minimum recall 0.947) | 3/3 (minimum recall 1.000) |

`recommended_exploratory_2.csv` contains only framework combinations WGGKGPYWMD and WGGMGPYWMD. Prefer the Lys option first: its original CDRH3 model affinity rank is higher and it avoids introducing another Met oxidation site. This is a qualitative selection, not a calibrated S-tier probability. Keep the His option in the investigation record, but do not prioritise its framework combination for scarce experimental slots.

The His original-framework result conflicts with the pre-existing untracked six-replica summary (6/6). Our new structures pass sequence assertions and the WT controls pass. This establishes variability under this sampling setup; do not silently merge datasets with unknown seed/sample provenance or claim the older result was erroneous. The matched 1/3 vs 1/3 result does not isolate a harmful effect of the framework edits.

Organiser clarification to obtain: “Are VH framework substitutions permitted alongside a novel CDRH3 in the optimisation track, and will the assay construct/expression system support N-linked glycosylation? Our proposed framework pair creates an NLT sequon; abguard accepts the combinations, but the published thermal benefit may not transfer without glycans.” No external message has been sent, and no designs have been submitted to the competition.

## Organiser clarification relayed by user

The user confirmed that mutations anywhere are permitted provided the design passes abguard screening. Framework eligibility is therefore resolved for this panel, whose three combinations already passed the live check.

The user also confirmed cell-free expression and said they think N-linked glycosylation is supported. Glycosylation support remains tentative, not verified. Even a capable system does not establish occupancy at the introduced NLT site or transfer of the published thermal benefit. These remain exploratory designs; no additional eligibility confirmation is needed for framework changes.
