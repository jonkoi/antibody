# One-shot de novo VHH hail Mary

Authorised by user, 14 Nov 2026. One 12-design batch, antibody-specific RFD4-Proteina checkpoint `antibody-2026-10-30`, seed 73, high/front_loaded, hotspot resampling off. No expansion into another generator campaign.

The prior BoltzGen campaign did not survive unguided epitope checks. This attempt changes the generator and uses explicit crystal hotspot atoms (Pro572 CD, Phe573 CD1, Ile591 CG2, Lys593 NZ), not sequence-only binding annotation. Conditions bias design; they do not establish competitive binding.

The shipped caplacizumab scaffold has expandable pseudo-residues rejected by the installed CLI (“Unresolved link valence”); no designs were generated from it. The corrected input uses the complete canonical-disulfide VHH from `runs/denovo-vhh/boltz2_unguided/cont_full_s1_r15.pdb` as framework geometry/sequence and crystal HER2 C489–580/D591–607 as antigen geometry/sequence. All three CDRs are redesigned; canonical framework Cys retained. The old VHH is a scaffold, not a positive HER2 epitope binder, and is independently refolded as a failed-route comparator.

`prepare.py` regenerates the corrected input. `design.yaml` specifies kept framework spans H1–22, H36–49, H60–96, H115–125 and redesigned loop lengths 13, 9–10, 12–21. Total chain remains under 250 aa. New sequences must independently pass abguard and be classified de_novo; VHH format alone does not establish track eligibility.

Independent check is unguided Boltz-2+MSA, 3 diffusion samples per design, seed 73. Compare actual 5 Å residue contacts against the 19 crystal epitope residues and four core residues; include paired WT positive control and old off-site VHH comparator. Three samples share a seed/MSA. No confidence-only acceptance and no claim of experimental binding.

Prespecified contact summaries: strict = ≥0.85 recall + 4/4 core; VHH overlap = ≥0.50 recall + 4/4 core (a smaller VHH footprint need not cover the whole Fab epitope). Both reported; the latter is a weaker structural compatibility criterion, not a predicted competition-assay pass. Require reproducible overlap before recommending any experimental candidate.

NanoMelt, if used, is VHH-only and licensed CC BY-NC-SA 4.0. Its apparent Tm predictions do not establish improvement relative to the competition's trastuzumab assay baseline. No measured affinity, thermostability or epitope competition is available for new designs.

## Generation and screening completed

12 unique designs generated, lengths 118–128 aa. All retain exactly two cysteines; generated S–S separations 2.214–2.377 Å and maximum adjacent peptide C–N distance 1.421–1.506 Å. This checks continuity/geometry, not experimental disulfide formation.

Live abguard with explicit `/home/hackathon/team15-claude/her2_novelty.duckdb`: 12/12 valid; 7 classified `de_novo`, 5 `novel` (optimisation). `screen_valid_denovo_7.csv` is a screening pool, not a recommended submission. Identity to trastuzumab for the de novo set ranges 0.6667–0.6953.

Independent fold job `14c0fd3d006eb13479a0407a7a013340`: only those 7, plus WT paired-Fv positive control and the prior off-site VHH comparator, 3 diffusion samples each. No further generation batch launched.

NanoMelt uses the existing conda environment `/home/hackathon/mamba/envs/nm`. The separate `nmelt/.venv` lacks pdbfixer and the conda binary initially lacked hmmscan on PATH; retry adds the conda bin directory. These are setup issues, not biological failures.

## Initial independent result

7 de novo candidates × 3 samples, plus two controls. Six new candidates have median epitope recall 0.0. Only `dn_vhh_her2_hail_mary_1` retains overlap: 10/19, 15/19, 14/19 crystal residues, core hits 3/4, 4/4, 4/4. Thus 2/3 meet the prespecified VHH-overlap gate and 0/3 meet the stricter Fab-footprint gate. WT control passes 3/3 with 19/19; prior off-site VHH comparator passes 0/3 with median recall 0.0. Initial fold charge $0.018392.

For candidate 1, both HER2 fragments jointly fit the crystal at 1.163–1.458 Å CA-RMSD; WT controls 0.937–1.110 Å. This checks target distortion rather than reading confidence alone. Candidate sequence verified unchanged in every independent structure.

Candidate 1: 128 aa, IMGT CDRH3 `AKWNRSASGWSSSNFWEAMDY`, live abguard identity to trastuzumab 0.6953125 (just below the de novo threshold), minimum reported CDRH3 distance 14. NanoMelt apparent Tm **69.68 °C**, not evidence of a gain over trastuzumab. Potential N-glycosylation motif **NRS at sequence positions 100–102 lies in CDRH3**: unlike the optimisation framework sequon, a glycan here could directly alter binding. Folding used unglycosylated sequence. Additional liabilities include Met/Trp oxidation sites. Abguard accepts the design, but these risks remain relevant.

One different-seed confirmation was submitted for candidate 1 plus WT only (no new sequence generation), job `f35e670f4de1e16c04091a6639b73aaa`, seed 1073, three diffusion samples each. This is a structural hail Mary, not an S-tier prediction. `exploratory_candidate_1.csv` and FASTA record the only candidate worth reviewing from this batch.

## Final decision

Different-seed confirmation (1073) again gives candidate 1 2/3 VHH-overlap passes, 0/3 strict Fab-footprint passes; recall 10/19–15/19 and minimum 3/4 core. Combined: **4/6 VHH-overlap passes across two seeds**, 0/6 strict passes. WT control passes 6/6 (19/19 in each). Candidate target joint CA-RMSD in confirmation 1.224–1.409 Å; the target fold remains near crystal. No affinity inference is made from this observation.

**One exploratory candidate is worth colleague review:** `dn_vhh_her2_hail_mary_1`, exact CSV and FASTA in `exploratory_candidate_1.*`. It offers substantially more independent epitope overlap than the failed earlier campaign, but is not a proven binder or a credible quantified S-tier prediction. A partial VHH footprint can be compatible with competition, which must be measured. The CDRH3 NRS sequon is a material binding risk in a glycosylating system; do not describe it as a stabilising framework glycan. Do not call 69.68 °C an improvement over trastuzumab without matched measurements.

Final hosted confirmation charge $0.004931; total hosted fold charge **$0.023323**, plus existing VM compute. NanoMelt used the existing VM install. One generation batch only; no ongoing design jobs and no competition submission.

Results: `epitope_samples.csv`, `confirmation_epitope_samples.csv`, their summaries, `target_geometry.csv`, live abguard reports, `nanomelt.csv`, and `decision.json`. Source structures and request/job/usage records are saved alongside them. All 33 independent fold structures have input sequences checked against submitted requests.
