# Team 05 — ranked HER2 antibody designs

Four-minute version · 7 slides · Oslo Antibodies × AI Hackathon · 14 Nov 2026

Sources: [team process overview](findings/team-process-overview.md) · [24-priority crosswalk](runs/reconciliation/priority24_crosswalk.csv) · [full ranked evidence](runs/reconciliation/ranked_with_evidence.csv) · [submission list](team15_01.csv) · [final SUBMISSION.csv](SUBMISSION.csv)

## 1 — The task and the real constraint (30 sec)

We redesigned trastuzumab’s CDRH3 to preserve HER2 binding while meeting novelty and developability screens. The measured library shows why this is hard: at three substitutions, 1,341 of 1,490 variants (90.0%) are still in the high FACS-binding class. That is a useful sequence-class signal, not an affinity measurement.

**The first generative batch did not find the conservative, novel region. That was a limitation of our search, not proof that the region was empty.**

## 2 — Search where the structure gives room (35 sec)

Crystal structure 1N8Z identifies four CDRH3 contact positions: W99, G103, F104 and Y105. We enumerated mutations across all ten positions, then used structure and sequence evidence to prioritize. The five F104P/A106W leaders preserve three of the four contact identities; they are hypotheses, not guaranteed binders.

The enumeration set contained 51,417 legal sequences, but we did not submit that pool wholesale. We collected designs that had actually been generated or shortlisted, screened them, and kept the evaluated candidates and backups.

## 3 — Use each model for what it measures (35 sec)

The team combined the measured FACS classification dataset, a separate small SPR affinity dataset, sequence models, structural refolds, StaB-ddG predictions and developability reports. These answer different questions. The colleague’s clonotype-aware FACS model reported test AUC 0.9842 and mixed-cluster macro-AUC 0.9306 across 44 clusters; FACS class does not equal KD.

The eight measured wild-type SPR values have median 0.1345 nM. Historical probability columns used a 0.146 nM threshold. A rough seven-fold contrast comes from this small, potentially mixed-batch sample; it is not a universal assay resolution floor. No new candidate has measured KD, competition or Tm.

## 4 — Screen first; interpret structures carefully (35 sec)

We deduplicated by format and full heavy/light sequences, then ran live abguard against the provided HER2 reference database and retained its track assignments. Of 653 unique candidates, 546 passed and 107 failed: 369 optimisation and 177 de novo. Screening establishes eligibility under this database and tool, not efficacy.

Structural checks compare epitope contacts with controls. Contact-residue identity, four-hotspot coverage, whole-epitope recall and one atom-pair distance are different measures. Fold confidence and predicted binding ΔΔG are useful evidence, not experimental binding results.

## 5 — Why the first 24 are ordered this way (40 sec)

We preserved the colleague’s 24 sequences and their order exactly, then appended all other validated candidates as a reserve. Ten of those priorities matched our earlier set; fourteen were newly shortlisted. The first 24 are a team decision combining different evidence and diversity, not a calibrated mathematical ranking.

The leading hypotheses include `WGGRGPYWMD` (low predicted binding penalty, with a liability trade-off), `SLGSADYAFD` (a distinct motif, favourable reported atom-pair geometry, five liabilities), and `WGHDGPYWMD` (strong historical sequence-affinity score but weaker latest matched contact replication). These are competing hypotheses; no single metric proves superiority. Framework variants test transfer of published S85N/R87T stability changes, but their Tm is unmeasured and the framework ddG is not a measured combined-construct result.

## 6 — Keep the de novo long shots honest (35 sec)

The 177 de novo-track candidates remain in the ranked list, with their track rank and evidence. A one-shot VHH design had partial HER2 overlap in 4/6 model samples across two seeds, but 0/6 met the stricter Fab-footprint criterion. Its predicted NanoMelt value is 69.68°C, not a measured Tm gain; its CDRH3 also contains an NRS glycosylation sequon. This is an exploratory candidate, not an S-tier claim.

We also appended two newer BoltzGen nanobody designs, `denovo_hotspot_01` and `denovo_hotspot_03`, generated with HER2 Pro572, Phe573, Ile591 and Lys593 set as binding residues. In the predicted structures all four hotspots have a binder atom within 3.3 Å and 4.6 Å respectively, but design→target ipTM is only 0.21–0.22 from one model and one seed, with no refold check, no thermostability score and 21–23 liability flags. Their sequences have no cysteines (the conserved Cys positions read as Glu), and they were not run through abguard. They are unvalidated long shots, not recommendations. An earlier unconstrained batch, with the epitope not enforced, was dropped in favour of these.

## 7 — What we learned and what happens next (30 sec)

We now have 546 screen-passing designs in the shared ranked CSV, including the colleague’s top 24 plus a broad reserve. `SUBMISSION.csv` holds 548 rows: those 546 plus the two unscreened hotspot-guided de novo VHHs at the end. We used literature to select tests and identify risks, while checking source data and keeping predictions separate from measurements. The most useful result is the traceable reasoning: why each sequence is present, what evidence supports it, and what could still be wrong.

**Next:** choose experimental slots from the first-ranked candidates and preserve representation across distinct hypotheses; then test binding, affinity and thermal stability directly.

### Sources for speaker Q&A

- [Process, reconciled reasoning and source map](findings/team-process-overview.md)
- [Original-priority IDs and new rank crosswalk](runs/reconciliation/priority24_crosswalk.csv)
- [All 546 candidates with evidence](runs/reconciliation/ranked_with_evidence.csv)
- [Full ranked inventory and limitations](runs/all-ranked/README.md)
- [De novo hotspot run outputs](runs/denovo/F_epitope.json) (hotspot distances for the 20 BoltzGen designs)
- [Original remote 24-row snapshot](runs/reconciliation/team15_01.remote.csv)
