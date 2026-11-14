# Correlated-risk review: F104P + A106W in the five enumeration track designs (enum_1 to enum_5)

Compiled 2026-11-14 from Amass BiomedCore plus geometry computed on
`inputs/12-her2/her2_fv_1n8z.cif` (chain B = heavy; CDRH3 = W99 G100 G101 D102 G103 F104
Y105 A106 M107 D108; sequential numbering, which matches the PDB auth numbering in this file).

## UPDATE 2026-11-14 (after the Boltz-2 refold of enum_1 to enum_5; supersedes the backbone argument in section 2)

Verified on `runs/enumeration track/boltz2/*/seed_42/*_model.cif` (same dihedral code):

| model | phi(101) | phi(104) |
| --- | --- | --- |
| 1N8Z crystal (wild-type) | +106.7 | -141.3 |
| Boltz-2 wild-type control | +77.6 | -99.4 |
| enum_1 / enum_2 / enum_3 / enum_4 / enum_5 | -68 / -93 / -106 / -54 / -71 | -64.7 / -64.1 / -78.4 / -59.4 / -59.1 |

- phi(104) is -59 to -78 in all five: proline-compatible. **The "forbidden by the crystal backbone"
  argument (section 2 and the enum_5 / G101P point in section 5) does not survive the refold and is
  withdrawn.** It rested on the wild-type conformation and was caveated as such.
- The model did **not** keep the crystal loop. phi(101) flips from positive to negative in every
  design, so the predicted CDRH3 is a rearranged loop, not wild-type with two edits.
- Weight to put on it: a co-folding model produces backbone angles consistent with the sequence it
  is given, so "phi is proline-compatible" is close to guaranteed and says little about binding.
  The informative readouts are the hotspot distances (reported 2.2 to 3.0 A vs 2.63 A for control)
  and per-design confidence (ipTM / iPAE vs the control). The Gly103 O to Lys593 NZ distance
  (2.42 A in the crystal) is the single best contact to check.
- Still standing: both positions pack against VL (section 2), A106W adds a second CDR3 Trp
  (section 4), and no experimental data exist for the pair (section 1).

## The five designs (runs/enumeration track/submission.csv)

| id | CDRH3 | substitutions vs trastuzumab | own predicted KD |
| --- | --- | --- | --- |
| enum_1 | WGHDGPYWMD | G101H, **F104P, A106W** | 0.50 nM |
| enum_2 | WGGKGPYWMD | D102K, **F104P, A106W** | 0.77 nM |
| enum_3 | WGGMGPYWMD | D102M, **F104P, A106W** | 0.89 nM |
| enum_4 | WGGRGPYWMD | D102R, **F104P, A106W** | 0.95 nM |
| enum_5 | WGPDGPYWMD | **G101P**, **F104P, A106W** | 1.52 nM |

All five pass the novelty gate (`is_novel_enough: true`, distance 3 from 524,346 references).
Each differs from the others by one residue, so this is **one design hypothesis tested five
times**, not five independent bets. If F104P/A106W fails, all five fail together.

## Bottom line

No Amass record reports F104P, A106W or any alanine scan of trastuzumab CDRH3. The literature
cannot validate this pair. The wild-type crystal backbone argued against the pair, but the Boltz-2 refold
(update above) removed that objection. Treat the submission as a correlated gamble and decorrelate it.

---

## 1. Mutagenesis / affinity-maturation of trastuzumab (hu4D5) CDRH3

| study | what it reports | bears on F104 / A106? |
| --- | --- | --- |
| Gerstner 2002, 20-aa saturation at selected hu4D5 CDR positions; one Trp-for-Asp change gave ~3-fold gain. `AMBC_MHGn8BBeOsu3pBUwK0R7QwVdBJS` PMID 12206766 doi:10.1016/s0022-2836(02)00677-0 | variability map vs earlier alanine scan | Abstract names no CDRH3 positions. Per-position data would be in the full text. |
| Moon 2016, CDR-H3/L3 randomised hu4D5 scFv; AH06 has 7.4-fold better KD (60 pM). `AMBC_TJrq1YVG3lzNqIHRh9sjhfQX2rs` PMID 26743905 doi:10.14348/molcells.2016.2235 | selected CDR-H3 residues are **W98, F100c, A101, L102** (Kabat); F100c, A101, L102 credited with heavy/light-chain hydrophobic interactions | **Yes, indirectly.** If I map F100c to Phe104 and A101 to Ala106 (my mapping, not the paper's) the winner **kept Phe and Ala** where you put Pro and Trp. The abstract does not say which positions differ from wild-type. |
| Nguyen 2018, CHO Fab display of a low-affinity hu4D5 variant; up to 400-fold gains. `AMBC_7KVnsst4I9HeHylLHqmlwLkT6ZT` PMID 29566240 doi:10.1093/protein/gzy004 | 3 of 4 improved clones carried free Cys in CDR-H3 | No. Shows CDR-H3 tolerates change when starting from a weak parent. |
| Balakrishnan 2021, ML on 193 trastuzumab single mutants; CDR-H3 mutations act partly through secondary-structure change. `AMBC_MR053BqgDwChZk2sMpMs9p875oL` PMID 34392800 doi:10.1080/07391102.2021.1961866 | a per-position dataset exists | Record does not list F104 or A106 effects. Worth pulling the dataset. |
| Eigenbrot 1993, X-ray of 4D5 variants; CDR-H3 differs up to 8 A between copies. `AMBC_6ozdxIDDNijfNoA8zrLh0filWP2` PMID 8095303 doi:10.1006/jmbi.1993.1099 | trastuzumab H3 is intrinsically variable in the free Fab | Supports conformational selection (section 3). |

Not trastuzumab: the shotgun alanine scan `AMBC_7vdLCwbOWsiqFHJ5dNI7tGd4SqK` (PMID 12079396) is Fab2C4 (pertuzumab).

## 2. Is Phe104 a hotspot? Structure says "buried VL-packing residue", not a direct hotspot

Computed from 1N8Z:

- **Phe104 nearest HER2 atom: 4.27 A** (backbone N to Lys593 NZ). It is a marginal HER2 contact.
  *(This corrects an earlier draft of `epitope-1n8z-verification.md` that listed 3.1 A; the
  table there is now fixed: Trp99 4.05, Phe104 4.27, Tyr105 3.32 A.)*
- Its side chain packs on the **light chain**: Tyr49 (3.4 A), Leu46 (3.8 A), Tyr55 (4.1 A);
  191 neighbour atoms within 6 A of its 7 side-chain atoms, i.e. buried.
- Ala106 touches no HER2 atom within 4.5 A (nearest 6.83 A) but its C-beta sits 3.9 to 4.1 A
  from light-chain Tyr49, Leu46 and Ala34 and 3.76 A from Tyr36. It also faces the VL.
- Moon 2016 independently credits F100c with heavy/light hydrophobic contacts.

So both positions are **VH/VL-interface packing residues**, not epitope hotspots. That is a
different risk from "we broke a contact": the swap can destabilise the paratope geometry.
The F104P + A106W pair removes a buried aromatic and puts a large indole where an Ala nests
against the VL. My inference: that clashes with Tyr36/Tyr49 unless the loop rearranges.
No loop model exists in the repo for these five, so this is untested.

**Backbone geometry (validated dihedral code; framework beta-strand gives phi -113 to -155):**

| residue | phi | psi | consequence |
| --- | --- | --- | --- |
| Gly101 | **+106.7** | +152.0 | positive phi: Gly-only region. **G101P (enum_5) is incompatible with this backbone.** |
| Asp102 | -88.0 | -7.3 | alpha-like; D102K/M/R are backbone-tolerant |
| Gly103 | -137.7 | +42.1 | carries the 2.42 A carbonyl H-bond to Lys593 |
| **Phe104** | **-141.3** | +93.7 | **proline requires phi about -65; this is about 75 degrees away** |
| Ala106 | -105.5 | +154.4 | beta; Trp is backbone-tolerant |

F104P cannot adopt the 1N8Z conformation as is. Any binding must come from a different loop
conformation, and the residue before the new Pro is Gly103, the carbonyl that anchors Lys593.
This is geometry on the wild-type crystal conformation. A designed loop may differ.

## 3. Proline in antibody CDR loops

No Amass record reports a Pro substitution in a CDR-H3 with binding and expression measured
together. Closest evidence:

- **Isomerisation heterogeneity:** a YPP motif in a CDR-H3 produced a two-peak SEC profile from
  proline cis/trans isomerisation; potency was not significantly affected.
  `AMBC_UWZc7Cq2njWUFJtySGquWRR2t7f` PMID 31791173 doi:10.1080/19420862.2019.1698128
- **Antigen-triggered CDR-H3 Pro isomerisation** (Pro103), with the adjacent Asp102 rotating.
  `AMBC_686XfOAdf4QSMTmCWKSrGR1GeZ6` PMID 29208911 doi:10.1038/s41598-017-16766-8
- **Position matters:** loop prolines in an adalimumab Fab shifted Tm by position (P41G lowers it,
  K209P raises it). Constant/framework loops, not CDRs.
  `AMBC_EAT3RHejXN7zVPO0mZdZlVA4sSF` PMID 42480179 doi:10.1016/j.bbrc.2026.154312
- **CDR-L3 cis-Pro** maintains stability in V-kappa (light chain). `AMBC_YfM2ibvLSIPp6E9iHChR16Qkyo4` PMID 20727359 doi:10.1016/j.jmb.2010.08.009
- **Rigidification does not reliably help.** Conflicting repertoire work: only a slight decrease in
  CDR-H3 flexibility on maturation (`AMBC_J27FcmcKP4CwVRNuTwOpuN7qBCb` PMID 29545810
  doi:10.3389/fimmu.2018.00413) vs substantial rigidification in 10 pairs
  (`AMBC_9dnqJvRzyIyliUSBgH11KzsDpTY` PMID 32850970 doi:10.3389/fmolb.2020.00182) and
  QSFR rigidity of CDR H3 (`AMBC_Rp9hiPl5X2FCEjbzv9OULpxsSAm` PMID 26132144 doi:10.1371/journal.pcbi.1004327).
- **CDR-H3 binds by conformational selection:** the binding-competent conformation pre-exists in the
  ensemble. A Pro that excludes it removes binding; one that enforces it can help.
  `AMBC_4FDmCzhkEU26tTlaKIbIXSM4FOS` PMID 31148507 doi:10.1080/19420862.2019.1618676
- Non-antibody analogues, opposite directions: Pro to Gly in an SH3 loop raised affinity 8-fold via
  flexibility (`AMBC_4L2PZxUfJtWW78Lq93JYLRAiMeP` PMID 17382937); Gly to Pro above a TCR CDR3 beta
  loop narrowed specificity with affinity unchanged (`AMBC_N8JKKqvVYcvUN1wmMQdke4MZLM6` PMID 38545094).

Usable tool: ITsFlexible classifies CDR3 loops as rigid/flexible
(`AMBC_3WNt7DR3hL8nb6QXrKszJ0nqZVI` PMID 41143207 doi:10.1038/s42256-025-01131-6).

## 4. Tryptophan in CDR3 as a developability liability

- CDR-H3 Trp oxidises when solvent-accessible, though susceptibility does not track accessibility
  linearly. `AMBC_Pu75OJ7Lt4EzKS5ig9TVylh7rhX` PMID 30608647 doi:10.1021/acs.analchem.8b04768
- Oxidation in or near CDRs raises variable-domain flexibility and disrupts binding.
  `AMBC_PvJNIa9tOf6oWK2D6JGd09vDbK6` PMID 30536043 doi:10.1007/s11095-018-2545-8
- Trp oxidation barely changes conformational stability but **promotes aggregation** over a wide
  range of conditions via covalent crosslinks. `AMBC_3YIace8LkcCtYFNxQqhaOPPrGAk` PMID 31678251 doi:10.1016/j.xphs.2019.10.051
- Heavy-chain CDR Trp100 oxidation drove CDR fragmentation (polysorbate-80 dependent).
  `AMBC_8Um38EyzvX5eA1fGOaE6dl2q0M4` PMID 35219863 doi:10.1016/j.ejpb.2022.02.015
- Susceptibility depends on solvent accessibility **and local electrostatic potential**
  (negative microenvironments raise risk); in one anti-CD33 antibody Trp to Phe abolished binding.
  `AMBC_CPdsaMpYQ1lPhZe2y8Mj1bc5iQ3` PMID 42590841 doi:10.1080/19420862.2026.2713231
- Only a subset of liability motifs actually modify, so experimental assessment is needed.
  `AMBC_Pw3yvGc8mFrVVEgL2RjiRB8hPBU` PMID 30118239 doi:10.1021/acs.molpharmaceut.8b00536
- Also: `AMBC_Sd5ej9sWSdBtbZfD4s4Pc1SVZDC` PMID 21390239; `AMBC_W0ka7h6RaNlpPCnFCcssmBJmXE8` PMID 27612038;
  `AMBC_QpARilhT4Alv1iMBbu9s5HsE5YR` PMID 21440514.

**Your own liability report** (`runs/enumeration track/all_liabilities.json`): wild-type has 2 TrpOx flags
(CDR3 W99, FR4 W103); each design has **3**, because A106W makes CDRH3 a two-Trp loop (medium
severity). W106 also sits two residues from Asp108; the electrostatics paper suggests negative
neighbours raise risk, which the authors' two-parameter model could test. enum_2 to enum_4 happen to
remove the CDR3 DG isomerisation motif (D102K/M/R), a real plus; enum_1 and enum_5 keep it.

## 5. What this means for the submission

1. **Do not submit five copies of one hypothesis.** The pair is unvalidated in the literature and
   the 1N8Z backbone argues against F104P.
2. **enum_5 is the weakest.** G101P (phi +107 for a Gly) plus F104P is two prolines in a 10-residue
   loop at geometries the crystal structure forbids.
3. **If you keep any, keep one**, ideally enum_1 (best own predicted KD) or one of enum_2 to enum_4
   (removes the DG motif), and fill the remaining slots with designs that retain Phe104.
4. **Checks that fit in the time:** refold the five with Boltz-2 (the repo already has the
   pipeline) and read off phi at 104 and the interface iPAE; run ITsFlexible; look at Moon 2016
   and Gerstner 2002 full text for the F100c / A101 positions.
5. The "predicted KD 0.5 to 1.5 nM" comes from your own model. It is the same model scoring all
   five, so its agreement across them is not independent evidence.

## Provenance

Geometry: heavy-atom distances at 4.5 A from the 1N8Z-derived file; dihedrals validated on a
framework beta-strand before use. Amass BiomedCore, 2026-11-14. No record I found reports
F104P or A106W on trastuzumab.
