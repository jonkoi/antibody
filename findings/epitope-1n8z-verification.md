# Was the 100-design batch aimed at the right surface patch?

**Verdict: yes.** All four hotspot atoms are genuine trastuzumab–HER2 contacts in
1N8Z, and all four are among only **seven** HER2 residues within 5 Å of CDRH3.
For a CDRH3 redesign this is a well-scoped target set, not a mis-aimed one.

Compiled 2026-11-14. Literature from Amass BiomedCore (amassId + PMID/DOI for every
record). Geometry computed directly from `inputs/12-her2/her2_fv_1n8z.cif`.

---

## 0. Method note — why this answer is stronger than a literature answer

Amass is a literature index, not a structural database, so it cannot tell you which
atoms touch which. It *can* tell you what the field asserts and what was measured
experimentally. The decisive evidence here came from computing contacts on the
structure file already in this repo. Chains: A = trastuzumab light (1–107),
B = heavy (1–120), C = HER2 489–580, D = HER2 591–607.

---

## 1. HER2 domain IV residues in the trastuzumab contact interface

**Computed, 4.5 Å heavy-atom cutoff, HER2 vs. whole Fv — 19 residues in three segments:**

| segment | HER2 residues in contact | contacted by |
| --- | --- | --- |
| loop 1 | Pro557, Glu558, Asp560, Gln561 | CDR-H2 (Arg50, Tyr52, Tyr57, Thr58, Arg59), Tyr33, light Thr94 |
| loop 2 | Lys569, Asp570, Pro571, **Pro572**, **Phe573**, Pro579 | **CDR-H3 (Trp99, Gly103, Tyr105)**, light His91/Tyr92/Thr93/Thr94, Tyr33, Arg50 |
| loop 3 | **Ile591**, **Lys593**, Asp596, Glu598, Ala600, Gln602, Pro603, Cys604, Pro605 | **CDR-H3 (Gly103, Phe104)**, light Asp28/Asn30/Thr31/Ser50/Phe53 |

This three-loop architecture is exactly what the primary structure paper describes
("Herceptin binds to the juxtamembrane region of HER2"):

- Cho HS, Mason K, Ramyar KX, Stanley AM, Gabelli SB, Denney DW, Leahy DJ.
  *Structure of the extracellular region of HER2 alone and in complex with the
  Herceptin Fab.* Nature 2003. **1,628 citations.** This is the source of 1N8Z.
  `AMBC_JeJdrXhjh7YpbBacDLhYNGDWz3d` · PMID 12610629 · doi:10.1038/nature01392
  *(Amass holds abstract only — `hasFulltext: false` — so the residue-level loop
  definitions are not retrievable from the record itself.)*

**Independent numbering cross-check.** An MD/MM-GBSA study that explicitly used
PDB 1N8Z reports the three binding loops as residues **991–995, 1004–1007 and
1027–1037**. Subtract a constant 434 and those become **557–561, 570–573 and
593–603** — i.e. the same three loops, renumbered by their simulation setup. The
offset is identical for all three loops, which confirms the canonical definitions.
- `AMBC_OCRs4Bt7QBxPiVt4ukowLu2fvEA` · doi:10.1096/fasebj.31.1_supplement.775.23
  (FASEB J conference abstract, no PMID, 0 citations — treat as corroborating
  arithmetic, not as an authority.)

Supporting structural records:
- Cryo-EM HER2–trastuzumab–pertuzumab ternary complex, 4.36 Å; confirms the two
  antibodies bind distinct subdomains with no cooperative interaction.
  `AMBC_EsuxDT0fSVGPZwXjYp1QuAAqryb` · PMID 31042744 · doi:10.1371/journal.pone.0216095
- Review placing trastuzumab's domain IV juxtamembrane epitope against the other
  four inhibitory ErbB antibody modes.
  `AMBC_3pEfFZH3UTILCbcVvpSyZxisuqx` · PMID 18992239 · doi:10.1016/j.yexcr.2008.10.008

---

## 2. Is Pro572 / Phe573 / Ile591 / Lys593 supported? — Yes, all four

**Computed geometry for the exact atoms the run targeted:**

| hotspot atom | nearest Fv atom | nearest CDR-H3 atom | verdict |
| --- | --- | --- | --- |
| `D/593/NZ` (LYS) | **2.42 Å** — Gly103 O | **2.42 Å** — Gly103 O | H-bond; **tightest contact in the entire interface** |
| `C/573/CD1` (PHE) | 3.42 Å — Tyr105 CZ | 3.44 Å — Tyr105 CZ | aromatic stacking with CDR-H3 Tyr105 |
| `D/591/CG2` (ILE) | 3.74 Å — Gly103 O | 3.74 Å — Gly103 O | genuine hydrophobic contact |
| `C/572/CD` (PRO) | 3.29 Å — light His91 O | 4.00 Å — Tyr105 CD1 | contact, but **primarily a light-chain contact** |

Only **seven** HER2 residues sit within 5 Å of CDRH3: Asp570, Pro571, Pro572,
Phe573, Pro579, Ile591, Lys593. The run picked four of those seven. Three of the
four (Lys593, Phe573, Ile591) have their closest approach to CDRH3 specifically.
**The hotspot set is correctly scoped for a CDRH3 redesign.**

### What the literature emphasises instead — and why that is a probe artifact

The only experimental residue-level epitope mapping of this complex uses mass
spectrometry, and it emphasises **Trp592 and Phe594**, not your four:

> "Statistically significant changes … in 5 residues on HER2 (H66, Y321, H327,
> **W592, and F594**)" … "peptides with the biggest change in deuteration are
> **587–594** and **590–597** … Additionally, a smaller change in HER2 protection
> after trastuzumab binding can also be observed within peptides 565 and 574,
> corresponding to the modification of **F573**."

- Fojtík L, Kalaninová Z, Fiala J, Halada P, Chmelík J, Man P, Kukačka Z, Novák P.
  *Structural Characterization of Monoclonal Antibodies and Epitope Mapping by FFAP
  Footprinting.* Anal Chem 2024. FFAP radical footprinting + HDX-MS, validated
  against 1N8Z. `AMBC_C5bjvuz7s6SFvT10zgLjEDLS61p` · PMID 38698660 · doi:10.1021/acs.analchem.3c04161

**Do not let this unsettle the design choice.** Computed from your own file:
Trp592 is **6.99 Å** and Phe594 is **6.68 Å** from the nearest Fv atom. Neither is
a contact. The explanation is the probe chemistry: **FFAP labels only aromatic and
Cys residues** (Trp/Tyr/Phe/His/Cys). Ile591 and Lys593 — the actual contacts — are
chemically invisible to it, so the method reported the flanking aromatics instead.
And the HDX peptides it names, 587–594 and 590–597, *contain* Ile591 and Lys593.
So the W592/F594 emphasis is the footprint of your own hotspot pair, read out at
peptide resolution through the nearest labelable side chains. It **corroborates**
the set rather than contradicting it. The same paper independently confirms Phe573.

### One real caveat on Lys593

Per-residue MM-GBSA decomposition found **Lys593 to be the least energetically
favourable residue in the whole complex** ("the least favourable residue in the
complex was Lys1027, the starting residue of the third binding loop" — 1027 − 434
= 593). Structurally it is the tightest contact; energetically a buried charged
Lys can carry a desolvation penalty that outweighs its H-bond.
`AMBC_OCRs4Bt7QBxPiVt4ukowLu2fvEA` · doi:10.1096/fasebj.31.1_supplement.775.23

**Implication:** keep Lys593 as a geometric anchor, but do **not** claim in the
pitch that it is the energetic driver of binding. The one published energy
decomposition says the opposite.

### What the set omits — and why the omission is correct

Loop 1 (Pro557/Glu558/Asp560/Gln561) and the Asp596–Pro605 stretch are absent from
the hotspot set. But the computed contact map shows those are engaged by **CDR-H2
(Arg50, Tyr52, Tyr57, Thr58, Arg59) and the light chain (Asp28, Asn30, Thr31,
Ser50, Phe53)** — not by CDRH3. In a CDRH3-only redesign they are not CDRH3's
business. Omitting them is correct by construction.

---

## 3. Mutational / alanine-scanning work on the trastuzumab epitope

**The key negative finding: there is no published alanine scan of the HER2 side of
the trastuzumab epitope in BiomedCore.** Epitope evidence is structural (Cho 2003),
MS-footprinting (Fojtík 2024) and computational (MM-GBSA). If the brief implies
"mutationally validated hotspots", that claim is not available.

What does exist:

**Antibody (paratope) side**
- **hu4D5 = trastuzumab.** Selected CDR positions substituted with all 20 amino
  acids in phage-displayed Fab libraries; the variability map is explicitly
  compared against "a map derived from previous alanine-scanning of the antibody".
  One single point mutation (Trp for Asp, at the most plastic position by
  Wu-Kabat coefficient) improved affinity ~3-fold.
  `AMBC_MHGn8BBeOsu3pBUwK0R7QwVdBJS` · PMID 12206766 · doi:10.1016/s0022-2836(02)00677-0
- **Careful — a common miscitation.** The comprehensive shotgun alanine/homolog
  scan of "an anti-ErbB2 antibody" is **Fab2C4, i.e. pertuzumab, not trastuzumab**.
  It reports a 610 Å² alanine-scan functional epitope shrinking to 369 Å² under
  homolog scanning, with most binding side chains in the heavy chain. Do not cite
  it as trastuzumab data.
  `AMBC_7vdLCwbOWsiqFHJ5dNI7tGd4SqK` · PMID 12079396 · doi:10.1016/S0022-2836(02)00264-4
- Herceptin variant bH1 ("two-in-one", HER2 + VEGF): crystallography plus
  mutagenesis showed distinct amino acids engage HER2 vs VEGF from a heavily
  overlapping surface. `AMBC_KUDgjIQIhLkWjHJImbvO7PsIQda` · PMID 19299620 · doi:10.1126/science.1165480
  Thermodynamic follow-up: `AMBC_FUPqy7uANxF0HHPcl9laZRH7aTW` · PMID 21526167 · doi:10.1371/journal.pone.0017887
- X-ray structures of three humanised 4D5 variants; notes CDR-H3 is the variable
  loop, "differing by up to 8 Å between molecules within an asymmetric unit".
  Directly relevant to how much CDRH3 conformational latitude to expect.
  `AMBC_6ozdxIDDNijfNoA8zrLh0filWP2` · PMID 8095303 · doi:10.1006/jmbi.1993.1099
- Alanine scan of an *anti-trastuzumab* dAb mask (PDB 7PKL) identified Tyr61 and
  Trp104 of the mask as the principal anchors against trastuzumab's CDRs; W104F
  gave a ~22-fold masked/unmasked window. Paratope-side, but a clean worked example
  of structure-guided interface tuning.
  `AMBC_3vl8zKeDUQa01XUiJPxtVpv9pC3` · PMID 42968171 · doi:10.1186/s13036-026-00746-6

**HER2 (antigen) side**
- Sequencing of the **ERBB2 juxtamembrane domain (the trastuzumab binding site)**
  in 54 ERBB2-overexpressing invasive breast cancers: **1 case (2%)** carried a
  missense mutation in exon 17 (reported as "H559A"). Conclusion: epitope mutation
  is a rare event and does not explain clinical trastuzumab resistance. Note the
  structure has ALA at 559, so that mutation label needs checking against whichever
  numbering the paper used before quoting it.
  `AMBC_YiWebLytBGYVokfSZBDGZEsF27j` · PMID 21499233 · doi:10.1038/modpathol.2011.64
- Saturation mutagenesis across 111,176 tumours of the HER2 **transmembrane and
  juxtamembrane** domains found activating G660D, R678Q, E693K, Q709L — all
  C-terminal to the ECD epitope, and still blocked by anti-HER2 antibodies.
  `AMBC_Nrrnw4mNNQYmyA9fR4jcGCTHjda` · PMID 30449325 · doi:10.1016/j.ccell.2018.09.010
- Exon-skipping survey: 29% of 824 HER2-high solid tumours showed domain IV exon
  skipping (exons 13–29), which would delete parts of the trastuzumab epitope.
  `AMBC_9VqHklvUK8nZKzrX5qBuPiC2sne` · doi:10.1200/jco.2024.42.16_suppl.1033 (ASCO abstract, no PMID)

---

## 4. Which CDRH3 residues of trastuzumab contact HER2

**Computed (4.5 Å): only four CDR-H3 residues touch HER2.** (Distances corrected 2026-11-14; an earlier draft listed uncomputed values for Trp99/Phe104/Tyr105.)

| CDR-H3 residue (chain B, sequential) | HER2 partners | closest |
| --- | --- | --- |
| **Trp99** | Pro572, Phe573 | 4.05 Å |
| **Gly103** | Asp570, Pro579, Ile591, **Lys593** | **2.42 Å** (O···Lys593 NZ) |
| **Phe104** | Lys593 | 4.27 Å |
| **Tyr105** | Asp570, Pro571, Pro572, Phe573 | 3.32 Å |

Trastuzumab CDR-H3 is `WGGDGFYAMDY`; the contacting subset is W99 … G103-F104-Y105.
Gly100, Gly101, Asp102, Ala106, Met107, Asp108, Tyr109 make no 4.5 Å contact.

**Literature corroboration, and it is good:**

- **The minimal motif is published.** CDR-H3-derived peptide **`GDGFYA`** (= G101
  D102 G103 F104 Y105 A106), grafted with a CDR-L3 peptide onto a fibronectin type
  III scaffold, binds HER2 at **KD = 58 nM**. That hexapeptide contains exactly the
  three contact residues G103/F104/Y105.
  `AMBC_2ZqITuPRlMpysYNr2IHgJsphec1` · PMID 32975036 · doi:10.1002/biot.202000078
- **CDR-H3 affinity maturation of hu4D5.** Randomising CDR-H3/L3 and panning gave
  variant AH06 with 7.4-fold better affinity (KD 60 pM) and 7.2-fold better
  anti-proliferative activity. Modelling attributed this to CDR-H3 **W98, F100c,
  A101, L102** (Kabat), with "direct hydrophobic interactions between the aromatic
  ring of W98 and the aliphatic group of I613 within HER2 domain IV". Kabat W98
  corresponds to our sequential Trp99 — the same residue our contact map flags.
  `AMBC_TJrq1YVG3lzNqIHRh9sjhfQX2rs` · PMID 26743905 · doi:10.14348/molcells.2016.2235
- FFAP footprinting found trastuzumab heavy-chain **Tyr109** significantly
  protected on complex formation, along with Phe27/Tyr33 (CDR-H1) and
  Trp47/Tyr52/Tyr57 (CDR-H2). It did not flag Trp99 or Phe104 — but again, note the
  method's aromatic-only chemistry and incomplete peptide coverage before reading
  that as absence of contact. PMID 38698660 (as above).
- ML model over 193 single-mutant affinity measurements of trastuzumab concluded
  that "secondary structure induced affinity change due to mutations in the CDR-H3
  is a novel insight", i.e. CDR-H3 mutations act partly through conformation rather
  than direct contact. Relevant to why a sequence-only novelty filter is a poor
  proxy for retained binding.
  `AMBC_MR053BqgDwChZk2sMpMs9p875oL` · PMID 34392800 · doi:10.1080/07391102.2021.1961866

---

## 5. Consequences for the batch

1. **The surface patch is right.** Do not re-run the batch on epitope grounds. The
   four hotspot atoms are real contacts and are the CDRH3-proximal subset of the
   interface. This closes open question 1 in `HANDOFF.md:111` and `STRATEGY.md:262`.
2. **The real fragility is on the antibody side, not the antigen side.** CDRH3's
   entire contribution to this interface runs through **four** residues — Trp99,
   Gly103, Phe104, Tyr105 — and `GDGFYA` is a published minimal binding motif
   (58 nM on a foreign scaffold). The question worth asking of the 100 designs is
   not "did we aim at the right patch" but **"how many designs preserve Trp99,
   Gly103, Phe104 and Tyr105?"** That is checkable from the design sequences
   immediately, and `PLAN.md:177` already independently identified Tyr105 as the
   dominant hotspot (0.866) — consistent with this contact map.
3. **Gly103 is the load-bearing residue and it is a glycine.** It makes the 2.42 Å
   H-bond to Lys593 *through its backbone carbonyl*, and also contacts Asp570,
   Pro579 and Ile591 — four of the seven CDRH3-proximal HER2 residues. A backbone
   contact from a Gly is a conformational constraint, not a side-chain one, so any
   substitution at 103 risks the contact via loop geometry rather than chemistry.
   Novelty pressure that mutates Gly103 is the highest-risk edit in the set.
4. **Re-label two things before the pitch.** `DECK.md:160` already carries the
   Leu573→Phe573 / Val591→Ile591 correction; that correction is confirmed here
   (573 = PHE, 591 = ILE in the file). Also drop any wording that implies the
   hotspots were mutationally validated — they were not (§3).
5. **Numbering discipline.** Three conventions are live in this literature:
   1N8Z/mature (what your file and the FFAP paper use), precursor/UniProt (+22,
   for constructs quoted as "HER2 [23–652]"), and ad-hoc simulation renumbering
   (the +434 offset in the FASEB abstract). Any residue number quoted in the
   submission should name its convention.

## 6. Provenance

Geometry: `inputs/12-her2/her2_fv_1n8z.cif`, heavy-atom distances, 4.5 Å cutoff for
the Fv interface and 5.0 Å for the CDRH3 shell; CDRH3 taken as chain B 99–109.
Script: `/tmp/claude-1001/-home-temp-user-antibody/77b38ae8-7662-45e5-9a48-5a5cf48255cc/scratchpad/contacts.py`
(scratchpad — copy it into `pipeline/` if this needs to be reproducible for judges).
Literature: Amass BiomedCore, 2026-11-14. Two cited records are conference
abstracts with no PMID and are cited by DOI; Cho 2003 and the FFAP paper differ in
that only the latter has fulltext in Amass.
