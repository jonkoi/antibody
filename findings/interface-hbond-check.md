# The contact that decides it — and two of our claims it overturns

Prompted by a reviewer's objection that our phi(104) result "says little about
binding", with the specific suggestion to check **Gly103 O → Lys593 NZ, 2.42 Å in
the crystal**: the tightest contact in the whole trastuzumab–HER2 interface.

They were right on both counts. Measuring it **inverted our submission**.

---

## The measurement

One Boltz-2 model per design, HER2 template, seed 42, with a **wild-type control
folded identically** so the model's own error is visible.

| design | res 103 | **O → Lys593 NZ** | Tyr105 → Phe573 | W99 | verdict |
| --- | --- | ---: | ---: | --- | --- |
| **her2_cdrh3_61** | ALA | **2.51 Å** | 3.22 Å | Ser | **best — tighter than the control** |
| *WT control* | GLY | *2.63 Å* | *3.39 Å* | *Trp* | *reference (crystal: 2.42 Å)* |
| cw_5 | GLY | 3.20 Å | 3.03 Å | Trp | holds |
| cw_4 | GLY | 3.48 Å | 3.30 Å | Trp | holds |
| cw_1 | GLY | 4.07 Å | 2.72 Å | Trp | marginal |
| her2_cdrh3_60 | TYR | 6.48 Å | 3.01 Å | Trp | **lost** |
| **cw_2** | GLY | **6.74 Å** | **11.11 Å** | Trp | **lost — paratope off the epitope** |
| her2_cdrh3_71 | TYR | 6.80 Å | 8.84 Å | Ser | **lost** |
| **cw_3** | GLY | **7.38 Å** | 7.83 Å | Trp | **lost** |
| her2_cdrh3_88 | GLY | 9.08 Å | 6.59 Å | Ser | **lost** |

The control folds to 2.63 Å against the crystal's 2.42 Å, so the model carries
roughly 0.2 Å of error on a case we know the answer to. Differences of 4–9 Å are
far outside that.

---

## Correction 1 — "no loop drift" was an artefact of a loose metric

We reported CDRH3 staying **2.2–3.0 Å** from the four epitope hotspots across all
five enumerated designs and concluded there was no drift. That metric was
**minimum distance from *any* CDRH3 atom to *any* atom of the four hotspot
residues** — far too permissive. A loop that has rearranged completely can still
put *something* near the epitope.

On the specific contact, **cw_2 and cw_3 have lost it entirely** (6.74 Å and
7.38 Å), and cw_2's Tyr105–Phe573 distance of **11.11 Å** says the paratope has
moved off the surface patch altogether. Those were our top two picks on ΔΔG.

A reviewer flagged that they had taken our 2.2–3.0 Å figures on trust without
recomputing them. That trust was misplaced, and the fault is ours.

## Correction 2 — the phi(104) result proved much less than we said

We presented phi(104) landing at −59° to −78° across all five as evidence that
`F104P` is geometrically fine. **A co-folding model will produce
proline-compatible backbone angles for any sequence containing a proline**, so
that observation is close to vacuous.

Worse, the same folds show **every enumerated design flips phi(101) negative**
where the wild-type control sits at +78°. The predicted loops are **rearranged**,
not wild-type plus two edits — which is consistent with what the contact distances
then revealed.

---

## What this means, and the irony in it

**The contact is made by a backbone carbonyl, and every residue has one.** That is
why `her2_cdrh3_61` holds the bond at 2.51 Å *despite* carrying Ala at position
103 — G103A keeps the carbonyl and the loop positions it correctly.

So the same property that made our alanine scan blind to G103's importance
(alanine retains the backbone carbonyl) is the property that lets a design mutate
G103 and keep the bond. **Side-chain identity at 103 is not what matters; loop
geometry is.** Designs that "preserve G103" by sequence can still lose the bond
(cw_2, cw_3, both with Gly at 103 and no contact), and designs that mutate it can
keep it (`_61`).

**Sequence-level contact preservation is therefore the wrong filter.** We had been
ranking on it. The right filter is the measured geometry of the refolded complex.

---

## The revised submission

Filtered on contact retention first — a design whose loop sits 7 Å off the
epitope cannot bind whatever its ΔΔG says — then ordered within the survivors.

| # | design | CDRH3 | H-bond | ΔΔG (fold) | liabilities |
| --- | --- | --- | ---: | ---: | ---: |
| 1 | **her2_cdrh3_61** | `SLGSADYAFD` | **2.51 Å** | 0.42 (2.0×) | **5** |
| 2 | **cw_4** | `WGGRGPYWMD` | 3.48 Å | **0.20 (1.4×)** | 8 |
| 3 | **cw_1** | `WGHDGPYWMD` | 4.07 Å | 0.50 (2.3×) | 9 |
| 4 | **cw_5** | `WGPDGPYWMD` | 3.20 Å | 0.74 (3.5×) | 9 |

`abguard team-validate --reference-db` → **4/4 valid**. All four sit inside the
~7× assay noise floor on ΔΔG, so contact retention was the binding constraint,
not energetics.

**Dropped:** cw_2 and cw_3 (contact lost), her2_cdrh3_71 (6.80 Å — it keeps Phe104
but loses the bond, so the F104 hedge we wanted it for is worthless), and
her2_cdrh3_88 and her2_cdrh3_60.

`cw_5` is included despite having been called the weakest of the enumerated set by
both us and the reviewer, on ΔΔG and on a now-withdrawn G101P objection. It holds
the contact at 3.20 Å. `cw_1` carries 9 liabilities and keeps wild-type's `DG`
isomerisation motif, which is the cost of its position here.

---

## What we still cannot say

- **One model, one seed, per design.** This is not a distribution, and Boltz-2 is
  the model whose *confidence* we showed is barely better than chance at
  separating real pairs from shuffled ones (ROC AUC 0.645). We are using its
  coordinates against a matched control rather than its confidence, which is a
  more defensible use — but it is still one model.
- **A refold is not a measurement.** None of these designs has been near an
  instrument. The honest claim is that four designs retain the interface
  hydrogen bond *in silico* under a template-guided co-fold.
- **The reviewer's two surviving objections are untouched by any of this:** F104
  packs against the light chain (Tyr49/Leu46/Tyr55), and `A106W` gives every
  enumerated design a second CDR3 tryptophan, a documented oxidation and
  aggregation liability (PMID 30608647, 31678251). Three of our four carry both.
  `her2_cdrh3_61` carries neither, which is a further reason it leads.
