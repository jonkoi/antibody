# Validating the enumerated design set — and correcting our headline claim

A second track in the team, working from the same data, produced five designs. We validated them
against three gates they had not been through, using methods neither team trained.
**They are better than ours, and in validating them we found that our central
claim was overstated.**

---

## 1. The claim we have to correct

We were going to present this as our headline:

> The novelty filter and binding affinity are structurally opposed. Acceptance is
> **0% at 2, 3 and 4 substitutions**; not one of our 13 conservative designs
> survived.

Those numbers are real, but **the conclusion drawn from them was wrong.** All five
enumeration track designs sit at **3 substitutions and all five pass the novelty gate**
(`abguard team-validate --reference-db`, 525,063 references, minimum distance 3;
the wild-type control failed as designed).

So 0% acceptance at ≤4 substitutions is a property of **our generator's output
distribution**, not of the design space. RFD4 placed its mutations across the
interface, which is exactly the region the 524,346-variant library samples
densely. Enumerating deliberately on the positions that are *not* measured finds
plenty of novel conservative designs.

**Corrected claim, which is more useful:** *generate-then-filter fails here.* The
conservative region is not exhausted — it is unreachable by a generative model
that does not know which positions it must leave alone. The fix is constraint-aware
enumeration, and the enumerated set demonstrates it.

---

## 2. Why their approach works: only 4 of 10 positions touch HER2

Computed from 1N8Z (`findings/epitope-1n8z-verification.md`): the CDRH3 residues
contacting HER2 are **W99, G103, F104, Y105**. The other six — G100, G101, D102,
A106, M107, D108 — make no HER2 contact.

So ≥3 mutations can be placed entirely on non-contact positions. There are
**2,091,995** such sequences at 3–4 substitutions.

| | our 14 (generative) | that track's 5 (enumerated) |
| --- | --- | --- |
| HER2 contacts preserved | **2 of 4** at best | **3 of 4** (all five) |
| substitutions | 5–6 | **3** |
| designs preserving ≥3 contacts, of 71 accepted | **0** | n/a |

Their five keep W99, G103 and Y105 and break only F104 (→Pro). Ours break two or
three contacts each.

**Note on our own ΔΔG scan's blind spot:** it rated G103A at 0.042 kcal/mol,
essentially zero — but G103's backbone carbonyl makes the **2.42 Å hydrogen bond
to Lys593, the tightest contact in the whole interface**. Alanine retains the
backbone carbonyl, so alanine scanning cannot see backbone-mediated contacts. We
would have kept mutating G103 on the strength of that number.

---

## 3. Independent validation — three gates

### Novelty: all five pass
`abguard team-validate --reference-db her2_novelty.duckdb` → 5 of 6 rows valid,
the only failure being the wild-type control. Minimum distance 3 for all five.

### Binding energetics: all five land within the assay noise floor
StaB-ddG, scored jointly as multi-site variants. **Neither team trained this
model**, and their predicted KD came from their own, so this is a genuine
independent check.

| design | ΔΔG kcal/mol | ±SE | fold worse than WT | their pred. KD | their P(beat) |
| --- | ---: | ---: | ---: | ---: | ---: |
| **enum_3** | **0.17** | 0.17 | **1.3×** | 0.89 nM | 0.153 |
| **enum_4** | **0.20** | 0.16 | **1.4×** | 0.95 nM | 0.142 |
| enum_2 | 0.42 | 0.17 | 2.0× | 0.77 nM | 0.179 |
| enum_1 | 0.50 | 0.18 | 2.3× | 0.50 nM | 0.244 |
| enum_5 | 0.74 | 0.18 | 3.5× | 1.52 nM | 0.098 |

Median ΔΔG **0.42** against our **1.38**. The 6.2× noise floor is ΔΔG 1.08, so
**all five are inside it** — where only 3 of our 14 were. **enum_3 and enum_4 are
statistically indistinguishable from wild-type** (0.17 ± 0.17, 0.20 ± 0.16), which
nothing of ours achieved.

The two models disagree on the *within-set* order — StaB-ddG puts enum_3 first,
their model puts enum_1 first — but agree on what matters: all five are good. Trust
the level, not the ranking.

### Developability: slightly worse than ours, and one real flag

| | that track's 5 | our 14 | wild-type |
| --- | ---: | ---: | ---: |
| liabilities | **8–9** | 5–7 | 8 |
| high-severity | 2–3 | 2 | 3 |
| CDR3 tryptophans | **2** (all five) | 0–1 | 1 |
| unpaired Cys / N-glyc sequon | **0 / 0** | 0 / 0 | 0 / 0 |

They pass the hard expression filters, so nothing is disqualifying. But:

- **Every design carries two CDR3 tryptophans** — the retained W99 plus `A106W`.
  A CDR3 Trp is a genuine oxidation and aggregation hotspot.
- **enum_1 and enum_5 carry 3 high-severity liabilities instead of 2**, because they
  keep both D102 and G103, which is wild-type's own `DG` isomerization motif.
  enum_2, enum_3 and enum_4 mutate D102 and lose it. *Preserving the contact
  reintroduces the parent's liability* — a real tension in the strategy.

### The correlated risk
**All five share `F104P` and `A106W`**, and differ at only one other position. So
they are five variations on one design, not five independent bets. F104 is a
**peripheral** to HER2, not a tight contact: its nearest HER2 atom is **4.27 Å**
(an earlier version of this file said 3.42 Å, which is actually the HER2
Phe573–CDRH3 **Tyr105** stack — a different pair). F104's real role is **VH/VL
packing**, against light-chain Tyr49/Leu46/Tyr55. The sharper objection to F104P
is geometric: **backbone phi at F104 is −141° in 1N8Z, and proline requires
about −65°**, so F104P is incompatible with the wild-type loop conformation. **If F104P is wrong,
all five fail together.**

---

## 4. Recommended submission

Take theirs, but hedge the correlated risk with one of ours.

| # | design | CDRH3 | why |
| --- | --- | --- | --- |
| 1 | **enum_3** | WGGMGPYWMD | best independent ΔΔG (0.17, WT-equivalent), 2 high-sev |
| 2 | **enum_2** | WGGKGPYWMD | best enumeration-track model KD among the low-liability three |
| 3 | **enum_4** | WGGRGPYWMD | ΔΔG 0.20, WT-equivalent, 2 high-sev |
| 4 | **her2_cdrh3_61** | SLGSADYAFD | **the hedge** — different mutation pattern, keeps F104 absent, cleanest liabilities (5), ΔΔG 0.42 tied with enum_2 |

Rationale for the hedge: enum_1 would be the fourth enumeration track pick on P(beat), but it
carries 3 high-severity liabilities *and* shares the F104P/A106W pattern. `_61`
costs little in predicted affinity and fails independently of them.

---

## 5. What this does to the pitch

It improves it. The arc becomes: we measured a trade-off, drew too strong a
conclusion from it, and a second method disproved our own conclusion — then both
routes converged on the same mechanism, that only 4 of 10 positions matter and the
design must avoid them. Our contribution is the mechanism and the negative
results; theirs is the design set that follows.

Two of our own claims get corrected on stage, which is the point: "0% at ≤4
substitutions" becomes a statement about *our generator*, and our ΔΔG scan is shown
to be blind to the interface's tightest contact.
