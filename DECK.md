# Team 15 — 4 minutes, 7 slides

Backing detail for Q&A: **[DECK-DETAIL.md](DECK-DETAIL.md)** ·
Submission: **[team15_01.csv](team15_01.csv)** · Evidence: **[findings/](findings/)**

Budget ~35 seconds per slide. One idea each. The bold line is the thing to say;
everything else is there to be pointed at, not read out.

---

## 1 — What we tried first, and it failed

Redesign trastuzumab's CDRH3 tip (10 residues) to bind HER2 better, while
differing from **525,063 known CDRH3s** by ≥3 mutations.

We generated 100 designs with RFD4-Proteina and ran the real novelty gate.

| substitutions from wild-type | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| passed the gate | **0%** | **0%** | **0%** | 33% | 82% | 74% | 93% |
| measured chance of binding | 98% | 91% | 73% | 59% | 42% | 32% | 29% |

> **Everything that could bind was rejected. Everything accepted was unlikely to
> bind. Two designs were rejected for being *identical* to sequences already
> measured.**

*Say: we thought this meant the competition was impossible. It doesn't.*

---

## 2 — Why it failed: only 4 of the 10 residues matter

We computed the actual HER2 contacts from the crystal structure, 1N8Z.

**Of the 10 designable residues, four touch HER2: W99, G103, F104, Y105.**
The other six touch nothing.

- The generator scattered its mutations **across the interface** — so it broke
  binding to earn novelty
- Our designs preserved **at most 2 of 4** contacts. Zero of 71 preserved three
- The measured library densely samples exactly the region we needed

> **The novelty filter isn't the problem. Generating-then-filtering is.**

*Say: the filter rewards mutations the assay has never seen — and the generator
put them where they do damage.*

---

## 3 — The fix: design only where it's free

**Fix the four contact residues. Enumerate the other six.**

- **2,091,995** sequences exist at 3–4 substitutions on non-contact positions
- Our second track did exactly that → five designs at **3 substitutions**
- **All five pass the same gate** that rejected 100% of our conservative
  generative designs

And we reordered the pipeline: **novelty gate *before* folding.** The default
ranking spent **30% of our GPU** on designs the gate later rejected.

> **The conservative region was never empty. It was unreachable by a generator
> that didn't know which positions to leave alone.**

---

## 4 — How we validated: four independent checks that disagree

| check | what it tests | result |
| --- | --- | --- |
| `abguard` + 525,063 refs | is it novel and submittable | **4/4 valid**, 0 hard errors |
| Epitope recall, Boltz-2 refold | does it still hit HER2 | **4/4 hotspots**, every design |
| StaB-ddG | binding energy vs wild-type | all within the assay's own noise |
| LAP liabilities | will it survive manufacture | **cleaner than trastuzumab** (5–10 vs 8) |

Each run against a **wild-type control put through the identical pipeline**.

> **The checks ranked the designs differently. So we submit one design leading on
> each, instead of pretending one metric wins.**

*Say: and fold confidence — ipTM — gave its lowest score to the design with the
best geometry. We stopped using it.*

---

## 5 — Literature as an instrument, not a footnote

We used **AMASS** (43M papers + 1.2M trials, dual-corpus) as a *measuring device*
in the loop. It changed the design strategy and **caught three of our errors.**

**It gave us the target.** Querying the trastuzumab–HER2 epitope, then computing
contacts on 1N8Z, is what produced **the four-residue answer on slide 2** — and
found that **Gly103's 2.42 Å hydrogen bond to Lys593 is the tightest contact in
the interface**, made by a *backbone carbonyl*. That is the single fact our whole
approach turned on, and it is why our alanine scan had been blind to it.

**It stopped three wrong claims reaching this talk:**

| we were about to say | it found |
| --- | --- |
| nanobody NB46 "has Tm 90 °C" | that figure is retention during radiolabeling — **no Tm for NB46 exists** |
| cite a nanobody Tm predictor as our number | all three are nanobody-specific, **none transfers** to our scaffold |
| thermostability implies better binding | the two are **explicitly uncorrelated** (PMID 39593207) |

**It priced our biggest risk.** Our top designs share `F104P` + `A106W`. AMASS
found **no** published data for either substitution — and documented CDR3
tryptophan as an oxidation and aggregation liability (PMID 30608647, 31678251).
So we labelled it a correlated bet instead of a result.

> **Every claim carries an amassId plus PMID/DOI. The useful part wasn't the
> citations — it was being told we were wrong three times before a judge did.**

*Say: this is also the one scoring axis we could not compute, so literature is
all we have — and we say so rather than substituting a proxy.*

---

## 6 — The designs, ranked

The rules changed: submit everything, ranked; the lab tests the top 4–14. **So the
ranking is the deliverable.** We submit **24, in two tiers — by strength of
evidence, not by score.**

| | ranks | what backs them |
| --- | --- | --- |
| **Tier 1** | **1–8** | measured on this complex: binding ΔΔG, Boltz-2 refold vs a wild-type control, liabilities, novelty |
| **Tier 2** | 9–24 | the FACS model only — 524,346 **measured** variants, test AUC 0.984 — plus novelty and liabilities |

**All 24: `abguard` valid, 0 hard errors, all novelty-legal at distance 3.**

| rank | CDRH3 | leads on |
| --- | --- | --- |
| 1 | `WGGRGPYWMD` | best measured ΔΔG, 0.20 — indistinguishable from wild-type |
| 2 | `SLGSADYAFD` | best interface H-bond; fewest liabilities; **fails independently** |
| 3–8 | — | remaining measured designs, incl. 2 framework thermostability tests |
| 9–24 | `WYGDGFYNML` … | **13 of 16 preserve all four HER2 contacts** |

**Two honest things about this ranking.**

Tier 1 is first because its **evidence** is stronger, not its scores — four
independent checks beats one model, even when the model prefers the other design.

And **the model saturates**: 24 candidates score exactly 1.000, and the top
scorers collapsed into one motif — 12 of our first 16 were `WYGDGFY**`
near-duplicates. Handing a lab 16 variants of one sequence wastes its slots. So
where the score can't discriminate, **we select for diversity instead** —
minimum 3 mutations apart.

> **A ranking we can't justify is worse than a short list. Every row says which
> tier it's in and what evidence is behind it.**

---

## 7 — What we got wrong

Nine claims we inherited or made ourselves, and disproved by checking:

- a dataset of **"454 measured non-binders" that does not exist** — we were
  calibrating against it
- **`abguard` doesn't check novelty by default.** Silently.
- two of four epitope residues **misnamed in our own notes** since commit one
- our alanine scan **missed the interface's tightest contact** — it's made by a
  backbone carbonyl, and alanine keeps it
- **"conservative designs are impossible"** — our own slide 1 conclusion, killed
  by our own slide 3
- and an hour ago, **we reversed this submission on one over-read measurement**,
  then reversed it back

**One axis we could not address at all:** nothing in either stack predicts
thermostability. We argue it from literature and say so.

> **We'd rather show you the nine things we caught than the one number we'd have
> reported if we hadn't.**

---

## Pre-flight, 17:40

- [ ] `team15_01.csv` attached, @Davide Campagna tagged, nothing else in the message
- [ ] Filename prefix confirmed with the organisers
- [ ] Slide 6 matches the CSV exactly — 4 rows, `fw_WGGMGPYWMD` in slot 3
- [ ] **No thermostability number spoken aloud** — literature rationale only
- [ ] No affinity claim below ~7× (the assay's own replicate spread)
- [ ] Chinery et al. (doi:10.1101/2024.03.26.586756) credited for the HER2 data
- [ ] Both presenters have read `findings/interface-hbond-check.md` — slide 6 changed late
- [ ] Read aloud once, timed. **If it runs over 4:00, cut slide 4's table to two rows.**
