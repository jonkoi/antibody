# Team 15 — 4 minutes, 6 slides

Backing detail for Q&A: **[DECK-DETAIL.md](DECK-DETAIL.md)** ·
Submission: **[team15_01.csv](team15_01.csv)** · Evidence: **[findings/](findings/)**

Budget ~40 seconds per slide. One idea each. The bold line is the thing to say;
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

## 5 — The four designs

| # | design | CDRH3 | what it leads on |
| --- | --- | --- | --- |
| 1 | `her2_cdrh3_61` | `SLGSADYAFD` | best interface H-bond, fewest liabilities, **fails independently of 2–4** |
| 2 | `enum_4` | `WGGRGPYWMD` | best predicted binding energy |
| 3 | `fw_WGGMGPYWMD` | `WGGMGPYWMD` | **+ framework S85N/R87T** — a thermostability experiment |
| 4 | `enum_1` | `WGHDGPYWMD` | best sequence-model score |

Slots 2–4 share two substitutions, so they rise and fall together — **slot 1 is
deliberately built to fail for different reasons.** Slot 3 is slot 2's CDRH3 plus
two framework changes, so it's a controlled test, not a fourth guess.

> **We are not claiming we beat trastuzumab. Nothing in the measured data beats
> it by a margin the assay can resolve.**

---

## 6 — What we got wrong

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
- [ ] Slide 5 matches the CSV exactly — 4 rows, `fw_WGGMGPYWMD` in slot 3
- [ ] **No thermostability number spoken aloud** — literature rationale only
- [ ] No affinity claim below ~7× (the assay's own replicate spread)
- [ ] Chinery et al. (doi:10.1101/2024.03.26.586756) credited for the HER2 data
- [ ] Both presenters have read `findings/interface-hbond-check.md` — slide 5 changed late
- [ ] Read aloud once, timed. **If it runs over 4:00, cut slide 4's table to two rows.**
