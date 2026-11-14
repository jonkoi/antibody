# Handoff — Team 15

State of the work as of **14 Nov 2026, ~11:00 CET**. Written for whoever picks this
up next: a teammate joining mid-day, or us after the VM is gone.

---

## TL;DR

The pipeline works end to end and is reproducible from this repo. One validation
run (4 designs) is complete and analysed. **No production batch has been run yet** —
that is the next action and it is the critical path.

The hackathon VM **shuts down at the end of Saturday 14 November**. Everything
needed to reproduce the work has already been pulled into this repo; nothing else
on that machine is load-bearing except the Docker images, which cannot be moved.

---

## Deadlines

| What | When |
| --- | --- |
| Design submission + pitch | **18:00 today** — *assumed, not confirmed with organisers. Verify this.* |
| VM shutdown | End of Saturday 14 November |
| Experimental results (Adaptyv) | Weeks later, out of our hands |

---

## Access and credentials

Four platforms. All four are connected and verified working. **No secret values are
in this repo** — the files below are gitignored and live only on the local machine.

| Platform | How it authenticates | Where the secret lives | Status |
| --- | --- | --- | --- |
| AMASS | MCP server, API key | `.env` → `AMASS_API_KEY` | verified, returns live data |
| JapanFold | MCP server over HTTP | `.env` → `JAPANFOLD_API_KEY` | verified, `/v1/health` ok |
| Nebius Scientific AI | `forge.py` + creds file | `credentials.private.json` | verified, 34/59 models available |
| Apheris Foundry | none — local Docker images | n/a, lives on the VM | verified, full pipeline ran |

The hackathon VM is reachable as the `hackathon` user; its IP and private key came
from the organisers' onboarding mail and are deliberately not recorded here. The key
is at `~/.ssh/team-15` locally. One H100 80GB, ~360GB free disk.

If `.env` or `credentials.private.json` are lost, they are re-issuable: AMASS from
`platform.amass.tech`, JapanFold from `japanfold.aiand.com/account/token`, Nebius
from the original `participant-12.zip`.

---

## What is done

- **Pipeline validated end to end** — design → refold-prepare → Boltz-2 → refold-score,
  ~4 min wall for 4 designs on one H100.
- **All run artefacts rescued from the VM** into `runs/her2-validation-seed1/`.
- **Batch driver written** (`pipeline/run_design.sh`) — generalises the organisers'
  hardcoded 4-design example; weights, design count, seed, sequence randomness and
  folding model are all environment knobs, plus a multi-model sweep.
- **Analysis tooling written** (`pipeline/analyze_designs.py`) — dependency-free mmCIF
  parsing, diffs designs against wild-type, joins to refold metrics.
- **`abguard` located and working** on the VM. It is *not* on `PATH` by default; it
  lives in a prebuilt uv venv and has been symlinked to `~/.local/bin/abguard`.

## What is not done

- **No production batch.** Everything so far is 4 designs at seed 1.
- **`abguard` has never been run on a design.** This is the highest-value cheap test
  outstanding — see Open question 1.
- **No VHH / de novo designs generated.** Only trastuzumab CDRH3 redesign so far.
- **No thermostability scoring.** Nebius `stab-ddg` and `protein-scoring` are available
  and unused.
- **Confidence cut-off not calibrated** against the 454 measured non-binders.

---

## Findings so far

**Confidence metrics are saturated.** The validation designs scored iptm 0.95–0.97 and
pLDDT 96–99. This is not evidence of quality — the refold is template-guided, the
framework and other CDRs are frozen, and only 10 residues ever move, so Boltz-2 is
largely reproducing a complex it was handed. **iptm and pLDDT will not rank a large
pool.** `min_ipae` showed some spread (0.310 vs 0.331) but two points is not a
distribution. Calibrate against the 454 SPR-measured non-binders in
`sigma2-ab-ag-champloo/` — 21 of them score ipTM ≥ 0.7.

**Conserved W99 tracks selection.** Designs that lost the Trp at the base of CDRH3
(designs 0, 1) were exactly the two `refold-prepare` did not rank into the top 2.
Possibly a cheap pre-filter before spending fold compute. Untested at n=4.

**design_2 sits on the novelty knife-edge.** Three substitutions is exactly abguard's
floor for the optimisation track, and 70% window identity is exactly the boundary the
task description words as "under 70%". Conservative designs must be checked with
abguard, not by arithmetic.

**The design window and the scored CDRH3 may not be the same residues.** The model
designs heavy chain 99–108; abguard scores novelty on the IMGT-defined CDRH3. If IMGT
extends beyond that window, designs can be scored on residues the model was never
allowed to change.

---

## Open questions

Three need a structural biologist. Full write-up with sequences and figures:
<https://claude.ai/artifact/K1Kd6wCkfJ2v9yUYExp6M1>

1. **Are the four hotspot atoms the right epitope definition?** The run aimed at HER2
   Lys593 NZ, Leu573 CD1, Val591 CG2, Pro572 CD (from 1N8Z). That set is the *entire*
   operational definition of "target the trastuzumab epitope". Everything downstream
   inherits it. **This blocks scaling to a real batch.**
2. **Did the designed loop stay on the epitope, or drift?** A folding model will return
   a confident complex with the paratope on the wrong surface patch, and ipTM cannot
   tell. Refolded structures are in `runs/her2-validation-seed1/refold_boltz2/`.
3. **Is design_2's double proline (GDG → PPY) clever or a liability?** Rigidifying the
   loop could help affinity and thermostability, or distort the paratope.

---

## Strategy, in one paragraph

Five prize pools: optimisation (€800/500/200), de novo (€800/500/200), best pitch
($500), AMASS best-use, and pitch placement — which converts 4 design slots into up to
16 and so multiplies the two experimental pools. **Three of the five are decided today
by the pitch, where Judgment alone is 40%.** S-tier requires improving affinity *and*
epitope preservation *and* thermostability simultaneously. Most teams will chase
affinity, where the 524k-variant FACS data is. **Thermostability is the bottleneck and
therefore the differentiator** — the HER2 weights only redesign CDRH3, CDRH3 edits
rarely improve Tm, and no thermostability predictor ships in the Apheris stack. Two
levers: format (VHH scaffolds are intrinsically more thermostable, and express far
better in Adaptyv's cell-free system — expression is an unscored gate that zeroes you
out), and explicit ddG scoring via Nebius `stab-ddg`. Supporting evidence from AMASS:
nanobody NB46 competes with trastuzumab for HER2 *and* retains affinity to 90 °C, so
two of the three axes are demonstrably satisfiable in VHH format at once.

---

## Next actions, in order

1. **Verify the 18:00 deadline** with the organisers.
2. **Answer open question 1** (epitope hotspots) — blocks everything downstream.
3. **Run `abguard` on the two validation designs.** Needs a table with columns
   `sequence_id, kind, heavy_sequence, light_sequence, design_description`. Both chains
   are in `inputs/12-her2/her2_fv_1n8z.cif` (A = light, B = heavy). Settles which
   novelty rule actually binds and whether the IMGT window differs. Cheap, high value.
4. **Launch the production batch** — `N=100 TOPN=20 bash pipeline/run_design.sh all`.
   Long-running; start it before the pitch work, not after.
5. **Score thermostability** on survivors with Nebius `stab-ddg`.
6. **Calibrate the confidence cut-off** against the 454 non-binders.
7. **Assemble up to 20 ranked designs**, split across whatever tracks abguard assigns.

---

## Gotchas

- **There is no Foundry Hub on the VM.** `foundry-toolkit module list` fails asking for
  credentials — this is expected and irrelevant. Use `docker run` directly; all images
  are already local. Documented in `~/rfd4/18-vm.md`.
- **`ssh host "cmd"` does not source `.bashrc`**, so `claude`, `uv` and others appear
  missing. Use `ssh host 'bash -lc "cmd"'`.
- **`abguard` ships its own reference** of 610,109 CDRH3s, used by default. Do not build
  one. `--no-novelty` is required to skip the check, so skipping is always deliberate.
- **abguard assigns the track, you do not.** ≥70% identity to nearest known antibody →
  optimisation (needs ≥3 CDRH3 substitutions); <70% → de novo (needs ≥2).
- **RFD4 HER2 weights only redesign CDRH3.** To touch other CDRs or framework, use
  `antibody-2026-10-30`.
- **The `pass` column is always empty** for the HER2 example — the scorer cannot measure
  absolute movement. Rank on `min_ipae`, `iptm`, `binder_plddt` instead.
- **RFD4-Proteina rejects pasted PDB** — mmCIF only.
- **Known binders fail the novelty filter.** Filter before ranking, never after.
