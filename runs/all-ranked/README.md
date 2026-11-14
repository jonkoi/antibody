> **Reconciled current state (14 Nov 2026):** the remote 24-row priority list is preserved in exact order as the first 24 entries. The combined unique input was 653; fresh live abguard accepted 546 and rejected 107 (369 optimisation, 177 de novo). `team15_01.csv` and `SUBMISSION.csv` contain the same 546 passing candidates. See [process overview](../../findings/team-process-overview.md), [priority crosswalk](../reconciliation/priority24_crosswalk.csv), and [full evidence CSV](../reconciliation/ranked_with_evidence.csv). The earlier 639/532 numbers below describe the pre-sync snapshot only.

# Consolidated ranked antibody designs

User requested all actual generated and shortlisted HER2 designs, including de novo long shots, with numbered ranks and selection/background notes. **The 51,417-sequence enumeration pool is excluded.** Calibration-panel antibodies, known measured binders used only as benchmarks, abguard test probes, and lysozyme design controls are excluded. Two WT-CDRH3 controls are saved separately.

## Files

- `all_candidates_ranked.csv`: every unique candidate, including screen-rejected rows marked as such, with numbered rank, original IDs, full heavy/light sequences, reason picked, background notes, risks and available evidence.
- `valid_candidates_ranked_with_notes.csv`: screen-passing candidates, with the same notes.
- `submission_all_valid_ranked.csv`: screen-passing candidates in ranked order, using the exact five-column abguard submission schema.
- `optimisation_ranked_with_notes.csv`, `de_novo_ranked_with_notes.csv`: separate track rankings with notes.
- `optimisation_submission_ranked.csv`, `de_novo_submission_ranked.csv`: per-track five-column submission tables.
- `controls_not_for_submission.csv`: controls, excluded from candidate count/submission.
- `summary.json`, live validator logs and JSON reports: screening counts and provenance.

## Scope and deduplication

639 unique candidates before screening: 335 paired designs, 304 VHHs. Exact deduplication uses format + full heavy sequence + full light sequence. Framework variants remain separate. Format is not track: VHHs may classify as optimisation, so final track comes from live abguard novelty-report.

Included: 18-design optimisation shortlist; three framework combinations and their matched originals; all 12 one-shot RFD4 VHHs; 85 earlier disulfide-restored VHH variants; original 100-design seed-1 and 100-design seed-7 production batches; pilot designs; the earlier 55-candidate winner-prior shortlist; 100 ML candidates actually evaluated with AF2-IG; all available results from previously submitted HER2 BoltzGen jobs, including uninspected later waves. This collection launches no new design jobs. Many source entries are exact duplicates and preserve all aliases after merging.

The 2,306 broad AF2-IG reserve sequences are not automatically treated as a curated shortlist. The 100 explicit ML rescoring candidates and 55 earlier shortlisted candidates are included. No full enumeration pool has been added.

## Ranking interpretation

This is a transparent **heuristic experimental priority order**, not a calibrated probability of affinity improvement or S-tier success. Use per-track rank when choosing within a track. The global order does not imply comparable measured performance across formats.

Priority: reproducible independent epitope contacts, followed by framework-stabilisation exploratory combinations; partially reproducible de novo epitope overlap; single-run favourable contacts; evaluated ML candidates; untested generated backups; known epitope concerns; screen failures last. All screen-passing long shots remain in the submission table as requested, including known off-site candidates, with concerns in the notes CSV.

Within appropriate paired/trastuzumab constructs, existing SPR/FACS model scores guide order. The model has wide uncertainty and random-split leakage concerns. CDR-only affinity scores are not transferred as predictions to changed frameworks or VHHs; a separately labelled parent score is shown for framework variants. The incumbent FACS log-odds is only a sequence safety score, not affinity. AF2-IG confidence is descriptive rather than treated as proof of binding. VHH Tm estimates only order otherwise untested VHHs; they are not directly comparable to a measured trastuzumab Tm baseline.

Full-sequence contact evidence from saved experimental batches overrides earlier single-run evidence. The untracked six-replica summary has uncertain provenance/possible overlapping samples and is not added to verified counts. For paired designs, strict contact success requires >=85% of crystal epitope residues plus all four core residues. The one-shot VHH overlap criterion is >=50% plus all four core residues; its strict Fab-footprint count is also disclosed. These are modelled residue contacts, not competitive-binding results.

## Screening and limitations

Fresh live abguard validation uses the explicit 525,063-provided-sequence HER2 reference DB on the VM (`/home/hackathon/team15-claude/her2_novelty.duckdb`). The JSON/log outputs are saved, and final submission rows must match the accepted IDs exactly. This confirms this tool/reference combination, not any additional unpublished organiser database check.

Potential N-glycosylation and absent canonical disulfides are flagged, not silently filtered by the local expression heuristic. Framework edits are permitted per the user's organiser clarification. Cell-free expression is confirmed; glycosylation capability/occupancy remains uncertain. All new designs lack experimental affinity, Tm and competitive-binding measurements. NanoMelt licence is CC BY-NC-SA 4.0.

Reproduce inventory and rank with `python3 runs/all-ranked/build_ranking.py --final` once the three live JSON screening reports are present. `collect_previous_jobs.py` only reads previously submitted HER2 jobs. No external competition submission has been made.

## Final export counts

Fresh screen: **532/639 pass, 107 rejected**. Passing tracks: **355 optimisation, 177 de novo**. `antibody_designs_ranked.csv` is the compact, numbered, colleague-facing CSV with reasons/background, risk notes and full sequences for all 532 passing designs. The all-639 audit preserves rejected candidates and their flags; they are excluded from submission files. All final submitted-table sequence/format fields were checked against the fresh accepted JSON, and the exported IDs equal the accepted set exactly. Exact duplicates removed; distinct framework/disulfide variants retained. Candidate numbers AB001–AB532 follow the current rank and can change if ranking changes; sequence_id and original_ids provide stable traceability.

Commit packaging: seed-7 full sequences and IDs are preserved in `source/seed7_candidates.csv`, with original generation metrics. The bulk structure directory/archive and unused historical 40k pass-list copy remain local and are excluded from git. The builder prefers the committed compact sequence table, preserving the same ranking without those bulk files.
