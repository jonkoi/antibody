# Enumeration track — raw artefacts

Candidates produced by enumerating substitutions on the six CDRH3 positions that
make no HER2 contact, rather than by generative design. Five survived the team's
own ranking and were validated here.

**ID mapping.** `SUBMISSION.csv`, `DECK.md` and the findings use `enum_N`. The raw
tool outputs in this folder and in `../stab-ddg/` retain the original `cw_N`
labels, because those are the identifiers actually submitted to `abguard` and to
StaB-ddG — renaming them would break traceability to the recorded results.

| submission ID | raw tool ID | CDRH3 |
| --- | --- | --- |
| enum_1 | cw_1 | WGHDGPYWMD |
| enum_2 | cw_2 | WGGKGPYWMD |
| enum_3 | cw_3 | WGGMGPYWMD |
| enum_4 | cw_4 | WGGRGPYWMD |
| enum_5 | cw_5 | WGPDGPYWMD |

Contents: `submission.csv` (the five plus a wild-type control),
`all_liabilities.json`, `boltz2/` (template-guided refolds, seed 42, with the
wild-type control folded identically), `queries.json` and `reseed_queries.json`.
