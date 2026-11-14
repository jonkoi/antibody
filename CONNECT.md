# Connecting the sources

Secrets live in `.env` (gitignored). Nothing secret belongs in a committed file.

## Load the keys, then start Claude Code

`.mcp.json` references `${AMASS_API_KEY}` rather than holding the key, so the key
has to be in the environment before Claude Code starts:

```bash
cd ~/antibody
set -a && . ./.env && set +a     # load keys into the environment
claude                            # restart so .mcp.json is picked up
```

MCP servers are read at startup. Editing `.mcp.json` mid-session does nothing
until you restart, and Claude Code will ask you to approve the server the first
time it sees it.

Verify with `/mcp` inside Claude Code, or `claude mcp list` from the shell.

## Status

| # | Source | Status | What is still needed |
| --- | --- | --- | --- |
| 1 | AMASS | **configured** — key verified against the live endpoint | restart to load it |
| 2 | JapanFold | not configured | a key from `japanfold.aiand.com/account/token` |
| 3 | Nebius Scientific AI | not configured | `credentials.private.json` from `participant-12.zip` |
| 4 | Apheris Foundry | **working** | nothing — local Docker on the VM, no auth |

### 1. AMASS — done

Endpoint `https://mcp.amass.tech/mcp`, streamable-http, bearer auth. Found via
`https://amass.tech/.well-known/ard.json`; note the documentation page at
`https://amass.tech/mcp` is **not** the transport endpoint and returns 405 to a
POST. Verified: `initialize` returns `serverInfo: Amass MCP 0.0.1` with six
cross-linked Cores (BiomedCore, TrialCore, DrugCore, GeneCore, PatentCore,
RegulatoryCore).

[summery.md](summery.md) is the AMASS evidence-scout skill, which refuses to run
without the four `*_biomedcore_*` / `*_trialcore_*` tools. Once the restart is
done, it is what builds the pitch's evidence section.

### 2. JapanFold — key missing

Get a token, then add it to `.env` as `JAPANFOLD_API_KEY` and add a server block
to `.mcp.json`. **Lowest priority** — five folding models already run on the VM.

### 3. Nebius Scientific AI — credentials missing

Needs `credentials.private.json` (gitignored) plus `forge.py`. Neither is on this
box or the VM, so they are on a teammate's machine or still in
`participant-12.zip`. **Chase this one** — Nebius `stab-ddg` is how we score
thermostability, which is the strategic differentiator.

### 4. Apheris Foundry — working

No credentials by design. Every model is a local Docker image on the VM; there is
no Foundry Hub, so `foundry-toolkit module list` fails asking for credentials and
that failure is expected and irrelevant. Use `docker run` directly.

```bash
ssh -i ~/.ssh/team-15 hackathon@89.169.96.144
```

## Where secrets live, and re-issue points

| Secret | Location | Re-issue |
| --- | --- | --- |
| AMASS API key | `.env` → `AMASS_API_KEY` | `platform.amass.tech` |
| JapanFold key | `.env` → `JAPANFOLD_API_KEY` | `japanfold.aiand.com/account/token` |
| Nebius creds | `credentials.private.json` | original `participant-12.zip` |
| VM SSH key | `~/.ssh/team-15` (mode 600) | organisers' onboarding mail |
| AWS creds | VM only, at `~/.aws/credentials` | organisers |

Treat the AMASS key and the VM SSH key as **exposed** — both were pasted into a
chat transcript. The VM dies tonight; rotate the AMASS key after the hackathon.
