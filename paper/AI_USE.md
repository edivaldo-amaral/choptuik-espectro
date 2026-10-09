# Use of AI assistants: dated account

Supplementary to Section 7.6 of the paper. The model identifiers are those recorded in the session logs of the
tools (Claude Code project logs and Codex session files); dates are in 2026.

| dates | tool | models | what was done |
|---|---|---|---|
| 10 Sep | OpenAI Codex agent | `gpt-5.6-luna`, `gpt-5.6-sol` | first exploratory session, from the author's request to study the unstable spectrum of the solution of Reiterer and Trubowitz |
| 11–23 Sep | Claude Code (Anthropic) | `claude-opus-5` | main assistant: formulation, code, certificates, documentation |
| 14–16 Sep | OpenAI Codex agent | `gpt-5.6-sol`, `gpt-6-astra`, `gpt-reserve` | work continued while the Claude quota was exhausted: notes and scripts on the formulation, the reduction of sector B, the reconstruction and the local gauge analysis; audited on 16 Sep (`docs/AUDITORIA_MATEMATICA_16SET.md`, `docs/AUDITORIA_FISICA_16SET.md`) and later re-derived or checked |
| 16–17 Sep | OpenAI Codex agent | `gpt-6-astra` | delegated tasks whose output is checked mechanically (`.codex-runs/2026-09-16-setorB-6x18/`) |
| 23 Sep – 8 Oct | Claude Code | `claude-opus-5-5` | main assistant: the counting certificates, the lemmas of Section 4, the reproduction package, the manuscript |
| 24 Sep | OpenAI Codex agent | `gpt-5.6-luna` | one check of the tiles, which turned out to be superficial and is not counted (`docs/AUDITORIA_S3A_TILES_24SET.md`) |
| 24 Sep – 9 Oct | Claude agents (clean context) | `claude-fable-5-1` | the twelve independent checks of Section 7.5 (`.codex-runs/*revisao*`), the last two on the whole manuscript (8 Oct) and on the certificate at s\* (9 Oct) |
| 8–9 Oct | ChatGPT (OpenAI) | Sol 6.1 | readings of the draft by the author's request; their comments were checked against the certificates and sources before being incorporated |

The author posed the problem, directed the work, chose among the approaches, ran the computations, and checked the
results.
