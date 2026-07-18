# `mirna_hallmark` — agent onboarding (pointer)

> **This is a SEPARATE subproject** inside the APM repo, with its **own** docs, analysis catalog, discovery
> registry, and run ledger under `mirna_hallmark/` — do **not** mix them with the main `analysis/` ones.

**There is ONE onboarding doc, and it is [`CLAUDE.md`](CLAUDE.md).** Read it before editing code or answering
questions about `mirna_hallmark/`. It is the canonical source for:

- **orientation** → it points you first to `docs/ARCHITECTURE.md` (the axes×models×analyses×results map;
  interactive view `docs/derived/architecture.html`) and `docs/STATE_OF_PLAY.md` (per-axis verdicts);
- the **reuse contract** (import from the parent repo, don't re-implement) and **essential commands / spine**;
- **⛔ the MANDATORY rules** — the documentation one-home protocol and the six working axioms. These are
  load-bearing; skipping them is how findings rot.

*Why this file exists:* Claude Code auto-loads `CLAUDE.md`; other tools (Cursor, …) read `AGENTS.md` and would
otherwise miss it. This pointer sends them there. All content lives in `CLAUDE.md` and the docs it links — this
file is deliberately thin so the two never drift.

*Parent conventions:* cohort / join rules in `analysis/docs/COHORT_AND_JOIN_CONVENTIONS.md`; WSL-from-Windows in
root `.cursor/rules/wsl-shell-bridge.mdc`.
