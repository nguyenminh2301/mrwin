# AGENTS.md

Working standards for **any** AI agent (tool-agnostic) on `mrwin` live in
[`CLAUDE.md`](CLAUDE.md). Read it before doing methodological, statistical, or
kernel work — it is the canonical academic-standards and process document for
this research programme.

Quick orientation:
- `ROADMAP.md` — the programme plan and what is shipped vs planned.
- `dev/README.md`, `dev/project-checkpoint.md` — process and current status.
- `inst/spec/validation-findings.md` — external-calibration lessons (ship with the package).
- `dev/findings-*.md` — banked research findings (incl. negative results).

The one rule everything else elaborates: **no documented claim before its test
passes**, and parity between two of your own implementations is not validation —
close the loop with an external ground truth.
