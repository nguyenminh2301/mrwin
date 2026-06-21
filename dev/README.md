# Development docs (not shipped with the package)

This directory holds the **process and planning** record for `mrwin`: the
work-package specs, roadmaps, the project checkpoint, the release checklist, and
the per-WP design/finding notes. It is excluded from the package build
(`.Rbuildignore`), so it never reaches the installed package — it is for
developers and the project history only.

Durable, user-facing reference material stays under `inst/spec/` and ships with
the package:

- `inst/spec/algorithm-spec.md` — the implementation specification.
- `inst/spec/validation-findings.md` — the external calibration findings
  (the R2 lesson, the Fieller fix, the WP19 grid, the doubly-ranked
  incompatibility).
- `inst/spec/benchmark-results.md` — the curated scaling/benchmark tables.

## Map of this directory

- `work-package-roadmap.md`, `acceleration-roadmap.md`, `phase2-direction.md`,
  `implementation-plan.md` — the programme of work (Phase I WP0–WP12, Phase II
  WP13–WP19).
- `project-checkpoint.md` — running status snapshot.
- `release-checklist.md` — pre-release QA gate.
- `wp1`–`wp8` — Phase I work packages (API, data validation, kernel, estimator,
  bootstrap, adjustment, SDPD, simulation engine).
- `wp13`–`wp19` — Phase II work packages (fast kernel, bootstrap acceleration,
  continuous ISG [negative finding], doubly-ranked strata [negative finding],
  analytic variance, theory/manuscript, scalability validation).

The reproducibility scripts behind the findings live in
`tools/validation-scripts/`; the manuscript draft is kept out of the repo
(`manuscript/`, gitignored).
