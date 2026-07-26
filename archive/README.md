# archive/ — historical artifacts (not built, not maintained)

This directory preserves earlier work that is **no longer part of the package**
but is kept for provenance and reproducibility of the project history. Nothing
here is shipped (`.Rbuildignore`d) or tested in CI; treat it as read-only.

## `python-prototype/` — the original P1 simulation engine (v5)

The method was first prototyped in Python (`p1_engine_v5/`) together with its
own test suite (`tests/`), packaging (`pyproject.toml`, `requirements.txt`) and
an end-to-end orchestrator (`run_all.sh`). That prototype established the
data-generating process, the hierarchical win/loss kernel, the multiplier
bootstrap, the pleiotropy-bounded CI, and the Q-statistic asymptotics.

It has since been **superseded by the R package** in `R/` + the compiled
`src/` kernel, which is the maintained implementation (subquadratic kernel,
analytic influence-function variance, Fieller intervals, ordinal IPTW, SDPD
diagnostics). The R results were differentially tested against this prototype
during development.

Kept because: it documents how the algorithm was first validated and is a useful
cross-language reference. It is **not** a dependency of the R package and need
not be installed to use `mrwin`.
