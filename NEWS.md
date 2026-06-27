# mrwin (development version)

## New features

* **Two-sample / summary-data win-ratio MR.** `mrwin_win_snp()` estimates a
  per-SNP win-odds coefficient (the log-win-odds slope on allele dosage) with a
  closed-form influence-function standard error; `mrwin_win_gwas()` runs it across
  a SNP panel (a "win-odds GWAS"); and `mrwin_twosample_ivw()` pools the per-SNP
  `(beta_GX, delta_winodds)` pairs by inverse-variance weighting — the win-ratio
  analogue of ordinary two-sample MR. The estimator is consistent for the causal
  win-odds gradient and its standard error is calibrated (SE/SD ≈ 1, 96% CI
  coverage; see `inst/spec/validation-findings.md`). A new unexported C++ helper
  (`mrwin_subject_win_loss_cpp`) supplies the per-subject win/loss counts for the
  subquadratic standard error.

## Validation

* **Paper-01 full Monte-Carlo grid** — type-I, power, weak-instrument, and
  pleiotropy/SDPD, recorded in `inst/spec/validation-findings.md`. Type-I is
  calibrated-to-conservative; weak-instrument coverage stays ~0.96 at every
  instrument strength (the method flags weak instruments and returns unbounded
  Fieller intervals rather than falsely excluding the truth); SDPD is calibrated
  under the null with a clean power curve against InSIDE-satisfying pleiotropy.

## Documentation & structure

Project organization and onboarding for a multi-paper, multi-feature package.

* **Getting-started guide for clinicians & epidemiologists** — a bilingual
  (Tiếng Việt + English) vignette covering installation of R/RStudio and the
  package, a first analysis, and how to read the output, for readers with no
  programming background.
* **Documentation website scaffolding** (`_pkgdown.yml`): the ~40 exported
  functions are organized into themed reference sections, with reserved
  placeholders for planned feature families.
* **`ROADMAP.md`** records shipped features, planned modules, and the papers in
  the programme.
* **Tidier repository**: the original Python prototype moved to
  `archive/python-prototype/`; non-English READMEs to `translations/`; paper
  drafts to private `papers/<slug>/` folders; per-WP dev notes grouped under
  `dev/work-packages/`. The R package core (`R/`, `src/`, `man/`,
  `tests/`) is unchanged.

# mrwin 0.1.0

First versioned release. Adds the scalability and inference work and records two
externally-validated negative findings.

## New features

* **Subquadratic, compiled win/loss kernel** (`backend = "fast"`). The
  hierarchical comparison is recast as a weighted multivariate dominance count
  (Fenwick/BIT for K = 1, 2; CDQ divide-and-conquer for K = 3) and ported to
  C++ via Rcpp. Bit-for-bit identical to the dense backend; brings the estimator
  to biobank scale (e.g. K = 3, N = 80,000 in under a second).
* **Analytic influence-function variance** (`inference = "analytic"`). A
  closed-form `cov_u` reproduces the multiplier bootstrap for the sampling
  component, with an exact Monte-Carlo term for GWAS-weight uncertainty and a
  first-order weighted form under ordinal IPTW.
* **Fieller confidence intervals are now the primary reported interval.**
  External calibration showed the original delta-method ratio interval is
  mis-calibrated (structurally too wide); the Fieller construction on the same
  covariance is correctly calibrated and is reported by
  `print`/`summary`/`tidy`/`mrwin_report`, with an explicit unbounded-interval
  report under a weak instrument. The delta-method interval is retained as a
  labelled reference.
* **Doubly-ranked stratification** (`mrwin_doubly_ranked_strata()`,
  `stratification = "doubly_ranked"`) is implemented and wired end to end.

## Validated negative findings

* **Continuous (bandwidth-controlled) ISG** reduces exactly to the decile
  estimator in the boxcar limit, but is noisier than the decile and does not
  reduce stratum-count sensitivity (the gradient is a ratio whose denominator
  shrinks under finer smoothing). The discrete estimator with a sensitivity
  analysis over the stratum count is recommended.
* **Doubly-ranked stratification is incompatible with the between-stratum
  DS-CWR estimand.** It balances the instrument across strata (its within-stratum
  LACE selling point), which makes the adjacent-stratum contrast
  confounder-driven (type-I error ~1.0 under confounding). `mrwin()` keeps it
  runnable but emits a `doubly_ranked_invalid` warning; `prs_rank` remains the
  default and only validated scheme.

See `inst/spec/validation-findings.md` for the calibration tables and diagnoses.

## Documentation and packaging

* Development and planning docs (work-package specs, roadmaps, checkpoint) moved
  to `dev/` and excluded from the package build; `inst/spec/` now ships only the
  durable references (`algorithm-spec.md`, `validation-findings.md`,
  `benchmark-results.md`).
* Reproducibility scripts collected under `tools/validation-scripts/`.
* `DESCRIPTION` updated (title, description, `URL`, `BugReports`).
