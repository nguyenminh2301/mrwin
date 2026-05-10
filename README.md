# mrwin: Causal Win Statistics for Hierarchical Composite Endpoints

[English](README.md) | [Tiếng Việt](README.vi.md) | [中文](README.zh.md) | [日本語](README.ja.md) | [한국어](README.ko.md) | [Français](README.fr.md) | [Deutsch](README.de.md)

## What This Package Does

`mrwin` is an R package that answers a specific question in epidemiology:

**Does a genetically predicted exposure shift people toward better or worse
overall outcome profiles, when outcomes are clinically ranked by severity?**

In many diseases, patients do not experience a single event. A person with
cardiorenal disease may die, be hospitalized for heart failure, or develop
renal decline. These events are not interchangeable: death is worse than
hospitalization, which is worse than a lab abnormality. Standard analyses
treat all events equally or analyze them one at a time, losing the clinical
ordering that matters most to patients and clinicians.

`mrwin` combines two statistical frameworks to address this:

1. **Win statistics** -- pairwise comparisons that respect clinical severity
   hierarchies (death beats hospitalization beats biomarker decline).
2. **Mendelian randomization** -- genetic instrumental variables that support
   causal claims even when confounders are unmeasured.

The result is a **dose-standardized causal win ratio (DS-CWR)**: a single
number that summarizes whether a genetically predicted exposure makes the
prioritized outcome profile better or worse.

---

## Mendelian Randomization: A Brief Introduction

### The Problem: Unmeasured Confounding

In observational epidemiology, we often want to know whether an exposure
(e.g., LDL cholesterol) causes an outcome (e.g., heart failure). The
difficulty is that people with high LDL also differ from people with low LDL
in many other ways -- diet, exercise, socioeconomic status, other medications.
These are confounders. Even after adjusting for measured confounders,
unmeasured confounders may remain.

### The Solution: Genetics as a Natural Experiment

At conception, each person randomly inherits one allele from each parent at
each genetic locus. This is called **Mendelian segregation**. If a genetic
variant affects an exposure (e.g., a variant in the LDLR gene raises LDL
cholesterol), then people carrying that variant are, on average, exposed to
higher LDL throughout life -- not because of their diet or lifestyle, but
because of their genotype.

This random assignment is the basis of **Mendelian randomization (MR)**. By
using genetic variants as **instrumental variables**, MR estimates the causal
effect of the exposure on the outcome, bypassing unmeasured confounding.

### Key Assumptions of MR

For the genetic instrument to be valid, three conditions must hold:

1. **Relevance**: The genetic variant must be associated with the exposure.
2. **Independence**: The genetic variant must not be associated with
   confounders (this follows from Mendelian segregation in well-designed
   studies).
3. **Exclusion restriction**: The genetic variant must affect the outcome
   only through the exposure, not through other pathways.

### Types of Data in MR

| Data type | What it is | Example |
|---|---|---|
| Individual-level | One row per person with genotype, exposure, outcome | A cohort of 10,000 people with SNP data, LDL measurements, and hospital records |
| Summary-level | Per-SNP effect estimates from GWAS | "SNP rs12345 has effect 0.05 (SE 0.01) on LDL" |
| Polygenic risk score (PRS) | Weighted sum of many SNPs | PRS for LDL cholesterol, computed from 100 SNPs |

`mrwin` uses **individual-level data** with a **PRS instrument**. The PRS
serves as a continuous instrument that captures the cumulative genetic
predisposition toward the exposure.

---

## Win Statistics: The Target Estimand

### What Is a Win Statistic?

A win statistic compares every pair of individuals in a study. For each pair
(A, B), the comparison follows a clinical hierarchy:

1. First, compare on the most important outcome (e.g., death). If A died and
   B did not, B wins.
2. If tied on the first outcome, move to the next (e.g., hospitalization).
3. Continue down the hierarchy until one person wins or all outcomes are tied.

The **win ratio** is the number of wins for the higher-exposed group divided
by the number of wins for the lower-exposed group.

### Why Win Statistics Matter

| Method | Treats events equally? | Respects clinical order? | Handles competing risks? |
|---|---|---|---|
| Time-to-first-event | Yes | No | No |
| Composite endpoint analysis | Yes | No | Partially |
| Win ratio / win statistic | No | Yes | Yes |

Win statistics were introduced by Pocock et al. (2012) and have been adopted
in over 36 randomized clinical trials between 2022 and 2024.

### The Causal Win Ratio (cCWR)

Standard win statistics require randomization or strong ignorability. In
observational epidemiology, this assumption is rarely defensible.

`mrwin` defines the **continuous Causal Win Ratio (cCWR)** as a gradient
estimand across quantiles of the genetic instrument:

- Divide the population into ordered strata by PRS (e.g., deciles).
- Within each adjacent stratum pair, compute the win ratio.
- Pool across strata using instrument-standardized gradients and GLS
  meta-analysis.

The result is a **marginal, population-level causal effect** that does not
require strong ignorability.

---

## Core Formulas: Plain-Language Explanation

### 1. The Hierarchical Comparison Kernel

For two individuals *i* and *j*, compare outcomes from highest to lowest
priority:

- If *j* had the event at priority *k* but *i* did not (or *i* survived
  longer): *i* wins (+1).
- If *i* had the event at priority *k* but *j* did not (or *j* survived
  longer): *i* loses (-1).
- If tied: move to the next priority.

**In words**: "Did this person have a better outcome profile than the other,
respecting the clinical severity order?"

### 2. The Stratum-Specific Win Ratio

Within a stratum (e.g., the 6th PRS decile vs. the 5th):

```
theta_d = (number of wins for higher stratum) / (number of losses)
log_theta_d = log(theta_d)
```

**In words**: "Among people with slightly higher genetic predisposition to
the exposure, did they tend to win more often than they lost?"

### 3. The Instrument-Standardized Gradient (ISG)

```
Delta_X_d = mean(exposure in stratum d) - mean(exposure in stratum d-1)
ISG_d = log_theta_d / Delta_X_d
```

**In words**: "How much does the log win ratio change per unit increase in
the exposure, as predicted by genetics?"

### 4. GLS Pooling

The ISG values across all stratum pairs are pooled using Generalized Least
Squares with a shrinkage covariance estimator:

```
delta_GLS = pooled ISG across all strata
DS-CWR = exp(delta_GLS)
```

**In words**: "What is the overall causal effect of the exposure on the
prioritized outcome profile?"

- **DS-CWR > 1**: The exposure is associated with a better outcome profile.
- **DS-CWR < 1**: The exposure is associated with a worse outcome profile.
- **DS-CWR = 1**: No evidence of a causal effect.

### 5. The Multiplier Bootstrap

Uncertainty is estimated by:
1. Perturbing the GWAS weights (external uncertainty).
2. Drawing random multiplier weights (internal uncertainty).
3. Recomputing the full pipeline for each bootstrap iteration.

This gives a confidence interval that accounts for both sources of
uncertainty.

### 6. The Step-Down Pleiotropy Diagnostic (SDPD)

The SDPD tests whether the genetic instrument affects the highest-priority
outcome through pathways other than the exposure (direct pleiotropy). This
is done by fitting an MR-Egger regression on the per-SNP exposure and
outcome effects.

**In words**: "Is the genetic instrument valid, or does it affect the
outcome through other pathways?"

---

## Clinical Examples

### Example 1: Cardiovascular Disease

**Research question**: Does genetically predicted LDL cholesterol worsen the
overall cardiorenal trajectory?

**Endpoint hierarchy** (highest to lowest priority):
1. Death
2. Heart failure hospitalization
3. Renal decline

**Interpretation**: If DS-CWR = 0.82 (95% CI: 0.70 to 0.96), this means
that per unit increase in genetically predicted LDL, the prioritized
cardiorenal outcome profile worsens by 18%. The CI does not cross 1,
suggesting a statistically significant harmful effect.

**Clinical meaning**: Higher LDL does not just cause heart attacks -- it
shifts the entire disease trajectory toward earlier death, more
hospitalizations, and faster kidney decline, in that clinical order.

### Example 2: Dementia

**Research question**: Does genetically predicted sleep disturbance shift
people toward earlier severe neurocognitive outcomes?

**Endpoint hierarchy**:
1. Death
2. Dementia diagnosis
3. Nursing-home admission

**Interpretation**: If DS-CWR = 0.91 (95% CI: 0.78 to 1.06), the point
estimate suggests harm, but the CI crosses 1. There is no strong evidence
of a causal effect.

**Why this matters**: In standard MR of sleep on dementia, people who die
before developing dementia are lost from the analysis. This is a form of
survivor bias: a harmful sleep variant may appear protective because its
carriers never live long enough to be diagnosed. `mrwin` avoids this by
treating death as the highest-priority event, so differential survival is
absorbed into the estimand rather than conditioned away.

### Example 3: Cancer Epidemiology

**Research question**: Does genetically predicted BMI worsen prioritized
cancer outcomes?

**Endpoint hierarchy**:
1. Cancer death
2. Progression or metastasis
3. Recurrence

**Interpretation**: If DS-CWR = 0.75 (95% CI: 0.60 to 0.94), this means
higher genetically predicted BMI shifts the cancer trajectory toward worse
outcomes at every level of the hierarchy.

**Clinical meaning**: BMI does not just increase cancer risk -- it worsens
the entire disease course from recurrence through death.

---

## Installation

```r
# From GitHub (development version)
# install.packages("devtools")
devtools::install_github("nguyenminh2301/mrwin")
```

---

## Quick Start

### Simulated Data (No Real Cohort Needed)

```r
library(mrwin)

# Simulate a cardiorenal dataset
cfg <- mrwin_config(n_outcome = 500, m_snps = 20, seed = 1)
dat <- mrwin_simulate(cfg)

# Fit the model
fit <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  controls = mrwin_controls(n_strata = 5, bootstrap = 100, seed = 2)
)

# Examine results
print(fit)
summary(fit)
plot(fit)
```

### Reading the Output

```
mrwin fit
  N: 500
  SNPs: 20
  Priorities: 3
  Strata: 5
  Bootstrap: 100 (100 valid)
  Backend: dense
  Adjustment: none
  delta_GLS: -0.2000
  DS-CWR: 0.8187
  95% CI: 0.7000 to 0.9600
  Q: 1.2345 ( df = 2 , p = 0.5394 )
```

| Field | What it means | What to look for |
|---|---|---|
| `delta_GLS` | Log-scale pooled effect | Negative = harmful, positive = beneficial |
| `DS-CWR` | Exponentiated effect (win ratio) | < 1 = harmful, > 1 = beneficial, 1 = no effect |
| `95% CI` | Bootstrap confidence interval | Does it cross 1? |
| `Q` | Heterogeneity statistic | Small p-value = effect varies across strata |
| `Warnings` | Structured caveats | Always check before interpreting |

### Interpreting the Result

If DS-CWR = 0.82 (95% CI: 0.70 to 0.96):

- The exposure is associated with an 18% worse prioritized outcome profile.
- The CI does not cross 1, so the effect is statistically significant.
- Check the Q p-value: if significant, the effect may vary across the
  exposure range.
- Check warnings: if `sdpd_rejected`, pleiotropy may invalidate the
  estimate.

---

## Real Cohort Template

```r
snp_cols <- paste0("rs", 1:20)

fit <- mrwin(
  data = dat,
  endpoint = mrwin_endpoint(
    time = c("t_death", "t_hf", "t_renal"),
    status = c("d_death", "d_hf", "d_renal"),
    priority = c("death", "heart_failure", "renal_decline")
  ),
  genotype = snp_cols,
  exposure = "ldl_cholesterol",
  gwas = mrwin_gwas(beta = beta_hat, se = se_beta, snp = snp_cols),
  covariates = c("age", "sex", "PC1", "PC2"),
  controls = mrwin_controls(
    n_strata = 10,
    bootstrap = 500,
    seed = 20260510,
    adjustment = "ordinal_iptw",
    sdpd_scale = "both"
  )
)

summary(fit)
plot(fit, type = "forest")
mrwin_report(fit, file = "analysis_report.md", format = "markdown")
tidy(fit)
```

---

## Output Formats

### Print

```r
print(fit)
```

Shows the main result: sample size, SNPs, strata, bootstrap validity,
DS-CWR, CI, Q statistic, and warnings.

### Summary

```r
summary(fit)
```

Produces a formatted report with:
- Main DS-CWR estimate, CI, p-value, and optional pleiotropy-bounded CI
- Adjacent stratum gradients (log-theta, Delta-X, ISG, CWR)
- Heterogeneity Q statistic
- Fieller sensitivity CI
- SDPD diagnostics (when enabled)
- Structured warnings

### Tidy (for tables and further analysis)

```r
tidy(fit)
#   term estimate    delta se_delta statistic p_value  ci_low ci_high ...
# DS-CWR    0.82 -0.2000   0.075    -2.667  0.0077  0.7000  0.9600 ...
```

### Report

```r
mrwin_report(fit, format = "markdown")
mrwin_report(fit, file = "report.md", format = "markdown")
```

### Plots

```r
plot(fit, type = "isg")       # Adjacent ISG gradients (default)
plot(fit, type = "forest")    # Forest plot of adjacent CWR
plot(fit, type = "bootstrap") # Bootstrap distribution of log-theta
```

---

## Covariate Adjustment

When confounders are measured (age, sex, principal components), use
ordinal IPTW adjustment:

```r
fit_adj <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  covariates = cbind(u_proxy = scale(dat$U)),
  controls = mrwin_controls(
    n_strata = 5,
    bootstrap = 100,
    seed = 2,
    adjustment = "ordinal_iptw",
    sdpd_scale = "both",
    pleiotropy_bias_radius = 0.05
  )
)
```

The IPTW adjustment improves precision by balancing covariates across PRS
strata. The SDPD tests for direct pleiotropy.

---

## Strengths and Limitations

### Strengths

- **Causal inference under unmeasured confounding**: Uses genetic
  instruments, not observed exposure alone.
- **Respects clinical hierarchy**: Death is prioritized over
  hospitalization, which is prioritized over biomarker decline.
- **Handles competing risks naturally**: Differential survival is absorbed
  into the estimand, not conditioned away.
- **Single summary measure**: One DS-CWR replaces multiple component-wise
  analyses.
- **Diagnostic framework**: SDPD tests the exclusion restriction; Q
  statistic tests heterogeneity; structured warnings flag problems.
- **Reproducible**: Bootstrap is seed-controlled; all outputs are
  deterministic given the seed.

### Limitations

- **Individual-level data required**: The package currently requires
  one-row-per-participant data with genotype, exposure, and outcomes.
  Summary-data-only MR is not supported.
- **Non-collapsibility**: The win ratio is not collapsible. The marginal
  DS-CWR differs from conditional win ratios. This is a property of the
  estimand, not a bug.
- **Linear PRS instrument**: The package assumes a linear polygenic risk
  score. Non-linear gene-gene interactions are not modeled.
- **Diagonal GWAS covariance**: The current bootstrap uses
  SNP-by-SNP standard errors, not the full LD covariance matrix.
  This may understate uncertainty when SNPs are in linkage
  disequilibrium.
- **Pleiotropy vulnerability**: The cCWR is sensitive to
  hierarchy-contaminating pleiotropy. Mortality-level pleiotropy as
  small as gamma = 0.05 can collapse coverage to 12%.
- **Sample size requirements**: Biobank-scale samples (N > 100,000) are
  a strict statistical prerequisite for reliable inference.

---

## Assumptions

The `mrwin` framework relies on the following assumptions:

1. **Relevance**: The PRS must be associated with the exposure. Weak
   instruments produce unstable estimates (check the `weak_instrument`
   warning).

2. **Independence**: The PRS must not be associated with confounders.
   This is satisfied by Mendelian segregation in homogeneous
   populations. Population stratification can violate this; use
   principal components as covariates.

3. **Exclusion restriction**: The PRS must affect the outcome only
   through the exposure. The SDPD tests this for the highest-priority
   endpoint. If rejected, the estimate may be invalid.

4. **Monotonicity**: The instrument must shift the exposure in the same
   direction for all individuals. Violations can bias the estimate.

5. **No effect modification by survival**: The causal effect should not
   vary between those who survive and those who do not. This is
   untestable.

6. **Correct priority ordering**: The clinical hierarchy must be
   specified before analysis. The package does not infer it from data.

---

## Warning Codes

| Code | Meaning | What to do |
|---|---|---|
| `weak_instrument` | Phenotypic shift near zero or Fieller CI unbounded | Treat estimate as unstable |
| `sdpd_rejected` | MR-Egger intercept suggests pleiotropy | Add sensitivity discussion |
| `sdpd_underpowered` | Too few SNPs for SDPD | Report power limitation |
| `positivity_failure` | IPTW ESS too small in a stratum | Inspect covariate overlap |
| `bridged_strata` | Strata skipped due to positivity failure | Report bridged schema |
| `discordant_components` | Component MR directions conflict | Do not overclaim global effect |

---

## Simulation Engine

Run scenario grids to validate the method:

```r
scenarios <- mrwin_scenarios(
  mrwin_config(n_outcome = 200, m_snps = 10, seed = 10),
  scenarios = c("A_null", "B_valid_IV", "C_pleiotropy", "D_hierarchy_discordant")
)

grid <- mrwin_run_simulation_grid(
  scenarios = scenarios, n_iter = 2,
  controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 11)
)

mrwin_simulation_summary(grid)
```

| Scenario | Meaning | Expected behavior |
|---|---|---|
| A_null | No causal effect | DS-CWR near 1, rejection rate near 5% |
| B_valid_IV | Valid causal effect | DS-CWR in expected direction |
| C_pleiotropy | Direct pleiotropy | SDPD rejection increases |
| D_hierarchy_discordant | Opposing component effects | Warning emitted |

---

## Citing This Package

If you use `mrwin` in your research, please cite:

> Nguyen Thien Minh, N. Ahmad Aziz. "Causal Win Statistics: Integrating
> Instrumental Variable Estimation with Hierarchical Composite Endpoints."
> *arXiv preprint*, 2026.

BibTeX:

```bibtex
@article{NguyenAziz2026,
  title={Causal Win Statistics: Integrating Instrumental Variable
         Estimation with Hierarchical Composite Endpoints},
  author={Nguyen Thien Minh and N. Ahmad Aziz},
  journal={arXiv preprint},
  year={2026}
}
```

**Authors:**
- Nguyen Thien Minh, University of Medicine and Pharmacy at Ho Chi Minh City,
  Vietnam (minhnt@ump.edu.vn)
- N. Ahmad Aziz (corresponding author), German Center for Neurodegenerative
  Diseases (DZNE) and University of Bonn, Germany (Ahmad.Aziz@dzne.de)

---

## Related References

**Win statistics:**
- Pocock SJ, et al. The win ratio: a new approach to the analysis of
  composite endpoints in clinical trials. *Eur Heart J*. 2012;33(14):1744-1749.
- Bebu I, Lachin JM. Large sample inference for a win ratio analysis of a
  composite endpoint based on prioritized components. *Biostatistics*.
  2016;17(1):178-191.
- Even Z, Josse A. Causal win ratio. *arXiv preprint*. 2025.

**Mendelian randomization:**
- Lawlor DA, et al. Mendelian randomization: using genes as instruments for
  making causal inferences in epidemiology. *Stat Med*. 2008;27(8):1133-1163.
- Davey Smith G, Hemani G. Mendelian randomization: genetic anchors for
  causal inference in epidemiological studies. *Hum Mol Genet*.
  2014;23(R1):R89-R98.
- Bowden J, et al. Mendelian randomization with invalid instruments: effect
  estimation and bias detection through Egger regression. *Int J Epidemiol*.
  2015;44(2):512-525.

**Polygenic risk scores:**
- Choi SW, et al. Tutorial: a guide to performing polygenic risk score
  analyses. *Nat Protoc*. 2020;15(9):2759-2772.

**Multi-state models:**
- Putter H, Fiocco M, Geskus RB. Tutorial in biostatistics: competing risks
  and multi-state models. *Stat Med*. 2007;26(11):2389-2430.

---

## Package Structure

```
mrwin/
  R/
    api.R              # High-level mrwin() workflow
    kernel.R           # Hierarchical comparison kernel
    estimate.R         # DS-CWR estimator, GLS pooling
    bootstrap.R        # Multiplier bootstrap inference
    adjustment.R       # Ordinal-IPTW, ESS, bridging
    sdpd.R             # Step-Down Pleiotropy Diagnostic
    simulate.R         # Simulation data-generating process
    simulation_engine.R # Scenario grid runner
    validation.R       # Input validation
    methods.R          # print, summary, plot, tidy, report
    benchmark.R        # Performance benchmarking
    backend_sparse.R   # Sparse backend
    config.R           # Simulation configuration
    strata.R           # PRS stratum assignment
  tests/
    testthat/          # 273 unit tests
  vignettes/           # 3 vignettes
  inst/spec/           # Implementation specifications
```

---

## License

MIT
