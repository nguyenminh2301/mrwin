# WP2 Data Validation Specification

Status: implemented as the public validation layer.

## Goal

WP2 adds a single validation gate before estimation. The high-level `mrwin()` workflow now calls `mrwin_validate_data()` before kernel construction, bootstrap, or SDPD diagnostics.

## Public Function

```r
validated <- mrwin_validate_data(
  data = cohort,
  endpoint = mrwin_endpoint(...),
  genotype = snp_columns,
  exposure = "X",
  gwas = mrwin_gwas(...),
  covariates = c("age", "sex", "PC1"),
  controls = mrwin_controls(...)
)
```

The return value has class:

```r
c("mrwin_validated_data", "list")
```

Returned fields:

- `time`
- `status`
- `G`
- `X`
- `covariates`
- `gwas`
- `endpoint`
- `controls`
- `strata`
- `rows_used`
- `issues`

## Error Checks

The validator stops before estimation on:

- endpoint `time`/`status` dimension mismatch.
- row-count mismatch across endpoint, genotype, exposure, or covariates.
- duplicate priority labels.
- genotype column count mismatch with GWAS beta length.
- GWAS SE length mismatch.
- missing/non-finite values unless `complete_cases = "drop"`.
- non-binary endpoint status.
- non-finite, negative, or disallowed zero event times.
- non-finite genotype, exposure, or covariate values.
- zero-variance genotype columns.
- non-PSD GWAS covariance.
- terminal priority observed before lower-priority observed events.
- more PRS strata than analysis rows.
- empty PRS strata.

Errors are structured internally as:

```r
list(code = "...", message = "...")
```

and are printed in the thrown error as `[code] message`.

## Warning Checks

The validator records warnings for:

- rows dropped under `complete_cases = "drop"`.
- all GWAS SEs equal to zero.
- small PRS strata with fewer than two observations.

Warnings flow into `fit$warnings` when validation is called through `mrwin()`.

## Current Boundaries

WP2 does not impute missing data and does not repair invalid endpoint ordering. It either drops incomplete rows when explicitly requested or stops. Production IPTW/GPS and ESS/bridging remain later work packages.

## Acceptance Criteria

WP2 is accepted when:

- `mrwin_validate_data()` is exported.
- `mrwin()` calls the validation layer.
- validation tests cover valid data, missing/drop behavior, endpoint errors, GWAS mismatch, zero-variance SNPs, and terminal-order violations.
