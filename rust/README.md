# mrwin Rust backend

Compiled implementation of the exact sub-quadratic right-censored hierarchical
win/loss algorithm (`inst/spec/fast-hierarchical-win-algorithm.md`). It is an
optional performance backend; the package builds, checks, and runs fully
without it using the `dense` and `sparse` R backends.

## Layout

```
rust/
  mrwinkernel/        # pure-Rust core: the algorithm. No external deps.
    src/lib.rs        #   dominance engine (Fenwick + value-pivot CDQ) + tests
    tests/parity.rs   #   cross-language parity vs the Python reference
    tests/parity_fixtures.txt
  mrwinrust/          # extendr bindings exposing the core to R (needs libR)
    src/lib.rs
  packaging/          # templates to wire extendr into the R package's src/
    Makevars
    Makevars.win
    entrypoint.c
```

## What is verified here

The **core crate is fully tested with no R and no network**:

```bash
cd rust/mrwinkernel
cargo test            # unit tests (fast vs brute, all K) + Python parity
```

The tests prove the Rust core reproduces the dense kernel exactly: bit-exact for
unweighted counts and within 1e-9 for weighted totals, for K = 1..4, including
heavy ties and right-censoring. The Python parity fixture is regenerated with:

```bash
PYTHONPATH=python python3 -m p1_engine_v5.gen_parity_fixtures \
    rust/mrwinkernel/tests/parity_fixtures.txt
```

## What still needs an R toolchain

The `mrwinrust` extendr crate and the R<->Rust linkage **require R (libR) plus a
Rust toolchain at build time** and are therefore not compiled in R-less CI. To
enable `mrwin(backend = "rust")` in an R environment:

1. Install Rust (https://rustup.rs) and the R package `rextendr`.
2. From the package root, scaffold the standard extendr `src/` layout
   (`rextendr::use_extendr()` or copy `rust/packaging/*` into `src/`, with the
   `mrwinrust` crate placed at `src/rust`). The crate path-depends on
   `mrwinkernel`, so vendor or co-locate both crates under `src/rust/`.
3. `rextendr::register_extendr()` to (re)generate `R/extendr-wrappers.R`, then
   `R CMD INSTALL .`.
4. Verify: `mrwin:::.mrwin_rust_available()` returns `TRUE`, and
   `tests/testthat/test-backend-rust.R` (skipped otherwise) checks parity
   against the `sparse` backend at 1e-10.

The R wrapper `R/backend_rust.R` already guards on availability, so installing
without Rust simply leaves the backend disabled with an actionable error.

## Status

| Component | State |
|---|---|
| `mrwinkernel` core algorithm | Implemented, tested (cargo test) |
| Cross-language parity vs Python | Verified (bit-exact / 1e-9) |
| `mrwinrust` extendr bindings | Source provided; build needs libR |
| R `src/` packaging | Templates in `packaging/`; maintainer step |
| `backend = "rust"` R wiring | Implemented + guarded in `R/backend_rust.R` |
