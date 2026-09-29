#!/usr/bin/env bash
#
# run-baseline.sh — reproducible local verification baseline for mrwin.
#
# Stage G0.1 of inst/spec/completion-roadmap-academic.vi.md. Every gate runs on
# the local machine; no cloud service is used. Record the summary it prints in
# inst/spec/baseline-verification.md after each run that should become the new
# reference.
#
# Usage (from anywhere inside the repository):
#   bash tools/baseline/run-baseline.sh [--skip-python] [--skip-coverage] [--offline]
#
#   --skip-python    do not run the Python oracle tests
#   --skip-coverage  do not run covr (the slowest gate)
#   --offline        skip the network-dependent CRAN-incoming checks of --as-cran
#
# Requirements: R >= 4.1 with Rcpp, testthat, knitr, rmarkdown, covr, xml2;
# pandoc; qpdf. A LaTeX installation is optional: without pdflatex the PDF
# manual is skipped (--no-manual) and the summary says so. Python >= 3.11 with
# numpy, scipy, statsmodels, pytest for the oracle tests.
#
# Output: out/baseline-<YYYYMMDD>/ (gitignored) with one log per gate and
# summary.txt. Exit status is non-zero when any gate fails.

set -uo pipefail

SKIP_PYTHON=0
SKIP_COVERAGE=0
OFFLINE=0
for arg in "$@"; do
  case "$arg" in
    --skip-python) SKIP_PYTHON=1 ;;
    --skip-coverage) SKIP_COVERAGE=1 ;;
    --offline) OFFLINE=1 ;;
    *) echo "Unknown option: $arg" >&2; exit 2 ;;
  esac
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$ROOT/out/baseline-$(date +%Y%m%d)"
mkdir -p "$OUT"
SUMMARY="$OUT/summary.txt"
: > "$SUMMARY"
FAILED=0

note() { echo "$*" | tee -a "$SUMMARY"; }

note "mrwin local baseline"
note "  date:   $(date -u +%Y-%m-%dT%H:%M:%SZ)"
DIRTY=""
if [[ -n "$(git -C "$ROOT" status --porcelain --untracked-files=no 2>/dev/null)" ]]; then
  DIRTY=" + uncommitted changes"
fi
note "  commit: $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)${DIRTY}"
note "  output: $OUT"
note ""

# ---------------------------------------------------------------------------
# 0. Environment
# ---------------------------------------------------------------------------
{
  Rscript -e 'print(sessionInfo()); cat("\n"); for (p in c("Rcpp","testthat","knitr","rmarkdown","covr","xml2")) cat(sprintf("%-10s %s\n", p, tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")))'
  echo
  echo "pandoc:   $(pandoc --version 2>/dev/null | head -1 || echo MISSING)"
  echo "qpdf:     $(qpdf --version 2>/dev/null | head -1 || echo MISSING)"
  echo "pdflatex: $(pdflatex --version 2>/dev/null | head -1 || echo MISSING)"
  echo "python:   $(python3 --version 2>/dev/null || echo MISSING)"
} > "$OUT/session-info.txt" 2>&1
note "[env] $(grep -m1 '^R version' "$OUT/session-info.txt")"

# ---------------------------------------------------------------------------
# 1. testthat (source tree)
# ---------------------------------------------------------------------------
(
  cd "$ROOT" && Rscript -e '
    res <- testthat::test_local(".", reporter = testthat::SummaryReporter$new(), stop_on_failure = FALSE)
    df <- as.data.frame(res)
    cat(sprintf("\nTESTTHAT files=%d blocks=%d expectations=%d failed=%d errors=%d skipped=%d warnings=%d\n",
      length(unique(df$file)), nrow(df), sum(df$nb), sum(df$failed), sum(df$error),
      sum(df$skipped), sum(df$warning)))
  '
) > "$OUT/testthat.log" 2>&1
TT_LINE="$(grep '^TESTTHAT' "$OUT/testthat.log" | tail -1)"
note "[testthat] ${TT_LINE:-no summary line (see testthat.log)}"
if [[ -z "$TT_LINE" ]] || ! grep -q 'failed=0 errors=0' <<<"$TT_LINE"; then FAILED=1; fi

# ---------------------------------------------------------------------------
# 2. R CMD build + R CMD check --as-cran
# ---------------------------------------------------------------------------
( cd "$OUT" && R CMD build "$ROOT" ) > "$OUT/build.log" 2>&1
TARBALL="$(ls -t "$OUT"/mrwin_*.tar.gz 2>/dev/null | head -1)"
if [[ -z "$TARBALL" ]]; then
  note "[build] FAILED (see build.log)"
  FAILED=1
else
  note "[build] $(basename "$TARBALL")"
  CHECK_ARGS=(--as-cran)
  if ! command -v pdflatex >/dev/null 2>&1; then
    CHECK_ARGS+=(--no-manual)
    note "[check] pdflatex not found: PDF manual skipped (--no-manual)"
  fi
  if [[ "$OFFLINE" == 1 ]]; then
    export _R_CHECK_CRAN_INCOMING_REMOTE_=false
    note "[check] offline: remote CRAN-incoming checks skipped"
  fi
  ( cd "$OUT" && R CMD check "${CHECK_ARGS[@]}" "$TARBALL" ) > "$OUT/check.log" 2>&1
  STATUS_LINE="$(grep -m1 '^Status:' "$OUT/check.log")"
  note "[check] ${STATUS_LINE:-no Status line (see check.log)}"
  grep -E '^\* checking .*\.\.\. (WARNING|NOTE|ERROR)' "$OUT/check.log" | sed 's/^/          /' | tee -a "$SUMMARY" >/dev/null
  if [[ -z "$STATUS_LINE" ]] || grep -q 'ERROR' <<<"$STATUS_LINE"; then FAILED=1; fi
fi

# ---------------------------------------------------------------------------
# 3. Python oracle tests
# ---------------------------------------------------------------------------
if [[ "$SKIP_PYTHON" == 1 ]]; then
  note "[pytest] skipped (--skip-python)"
else
  ( cd "$ROOT" && PYTHONPATH="$ROOT/python${PYTHONPATH:+:$PYTHONPATH}" python3 -m pytest tests/python -q ) > "$OUT/pytest.log" 2>&1
  PY_EXIT=$?
  note "[pytest] $(tail -1 "$OUT/pytest.log")"
  if [[ "$PY_EXIT" != 0 ]]; then FAILED=1; fi
fi

# ---------------------------------------------------------------------------
# 4. Coverage
# ---------------------------------------------------------------------------
if [[ "$SKIP_COVERAGE" == 1 ]]; then
  note "[covr] skipped (--skip-coverage)"
else
  (
    cd "$ROOT" && Rscript -e '
      cov <- covr::package_coverage(".")
      print(cov, group = "filename")
      cat(sprintf("\nCOVERAGE total=%.2f%%\n", covr::percent_coverage(cov)))
    '
  ) > "$OUT/coverage.txt" 2>&1
  COV_LINE="$(grep '^COVERAGE' "$OUT/coverage.txt" | tail -1)"
  note "[covr] ${COV_LINE:-no summary line (see coverage.txt)}"
  if [[ -z "$COV_LINE" ]]; then FAILED=1; fi
fi

note ""
if [[ "$FAILED" == 0 ]]; then note "BASELINE: all gates passed"; else note "BASELINE: at least one gate failed"; fi
exit "$FAILED"
