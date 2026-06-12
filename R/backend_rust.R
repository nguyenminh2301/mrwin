#' Rust performance backend (extendr)
#'
#' These functions route the hot win/loss computation through the compiled
#' `mrwinrust` extendr library (a binding over the dependency-free `mrwinkernel`
#' crate, which implements the exact sub-quadratic algorithm specified in
#' `inst/spec/fast-hierarchical-win-algorithm.md`). The Rust backend reuses the
#' entire sparse-backend pipeline (stratification, IPTW, GLS, Fieller); it only
#' replaces the per-stratum-pair win/loss kernel.
#'
#' The compiled library is optional. It is built only when the package is
#' installed with a Rust toolchain present (see `rust/README.md`). When it
#' is unavailable, `mrwin(backend = "rust")` fails with an actionable message
#' and the `dense`/`sparse` backends remain fully functional.

#' Report whether the compiled Rust kernel is available in this installation.
#' @return `TRUE` if the extendr symbol is loaded, otherwise `FALSE`.
#' @keywords internal
.mrwin_rust_available <- function() {
  tryCatch(isTRUE(is.loaded("wrap__rust_pair_win_loss")), error = function(e) FALSE)
}

#' Win/loss for one stratum pair via the compiled Rust kernel.
#'
#' Drop-in replacement for [mrwin_pair_win_loss()] with the same signature and
#' return contract (`c(wins, losses, total)`). Matrices are flattened row-major
#' for the Rust ABI; `k` is the number of priority columns.
#' @keywords internal
.mrwin_rust_pair_win_loss <- function(
    time_high, status_high, time_low, status_low,
    weights_high = NULL, weights_low = NULL) {
  th <- as.matrix(time_high)
  sh <- as.matrix(status_high)
  tl <- as.matrix(time_low)
  sl <- as.matrix(status_low)
  k <- ncol(th)
  n_high <- nrow(th)
  n_low <- nrow(tl)
  if (is.null(weights_high)) weights_high <- rep(1, n_high)
  if (is.null(weights_low)) weights_low <- rep(1, n_low)

  # row-major flattening: as.double(t(M))
  res <- .Call(
    "wrap__rust_pair_win_loss",
    as.double(t(th)), as.double(t(sh)), as.double(weights_high),
    as.double(t(tl)), as.double(t(sl)), as.double(weights_low),
    as.integer(k),
    PACKAGE = "mrwin"
  )
  c(
    wins = as.numeric(res[1]),
    losses = as.numeric(res[2]),
    total = sum(weights_high) * sum(weights_low)
  )
}

#' Multiplier bootstrap using the Rust win/loss kernel
#'
#' Identical contract to [mrwin_sparse_bootstrap()]; it delegates to that
#' function with the compiled Rust kernel injected as `pair_fun`, so all
#' downstream estimation is shared and tested. Requires the compiled extendr
#' library (see `rust/README.md`).
#'
#' @param ... Arguments forwarded to [mrwin_sparse_bootstrap()].
#' @return An object of class `mrwin_bootstrap`, as from
#'   [mrwin_sparse_bootstrap()].
#' @export
mrwin_rust_bootstrap <- function(...) {
  if (!.mrwin_rust_available()) {
    stop(
      "backend = 'rust' requires the compiled extendr library, which was not ",
      "found in this installation. Install the package with a Rust toolchain ",
      "available (see rust/README.md), or use backend = 'dense' / 'sparse'.",
      call. = FALSE
    )
  }
  mrwin_sparse_bootstrap(..., pair_fun = .mrwin_rust_pair_win_loss)
}
