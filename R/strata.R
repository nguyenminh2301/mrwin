mrwin_prs_strata <- function(G, beta, n_strata = 10L) {
  G <- as.matrix(G)
  beta <- as.numeric(beta)
  if (ncol(G) != length(beta)) {
    stop("`G` columns must match length of `beta`.", call. = FALSE)
  }
  if (n_strata < 2L) {
    stop("`n_strata` must be at least 2.", call. = FALSE)
  }

  score <- drop(G %*% beta)
  strata <- ceiling(rank(score, ties.method = "first") * n_strata / length(score))
  strata <- pmin(pmax(strata, 1L), n_strata)
  list(score = score, strata = as.integer(strata))
}
