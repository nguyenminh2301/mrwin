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

# WP16 / M2: doubly-ranked stratification (Tian, Burgess et al., 2023).
# Form pre-strata of size `n_strata` along the instrument (PRS) ranking, then
# within each pre-stratum rank by the exposure and assign the exposure rank as
# the final stratum. This yields strata that differ in mean exposure while
# preserving the instrument ordering, valid under a weaker rank-preservation
# assumption than the linearity/homogeneity needed by raw PRS-rank ("residual")
# stratification. Deterministic ties: instrument order breaks ties by original
# index; within-block exposure order breaks ties by PRS then index. A trailing
# incomplete pre-stratum (when N is not a multiple of n_strata) is spread across
# 1..n_strata by scaled exposure rank.
mrwin_doubly_ranked_strata <- function(prs, exposure, n_strata = 10L) {
  prs <- as.numeric(prs)
  exposure <- as.numeric(exposure)
  if (length(prs) != length(exposure)) {
    stop("`prs` and `exposure` must have the same length.", call. = FALSE)
  }
  if (n_strata < 2L) {
    stop("`n_strata` must be at least 2.", call. = FALSE)
  }
  n <- length(prs)
  if (n < n_strata) {
    stop("`n_strata` cannot exceed the number of individuals.", call. = FALSE)
  }
  if (any(!is.finite(prs)) || any(!is.finite(exposure))) {
    stop("`prs` and `exposure` must be finite.", call. = FALSE)
  }
  n_strata <- as.integer(n_strata)

  # instrument ranking (deterministic: ties broken by original index)
  ord_prs <- order(prs)
  r_prs <- integer(n)
  r_prs[ord_prs] <- seq_len(n)
  pre_id <- (r_prs - 1L) %/% n_strata

  strata <- integer(n)
  for (blk in sort(unique(pre_id))) {
    idx <- which(pre_id == blk)
    bs <- length(idx)
    # within-block exposure ranking (ties broken by PRS then index)
    o <- order(exposure[idx], prs[idx], idx)
    rank_in_block <- integer(bs)
    rank_in_block[o] <- seq_len(bs)
    if (bs == n_strata) {
      strata[idx] <- rank_in_block
    } else {
      strata[idx] <- pmin(pmax(ceiling(rank_in_block * n_strata / bs), 1L), n_strata)
    }
  }

  list(score = prs, strata = as.integer(strata))
}

# Internal stratification dispatch used across the estimator/bootstrap paths.
# "prs_rank" = rank-bin the PRS (the default). "doubly_ranked" = Tian/Burgess
# doubly-ranked strata (needs the exposure X). Returns list(score, strata),
# matching the mrwin_prs_strata contract.
.mrwin_assign_strata <- function(G, beta, X, n_strata,
                                 stratification = c("prs_rank", "doubly_ranked")) {
  stratification <- match.arg(stratification)
  score <- drop(as.matrix(G) %*% as.numeric(beta))
  if (stratification == "doubly_ranked") {
    if (is.null(X)) stop("`X` is required for doubly-ranked stratification.", call. = FALSE)
    return(mrwin_doubly_ranked_strata(score, X, n_strata = n_strata))
  }
  strata <- ceiling(rank(score, ties.method = "first") * n_strata / length(score))
  strata <- pmin(pmax(strata, 1L), n_strata)
  list(score = score, strata = as.integer(strata))
}
