#!/usr/bin/env Rscript
# ============================================================================
# PROBE #2, attempt 3 -- lift the AR-style tau test to the setting where the
# portability problem ACTUALLY LIVES: TWO-SAMPLE, MULTI-SNP, within-family.
#
# WHY. p2-ar-style-tau-test.R built a rigorous AR-type test for tau but only in
# the ONE-SAMPLE, SCALAR-instrument case, where it turned out statistically
# equivalent to the naive ar$ci/D_hat plug-in and bounded only 2-7% of the time
# under a weak instrument. That is the wrong arena on two counts:
#   (a) the portability question ("same causal architecture, different cohort
#       follow-up/censoring") is inherently a TWO-SAMPLE / cross-biobank problem;
#   (b) a single SNP cannot pin tau down, so the test had almost no power to
#       return a bounded set. Paper 04's deployable pipeline has L=25 SNPs.
# With L instruments the statistic is chi^2_L rather than chi^2_1, so the set
# should be bounded far more often -- the regime where an AR-type construction
# can actually beat a plug-in rather than merely match it.
#
# DERIVATION. Per-SNP two-sample moments: delta_l = gamma*beta_l with
# gamma = tau*D. Under H0: tau = t, the residual is  r_l(t) = delta_l - t*D*beta_l.
# Propagating BOTH the per-SNP SEs and D_hat's own sampling error (delta method;
# D_hat comes from the outcome sample, beta_l from the exposure regression, so
# they are treated as independent):
#     Var(t*D_hat*beta_l) ~= t^2 * [ D_hat^2 * se_beta_l^2 + beta_l^2 * se_D^2 ]
#     AR_tau(t) = SUM_l  r_l(t)^2 / ( se_delta_l^2 + t^2*(D_hat^2*se_beta_l^2
#                                                        + beta_l^2*se_D^2) )
# which is chi^2_L at the true tau regardless of instrument strength (it never
# divides by beta or by D), and min_t AR_tau(t) ~ chi^2_{L-1} is an over-ID /
# pleiotropy test on the tau scale. Confidence set = {t : AR_tau(t) <= q}.
#
# se_D is CLUSTER-ROBUST by family: D_hat = mean(d_i) is a degree-2 U-statistic
# with Hajek influence IF_D_i = 2*(d_i - D_hat), and sibs within a family are
# correlated, so  Var(D_hat) = (1/n^2) * SUM_families ( SUM_{i in fam} IF_D_i )^2.
#
# TEST: coverage of the portable target tau* across a censoring ladder (cohorts
# differing ONLY in follow-up), vs. the naive plug-in, plus the bounded-set rate
# and the over-ID type-I. Public sim + package only.
# ============================================================================
suppressMessages(library(mrwin)); ncores <- max(1L, min(4L, parallel::detectCores())); cfg <- mrwin_config()
subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")

outc <- function(X, u, w, cn, uf) {
  F <- length(X); et <- matrix(NA_real_, F, 3); lw <- log(w)
  for (j in 1:3) { lp <- cfg$alpha_x[j] * X + cfg$nu_u[j] * u + lw
    et[, j] <- (-log(uf[, j]) / (cfg$baseline_haz[j] * exp(lp)))^(1 / cfg$shape_weibull) }
  td <- et[, 1]; th <- et[, 2]; tr <- et[, 3]
  list(time = cbind(pmin(td, cn), pmin(th, cn, td), pmin(tr, cn, td)),
       status = cbind(as.integer(td <= cn), as.integer(th <= cn & th <= td), as.integer(tr <= cn & tr <= td)))
}
tr1 <- function(G) ifelse(G == 2, 1L, ifelse(G == 0, 0L, rbinom(length(G), 1, 0.5)))
mkarch <- function(seed = 11L, M = 25L) { set.seed(seed); list(base = runif(M, .15, .35), betas = rnorm(M, 0, .15), M = M) }
genSib <- function(Fn, seed, arch, cens) {
  set.seed(seed); M <- arch$M; base <- arch$base; betas <- arch$betas
  pp <- pmin(pmax(matrix(base, Fn, M, byrow = TRUE), .02), .95)
  Gf <- matrix(rbinom(Fn * M, 2, pp), Fn, M); Gm <- matrix(rbinom(Fn * M, 2, pp), Fn, M)
  G <- NULL; fam <- NULL; X <- NULL; U <- NULL; Wf <- NULL; CN <- NULL; UF <- NULL
  for (k in 1:2) {
    Gs <- matrix(tr1(Gf), Fn, M) + matrix(tr1(Gm), Fn, M); Z <- as.numeric(Gs %*% betas)
    u <- rnorm(Fn); x <- cfg$alpha_s * Z + cfg$alpha_u * u + rnorm(Fn)
    wf <- rgamma(Fn, 1 / cfg$theta_f, scale = cfg$theta_f)
    cn <- pmin(rexp(Fn, cens), cfg$max_follow_up); uf <- matrix(runif(Fn * 3), Fn, 3)
    G <- rbind(G, Gs); fam <- c(fam, 1:Fn); X <- c(X, x); U <- c(U, u); Wf <- c(Wf, wf); CN <- c(CN, cn); UF <- rbind(UF, uf)
  }
  oc <- outc(X, U, Wf, CN, UF)
  list(fam = fam, G = G, X = X, u = U, w = Wf, cn = CN, uf = UF, time = oc$time, status = oc$status, M = M)
}
oracle <- function(d, eps = 0.25) {
  np <- length(d$X)
  Bp <- outc(d$X + eps, d$u, d$w, d$cn, d$uf); Bm <- outc(d$X - eps, d$u, d$w, d$cn, d$uf)
  v <- mrwin_fast_pair_win_loss(Bp$time, Bp$status, Bm$time, Bm$status)
  (as.numeric(v[["wins"]]) - as.numeric(v[["losses"]])) / (as.numeric(np)^2) / (2 * eps)
}
winscore <- function(d) { c <- subj(d$time, d$status, d$time, d$status, rep(1, length(d$X))); (c[, 1] - c[, 2]) / (length(d$X) - 1) }
# D_hat plus a CLUSTER-ROBUST (by family) SE via the Hajek influence function
D_hat_se <- function(d) {
  n <- length(d$X); cw <- subj(d$time, d$status, d$time, d$status, rep(1, n))
  d_i <- (cw[, 1] + cw[, 2]) / (n - 1); Dh <- mean(d_i)
  IF <- 2 * (d_i - Dh)
  S <- tapply(IF, d$fam, sum)
  list(D = Dh, se = sqrt(sum(S^2)) / n)
}
gwas_fam <- function(d, wv) {
  N <- length(d$X); G <- d$G; M <- d$M; cnt <- as.vector(table(d$fam))
  cX <- as.numeric(d$X - tapply(d$X, d$fam, mean)[d$fam]); cW <- as.numeric(wv - tapply(wv, d$fam, mean)[d$fam])
  fmG <- rowsum(G, d$fam) / cnt; cG <- G - fmG[d$fam, ]; df <- N - length(cnt) - 1L
  den <- colSums(cG^2); numW <- colSums(cG * cW); numX <- colSums(cG * cX)
  rssW <- sum(cW^2) - numW^2 / den; rssX <- sum(cX^2) - numX^2 / den
  list(beta = numX / den, se_beta = sqrt((rssX / df) / den), delta = numW / den, se_delta = sqrt((rssW / df) / den))
}

# multi-SNP two-sample AR test on the TAU scale
ar_tau_twosample <- function(gf, Dh, seD, alpha = 0.05, n_grid = 6000L, span = 40) {
  L <- length(gf$beta)
  ivw <- sum(gf$beta * gf$delta / gf$se_delta^2) / sum(gf$beta^2 / gf$se_delta^2)
  se_ivw <- sqrt(1 / sum(gf$beta^2 / gf$se_delta^2))
  tau_c <- ivw / Dh
  hw <- span * max(se_ivw / max(Dh, 1e-8), abs(tau_c) * 0.5, 1e-3)
  grid <- seq(tau_c - hw, tau_c + hw, length.out = n_grid)
  q <- stats::qchisq(1 - alpha, L)
  ARv <- vapply(grid, function(t) {
    num <- (gf$delta - t * Dh * gf$beta)^2
    den <- gf$se_delta^2 + t^2 * (Dh^2 * gf$se_beta^2 + gf$beta^2 * seD^2)
    sum(num / den)
  }, numeric(1))
  acc <- grid[ARv <= q]
  Q <- min(ARv); Qp <- stats::pchisq(Q, L - 1L, lower.tail = FALSE)
  if (!length(acc)) return(list(ci = c(NA, NA), bounded = FALSE, empty = TRUE, tau_point = grid[which.min(ARv)], Q = Q, Q_p = Qp))
  lo <- min(acc); hi <- max(acc)
  bounded <- (lo > grid[1] + 1e-12) && (hi < grid[n_grid] - 1e-12)
  list(ci = c(lo, hi), bounded = bounded, empty = FALSE, tau_point = grid[which.min(ARv)], Q = Q, Q_p = Qp)
}

arch <- mkarch(); Fn <- 1500L; R <- 60L
cens_grid <- c(0.02, 0.1, 0.3, 0.8, 2.0, 5.0)

cat("== PROBE 2 attempt 3: TWO-SAMPLE MULTI-SNP (L=25) AR test on the tau scale ==\n\n")
dref <- genSib(30000L, 999L, arch, cens_grid[1])
Dref <- D_hat_se(dref); tau_star <- oracle(dref) / Dref$D
cat(sprintf("Reference: oracle_ref=%.4f  D_ref=%.4f  tau*=%.4f\n\n", oracle(dref), Dref$D, tau_star))

cat(sprintf("F=%d fam x2 sibs, M=25 SNPs, R=%d; cohorts differ ONLY in censoring\n\n", Fn, R))
cat("APPLES-TO-APPLES NOTE: an UNBOUNDED confidence set trivially covers, so raw\n")
cat("coverage is only interpretable alongside the bounded rate. Naive coverage is\n")
cat("scored GENEROUSLY (range as a conservative superset, contiguity NOT required),\n")
cat("so any remaining AR-tau advantage is not a counting artifact.\n\n")
cat(sprintf("%-7s %-8s %-9s %-9s %-9s %-9s %-9s %-9s\n",
            "cens", "D_hat", "ARt bnd%", "nv bnd%", "ARt cov", "nv cov", "ARt c|bnd", "nv c|bnd"))
rows <- list()
for (cens in cens_grid) {
  rr <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
    d <- genSib(Fn, 7000L + round(cens * 1000) + r, arch, cens)
    wv <- winscore(d); ds <- D_hat_se(d); gf <- gwas_fam(d, wv)
    at <- ar_tau_twosample(gf, ds$D, ds$se)
    nv <- mrwin_winmr_ar(gf$beta, gf$delta, gf$se_delta, se_gx = gf$se_beta)
    nci <- sort(nv$ci / ds$D)   # D_hat>0 so order preserved; sort() guards sign conventions
    nv_bnd <- as.numeric(isTRUE(nv$bounded))
    c(D = ds$D, seD = ds$se,
      bnd = as.numeric(at$bounded), nbnd = nv_bnd,
      cov = as.numeric(!at$empty && tau_star >= at$ci[1] && tau_star <= at$ci[2]),
      ncov = as.numeric(tau_star >= nci[1] && tau_star <= nci[2]),  # generous: no contiguity req.
      qp = as.numeric(at$Q_p < 0.05),
      wid = if (at$bounded) at$ci[2] - at$ci[1] else NA_real_,
      nwid = if (nv_bnd == 1) nci[2] - nci[1] else NA_real_)
  }, mc.cores = ncores))
  rr <- rr[stats::complete.cases(rr[, 1:7]), , drop = FALSE]
  m <- colMeans(rr); rows[[as.character(cens)]] <- rr
  cb <- if (sum(rr[, "bnd"] == 1) > 0) mean(rr[rr[, "bnd"] == 1, "cov"]) else NA
  nb <- if (sum(rr[, "nbnd"] == 1) > 0) mean(rr[rr[, "nbnd"] == 1, "ncov"]) else NA
  cat(sprintf("%-7.2f %-8.4f %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f\n",
              cens, m["D"], m["bnd"], m["nbnd"], m["cov"], m["ncov"], cb, nb))
}

cat("\n== Head-to-head on the SAME replicates where BOTH sets are bounded ==\n")
cat(sprintf("%-7s %-6s %-11s %-11s %-11s %-11s %-9s\n",
            "cens", "n_both", "ARt cov", "naive cov", "ARt wid", "naive wid", "wid ratio"))
for (cens in cens_grid) {
  rr <- rows[[as.character(cens)]]
  ok <- rr[, "bnd"] == 1 & rr[, "nbnd"] == 1 & !is.na(rr[, "wid"]) & !is.na(rr[, "nwid"])
  if (sum(ok) < 3) { cat(sprintf("%-7.2f %-6d %-11s\n", cens, sum(ok), "(too few)")); next }
  s <- rr[ok, , drop = FALSE]
  cat(sprintf("%-7.2f %-6d %-11.3f %-11.3f %-11.4f %-11.4f %-9.3f\n",
              cens, sum(ok), mean(s[, "cov"]), mean(s[, "ncov"]),
              median(s[, "wid"]), median(s[, "nwid"]),
              median(s[, "wid"]) / median(s[, "nwid"])))
}
cat("\n== over-ID (pleiotropy) test on the tau scale: type-I under a VALID instrument ==\n")
for (cens in cens_grid) {
  cat(sprintf("  cens=%.2f  P(Q_p<0.05) = %.3f\n", cens, mean(rows[[as.character(cens)]][, "qp"])))
}
cat("DONE-P2-TWOSAMPLE-ARTAU\n")
