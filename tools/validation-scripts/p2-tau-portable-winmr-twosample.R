#!/usr/bin/env Rscript
# ============================================================================
# CROSS-CUTTING PROBE #2: tau = D-corrected win-MR gradient AS AN MR ESTIMAND,
# chained into the ACTUAL two-sample within-family estimation pipeline (not
# just the population oracle sweep of p5-collapsible-transport-hierarchical.R).
#
# MOTIVATION. dev/findings-transport-collapsibility.md shows the raw win effect
# gradient is SQUASHED by censoring through the decidability functional D (a
# non-collapsibility phenomenon: two cohorts with the SAME causal architecture
# but different follow-up/censoring have DIFFERENT raw win-MR targets, because
# the oracle win gradient itself depends on D). Two-sample win-ratio MR is
# routinely applied ACROSS cohorts with different follow-up length (different
# biobanks, different censoring rates) -- so the raw win-MR estimate is not
# directly comparable/poolable across such cohorts. This script asks: does
# dividing the AR/IVW point estimate by a CHEAP, data-estimable decidability
# D_hat (fraction of decided pairs in the pooled outcome margin -- already
# validated as an (almost) free byproduct of the win-score kernel) restore
# portability -- i.e. does tau_hat = gamma_hat / D_hat stay close to a common
# reference across a censoring ladder, while raw gamma_hat does not?
#
# Design: SAME causal architecture (mode="none", fixed SNP effects), sweep
# ONLY the censoring rate across "cohorts" (as if from biobanks with
# different follow-up). Ground truth = each cohort's own interventional
# oracle (do(X+-eps)); the reference tau* is oracle_ref/D_ref from a
# near-uncensored, large-N reference cohort. Public sim + package only.
# ============================================================================
suppressMessages(library(mrwin)); ncores <- max(1L, min(4L, parallel::detectCores())); cfg <- mrwin_config()
subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")

outc <- function(X, u, w, cn, uf, cens_shape = NULL) {
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
    wf <- rgamma(Fn, 1 / cfg$theta_f, scale = cfg$theta_f); cn <- pmin(rexp(Fn, cens), cfg$max_follow_up); uf <- matrix(runif(Fn * 3), Fn, 3)
    G <- rbind(G, Gs); fam <- c(fam, 1:Fn); X <- c(X, x); U <- c(U, u); Wf <- c(Wf, wf); CN <- c(CN, cn); UF <- rbind(UF, uf)
  }
  oc <- outc(X, U, Wf, CN, UF); list(fam = fam, G = G, X = X, u = U, w = Wf, cn = CN, uf = UF, time = oc$time, status = oc$status, M = M)
}
oracle <- function(d, eps = 0.25) {
  np <- length(d$X)
  Bp <- outc(d$X + eps, d$u, d$w, d$cn, d$uf); Bm <- outc(d$X - eps, d$u, d$w, d$cn, d$uf)
  v <- mrwin_fast_pair_win_loss(Bp$time, Bp$status, Bm$time, Bm$status)
  (as.numeric(v[["wins"]]) - as.numeric(v[["losses"]])) / (as.numeric(np)^2) / (2 * eps)
}
winscore <- function(d) { c <- subj(d$time, d$status, d$time, d$status, rep(1, length(d$X))); (c[, 1] - c[, 2]) / (length(d$X) - 1) }
D_hat_fn <- function(d) { n <- length(d$X); cw <- subj(d$time, d$status, d$time, d$status, rep(1, n)); mean((cw[, 1] + cw[, 2]) / (n - 1)) }
gwas_fam <- function(d, wv) {
  N <- length(d$X); G <- d$G; M <- d$M; cnt <- as.vector(table(d$fam))
  cX <- as.numeric(d$X - tapply(d$X, d$fam, mean)[d$fam]); cW <- as.numeric(wv - tapply(wv, d$fam, mean)[d$fam])
  fmG <- rowsum(G, d$fam) / cnt; cG <- G - fmG[d$fam, ]; df <- N - length(cnt) - 1L
  den <- colSums(cG^2); numW <- colSums(cG * cW); numX <- colSums(cG * cX)
  rssW <- sum(cW^2) - numW^2 / den; rssX <- sum(cX^2) - numX^2 / den
  list(beta = numX / den, se_beta = sqrt((rssX / df) / den), delta = numW / den, se_delta = sqrt((rssW / df) / den))
}

arch <- mkarch(); Fn <- 1500L; R <- 80L
cens_grid <- c(0.02, 0.1, 0.3, 0.8, 2.0, 5.0)

cat("== Reference (near-uncensored, large N): oracle_ref, D_ref, tau* = oracle_ref/D_ref ==\n")
dref <- genSib(30000L, 999L, arch, cens_grid[1])
oracle_ref <- oracle(dref); D_ref <- D_hat_fn(dref); tau_star <- oracle_ref / D_ref
cat(sprintf("oracle_ref=%.4f  D_ref=%.4f  tau*=%.4f\n\n", oracle_ref, D_ref, tau_star))

cat("== Censoring ladder (SAME architecture, F=1500 fam x2, M=25 SNPs, R=80) ==\n")
cat(sprintf("%-8s %-10s %-10s %-10s %-20s %-20s %-10s %-10s %-10s\n",
            "cens", "oracle(c)", "D_hat(c)", "orac/D", "gamma_AR (bias)", "tau_hat=g/D (bias v tau*)",
            "AR-cov(own)", "tau-cov(tau*)", "rawCI~tau*"))
rows <- list(); raw_all <- list(); tau_all <- list()
for (cens in cens_grid) {
  rr <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
    d <- genSib(Fn, 6000L + round(cens * 1000) + r, arch, cens)
    o <- oracle(d); wv <- winscore(d); Dh <- D_hat_fn(d)
    gf <- gwas_fam(d, wv)
    ar <- mrwin_winmr_ar(gf$beta, gf$delta, gf$se_delta, se_gx = gf$se_beta)
    tau_hat <- ar$gamma / Dh; tau_ci <- ar$ci / Dh
    c(o = o, D = Dh, g = ar$gamma, tau = tau_hat,
      cov_own = as.numeric(ar$contiguous && o >= ar$ci[1] && o <= ar$ci[2]),
      cov_tau = as.numeric(ar$contiguous && tau_star >= tau_ci[1] && tau_star <= tau_ci[2]),
      rawci_covers_taustar = as.numeric(ar$contiguous && tau_star >= ar$ci[1] && tau_star <= ar$ci[2]),
      ci_lo = ar$ci[1], ci_hi = ar$ci[2])
  }, mc.cores = ncores))
  m <- colMeans(rr, na.rm = TRUE); med <- apply(rr, 2, median, na.rm = TRUE)
  rows[[as.character(cens)]] <- m
  cat(sprintf("%-8.2f %-10.4f %-10.4f %-10.4f %-20s %-20s %-10.3f %-10.3f %-10.3f\n",
              cens, m["o"], m["D"], m["o"] / m["D"],
              sprintf("%.4f/%.4f", m["g"], med["g"]),
              sprintf("%.4f/%.4f", m["tau"], med["tau"]), m["cov_own"], m["cov_tau"], m["rawci_covers_taustar"]))
}
ovec <- sapply(rows, function(m) m["o"]); gvec <- sapply(rows, function(m) m["g"]); tvec <- sapply(rows, function(m) m["tau"])
odvec <- ovec / sapply(rows, function(m) m["D"])
cat(sprintf("\nPortability check -- coefficient of variation across the censoring ladder:\n"))
cat(sprintf("  raw oracle(c) [population, no estimation noise]:        CV = %.3f\n", sd(ovec) / mean(ovec)))
cat(sprintf("  oracle(c)/D_hat(c) [population, no estimation noise]:   CV = %.3f\n", sd(odvec) / mean(odvec)))
cat(sprintf("  raw gamma_AR_hat(c) [estimated, weak-IV point]:         CV = %.3f\n", sd(gvec) / mean(gvec)))
cat(sprintf("  tau_hat = gamma_AR_hat/D_hat(c) [estimated point]:      CV = %.3f\n", sd(tvec) / mean(tvec)))
cat("DONE-TAU-PORTABLE\n")
