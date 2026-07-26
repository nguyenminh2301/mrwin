#!/usr/bin/env Rscript
# Probe §8 for Paper 03 (dev/p3-weak-iv-robust-winmr.md): identification-robust
# (Anderson-Rubin) inference for win-ratio MR on the pairwise U-statistic moment.
# Three falsification probes against the interventional oracle, all reusing the
# shipped per-subject win-score primitive (mrwin_subject_win_loss_cpp).
#
#   8.1  AR vs Wald coverage across instrument strength (scalar Z = PRS).
#   8.2  degeneracy stress: naive chi2_1 AR coverage under weak instrument + heavy
#        censoring (does the U-statistic degeneracy bite?).
#   8.3  over-identified min-beta AR: a weak-IV-robust pleiotropy (Sargan/Q) test.
#
# RESULT (recorded in dev/p3-weak-iv-robust-winmr.md §8): AR is calibrated ~0.95
# uniformly while Wald over-covers (mis-calibrated) off the strong regime;
# degeneracy is real but MILD (naive-AR holds ~0.94 under extreme stress); the
# over-ID AR pleiotropy test is calibrated, correctly blind to InSIDE-violating
# pleiotropy, and powerful against InSIDE-satisfying pleiotropy. Public sim only.
suppressMessages(library(mrwin))
ncores <- max(1L, min(4L, parallel::detectCores()))
subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")
cfg0 <- mrwin_config(m_snps = 40L)

gen_out <- function(X, u, w, s, cn, uf, cfg) {
  n <- length(X); et <- matrix(NA_real_, n, 3); lw <- log(w)
  for (j in 1:3) { lp <- cfg$alpha_x[j]*X + cfg$nu_u[j]*u + cfg$gamma_direct[j]*s + lw
    et[, j] <- (-log(uf[, j]) / (cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull) }
  td <- et[,1]; th <- et[,2]; tr <- et[,3]
  list(time = cbind(pmin(td,cn), pmin(th,cn,td), pmin(tr,cn,td)),
       status = cbind(as.integer(td<=cn), as.integer(th<=cn & th<=td), as.integer(tr<=cn & tr<=td)))
}
gen_pop <- function(n, seed, mafs, betas, alpha_s, cens = NULL, gd = 0) {
  set.seed(seed); cfg <- cfg0; cfg$gamma_direct <- rep(gd, 3)
  graw <- sapply(mafs, function(p) rbinom(n, 2, p)); g <- scale(graw); storage.mode(g) <- "double"
  s <- drop(g %*% betas); u <- rnorm(n); w <- rgamma(n, 1/cfg$theta_f, scale = cfg$theta_f)
  cr <- if (is.null(cens)) cfg$censoring_rate else cens
  cn <- pmin(rexp(n, cr), cfg$max_follow_up); uf <- matrix(runif(n*3), n, 3)
  X <- alpha_s*s + cfg$alpha_u*u + rnorm(n); oc <- gen_out(X, u, w, s, cn, uf, cfg)
  list(time = oc$time, status = oc$status, gstd = g, X = X, s = as.numeric(s),
       u = u, w = w, cn = cn, uf = uf)
}
swin <- function(time, status, n) { cu <- subj(time, status, time, status, rep(1, n)); cu[,1] - cu[,2] }

# scalar AR statistic n*U(beta)^2 / (4*zeta1(beta));  Sigma = 4*Var(projection) (NOT /n)
ar_scalar <- function(time, status, X, Z) {
  n <- length(X); s_i <- swin(time, status, n)
  cz <- subj(time, status, time, status, Z); r_i <- cz[,1] - cz[,2]
  a <- 2*sum(Z*s_i)/(n*(n-1)); b <- 2*(n*sum(X*Z) - sum(X)*sum(Z))/(n*(n-1))
  gh <- (Z*s_i - r_i)/(n-1); SX <- sum(X); SZ <- sum(Z); SXZ <- sum(X*Z)
  gx <- ((n-1)*X*Z - X*(SZ-Z) - Z*(SX-X) + (SXZ - X*Z))/(n-1)
  list(a = a, b = b, c0 = var(gh), c1 = cov(gh, gx), c2 = var(gx), n = n)
}
ARstat <- function(m, beta) m$n*(m$a - m$b*beta)^2 / (4*(m$c0 - 2*beta*m$c1 + beta^2*m$c2))
ar_ci <- function(m, bstar, q = 3.841) {
  A <- m$n*m$b^2 - 4*q*m$c2; B <- -2*m$n*m$a*m$b + 8*q*m$c1; C <- m$n*m$a^2 - 4*q*m$c0
  disc <- B^2 - 4*A*C
  list(covers = ARstat(m, bstar) <= q, unb = (A <= 0),
       width = if (A > 0 && disc >= 0) sqrt(disc)/A else Inf)
}
wald_ci <- function(m, bstar) {
  bh <- m$a/m$b; z1 <- max(m$c0 - 2*bh*m$c1 + bh^2*m$c2, 0); se <- 2*sqrt(z1/m$n)/abs(m$b)
  list(covers = abs(bh - bstar) <= 1.96*se, width = 2*1.96*se)
}
# over-identified min-beta AR (vector Z = SNP dosages) -> chi2_{L-1}; Sigma = 4*Cov(proj)
overid <- function(p) {
  n <- length(p$X); G <- p$gstd; L <- ncol(G); s_i <- swin(p$time, p$status, n)
  Uh <- 2*as.numeric(crossprod(G, s_i))/(n*(n-1))
  Ux <- 2*(n*as.numeric(crossprod(G, p$X)) - colSums(G)*sum(p$X))/(n*(n-1))
  GH <- matrix(0, n, L)
  for (l in 1:L) { cz <- subj(p$time, p$status, p$time, p$status, G[, l]); r <- cz[,1] - cz[,2]; GH[, l] <- (G[,l]*s_i - r)/(n-1) }
  SX <- sum(p$X); GX <- matrix(0, n, L)
  for (l in 1:L) { Z <- G[,l]; SZ <- sum(Z); SXZ <- sum(p$X*Z); GX[, l] <- ((n-1)*p$X*Z - p$X*(SZ-Z) - Z*(SX-p$X) + (SXZ - p$X*Z))/(n-1) }
  b0 <- sum(Ux*Uh)/sum(Ux*Ux); Gp <- GH - b0*GX
  Sig <- 4*cov(Gp)                                              # 4*Cov(g)  (no /n)
  Si <- tryCatch(solve(Sig + diag(1e-8*mean(diag(Sig)), L)), error = function(e) NULL); if (is.null(Si)) return(NA)
  num <- as.numeric(t(Ux) %*% Si %*% Uh); den <- as.numeric(t(Ux) %*% Si %*% Ux)
  n*(as.numeric(t(Uh) %*% Si %*% Uh) - num^2/den)              # -> chi2_{L-1}
}

set.seed(7); mafs <- runif(40, cfg0$maf_low, cfg0$maf_high); betas <- rnorm(40, 0, cfg0$sigma_beta)
# beta* = net-benefit gradient: interventional oracle + moment-ratio plim (agree)
eps <- 0.25; po <- gen_pop(40000L, 99L, mafs, betas, alpha_s = 0.4)
Ap <- gen_out(po$X+eps, po$u,po$w,po$s,po$cn,po$uf, cfg0); Am <- gen_out(po$X-eps, po$u,po$w,po$s,po$cn,po$uf, cfg0)
nb <- mrwin_fast_pair_win_loss(Ap$time, Ap$status, Am$time, Am$status)
nb_oracle <- (as.numeric(nb[["wins"]]) - as.numeric(nb[["losses"]])) / (as.numeric(nrow(Ap$time))^2) / (2*eps)
bm <- unlist(parallel::mclapply(1:8, function(r) { p <- gen_pop(20000L, 300L+r, mafs, betas, 0.4)
  m <- ar_scalar(p$time, p$status, p$X, p$s); m$a/m$b }, mc.cores = ncores))
bstar <- mean(bm)
cat(sprintf("beta* net-benefit gradient: interventional=%.4f  moment-ratio=%.4f -> using %.4f\n\n", nb_oracle, bstar, bstar))

cat("== 8.1  AR vs Wald coverage across instrument strength (N=4000, R=300) ==\n")
for (as_ in c(0.10, 0.20, 0.30, 0.40)) {
  res <- parallel::mclapply(1:300, function(r) { p <- gen_pop(4000L, as.integer(as_*1e4)+r, mafs, betas, as_)
    m <- ar_scalar(p$time, p$status, p$X, p$s); list(ar = ar_ci(m, bstar), wd = wald_ci(m, bstar)) }, mc.cores = ncores)
  arc <- mean(sapply(res, function(x) x$ar$covers)); aru <- mean(sapply(res, function(x) x$ar$unb))
  arw <- median(Filter(is.finite, sapply(res, function(x) x$ar$width)))
  wdc <- mean(sapply(res, function(x) x$wd$covers)); wdw <- median(sapply(res, function(x) x$wd$width))
  cat(sprintf("  alpha_s=%.2f | AR cover=%.3f unbounded=%.2f medW=%.3f | Wald cover=%.3f medW=%.3f\n",
              as_, arc, aru, arw, wdc, wdw))
}

cat("\n== 8.2  degeneracy stress: naive chi2_1 AR coverage (R=400) ==\n")
for (cd in list(c(0.05,0.05), c(0.02,0.05), c(0.05,0.40), c(0.02,0.40), c(0.05,1.00))) {
  res <- unlist(parallel::mclapply(1:400, function(r) { p <- gen_pop(4000L, 7000L+r, mafs, betas, cd[1], cens = cd[2])
    ARstat(ar_scalar(p$time, p$status, p$X, p$s), bstar) }, mc.cores = ncores))
  cat(sprintf("  alpha_s=%.2f cens=%.2f | naive-AR coverage=%.3f\n", cd[1], cd[2], mean(res <= 3.841, na.rm = TRUE)))
}

cat("\n== 8.3  over-ID min-beta AR pleiotropy test (L=40, N=2500, R=120; chi2_39 crit=54.6) ==\n")
run83 <- function(gd, frac) unlist(parallel::mclapply(1:120, function(r) {
  if (frac > 0) { p <- gen_pop(2500L, 9000L+r, mafs, betas, 0.4); set.seed(9000L+r+99L)
    pi <- numeric(40); pi[sample(40, 12)] <- frac; lp_add <- drop(p$gstd %*% pi)
    et <- matrix(NA_real_, 2500, 3); lw <- log(p$w)
    for (j in 1:3) { lp <- cfg0$alpha_x[j]*p$X + cfg0$nu_u[j]*p$u + lp_add + lw; et[, j] <- (-log(p$uf[,j])/(cfg0$baseline_haz[j]*exp(lp)))^(1/cfg0$shape_weibull) }
    td <- et[,1]; th <- et[,2]; tr <- et[,3]
    p$time <- cbind(pmin(td,p$cn), pmin(th,p$cn,td), pmin(tr,p$cn,td))
    p$status <- cbind(as.integer(td<=p$cn), as.integer(th<=p$cn & th<=td), as.integer(tr<=p$cn & tr<=td))
  } else p <- gen_pop(2500L, 9000L+r, mafs, betas, 0.4, gd = gd)
  overid(p) }, mc.cores = ncores))
q39 <- qchisq(0.95, 39)
cat(sprintf("  null (no pleio)                  reject=%.3f (target 0.05)\n", mean(run83(0, 0) > q39, na.rm = TRUE)))
cat(sprintf("  InSIDE-violating (gd=0.3, prop)  reject=%.3f (correctly ~null: absorbed into beta-hat)\n", mean(run83(0.3, 0) > q39, na.rm = TRUE)))
cat(sprintf("  InSIDE-satisfying (tau=0.06,30%%) reject=%.3f (power)\n", mean(run83(0, 0.06) > q39, na.rm = TRUE)))
