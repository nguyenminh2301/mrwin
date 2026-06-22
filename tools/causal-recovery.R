#!/usr/bin/env Rscript
# Ground-truth causal-recovery validation for the mrwin methods write-up.
#
# Question answered: does DS-CWR recover the *causal* win-odds gradient that a
# randomized experiment would measure -- and does the naive (non-instrumental)
# alternative fail? We build an interventional oracle gamma* by simulating
# do(X = x_lo) vs do(X = x_hi) on a paired population (same confounder/frailty
# draws, exposure set exogenously), then show across growing N that DS-CWR
# (PRS as instrument) -> gamma* while the naive estimator (observed exposure as
# the stratifier) converges to a confounded, biased value.
#
# Entirely from the public package + simulator; no confidential data.
#   Rscript tools/causal-recovery.R          # writes manuscript/figures/fig-recovery.pdf
#   MRWIN_FIG_DIR=some/dir Rscript tools/causal-recovery.R
suppressMessages(library(mrwin))
set.seed(1)
outdir <- Sys.getenv("MRWIN_FIG_DIR", unset = file.path("manuscript", "figures"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- DGP replica that exposes latent draws so X can be set exogenously -------
gen_outcomes <- function(X, u, w, s, cn, uf, cfg) {
  n <- length(X); et <- matrix(NA_real_, n, 3); lw <- log(w)
  for (j in 1:3) {
    lp <- cfg$alpha_x[j]*X + cfg$nu_u[j]*u + cfg$gamma_direct[j]*s + lw
    et[, j] <- (-log(uf[, j]) / (cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)
  }
  td <- et[,1]; th <- et[,2]; tr <- et[,3]
  time <- cbind(pmin(td,cn), pmin(th,cn,td), pmin(tr,cn,td))
  st <- cbind(as.integer(td<=cn), as.integer(th<=cn & th<=td), as.integer(tr<=cn & tr<=td))
  colnames(time) <- colnames(st) <- c("death","hf","renal"); list(time=time, status=st)
}
gen_pop <- function(n, cfg, seed) {
  set.seed(seed); m <- cfg$m_snps
  mafs <- runif(m, cfg$maf_low, cfg$maf_high)
  g <- scale(sapply(mafs, function(p) rbinom(n,2,p))); storage.mode(g) <- "double"
  b <- rnorm(m, 0, cfg$sigma_beta); s <- drop(g %*% b)
  u <- rnorm(n); w <- rgamma(n, 1/cfg$theta_f, scale=cfg$theta_f)
  cn <- pmin(rexp(n, cfg$censoring_rate), cfg$max_follow_up); uf <- matrix(runif(n*3), n, 3)
  X <- cfg$alpha_s*s + cfg$alpha_u*u + rnorm(n)
  oc <- gen_outcomes(X, u, w, s, cn, uf, cfg)
  list(time=oc$time, status=oc$status, G=g, X=X, beta=b, s=s, u=u, w=w, cn=cn, uf=uf)
}
WL <- function(th,sh,tl,sl) {
  v <- mrwin_fast_pair_win_loss(th,sh,tl,sl)
  log(as.numeric(v[["wins"]]) / as.numeric(v[["losses"]]))
}
qs <- function(z, D) findInterval(z, quantile(z, seq(0,1,length.out=D+1)),
                                  rightmost.closed=TRUE, all.inside=TRUE)
delta_gls <- function(p, G, beta, se, D)
  mrwin(endpoint = mrwin_endpoint(p$time, p$status, colnames(p$time)),
        genotype = G, exposure = p$X, gwas = mrwin_gwas(beta, se),
        controls = mrwin_controls(n_strata=D, bootstrap=2L, seed=1L,
                                  backend="fast", run_sdpd=FALSE))$point$delta_gls

cfg <- mrwin_config(m_snps = 40L)   # alpha_x=-0.4 (protective), alpha_u=0.8 (confounding)
D <- 5L

## ---- interventional oracle: gamma* = causal log-win-odds gradient -----------
ref <- gen_pop(20000L, cfg, 1L); sref <- qs(ref$s, D)
xlo <- mean(ref$X[sref == 1]); xhi <- mean(ref$X[sref == D])
po <- gen_pop(200000L, cfg, 99L)
A <- gen_outcomes(rep(xlo, 200000L), po$u, po$w, po$s, po$cn, po$uf, cfg)  # do(X=xlo)
B <- gen_outcomes(rep(xhi, 200000L), po$u, po$w, po$s, po$cn, po$uf, cfg)  # do(X=xhi)
gstar <- WL(B$time, B$status, A$time, A$status) / (xhi - xlo)
cat(sprintf("oracle gamma* = %.4f  (causal, over X-span %.3f)\n\n", gstar, xhi - xlo))

## ---- recovery across N: DS-CWR (PRS instr) vs naive (observed-X instr) ------
Ns <- c(4000, 8000, 16000, 32000, 64000); R <- 30L
ds_m <- ds_se <- nv_m <- nv_se <- numeric(length(Ns))
for (i in seq_along(Ns)) {
  N <- Ns[i]; ds <- nv <- numeric(R)
  for (r in 1:R) {
    p <- gen_pop(as.integer(N), cfg, 1000L + i*1000L + r)
    ds[r] <- delta_gls(p, p$G, p$beta, rep(0.01, length(p$beta)), D)  # instrument = PRS
    nv[r] <- delta_gls(p, matrix(p$X, N, 1), 1, 0, D)                 # instrument = observed X
  }
  ds_m[i] <- mean(ds); ds_se[i] <- sd(ds)/sqrt(R)
  nv_m[i] <- mean(nv); nv_se[i] <- sd(nv)/sqrt(R)
  cat(sprintf("N=%6d  DS-CWR %.3f (se %.3f) bias %+.3f | naive %.3f (se %.3f) bias %+.3f\n",
              N, ds_m[i], ds_se[i], ds_m[i]-gstar, nv_m[i], nv_se[i], nv_m[i]-gstar))
}

## ---- figure: estimates +/- 2 se vs N, with the causal target gamma* --------
pdf(file.path(outdir, "fig-recovery.pdf"), width = 6, height = 4.2)
yl <- range(c(ds_m-2*ds_se, ds_m+2*ds_se, nv_m-2*nv_se, gstar, 0))
plot(Ns, ds_m, log = "x", type = "b", pch = 19, col = "steelblue", ylim = yl,
     xlab = "sample size N", ylab = expression(hat(delta)[GLS]),
     main = "Causal recovery: DS-CWR vs naive")
arrows(Ns, ds_m-2*ds_se, Ns, ds_m+2*ds_se, angle=90, code=3, length=0.03, col="steelblue")
lines(Ns, nv_m, type = "b", pch = 17, col = "firebrick")
arrows(Ns, nv_m-2*nv_se, Ns, nv_m+2*nv_se, angle=90, code=3, length=0.03, col="firebrick")
abline(h = gstar, lty = 2, col = "grey30")
text(Ns[1], gstar, expression(gamma*"* (causal truth)"), pos = 3, col = "grey30", cex = 0.9)
legend("topright", bty = "n", pch = c(19,17), col = c("steelblue","firebrick"),
       legend = c("DS-CWR (PRS instrument) -> gamma*",
                  "naive (observed X) -> confounded"))
dev.off()
cat(sprintf("\nfigure written to %s\n", normalizePath(file.path(outdir, "fig-recovery.pdf"))))
