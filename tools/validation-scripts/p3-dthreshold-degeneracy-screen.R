#!/usr/bin/env Rscript
# ============================================================================
# CROSS-CUTTING PROBE #3: a cheap decidability (D) pre-screen for Paper 03's
# degeneracy flag, motivated by the transport-collapsibility finding that the
# SAME decidability functional D = P(pair decided, not tied) that squashes the
# win effect (non-collapsibility) also tracks the AR test's degenerate
# boundary (dev/findings-transport-collapsibility.md, Claim: Spearman(D,Var(w))
# = 1.0000 in the hierarchical censoring sweep).
#
# mrwin_ar_onesample()'s degeneracy flag currently costs boot_reps (200 default)
# full re-computations of .mrwin_ar1_blocks() -- each an O(n log n / n^2)
# per-subject kernel call -- to bootstrap a CI for c^2 = n*zeta1(beta_hat),
# because "a point estimate alone is unreliable near the boundary" (the
# existing docstring's stated reason for not just thresholding c2_hat).
#
# QUESTION THIS SCRIPT ANSWERS: can a functional that is available FOR FREE
# from the *already-computed* .mrwin_ar1_blocks() call -- either the point
# estimate c2_hat itself, or the raw per-subject decidability D_hat = mean
# fraction of decided pairs (from the `cw` win/loss counts already computed
# inside .mrwin_ar1_blocks, zero extra cost) -- serve as a conservative
# pre-screen: "when the screen clears a validated threshold, skip the
# bootstrap and report degenerate=FALSE; only pay for the bootstrap when the
# screen is ambiguous"? This is a two-stage design, not a replacement: the
# bootstrap remains the ground truth in the region the screen flags.
#
# We need the ONE-SIDED guarantee: screen-clears-threshold must essentially
# NEVER co-occur with the (expensive, ground-truth) bootstrap flagging
# degenerate=TRUE, across a wide censoring/N/instrument-strength grid.
#
# Ground truth = mrwin_ar_onesample()'s own bootstrap flag (boot_reps=200,
# unchanged, validated machinery -- lesson G). Public sim + package only.
# ============================================================================
suppressMessages(library(mrwin))
ncores <- max(1L, min(4L, parallel::detectCores()))
blocks <- getFromNamespace(".mrwin_ar1_blocks", "mrwin")
subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")

cfg0 <- mrwin_config(m_snps = 1L, n_outcome = 1200L)
set.seed(3); maf <- runif(1, cfg0$maf_low, cfg0$maf_high); beta1 <- rnorm(1, 0, cfg0$sigma_beta)

gen_one <- function(n, seed, alpha_s, cens) {
  set.seed(seed)
  g <- scale(rbinom(n, 2, maf)); storage.mode(g) <- "double"
  Z <- as.numeric(g) * beta1
  u <- rnorm(n); w <- rgamma(n, 1 / cfg0$theta_f, scale = cfg0$theta_f)
  cn <- pmin(rexp(n, cens), cfg0$max_follow_up); uf <- matrix(runif(n * 3), n, 3)
  X <- alpha_s * Z + cfg0$alpha_u * u + rnorm(n)
  et <- matrix(NA_real_, n, 3); lw <- log(w)
  for (j in 1:3) { lp <- cfg0$alpha_x[j] * X + cfg0$nu_u[j] * u + lw
    et[, j] <- (-log(uf[, j]) / (cfg0$baseline_haz[j] * exp(lp)))^(1 / cfg0$shape_weibull) }
  td <- et[, 1]; th <- et[, 2]; tr <- et[, 3]
  list(time = cbind(pmin(td, cn), pmin(th, cn, td), pmin(tr, cn, td)),
       status = cbind(as.integer(td <= cn), as.integer(th <= cn & th <= td), as.integer(tr <= cn & tr <= td)),
       X = X, Z = Z)
}

# D_hat: cheap, per-subject decided-fraction from the SAME cw the blocks fn computes internally.
d_hat_cheap <- function(time, status) {
  n <- nrow(time)
  cw <- subj(time, status, time, status, rep(1, n))
  mean((cw[, 1] + cw[, 2]) / (n - 1))
}

one_rep <- function(n, seed, alpha_s, cens) {
  d <- gen_one(n, seed, alpha_s, cens)
  t0 <- proc.time()[["elapsed"]]
  m <- blocks(d$time, d$status, d$X, d$Z)
  gamma <- m$a / m$b
  c2_hat_cheap <- max(m$n * (m$c0 - 2 * gamma * m$c1 + gamma^2 * m$c2), 0)
  D_hat <- d_hat_cheap(d$time, d$status)
  t_cheap <- proc.time()[["elapsed"]] - t0
  t1 <- proc.time()[["elapsed"]]
  full <- mrwin_ar_onesample(d$time, d$status, d$X, d$Z, degeneracy_check = TRUE,
                              boot_reps = 200L, seed = seed)
  t_full <- proc.time()[["elapsed"]] - t1
  c(D_hat = D_hat, c2_hat_cheap = c2_hat_cheap, c2_hat_full = full$c2_hat,
    c2_lo = full$c2_boot_ci[1], degenerate = as.numeric(full$degenerate),
    t_cheap = t_cheap, t_full = t_full)
}

cens_grid <- c(0.01, 0.05, 0.15, 0.4, 1.0, 2.5, 6.0)
N_grid <- c(800L, 2500L)
R <- 30L
cat("== D_hat / c2_hat (cheap) vs bootstrap degenerate flag (expensive, ground truth) ==\n")
cat(sprintf("%-6s %-6s %-10s %-12s %-12s %-10s %-12s\n",
            "N", "cens", "D_hat", "c2_hat", "c2_boot_lo", "P(degen)", "speedup"))
all_rows <- list()
for (n in N_grid) for (cens in cens_grid) {
  rr <- do.call(rbind, parallel::mclapply(seq_len(R), function(r) {
    tryCatch(one_rep(n, 10000L * n + 7L * round(cens * 100) + r, alpha_s = 0.05, cens = cens),
             error = function(e) NULL)
  }, mc.cores = ncores))
  rr <- rr[stats::complete.cases(rr), , drop = FALSE]
  all_rows[[paste(n, cens)]] <- data.frame(N = n, cens = cens, rr)
  m <- colMeans(rr)
  cat(sprintf("%-6d %-6.2f %-10.4f %-12.4f %-12.4f %-10.3f %-12.1f\n",
              n, cens, m["D_hat"], m["c2_hat_cheap"], m["c2_lo"], m["degenerate"],
              m["t_full"] / m["t_cheap"]))
}
dat <- do.call(rbind, all_rows)
cat(sprintf("\nTotal rows: %d\n", nrow(dat)))
cat(sprintf("Spearman(D_hat, c2_hat_cheap) = %.4f\n", cor(dat$D_hat, dat$c2_hat_cheap, method = "spearman")))
cat(sprintf("Spearman(D_hat, c2_boot_lo)   = %.4f\n", cor(dat$D_hat, dat$c2_lo, method = "spearman")))
cat(sprintf("Spearman(D_hat, degenerate)   = %.4f\n", cor(dat$D_hat, dat$degenerate, method = "spearman")))

cat("\n== Screen calibration: does D_hat > tau ever co-occur with degenerate==1? ==\n")
cat(sprintf("%-8s %-14s %-16s %-10s\n", "tau", "n_above_tau", "n_degen_above", "coverage(skip%)"))
for (tau in c(0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30)) {
  above <- dat$D_hat > tau
  n_above <- sum(above); n_bad <- sum(above & dat$degenerate == 1)
  cat(sprintf("%-8.2f %-14d %-16d %-10.1f\n", tau, n_above, n_bad, 100 * n_above / nrow(dat)))
}
cat("\n== Same screen using c2_hat_cheap (point estimate, no bootstrap) instead of D_hat ==\n")
cat(sprintf("%-8s %-14s %-16s %-10s\n", "tau", "n_above_tau", "n_degen_above", "coverage(skip%)"))
for (tau in c(1, 2, 3, 5, 8, 12, 20)) {
  above <- dat$c2_hat_cheap > tau
  n_above <- sum(above); n_bad <- sum(above & dat$degenerate == 1)
  cat(sprintf("%-8.1f %-14d %-16d %-10.1f\n", tau, n_above, n_bad, 100 * n_above / nrow(dat)))
}
cat("\nPer-row scatter (D_hat, c2_hat_cheap, degenerate) for the borderline region (D_hat<0.15):\n")
sub <- dat[dat$D_hat < 0.15, c("N", "cens", "D_hat", "c2_hat_cheap", "c2_lo", "degenerate")]
print(sub[order(sub$D_hat), ], row.names = FALSE)
cat("DONE-DTHRESHOLD\n")
