#!/usr/bin/env Rscript
# ============================================================================
# THE TRANSPORT-VELOCITY STRUCTURE, LIFTED TO THE HIERARCHICAL CENSORED OUTCOME.
#
# Scalar theory (p5-collapsible-transport-scalar.R) proved the interventional win
# gradient = 4*INT f^2 v = (density-overlap) x (transport velocity), so the win
# statistic is the causal effect SQUASHED by the outcome-law overlap INT f^2, and
# dividing by the (estimable) overlap recovers the collapsible effect.
#
# This script tests the two claims that this mechanism LIFTS to the K=3 censored
# hierarchical composite the mrwin programme actually uses, and UNIFIES two
# previously-separate phenomena:
#
#  (H1) CENSORING = SQUASH. Censoring turns informative comparisons into ties, so
#       it plays the role of "increasing sigma" in the scalar theory: it shrinks the
#       decidability overlap D = P(a random pair is decided, not a tie) and squashes
#       the win gradient. Ground truth = the near-uncensored gradient g0. Claim: the
#       decidability-corrected gradient g(c)/D(c) is far more censoring-STABLE than
#       the raw g(c) (recovers most of the squash), the hierarchical analogue of the
#       scalar collapsibility correction.
#
#  (H2) THE SQUASH FACTOR IS THE PAPER-03 DEGENERACY FACTOR. The same decidability
#       overlap that squashes the win gradient (H1, non-collapsibility, this work)
#       is what drives the win-MR U-statistic's first-projection variance zeta1
#       (Paper 03's degeneracy boundary) to zero. So the object that makes the win
#       effect non-collapsible and the object that makes Paper-03 inference degenerate
#       are the SAME functional: max censoring -> min decidability -> simultaneously
#       max non-collapsibility squash AND the zeta1->0 degenerate boundary. Claim:
#       zeta1_proxy(c) and D(c) fall together (strongly rank-correlated), and the
#       win gradient scales with D.
#
# Ground truth: the mrwin interventional oracle do(X+-eps). Public sim only.
# ============================================================================
suppressMessages(library(mrwin))
ncores <- max(1L, min(4L, parallel::detectCores()))
subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")

mk_cfg <- function(maxfu) { c <- mrwin_config(m_snps = 40L); c$max_follow_up <- maxfu; c }
gen_out <- function(X, u, w, cn, uf, cfg) {
  n <- length(X); et <- matrix(NA_real_, n, 3); lw <- log(w)
  for (j in 1:3) { lp <- cfg$alpha_x[j]*X + cfg$nu_u[j]*u + lw
    et[, j] <- (-log(uf[, j]) / (cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull) }
  td <- et[,1]; th <- et[,2]; tr <- et[,3]
  list(time = cbind(pmin(td,cn), pmin(th,cn,td), pmin(tr,cn,td)),
       status = cbind(as.integer(td<=cn), as.integer(th<=cn & th<=td), as.integer(tr<=cn & tr<=td)))
}
gen_pop <- function(n, seed, mafs, betas, alpha_s, cens, cfg) {
  set.seed(seed)
  g <- scale(sapply(mafs, function(p) rbinom(n, 2, p))); storage.mode(g) <- "double"
  s <- drop(g %*% betas); u <- rnorm(n); w <- rgamma(n, 1/cfg$theta_f, scale = cfg$theta_f)
  cn <- pmin(rexp(n, cens), cfg$max_follow_up); uf <- matrix(runif(n*3), n, 3)
  X <- alpha_s*s + cfg$alpha_u*u + rnorm(n); oc <- gen_out(X, u, w, cn, uf, cfg)
  list(time = oc$time, status = oc$status, X = X, s = as.numeric(s), u = u, w = w, cn = cn, uf = uf, cfg = cfg)
}
# oracle win gradient (codebase convention: /(2 eps)); also returns decided fraction of the +/- comparison
oracle_grad <- function(d, eps = 0.25) {
  cfg <- d$cfg; np <- length(d$X)
  Bp <- gen_out(d$X + eps, d$u, d$w, d$cn, d$uf, cfg)
  Bm <- gen_out(d$X - eps, d$u, d$w, d$cn, d$uf, cfg)
  v <- mrwin_fast_pair_win_loss(Bp$time, Bp$status, Bm$time, Bm$status)
  wins <- as.numeric(v[["wins"]]); loss <- as.numeric(v[["losses"]]); tot <- as.numeric(np)^2
  list(grad = (wins - loss) / tot / (2 * eps), decided = (wins + loss) / tot)
}
# decidability of the BASE population (self-pairs) + outcome-side zeta1 proxy Var(win score)
base_stats <- function(d) {
  n <- length(d$X)
  cu <- subj(d$time, d$status, d$time, d$status, rep(1, n))
  wins <- sum(cu[,1]); loss <- sum(cu[,2])
  Dbase <- (wins + loss) / (as.numeric(n)^2)          # P(pair decided) at base pop
  ws <- (cu[,1] - cu[,2]) / (n - 1)                    # per-subject net win score
  list(D = Dbase, var_w = stats::var(ws))
}

cfgH <- mk_cfg(maxfu = 300)                            # long follow-up so low-rate censoring is near-clean
set.seed(11); mafs <- runif(40, cfgH$maf_low, cfgH$maf_high); betas <- rnorm(40, 0, cfgH$sigma_beta)
N <- 6000L; R <- 24L; alpha_s <- 0.4
cens_grid <- c(0.002, 0.05, 0.15, 0.35, 0.8, 1.6)

cat("== (H1)+(H2) censoring sweep: win gradient, decidability D, zeta1-proxy Var(w) ==\n")
cat(sprintf("K=3 hierarchical, N=%d, R=%d, alpha_s=%.2f; ground truth g0 = gradient at lowest censoring\n\n", N, R, alpha_s))
cat(sprintf("%-8s %-16s %-10s %-12s %-16s %-14s\n",
            "cens", "win grad g(c)", "decided D", "Var(w)", "corrected g/D", "g/g0 vs D/D0"))
res <- list()
for (cr in cens_grid) {
  rr <- do.call(rbind, parallel::mclapply(1:R, function(r) {
    d <- gen_pop(N, 2000L + r, mafs, betas, alpha_s, cr, cfgH)
    og <- oracle_grad(d); bs <- base_stats(d)
    c(g = og$grad, D = bs$D, vw = bs$var_w)
  }, mc.cores = ncores))
  m <- colMeans(rr); se <- apply(rr, 2, sd) / sqrt(R)
  res[[as.character(cr)]] <- m
  cat(sprintf("%-8.3f %-16s %-10.4f %-12.3e %-16s\n",
              cr, sprintf("%.4f(%.4f)", m["g"], se["g"]), m["D"], m["vw"],
              sprintf("%.4f", m["g"] / m["D"])))
}
g0 <- res[[1]]["g"]; D0 <- res[[1]]["D"]; vw0 <- res[[1]]["vw"]
cat(sprintf("\n%-8s %-14s %-14s %-14s\n", "cens", "g/g0", "D/D0", "Var(w)/Var(w)0"))
for (cr in cens_grid) { m <- res[[as.character(cr)]]
  cat(sprintf("%-8.3f %-14.4f %-14.4f %-14.4f\n", cr, m["g"]/g0, m["D"]/D0, m["vw"]/vw0)) }

# (H1) stability: SD of raw g vs SD of corrected g/D across the censoring sweep
gvec <- sapply(res, function(m) m["g"]); Dvec <- sapply(res, function(m) m["D"]); vwvec <- sapply(res, function(m) m["vw"])
cat(sprintf("\n(H1) coefficient of variation across censoring sweep:  raw g = %.3f   corrected g/D = %.3f\n",
            sd(gvec)/mean(gvec), sd(gvec/Dvec)/mean(gvec/Dvec)))
# (H2) do decidability and zeta1-proxy fall together?
cat(sprintf("(H2) Spearman corr( decidability D , zeta1-proxy Var(w) ) across sweep = %.4f\n",
            cor(Dvec, vwvec, method = "spearman")))
cat(sprintf("(H2) Var(w)/D across sweep (should be ~constant if zeta1 ∝ D): %s\n",
            paste(sprintf("%.4f", vwvec / Dvec), collapse = " ")))
cat("DONE-HIER\n")
