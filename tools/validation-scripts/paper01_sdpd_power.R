#!/usr/bin/env Rscript
# SDPD power against DETECTABLE (InSIDE-satisfying) pleiotropy — the supplement
# the main grid (paper01_mc_grid.R) could not provide. The simulator's
# `gamma_direct` enters the outcome as gamma_direct * s_true (s_true = G %*% betas),
# i.e. pleiotropy PROPORTIONAL to each SNP's instrument strength: that is the
# InSIDE-violating case, observationally equivalent to a stronger causal effect,
# which NO MR-Egger / SDPD test can detect (its ~null power there is correct, not
# a defect). To measure real SDPD power we inject per-SNP direct effects pi_k that
# are INDEPENDENT of betas (InSIDE-satisfying, directional) on a fraction of SNPs,
# faithfully regenerating the package outcome model otherwise, and sweep their
# magnitude. tau = 0 reproduces the null (cross-check vs the main grid's 0.047).
# Public package + simulator only; parallel reps.
suppressMessages(library(mrwin))
ncores <- max(1L, min(4L, parallel::detectCores()))

weibull_inv <- function(rate, shape, lp) {           # mirrors .mrwin_weibull_inv
  (-log(stats::runif(length(lp))) / (rate * exp(lp)))^(1 / shape)
}
# regenerate the 3-endpoint outcome with gamma_direct=0 + per-SNP pleiotropy G%*%pi
outcome_pleio <- function(cfg, X, U, W, G, pi) {
  n <- length(X); pl <- drop(G %*% pi); log_w <- log(W)
  et <- matrix(NA_real_, n, 3)
  for (j in 1:3) {
    lp <- cfg$alpha_x[j]*X + cfg$nu_u[j]*U + pl + log_w
    et[, j] <- weibull_inv(cfg$baseline_haz[j], cfg$shape_weibull, lp)
  }
  cn <- pmin(stats::rexp(n, cfg$censoring_rate), cfg$max_follow_up)
  td <- et[,1]; th <- et[,2]; tr <- et[,3]
  time <- cbind(pmin(td,cn), pmin(th,cn,td), pmin(tr,cn,td))
  status <- cbind(as.integer(td<=cn), as.integer(th<=cn & th<=td), as.integer(tr<=cn & tr<=td))
  colnames(time) <- colnames(status) <- c("death","hf","renal")
  list(time = time, status = status)
}

# one rep: invalid_frac of SNPs get a directional pi_k = tau (independent of beta)
rep_sdpd <- function(seed, N, m_snps, tau, invalid_frac = 0.3) {
  cfg <- mrwin_config(n_outcome = N, m_snps = m_snps, alpha_x = rep(-0.4, 3),
                      gamma_direct = c(0, 0, 0), seed = seed)
  dat <- mrwin_simulate(cfg, seed = seed)
  set.seed(seed + 7L)
  pi <- numeric(m_snps)
  inv <- sample.int(m_snps, max(1L, round(invalid_frac * m_snps)))
  pi[inv] <- tau                                   # directional, pi ⟂ betas (InSIDE)
  oc <- outcome_pleio(cfg, dat$X, dat$U, dat$W, dat$G, pi)
  ep <- mrwin_endpoint(oc$time, oc$status, colnames(oc$time))
  gw <- mrwin_gwas(dat$true_betas, rep(0, length(dat$true_betas)))
  ctrl <- mrwin_controls(n_strata = 5L, inference = "analytic", backend = "dense",
                         run_sdpd = TRUE, seed = seed)
  fit <- tryCatch(mrwin(endpoint = ep, genotype = dat$G, exposure = dat$X,
                        gwas = gw, controls = ctrl), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  c(rej = isTRUE(fit$sdpd$rejected), up = isTRUE(fit$sdpd$underpowered))
}

cat(sprintf("SDPD power vs InSIDE-satisfying per-SNP pleiotropy | cores=%d, N=4000, 30%% invalid\n", ncores))
for (tau in c(0.0, 0.03, 0.06, 0.10, 0.15)) {
  M <- 300L
  res <- parallel::mclapply(seq_len(M), function(m) rep_sdpd(40000L + as.integer(tau*1000) + m, 4000L, 100L, tau),
                            mc.cores = ncores)
  res <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
  rej <- mean(res[, "rej"]); n <- nrow(res); se <- sqrt(rej*(1-rej)/max(n,1))
  lbl <- if (tau == 0) "type-I" else "power "
  cat(sprintf("  tau=%.2f  SDPD reject (%s)=%.3f (SE %.3f)  underpowered=%.2f (n=%d)\n",
              tau, lbl, rej, se, mean(res[, "up"]), n))
  flush.console()
}
cat("DONE\n")
