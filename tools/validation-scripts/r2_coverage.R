suppressMessages(library(mrwin))
# 1) Estimate the estimand "truth" = large-N limit of delta_GLS under valid IV.
truth_fit <- function(N, seed){
  cfg <- mrwin_config(n_outcome=N, m_snps=15, alpha_x=c(-0.4,-0.4,-0.4), gamma_direct=c(0,0,0), seed=seed)
  dat <- mrwin_simulate(cfg, seed=seed)
  ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
  gw <- mrwin_gwas(dat$true_betas, rep(0, length(dat$true_betas)))
  fit <- mrwin(endpoint=ep, genotype=dat$G, exposure=dat$X, gwas=gw,
    controls=mrwin_controls(n_strata=5, bootstrap=50, seed=seed, inference="analytic", run_sdpd=FALSE))
  fit$point$delta_gls
}
truths <- vapply(1:8, function(i) truth_fit(20000, 90000+i), numeric(1))
delta_true <- mean(truths)
cat(sprintf("estimand truth (large-N delta_GLS) = %.4f  (sd across reps %.4f)\n", delta_true, sd(truths)))

cover <- function(M, N=2000, D=5, seed0=8000){
  cov <- rep(NA, M)
  for (m in seq_len(M)){
    s <- seed0+m
    cfg <- mrwin_config(n_outcome=N, m_snps=15, alpha_x=c(-0.4,-0.4,-0.4), gamma_direct=c(0,0,0), seed=s)
    dat <- mrwin_simulate(cfg, seed=s)
    ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
    gw <- mrwin_gwas(dat$true_betas, rep(0, length(dat$true_betas)))
    fit <- tryCatch(mrwin(endpoint=ep, genotype=dat$G, exposure=dat$X, gwas=gw,
      controls=mrwin_controls(n_strata=D, bootstrap=200, seed=s, inference="analytic", run_sdpd=FALSE)),
      error=function(e) NULL)
    if (is.null(fit)) next
    ci <- fit$inference$ci95_delta
    cov[m] <- (ci[[1]] <= delta_true) && (delta_true <= ci[[2]])
  }
  ok <- !is.na(cov); n <- sum(ok)
  cat(sprintf("coverage of true delta=%.3f: M=%d valid=%d  coverage=%.3f (target 0.95, MC SE %.3f)\n",
      delta_true, M, n, mean(cov[ok]), sqrt(0.95*0.05/n)))
}
cover(400)
cat("DONE\n")
