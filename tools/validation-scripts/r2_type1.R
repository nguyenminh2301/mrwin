suppressMessages(library(mrwin))
type1 <- function(M, N=2000, D=5, seed0=7000, inf="analytic"){
  rej <- rep(NA, M); dg <- rep(NA, M); weak <- rep(FALSE, M)
  for (m in seq_len(M)){
    s <- seed0 + m
    cfg <- mrwin_config(n_outcome=N, m_snps=15, alpha_x=c(0,0,0), gamma_direct=c(0,0,0), seed=s)
    dat <- mrwin_simulate(cfg, seed=s)
    ep <- mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
    gw <- mrwin_gwas(dat$true_betas, rep(0, length(dat$true_betas)))  # oracle weights, sigma_beta=0
    fit <- tryCatch(mrwin(endpoint=ep, genotype=dat$G, exposure=dat$X, gwas=gw,
      controls=mrwin_controls(n_strata=D, bootstrap=200, seed=s, inference=inf, run_sdpd=FALSE)),
      error=function(e) NULL)
    if (is.null(fit)) next
    ci <- fit$inference$ci95_delta
    rej[m] <- (ci[[1]] > 0) || (ci[[2]] < 0)
    dg[m] <- fit$point$delta_gls
    weak[m] <- any(vapply(fit$warnings, function(w) w$code=="weak_instrument", logical(1)))
  }
  ok <- !is.na(rej); n <- sum(ok)
  cat(sprintf("[%s] type-I: M=%d valid=%d  rejection rate=%.3f  (target 0.05, MC SE %.3f)\n",
      inf, M, n, mean(rej[ok]), sqrt(0.05*0.95/n)))
  cat(sprintf("        mean delta_GLS=%.4f  sd=%.4f  (bias check; should be ~0)  weak-instr frac=%.2f\n",
      mean(dg[ok]), sd(dg[ok]), mean(weak[ok])))
}
type1(500)
cat("DONE\n")
