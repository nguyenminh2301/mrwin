suppressMessages(library(mrwin))
cmp <- function(N, sb, D=5, B=4000, Bg=3000, seed=21){
  cfg <- mrwin_config(n_outcome=N, m_snps=15, seed=seed); dat <- mrwin_simulate(cfg, seed=seed)
  est <- mrwin_estimate(dat$time, dat$status, dat$G, dat$X, dat$true_betas, n_strata=D)
  ana <- mrwin_analytic_inference(est, dat$X, G=dat$G, beta_hat=dat$true_betas,
           sigma_beta=sb, n_strata=D, B_gwas=Bg, seed=seed)
  boot <- mrwin_multiplier_bootstrap(time=dat$time, status=dat$status, G=dat$G, X=dat$X,
           beta_hat=dat$true_betas, sigma_beta=sb, n_strata=D, B=B, seed=seed, adjustment="none")
  gfrac <- sum(diag(ana$cov_gwas))/sum(diag(ana$cov_u))   # share of variance from GWAS
  cat(sprintf("N=%5d sb=%.2f  gwasVarShare=%.2f  se_ana=%.4f se_boot=%.4f ratio=%.3f\n",
      N, sb, gfrac, ana$se_delta_gls, boot$se_delta_gls, ana$se_delta_gls/boot$se_delta_gls))
}
for (sb in c(0.05, 0.15, 0.30)) cmp(800, sb)
cmp(1500, 0.15)
cat("DONE\n")
