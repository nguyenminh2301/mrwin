suppressMessages(library(mrwin))
measure <- function(N, D=5, B=1500, seed=31){
  cfg <- mrwin_config(n_outcome=N, m_snps=12, seed=seed); dat <- mrwin_simulate(cfg, seed=seed)
  Z <- cbind(u = as.numeric(scale(dat$U)))
  boot <- mrwin_multiplier_bootstrap(time=dat$time, status=dat$status, G=dat$G, X=dat$X,
            beta_hat=dat$true_betas, sigma_beta=0, n_strata=D, B=B, seed=seed,
            adjustment="ordinal_iptw", covariates=Z)
  om <- boot$adjustment$weights            # point IPTW weights (xi=1)
  acov <- mrwin:::mrwin_analytic_covariance(boot$kernel, boot$strata, dat$X,
            boot$contrast_plan, weights=om)
  bcov <- boot$cov_u
  # downstream se from each
  sig_a <- mrwin:::.mrwin_isg_covariance(acov, boot$point_log_theta, boot$point_delta_x)
  pa <- mrwin_gls_pool(boot$delta_isg, sig_a, shrink=TRUE)
  fro <- sqrt(sum((acov-bcov)^2))/sqrt(sum(bcov^2))
  dgl <- mean(abs(diag(acov)-diag(bcov))/pmax(abs(diag(bcov)),1e-12))
  cat(sprintf("N=%4d  relFrob=%.3f  diagRelErr=%.3f  se_wfix=%.4f se_boot=%.4f ratio=%.3f\n",
      N, fro, dgl, pa$se_delta_gls, boot$se_delta_gls, pa$se_delta_gls/boot$se_delta_gls))
}
for (N in c(400, 800)) measure(N)
cat("DONE\n")
