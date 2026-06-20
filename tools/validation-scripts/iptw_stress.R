suppressMessages(library(mrwin))
# Adversarial: covariate Z strongly correlated with the PRS score -> informative
# IPTW weights -> maximal estimated-weights correction.
stress <- function(N, rho, D=5, B=1500, seed=51){
  cfg <- mrwin_config(n_outcome=N, m_snps=12, seed=seed); dat <- mrwin_simulate(cfg, seed=seed)
  score <- as.numeric(dat$G %*% dat$true_betas)
  Z <- cbind(z = as.numeric(scale(rho*scale(score) + sqrt(1-rho^2)*rnorm(N))))
  boot <- mrwin_multiplier_bootstrap(time=dat$time, status=dat$status, G=dat$G, X=dat$X,
            beta_hat=dat$true_betas, sigma_beta=0, n_strata=D, B=B, seed=seed,
            adjustment="ordinal_iptw", covariates=Z)
  om <- boot$adjustment$weights
  acov <- mrwin:::mrwin_analytic_covariance(boot$kernel, boot$strata, dat$X, boot$contrast_plan, weights=om)
  sig_a <- mrwin:::.mrwin_isg_covariance(acov, boot$point_log_theta, boot$point_delta_x)
  pa <- mrwin_gls_pool(boot$delta_isg, sig_a, shrink=TRUE)
  wvar <- sd(om)/mean(om)
  cat(sprintf("N=%4d rho=%.1f  weightCV=%.2f  se_wfix=%.4f se_boot=%.4f ratio=%.3f\n",
      N, rho, wvar, pa$se_delta_gls, boot$se_delta_gls, pa$se_delta_gls/boot$se_delta_gls))
}
for (rho in c(0.5, 0.8, 0.95)) stress(500, rho)
stress(300, 0.9)
cat("DONE\n")
