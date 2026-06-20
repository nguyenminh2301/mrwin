suppressMessages(library(mrwin))
cmp <- function(N, D=5, B=4000, seed=11){
  cfg <- mrwin_config(n_outcome=N, m_snps=15, seed=seed); dat <- mrwin_simulate(cfg, seed=seed)
  kernel <- mrwin_kernel(dat$time, dat$status)
  boot <- mrwin_sparse_bootstrap(time=dat$time, status=dat$status, G=dat$G, X=dat$X,
            beta_hat=dat$true_betas, sigma_beta=0, n_strata=D, B=B, seed=seed,
            adjustment="none", fast=TRUE)
  cp <- boot$contrast_plan; strata <- boot$strata
  acov <- mrwin:::mrwin_analytic_covariance(kernel, strata, dat$X, cp)
  bcov <- boot$cov_u
  fro <- sqrt(sum((acov-bcov)^2))/sqrt(sum(bcov^2))
  dgl <- mean(abs(diag(acov)-diag(bcov))/pmax(abs(diag(bcov)),1e-12))
  sig_a <- mrwin:::.mrwin_isg_covariance(acov, boot$point_log_theta, boot$point_delta_x)
  pooled_a <- mrwin_gls_pool(boot$delta_isg, sig_a, shrink=TRUE)
  cat(sprintf("N=%5d B=%d  relFrob=%.3f  diagRelErr=%.3f  se_analytic=%.4f se_boot=%.4f ratio=%.3f\n",
      N, B, fro, dgl, pooled_a$se_delta_gls, boot$se_delta_gls, pooled_a$se_delta_gls/boot$se_delta_gls))
}
for (N in c(500, 1500, 4000)) cmp(N)
cat("DONE\n")
