suppressMessages(library(mrwin))
M<-150; N<-6000; msnp<-80; D<-5; seed0<-42000
dg<-se_s<-se_ns<-rep(NA,M)
for(m in seq_len(M)){
  s<-seed0+m
  cfg<-mrwin_config(n_outcome=N, m_snps=msnp, alpha_x=c(0,0,0), gamma_direct=c(0,0,0), seed=s)
  dat<-mrwin_simulate(cfg, seed=s)
  ep<-mrwin_endpoint(dat$time, dat$status, colnames(dat$time))
  gw<-mrwin_gwas(dat$true_betas, rep(0, length(dat$true_betas)))
  fit<-tryCatch(mrwin(endpoint=ep, genotype=dat$G, exposure=dat$X, gwas=gw,
    controls=mrwin_controls(n_strata=D, bootstrap=200, seed=s, inference="analytic", run_sdpd=FALSE)),
    error=function(e) NULL)
  if(is.null(fit)) next
  dg[m]<-fit$point$delta_gls; se_s[m]<-fit$inference$se_delta_gls
  # recompute se WITHOUT Ledoit-Wolf shrinkage
  ns<-tryCatch(mrwin_gls_pool(fit$bootstrap$delta_isg, fit$bootstrap$sigma_isg, shrink=FALSE), error=function(e) NULL)
  if(!is.null(ns)) se_ns[m]<-ns$se_delta_gls
}
ok<-!is.na(dg)
cat(sprintf("N=%d msnp=%d M=%d: empirical sd(delta_gls)=%.4f\n", N,msnp,sum(ok),sd(dg[ok])))
cat(sprintf("  mean reported se (Ledoit-Wolf shrink) = %.4f  ->  ratio se/sd = %.2f\n", mean(se_s[ok]), mean(se_s[ok])/sd(dg[ok])))
cat(sprintf("  mean reported se (NO shrink)          = %.4f  ->  ratio se/sd = %.2f\n", mean(se_ns[ok],na.rm=TRUE), mean(se_ns[ok],na.rm=TRUE)/sd(dg[ok])))
cat("DONE\n")
