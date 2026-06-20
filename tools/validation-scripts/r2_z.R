suppressMessages(library(mrwin))
M<-200; N<-6000; msnp<-80; D<-5; seed0<-13000
dg<-se<-rep(NA,M)
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
  dg[m]<-fit$point$delta_gls; se[m]<-fit$inference$se_delta_gls
}
ok<-!is.na(dg); z<-abs(dg[ok]/se[ok])
cat(sprintf("N=%d msnp=%d valid=%d\n", N,msnp,sum(ok)))
cat(sprintf("  delta_gls: median=%.4f  MAD=%.4f  (robust spread)\n", median(dg[ok]), mad(dg[ok])))
cat(sprintf("  se:        median=%.4f  (vs MAD above)  -> typical se/spread = %.2f\n", median(se[ok]), median(se[ok])/mad(dg[ok])))
cat(sprintf("  |z|=|delta/se|: median=%.3f  q90=%.3f  q95=%.3f  max=%.3f\n", median(z), quantile(z,.9), quantile(z,.95), max(z)))
cat(sprintf("  rejections (|z|>1.96) = %d/%d = %.3f   (calibrated -> ~0.05)\n", sum(z>1.96), sum(ok), mean(z>1.96)))
cat("DONE\n")
