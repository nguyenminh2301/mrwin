suppressMessages(library(mrwin))
M<-200; N<-6000; msnp<-80; D<-5; seed0<-13000   # same seeds as z-diagnostic
rej_bd<-rej_fi<-fi_unb<-rep(NA,M)
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
  bd<-fit$inference$ci95_delta
  rej_bd[m]<-(bd[[1]]>0)||(bd[[2]]<0)
  fi<-fit$inference$ci95_delta_fieller; fi_unb[m]<-fit$bootstrap$fieller_unbounded
  rej_fi[m]<-is.finite(fi[[1]])&&is.finite(fi[[2]])&&((fi[[1]]>0)||(fi[[2]]<0))
}
ok<-!is.na(rej_bd); n<-sum(ok)
cat(sprintf("N=%d msnp=%d valid=%d\n",N,msnp,n))
cat(sprintf("  bivariate-Delta type-I = %.3f (SE %.3f)\n", mean(rej_bd[ok]), sqrt(.05*.95/n)))
cat(sprintf("  Fieller         type-I = %.3f   (Fieller unbounded frac = %.2f)\n", mean(rej_fi[ok],na.rm=TRUE), mean(fi_unb[ok],na.rm=TRUE)))
cat("DONE\n")
