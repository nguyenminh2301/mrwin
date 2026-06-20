suppressMessages(library(mrwin))
# truth = large-N pooled delta under valid IV
truth1<-function(N,seed){
  cfg<-mrwin_config(n_outcome=N,m_snps=80,alpha_x=c(-0.4,-0.4,-0.4),gamma_direct=c(0,0,0),seed=seed)
  dat<-mrwin_simulate(cfg,seed=seed); ep<-mrwin_endpoint(dat$time,dat$status,colnames(dat$time))
  gw<-mrwin_gwas(dat$true_betas,rep(0,length(dat$true_betas)))
  mrwin(endpoint=ep,genotype=dat$G,exposure=dat$X,gwas=gw,
    controls=mrwin_controls(n_strata=5,bootstrap=50,seed=seed,inference="analytic",run_sdpd=FALSE))$point$delta_gls
}
tr<-vapply(1:8,function(i) truth1(20000,90000+i),numeric(1)); delta_true<-mean(tr)
cat(sprintf("valid-IV truth (large-N delta) = %.4f (rep sd %.4f)\n",delta_true,sd(tr)))
M<-200; N<-6000
cov_bd<-cov_fi<-rep(NA,M)
for(m in seq_len(M)){
  s<-22000+m
  cfg<-mrwin_config(n_outcome=N,m_snps=80,alpha_x=c(-0.4,-0.4,-0.4),gamma_direct=c(0,0,0),seed=s)
  dat<-mrwin_simulate(cfg,seed=s); ep<-mrwin_endpoint(dat$time,dat$status,colnames(dat$time))
  gw<-mrwin_gwas(dat$true_betas,rep(0,length(dat$true_betas)))
  fit<-tryCatch(mrwin(endpoint=ep,genotype=dat$G,exposure=dat$X,gwas=gw,
    controls=mrwin_controls(n_strata=5,bootstrap=200,seed=s,inference="analytic",run_sdpd=FALSE)),error=function(e)NULL)
  if(is.null(fit)) next
  bd<-fit$inference$ci95_delta; fi<-fit$inference$ci95_delta_fieller
  cov_bd[m]<-(bd[[1]]<=delta_true)&&(delta_true<=bd[[2]])
  cov_fi[m]<-is.finite(fi[[1]])&&is.finite(fi[[2]])&&(fi[[1]]<=delta_true)&&(delta_true<=fi[[2]])
}
ok<-!is.na(cov_bd);n<-sum(ok)
cat(sprintf("coverage of truth=%.3f (M=%d valid=%d):\n",delta_true,M,n))
cat(sprintf("  bivariate-Delta = %.3f   Fieller = %.3f   (target 0.95, SE %.3f)\n",
    mean(cov_bd[ok]),mean(cov_fi[ok],na.rm=TRUE),sqrt(.95*.05/n)))
cat("DONE\n")
