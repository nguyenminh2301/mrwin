suppressMessages(library(mrwin))
type1 <- function(M, N, msnp, D, seed0=5000){
  rej<-rep(NA,M); dg<-rep(NA,M); weak<-rep(NA,M)
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
    ci<-fit$inference$ci95_delta
    rej[m]<-(ci[[1]]>0)||(ci[[2]]<0); dg[m]<-fit$point$delta_gls
    weak[m]<-any(vapply(fit$warnings, function(w) w$code=="weak_instrument", logical(1)))
  }
  ok<-!is.na(rej); n<-sum(ok)
  cat(sprintf("N=%d msnp=%d D=%d : valid=%d  type-I=%.3f (SE %.3f)  meanDg=%.3f medianDg=%.3f sdDg=%.3f  weak=%.2f\n",
      N,msnp,D,n,mean(rej[ok]),sqrt(.05*.95/n),mean(dg[ok]),median(dg[ok]),sd(dg[ok]),mean(weak[ok])))
}
type1(250, 5000, 60, 5)
type1(250, 8000, 100, 5)
cat("DONE\n")
