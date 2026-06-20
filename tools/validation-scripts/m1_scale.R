suppressMessages(library(mrwin))
CISG<-mrwin:::mrwin_continuous_isg
Ds<-c(4,6,8,10,14); bs<-c(0.06,0.10,0.14,0.18,0.25); seeds<-1:3; N<-12000; MS<-100
DA<-matrix(NA,length(seeds),length(Ds)); CA<-matrix(NA,length(seeds),length(bs))
for(si in seq_along(seeds)){ s<-seeds[si]
  cfg<-mrwin_config(n_outcome=N,m_snps=MS,alpha_x=c(-0.4,-0.4,-0.4),gamma_direct=c(0,0,0),seed=s)
  dat<-mrwin_simulate(cfg,seed=s); kern<-mrwin_kernel(dat$time,dat$status); sc<-as.numeric(dat$G%*%dat$true_betas)
  for(i in seq_along(Ds)){est<-mrwin_estimate(dat$time,dat$status,dat$G,dat$X,dat$true_betas,n_strata=Ds[i],kernel=kern)
    DA[si,i]<-tryCatch(mrwin_analytic_inference(est,dat$X)$delta_gls,error=function(e)NA)}
  for(i in seq_along(bs)) CA[si,i]<-tryCatch(CISG(kern,sc,dat$X,bandwidth=bs[i],n_centers=10,kernel_fn="epanechnikov")$delta_gls,error=function(e)NA)
}
cat(sprintf("LARGE-N (N=%d, %d cohorts) delta_gls -- plateau check:\n",N,length(seeds)))
cat("decile vs D:    "); for(i in seq_along(Ds)) cat(sprintf("D%d=%.3f ",Ds[i],mean(DA[,i],na.rm=TRUE))); cat("\n")
cat("continuous vs b:"); for(i in seq_along(bs)) cat(sprintf("b%.2f=%.3f ",bs[i],mean(CA[,i],na.rm=TRUE))); cat("\n")
cat(sprintf("range decile=%.3f  range continuous=%.3f  (per-cohort sd across cohorts: dec %.3f con %.3f)\n",
   diff(range(colMeans(DA,na.rm=TRUE))), diff(range(colMeans(CA,na.rm=TRUE))),
   mean(apply(DA,2,sd,na.rm=TRUE)), mean(apply(CA,2,sd,na.rm=TRUE))))
cat("DONE\n")
