suppressMessages(library(mrwin))
CISG <- mrwin:::mrwin_continuous_isg
Ds <- c(3,4,5,6,8,10,12,15); bs <- c(0.05,0.07,0.09,0.11,0.13,0.15,0.18,0.22)
seeds <- 1:8; N<-6000; MS<-90
dec_cv <- con_cv <- numeric(length(seeds)); dec_all <- con_all <- list()
for (si in seq_along(seeds)) {
  s<-seeds[si]
  cfg<-mrwin_config(n_outcome=N,m_snps=MS,alpha_x=c(-0.4,-0.4,-0.4),gamma_direct=c(0,0,0),seed=s)
  dat<-mrwin_simulate(cfg,seed=s); kern<-mrwin_kernel(dat$time,dat$status)
  sc<-as.numeric(dat$G%*%dat$true_betas)
  dvals<-sapply(Ds,function(D){est<-mrwin_estimate(dat$time,dat$status,dat$G,dat$X,dat$true_betas,n_strata=D,kernel=kern)
     tryCatch(mrwin_analytic_inference(est,dat$X)$delta_gls,error=function(e)NA)})
  cvals<-sapply(bs,function(b){tryCatch(CISG(kern,sc,dat$X,bandwidth=b,n_centers=12,kernel_fn="epanechnikov")$delta_gls,error=function(e)NA)})
  dec_all[[si]]<-dvals; con_all[[si]]<-cvals
  dec_cv[si]<-sd(dvals,na.rm=TRUE)/abs(mean(dvals,na.rm=TRUE))
  con_cv[si]<-sd(cvals,na.rm=TRUE)/abs(mean(cvals,na.rm=TRUE))
}
DM<-rowMeans(sapply(dec_all,identity),na.rm=TRUE); CM<-rowMeans(sapply(con_all,identity),na.rm=TRUE)
cat("decile delta_gls vs D (mean over seeds):\n"); for(i in seq_along(Ds)) cat(sprintf("  D=%2d: %.4f\n",Ds[i],DM[i]))
cat("continuous delta_gls vs b (mean over seeds):\n"); for(i in seq_along(bs)) cat(sprintf("  b=%.2f: %.4f\n",bs[i],CM[i]))
cat(sprintf("\nTUNING SENSITIVITY (CV across tuning, lower=less sensitive, mean over %d cohorts):\n", length(seeds)))
cat(sprintf("  decile (vary D):     CV = %.3f\n", mean(dec_cv,na.rm=TRUE)))
cat(sprintf("  continuous (vary b): CV = %.3f\n", mean(con_cv,na.rm=TRUE)))
cat("DONE\n")
