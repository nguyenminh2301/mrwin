suppressMessages(library(mrwin)); CISG<-mrwin:::mrwin_continuous_isg
seeds<-1:8; N<-8000; MS<-100
dec<-cgls<-cwald<-numeric(length(seeds)); condn<-numeric(length(seeds))
for(si in seq_along(seeds)){ s<-seeds[si]
  cfg<-mrwin_config(n_outcome=N,m_snps=MS,alpha_x=c(-0.4,-0.4,-0.4),gamma_direct=c(0,0,0),seed=s)
  dat<-mrwin_simulate(cfg,seed=s); kern<-mrwin_kernel(dat$time,dat$status); sc<-as.numeric(dat$G%*%dat$true_betas)
  est<-mrwin_estimate(dat$time,dat$status,dat$G,dat$X,dat$true_betas,n_strata=5,kernel=kern)
  dec[si]<-mrwin_analytic_inference(est,dat$X)$delta_gls
  cc<-CISG(kern,sc,dat$X,bandwidth=0.10,n_centers=10,kernel_fn="epanechnikov")
  cgls[si]<-cc$delta_gls
  cwald[si]<-sum(cc$log_theta)/sum(cc$delta_x)         # simple Wald ratio (no GLS-of-ratios)
  condn[si]<-kappa(cc$sigma_isg)                        # ISG covariance condition number
}
f<-function(x) sprintf("mean=%.3f sd=%.3f", mean(x,na.rm=TRUE), sd(x,na.rm=TRUE))
cat(sprintf("N=%d, %d cohorts (matched scale: decile D=5 ~ continuous b=0.10):\n",N,length(seeds)))
cat("  decile (Fieller pool):       ", f(dec),"\n")
cat("  continuous (GLS-of-ratios):  ", f(cgls),"\n")
cat("  continuous (simple Wald):    ", f(cwald),"\n")
cat(sprintf("  continuous sigma_isg condition number: median %.1e (near-singular if huge)\n", median(condn)))
cat("DONE\n")
