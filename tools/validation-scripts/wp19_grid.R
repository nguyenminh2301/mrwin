suppressMessages(library(mrwin))
V <- mrwin:::.mrwin_validate_calibration
pr <- function(r) cat(sprintf("[%-8s %-9s sb=%-4g %-12s] rate=%.3f (SE %.3f, n=%d, weak=%.2f)\n",
   r$scenario, r$inference, r$sigma_beta, r$adjustment, r$rate, r$se, r$valid, r$weak_frac))
N<-3000; MS<-100; D<-5
# ---- TYPE-I (null), Fieller ----
pr(V(M=200,N=N,m_snps=MS,n_strata=D, scenario="null", inference="analytic",  sigma_beta=0,    adjustment="none",        seed0=1000))
pr(V(M=150,N=N,m_snps=MS,n_strata=D, scenario="null", inference="analytic",  sigma_beta=0.01, adjustment="none",  bootstrap=150, seed0=2000))
pr(V(M=150,N=N,m_snps=MS,n_strata=D, scenario="null", inference="bootstrap", sigma_beta=0,    adjustment="none",  bootstrap=150, seed0=3000))
pr(V(M=150,N=N,m_snps=MS,n_strata=D, scenario="null", inference="bootstrap", sigma_beta=0.01, adjustment="none",  bootstrap=150, seed0=4000))
pr(V(M=150,N=N,m_snps=MS,n_strata=D, scenario="null", inference="analytic",  sigma_beta=0,    adjustment="ordinal_iptw", seed0=5000))
# ---- COVERAGE (valid_iv), Fieller ----
tr <- vapply(1:6, function(i){
  cfg<-mrwin_config(n_outcome=12000,m_snps=MS,alpha_x=c(-0.4,-0.4,-0.4),gamma_direct=c(0,0,0),seed=70000+i)
  dat<-mrwin_simulate(cfg,seed=70000+i); ep<-mrwin_endpoint(dat$time,dat$status,colnames(dat$time))
  gw<-mrwin_gwas(dat$true_betas,rep(0,length(dat$true_betas)))
  mrwin(endpoint=ep,genotype=dat$G,exposure=dat$X,gwas=gw,controls=mrwin_controls(n_strata=D,bootstrap=50,seed=70000+i,inference="analytic",run_sdpd=FALSE))$point$delta_gls
}, numeric(1))
truth<-mean(tr); cat(sprintf("valid_iv truth=%.4f (rep sd %.4f)\n", truth, sd(tr)))
pr(V(M=150,N=N,m_snps=MS,n_strata=D, scenario="valid_iv", inference="bootstrap", sigma_beta=0.01, adjustment="none", bootstrap=150, truth=truth, seed0=6000))
cat("DONE\n")
