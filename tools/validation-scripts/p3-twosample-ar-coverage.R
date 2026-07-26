#!/usr/bin/env Rscript
# Two-sample mrwin_winmr_ar() coverage check (Paper 03). SHARED genetic
# architecture across the exposure and outcome cohorts (only individuals differ) --
# the correct two-sample setup. Result: AR point plim = 0.307 (= interventional
# oracle ~0.30, consistent), AR-CI coverage of the plim ~0.90-0.92 at N_out=8k-16k
# (approaching nominal; matches the one-sample AR calibration). A non-shared-betas
# version spuriously shows ~0.1 coverage -- a reminder to fix the genetic
# architecture across two-sample cohorts. Public package + simulator only.
suppressMessages(library(mrwin)); ncores<-min(4L,parallel::detectCores())
cfg<-mrwin_config(m_snps=40L)
# SHARED genetic architecture across cohorts (only individuals differ) -- the
# correct two-sample setup; the earlier check drew fresh betas per cohort (bug).
set.seed(7); mafs<-runif(40,cfg$maf_low,cfg$maf_high); betas<-rnorm(40,0,cfg$sigma_beta)
gen<-function(n,seed,as_){set.seed(seed)
 g<-scale(sapply(mafs,function(p)rbinom(n,2,p)));storage.mode(g)<-"double"
 s<-drop(g%*%betas);u<-rnorm(n);w<-rgamma(n,1/cfg$theta_f,scale=cfg$theta_f)
 cn<-pmin(rexp(n,cfg$censoring_rate),cfg$max_follow_up);uf<-matrix(runif(n*3),n,3)
 X<-as_*s+cfg$alpha_u*u+rnorm(n)
 et<-matrix(NA,n,3);lw<-log(w);for(j in 1:3){lp<-cfg$alpha_x[j]*X+cfg$nu_u[j]*u+lw
  et[,j]<-(-log(uf[,j])/(cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)}
 td<-et[,1];th<-et[,2];tr<-et[,3]
 list(time=cbind(pmin(td,cn),pmin(th,cn,td),pmin(tr,cn,td)),
   status=cbind(as.integer(td<=cn),as.integer(th<=cn&th<=td),as.integer(tr<=cn&tr<=td)),G=g,X=X)}
expg<-function(G,X){n<-length(X);apply(G,2,function(d){dc<-d-mean(d);b<-sum(dc*X)/sum(dc^2)
  res<-X-mean(X)-b*dc;c(b=b,se=sqrt(sum(res^2)/((n-2)*sum(dc^2))))})}
# plim of the AR point estimate (shared arch), large outcome N
tr<-unlist(parallel::mclapply(1:4,function(r){pe<-gen(60000L,100L+r,0.4);eg<-expg(pe$G,pe$X)
  po<-gen(30000L,200L+r,0.4);wg<-mrwin_win_gwas(po$time,po$status,po$G,max_subjects=600L,seed=r)
  mrwin_winmr_ar(eg["b",],wg$delta,wg$se,se_gx=eg["se",])$gamma},mc.cores=ncores))
plim<-mean(tr); cat(sprintf("AR point plim (N_out=30k, shared arch) = %.3f (sd %.3f)\n",plim,sd(tr)))
# coverage of the plim at moderate N (calibration of the AR CI itself)
for(Nout in c(8000L,16000L)){hit<-unlist(parallel::mclapply(1:50,function(r){
  pe<-gen(60000L,5e4+r,0.4);eg<-expg(pe$G,pe$X)
  po<-gen(Nout,6e4+r,0.4);wg<-mrwin_win_gwas(po$time,po$status,po$G,max_subjects=500L,seed=r)
  ar<-mrwin_winmr_ar(eg["b",],wg$delta,wg$se,se_gx=eg["se",])
  ok<-ar$ci[1]<=plim && plim<=ar$ci[2]; if(is.na(ok))FALSE else ok},mc.cores=ncores))
 cat(sprintf("N_out=%5d  AR-CI coverage of plim=%.3f : %.2f (n=%d)\n",Nout,plim,mean(hit),length(hit)))}
cat("DONE\n")
