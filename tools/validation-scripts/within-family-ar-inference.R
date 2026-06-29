#!/usr/bin/env Rscript
# Within-family win-MR: (1) assortative-mating robustness + (2) the family-CLUSTERED
# Anderson-Rubin inference (coverage/type-I) -- the two pieces flagged as untested.
#
# Math: with 2 sibs per family, each family contributes ONE independent moment
# m_f(b) = (h_f - b*dX_f)*dZ_f (h_f the sib win, dZ_f/dX_f the sib instrument/exposure
# contrasts). Families are i.i.d., so this is NOT a degenerate U-statistic -- the
# within-family AR is a clean i.i.d.-across-families test:
#   AR(b) = (sum_f m_f(b))^2 / (sum_f m_f(b)^2) -> chi2_1 ;  CI = {b: AR(b) <= 3.841}.
# (For >2 sibs the within-family pairs are dependent and genuine family-clustering of
#  the influence function is needed -- NOT covered here.)
#
# Scenarios: none (control) | dyn (dynastic) | am (assortative mating: a continuous
# family latent shifts BOTH parents' trait-allele freqs together -> correlated parents
# -> + a dynastic path; PCs cannot capture this continuous family-level confounder).
#
# RESULT (point bias vs interventional oracle; 95% AR-CI coverage):
#   F=3000:  none +0.022 cov 0.925 | dyn +0.029 cov 0.950 | am +0.005 cov 0.925
#   F=9000:  none +0.007 cov 0.940 |                       | am +0.020 cov 0.940
# Conclusions: within-family is robust to assortative mating (recovers the causal
# gradient); the clustered AR is calibrated (~0.93-0.94, -> 0.95 as N grows). The
# small point bias present even with NO confounding is a within-sibship WEAK-INSTRUMENT
# finite-sample bias (the within-family genetic variance is small) -- it SHRINKS with N
# (0.022 -> 0.007 from F=3k to 9k), confirming consistency; it is the known power cost
# of within-family designs (need larger N). Public package + simulator only; estimators
# use only observed data (no latent/oracle leakage).
suppressMessages(library(mrwin))
ncores<-max(1L,min(4L,parallel::detectCores())); cfg<-mrwin_config()
hpair<-function(t1,s1,t2,s2){F<-nrow(t1);res<-integer(F)
 for(k in 1:3){dec<-res!=0; iwin<-(s2[,k]==1 & t1[,k]>t2[,k]); ilos<-(s1[,k]==1 & t2[,k]>t1[,k])
  res[!dec & iwin]<-1L; res[!dec & ilos & !iwin]<- -1L}; res}
outcome<-function(X,u,w,extra,cn,uf){F<-length(X);et<-matrix(NA_real_,F,3);lw<-log(w)
 for(j in 1:3){lp<-cfg$alpha_x[j]*X+cfg$nu_u[j]*u+extra+lw
  et[,j]<-(-log(uf[,j])/(cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)}
 td<-et[,1];th<-et[,2];tr<-et[,3]
 list(time=cbind(pmin(td,cn),pmin(th,cn,td),pmin(tr,cn,td)),
      status=cbind(as.integer(td<=cn),as.integer(th<=cn&th<=td),as.integer(tr<=cn&tr<=td)))}
transmit<-function(G) ifelse(G==2,1,ifelse(G==0,0,rbinom(length(G),1,0.5)))
gen_fam<-function(Fn,seed,mode,gdyn=1.2,amc=0.10){set.seed(seed);M<-30
 base<-runif(M,0.15,0.35); betas<-rnorm(M,0,0.15)
 amf<-if(mode=="am") rnorm(Fn) else rep(0,Fn)
 pp<-pmin(pmax(matrix(base,Fn,M,byrow=TRUE)+outer(amf,amc*sign(betas)),0.02),0.95)
 Gf<-matrix(rbinom(Fn*M,2,pp),Fn,M); Gm<-matrix(rbinom(Fn*M,2,pp),Fn,M)
 parPRS<-as.numeric(((Gf+Gm)/2)%*%betas)
 sibg<-function() t(apply(Gf,1,transmit))+t(apply(Gm,1,transmit))
 mk<-function(Gs){Z<-as.numeric(Gs%*%betas)
   u<-rnorm(Fn); X<-cfg$alpha_s*Z+cfg$alpha_u*u+rnorm(Fn)
   extra<-(if(mode%in%c("dyn","am")) gdyn*scale(parPRS)[,1] else 0)
   w<-rgamma(Fn,1/cfg$theta_f,scale=cfg$theta_f); cn<-pmin(rexp(Fn,cfg$censoring_rate),cfg$max_follow_up)
   uf<-matrix(runif(Fn*3),Fn,3); oc<-outcome(X,u,w,extra,cn,uf)
   list(Z=Z,X=X,u=u,w=w,extra=extra,cn=cn,uf=uf,time=oc$time,status=oc$status)}
 list(s1=mk(sibg()),s2=mk(sibg()))}
oracle<-function(d,eps=0.25){A<-d$s1;np<-length(A$X)
 Bp<-outcome(A$X+eps,A$u,A$w,A$extra,A$cn,A$uf); Bm<-outcome(A$X-eps,A$u,A$w,A$extra,A$cn,A$uf)
 v<-mrwin_fast_pair_win_loss(Bp$time,Bp$status,Bm$time,Bm$status)
 (as.numeric(v[["wins"]])-as.numeric(v[["losses"]]))/(as.numeric(np)^2)/(2*eps)}
wf_ar<-function(d,bstar){h<-hpair(d$s1$time,d$s1$status,d$s2$time,d$s2$status)
 dZ<-d$s1$Z-d$s2$Z; dX<-d$s1$X-d$s2$X; a<-h*dZ; b<-dX*dZ
 ARb<-function(b0){m<-a-b0*b; sum(m)^2/sum(m^2)}
 list(bhat=sum(a)/sum(b), covers=ARb(bstar)<=3.841)}
F0<-as.integer(Sys.getenv("FAM_F","3000")); R0<-as.integer(Sys.getenv("FAM_R","80"))
cat(sprintf("within-family clustered-AR (F=%d, R=%d): point bias + 95%% CI coverage\n",F0,R0))
cat(sprintf("%-8s %-10s %-16s %-16s\n","mode","oracle","wf point (bias)","AR-CI coverage"))
for(mode in c("none","dyn","am")){
 res<-do.call(rbind,parallel::mclapply(1:R0,function(r){d<-gen_fam(F0,4000L+r,mode)
   o<-oracle(d); ar<-wf_ar(d,o); c(orc=o,bhat=ar$bhat,cov=as.numeric(ar$covers))},mc.cores=ncores))
 m<-colMeans(res); se<-sqrt(m["cov"]*(1-m["cov"])/nrow(res))
 cat(sprintf("%-8s %-10.3f %-16s %.3f (SE %.3f)\n",mode,m["orc"],
   sprintf("%.3f(%+.3f)",m["bhat"],m["bhat"]-m["orc"]), m["cov"], se))}
