#!/usr/bin/env Rscript
# RIGOROUS re-check of within-family / sibling win-ratio MR (Paper 04 candidate).
# Compares FOUR estimators across THREE confounding scenarios, against the
# interventional oracle. Crucially it gives population MR its standard defence --
# PRINCIPAL-COMPONENT adjustment -- so the comparison is NOT a strawman.
#
# Estimators (all use ONLY observed data: Z = offspring PRS, X, win-score from
# time/status, and the genotype matrix for PCs; NEVER the latent stratum, the
# confounder u, the parental PRS, or the oracle -- no information leakage):
#   naive-pop  : Cov(Z, win-score)/Cov(Z, X)                    (no stratification control)
#   PC-adj-pop : same on residuals after regressing out 10 genotype PCs (standard MR control)
#   within-fam : sum_sib h(O_i,O_j) dZ / sum_sib dX dZ          (sibling-pair restriction)
#
# Scenarios:
#   none  : no confounding (positive control -- all four must agree)
#   strat : population stratification (stratum shifts PRS and worsens outcome)
#   dyn   : dynastic / indirect genetic effect (parental PRS -> offspring outcome)
#
# RESULT (R=12, F=4000), bias vs oracle in parentheses:
#   none  : naive +0.01 | PC +0.01 | within -0.01   (all correct)
#   strat : naive -0.89 | PC -0.14 | within +0.02   (PC fixes MOST of stratification)
#   dyn   : naive -1.03 | PC -1.06 | within -0.01   (PC FAILS; only within-family works)
#
# HONEST CONCLUSION: within-family's unique value is NOT stratification (PC adjustment
# already handles that) but DYNASTIC / indirect-genetic (demonstrated) and, by the same
# Mendelian logic, ASSORTATIVE MATING (NOT yet tested here) -- exactly the confounders
# PC adjustment cannot remove. Still UNTESTED: assortative mating; the family-clustered
# AR inference (coverage/type-I); >2 sibs / trios; the two-sample within-sibship form.
suppressMessages(library(mrwin))
ncores<-max(1L,min(4L,parallel::detectCores())); cfg<-mrwin_config()
subj<-getFromNamespace("mrwin_subject_win_loss_cpp","mrwin")
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
gen_fam<-function(Fn,seed,mode,gstrat=2.5,fgap=0.18,gdyn=1.2){set.seed(seed);M<-30
 base<-runif(M,0.15,0.35); betas<-rnorm(M,0,0.15); strat<-rbinom(Fn,1,0.5)
 fg<-if(mode=="strat") fgap else 0
 p<-matrix(base,Fn,M,byrow=TRUE)+outer(strat,fg*sign(betas)); p<-pmin(pmax(p,0.02),0.95)
 Gf<-matrix(rbinom(Fn*M,2,p),Fn,M); Gm<-matrix(rbinom(Fn*M,2,p),Fn,M)
 parPRS<-as.numeric(((Gf+Gm)/2)%*%betas)
 sibg<-function() t(apply(Gf,1,transmit))+t(apply(Gm,1,transmit))
 mk<-function(Gs){Z<-as.numeric(Gs%*%betas)
   u<-rnorm(Fn); X<-cfg$alpha_s*Z+cfg$alpha_u*u+rnorm(Fn)
   extra<-(if(mode=="strat") gstrat*strat else 0)+(if(mode=="dyn") gdyn*scale(parPRS)[,1] else 0)
   w<-rgamma(Fn,1/cfg$theta_f,scale=cfg$theta_f); cn<-pmin(rexp(Fn,cfg$censoring_rate),cfg$max_follow_up)
   uf<-matrix(runif(Fn*3),Fn,3); oc<-outcome(X,u,w,extra,cn,uf)
   list(Z=Z,X=X,u=u,w=w,extra=extra,cn=cn,uf=uf,G=Gs,time=oc$time,status=oc$status)}
 list(s1=mk(sibg()),s2=mk(sibg()),strat=strat)}
winscore<-function(d){ti<-rbind(d$s1$time,d$s2$time);si<-rbind(d$s1$status,d$s2$status)
 cc<-subj(ti,si,ti,si,rep(1,nrow(ti))); (cc[,1]-cc[,2])/(nrow(ti)-1)}
naive_pop<-function(d,wsc){Z<-c(d$s1$Z,d$s2$Z);X<-c(d$s1$X,d$s2$X); cov(Z,wsc)/cov(Z,X)}
pc_pop<-function(d,wsc,k=10){G<-rbind(d$s1$G,d$s2$G); Z<-c(d$s1$Z,d$s2$Z);X<-c(d$s1$X,d$s2$X)
 pc<-prcomp(G,center=TRUE,scale.=FALSE)$x[,1:k]
 rz<-resid(lm(Z~pc)); rw<-resid(lm(wsc~pc)); rx<-resid(lm(X~pc)); cov(rz,rw)/cov(rz,rx)}
wf<-function(d){h<-hpair(d$s1$time,d$s1$status,d$s2$time,d$s2$status)
 dZ<-d$s1$Z-d$s2$Z; dX<-d$s1$X-d$s2$X; sum(h*dZ)/sum(dX*dZ)}
oracle<-function(d,eps=0.25){A<-d$s1;np<-length(A$X)
 Bp<-outcome(A$X+eps,A$u,A$w,A$extra,A$cn,A$uf); Bm<-outcome(A$X-eps,A$u,A$w,A$extra,A$cn,A$uf)
 v<-mrwin_fast_pair_win_loss(Bp$time,Bp$status,Bm$time,Bm$status)
 (as.numeric(v[["wins"]])-as.numeric(v[["losses"]]))/(as.numeric(np)^2)/(2*eps)}
cat(sprintf("%-8s %-8s %-14s %-14s %-14s\n","scenario","oracle","naive-pop","PC-adj-pop","within-fam"))
for(mode in c("none","strat","dyn")){
 res<-do.call(rbind,parallel::mclapply(1:12,function(r){d<-gen_fam(4000L,1000L+r,mode)
   ws<-winscore(d); c(orc=oracle(d),nv=naive_pop(d,ws),pc=pc_pop(d,ws),wf=wf(d))},mc.cores=ncores))
 m<-colMeans(res)
 cat(sprintf("%-8s %-8.3f %-14s %-14s %-14s\n",mode,m["orc"],
   sprintf("%.3f(%+.2f)",m["nv"],m["nv"]-m["orc"]),
   sprintf("%.3f(%+.2f)",m["pc"],m["pc"]-m["orc"]),
   sprintf("%.3f(%+.2f)",m["wf"],m["wf"]-m["orc"])))}
