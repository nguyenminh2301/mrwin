#!/usr/bin/env Rscript
# Efficient (optimally-weighted) within-family win-MR -- NEGATIVE FINDING.
# Each family gives a P=C(s,2) pair-moment vector; the GMM-efficient weight is
# Sigma_g^-1. RESULT (F=3000,R=90,dynastic): the efficient (CUE) point estimator
# has the SAME sd as the equal-weight one (s=2: .076/.076; s=3: .057/.057; s=4:
# .043/.044 -- efficient slightly WORSE); coverage ~0.93-0.96 both.
# WHY (theory): sib pairs within a family are EXCHANGEABLE -> Sigma_g is compound-
# symmetric and every pair loads equally on beta -> the optimal weights ARE equal.
# Estimating the PxP weight only adds finite-sample noise. So equal-weight clustered
# AR is already (essentially) efficient for full-sibs; optimal within-family weighting
# is a dead end. The genuine efficiency question is BETWEEN families of different
# SIZES (equal-pair-weighting over-weights large families) = the mixed-size problem.
# Efficient (optimally-weighted) within-family win-MR for s sibs. Each family gives a
# P=C(s,2) moment vector g_f(b)=a_f-b*b_f. The within-family pairs are correlated, so
# equal-weighting (summing pairs) is sub-optimal; the efficient GMM weight is Sigma_g^-1.
# Efficient AR (continuous-updating): AR(b)=(colSums G)' (G'G)^-1 (colSums G) -> chi2_P,
# G=am-b*bm (F x P). Compare empirical sd of the EQUAL-weight vs EFFICIENT point
# estimator (the efficiency gain) and the efficient-AR coverage. Observed data only.
suppressMessages(library(mrwin)); ncores<-max(1L,min(4L,parallel::detectCores())); cfg<-mrwin_config()
hpair<-function(t1,s1,t2,s2){F<-nrow(t1);r<-integer(F)
 for(k in 1:3){d<-r!=0; iw<-(s2[,k]==1&t1[,k]>t2[,k]); il<-(s1[,k]==1&t2[,k]>t1[,k])
  r[!d&iw]<-1L; r[!d&il&!iw]<- -1L}; r}
outc<-function(X,u,w,ex,cn,uf){F<-length(X);et<-matrix(NA_real_,F,3);lw<-log(w)
 for(j in 1:3){lp<-cfg$alpha_x[j]*X+cfg$nu_u[j]*u+ex+lw; et[,j]<-(-log(uf[,j])/(cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)}
 td<-et[,1];th<-et[,2];tr<-et[,3]
 list(time=cbind(pmin(td,cn),pmin(th,cn,td),pmin(tr,cn,td)),status=cbind(as.integer(td<=cn),as.integer(th<=cn&th<=td),as.integer(tr<=cn&tr<=td)))}
tr1<-function(G) ifelse(G==2,1L,ifelse(G==0,0L,rbinom(length(G),1,0.5)))
genS<-function(Fn,seed,s,gdyn=1.2){set.seed(seed);M<-30;base<-runif(M,.15,.35);betas<-rnorm(M,0,.15)
 cc<-col(matrix(0,Fn,M)); Gf<-matrix(rbinom(Fn*M,2,base[cc]),Fn,M); Gm<-matrix(rbinom(Fn*M,2,base[cc]),Fn,M)
 pp<-scale(as.numeric(((Gf+Gm)/2)%*%betas))[,1]; sib<-vector("list",s); o1<-NULL
 for(k in 1:s){Gs<-matrix(tr1(Gf),Fn,M)+matrix(tr1(Gm),Fn,M); Z<-as.numeric(Gs%*%betas)
  u<-rnorm(Fn); X<-cfg$alpha_s*Z+cfg$alpha_u*u+rnorm(Fn); ex<-gdyn*pp
  w<-rgamma(Fn,1/cfg$theta_f,scale=cfg$theta_f); cn<-pmin(rexp(Fn,cfg$censoring_rate),cfg$max_follow_up); uf<-matrix(runif(Fn*3),Fn,3)
  oc<-outc(X,u,w,ex,cn,uf); sib[[k]]<-list(Z=Z,X=X,time=oc$time,status=oc$status); if(k==1) o1<-list(u=u,w=w,ex=ex,cn=cn,uf=uf,X=X)}
 list(sib=sib,o1=o1)}
oracle<-function(d,eps=0.25){A<-d$o1;np<-length(A$X)
 Bp<-outc(A$X+eps,A$u,A$w,A$ex,A$cn,A$uf); Bm<-outc(A$X-eps,A$u,A$w,A$ex,A$cn,A$uf)
 v<-mrwin_fast_pair_win_loss(Bp$time,Bp$status,Bm$time,Bm$status)
 (as.numeric(v[["wins"]])-as.numeric(v[["losses"]]))/(as.numeric(np)^2)/(2*eps)}
pairmats<-function(d){s<-length(d$sib);Fn<-length(d$sib[[1]]$Z);np<-s*(s-1)/2;am<-matrix(0,Fn,np);bm<-matrix(0,Fn,np);c0<-0
 for(i in 1:(s-1)) for(j in (i+1):s){c0<-c0+1;si<-d$sib[[i]];sj<-d$sib[[j]]
   h<-hpair(si$time,si$status,sj$time,sj$status); dZ<-si$Z-sj$Z; dX<-si$X-sj$X; am[,c0]<-h*dZ; bm[,c0]<-dX*dZ}
 list(am=am,bm=bm)}
ARv<-function(am,bm,b0){G<-am-b0*bm; cs<-colSums(G); P<-ncol(G); M<-crossprod(G)+diag(1e-8*mean(diag(crossprod(G)))+1e-12,P)
 as.numeric(t(cs)%*%solve(M,cs))}
cat(sprintf("efficient (optimally-weighted) within-family win-MR; F=3000,R=90, dynastic\n%-4s %-8s %-13s %-13s %-9s %-9s\n","s","oracle","bias(eq)","bias(eff)","sd_eq","sd_eff/cov"))
for(s in 2:4){P<-s*(s-1)/2; cr<-qchisq(.95,P)
 res<-do.call(rbind,parallel::mclapply(1:90,function(r){d<-genS(3000L,7000L+r,s);o<-oracle(d);pm<-pairmats(d)
   beq<-sum(pm$am)/sum(pm$bm)
   beff<-optimize(function(b)ARv(pm$am,pm$bm,b),c(beq-0.4,beq+0.4))$minimum
   c(o=o,beq=beq,beff=beff,cov=as.numeric(ARv(pm$am,pm$bm,o)<=cr))},mc.cores=ncores))
 sdq<-sd(res[,"beq"]);sde<-sd(res[,"beff"]);m<-colMeans(res)
 cat(sprintf("%-4d %-8.3f %-13s %-13s %-9.3f sd_eff=%.3f cov=%.3f\n",s,m["o"],
   sprintf("%+.3f",m["beq"]-m["o"]),sprintf("%+.3f",m["beff"]-m["o"]),sdq,sde,m["cov"]))}
cat("DONEEFF\n")
