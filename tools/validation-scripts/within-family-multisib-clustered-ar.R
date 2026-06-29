#!/usr/bin/env Rscript
# >2 sibs family-clustered AR for within-family win-MR -- proves the clustering is
# NECESSARY and SUFFICIENT. With s sibs the C(s,2) within-family pairs are dependent;
# aggregating them to the family level m_f=sum_{i<j in f}(h_ij-b*dX_ij)dZ_ij gives an
# i.i.d.-across-families moment, so AR=(sum_f m_f)^2/sum_f m_f^2 -> chi2_1.
# RESULT (F=3000,R=80, dynastic): point bias stable +0.012..+0.014 (within-sibship
# weak-IV finite-sample bias, not growing with s); 95% CI coverage:
#   s=2: clustered 0.950 = naive 0.950   (1 pair/family -- identical)
#   s=3: clustered 0.925 , naive 0.887
#   s=4: clustered 0.950 , naive 0.850   <- naive DEGRADES, clustered HOLDS
# => within-family pair dependence is real (naive under-covers as s grows) and the
# family-clustered AR fixes it. Public package + simulator; observed data only.
# Prove the family-clustering is NECESSARY, not just sufficient: compare
#   clustered : AR=(sum_f m_f)^2 / sum_f m_f^2     (m_f = within-family pair sum)
#   naive     : AR=(sum_p m_p)^2 / sum_p m_p^2     (treats C(s,2) pairs as independent)
# same numerator; denominators differ for s>2. If naive under-covers as s grows while
# clustered holds ~0.95, the within-family pair dependence is real and clustering fixes it.
suppressMessages(library(mrwin)); ncores<-max(1L,min(4L,parallel::detectCores())); cfg<-mrwin_config()
hpair<-function(t1,s1,t2,s2){F<-nrow(t1);r<-integer(F)
 for(k in 1:3){d<-r!=0; iw<-(s2[,k]==1&t1[,k]>t2[,k]); il<-(s1[,k]==1&t2[,k]>t1[,k])
  r[!d&iw]<-1L; r[!d&il&!iw]<- -1L}; r}
outc<-function(X,u,w,ex,cn,uf){F<-length(X);et<-matrix(NA_real_,F,3);lw<-log(w)
 for(j in 1:3){lp<-cfg$alpha_x[j]*X+cfg$nu_u[j]*u+ex+lw; et[,j]<-(-log(uf[,j])/(cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)}
 td<-et[,1];th<-et[,2];tr<-et[,3]
 list(time=cbind(pmin(td,cn),pmin(th,cn,td),pmin(tr,cn,td)),
      status=cbind(as.integer(td<=cn),as.integer(th<=cn&th<=td),as.integer(tr<=cn&tr<=td)))}
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
estim<-function(d,bs){s<-length(d$sib);Fn<-length(d$sib[[1]]$Z);np<-s*(s-1)/2
 am<-matrix(0,Fn,np);bm<-matrix(0,Fn,np);c0<-0
 for(i in 1:(s-1)) for(j in (i+1):s){c0<-c0+1;si<-d$sib[[i]];sj<-d$sib[[j]]
   h<-hpair(si$time,si$status,sj$time,sj$status); dZ<-si$Z-sj$Z; dX<-si$X-sj$X
   am[,c0]<-h*dZ; bm[,c0]<-dX*dZ}
 Af<-rowSums(am);Bf<-rowSums(bm); N<-sum(Af)-bs*sum(Bf)
 clus<-N^2/sum((Af-bs*Bf)^2); naive<-N^2/sum((am-bs*bm)^2)
 c(bhat=sum(Af)/sum(Bf), cl=as.numeric(clus<=3.841), nv=as.numeric(naive<=3.841))}
cat(sprintf("%-4s %-8s %-13s %-13s %-13s\n","s","oracle","bias","cov-CLUSTERED","cov-naive"))
for(s in 2:4){res<-do.call(rbind,parallel::mclapply(1:80,function(r){d<-genS(3000L,5000L+r,s)
   o<-oracle(d); e<-estim(d,o); c(o=o,e)},mc.cores=ncores))
 m<-colMeans(res);se<-function(p)sqrt(p*(1-p)/nrow(res))
 cat(sprintf("%-4d %-8.3f %-13s %.3f(SE%.3f)  %.3f(SE%.3f)\n",s,m["o"],sprintf("%+.3f",m["bhat"]-m["o"]),m["cl"],se(m["cl"]),m["nv"],se(m["nv"])))}
cat("DONES2\n")
