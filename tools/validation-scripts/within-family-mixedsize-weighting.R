#!/usr/bin/env Rscript
# MIXED FAMILY SIZES -- is pair-pooling efficient? (the genuine remaining efficiency
# question after exchangeability killed WITHIN-family optimal weighting). Families have
# different sizes (1/3 each s=2,3,4). Each family f gives A_f=sum_{i<j in f} h_ij*dZ_ij,
# B_f=sum_{i<j in f} dX_ij*dZ_ij; m_f(b)=A_f-b*B_f is i.i.d. across families.
#
# Pair-pooling beta_hat = sum_f A_f / sum_f B_f gives each family weight w_f=1, so its
# influence scales with C(s,2) -> it "over-weights" big families. The GMM-efficient
# weight for this stratified ratio is w_s prop b_s/v_s, b_s=E[B_f|s], v_s=Var(m_f|s).
# Three estimators compared by EMPIRICAL sd over reps vs the interventional oracle:
#   PAIR  w=1            (pool raw pairs; "over-weights" big families)
#   FAMEQ w=1/C(s,2)     (each family equal; the intuitive "fix")
#   EFF   w=b_s/v_s      (two-step optimal stratified weighting)
#
# RESULT (F=3000, R=120, dynastic):
#   PAIR  bias +0.0136  sd 0.0422
#   FAMEQ bias +0.0166  sd 0.0544
#   EFF   bias +0.0119  sd 0.0423
#   rel-eff sd(PAIR)/sd(EFF)=0.997 ; sd(FAMEQ)/sd(EFF)=1.285 ; AR coverage 0.975
#
# CONCLUSION (positive + slightly counterintuitive): PAIR-POOLING IS ALREADY EFFICIENT
# even with mixed sizes -- optimal weighting buys nothing (rel-eff 1.00). The closed form
# explains it: rel-eff = [sum C(s,2)g_s * sum C(s,2)/g_s]/(sum C(s,2))^2 with
# g_s = nu + 2(s-2)c  (nu = per-pair var, c = covariance of two within-family pairs that
# share a sib). This is 1 iff c=0; the data say rel-eff~1, i.e. c << nu empirically (the
# per-pair win moment's idiosyncratic variance dominates the one-shared-sib covariance),
# so pooling raw pairs is ~optimal. The ONLY mistake is FAMEQ (normalizing per family):
# predicted sd penalty sqrt((1/6)/(1/10))=1.29, matching the observed 1.285. So the simple
# estimator (pool all within-family pairs, family-cluster the AR variance) is BOTH simplest
# AND efficient under mixed family sizes -- do NOT down-weight large families. Together with
# within-family-efficient-weighting.R (exchangeability => equal pairs optimal WITHIN a
# family), the efficiency question is fully resolved: equal-pair-weight clustered AR is the
# efficient estimator. Public package + simulator only; estimators use observed data.
suppressMessages(library(mrwin)); ncores<-max(1L,min(4L,parallel::detectCores())); cfg<-mrwin_config()
hpair<-function(t1,s1,t2,s2){F<-nrow(t1);r<-integer(F)
 for(k in 1:3){d<-r!=0; iw<-(s2[,k]==1&t1[,k]>t2[,k]); il<-(s1[,k]==1&t2[,k]>t1[,k])
  r[!d&iw]<-1L; r[!d&il&!iw]<- -1L}; r}
outc<-function(X,u,w,ex,cn,uf){F<-length(X);et<-matrix(NA_real_,F,3);lw<-log(w)
 for(j in 1:3){lp<-cfg$alpha_x[j]*X+cfg$nu_u[j]*u+ex+lw; et[,j]<-(-log(uf[,j])/(cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)}
 td<-et[,1];th<-et[,2];tr<-et[,3]
 list(time=cbind(pmin(td,cn),pmin(th,cn,td),pmin(tr,cn,td)),status=cbind(as.integer(td<=cn),as.integer(th<=cn&th<=td),as.integer(tr<=cn&tr<=td)))}
tr1<-function(G) ifelse(G==2,1L,ifelse(G==0,0L,rbinom(length(G),1,0.5)))
genS<-function(Fn,seed,smax=4L,gdyn=1.2){set.seed(seed);M<-30;base<-runif(M,.15,.35);betas<-rnorm(M,0,.15)
 cc<-col(matrix(0,Fn,M)); Gf<-matrix(rbinom(Fn*M,2,base[cc]),Fn,M); Gm<-matrix(rbinom(Fn*M,2,base[cc]),Fn,M)
 pp<-scale(as.numeric(((Gf+Gm)/2)%*%betas))[,1]; sib<-vector("list",smax); o1<-NULL
 for(k in 1:smax){Gs<-matrix(tr1(Gf),Fn,M)+matrix(tr1(Gm),Fn,M); Z<-as.numeric(Gs%*%betas)
  u<-rnorm(Fn); X<-cfg$alpha_s*Z+cfg$alpha_u*u+rnorm(Fn); ex<-gdyn*pp
  w<-rgamma(Fn,1/cfg$theta_f,scale=cfg$theta_f); cn<-pmin(rexp(Fn,cfg$censoring_rate),cfg$max_follow_up); uf<-matrix(runif(Fn*3),Fn,3)
  oc<-outc(X,u,w,ex,cn,uf); sib[[k]]<-list(Z=Z,X=X,time=oc$time,status=oc$status); if(k==1) o1<-list(u=u,w=w,ex=ex,cn=cn,uf=uf,X=X)}
 list(sib=sib,o1=o1)}
oracle<-function(d,eps=0.25){A<-d$o1;np<-length(A$X)
 Bp<-outc(A$X+eps,A$u,A$w,A$ex,A$cn,A$uf); Bm<-outc(A$X-eps,A$u,A$w,A$ex,A$cn,A$uf)
 v<-mrwin_fast_pair_win_loss(Bp$time,Bp$status,Bm$time,Bm$status)
 (as.numeric(v[["wins"]])-as.numeric(v[["losses"]]))/(as.numeric(np)^2)/(2*eps)}
pidx<-rbind(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))            # all 6 pairs of 4 sibs
actcols<-list(`2`=1L,`3`=c(1L,2L,4L),`4`=1:6)                     # pairs present at each size
allpairs<-function(d){Fn<-length(d$sib[[1]]$Z);am<-matrix(0,Fn,6);bm<-matrix(0,Fn,6)
 for(p in 1:6){i<-pidx[p,1];j<-pidx[p,2];si<-d$sib[[i]];sj<-d$sib[[j]]
   h<-hpair(si$time,si$status,sj$time,sj$status);dZ<-si$Z-sj$Z;dX<-si$X-sj$X;am[,p]<-h*dZ;bm[,p]<-dX*dZ}
 list(am=am,bm=bm)}
famAB<-function(pm,sizes){Fn<-nrow(pm$am);A<-numeric(Fn);B<-numeric(Fn)
 for(s in c(2,3,4)){idx<-which(sizes==s);cols<-actcols[[as.character(s)]]
   A[idx]<-rowSums(pm$am[idx,cols,drop=FALSE]);B[idx]<-rowSums(pm$bm[idx,cols,drop=FALSE])}
 list(A=A,B=B)}
run<-function(d,sizes,bstar){pm<-allpairs(d);ab<-famAB(pm,sizes);A<-ab$A;B<-ab$B
 b_pair<-sum(A)/sum(B)
 we<-1/choose(sizes,2); b_feq<-sum(we*A)/sum(we*B)
 m<-A-b_pair*B; w<-numeric(length(A))
 for(s in c(2,3,4)){idx<-which(sizes==s);w[idx]<-mean(B[idx])/mean(m[idx]^2)}
 b_eff<-sum(w*A)/sum(w*B)
 covp<-(sum(A-bstar*B)^2/sum((A-bstar*B)^2))<=3.841
 c(o=bstar,b_pair=b_pair,b_feq=b_feq,b_eff=b_eff,covp=as.numeric(covp))}
Fn<-as.integer(Sys.getenv("MIX_F","3000")); R<-as.integer(Sys.getenv("MIX_R","120"))
cat(sprintf("MIXED family sizes (1/3 each s=2,3,4; F=%d, R=%d, dynastic)\n%-13s %-9s %-9s\n",Fn,R,"estimator","bias","sd"))
res<-do.call(rbind,parallel::mclapply(1:R,function(r){d<-genS(Fn,9000L+r);o<-oracle(d)
   sizes<-rep(c(2L,3L,4L),length.out=Fn);run(d,sizes,o)},mc.cores=ncores))
m<-colMeans(res)
for(nm in c("b_pair","b_feq","b_eff")){lab<-c(b_pair="PAIR(w=1)",b_feq="FAMEQ(1/Cs2)",b_eff="EFF(b_s/v_s)")[nm]
 cat(sprintf("%-13s %+9.4f %9.4f\n",lab,m[nm]-m["o"],sd(res[,nm])))}
cat(sprintf("clustered-AR coverage (pair-pooled) = %.3f\n",m["covp"]))
cat(sprintf("rel-eff sd(PAIR)/sd(EFF)=%.3f ; sd(FAMEQ)/sd(EFF)=%.3f\n",
  sd(res[,"b_pair"])/sd(res[,"b_eff"]), sd(res[,"b_feq"])/sd(res[,"b_eff"])))
cat("DONEMIX\n")
