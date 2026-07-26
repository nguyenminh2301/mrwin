#!/usr/bin/env Rscript
# P3 section 6: Andrews-Cheng sup-over-CI degeneracy-robust AR test, implemented and
# probed on the synthetic kernel phi=eps*(Ai+Aj)+Ai*Aj (zeta1=eps^2, rank-1
# degenerate part). Construction:
#   - the limit quantile q(c^2) = quantile of R(c)^2 is MONOTONE INCREASING in
#     c^2 = n*zeta1 (the Gaussian part has variance 4c^2), so the sup over a CI for
#     c^2 is attained at its UPPER endpoint c2_hi;
#   - build c2_hi by BOOTSTRAP (the analytic iid se of zeta1_hat underestimates badly
#     because the empirical projections g_hat_i are dependent through the degenerate
#     common mode); critical value = (1-a2) quantile of R(c2_hi, lam_hat)^2;
#   - reject (n*U_n)^2 > crit, with a1 + a2 = 0.05 (Bonferroni).
#
# FINDING (dev/p3-weak-iv-robust-winmr.md section 6): the sup-CI test is uniformly
# VALID (size <= 0.05 everywhere, fixing naive AR's ~0.18-0.20 over-rejection under
# degeneracy) but CONSERVATIVE near the boundary, and this is FUNDAMENTAL: c^2=n*zeta1
# is weakly identified there (its estimation error is O(1) relative to its value), so
# the valid upper CI is necessarily wide. The oracle (known c^2) achieves exact 0.05;
# no data-driven calibration can, near the boundary. -> ship basic AR + degeneracy
# flag; report the (honestly wide) sup-CI interval only when flagged.
ncores<-max(1L,min(4L,parallel::detectCores()))
a1<-0.025; a2<-0.025; z1<-qnorm(1-a1)
ghat<-function(A,eps){n<-length(A);SA<-sum(A);(eps*((n-2)*A+SA)+A*(SA-A))/(n-1)}
c2point<-function(A,eps){n<-length(A);g<-ghat(A,eps);m2<-mean((g-mean(g))^2)
  d2<-mean(A^2)^2; max(n*(m2-d2/(n-1)),0)}
one<-function(A,eps,B=4000,Bb=300,r=6){n<-length(A)
 P<-outer(A,A)+eps*(outer(A,rep(1,n))+outer(rep(1,n),A)); diag(P)<-0
 Un<-sum(P[upper.tri(P)])/choose(n,2); gi<-rowSums(P)/(n-1); m2<-mean((gi-mean(gi))^2)
 Pt<-P-outer(gi,rep(1,n))-outer(rep(1,n),gi)+Un; diag(Pt)<-0
 ev<-eigen(Pt,symmetric=TRUE,only.values=TRUE)$values/n; lam<-ev[order(abs(ev),decreasing=TRUE)][1:r]
 cb<-replicate(Bb,c2point(A[sample(n,n,replace=TRUE)],eps)); c2hi<-as.numeric(quantile(cb,1-a1))
 simq<-function(cc,ll,q){Z0<-rnorm(B);Zk<-matrix(rnorm(B*length(ll)),B,length(ll))
   as.numeric(quantile((2*cc*Z0+as.numeric((Zk^2-1)%*%ll))^2,q))}
 nUn2<-(n*Un)^2
 c(naive=as.numeric(n*Un^2/(4*m2)>3.841),
   oracle=as.numeric(nUn2>simq(sqrt(n)*eps,1,0.95)),
   supCI =as.numeric(nUn2>simq(sqrt(c2hi),lam,1-a2)))}
cat(sprintf("Andrews-Cheng sup-CI (percentile bootstrap c2-limit, a1=%.3f,a2=%.3f); H0 size target 0.05, n=300, R=1200\n",a1,a2))
cat(sprintf("%-7s %-9s %-9s %-9s %-9s\n","eps","zeta1","naive","oracle","sup-CI"))
for(eps in c(1.0,0.3,0.1,0.05,0.02,0.0)){
 res<-do.call(rbind,parallel::mclapply(1:1200,function(r){set.seed(20000L+as.integer(eps*1000)+r);one(rnorm(300),eps)},mc.cores=ncores))
 cat(sprintf("%-7.2f %-9.4f %-9.3f %-9.3f %-9.3f\n",eps,eps^2,mean(res[,"naive"]),mean(res[,"oracle"]),mean(res[,"supCI"])))}
