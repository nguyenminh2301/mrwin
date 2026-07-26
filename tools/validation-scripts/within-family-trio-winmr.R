#!/usr/bin/env Rscript
# PARENT-OFFSPRING TRIO win-ratio MR -- a genuinely DIFFERENT structure from within-sibship.
# =============================================================================================
# STRUCTURE. Sibship: contrast is WITHIN-family (sib i vs sib j), dZ=Z_i-Z_j -> a design-
# RESTRICTED (incomplete) U-statistic over sib pairs only. Trio: the Mendelian instrument is the
# per-offspring MID-PARENT RESIDUAL g_o = Z_o - (Z_f+Z_m)/2 (pure meiotic segregation; mean-zero
# in every family AND every subpopulation, because E[Z_o|parents] = mid-parent PRS). g_o is thus
# clean of dynastic, assortative-mating, and population-stratification confounding -- all of which
# act through the PARENTAL background. Because the instrument is cleaned at the INDIVIDUAL level,
# the win comparison runs over ALL offspring pairs -> a COMPLETE U-statistic (more pairs than the
# sibship's within-family-only pairs).
#
# ESTIMATOR (O(N^2 K), reuses the validated per-subject kernel). With antisymmetric h,
#   sum_{a<b} h(O_a,O_b)(g_a-g_b) = sum_a g_a * W_a,   W_a = net win count of offspring a vs all.
# So beta_trio = Cov(g, w(O)) / Cov(g, X),  w(O)_a = W_a/(N-1) (per-subject win-loss / kernel).
# Naive comparator = same complete-U but instrument = offspring PRS Z_o (NO mid-parent
# correction) -> confounded by dynastic/AM/stratification. Oracle = interventional do(X+-eps).
#
# INFERENCE (complete-U influence function; no bootstrap). Perturbing obs a changes both
# g_a*w(O_a) and every w(O_b) averaging over a, giving
#   IF_num,a = g_a w(O_a) - r(O_a) - 2 Cov(g,w),  r(O_a)=E_b[g_b h(O_a,O_b)].
# r(O_a) = R_a/(N-1) with R_a = sum_b g_b h_ab is a SECOND kernel call with weights = g (signed).
# Ratio IF per offspring: phi_a = (1/D)[(psi_num,a-mean) - beta(g_a X_a - mean)],
# psi_num,a = g_a w(O_a) - r(O_a), D = Cov(g,X); Var(beta) = sum phi_a^2 / N^2; CI = beta +-1.96 se.
#
# RESULTS (R 4.3.3; public package + simulator only; estimators use observed data).
# (A) Robustness, F=3000 trios, R=80 -- bias vs interventional oracle:
#       mode    oracle  trio g=Zo-midpar (bias)   naive Zo complete-U (bias)
#       none    0.112   +0.007 (SE .008)          +0.002
#       dyn     0.092   +0.002                    -1.032
#       am      0.094   +0.001                    -0.722
#       strat   0.107   +0.009                    -0.353
#     => trio mid-parent estimator is UNBIASED under dynastic/AM/stratification; the naive
#        offspring-PRS complete-U is catastrophically confounded. (PCs would handle 'strat' but
#        NOT dyn/am -- the PC-irreducible confounders are the trio's unique value, mirroring the
#        sibship scope correction in findings-paper04-litscan.md / research-frontier-roadmap.md.)
# (B) Consistency (convergence-in-N ladder, mode=dyn): bias 1500->3000->6000 = +0.013->+0.011->
#     +0.007 -> shrinks toward 0 (within-trio weak-IV finite-sample bias; consistent, lesson B).
# (C) Calibrated inference (F=3000, R=90): complete-U IF Wald-CI coverage
#       none 0.978 | dyn 0.944 | am 0.967 | strat 0.944  (mean se 0.067-0.072 vs emp sd
#     0.071-0.078). The IF se slightly under-estimates in confounded scenarios (ratio ~0.90),
#     attributable to the omitted higher-order ZETA_2 U-statistic term; coverage holds at nominal.
# =============================================================================================
suppressMessages(library(mrwin)); ncores<-max(1L,min(4L,parallel::detectCores())); cfg<-mrwin_config()
subj<-getFromNamespace("mrwin_subject_win_loss_cpp","mrwin")
outc<-function(X,u,w,ex,cn,uf){F<-length(X);et<-matrix(NA_real_,F,3);lw<-log(w)
 for(j in 1:3){lp<-cfg$alpha_x[j]*X+cfg$nu_u[j]*u+ex+lw; et[,j]<-(-log(uf[,j])/(cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)}
 td<-et[,1];th<-et[,2];tr<-et[,3]
 list(time=cbind(pmin(td,cn),pmin(th,cn,td),pmin(tr,cn,td)),status=cbind(as.integer(td<=cn),as.integer(th<=cn&th<=td),as.integer(tr<=cn&tr<=td)))}
tr1<-function(G) ifelse(G==2,1L,ifelse(G==0,0L,rbinom(length(G),1,0.5)))   # transmit one allele
genTrio<-function(Fn,seed,mode,gdyn=1.2,amc=0.12,fdiff=0.12,strat=1.0){set.seed(seed);M<-30
 base<-runif(M,.15,.35); betas<-rnorm(M,0,.15)
 amf<-if(mode=="am") rnorm(Fn) else rep(0,Fn); sub<-if(mode=="strat") rep(c(0,1),length.out=Fn) else rep(0,Fn)
 pp<-pmin(pmax(matrix(base,Fn,M,byrow=TRUE)+outer(amf,amc*sign(betas))+outer(sub,fdiff*sign(betas)),.02),.95)
 Gf<-matrix(rbinom(Fn*M,2,pp),Fn,M); Gm<-matrix(rbinom(Fn*M,2,pp),Fn,M)
 Zf<-as.numeric(Gf%*%betas); Zm<-as.numeric(Gm%*%betas); mp<-(Zf+Zm)/2
 Go<-matrix(tr1(Gf),Fn,M)+matrix(tr1(Gm),Fn,M); Zo<-as.numeric(Go%*%betas); g<-Zo-mp
 u<-rnorm(Fn); X<-cfg$alpha_s*Zo+cfg$alpha_u*u+rnorm(Fn)
 ex<-(if(mode%in%c("dyn","am")) gdyn*scale(mp)[,1] else 0)+(if(mode=="strat") strat*(sub-mean(sub)) else 0)
 w<-rgamma(Fn,1/cfg$theta_f,scale=cfg$theta_f); cn<-pmin(rexp(Fn,cfg$censoring_rate),cfg$max_follow_up); uf<-matrix(runif(Fn*3),Fn,3)
 oc<-outc(X,u,w,ex,cn,uf); list(Zo=Zo,g=g,X=X,u=u,w=w,ex=ex,cn=cn,uf=uf,time=oc$time,status=oc$status)}
oracle<-function(d,eps=0.25){np<-length(d$X)
 Bp<-outc(d$X+eps,d$u,d$w,d$ex,d$cn,d$uf); Bm<-outc(d$X-eps,d$u,d$w,d$ex,d$cn,d$uf)
 v<-mrwin_fast_pair_win_loss(Bp$time,Bp$status,Bm$time,Bm$status)
 (as.numeric(v[["wins"]])-as.numeric(v[["losses"]]))/(as.numeric(np)^2)/(2*eps)}
wsc<-function(d){c<-subj(d$time,d$status,d$time,d$status,rep(1,nrow(d$time))); (c[,1]-c[,2])/(nrow(d$time)-1)}
est<-function(d,iv) as.numeric(cov(iv,wsc(d))/cov(iv,d$X))                     # point (any instrument)
estIF<-function(d){N<-length(d$X);g<-d$g;X<-d$X                                # point + complete-U IF se
 cu<-subj(d$time,d$status,d$time,d$status,rep(1,N)); wsc<-(cu[,1]-cu[,2])/(N-1)
 cg<-subj(d$time,d$status,d$time,d$status,g);        r  <-(cg[,1]-cg[,2])/(N-1)
 Nnum<-mean(g*wsc)-mean(g)*mean(wsc); Dden<-mean(g*X)-mean(g)*mean(X); beta<-Nnum/Dden
 psin<-g*wsc-r; psid<-g*X; phi<-((psin-mean(psin))-beta*(psid-mean(psid)))/Dden
 list(beta=beta, se=sqrt(sum(phi^2)/N^2))}

cat("== (A) ROBUSTNESS: F=3000 trios, R=80 -- bias vs interventional oracle ==\n")
cat(sprintf("%-8s %-9s %-22s %-22s\n","mode","oracle","trio g=Zo-midpar (bias)","naive Zo complete-U (bias)"))
for(mode in c("none","dyn","am","strat")){R<-80L
 res<-do.call(rbind,parallel::mclapply(1:R,function(r){d<-genTrio(3000L,3000L+r,mode)
   c(o=oracle(d),bg=est(d,d$g),bz=est(d,d$Zo))},mc.cores=ncores))
 m<-colMeans(res); seg<-sd(res[,"bg"])/sqrt(R)
 cat(sprintf("%-8s %-9.3f %-22s %-22s\n",mode,m["o"],sprintf("%.3f (%+.3f, SE%.3f)",m["bg"],m["bg"]-m["o"],seg),
   sprintf("%.3f (%+.3f)",m["bz"],m["bz"]-m["o"])))}

cat("\n== (B) CONSISTENCY: convergence-in-N ladder, mode=dyn ==\n")
cat(sprintf("%-7s %-9s %-20s\n","F","oracle","trio-g point (bias)"))
for(Fn in c(1500L,3000L,6000L)){R<-60L
 res<-do.call(rbind,parallel::mclapply(1:R,function(r){d<-genTrio(Fn,3000L+r,"dyn"); c(o=oracle(d),b=estIF(d)$beta)},mc.cores=ncores))
 m<-colMeans(res); se<-sd(res[,"b"])/sqrt(R)
 cat(sprintf("%-7d %-9.3f %.3f (%+.3f, SE %.3f)\n",Fn,m["o"],m["b"],m["b"]-m["o"],se))}

cat("\n== (C) CALIBRATED INFERENCE: complete-U IF-variance Wald-CI coverage (F=3000,R=90) ==\n")
cat(sprintf("%-8s %-9s %-10s %-10s %-14s\n","mode","oracle","mean se","emp sd","95% coverage"))
for(mode in c("none","dyn","am","strat")){R<-90L
 res<-do.call(rbind,parallel::mclapply(1:R,function(r){d<-genTrio(3000L,4200L+r,mode)
   o<-oracle(d); e<-estIF(d); c(o=o,b=e$beta,se=e$se,cov=as.numeric(abs(e$beta-o)<=1.96*e$se))},mc.cores=ncores))
 m<-colMeans(res); cse<-sqrt(m["cov"]*(1-m["cov"])/R)
 cat(sprintf("%-8s %-9.3f %-10.4f %-10.4f %.3f (SE %.3f)\n",mode,m["o"],m["se"],sd(res[,"b"]),m["cov"],cse))}
cat("DONETRIO\n")
