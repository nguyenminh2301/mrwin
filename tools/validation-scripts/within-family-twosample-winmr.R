#!/usr/bin/env Rscript
# TWO-SAMPLE / SUMMARY-DATA WITHIN-FAMILY win-ratio MR -- the programme capstone, chaining
# Papers 02 (IVW) + 03 (weak-IV-robust AR) into the within-family (Paper 04) design.
# =============================================================================================
# IDEA. The family design changes ONLY how the per-SNP summary stats are computed; the pooling
# is the EXISTING mrwin_twosample_ivw / mrwin_winmr_ar. Per SNP l:
#   beta_GX,l^wf = WITHIN-FAMILY (sib-mean-centered) slope of exposure X on dosage G_l
#                  (= exactly what a within-sibship GWAS, e.g. Howe 2022, produces)
#   delta_l^wf   = WITHIN-FAMILY slope of the marginal WIN-SCORE phenotype w(O_i)=E_j h(O_i,O_j)
#                  on dosage G_l.
# The win-score is a per-subject phenotype computed ONCE over the whole sample (per-subject
# kernel) then GWAS'd within family -> deployable on standard within-sibship GWAS pipelines.
# Family-mean-centering purges dynastic/AM (family-level) confounding; the within-family slope
# isolates the Mendelian path G_l -> X -> w. Feed (beta^wf, delta^wf, se) to IVW and AR.
#
# RESULTS (R 4.3.3; public package + simulator only; estimators use observed data).
# (A) Main, F=2500 families x 2 sibs = 5000, M=30 SNPs, R=80 -- vs interventional oracle:
#       mode  oracle | within-fam IVW (bias)   IVWcov | AR cov  Q-typeI | NAIVE pop IVW (bias)
#       none  0.114  | 0.043 (-0.071)          0.825  | 0.963   0.000   | 0.074 (-0.040)
#       dyn   0.093  | 0.037 (-0.056)          0.850  | 0.912   0.013   | -0.514 (-0.607)
#       am    0.094  | 0.044 (-0.050)          0.900  | 0.975   0.000   | -0.803 (-0.897)
#     => the WITHIN-FAMILY design removes dynastic/AM confounding (naive IVW is catastrophically
#        biased, -0.61 / -0.90; within-fam is not). But the within-family IVW POINT is attenuated
#        and its Wald CI UNDER-covers (0.82-0.90) -- the classic WEAK-INSTRUMENT failure of IVW.
#        The AR (Paper 03) is calibrated (0.91-0.98) and is the recommended inference; the over-ID
#        Q type-I is conservative (0.00-0.013, also a weak-IV effect).
# (B) F-ladder (mode=none) -- the attenuation is finite-sample weak-IV, not structural:
#       F      mean per-SNP F-stat   within-fam IVW bias   AR cover
#       2500   1.58                  -0.067                0.960
#       6000   2.24                  -0.047                1.000
#       12000  4.02                  -0.040                1.000
#     => as F grows the per-SNP F-stat grows and the IVW bias shrinks toward 0 (consistency);
#        within-family instruments are intrinsically WEAK (F-stat << 10 even at 24k subjects),
#        which is precisely WHY within-family win-MR needs the AR. AR coverage 0.96 -> 1.00
#        (valid, conservative as the AR set widens honestly under weak IV).
#
# HEADLINE: the two-sample within-family win-MR form is deployable on within-sibship GWAS
# summary statistics and unconfounded by dynastic/AM; report the AR confidence set (NOT the
# weak-IV-attenuated IVW point). Ties Papers 02 + 03 + 04 into one usable method.
# =============================================================================================
suppressMessages(library(mrwin)); ncores<-max(1L,min(4L,parallel::detectCores())); cfg<-mrwin_config()
subj<-getFromNamespace("mrwin_subject_win_loss_cpp","mrwin")
outc<-function(X,u,w,ex,cn,uf){F<-length(X);et<-matrix(NA_real_,F,3);lw<-log(w)
 for(j in 1:3){lp<-cfg$alpha_x[j]*X+cfg$nu_u[j]*u+ex+lw; et[,j]<-(-log(uf[,j])/(cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)}
 td<-et[,1];th<-et[,2];tr<-et[,3]
 list(time=cbind(pmin(td,cn),pmin(th,cn,td),pmin(tr,cn,td)),status=cbind(as.integer(td<=cn),as.integer(th<=cn&th<=td),as.integer(tr<=cn&tr<=td)))}
tr1<-function(G) ifelse(G==2,1L,ifelse(G==0,0L,rbinom(length(G),1,0.5)))
mkarch<-function(seed=11L,M=30L){set.seed(seed); list(base=runif(M,.15,.35),betas=rnorm(M,0,.15),M=M)}
genSib<-function(Fn,seed,mode,arch,gdyn=1.2,amc=0.12){set.seed(seed);M<-arch$M;base<-arch$base;betas<-arch$betas
 amf<-if(mode=="am") rnorm(Fn) else rep(0,Fn)
 pp<-pmin(pmax(matrix(base,Fn,M,byrow=TRUE)+outer(amf,amc*sign(betas)),.02),.95)
 Gf<-matrix(rbinom(Fn*M,2,pp),Fn,M); Gm<-matrix(rbinom(Fn*M,2,pp),Fn,M); mp<-as.numeric(((Gf+Gm)/2)%*%betas)
 G<-NULL;fam<-NULL;X<-NULL;U<-NULL;Wf<-NULL;EX<-NULL;CN<-NULL;UF<-NULL
 for(k in 1:2){Gs<-matrix(tr1(Gf),Fn,M)+matrix(tr1(Gm),Fn,M); Z<-as.numeric(Gs%*%betas)
  u<-rnorm(Fn); x<-cfg$alpha_s*Z+cfg$alpha_u*u+rnorm(Fn); ex<-if(mode%in%c("dyn","am")) gdyn*scale(mp)[,1] else rep(0,Fn)
  wf<-rgamma(Fn,1/cfg$theta_f,scale=cfg$theta_f); cn<-pmin(rexp(Fn,cfg$censoring_rate),cfg$max_follow_up); uf<-matrix(runif(Fn*3),Fn,3)
  G<-rbind(G,Gs);fam<-c(fam,1:Fn);X<-c(X,x);U<-c(U,u);Wf<-c(Wf,wf);EX<-c(EX,ex);CN<-c(CN,cn);UF<-rbind(UF,uf)}
 oc<-outc(X,U,Wf,EX,CN,UF); list(fam=fam,G=G,X=X,u=U,w=Wf,ex=EX,cn=CN,uf=UF,time=oc$time,status=oc$status,M=M)}
oracle<-function(d,eps=0.25){np<-length(d$X)
 Bp<-outc(d$X+eps,d$u,d$w,d$ex,d$cn,d$uf); Bm<-outc(d$X-eps,d$u,d$w,d$ex,d$cn,d$uf)
 v<-mrwin_fast_pair_win_loss(Bp$time,Bp$status,Bm$time,Bm$status)
 (as.numeric(v[["wins"]])-as.numeric(v[["losses"]]))/(as.numeric(np)^2)/(2*eps)}
winscore<-function(d){c<-subj(d$time,d$status,d$time,d$status,rep(1,length(d$X))); (c[,1]-c[,2])/(length(d$X)-1)}
# per-SNP within-family ("fam") or naive-population ("glob") slopes of (X, win-score) on each SNP
gwas<-function(d,wv,ctr){N<-length(d$X);G<-d$G;M<-d$M
 if(ctr=="fam"){cnt<-as.vector(table(d$fam))
   cX<-as.numeric(d$X-tapply(d$X,d$fam,mean)[d$fam]); cW<-as.numeric(wv-tapply(wv,d$fam,mean)[d$fam])
   fmG<-rowsum(G,d$fam)/cnt; cG<-G-fmG[d$fam,]; df<-N-length(cnt)-1L
 } else {cX<-as.numeric(d$X-mean(d$X)); cW<-as.numeric(wv-mean(wv)); cG<-G-matrix(colMeans(G),N,M,byrow=TRUE); df<-N-2L}
 den<-colSums(cG^2); numW<-colSums(cG*cW); numX<-colSums(cG*cX)
 rssW<-sum(cW^2)-numW^2/den; rssX<-sum(cX^2)-numX^2/den
 list(beta=numX/den,se_beta=sqrt((rssX/df)/den),delta=numW/den,se_delta=sqrt((rssW/df)/den))}

cat("== (A) MAIN: F=2500 fam x2 sibs=5000, M=30 SNPs, R=80 -- vs interventional oracle ==\n")
cat(sprintf("%-6s %-8s | %-24s %-7s | %-9s %-8s | %-18s\n","mode","oracle","within-fam IVW (bias)","IVWcov","AR cov","Q-typeI","NAIVE pop IVW (bias)"))
for(mode in c("none","dyn","am")){arch<-mkarch(); R<-80L
 res<-do.call(rbind,parallel::mclapply(1:R,function(r){d<-genSib(2500L,6000L+r,mode,arch); o<-oracle(d); wv<-winscore(d)
   gf<-gwas(d,wv,"fam"); gn<-gwas(d,wv,"glob")
   iv<-mrwin_twosample_ivw(gf$beta,gf$delta,gf$se_delta,se_gx=gf$se_beta)
   ar<-mrwin_winmr_ar(gf$beta,gf$delta,gf$se_delta,se_gx=gf$se_beta)
   ivn<-mrwin_twosample_ivw(gn$beta,gn$delta,gn$se_delta,se_gx=gn$se_beta)
   c(o=o,bf=iv$gamma,cf=as.numeric(o>=iv$ci95[1]&&o<=iv$ci95[2]),
     arc=as.numeric(ar$contiguous&&o>=ar$ci[1]&&o<=ar$ci[2]),qp=as.numeric(ar$Q_p<0.05),bn=ivn$gamma)},mc.cores=ncores))
 m<-colMeans(res)
 cat(sprintf("%-6s %-8.3f | %-24s %.3f   | %.3f     %.3f    | %.3f (%+.3f)\n",mode,m["o"],
   sprintf("%.3f (%+.3f)",m["bf"],m["bf"]-m["o"]),m["cf"],m["arc"],m["qp"],m["bn"],m["bn"]-m["o"]))}

cat("\n== (B) F-LADDER (mode=none): IVW attenuation is finite-sample weak-IV ==\n")
cat(sprintf("%-7s %-8s %-22s %-10s %-10s\n","F","oracle","within-fam IVW (bias)","AR cover","mean F-stat"))
for(spec in list(c(2500,50),c(6000,40),c(12000,25))){Fn<-as.integer(spec[1]);R<-as.integer(spec[2]);arch<-mkarch()
 res<-do.call(rbind,parallel::mclapply(1:R,function(r){d<-genSib(Fn,6000L+r,"none",arch); o<-oracle(d); wv<-winscore(d)
   g<-gwas(d,wv,"fam"); iv<-mrwin_twosample_ivw(g$beta,g$delta,g$se_delta,se_gx=g$se_beta)
   ar<-mrwin_winmr_ar(g$beta,g$delta,g$se_delta,se_gx=g$se_beta)
   c(o=o,b=iv$gamma,arc=as.numeric(ar$contiguous&&o>=ar$ci[1]&&o<=ar$ci[2]),fs=mean(g$beta^2/g$se_beta^2))},mc.cores=ncores))
 m<-colMeans(res)
 cat(sprintf("%-7d %-8.3f %-22s %.3f      %.2f\n",Fn,m["o"],sprintf("%.3f (%+.3f)",m["b"],m["b"]-m["o"]),m["arc"],m["fs"]))}
cat("DONE2SWF\n")
