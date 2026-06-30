#!/usr/bin/env Rscript
# Reproduces the figures for the Paper 04 (within-family win-ratio MR) manuscript,
# entirely from the public package + simulator -- no confidential data. Mirrors the
# tools/paper-figures.R convention (Paper 01). Run from the repository root after
# installing mrwin:
#
#     Rscript tools/paper04-figures.R
#
# By default writes vector PDFs to papers/04-within-family-winmr/figures/ (gitignored
# until submission). Override with MRWIN_FIG_DIR=some/dir. Select a subset of figures
# (useful for staged runs) with MRWIN_FIG_WHICH="1,2,3" (1=confounder bars, 2=C^2
# portrait, 3=weak-IV ladder; default = all). Seeds are fixed for bit-reproducibility.
#
# Each block REUSES the exact validated data-generating processes from the committed
# probes (lesson G: read/reuse, do not refork validated machinery):
#   fig 1 <- merges tools/validation-scripts/within-family-winmr-probe.R (naive/PC/
#            within-family estimators, stratification) with the AM generator of
#            within-family-ar-inference.R (continuous family-level mating latent).
#   fig 2 <- tools/figures/c2-portrait.R (PNG prototype), re-targeted to a vector PDF.
#   fig 3 <- the F-ladder of tools/validation-scripts/within-family-twosample-winmr.R
#            (part B), re-run live and plotted with MC-error bars.
suppressMessages(library(mrwin))
outdir <- Sys.getenv("MRWIN_FIG_DIR", unset = file.path("papers", "04-within-family-winmr", "figures"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pdf_open <- function(f, w = 7, h = 5) pdf(file.path(outdir, f), width = w, height = h)
which_figs <- Sys.getenv("MRWIN_FIG_WHICH", unset = "1,2,3")
which_figs <- as.integer(strsplit(which_figs, ",")[[1]])
ncores <- max(1L, min(4L, parallel::detectCores()))
cfg <- mrwin_config()
subj <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")

hpair <- function(t1,s1,t2,s2){F<-nrow(t1);res<-integer(F)
 for(k in 1:3){dec<-res!=0; iwin<-(s2[,k]==1 & t1[,k]>t2[,k]); ilos<-(s1[,k]==1 & t2[,k]>t1[,k])
  res[!dec & iwin]<-1L; res[!dec & ilos & !iwin]<- -1L}; res}
outcome <- function(X,u,w,extra,cn,uf){F<-length(X);et<-matrix(NA_real_,F,3);lw<-log(w)
 for(j in 1:3){lp<-cfg$alpha_x[j]*X+cfg$nu_u[j]*u+extra+lw
  et[,j]<-(-log(uf[,j])/(cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)}
 td<-et[,1];th<-et[,2];tr<-et[,3]
 list(time=cbind(pmin(td,cn),pmin(th,cn,td),pmin(tr,cn,td)),
      status=cbind(as.integer(td<=cn),as.integer(th<=cn&th<=td),as.integer(tr<=cn&tr<=td)))}
transmit <- function(G) ifelse(G==2,1,ifelse(G==0,0,rbinom(length(G),1,0.5)))

# ============================================================================
# FIGURE 1: confounder-robustness bars (naive-pop / PC-adj-pop / within-family
# x none / stratification / dynastic / assortative-mating), with MC-SE error bars.
# Generator MERGES the validated stratification DGP (within-family-winmr-probe.R)
# with the validated AM DGP (within-family-ar-inference.R) -- same building blocks,
# not a new model.
# ============================================================================
if (1L %in% which_figs) {
gen_fam <- function(Fn,seed,mode,gstrat=2.5,fgap=0.18,gdyn=1.2,amc=0.10){set.seed(seed);M<-30
 base<-runif(M,0.15,0.35); betas<-rnorm(M,0,0.15)
 strat<-if(mode=="strat") rbinom(Fn,1,0.5) else rep(0,Fn)
 amf<-if(mode=="am") rnorm(Fn) else rep(0,Fn)
 fg<-if(mode=="strat") fgap else 0
 p<-matrix(base,Fn,M,byrow=TRUE)+outer(strat,fg*sign(betas))+outer(amf,amc*sign(betas))
 p<-pmin(pmax(p,0.02),0.95)
 Gf<-matrix(rbinom(Fn*M,2,p),Fn,M); Gm<-matrix(rbinom(Fn*M,2,p),Fn,M)
 parPRS<-as.numeric(((Gf+Gm)/2)%*%betas)
 sibg<-function() t(apply(Gf,1,transmit))+t(apply(Gm,1,transmit))
 mk<-function(Gs){Z<-as.numeric(Gs%*%betas)
   u<-rnorm(Fn); X<-cfg$alpha_s*Z+cfg$alpha_u*u+rnorm(Fn)
   extra<-(if(mode=="strat") gstrat*strat else 0)+(if(mode%in%c("dyn","am")) gdyn*scale(parPRS)[,1] else 0)
   w<-rgamma(Fn,1/cfg$theta_f,scale=cfg$theta_f); cn<-pmin(rexp(Fn,cfg$censoring_rate),cfg$max_follow_up)
   uf<-matrix(runif(Fn*3),Fn,3); oc<-outcome(X,u,w,extra,cn,uf)
   list(Z=Z,X=X,u=u,w=w,extra=extra,cn=cn,uf=uf,G=Gs,time=oc$time,status=oc$status)}
 list(s1=mk(sibg()),s2=mk(sibg()))}
winscore <- function(d){ti<-rbind(d$s1$time,d$s2$time);si<-rbind(d$s1$status,d$s2$status)
 cc<-subj(ti,si,ti,si,rep(1,nrow(ti))); (cc[,1]-cc[,2])/(nrow(ti)-1)}
naive_pop <- function(d,wsc){Z<-c(d$s1$Z,d$s2$Z);X<-c(d$s1$X,d$s2$X); cov(Z,wsc)/cov(Z,X)}
pc_pop <- function(d,wsc,k=10){G<-rbind(d$s1$G,d$s2$G); Z<-c(d$s1$Z,d$s2$Z);X<-c(d$s1$X,d$s2$X)
 pc<-prcomp(G,center=TRUE,scale.=FALSE)$x[,1:k]
 rz<-resid(lm(Z~pc)); rw<-resid(lm(wsc~pc)); rx<-resid(lm(X~pc)); cov(rz,rw)/cov(rz,rx)}
wf_est <- function(d){h<-hpair(d$s1$time,d$s1$status,d$s2$time,d$s2$status)
 dZ<-d$s1$Z-d$s2$Z; dX<-d$s1$X-d$s2$X; sum(h*dZ)/sum(dX*dZ)}
oracle1 <- function(d,eps=0.25){A<-d$s1;np<-length(A$X)
 Bp<-outcome(A$X+eps,A$u,A$w,A$extra,A$cn,A$uf); Bm<-outcome(A$X-eps,A$u,A$w,A$extra,A$cn,A$uf)
 v<-mrwin_fast_pair_win_loss(Bp$time,Bp$status,Bm$time,Bm$status)
 (as.numeric(v[["wins"]])-as.numeric(v[["losses"]]))/(as.numeric(np)^2)/(2*eps)}

modes <- c("none","strat","dyn","am"); R1 <- 30L; Fn1 <- 4000L
cat(sprintf("[fig1] confounder-robustness bars: F=%d, R=%d per scenario\n", Fn1, R1))
bias_m <- bias_se <- matrix(NA_real_, 3, length(modes), dimnames=list(c("naive","pc","wf"), modes))
for (mi in seq_along(modes)) {
  mode <- modes[mi]
  res <- do.call(rbind, parallel::mclapply(1:R1, function(r) {
    d <- gen_fam(Fn1, 1000L + r, mode); ws <- winscore(d)
    c(orc=oracle1(d), nv=naive_pop(d,ws), pc=pc_pop(d,ws), wf=wf_est(d))
  }, mc.cores = ncores))
  b_nv <- res[,"nv"]-res[,"orc"]; b_pc <- res[,"pc"]-res[,"orc"]; b_wf <- res[,"wf"]-res[,"orc"]
  bias_m["naive",mi] <- mean(b_nv); bias_m["pc",mi] <- mean(b_pc); bias_m["wf",mi] <- mean(b_wf)
  bias_se["naive",mi] <- sd(b_nv)/sqrt(R1); bias_se["pc",mi] <- sd(b_pc)/sqrt(R1); bias_se["wf",mi] <- sd(b_wf)/sqrt(R1)
  cat(sprintf("  %-6s naive %+.3f(SE%.3f)  PC-adj %+.3f(SE%.3f)  within-fam %+.3f(SE%.3f)\n",
    mode, bias_m["naive",mi], bias_se["naive",mi], bias_m["pc",mi], bias_se["pc",mi], bias_m["wf",mi], bias_se["wf",mi]))
}
pdf_open("fig-confounder-bars.pdf", w = 7.5, h = 5)
par(mar=c(4.5,4.6,3,1))
cols <- c(naive="#C0392B", pc="#E08E45", wf="#2A6FB0")
bp <- barplot(bias_m, beside=TRUE, col=cols[rownames(bias_m)], ylim=c(min(bias_m-2*bias_se,-0.05)*1.1, max(bias_m+2*bias_se,0.05)*1.1+0.05),
  names.arg=c("none\n(control)","stratification","dynastic","assortative\nmating"),
  ylab="bias vs. interventional oracle  (estimate - truth)",
  main="Confounder robustness: naive vs. PC-adjusted vs. within-family win-MR", cex.main=1.0, cex.names=0.85)
arrows(bp, bias_m-bias_se, bp, bias_m+bias_se, angle=90, code=3, length=0.04, lwd=1.3)
abline(h=0, col="grey40", lty=2)
legend("bottomleft", bty="n", fill=cols, legend=c("naive population (no adjustment)","PC-adjusted population (standard MR control)","within-family (this paper)"), cex=0.82)
dev.off()
cat("[fig1] wrote fig-confounder-bars.pdf\n")
}

# ============================================================================
# FIGURE 2: Concordance-Contrast (C^2) portrait -- flagship figure, vector PDF.
# Identical construction to tools/figures/c2-portrait.R (PNG prototype).
# ============================================================================
if (2L %in% which_figs) {
hpair_tier <- function(t1,s1,t2,s2){F<-nrow(t1);r<-integer(F);tier<-integer(F)
 for(k in 1:3){iw<-tier==0 & s2[,k]==1 & t1[,k]>t2[,k]; il<-tier==0 & !iw & s1[,k]==1 & t2[,k]>t1[,k]
  r[iw]<-1L; tier[iw]<-k; r[il]<- -1L; tier[il]<-k}; list(h=r,tier=tier)}
genSib2 <- function(Fn,seed,gdyn=1.4){set.seed(seed);M<-30;base<-runif(M,.15,.35);betas<-rnorm(M,0,.15)
 Gf<-matrix(rbinom(Fn*M,2,base[col(matrix(0,Fn,M))]),Fn,M); Gm<-matrix(rbinom(Fn*M,2,base[col(matrix(0,Fn,M))]),Fn,M)
 mp<-as.numeric(((Gf+Gm)/2)%*%betas); sib<-vector("list",2); o1<-NULL
 for(k in 1:2){Gs<-matrix(transmit(Gf),Fn,M)+matrix(transmit(Gm),Fn,M); Z<-as.numeric(Gs%*%betas)
  u<-rnorm(Fn); X<-cfg$alpha_s*Z+cfg$alpha_u*u+rnorm(Fn); ex<-gdyn*scale(mp)[,1]
  w<-rgamma(Fn,1/cfg$theta_f,scale=cfg$theta_f); cn<-pmin(rexp(Fn,cfg$censoring_rate),cfg$max_follow_up); uf<-matrix(runif(Fn*3),Fn,3)
  oc<-outcome(X,u,w,ex,cn,uf); sib[[k]]<-list(Z=Z,X=X,time=oc$time,status=oc$status); if(k==1) o1<-list(u=u,w=w,ex=ex,cn=cn,uf=uf,X=X)}
 list(sib=sib,o1=o1)}
oracle2 <- function(d,eps=0.25){A<-d$o1;np<-length(A$X)
 Bp<-outcome(A$X+eps,A$u,A$w,A$ex,A$cn,A$uf); Bm<-outcome(A$X-eps,A$u,A$w,A$ex,A$cn,A$uf)
 v<-mrwin_fast_pair_win_loss(Bp$time,Bp$status,Bm$time,Bm$status)
 (as.numeric(v[["wins"]])-as.numeric(v[["losses"]]))/(as.numeric(np)^2)/(2*eps)}

Fn2 <- 12000L; d2 <- genSib2(Fn2, 2026L); s1 <- d2$sib[[1]]; s2 <- d2$sib[[2]]
ht <- hpair_tier(s1$time,s1$status,s2$time,s2$status)
dZ <- s1$Z-s2$Z; dX <- s1$X-s2$X; h <- ht$h; tier <- ht$tier
sh <- c(2:Fn2,1L); htn <- hpair_tier(s1$time,s1$status,s1$time[sh,],s1$status[sh,]); dZn <- s1$Z-s1$Z[sh]; hn <- htn$h
sX <- sum(dZ*dX)/sum(dZ^2); sH <- sum(dZ*h)/sum(dZ^2); bhat <- sH/sX
ARv2 <- function(b){m<-(h-b*dX)*dZ; sum(m)^2/sum(m^2)}
span0 <- 6/sqrt(sum((dX*dZ)^2)/sum(((h-bhat*dX)*dZ)^2)+1e-9)
grid2 <- seq(bhat-span0, bhat+span0, length.out=4000)
acc2 <- grid2[vapply(grid2,ARv2,1)<=3.841]; blo<-min(acc2); bhi<-max(acc2)
borc <- oracle2(d2)
brk <- seq(-quantile(abs(dZn),.97), quantile(abs(dZn),.97), length.out=19); mid <- (brk[-1]+brk[-length(brk)])/2
binmean <- function(x,y){i<-findInterval(x,brk); m<-tapply(y,factor(i,levels=1:(length(brk)-1)),mean); n<-tapply(y,factor(i,levels=1:(length(brk)-1)),length); list(m=as.numeric(m),n=as.numeric(n))}
wfb <- binmean(dZ,h); nvb <- binmean(dZn,hn)
dec <- tier>0; ti <- findInterval(dZ,brk)
tshare <- sapply(1:(length(brk)-1), function(b){idx<-dec&ti==b; if(sum(idx)<5) return(c(NA,NA,NA)); tabulate(tier[idx],3)/sum(idx)})
cat(sprintf("[fig2] within-fam gradient %.3f, AR [%.2f,%.2f], oracle %.3f, naive slope %.3f\n",
  bhat, blo, bhi, borc, sum(dZn*hn)/sum(dZn^2)/sX))

pdf_open("fig-c2-portrait.pdf", w = 7.2, h = 7.5)
layout(matrix(c(1,2),2,1),heights=c(2.3,1)); par(mar=c(4.2,4.6,3.2,1.4))
xr <- range(brk)
plot(NA,xlim=xr,ylim=c(-.6,.6),xlab=expression(paste("instrument contrast  ",Delta,"Z (genetic)")),
     ylab=expression(paste("reduced-form win  ",hat(E),"[ h | ",Delta,"Z ]")),
     main=expression("Concordance-Contrast ("*C^2*") portrait  -  within-family win-MR (dynastic confounding)"),cex.main=1.0)
abline(h=0,col="grey80"); abline(v=0,col="grey80")
xx <- c(xr[1],xr[2]); polygon(c(xx,rev(xx)),c(blo*sX*xx,rev(bhi*sX*xx)),col=rgb(.20,.45,.80,.18),border=NA)
lines(xx,borc*sX*xx,lty=2,lwd=2,col="black")
lines(xx,sH*xx,lwd=2.4,col=rgb(.13,.35,.75))
pts <- function(b,col){ok<-is.finite(b$m); points(mid[ok],b$m[ok],pch=19,col=col,cex=pmin(2.4,.5+b$n[ok]/max(b$n,na.rm=TRUE)*2.2))}
pts(nvb,rgb(.85,.30,.25,.9)); lines(mid,nvb$m,col=rgb(.85,.30,.25,.7),lwd=1.5)
pts(wfb,rgb(.13,.35,.75,.95))
legend("topleft",bty="n",cex=.85,
  legend=c(sprintf("within-family spine  (gradient %.2f)",bhat),
           sprintf("AR 95%% fan  [%.2f, %.2f]  -> %s instrument",blo,bhi,ifelse(bhi-blo>2*abs(bhat),"WEAK","strong")),
           sprintf("oracle gradient  %.2f  (dashed)",borc), "naive between-family (confounded)"),
  pch=c(NA,15,NA,19),lty=c(1,NA,2,NA),lwd=c(2.4,NA,2,NA),
  col=c(rgb(.13,.35,.75),rgb(.20,.45,.80,.5),"black",rgb(.85,.30,.25)))
par(mar=c(4.4,4.6,1.2,1.4))
cols3 <- c(rgb(.20,.20,.25),rgb(.45,.55,.70),rgb(.75,.80,.90))
plot(NA,xlim=xr,ylim=c(0,1),xlab=expression(paste("instrument contrast  ",Delta,"Z")),ylab="share of decided wins")
cum <- rep(0,ncol(tshare))
for(k in 1:3){top<-cum+ifelse(is.na(tshare[k,]),0,tshare[k,]); ok<-is.finite(top)
  polygon(c(mid[ok],rev(mid[ok])),c(cum[ok],rev(top[ok])),col=cols3[k],border=NA); cum<-top}
legend("top",horiz=TRUE,bty="n",cex=.85,fill=cols3,legend=c("tier 1: death","tier 2: hospitalization","tier 3: recurrence"))
dev.off()
cat("[fig2] wrote fig-c2-portrait.pdf\n")
}

# ============================================================================
# FIGURE 3: weak-instrument ladder for the two-sample within-family form --
# re-runs the F-ladder of within-family-twosample-winmr.R part B, plotted with
# MC-error bars (IVW bias) and AR-vs-Wald coverage.
# ============================================================================
if (3L %in% which_figs) {
mkarch <- function(seed=11L,M=30L){set.seed(seed); list(base=runif(M,.15,.35),betas=rnorm(M,0,.15),M=M)}
genSib3 <- function(Fn,seed,mode,arch,gdyn=1.2,amc=0.12){set.seed(seed);M<-arch$M;base<-arch$base;betas<-arch$betas
 amf<-if(mode=="am") rnorm(Fn) else rep(0,Fn)
 pp<-pmin(pmax(matrix(base,Fn,M,byrow=TRUE)+outer(amf,amc*sign(betas)),.02),.95)
 Gf<-matrix(rbinom(Fn*M,2,pp),Fn,M); Gm<-matrix(rbinom(Fn*M,2,pp),Fn,M); mp<-as.numeric(((Gf+Gm)/2)%*%betas)
 G<-NULL;fam<-NULL;X<-NULL;U<-NULL;Wf<-NULL;EX<-NULL;CN<-NULL;UF<-NULL
 for(k in 1:2){Gs<-matrix(transmit(Gf),Fn,M)+matrix(transmit(Gm),Fn,M); Z<-as.numeric(Gs%*%betas)
  u<-rnorm(Fn); x<-cfg$alpha_s*Z+cfg$alpha_u*u+rnorm(Fn); ex<-if(mode%in%c("dyn","am")) gdyn*scale(mp)[,1] else rep(0,Fn)
  wf<-rgamma(Fn,1/cfg$theta_f,scale=cfg$theta_f); cn<-pmin(rexp(Fn,cfg$censoring_rate),cfg$max_follow_up); uf<-matrix(runif(Fn*3),Fn,3)
  G<-rbind(G,Gs);fam<-c(fam,1:Fn);X<-c(X,x);U<-c(U,u);Wf<-c(Wf,wf);EX<-c(EX,ex);CN<-c(CN,cn);UF<-rbind(UF,uf)}
 oc<-outcome(X,U,Wf,EX,CN,UF); list(fam=fam,G=G,X=X,u=U,w=Wf,ex=EX,cn=CN,uf=UF,time=oc$time,status=oc$status,M=M)}
oracle3 <- function(d,eps=0.25){np<-length(d$X)
 Bp<-outcome(d$X+eps,d$u,d$w,d$ex,d$cn,d$uf); Bm<-outcome(d$X-eps,d$u,d$w,d$ex,d$cn,d$uf)
 v<-mrwin_fast_pair_win_loss(Bp$time,Bp$status,Bm$time,Bm$status)
 (as.numeric(v[["wins"]])-as.numeric(v[["losses"]]))/(as.numeric(np)^2)/(2*eps)}
winscore3 <- function(d){c<-subj(d$time,d$status,d$time,d$status,rep(1,length(d$X))); (c[,1]-c[,2])/(length(d$X)-1)}
gwas3 <- function(d,wv){N<-length(d$X);G<-d$G;M<-d$M; cnt<-as.vector(table(d$fam))
 cX<-as.numeric(d$X-tapply(d$X,d$fam,mean)[d$fam]); cW<-as.numeric(wv-tapply(wv,d$fam,mean)[d$fam])
 fmG<-rowsum(G,d$fam)/cnt; cG<-G-fmG[d$fam,]; df<-N-length(cnt)-1L
 den<-colSums(cG^2); numW<-colSums(cG*cW); numX<-colSums(cG*cX)
 rssW<-sum(cW^2)-numW^2/den; rssX<-sum(cX^2)-numX^2/den
 list(beta=numX/den,se_beta=sqrt((rssX/df)/den),delta=numW/den,se_delta=sqrt((rssW/df)/den))}

Fspecs <- list(c(2500,40),c(6000,30),c(12000,20)); ladder <- NULL
cat("[fig3] weak-IV ladder (mode=none)\n")
for (spec in Fspecs) {
  Fn3 <- as.integer(spec[1]); R3 <- as.integer(spec[2]); arch <- mkarch()
  res <- do.call(rbind, parallel::mclapply(1:R3, function(r) {
    d <- genSib3(Fn3, 6000L+r, "none", arch); o <- oracle3(d); wv <- winscore3(d)
    g <- gwas3(d, wv); iv <- mrwin_twosample_ivw(g$beta,g$delta,g$se_delta,se_gx=g$se_beta)
    ar <- mrwin_winmr_ar(g$beta,g$delta,g$se_delta,se_gx=g$se_beta)
    c(o=o, b=iv$gamma, wc=as.numeric(o>=iv$ci95[1]&&o<=iv$ci95[2]),
      arc=as.numeric(ar$contiguous&&o>=ar$ci[1]&&o<=ar$ci[2]), fs=mean(g$beta^2/g$se_beta^2))
  }, mc.cores=ncores))
  bias <- res[,"b"]-res[,"o"]
  ladder <- rbind(ladder, c(F=Fn3, fstat=mean(res[,"fs"]), bias=mean(bias), bias_se=sd(bias)/sqrt(R3),
                             wald_cov=mean(res[,"wc"]), ar_cov=mean(res[,"arc"])))
  cat(sprintf("  F=%-6d Fstat=%.2f  bias=%+.4f(SE%.4f)  Wald-cov=%.3f  AR-cov=%.3f\n",
    Fn3, ladder[nrow(ladder),"fstat"], ladder[nrow(ladder),"bias"], ladder[nrow(ladder),"bias_se"],
    ladder[nrow(ladder),"wald_cov"], ladder[nrow(ladder),"ar_cov"]))
}
ladder <- as.data.frame(ladder)
pdf_open("fig-weakiv-ladder.pdf", w = 6.8, h = 7)
layout(matrix(1:2,2,1)); par(mar=c(2.2,4.8,3,2.0))
xpad <- range(ladder$fstat) * c(0.85, 1.18)
ylo <- min(ladder$bias-2*ladder$bias_se) - 0.022
plot(ladder$fstat, ladder$bias, type="b", pch=19, col="#2A6FB0", log="x", xlim=xpad,
  xlab="", ylab="IVW point bias (vs. oracle)",
  ylim=c(ylo, max(0.012, ladder$bias+2*ladder$bias_se)),
  main="Two-sample within-family win-MR: weak-instrument ladder")
arrows(ladder$fstat, ladder$bias-ladder$bias_se, ladder$fstat, ladder$bias+ladder$bias_se, angle=90, code=3, length=0.05)
abline(h=0, col="grey50", lty=2)
text(ladder$fstat, rep(ylo+0.006,nrow(ladder)), labels=sprintf("F=%d",ladder$F), cex=0.72, col="grey30")
par(mar=c(4.6,4.8,1.6,2.0))
plot(ladder$fstat, ladder$ar_cov, type="b", pch=19, col="#2A6FB0", log="x", ylim=c(0.45,1.02), xlim=xpad,
  xlab="mean per-SNP instrument strength (F-statistic)", ylab="95% CI coverage of true gradient")
lines(ladder$fstat, ladder$wald_cov, type="b", pch=17, col="#C0392B")
abline(h=0.95, col="grey50", lty=2)
legend("bottomright", bty="n", pch=c(19,17), col=c("#2A6FB0","#C0392B"),
  legend=c("Anderson-Rubin (Paper 03)","IVW Wald"))
dev.off()
cat("[fig3] wrote fig-weakiv-ladder.pdf\n")
}
cat("paper04-figures.R done.\n")
