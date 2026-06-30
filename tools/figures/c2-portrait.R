#!/usr/bin/env Rscript
# ============================================================================================
# THE CONCORDANCE-CONTRAST (C^2) PORTRAIT -- a NEW figure type for the mrwin programme.
# One glyph for "an instrument identifying a pairwise WIN on a HIERARCHICAL outcome." It is the
# reduced form of a U-statistic-valued IV, with three overlays no MR/win figure offers:
#   (1) identification SPINE + AR FAN   -- slope through origin = causal win gradient; the fan
#       of admissible slopes = the Anderson-Rubin confidence set mapped into this plane, so its
#       angular WIDTH reads off instrument strength (pencil=strong, wedge=weak, all-angles=degen).
#   (2) DE-CONFOUNDING overlay          -- within-family (de-confounded) vs naive between-family
#       reduced forms; the gap IS the dynastic/AM confounding the design removes.
#   (3) hierarchical RESOLUTION ribbon  -- share of decided wins resolved at tier1(death)/
#       tier2(hospitalization)/tier3(recurrence) across the instrument-contrast axis: WHICH
#       endpoint drives the genetic win gradient.
# Public package + simulator only; fixed seed; base R + cairo. Usage: Rscript c2-portrait.R [out.png]
# ============================================================================================
suppressMessages(library(mrwin)); cfg<-mrwin_config()
args<-commandArgs(trailingOnly=TRUE); out<-if(length(args)>=1) args[1] else "c2-portrait.png"
outc<-function(X,u,w,ex,cn,uf){F<-length(X);et<-matrix(NA_real_,F,3);lw<-log(w)
 for(j in 1:3){lp<-cfg$alpha_x[j]*X+cfg$nu_u[j]*u+ex+lw; et[,j]<-(-log(uf[,j])/(cfg$baseline_haz[j]*exp(lp)))^(1/cfg$shape_weibull)}
 td<-et[,1];th<-et[,2];tr<-et[,3]
 list(time=cbind(pmin(td,cn),pmin(th,cn,td),pmin(tr,cn,td)),status=cbind(as.integer(td<=cn),as.integer(th<=cn&th<=td),as.integer(tr<=cn&tr<=td)))}
tr1<-function(G) ifelse(G==2,1L,ifelse(G==0,0L,rbinom(length(G),1,0.5)))
# pairwise win h in {-1,0,1} AND the deciding priority tier (0 = tie)
hpair_tier<-function(t1,s1,t2,s2){F<-nrow(t1);r<-integer(F);tier<-integer(F)
 for(k in 1:3){iw<-tier==0 & s2[,k]==1 & t1[,k]>t2[,k]; il<-tier==0 & !iw & s1[,k]==1 & t2[,k]>t1[,k]
  r[iw]<-1L; tier[iw]<-k; r[il]<- -1L; tier[il]<-k}; list(h=r,tier=tier)}
genSib<-function(Fn,seed,gdyn=1.4){set.seed(seed);M<-30;base<-runif(M,.15,.35);betas<-rnorm(M,0,.15)
 Gf<-matrix(rbinom(Fn*M,2,base[col(matrix(0,Fn,M))]),Fn,M); Gm<-matrix(rbinom(Fn*M,2,base[col(matrix(0,Fn,M))]),Fn,M)
 mp<-as.numeric(((Gf+Gm)/2)%*%betas); sib<-vector("list",2); o1<-NULL
 for(k in 1:2){Gs<-matrix(tr1(Gf),Fn,M)+matrix(tr1(Gm),Fn,M); Z<-as.numeric(Gs%*%betas)
  u<-rnorm(Fn); X<-cfg$alpha_s*Z+cfg$alpha_u*u+rnorm(Fn); ex<-gdyn*scale(mp)[,1]
  w<-rgamma(Fn,1/cfg$theta_f,scale=cfg$theta_f); cn<-pmin(rexp(Fn,cfg$censoring_rate),cfg$max_follow_up); uf<-matrix(runif(Fn*3),Fn,3)
  oc<-outc(X,u,w,ex,cn,uf); sib[[k]]<-list(Z=Z,X=X,time=oc$time,status=oc$status); if(k==1) o1<-list(u=u,w=w,ex=ex,cn=cn,uf=uf,X=X)}
 list(sib=sib,o1=o1)}
oracle<-function(d,eps=0.25){A<-d$o1;np<-length(A$X)
 Bp<-outc(A$X+eps,A$u,A$w,A$ex,A$cn,A$uf); Bm<-outc(A$X-eps,A$u,A$w,A$ex,A$cn,A$uf)
 v<-mrwin_fast_pair_win_loss(Bp$time,Bp$status,Bm$time,Bm$status)
 (as.numeric(v[["wins"]])-as.numeric(v[["losses"]]))/(as.numeric(np)^2)/(2*eps)}

Fn<-12000L; d<-genSib(Fn,2026L); s1<-d$sib[[1]]; s2<-d$sib[[2]]
# WITHIN-FAMILY pairs (the design pairs): contrast, win, exposure-contrast, deciding tier
ht<-hpair_tier(s1$time,s1$status,s2$time,s2$status)
dZ<-s1$Z-s2$Z; dX<-s1$X-s2$X; h<-ht$h; tier<-ht$tier
# NAIVE between-family pairs: pair family f's sib1 with family (f+1)'s sib1 (different families)
sh<-c(2:Fn,1L); htn<-hpair_tier(s1$time,s1$status,s1$time[sh,],s1$status[sh,]); dZn<-s1$Z-s1$Z[sh]; hn<-htn$h
# reduced-form slopes + clustered-AR set (within-family) + oracle
sX<-sum(dZ*dX)/sum(dZ^2); sH<-sum(dZ*h)/sum(dZ^2); bhat<-sH/sX               # spine slope sH; gradient bhat
ARv<-function(b){m<-(h-b*dX)*dZ; sum(m)^2/sum(m^2)}
grid<-seq(bhat-6/sqrt(sum((dX*dZ)^2)/sum(((h-bhat*dX)*dZ)^2)+1e-9), bhat+6/sqrt(sum((dX*dZ)^2)/sum(((h-bhat*dX)*dZ)^2)+1e-9), length.out=4000)
acc<-grid[vapply(grid,ARv,1)<=3.841]; blo<-min(acc); bhi<-max(acc)           # AR confidence set on the gradient
borc<-oracle(d)
# bin the reduced forms over a common contrast range
brk<-seq(-quantile(abs(dZn),.97),quantile(abs(dZn),.97),length.out=19); mid<-(brk[-1]+brk[-length(brk)])/2
binmean<-function(x,y){i<-findInterval(x,brk); m<-tapply(y,factor(i,levels=1:(length(brk)-1)),mean); n<-tapply(y,factor(i,levels=1:(length(brk)-1)),length); list(m=as.numeric(m),n=as.numeric(n))}
wf<-binmean(dZ,h); nv<-binmean(dZn,hn)
# tier shares across contrast bins (within-family decided pairs)
dec<-tier>0; ti<-findInterval(dZ,brk)
tshare<-sapply(1:(length(brk)-1),function(b){idx<-dec&ti==b; if(sum(idx)<5) return(c(NA,NA,NA)); tabulate(tier[idx],3)/sum(idx)})

png(out,width=1500,height=1500,res=160,type="cairo")
layout(matrix(c(1,2),2,1),heights=c(2.3,1)); par(mar=c(4.2,4.6,3.2,1.4),xpd=FALSE)
xr<-range(brk)
plot(NA,xlim=xr,ylim=c(-.6,.6),xlab=expression(paste("instrument contrast  ",Delta,"Z (genetic)")),
     ylab=expression(paste("reduced-form win  ",hat(E),"[ h | ",Delta,"Z ]")),
     main=expression("Concordance-Contrast ("*C^2*") portrait  -  within-family win-MR (dynastic confounding)"),cex.main=1.02)
abline(h=0,col="grey80"); abline(v=0,col="grey80")
# (1b) AR fan = admissible reduced-form slopes (gradient set x exposure slope), shaded wedge
xx<-c(xr[1],xr[2]); polygon(c(xx,rev(xx)),c(blo*sX*xx,rev(bhi*sX*xx)),col=rgb(.20,.45,.80,.18),border=NA)
# oracle slope (dashed) and the identification spine (within-family point estimate)
lines(xx,borc*sX*xx,lty=2,lwd=2,col="black")
lines(xx,sH*xx,lwd=2.4,col=rgb(.13,.35,.75))
# (2) de-confounding: within-family (blue) vs naive between-family (red) binned reduced forms
pts<-function(b,col){ok<-is.finite(b$m); points(mid[ok],b$m[ok],pch=19,col=col,cex=pmin(2.4,.5+b$n[ok]/max(b$n,na.rm=TRUE)*2.2))}
pts(nv,rgb(.85,.30,.25,.9)); lines(mid,nv$m,col=rgb(.85,.30,.25,.7),lwd=1.5)
pts(wf,rgb(.13,.35,.75,.95))
legend("topleft",bty="n",cex=.92,
  legend=c(sprintf("within-family spine  (gradient %.2f)",bhat),
           sprintf("AR 95%% fan  [%.2f, %.2f]  -> %s instrument",blo,bhi,ifelse(bhi-blo>2*abs(bhat),"WEAK","strong")),
           sprintf("oracle gradient  %.2f  (dashed)",borc),
           "naive between-family (confounded)"),
  pch=c(NA,15,NA,19),lty=c(1,NA,2,NA),lwd=c(2.4,NA,2,NA),
  col=c(rgb(.13,.35,.75),rgb(.20,.45,.80,.5),"black",rgb(.85,.30,.25)))
# (3) hierarchical resolution ribbon
par(mar=c(4.4,4.6,1.2,1.4))
cols<-c(rgb(.20,.20,.25),rgb(.45,.55,.70),rgb(.75,.80,.90))
plot(NA,xlim=xr,ylim=c(0,1),xlab=expression(paste("instrument contrast  ",Delta,"Z")),ylab="share of decided wins",
     main="",cex.lab=1)
cum<-rep(0,ncol(tshare))
for(k in 1:3){top<-cum+ifelse(is.na(tshare[k,]),0,tshare[k,]); ok<-is.finite(top)
  polygon(c(mid[ok],rev(mid[ok])),c(cum[ok],rev(top[ok])),col=cols[k],border=NA); cum<-top}
legend("top",horiz=TRUE,bty="n",cex=.9,fill=cols,legend=c("tier 1: death","tier 2: hospitalization","tier 3: recurrence"))
dev.off()
cat(sprintf("wrote %s | within-fam gradient %.3f, AR [%.2f,%.2f], oracle %.3f, naive slope %.3f\n",
  out,bhat,blo,bhi,borc,sum(dZn*hn)/sum(dZn^2)/sX))
