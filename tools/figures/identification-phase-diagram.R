#!/usr/bin/env Rscript
# ============================================================================================
# THE IDENTIFICATION PHASE DIAGRAM -- a NEW figure type for Paper 03 (weak-IV-robust win-MR).
# Turns the "uniform validity across the strong->weak->degenerate continuum" theorem into a
# single phase portrait. Two INDEPENDENT axes of difficulty:
#   x = INSTRUMENT STRENGTH  (per-SNP exposure F-stat = (beta_GX/se_gx)^2; low = weak IV)
#   y = DEGENERACY           (the win U-statistic's first projection zeta1 = eps^2; eps->0 =
#                             degenerate, the moment's limit leaves Gaussian for a chi^2 mixture)
# Each cell: Monte-Carlo coverage of the true gradient by the SHIPPED basic Anderson-Rubin set
# vs the IVW WALD CI. Reuses the VALIDATED degeneracy kernel phi=eps*(Ai+Aj)+Ai*Aj from
# tools/validation-scripts/p3-degeneracy-supci-probe.R (closed-form O(n) projection + zeta1).
#
# READS AS THREE PHASES:
#   - strong & non-degenerate (upper right): Wald and AR both ~0.95.
#   - weak edge (left):  Wald UNDER-covers (weak-IV); AR HOLDS -> "AR required".
#   - degenerate edge (bottom): the plug-in se (4*zeta1/n) collapses while the true spread does
#     not, so BOTH Wald and basic AR under-cover -> "degeneracy-robust (sup-AR) required" (the
#     Andrews-Cheng sup-CI of P3 section 6 restores validity here; probed separately).
# Public simulator only; fixed seeds; base R + cairo. Usage: Rscript identification-phase-diagram.R [out.png]
# ============================================================================================
suppressMessages(library(mrwin)); ncores<-max(1L,min(4L,parallel::detectCores()))
args<-commandArgs(trailingOnly=TRUE); out<-if(length(args)>=1) args[1] else "identification-phase-diagram.png"
b0<-1.0; L<-6L; n<-150L; R<-600L
set.seed(7L); betaGX<-runif(L,0.40,0.60)                         # fixed instrument-exposure architecture
Fgrid<-c(1,5,25,125,625,3125); Egrid<-c(1.0,0.3,0.1,0.03,0.01,0.003)  # x = F-stat (log) ; y = degeneracy eps
# one SNP's win-odds estimate from the validated degenerate kernel: signal b0*betaGX + U-stat noise,
# with the plug-in se = sqrt(4*zeta1hat/n) (this is what UNDER-estimates the spread as eps->0).
usnp<-function(A,eps){nn<-length(A); SA<-sum(A); SS<-sum(A^2)
 Un<-(eps*(nn-1)*SA + (SA^2-SS)/2)/choose(nn,2)                  # U_n = mean_{i<j} phi(A_i,A_j)
 g<-(eps*((nn-2)*A+SA)+A*(SA-A))/(nn-1)                          # Hajek projection g_i (closed form)
 z1<-max(mean((g-mean(g))^2)-mean(A^2)^2/(nn-1),0)              # zeta1 estimate (>=0)
 c(Un=Un, se=2*sqrt(z1/nn))}
onecell<-function(Fstr,eps,seed){set.seed(seed)
 se_gx<-betaGX/sqrt(Fstr)
 hit<-replicate(R,{
   uu<-vapply(1:L,function(l) usnp(rnorm(n),eps),numeric(2))    # 2 x L : (Un, se) per SNP
   bhatGX<-betaGX+rnorm(L,0,se_gx)                              # NOISY exposure estimate -> weak-IV bias
   dhat<-b0*betaGX+uu[1,]; segy<-pmax(uu[2,],1e-6)              # outcome = true signal + U-stat noise
   iv<-tryCatch(mrwin_twosample_ivw(bhatGX,dhat,segy,se_gx=se_gx),error=function(e)NULL)
   ar<-tryCatch(mrwin_winmr_ar(bhatGX,dhat,segy,se_gx=se_gx),error=function(e)NULL)
   wc<-if(is.null(iv)) NA else as.numeric(b0>=iv$ci95[1] && b0<=iv$ci95[2])
   ac<-if(is.null(ar)||any(is.na(ar$ci))) NA else as.numeric(b0>=ar$ci[1] && b0<=ar$ci[2])  # hull contains truth
   bd<-if(is.null(ar)) NA else as.numeric(isTRUE(ar$bounded))
   c(wc,ac,bd)})
 rowMeans(hit,na.rm=TRUE)}
cat(sprintf("Identification phase diagram: %dx%d grid, L=%d SNPs, n=%d, R=%d reps/cell\n",length(Fgrid),length(Egrid),L,n,R))
ARm<-Wm<-Bm<-matrix(NA_real_,length(Fgrid),length(Egrid))
cells<-expand.grid(i=seq_along(Fgrid),j=seq_along(Egrid))
vals<-parallel::mclapply(1:nrow(cells),function(k){i<-cells$i[k];j<-cells$j[k]
  onecell(Fgrid[i],Egrid[j],1000L+k)},mc.cores=ncores)
for(k in 1:nrow(cells)){i<-cells$i[k];j<-cells$j[k];v<-vals[[k]];Wm[i,j]<-v[1];ARm[i,j]<-v[2];Bm[i,j]<-v[3]}
cat("AR coverage:\n"); print(round(t(ARm[,length(Egrid):1]),3))
cat("Wald coverage:\n"); print(round(t(Wm[,length(Egrid):1]),3))

# ---- render: two panels (AR | Wald), shared diverging palette centred at nominal 0.95 ----
# y orientation: non-degenerate (eps=1) at BOTTOM, degenerate (eps->0) at TOP (Egrid order).
pal<-colorRampPalette(c("#7F0000","#C13639","#E8896A","#F7F0EA","#9EC6E0","#2F6FB0"))(100)
zlim<-c(0.40,1.00); zmap<-function(z) pmin(pmax(z,zlim[1]),zlim[2])
xs<-seq_along(Fgrid); ys<-seq_along(Egrid)
png(out,width=1750,height=980,res=160,type="cairo")
layout(matrix(c(1,2,3),1,3),widths=c(1,1,0.20)); par(mar=c(4.8,5.0,3.6,1.0))
panel<-function(M,ttl){
 image(xs,ys,zmap(M),col=pal,zlim=zlim,axes=FALSE,xlab="instrument strength   (per-SNP F-statistic)",
       ylab=expression(paste("degeneracy   ",zeta[1]==epsilon^2,"   (up = ",epsilon%->%0,", degenerate)")),main=ttl,cex.main=1.05)
 axis(1,at=xs,labels=Fgrid); axis(2,at=ys,labels=Egrid,las=1); box(lwd=1.2)
 for(i in xs) for(j in ys) text(i,j,sprintf("%.2f",M[i,j]),cex=.74,col=ifelse(M[i,j]<0.72,"white","grey12"))}
# left: basic AR
panel(ARm,"Basic Anderson-Rubin set  -  coverage of true gradient")
rect(0.5,0.5,1.5,length(ys)+0.5,border="grey15",lwd=2,lty=1)
text(1.0,length(ys)+0.62,"weak IV:\nAR HOLDS",cex=.74,font=2,xpd=NA)
rect(2.95,3.95,5.05,4.75,col=rgb(1,1,1,.85),border="grey40")    # boxed crossover note
text(4.0,4.35,"near-degenerate crossover:\nbasic AR degrades  ->  sup-AR (P3 §6)",cex=.66,font=2,col="grey15")
arrows(5.0,4.2,5.92,3.05,length=.06,lwd=1.4,col="grey20")       # -> (F=3125, eps=0.1)
arrows(4.9,3.98,5.92,2.05,length=.06,lwd=1.4,col="grey20")      # -> (F=3125, eps=0.3)
# right: Wald
panel(Wm,"IVW Wald CI  -  coverage of true gradient")
rect(0.5,0.5,2.5,length(ys)+0.5,border="white",lwd=2,lty=1)
text(1.5,length(ys)+0.62,"weak IV:\nWald COLLAPSES",cex=.74,font=2,col="#7F0000",xpd=NA)
# colourbar
par(mar=c(4.8,0.6,3.6,3.6)); image(1,seq(zlim[1],zlim[2],length.out=100),matrix(seq(zlim[1],zlim[2],length.out=100),1),
      col=pal,zlim=zlim,axes=FALSE,xlab="",ylab=""); axis(4,at=c(.4,.6,.8,.9,.95,1),las=1); box()
mtext("coverage of the true gradient",4,line=2.4,cex=.85); abline(h=0.95,lwd=2.5)
dev.off()
cat(sprintf("wrote %s\n",out))
