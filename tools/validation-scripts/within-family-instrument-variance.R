#!/usr/bin/env Rscript
# Why are within-family instruments weak? A mechanistic decomposition (Paper 04,
# Simulation study, "why these patterns arise"), verified by direct simulation
# rather than asserted from algebra alone.
#
# CLAIM 1 (variance is NOT smaller per genotyped subject). For a single unlinked,
# additive locus with allele frequency p, population dosage variance = 2p(1-p)
# (Hardy-Weinberg). Conditional on parental genotypes, two siblings' transmitted
# alleles at that locus are independent fair-coin draws (Mendel's first law), so
# the sib difference's conditional variance is 1/4 if exactly one parent is
# heterozygous (probability 2p(1-p) under random mating) and 0 otherwise;
# averaging over random-mating parental genotypes gives an UNCONDITIONAL sib-
# difference variance of exactly 2p(1-p) -- the SAME as the population dosage
# variance. Summed over unlinked loci: Var(dZ) = Var(Z) exactly.
#
# CLAIM 2 (the DESIGN, not the instrument, is weaker). Comparing dX to dZ uses F
# family contrasts instead of 2F individual observations (half the effective N),
# and differencing two siblings' independent idiosyncratic exposure noise DOUBLES
# the first-stage residual variance relative to a population regression of X on
# Z. Both effects reduce the first-stage F-statistic and compound multiplicatively
# (2x2=4x), even though Var(dZ)=Var(Z).
#
# RESULT (200 replications, F=2500 families = 2F=5000 individuals):
#   mean Var(Z) = 0.2494   mean Var(dZ) = 0.2494   ratio = 1.000  (confirms claim 1)
#   mean population first-stage F  (N=5000 individuals) = 120.98 (SE 2.68)
#   mean within-family first-stage F (F=2500 families)   =  30.36 (SE 0.95)
#   ratio of means = 3.99  (confirms claim 2: ~4x, i.e. 2x from halved N times 2x
#   from doubled residual noise)
#
# This decomposition is specific to this paper (not a quotation from the within-
# family MR literature, which documents the qualitative power cost without, to our
# knowledge, decomposing it this way). Public package + simulator only.
suppressMessages(library(mrwin)); cfg <- mrwin_config()
transmit <- function(G) ifelse(G==2,1,ifelse(G==0,0,rbinom(length(G),1,0.5)))
one_rep <- function(seed, Fn=2500L, M=30L) {
  set.seed(seed); base <- runif(M,.15,.35); betas <- rnorm(M,0,.15)
  Gf <- matrix(rbinom(Fn*M,2,base[col(matrix(0,Fn,M))]),Fn,M)
  Gm <- matrix(rbinom(Fn*M,2,base[col(matrix(0,Fn,M))]),Fn,M)
  sib <- function() t(apply(Gf,1,transmit)) + t(apply(Gm,1,transmit))
  Z1 <- as.numeric(sib() %*% betas); Z2 <- as.numeric(sib() %*% betas)
  u1 <- rnorm(Fn); u2 <- rnorm(Fn)
  X1 <- cfg$alpha_s*Z1 + cfg$alpha_u*u1 + rnorm(Fn)
  X2 <- cfg$alpha_s*Z2 + cfg$alpha_u*u2 + rnorm(Fn)
  Zpop <- c(Z1,Z2); Xpop <- c(X1,X2)
  dZ <- Z1-Z2; dX <- X1-X2
  c(varZ=var(Zpop), vardZ=var(dZ),
    fpop=unname(summary(lm(Xpop~Zpop))$fstatistic[1]),
    fwf=unname(summary(lm(dX~dZ))$fstatistic[1]))
}
Fn <- as.integer(Sys.getenv("IVAR_F","2500")); R0 <- as.integer(Sys.getenv("IVAR_R","200"))
ncores <- max(1L, min(4L, parallel::detectCores()))
cat(sprintf("Within-family instrument-variance decomposition: F=%d families, R=%d reps\n", Fn, R0))
res <- do.call(rbind, parallel::mclapply(1:R0, one_rep, Fn=Fn, mc.cores=ncores))
cat(sprintf("mean Var(Z)=%.4f  mean Var(dZ)=%.4f  ratio=%.3f\n",
  mean(res[,"varZ"]), mean(res[,"vardZ"]), mean(res[,"vardZ"])/mean(res[,"varZ"])))
cat(sprintf("mean population F (N=%d individuals) = %.2f (SE %.2f)\n",
  2*Fn, mean(res[,"fpop"]), sd(res[,"fpop"])/sqrt(R0)))
cat(sprintf("mean within-family F (F=%d families) = %.2f (SE %.2f)\n",
  Fn, mean(res[,"fwf"]), sd(res[,"fwf"])/sqrt(R0)))
cat(sprintf("ratio of means = %.2f\n", mean(res[,"fpop"])/mean(res[,"fwf"])))
cat("DONEIVAR\n")
