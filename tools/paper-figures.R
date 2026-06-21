#!/usr/bin/env Rscript
# Reproduces the figures used in the mrwin methods write-up, entirely from the
# public package (simulation + estimator) -- no confidential data. Run from the
# repository root after installing mrwin:
#
#     Rscript tools/paper-figures.R
#
# By default it writes vector PDFs to manuscript/figures/ (the manuscript folder
# is gitignored / kept private until submission). Override the destination with
#     MRWIN_FIG_DIR=some/dir Rscript tools/paper-figures.R
# Seeds are fixed so every figure is bit-reproducible.

suppressMessages(library(mrwin))
outdir <- Sys.getenv("MRWIN_FIG_DIR", unset = file.path("manuscript", "figures"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pdf_open <- function(f, w = 5, h = 4) pdf(file.path(outdir, f), width = w, height = h)

## ---- Figure 1: wall-clock scaling, dense vs compiled fast (log-log) ----------
time_fit <- function(n, backend) {
  cfg <- mrwin_config(n_outcome = n, m_snps = 40L, seed = 909L)
  d <- mrwin_simulate(cfg, seed = 909L)
  ep <- mrwin_endpoint(d$time, d$status, colnames(d$time))
  gw <- mrwin_gwas(d$true_betas, rep(0.01, length(d$true_betas)))
  ctrl <- mrwin_controls(n_strata = 5L, bootstrap = 2L, seed = 1L,
                         backend = backend, run_sdpd = FALSE)
  as.numeric(system.time(
    mrwin(endpoint = ep, genotype = d$G, exposure = d$X, gwas = gw, controls = ctrl)
  )[["elapsed"]])
}
n_fast  <- c(1000, 2000, 4000, 8000, 16000)
n_dense <- c(1000, 2000, 4000)
t_fast  <- vapply(n_fast,  time_fit, numeric(1), backend = "fast")
t_dense <- vapply(n_dense, time_fit, numeric(1), backend = "dense")
sl_fast  <- coef(lm(log(t_fast)  ~ log(n_fast)))[2]
sl_dense <- coef(lm(log(t_dense) ~ log(n_dense)))[2]
pdf_open("fig-scaling.pdf")
plot(n_dense, t_dense, log = "xy", type = "b", pch = 19, col = "firebrick",
     xlim = range(c(n_fast, n_dense)), ylim = range(c(t_fast, t_dense)),
     xlab = "sample size N", ylab = "wall-clock (s)",
     main = "Kernel scaling (K = 3)")
lines(n_fast, t_fast, type = "b", pch = 17, col = "steelblue")
legend("topleft", bty = "n", pch = c(19, 17), col = c("firebrick", "steelblue"),
       legend = c(sprintf("dense  (slope %.2f)", sl_dense),
                  sprintf("fast   (slope %.2f)", sl_fast)))
dev.off()
cat(sprintf("[fig1] slopes  dense=%.2f  fast=%.2f\n", sl_dense, sl_fast))

## ---- Figure 2: delta-method SE is driven by the smallest exposure gap --------
nullcfg <- function(s) mrwin_config(n_outcome = 1500L, m_snps = 40L,
                                    alpha_x = c(0,0,0), gamma_direct = c(0,0,0), seed = s)
M <- 60L
se <- minDX <- numeric(M); unb <- logical(M)
for (i in seq_len(M)) {
  d <- mrwin_simulate(nullcfg(5000L + i), seed = 5000L + i)
  e <- mrwin_estimate(d$time, d$status, d$G, d$X, d$true_betas, n_strata = 5L)
  a <- mrwin_analytic_inference(e, d$X)
  se[i] <- a$se_delta_gls
  mX <- tapply(d$X, e$strata, mean); mX <- mX[order(as.integer(names(mX)))]
  minDX[i] <- min(abs(diff(mX)))
  unb[i] <- isTRUE(a$fieller_unbounded)
}
rr <- cor(log(se), log(minDX))
pdf_open("fig-denominator.pdf")
plot(minDX, se, log = "xy", pch = ifelse(unb, 1, 19),
     col = ifelse(unb, "darkorange", "grey30"),
     xlab = expression(min[d] ~ "|" * Delta * hat(X)[d] * "|  (smallest adjacent exposure gap)"),
     ylab = expression("delta-method  se(" * delta[GLS] * ")"),
     main = "Why the delta-method interval fails")
abline(lm(log(se) ~ log(minDX)), col = "steelblue", lwd = 2)
legend("topright", bty = "n",
       legend = c(sprintf("cor(log se, log min|dX|) = %.2f", rr),
                  "open = Fieller reports unbounded"))
dev.off()
cat(sprintf("[fig2] cor(log se, log minDX)=%.3f  se q50=%.2f q90=%.2f\n",
            rr, quantile(se, .5), quantile(se, .9)))

## ---- Figure 3: doubly-ranked balances the instrument across strata -----------
d <- mrwin_simulate(mrwin_config(n_outcome = 4000L, m_snps = 40L, seed = 7L), seed = 7L)
score <- as.numeric(d$G %*% d$true_betas)
e <- mrwin_estimate(d$time, d$status, d$G, d$X, d$true_betas, n_strata = 5L)
s_dr <- mrwin_doubly_ranked_strata(score, d$X, n_strata = 5L)$strata
ms <- function(s, v) { m <- tapply(v, s, mean); m[order(as.integer(names(m)))] }
pdf_open("fig-doubly-ranked.pdf", w = 5.5)
plot(1:5, ms(e$strata, score), type = "b", pch = 19, col = "firebrick", ylim = c(-1.6, 1.6),
     xlab = "stratum (low -> high)", ylab = "stratum mean",
     main = "PRS-rank vs doubly-ranked strata")
lines(1:5, ms(e$strata, d$X), type = "b", pch = 1, col = "firebrick", lty = 2)
lines(1:5, ms(s_dr, score), type = "b", pch = 17, col = "steelblue")
lines(1:5, ms(s_dr, d$X),  type = "b", pch = 2, col = "steelblue", lty = 2)
abline(h = 0, col = "grey80")
legend("topleft", bty = "n", cex = 0.85,
       legend = c("PRS-rank: instrument", "PRS-rank: exposure",
                  "doubly-ranked: instrument (flat ~ 0)", "doubly-ranked: exposure"),
       pch = c(19, 1, 17, 2), lty = c(1, 2, 1, 2),
       col = c("firebrick", "firebrick", "steelblue", "steelblue"))
dev.off()
cat(sprintf("[fig3] instrument spread prs=%.3f doubly=%.3f\n",
            diff(range(ms(e$strata, score))), diff(range(ms(s_dr, score)))))

## ---- Figure 4: example adjacent-stratum ISG forest ---------------------------
dd <- mrwin_simulate(mrwin_config(n_outcome = 3000L, m_snps = 50L, seed = 3L), seed = 3L)
fit <- mrwin(endpoint = mrwin_endpoint(dd$time, dd$status, colnames(dd$time)),
             genotype = dd$G, exposure = dd$X,
             gwas = mrwin_gwas(dd$true_betas, rep(0.01, length(dd$true_betas))),
             controls = mrwin_controls(n_strata = 5L, bootstrap = 200L, seed = 3L,
                                       run_sdpd = FALSE))
pdf_open("fig-isg-forest.pdf")
plot(fit, type = "isg")
dev.off()
cat("[fig4] ISG forest written\n")

cat("figures written to ", normalizePath(outdir), "\n")
