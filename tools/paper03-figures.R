#!/usr/bin/env Rscript
# Reproduces the figures for the Paper 03 (identification-robust win-ratio MR)
# manuscript, entirely from the public package + simulator -- no confidential
# data. Mirrors the tools/paper-figures.R / tools/paper04-figures.R convention.
#
#     Rscript tools/paper03-figures.R
#
# Writes to papers/03-weak-iv-robust-winmr/figures/ by default (gitignored
# private manuscript dir); override with MRWIN_FIG_DIR=some/dir. Seeds fixed
# for bit-reproducibility.
#
# fig 1 <- the Identification Phase Diagram (previously a PNG-only prototype in
#          tools/figures/identification-phase-diagram.R), re-targeted to a
#          vector PDF and inserted into the manuscript as Figure 1.
suppressMessages(library(mrwin))
outdir <- Sys.getenv("MRWIN_FIG_DIR", unset = file.path("papers", "03-weak-iv-robust-winmr", "figures"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ncores <- max(1L, min(4L, parallel::detectCores()))

## ---- Figure 1: Identification Phase Diagram --------------------------------
b0 <- 1.0; L <- 6L; n <- 150L; R <- 600L
set.seed(7L); betaGX <- runif(L, 0.40, 0.60)
Fgrid <- c(1, 5, 25, 125, 625, 3125); Egrid <- c(1.0, 0.3, 0.1, 0.03, 0.01, 0.003)
usnp <- function(A, eps) {
  nn <- length(A); SA <- sum(A); SS <- sum(A^2)
  Un <- (eps * (nn - 1) * SA + (SA^2 - SS) / 2) / choose(nn, 2)
  g <- (eps * ((nn - 2) * A + SA) + A * (SA - A)) / (nn - 1)
  z1 <- max(mean((g - mean(g))^2) - mean(A^2)^2 / (nn - 1), 0)
  c(Un = Un, se = 2 * sqrt(z1 / nn))
}
onecell <- function(Fstr, eps, seed) {
  set.seed(seed); se_gx <- betaGX / sqrt(Fstr)
  hit <- replicate(R, {
    uu <- vapply(1:L, function(l) usnp(rnorm(n), eps), numeric(2))
    bhatGX <- betaGX + rnorm(L, 0, se_gx)
    dhat <- b0 * betaGX + uu[1, ]; segy <- pmax(uu[2, ], 1e-6)
    iv <- tryCatch(mrwin_twosample_ivw(bhatGX, dhat, segy, se_gx = se_gx), error = function(e) NULL)
    ar <- tryCatch(mrwin_winmr_ar(bhatGX, dhat, segy, se_gx = se_gx), error = function(e) NULL)
    wc <- if (is.null(iv)) NA else as.numeric(b0 >= iv$ci95[1] && b0 <= iv$ci95[2])
    ac <- if (is.null(ar) || any(is.na(ar$ci))) NA else as.numeric(b0 >= ar$ci[1] && b0 <= ar$ci[2])
    c(wc, ac)
  })
  rowMeans(hit, na.rm = TRUE)
}
cat(sprintf("[fig1] Identification phase diagram: %dx%d grid, L=%d SNPs, n=%d, R=%d reps/cell\n",
            length(Fgrid), length(Egrid), L, n, R))
ARm <- Wm <- matrix(NA_real_, length(Fgrid), length(Egrid))
cells <- expand.grid(i = seq_along(Fgrid), j = seq_along(Egrid))
vals <- parallel::mclapply(1:nrow(cells), function(k) {
  i <- cells$i[k]; j <- cells$j[k]; onecell(Fgrid[i], Egrid[j], 1000L + k)
}, mc.cores = ncores)
for (k in 1:nrow(cells)) { i <- cells$i[k]; j <- cells$j[k]; v <- vals[[k]]; Wm[i, j] <- v[1]; ARm[i, j] <- v[2] }
cat("[fig1] AR coverage:\n"); print(round(t(ARm[, length(Egrid):1]), 3))
cat("[fig1] Wald coverage:\n"); print(round(t(Wm[, length(Egrid):1]), 3))

pal <- colorRampPalette(c("#7F0000", "#C13639", "#E8896A", "#F7F0EA", "#9EC6E0", "#2F6FB0"))(100)
zlim <- c(0.40, 1.00); zmap <- function(z) pmin(pmax(z, zlim[1]), zlim[2])
xs <- seq_along(Fgrid); ys <- seq_along(Egrid)
pdf(file.path(outdir, "fig-identification-phase-diagram.pdf"), width = 11.5, height = 6.2)
layout(matrix(c(1, 2, 3), 1, 3), widths = c(1, 1, 0.20)); par(mar = c(4.8, 5.0, 3.6, 1.0))
panel <- function(M, ttl) {
  image(xs, ys, zmap(M), col = pal, zlim = zlim, axes = FALSE,
        xlab = "instrument strength   (per-SNP F-statistic)",
        ylab = expression(paste("degeneracy   ", zeta[1] == epsilon^2, "   (up = ", epsilon %->% 0, ", degenerate)")),
        main = ttl, cex.main = 1.05)
  axis(1, at = xs, labels = Fgrid); axis(2, at = ys, labels = Egrid, las = 1); box(lwd = 1.2)
  for (i in xs) for (j in ys) text(i, j, sprintf("%.2f", M[i, j]), cex = .74, col = ifelse(M[i, j] < 0.72, "white", "grey12"))
}
panel(ARm, "Basic Anderson-Rubin set  -  coverage of true gradient")
rect(0.5, 0.5, 1.5, length(ys) + 0.5, border = "grey15", lwd = 2, lty = 1)
text(1.0, length(ys) + 0.62, "weak IV:\nAR HOLDS", cex = .74, font = 2, xpd = NA)
rect(2.95, 3.95, 5.05, 4.75, col = rgb(1, 1, 1, .85), border = "grey40")
text(4.0, 4.35, "near-degenerate crossover:\nbasic AR degrades  ->  sup-CI (Sec. 6)", cex = .66, font = 2, col = "grey15")
arrows(5.0, 4.2, 5.92, 3.05, length = .06, lwd = 1.4, col = "grey20")
arrows(4.9, 3.98, 5.92, 2.05, length = .06, lwd = 1.4, col = "grey20")
panel(Wm, "IVW Wald CI  -  coverage of true gradient")
rect(0.5, 0.5, 2.5, length(ys) + 0.5, border = "white", lwd = 2, lty = 1)
text(1.5, length(ys) + 0.62, "weak IV:\nWald COLLAPSES", cex = .74, font = 2, col = "#7F0000", xpd = NA)
par(mar = c(4.8, 0.6, 3.6, 3.6))
image(1, seq(zlim[1], zlim[2], length.out = 100), matrix(seq(zlim[1], zlim[2], length.out = 100), 1),
      col = pal, zlim = zlim, axes = FALSE, xlab = "", ylab = ""); axis(4, at = c(.4, .6, .8, .9, .95, 1), las = 1); box()
mtext("coverage of the true gradient", 4, line = 2.4, cex = .85); abline(h = 0.95, lwd = 2.5)
dev.off()
cat("[fig1] wrote fig-identification-phase-diagram.pdf\n")
cat("paper03-figures.R done.\n")
