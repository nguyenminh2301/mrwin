# Two-sample / summary-data win-ratio MR: per-SNP win-odds estimator + SE, the
# win-odds GWAS, and the IVW pooling.

# brute-force per-subject win/loss counts of A vs B (reference for the C++ kernel)
.brute_counts <- function(At, As, Bt, Bs) {
  na <- nrow(At); nb <- nrow(Bt); K <- ncol(At)
  out <- matrix(0, na, 2)
  for (i in seq_len(na)) for (j in seq_len(nb)) {
    res <- 0L
    for (p in seq_len(K)) {
      if (Bs[j, p] == 1L && At[i, p] > Bt[j, p]) { res <- 1L; break }
      if (As[i, p] == 1L && Bt[j, p] > At[i, p]) { res <- -1L; break }
    }
    if (res == 1L) out[i, 1] <- out[i, 1] + 1 else if (res == -1L) out[i, 2] <- out[i, 2] + 1
  }
  out
}

test_that("mrwin_subject_win_loss_cpp matches brute-force per-subject counts (K=3, ties/censoring)", {
  cfg <- mrwin_config(n_outcome = 200, m_snps = 4, seed = 11)
  dat <- mrwin_simulate(cfg, seed = 11)
  d <- dat$G[, 1]
  hi <- which(d > stats::median(d)); lo <- which(d <= stats::median(d))
  At <- dat$time[hi, , drop = FALSE]; As <- dat$status[hi, , drop = FALSE]
  Bt <- dat$time[lo, , drop = FALSE]; Bs <- dat$status[lo, , drop = FALSE]
  f <- getFromNamespace("mrwin_subject_win_loss_cpp", "mrwin")
  cpp <- f(At, As, Bt, Bs, rep(1, nrow(Bt)))
  br <- .brute_counts(At, As, Bt, Bs)
  expect_equal(cpp[, 1], br[, 1])
  expect_equal(cpp[, 2], br[, 2])
  # totals reconcile with the production kernel
  wl <- mrwin_fast_pair_win_loss(At, As, Bt, Bs)
  expect_equal(sum(cpp[, 1]), as.numeric(wl[["wins"]]))
  expect_equal(sum(cpp[, 2]), as.numeric(wl[["losses"]]))
})

test_that("mrwin_win_snp SE equals the dense influence-function variance (exact, no subsampling)", {
  cfg <- mrwin_config(n_outcome = 1200, m_snps = 6, seed = 22)
  dat <- mrwin_simulate(cfg, seed = 22)
  d <- dat$G[, 2]
  cut <- stats::median(d); hi <- which(d > cut); lo <- which(d <= cut)
  # dense influence-function SE of log-theta for this contrast
  tt <- rbind(dat$time[hi, , drop = FALSE], dat$time[lo, , drop = FALSE])
  ss <- rbind(dat$status[hi, , drop = FALSE], dat$status[lo, , drop = FALSE])
  Kr <- mrwin_kernel(tt, ss)
  iH <- seq_along(hi); iL <- length(hi) + seq_along(lo)
  win <- (Kr[iH, iL, drop = FALSE] == 1) + 0; loss <- (Kr[iH, iL, drop = FALSE] == -1) + 0
  W0 <- sum(win); L0 <- sum(loss)
  cH <- rowSums(win) / W0 - rowSums(loss) / L0
  cL <- colSums(win) / W0 - colSums(loss) / L0
  se_dense <- sqrt(sum(cH^2) + sum(cL^2))

  # max_subjects above group sizes -> no subsampling -> must match exactly
  r <- mrwin_win_snp(dat$time, dat$status, d, max_subjects = 10000L)
  expect_equal(r$se_log_theta, se_dense, tolerance = 1e-6)
  expect_equal(r$log_theta, log(W0 / L0), tolerance = 1e-9)
  expect_equal(r$delta, r$log_theta / r$d_diff, tolerance = 1e-12)
})

test_that("mrwin_win_gwas returns one row per SNP with finite estimates", {
  cfg <- mrwin_config(n_outcome = 800, m_snps = 8, seed = 33)
  dat <- mrwin_simulate(cfg, seed = 33)
  wg <- mrwin_win_gwas(dat$time, dat$status, dat$G, max_subjects = 500L, seed = 1)
  expect_equal(nrow(wg), ncol(dat$G))
  expect_true(all(is.finite(wg$delta)))
  expect_true(all(wg$se > 0))
})

test_that("mrwin_twosample_ivw reduces to the Wald ratio for one SNP and weights by precision", {
  # single SNP: gamma = delta/beta, se = se_gy/|beta|
  one <- mrwin_twosample_ivw(beta_gx = 0.4, delta_gy = 0.12, se_gy = 0.05)
  expect_equal(one$gamma, 0.12 / 0.4, tolerance = 1e-12)
  expect_equal(one$se, 0.05 / 0.4, tolerance = 1e-12)
  # two SNPs, same ratio -> gamma unchanged, tighter SE; Q ~ 0
  two <- mrwin_twosample_ivw(c(0.4, 0.4), c(0.12, 0.12), c(0.05, 0.05))
  expect_equal(two$gamma, 0.3, tolerance = 1e-12)
  expect_lt(two$se, one$se)
  expect_lt(two$Q, 1e-8)
})
