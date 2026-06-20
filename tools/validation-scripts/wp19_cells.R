suppressMessages(library(mrwin))
vc <- getFromNamespace(".mrwin_validate_calibration", "mrwin")
out <- "/tmp/wp19_cells.out"
cat("WP19 additional cells (strong instrument N=3000, m=100, D=5)\n", file = out)
log_cell <- function(tag, r) {
  line <- sprintf("%-40s rate=%.3f  (M=%d valid=%d, MC_SE=%.3f, weak=%.2f)\n",
                  tag, r$rate, r$M, r$valid, sqrt(0.05 * 0.95 / max(r$valid, 1)), r$weak_frac)
  cat(line); cat(line, file = out, append = TRUE)
}

# Cell 1: documented gap -- bootstrap x IPTW (null type-I)
t1 <- system.time(r1 <- vc(M = 120L, N = 3000L, m_snps = 100L, n_strata = 5L,
  scenario = "null", inference = "bootstrap", sigma_beta = 0,
  adjustment = "ordinal_iptw", stratification = "prs_rank",
  bootstrap = 200L, seed0 = 2000L))
log_cell("null | bootstrap | IPTW | prs_rank", r1)
cat(sprintf("  [%.0fs]\n", t1[["elapsed"]]), file = out, append = TRUE)

# Cell 2: NEW M2 path -- doubly-ranked, bootstrap (null type-I)
t2 <- system.time(r2 <- vc(M = 150L, N = 3000L, m_snps = 100L, n_strata = 5L,
  scenario = "null", inference = "bootstrap", sigma_beta = 0,
  adjustment = "none", stratification = "doubly_ranked",
  bootstrap = 200L, seed0 = 3000L))
log_cell("null | bootstrap | none | doubly_ranked", r2)
cat(sprintf("  [%.0fs]\n", t2[["elapsed"]]), file = out, append = TRUE)

# Cell 3: NEW M2 path -- doubly-ranked, analytic (null type-I)
t3 <- system.time(r3 <- vc(M = 150L, N = 3000L, m_snps = 100L, n_strata = 5L,
  scenario = "null", inference = "analytic", sigma_beta = 0,
  adjustment = "none", stratification = "doubly_ranked",
  bootstrap = 200L, seed0 = 4000L))
log_cell("null | analytic | none | doubly_ranked", r3)
cat(sprintf("  [%.0fs]\n", t3[["elapsed"]]), file = out, append = TRUE)

cat("WP19_CELLS_DONE\n", file = out, append = TRUE)
cat("WP19_CELLS_DONE\n")
