# Generated from python.p1_engine_v5.kernel using fixed deterministic arrays.
# This fixture verifies WP4 point DS-CWR estimator parity independent of R RNG.
wp4_estimator_small <- list(
  n_strata = 4L,
  time = matrix(c(
    5, 7, 9,
    3, 8, 6,
    6, 4, 5,
    8, 6, 7,
    2, 9, 8,
    9, 3, 4,
    4, 5, 8,
    7, 2, 6
  ), ncol = 3, byrow = TRUE),
  status = matrix(c(
    0, 1, 1,
    1, 0, 1,
    0, 1, 0,
    1, 1, 1,
    1, 0, 0,
    0, 1, 1,
    1, 1, 0,
    0, 1, 1
  ), ncol = 3, byrow = TRUE),
  G = matrix(c(
    0, 2, 0,
    0, 1, 0,
    0, 0, 0,
    1, 1, 0,
    1, 0, 0,
    1, 0, 1,
    2, 0, 0,
    2, 0, 1
  ), ncol = 3, byrow = TRUE),
  X = c(0.0, 0.2, 0.5, 0.7, 1.2, 1.4, 2.0, 2.2),
  beta_hat = c(0.2, -0.1, 0.15),
  expected = list(
    score = c(-0.2, -0.1, 0.0, 0.1, 0.2, 0.35, 0.4, 0.55),
    strata = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L),
    wins = c(s2_vs_s1 = 2, s3_vs_s2 = 1, s4_vs_s3 = 2),
    losses = c(s2_vs_s1 = 2, s3_vs_s2 = 3, s4_vs_s3 = 2),
    total = c(s2_vs_s1 = 4, s3_vs_s2 = 4, s4_vs_s3 = 4),
    log_theta = c(
      s2_vs_s1 = 0.0,
      s3_vs_s2 = -1.0986122886681098,
      s4_vs_s3 = 0.0
    ),
    cwr = c(s2_vs_s1 = 1.0, s3_vs_s2 = 0.3333333333333333, s4_vs_s3 = 1.0),
    delta_x = c(
      s2_vs_s1 = 0.5,
      s3_vs_s2 = 0.6999999999999998,
      s4_vs_s3 = 0.8000000000000003
    ),
    delta_isg = c(
      s2_vs_s1 = 0.0,
      s3_vs_s2 = -1.5694461266687285,
      s4_vs_s3 = 0.0
    )
  ),
  sigma_isg = matrix(c(
    0.25, 0.03, 0.02,
    0.03, 0.20, 0.04,
    0.02, 0.04, 0.30
  ), ncol = 3, byrow = TRUE),
  gls_no_shrink = list(
    delta_gls = -0.6331099024778309,
    se_delta_gls = 0.31810212143319205,
    dscwr = 0.5309380652565682,
    q = 8.890390542837165,
    q_df = 2L
  )
)
