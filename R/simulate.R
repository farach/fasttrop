#' Simulate a panel with a factor structure
#'
#' Generates a balanced panel with unit and time effects, a low-rank
#' interactive component, and Gaussian noise. Treatment selection can depend
#' on the factor loadings, which biases two-way fixed effects and is the
#' setting TROP targets.
#'
#' @param n_units,n_periods Panel dimensions.
#' @param n_treated Number of treated units.
#' @param n_post Number of treated periods at the end of the panel (block
#'   design) or the latest adoption window (staggered design).
#' @param effect Constant treatment effect on treated cells.
#' @param rank Rank of the interactive component.
#' @param sigma Noise standard deviation.
#' @param selection Strength of selection on the first factor loading.
#' @param design `"block"` or `"staggered"`.
#' @param seed Optional seed; the global random state is restored afterwards.
#'
#' @return A data frame with columns `unit`, `time`, `y`, `d`, and `y0`
#'   (the untreated potential outcome).
#' @examples
#' head(trop_simulate(seed = 1))
#' @export
trop_simulate <- function(n_units = 40L, n_periods = 30L, n_treated = 8L, n_post = 5L,
                          effect = 1, rank = 2L, sigma = 0.5, selection = 1,
                          design = c("block", "staggered"), seed = NULL) {
  design <- match.arg(design)
  check_count(n_units, "n_units")
  check_count(n_periods, "n_periods")
  check_count(n_treated, "n_treated")
  check_count(n_post, "n_post")
  if (n_treated >= n_units) stop("`n_treated` must be smaller than `n_units`.", call. = FALSE)
  if (n_post >= n_periods) stop("`n_post` must be smaller than `n_periods`.", call. = FALSE)
  rng <- seed_scope(seed)
  on.exit(rng$restore(), add = TRUE)
  N <- as.integer(n_units)
  TT <- as.integer(n_periods)
  a <- stats::rnorm(N)
  b <- cumsum(stats::rnorm(TT, 0, 0.3))
  G <- matrix(stats::rnorm(N * rank), N, rank)
  F <- matrix(stats::rnorm(TT * rank), TT, rank)
  Y0 <- outer(a, rep(1, TT)) + outer(rep(1, N), b) + G %*% t(F) +
    matrix(stats::rnorm(N * TT, 0, sigma), N, TT)
  score <- selection * G[, 1] + stats::rnorm(N)
  treated <- order(score, decreasing = TRUE)[seq_len(n_treated)]
  W <- matrix(0L, N, TT)
  if (design == "block") {
    W[treated, (TT - n_post + 1L):TT] <- 1L
  } else {
    starts <- sample((TT - n_post + 1L):TT, n_treated, replace = TRUE)
    for (m in seq_along(treated)) W[treated[m], starts[m]:TT] <- 1L
  }
  Y <- Y0 + effect * W
  data.frame(unit = rep(seq_len(N), times = TT), time = rep(seq_len(TT), each = N),
             y = as.vector(Y), d = as.vector(W), y0 = as.vector(Y0))
}
