#' Simulated panel with a factor structure
#'
#' @description
#' `trop_simulate()` makes a balanced panel with unit effects, time effects, a
#' low-rank interactive term, and Gaussian noise. Each unit gets a score, which
#' is `selection` times its first factor loading plus random noise. The units
#' with the highest scores get the treatment.
#'
#' With a positive value of `selection`, a unit with a larger loading has a
#' higher probability of treatment. With a negative value, a unit with a smaller
#' loading has a higher probability. With zero, the loading has no effect on
#' the selection. Imputation with two-way fixed effects (TWFE imputation) does
#' not model this factor, so its estimate for one panel can have a large error.
#' TROP can reduce such errors, but some error can stay.
#'
#' @param n_units,n_periods The number of units and the number of periods.
#' @param n_treated The number of treated units.
#' @param n_post For a block design, the number of treated periods at the end of
#'   the panel. For a staggered design, the treatment starts in one of the last
#'   `n_post` periods.
#' @param effect The constant treatment effect on each treated cell.
#' @param rank The rank of the interactive term.
#' @param sigma The standard deviation of the noise.
#' @param selection The weight of the first factor loading in the selection of
#'   the treated units. With zero, the selection is random.
#' @param design Use `"block"` or `"staggered"`.
#' @param seed An optional seed. The function does not change the state of the
#'   random number generator in the R session.
#'
#' @return A data frame with the columns `unit`, `time`, `y`, `d`, and `y0`.
#'   The column `y0` contains the outcome without the treatment.
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
