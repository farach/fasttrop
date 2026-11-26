#' Make bootstrap rows
#'
#' `trop_boot_rows()` makes the matrix of unit indices for the unit bootstrap
#' of Algorithm 3. Each row of the matrix is one bootstrap draw. You can give
#' the matrix to the argument `boot_rows` of [trop()] and [trop_matrix()]. When
#' two or more fits use the same matrix, the draws have the same units. Then
#' the fits have a joint bootstrap distribution. Use this for a test of many
#' outcomes, for example with [trop_stepdown()].
#'
#' A unit is a treated unit if its row of `W` has at least one value of 1. The
#' function draws the treated units and the other units separately. With
#' `strata`, it draws the units of each stratum separately. Thus, each draw has
#' the same number of treated units and of other units in each stratum as the
#' data.
#'
#' A fit with fixed tuning parameters makes the same draws from the same
#' `seed`. Cross-validation with `cv_cells` uses random numbers before the
#' bootstrap. Then the draws are different.
#'
#' A fit uses a treated unit only if the unit has a treated cell with an
#' outcome. If all treated cells of a unit have an `NA` outcome, the fit
#' draws the unit with the other units. Use `W` with such units set to 0 to get
#' the same draws.
#'
#' @param W A matrix of treatment indicators with the values 0 and 1. Units are
#'   in rows.
#' @param n_boot The number of bootstrap draws.
#' @param strata An optional vector with one value for each unit. The values
#'   must not be `NA`.
#' @param resample Use `"stratified"` or `"controls"`. With `"controls"`, the
#'   treated units do not change in a draw. Refer to the section Inference of
#'   [trop()].
#' @param seed An optional integer seed. The function does not change the state
#'   of the random number generator in the R session.
#'
#' @return An integer matrix with `n_boot` rows and one column for each unit.
#'   Each value is a row number of `W`.
#' @seealso [trop_effects()], [trop_stepdown()]
#' @examples
#' W <- matrix(0L, 8, 5)
#' W[1:3, 4:5] <- 1L
#' trop_boot_rows(W, n_boot = 3, seed = 1)
#' trop_boot_rows(W, n_boot = 3, strata = rep(c("a", "b"), 4), seed = 1)
#' @export
trop_boot_rows <- function(W, n_boot = 200L, strata = NULL,
                           resample = c("stratified", "controls"), seed = NULL) {
  resample <- match.arg(resample)
  if (!is.matrix(W) || anyNA(W) || !all(W %in% c(0, 1))) {
    stop("`W` must be a matrix with only 0 and 1.", call. = FALSE)
  }
  check_count(n_boot, "n_boot")
  N <- nrow(W)
  if (!is.null(strata)) strata <- check_strata(strata, N)
  treated <- which(rowSums(W == 1) > 0L)
  never <- setdiff(seq_len(N), treated)
  rng <- seed_scope(seed)
  on.exit(rng$restore(), add = TRUE)
  boot_rows_matrix(treated, never, as.integer(n_boot), resample, strata)
}

# Draw units with replacement. The result has one row for each draw and the
# row numbers start at 1. Without strata, the random number stream is the
# stream of the unit bootstrap in earlier versions.
boot_rows_matrix <- function(treated, never, n_boot, resample = "stratified", strata = NULL) {
  N <- length(treated) + length(never)
  split_by <- function(idx) {
    if (is.null(strata) || length(idx) == 0L) return(list(idx))
    unname(split(idx, factor(strata[idx])))
  }
  draw <- function(groups) {
    unlist(lapply(groups, function(g) g[sample.int(length(g), length(g), replace = TRUE)]),
           use.names = FALSE)
  }
  tg <- split_by(treated)
  ng <- split_by(never)
  rows <- matrix(0L, n_boot, N)
  for (b in seq_len(n_boot)) {
    a <- if (resample == "stratified") draw(tg) else treated
    rows[b, ] <- c(a, draw(ng))
  }
  rows
}

check_strata <- function(strata, N) {
  if (!is.atomic(strata) || length(strata) != N || anyNA(strata)) {
    stop("`strata` must have one value for each unit, and no `NA` values.", call. = FALSE)
  }
  strata
}

check_boot_rows <- function(rows, N) {
  ok <- is.matrix(rows) && is.numeric(rows) && nrow(rows) >= 2L && ncol(rows) == N &&
    !anyNA(rows) && all(rows == round(rows)) && all(rows >= 1) && all(rows <= N)
  if (!ok) {
    stop("`boot_rows` must be a matrix with at least 2 rows and one column for each unit. ",
         "The values must be row numbers of the panel, from 1 to the number of units. ",
         "Use `trop_boot_rows()` to make such a matrix.", call. = FALSE)
  }
  storage.mode(rows) <- "integer"
  rows
}
