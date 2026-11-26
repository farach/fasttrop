#' Scale the tuning grid to the data
#'
#' `trop_grid_scaled()` changes the values of a tuning grid so that they fit
#' the scale of the outcome. The default grid of [trop_grid()] uses fixed
#' values. These values fit an outcome with a standard deviation of about 1.
#' For an outcome with other units, the best value of a parameter can be
#' outside of the grid. Then cross-validation chooses the edge of the grid.
#'
#' The function changes the grid in two steps:
#' * The unit weights decrease with the root mean square gap between the
#'   outcome histories of two units. The function divides the values of `unit`
#'   by the median of this gap. The gap uses the untreated cells that both units
#'   have. For a panel with more than `max_units` units, the median uses
#'   `max_units` units that are the same distance apart in the order of the
#'   rows.
#' * The nuclear-norm penalty depends on the scale of the noise. The function
#'   multiplies the finite values of `nn` by a scale. The scale is the root
#'   mean square of the residuals of a two-way fixed effects fit to the
#'   untreated cells.
#'
#' The values of `time` and the value `Inf` do not change. The scales are in
#' the attribute `scales` of the result.
#'
#' With `grid = "auto"`, [trop()] and [trop_matrix()] use this function with
#' the default grid, with the same `demean` and the same `target`.
#'
#' @param Y A numeric matrix of outcomes with units in rows. `NA` shows a cell
#'   without an outcome.
#' @param W A matrix of treatment indicators with the values 0 and 1. It must
#'   have the same dimensions as `Y`.
#' @param grid The tuning grid from [trop_grid()] to scale.
#' @param demean Logical. If `TRUE`, the function subtracts the mean of the
#'   untreated outcomes of each unit before it calculates the scales. Refer to
#'   the argument `demean` of [trop()].
#' @param max_units The maximum number of units for the median gap.
#'
#' @return A list of class `trop_grid` with the attribute `scales`. The
#'   attribute is a list with the elements `unit` (the median gap) and `nn` (the
#'   residual scale).
#' @examples
#' sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5,
#'                      n_post = 3, effect = 1, seed = 1)
#' pm <- tapply(sim$y, list(sim$unit, sim$time), sum)
#' W <- tapply(sim$d, list(sim$unit, sim$time), sum)
#' g <- trop_grid_scaled(pm * 100, W)
#' attr(g, "scales")
#' @export
trop_grid_scaled <- function(Y, W, grid = trop_grid(), demean = FALSE, max_units = 500L) {
  if (!inherits(grid, "trop_grid")) grid <- do.call(trop_grid, as.list(grid))
  check_flag(demean, "demean")
  check_count(max_units, "max_units")
  p <- as_panel_inputs(Y, W)
  if (demean) p$Y <- demean_units(p$Y, p$W, p$obs)
  sc <- grid_scales(p$Y, p$W, p$obs, max_units)
  out <- trop_grid(time = grid$time, unit = grid$unit / sc$unit,
                   nn = ifelse(is.finite(grid$nn), grid$nn * sc$nn, Inf))
  attr(out, "scales") <- sc
  out
}

grid_scales <- function(Y, W, obs, max_units = 500L) {
  M <- (W == 0L & obs == 1L)
  N <- nrow(Y)
  keep <- if (N > max_units) unique(round(seq(1, N, length.out = max_units))) else seq_len(N)
  Mk <- M[keep, , drop = FALSE] * 1
  Y0 <- Y[keep, , drop = FALSE]
  Y0[Mk == 0] <- 0
  S <- tcrossprod(Y0^2, Mk) + tcrossprod(Mk, Y0^2) - 2 * tcrossprod(Y0)
  n <- tcrossprod(Mk)
  up <- upper.tri(S)
  ok <- up & n > 0
  gap <- sqrt(pmax(S[ok], 0) / n[ok])
  gap <- gap[is.finite(gap) & gap > 0]
  unit <- if (length(gap) > 0L) stats::median(gap) else 1
  nn <- twfe_residual_scale(Y, M)
  list(unit = unit, nn = nn)
}

# Root mean square of the residuals of a two-way fixed effects fit to the cells
# in the mask M. Alternating means give the fit.
twfe_residual_scale <- function(Y, M) {
  r <- Y
  r[!M] <- NA_real_
  for (k in seq_len(100L)) {
    ru <- rowMeans(r, na.rm = TRUE)
    ru[is.nan(ru)] <- 0
    r <- r - ru
    rt <- colMeans(r, na.rm = TRUE)
    rt[is.nan(rt)] <- 0
    r <- sweep(r, 2L, rt)
    if (max(abs(ru), abs(rt)) < 1e-10) break
  }
  s <- sqrt(mean(r^2, na.rm = TRUE))
  if (is.finite(s) && s > 0) s else 1
}

# Subtract the mean of the untreated outcomes of each unit. A unit effect
# absorbs a unit constant. Thus, the fit and the effects do not change. Only
# the distances between units change.
demean_units <- function(Y, W, obs) {
  M <- (W == 0L & obs == 1L)
  Y0 <- Y
  Y0[!M] <- 0
  n <- rowSums(M)
  mu <- ifelse(n > 0L, rowSums(Y0) / pmax(n, 1L), 0)
  Y - mu
}

# Find the tuning values at the edge of the grid.
grid_boundary <- function(lambda, grid) {
  tol <- function(a, b) is.finite(b) && abs(a - b) <= 1e-9 * max(1, abs(b))
  one <- function(name) {
    val <- lambda[[name]]
    vec <- grid[[name]]
    if (length(unique(vec)) < 2L || !is.finite(val)) return(NULL)
    if (tol(val, max(vec))) return(paste0("largest value of `", name, "` (", format(val, digits = 4L), ")"))
    if (val > 0 && tol(val, min(vec))) {
      return(paste0("smallest value of `", name, "` (", format(val, digits = 4L), ")"))
    }
    NULL
  }
  unlist(lapply(c("time", "unit", "nn"), one), use.names = FALSE)
}

package_info <- function() {
  d <- tryCatch(
    suppressWarnings(read.dcf(system.file("DESCRIPTION", package = "fasttrop"),
                              fields = c("Version", "RemoteSha", "GithubSHA1"))),
    error = function(e) NULL
  )
  if (is.null(d)) return(list(version = NA_character_, sha = NA_character_))
  sha <- d[1L, c("RemoteSha", "GithubSHA1")]
  sha <- sha[!is.na(sha)]
  list(version = unname(d[1L, "Version"]), sha = if (length(sha) > 0L) unname(sha[1L]) else NA_character_)
}

check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  invisible(x)
}
