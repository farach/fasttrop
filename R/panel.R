#' Make outcome and treatment matrices from long panel data
#'
#' @param data A data frame with one row for each unit and period.
#' @param outcome,treatment,unit,time Column names.
#' @return A list with `Y` (double, with `NA` for a cell without an outcome),
#'   `W` (integer, 0 or 1), and the sorted `units` and `times`.
#' @noRd
panel_matrices <- function(data, outcome, treatment, unit, time) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  for (nm in c("outcome", "treatment", "unit", "time")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L) {
      stop(sprintf("`%s` must be one column name.", nm), call. = FALSE)
    }
    if (!v %in% names(data)) stop(sprintf("Column `%s` is not in `data`.", v), call. = FALSE)
  }
  y <- data[[outcome]]
  d <- data[[treatment]]
  if (!is.numeric(y)) stop(sprintf("Outcome `%s` must be numeric.", outcome), call. = FALSE)
  if (is.logical(d)) d <- as.integer(d)
  if (!is.numeric(d) || anyNA(d) || !all(d %in% c(0, 1))) {
    stop(sprintf("Treatment `%s` must contain only 0 and 1, or TRUE and FALSE, and no `NA` values.",
                 treatment), call. = FALSE)
  }
  if (anyNA(data[[unit]]) || anyNA(data[[time]])) {
    stop("Unit and time identifiers must not be `NA`.", call. = FALSE)
  }
  units <- sort(unique(data[[unit]]))
  times <- sort(unique(data[[time]]))
  i <- match(data[[unit]], units)
  t <- match(data[[time]], times)
  if (anyDuplicated(i + length(units) * (t - 1L))) {
    stop("`data` has two or more rows for the same unit and period.", call. = FALSE)
  }
  N <- length(units)
  TT <- length(times)
  Y <- matrix(NA_real_, N, TT, dimnames = list(as.character(units), as.character(times)))
  W <- matrix(0L, N, TT, dimnames = dimnames(Y))
  idx <- cbind(i, t)
  Y[idx] <- as.double(y)
  W[idx] <- as.integer(d)
  list(Y = Y, W = W, units = units, times = times)
}

#' Find the assignment design
#'
#' In a block design, all treated units start in the same period and stay
#' treated to the end of the panel. In a staggered design, the treatment does
#' not stop, but the start periods are different. In a general design, the
#' treatment stops for one or more units.
#' @noRd
trop_design <- function(W, obs) {
  N <- nrow(W)
  TT <- ncol(W)
  treated <- which(rowSums(W == 1L & obs == 1L) > 0L)
  if (length(treated) == 0L) stop("The panel has no treated cells.", call. = FALSE)
  if (!any(W == 0L & obs == 1L)) stop("The panel has no untreated cells with an outcome.", call. = FALSE)
  first <- vapply(treated, function(i) which(W[i, ] == 1L)[1L], integer(1))
  absorbing <- all(vapply(seq_along(treated), function(m) {
    all(W[treated[m], first[m]:TT] == 1L)
  }, logical(1)))
  type <- if (!absorbing) {
    "general"
  } else if (length(unique(first)) == 1L) {
    "block"
  } else {
    "staggered"
  }
  k <- if (type == "block") TT - first[1L] + 1L else NA_integer_
  list(type = type, treated = treated, first = first, k = k,
       n_treated_units = length(treated), n_treated_cells = sum(W == 1L & obs == 1L))
}

as_panel_inputs <- function(Y, W) {
  if (!is.matrix(Y) || !is.numeric(Y)) stop("`Y` must be a numeric matrix with units in rows and periods in columns.", call. = FALSE)
  if (!is.matrix(W) || !all(dim(W) == dim(Y))) stop("`W` must be a matrix with the same dimensions as `Y`.", call. = FALSE)
  if (is.logical(W)) W <- W * 1L
  if (anyNA(W) || !all(W %in% c(0, 1))) stop("`W` must contain only 0 and 1.", call. = FALSE)
  if (nrow(Y) < 3L || ncol(Y) < 3L) stop("The panel must have at least 3 units and 3 periods.", call. = FALSE)
  obs <- matrix(as.integer(is.finite(Y)), nrow(Y), ncol(Y))
  Y0 <- Y
  storage.mode(Y0) <- "double"
  Wi <- matrix(as.integer(W), nrow(W), ncol(W))
  list(Y = Y0, W = Wi, obs = obs)
}
