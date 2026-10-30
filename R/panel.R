#' Build outcome and treatment matrices from long panel data
#'
#' @param data A data frame with one row per unit-period.
#' @param outcome,treatment,unit,time Column names.
#' @return A list with `Y` (double, NA where unobserved), `W` (integer 0/1),
#'   `obs` (integer 0/1), and the sorted `units` and `times`.
#' @noRd
panel_matrices <- function(data, outcome, treatment, unit, time) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  for (nm in c("outcome", "treatment", "unit", "time")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L) {
      stop(sprintf("`%s` must be a single column name.", nm), call. = FALSE)
    }
    if (!v %in% names(data)) stop(sprintf("Column `%s` not found in `data`.", v), call. = FALSE)
  }
  y <- data[[outcome]]
  d <- data[[treatment]]
  if (!is.numeric(y)) stop(sprintf("Outcome `%s` must be numeric.", outcome), call. = FALSE)
  if (is.logical(d)) d <- as.integer(d)
  if (!is.numeric(d) || anyNA(d) || !all(d %in% c(0, 1))) {
    stop(sprintf("Treatment `%s` must be 0/1 (or logical) without missing values.", treatment),
         call. = FALSE)
  }
  if (anyNA(data[[unit]]) || anyNA(data[[time]])) {
    stop("Unit and time identifiers must not be missing.", call. = FALSE)
  }
  units <- sort(unique(data[[unit]]))
  times <- sort(unique(data[[time]]))
  i <- match(data[[unit]], units)
  t <- match(data[[time]], times)
  if (anyDuplicated(i + length(units) * (t - 1L))) {
    stop("`data` has duplicated unit-time rows.", call. = FALSE)
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

#' Classify the assignment design
#'
#' `block`: all treated units start in the same period and stay treated to the
#' end of the panel. `staggered`: absorbing treatment with different start
#' periods. `general`: treatment switches off for some unit.
#' @noRd
trop_design <- function(W, obs) {
  N <- nrow(W)
  TT <- ncol(W)
  treated <- which(rowSums(W == 1L & obs == 1L) > 0L)
  if (length(treated) == 0L) stop("No treated observations.", call. = FALSE)
  if (!any(W == 0L & obs == 1L)) stop("No untreated observed cells.", call. = FALSE)
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
  if (!is.matrix(Y) || !is.numeric(Y)) stop("`Y` must be a numeric matrix (units x periods).", call. = FALSE)
  if (!is.matrix(W) || !all(dim(W) == dim(Y))) stop("`W` must be a matrix with the same dimensions as `Y`.", call. = FALSE)
  if (is.logical(W)) W <- W * 1L
  if (anyNA(W) || !all(W %in% c(0, 1))) stop("`W` must contain only 0 and 1.", call. = FALSE)
  if (nrow(Y) < 3L || ncol(Y) < 3L) stop("The panel needs at least 3 units and 3 periods.", call. = FALSE)
  obs <- matrix(as.integer(is.finite(Y)), nrow(Y), ncol(Y))
  Y0 <- Y
  storage.mode(Y0) <- "double"
  Wi <- matrix(as.integer(W), nrow(W), ncol(W))
  list(Y = Y0, W = Wi, obs = obs)
}
