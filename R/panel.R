#' Make outcome and treatment matrices from long panel data
#'
#' @param data A data frame with one row for each unit and period.
#' @param outcome,treatment,unit,time Column names.
#' @param covariates Optional covariate column names.
#' @return A list with `Y` (double, with `NA` for a cell without an outcome),
#'   `W` (integer, 0 or 1), `X` (an N x T x p array or `NULL`), and the sorted
#'   `units` and `times`. A cell without a row gets a treatment from the rows
#'   of its unit. Refer to fill_missing_treatment(). The element `ambiguous`
#'   has the units with an uncertain start.
#' @noRd
panel_matrices <- function(data, outcome, treatment, unit, time, covariates = NULL,
                           target = NULL, strata = NULL, event_time = NULL) {
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
  present <- matrix(FALSE, N, TT)
  present[idx] <- TRUE
  event_m <- NULL
  if (!is.null(event_time)) {
    if (!is.character(event_time) || length(event_time) != 1L || !event_time %in% names(data)) {
      stop("`event_time` must be the name of one column in `data`.", call. = FALSE)
    }
    ev <- data[[event_time]]
    if (!is.numeric(ev) && !all(is.na(ev))) {
      stop(sprintf("Column `%s` must be numeric.", event_time), call. = FALSE)
    }
    event_m <- matrix(NA_real_, N, TT, dimnames = dimnames(Y))
    event_m[idx] <- as.double(ev)
  }
  fill <- fill_missing_treatment(W, present, event_m)
  W <- fill$W
  X <- NULL
  if (!is.null(covariates)) {
    X <- covariate_array(data, covariates, i, t, N, TT, units, times, is.finite(y))
  }
  target_m <- NULL
  if (!is.null(target)) {
    if (!is.character(target) || length(target) != 1L || !target %in% names(data)) {
      stop("`target` must be the name of one column in `data`.", call. = FALSE)
    }
    tv <- data[[target]]
    if (is.numeric(tv) && !anyNA(tv) && all(tv %in% c(0, 1))) tv <- tv == 1
    if (!is.logical(tv) || anyNA(tv)) {
      stop(sprintf("Column `%s` must contain only TRUE and FALSE, or 1 and 0, and no `NA` values.",
                   target), call. = FALSE)
    }
    target_m <- matrix(TRUE, N, TT, dimnames = dimnames(Y))
    target_m[idx] <- tv
  }
  strata_v <- NULL
  if (!is.null(strata)) {
    if (!is.character(strata) || length(strata) != 1L || !strata %in% names(data)) {
      stop("`strata` must be the name of one column in `data`.", call. = FALSE)
    }
    sv <- data[[strata]]
    if (anyNA(sv)) stop(sprintf("Column `%s` must not have `NA` values.", strata), call. = FALSE)
    strata_v <- sv[match(seq_len(N), i)]
    if (any(sv != strata_v[i])) {
      stop(sprintf("Column `%s` must have one value for each unit.", strata), call. = FALSE)
    }
  }
  list(Y = Y, W = W, X = X, units = units, times = times, target = target_m, strata = strata_v,
       event = event_m, ambiguous = units[fill$ambiguous])
}

# The treatment of the cells without a row. A cell gets the value 1 if the
# nearest row before it has the value 1, and the nearest row after it has the
# value 1 or does not exist. Otherwise, it gets 0. Thus, a missing row does not
# stop a treatment run. With event times, the treated rows after the cell can
# show that their run starts after the cell. Then the treatment stopped and
# started again, and the cell gets 0. The start of a run is not certain if the
# cells just before its first treated row have no row. The element `ambiguous`
# has the numbers of these units.
fill_missing_treatment <- function(W, present, event = NULL) {
  if (all(present)) return(list(W = W, ambiguous = integer()))
  N <- nrow(W)
  TT <- ncol(W)
  prev_w <- next_w <- matrix(NA_integer_, N, TT)
  last <- rep(NA_integer_, N)
  for (t in seq_len(TT)) {
    prev_w[, t] <- last
    last[present[, t]] <- W[present[, t], t]
  }
  last <- rep(NA_integer_, N)
  for (t in rev(seq_len(TT))) {
    next_w[, t] <- last
    last[present[, t]] <- W[present[, t], t]
  }
  miss <- !present
  one <- miss & !is.na(prev_w) & prev_w == 1L & (is.na(next_w) | next_w == 1L)
  if (!is.null(event) && any(one)) {
    later <- next_run_start(W, present, event) > col(W)
    one <- one & !(!is.na(later) & later)
  }
  W[one] <- 1L
  amb <- miss & (is.na(prev_w) | prev_w == 0L) & !is.na(next_w) & next_w == 1L
  list(W = W, ambiguous = which(rowSums(amb) > 0L))
}

# For each cell, the start (period minus event time) of the treated rows that
# come next. The treated rows next to each other form a segment. A segment
# gets the start of any of its rows with an event time. A cell without a row
# gets the start of the next segment, or NA.
next_run_start <- function(W, present, event) {
  N <- nrow(W)
  TT <- ncol(W)
  tr <- present & W == 1L
  first <- tr & cbind(FALSE, !tr[, -TT, drop = FALSE])
  seg <- t(apply(first, 1L, cumsum))
  if (TT == 1L) seg <- matrix(seg, N, 1L)
  key <- (row(W) - 1L) * (TT + 1L) + seg
  key[!tr] <- NA
  start <- col(W) - event
  known <- tr & !is.na(start)
  seg_start <- start[known][match(key, key[known])]
  dim(seg_start) <- dim(W)
  out <- matrix(NA_real_, N, TT)
  last <- rep(NA_real_, N)
  for (t in rev(seq_len(TT))) {
    out[, t] <- last
    on <- present[, t]
    last[on] <- seg_start[on, t]
  }
  out
}

# The event times of the cells, as an integer matrix. Each treated cell with an
# outcome must have a value. In a run of treated periods, the event time must
# increase by 1 from one period to the next. Thus, all cells of a run have the
# same start. A new run can have a different start.
as_event_time <- function(event_time, W, obs, units = NULL) {
  if (is.null(event_time)) return(NULL)
  if (!is.matrix(event_time) || !(is.numeric(event_time) || all(is.na(event_time))) ||
      !all(dim(event_time) == dim(W))) {
    stop("`event_time` must be a numeric matrix with the dimensions of `Y`.", call. = FALSE)
  }
  ok <- is.finite(event_time)
  if (any(event_time[ok] != round(event_time[ok])) || any(abs(event_time[ok]) > 1e6)) {
    stop("The values of `event_time` must be whole numbers from -1e6 to 1e6.", call. = FALSE)
  }
  need <- W == 1L & obs == 1L & !ok
  if (any(need)) {
    stop(sprintf("`event_time` has no value for %d treated cell(s) with an outcome. Give a value for each treated cell.",
                 sum(need)), call. = FALSE)
  }
  out <- matrix(NA_integer_, nrow(W), ncol(W))
  out[ok] <- as.integer(event_time[ok])
  sel <- W == 1L & ok
  if (any(sel)) {
    run <- ((row(W) - 1L) * ncol(W) + run_start(W))[sel]
    start <- (col(W) - out)[sel]
    bad <- start != start[match(run, run)]
    if (any(bad)) {
      rows <- unique(row(W)[sel][bad])
      who <- if (is.null(units)) rows else units[rows]
      stop(sprintf("In %d unit(s), the event times do not increase by 1 from one period to the next in a run of treated periods. The event time counts the periods of the panel. A period without a row for any unit is not a period of the panel. Units: %s.",
                   length(rows), paste(utils::head(who, 5L), collapse = ", ")), call. = FALSE)
    }
  }
  out
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

as_target <- function(target, W) {
  if (is.numeric(target) && !anyNA(target) && all(target %in% c(0, 1))) target <- target == 1
  if (!is.matrix(target) || !is.logical(target) || !all(dim(target) == dim(W)) || anyNA(target)) {
    stop("`target` must be a logical matrix with the dimensions of `Y`, and no `NA` values.",
         call. = FALSE)
  }
  target
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
