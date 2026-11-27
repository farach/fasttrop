#' Linear combinations of group effects
#'
#' `trop_contrast()` calculates linear combinations of the estimates of
#' [trop_effects()] or [trop()] objects. Examples are the mean effect of event
#' times 0 to 3, or the difference between two cohorts. The function uses the
#' bootstrap draws of the objects. It does not calculate the estimator again.
#'
#' For each combination, the estimate is the sum of the coefficients times the
#' estimates. Each draw gets the same sum of its draws. A draw without a value
#' for an estimate with a coefficient that is not 0 has no value for the
#' combination. The standard error is the standard deviation of the draws.
#'
#' All objects must have the same bootstrap draws. Draw b of each object must
#' use the same units. Refer to the section Joint draws of [trop_stepdown()].
#' The names of the estimates
#' are the names of [trop_stepdown()]. For one `trop_effects` object, the names
#' are the groups, for example `"0"` and `"1"` for event times. For a named
#' list, the names have the form `"event: 0"`.
#'
#' The result is a `trop_effects` object. Thus, [trop_stepdown()] can adjust
#' the p-values of the combinations, and `vcov()` gives their covariance.
#'
#' @param x A `trop_effects` object, a `trop` object, or a named list of such
#'   objects, all with a bootstrap.
#' @param contrasts A named numeric vector of coefficients, with the names of
#'   the estimates. A named list of such vectors gives one combination for each
#'   element. An estimate that is not in a vector has the coefficient 0.
#' @param conf_level The confidence level. The default is the level of the
#'   first object.
#' @param ci Use `"normal"` or `"percentile"`. Refer to [trop()].
#'
#' @return A data frame of class `trop_effects` with one row for each
#'   combination and the columns `contrast`, `estimate`, `std.error`,
#'   `conf.low`, `conf.high`, `band.low`, `band.high`, and `n_draws`. The
#'   attribute `coefficients` is the matrix of coefficients.
#' @seealso [trop_effects()], [trop_stepdown()]
#' @examples
#' sim <- trop_simulate(n_units = 20, n_periods = 12, n_treated = 5,
#'                      n_post = 3, effect = 1, seed = 1)
#' fit <- trop(sim, "y", "d", "unit", "time",
#'             lambda = c(time = 0.2, unit = 0.5, nn = 0.5),
#'             n_boot = 50, seed = 1, control = trop_control(threads = 1))
#' ev <- trop_effects(fit, by = "event_time")
#' trop_contrast(ev, list(mean_0_2 = c("0" = 1, "1" = 1, "2" = 1) / 3,
#'                        last_minus_first = c("2" = 1, "0" = -1)))
#' @export
trop_contrast <- function(x, contrasts, conf_level = NULL, ci = c("normal", "percentile")) {
  ci <- match.arg(ci)
  st <- stack_draws(x, "trop_contrast")
  if (is.null(conf_level)) conf_level <- st$conf_level
  if (!is.numeric(conf_level) || length(conf_level) != 1L || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be more than 0 and less than 1.", call. = FALSE)
  }
  L <- contrast_matrix(contrasts, st$labels)
  # Each combination uses only its own estimates. Thus, a missing estimate of
  # another combination has no effect on it.
  est <- vapply(seq_len(nrow(L)), function(k) {
    j <- which(L[k, ] != 0)
    if (any(!is.finite(st$estimate[j]))) NA_real_ else sum(L[k, j] * st$estimate[j])
  }, double(1))
  D <- st$draws
  draws <- matrix(NA_real_, nrow(D), nrow(L), dimnames = list(NULL, rownames(L)))
  for (k in seq_len(nrow(L))) {
    j <- which(L[k, ] != 0)
    ok <- rowSums(!is.finite(D[, j, drop = FALSE])) == 0L
    draws[ok, k] <- D[ok, j, drop = FALSE] %*% L[k, j]
  }
  out <- data.frame(contrast = rownames(L), estimate = unname(est))
  out <- add_boot_columns(out, unname(est), draws, st$rows, conf_level, ci)
  attr(out, "by") <- "contrast"
  attr(out, "conf.level") <- conf_level
  attr(out, "contrast_note") <- "Each row is a linear combination of estimates."
  attr(out, "coefficients") <- L
  attr(out, "units") <- st$units
  class(out) <- c("trop_effects", "data.frame")
  out
}

contrast_matrix <- function(contrasts, labels) {
  if (is.numeric(contrasts) && !is.list(contrasts)) contrasts <- list(contrast = contrasts)
  if (!is.list(contrasts) || length(contrasts) == 0L) {
    stop("`contrasts` must be a named numeric vector or a named list of such vectors.", call. = FALSE)
  }
  nm <- names(contrasts)
  if (is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop("Each element of `contrasts` must have a name, and the names must be different.",
         call. = FALSE)
  }
  L <- matrix(0, length(contrasts), length(labels), dimnames = list(nm, labels))
  for (k in seq_along(contrasts)) {
    v <- contrasts[[k]]
    if (!is.numeric(v) || length(v) == 0L || is.null(names(v)) || anyNA(v) || !all(is.finite(v))) {
      stop(sprintf("Contrast `%s` must be a numeric vector with names and finite values.", nm[k]),
           call. = FALSE)
    }
    if (anyDuplicated(names(v))) {
      stop(sprintf("Contrast `%s` has a name two or more times.", nm[k]), call. = FALSE)
    }
    miss <- setdiff(names(v), labels)
    if (length(miss) > 0L) {
      stop(sprintf("Contrast `%s` has names that are not estimates: %s. The estimates are: %s.",
                   nm[k], paste(miss, collapse = ", "),
                   paste(utils::head(labels, 12L), collapse = ", ")), call. = FALSE)
    }
    if (!any(v != 0)) stop(sprintf("Contrast `%s` has only coefficients of 0.", nm[k]), call. = FALSE)
    L[k, names(v)] <- v
  }
  L
}

#' Inputs for a sensitivity analysis of pre-trends
#'
#' `trop_honest_inputs()` makes the inputs of the 'HonestDiD' package from
#' event-time effects. The package gives confidence sets for the effects when
#' parallel trends is not exactly true.
#'
#' 'HonestDiD' uses estimates that are relative to a reference period. The
#' function subtracts the estimate of the `reference` event time from the
#' other estimates. It also calculates their covariance from the bootstrap
#' draws. The covariance uses only the draws with a value for all event times.
#'
#' The negative event times must be placebo estimates. To get them, set the
#' periods before the start to `W = 1` and keep them in the target. The
#' function [trop_rollout()] makes such a design with a negative `window[1]`. Give its
#' event times to [trop()] with `event_time`. Without `event_time`, give the
#' number of these periods to `anticipation` of [trop_effects()].
#'
#' In this design, the fit does not use these periods as untreated cells. Thus, the
#' effects after the start can be different from the effects of a fit that
#' uses these periods.
#'
#' 'HonestDiD' treats the estimates as consecutive periods. Thus, the event
#' times must have no gaps. The event times before the reference are the
#' periods before the treatment. The event times after the reference are the
#' periods after the treatment, with the anticipation periods after the
#' reference.
#'
#' By default, 'HonestDiD' gives a confidence set for the first period after
#' the reference. With `reference = -2`, this period is the event time -1, not
#' the event time 0. Thus, give `l_vec = inp$l_vec` to the functions of
#' 'HonestDiD'. The vector `l_vec` selects the mean of the event times in
#' `event_times`.
#'
#' @param x A `trop_effects` object with `by = "event_time"` and a bootstrap.
#' @param reference The event time of the reference period. It must be less
#'   than 0.
#' @param event_times The event times of the parameter of interest. They must
#'   be after the reference. The element `l_vec` of the result gives their
#'   mean.
#'
#' @return A list with `betahat`, `sigma`, `numPrePeriods`, and
#'   `numPostPeriods`, the arguments of the functions of 'HonestDiD'. It also
#'   has `l_vec` (the weights of the post periods for `event_times`),
#'   `event_time` (the event time of each element of `betahat`), `reference`,
#'   and `n_draws` (the number of draws in `sigma`).
#' @seealso [trop_effects()], [trop_rollout()]
#' @examples
#' sim <- trop_simulate(n_units = 24, n_periods = 10, n_treated = 6,
#'                      n_post = 3, effect = 1, seed = 1)
#' # The start of each treated unit is its first treated period.
#' start <- ave(ifelse(sim$d == 1, sim$time, NA), sim$unit,
#'              FUN = function(v) if (all(is.na(v))) NA else min(v, na.rm = TRUE))
#' # Two placebo periods before the start, in the target.
#' sim <- cbind(sim, trop_rollout(sim$unit, sim$time, start, anticipation = 2,
#'                                window = c(-2, Inf)))
#' fit <- trop(sim, "y", "treated", "unit", "time", target = "target",
#'             event_time = "event_time", lambda = c(time = 0.2, unit = 0.5, nn = 1),
#'             se = "none", control = trop_control(threads = 1))
#' ev <- trop_effects(fit, by = "event_time", n_boot = 30, seed = 1)
#' inp <- trop_honest_inputs(ev, reference = -1)
#' inp$betahat
#' inp$l_vec
#' @export
trop_honest_inputs <- function(x, reference = -1L, event_times = 0L) {
  if (!inherits(x, "trop_effects") || !identical(attr(x, "by"), "event_time")) {
    stop("`x` must be a `trop_effects` object with `by = \"event_time\"`.", call. = FALSE)
  }
  D <- effect_draws(x)
  if (is.null(D)) stop("`x` must have bootstrap draws.", call. = FALSE)
  if (!is.numeric(reference) || length(reference) != 1L || is.na(reference) ||
      reference != round(reference) || reference >= 0) {
    stop("`reference` must be one integer less than 0.", call. = FALSE)
  }
  e <- x$event_time
  if (!reference %in% e) {
    stop(sprintf("The reference event time %d is not in `x`. The event times are %s.",
                 as.integer(reference), paste(e, collapse = ", ")), call. = FALSE)
  }
  gaps <- setdiff(seq(min(e), max(e)), e)
  if (length(gaps) > 0L) {
    stop(sprintf("The event times of `x` must have no gaps. These event times are not in `x`: %s. 'HonestDiD' treats the estimates as consecutive periods. Use `window` in trop_effects() to remove the event times on one side of a gap.",
                 paste(gaps, collapse = ", ")), call. = FALSE)
  }
  pre <- e[e < reference]
  post <- e[e > reference]
  if (length(pre) == 0L) {
    stop("`x` must have at least one event time before the reference.", call. = FALSE)
  }
  if (length(post) == 0L) {
    stop("`x` must have at least one event time after the reference.", call. = FALSE)
  }
  if (!is.numeric(event_times) || length(event_times) == 0L || anyNA(event_times) ||
      any(event_times != round(event_times)) || anyDuplicated(event_times)) {
    stop("`event_times` must be different whole numbers.", call. = FALSE)
  }
  miss <- setdiff(event_times, post)
  if (length(miss) > 0L) {
    stop(sprintf("Each value of `event_times` must be an event time of `x` after the reference %d. These values are not: %s. The event times of `x` after the reference are: %s.",
                 as.integer(reference), paste(miss, collapse = ", "), paste(post, collapse = ", ")),
         call. = FALSE)
  }
  l_vec <- as.double(post %in% event_times) / length(event_times)
  names(l_vec) <- as.character(post)
  r <- which(e == reference)
  keep <- rowSums(!is.finite(D)) == 0L
  if (sum(keep) < 2L) {
    stop("Fewer than 2 draws have a value for all event times.", call. = FALSE)
  }
  other <- which(e != reference)
  Dd <- D[keep, other, drop = FALSE] - D[keep, r]
  beta <- x$estimate[other] - x$estimate[r]
  lab <- as.character(e[other])
  sigma <- stats::cov(Dd)
  dimnames(sigma) <- list(lab, lab)
  list(betahat = stats::setNames(beta, lab), sigma = sigma, numPrePeriods = length(pre),
       numPostPeriods = length(post), l_vec = l_vec, event_time = e[other],
       reference = as.integer(reference), n_draws = sum(keep))
}
