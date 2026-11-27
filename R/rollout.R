#' Make the treatment and the target from start dates
#'
#' `trop_rollout()` makes the treatment indicator, the target, and the event
#' time for a panel in long format. Each unit has one start: the first period
#' of the treatment. A unit without a start is never treated. The treatment
#' does not stop after it starts.
#'
#' With `anticipation = k`, the `k` periods before the start are also treated.
#' Thus, the fit does not use them as untreated cells. The event time is the
#' number of periods since the true start. Thus, the anticipation periods have
#' the event times `-k` to `-1`.
#'
#' The target contains the treated cells with an event time in `window`. With
#' the default window, the anticipation periods are not in the target. Thus,
#' they are not in the estimate. With a negative `window[1]`, the anticipation
#' periods in the window are in the target. Use this for placebo estimates of
#' the periods before the start.
#'
#' Give the columns to [trop()] with the arguments `treatment`, `target`, and
#' `event_time`. With `event_time`, [trop_effects()] and [trop_cells()] use
#' these event times. Then the groups are correct also when the panel has no
#' row for some periods, or when the anticipation periods start before the
#' panel. Do not give `anticipation` to these functions, because the event
#' times already count the anticipation periods. If you do not give
#' `event_time` to [trop()], give the same `anticipation` to [trop_effects()]
#' and [trop_cells()].
#'
#' The periods are the sorted values of `time`. The event time counts these
#' periods, not the difference of the time values. Thus, a panel with one row
#' for each month gives the event time in months.
#'
#' @param unit,time Vectors with one value for each row of the panel. A unit
#'   must have only one row for each period.
#' @param start A vector with the start of each row, in the units of `time`.
#'   The value must be the same for all rows of a unit. Use `NA` for a unit that
#'   is never treated, or that starts after the last period of the panel.
#'   Each other value must be one of the values of `time`.
#' @param anticipation The number of periods before the start that are treated.
#' @param window Two event times. The target contains the treated cells with
#'   an event time from `window[1]` to `window[2]`. The default contains all
#'   treated periods from the start.
#'
#' @return A data frame with one row for each row of the input and these
#'   columns: `treated` (0 or 1), `target` (`TRUE` or `FALSE`), and
#'   `event_time` (the number of periods since the start, or `NA` for a unit
#'   that is never treated).
#' @seealso [trop()], [trop_effects()], [trop_cells()]
#' @examples
#' panel <- expand.grid(unit = 1:4, month = 1:6)
#' starts <- c(3, 5, NA, NA)
#' panel$start <- starts[panel$unit]
#' r <- trop_rollout(panel$unit, panel$month, panel$start, anticipation = 1,
#'                   window = c(0, 1))
#' cbind(panel, r)[panel$unit == 1, ]
#' # Then, for example:
#' # trop(data, "y", "treated", "unit", "month", target = "target",
#' #      event_time = "event_time")
#' @export
trop_rollout <- function(unit, time, start, anticipation = 0L, window = c(0, Inf)) {
  n <- length(unit)
  if (length(time) != n || length(start) != n) {
    stop("`unit`, `time`, and `start` must have the same length.", call. = FALSE)
  }
  if (n == 0L) stop("`unit` must have at least one value.", call. = FALSE)
  if (anyNA(unit) || anyNA(time)) stop("`unit` and `time` must not have `NA` values.", call. = FALSE)
  anticipation <- check_anticipation(anticipation)
  if (!is.numeric(window) || length(window) != 2L || anyNA(window) || window[1L] > window[2L]) {
    stop("`window` must be 2 numbers, from the first to the last event time.", call. = FALSE)
  }
  times <- sort(unique(time))
  ti <- match(time, times)
  ui <- match(unit, unique(unit))
  if (anyDuplicated(cbind(ui, ti))) {
    stop("A unit has two or more rows for the same period.", call. = FALSE)
  }
  per_unit <- split(start, ui)
  mixed <- vapply(per_unit, function(s) length(unique(s)) > 1L, logical(1))
  if (any(mixed)) {
    bad <- unique(unit)[as.integer(names(per_unit)[mixed])]
    stop(sprintf("`start` must have one value for each unit. These units have more than one value: %s.",
                 paste(utils::head(bad, 3L), collapse = ", ")), call. = FALSE)
  }
  si <- match(start, times)
  unknown <- !is.na(start) & is.na(si)
  if (any(unknown)) {
    stop(sprintf("Each value of `start` must be `NA` or a value of `time`. This value is not: %s. For a unit that starts after the last period, use `NA`.",
                 format(start[unknown][1L])), call. = FALSE)
  }
  event <- ti - si
  treated <- !is.na(event) & event >= -anticipation
  target <- treated & event >= window[1L] & event <= window[2L]
  target[is.na(target)] <- FALSE
  data.frame(treated = as.integer(treated), target = target, event_time = as.integer(event))
}
