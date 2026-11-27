#' TROP estimator for outcome and treatment matrices
#'
#' `trop_matrix()` is the matrix version of [trop()]. The rows are units and
#' the columns are periods.
#'
#' @param Y A numeric matrix of outcomes with one row for each unit and one
#'   column for each period. `NA` shows a cell without an outcome.
#' @param W A matrix of treatment indicators with the values 0 and 1. It must
#'   have the same dimensions as `Y`.
#' @param X Optional covariates. You can use one of three types. An array has
#'   the dimensions of `Y` and one more dimension for the covariates. A matrix
#'   with the dimensions of `Y` gives one covariate. A named list can contain
#'   such matrices. Refer to the section Covariates of [trop()].
#' @param event_time An optional numeric matrix with the dimensions of `Y`. It
#'   gives the event time of each treated cell. Refer to the section Event
#'   times of [trop()].
#' @inheritParams trop
#' @return An object of class `trop`. Refer to [trop()].
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(15 * 10), 15, 10)
#' W <- matrix(0L, 15, 10)
#' W[1:3, 8:10] <- 1L
#' fit <- trop_matrix(Y, W, lambda = c(time = 0, unit = 0, nn = Inf), se = "none",
#'                    control = trop_control(threads = 1))
#' coef(fit)
#' @export
trop_matrix <- function(Y, W, method = c("local", "global"), lambda = NULL, grid = trop_grid(),
                        cv_cells = NULL, n_placebo = 200L, placebo_size = 1L, max_cycles = 10L,
                        se = c("bootstrap", "none"), n_boot = 200L,
                        resample = c("stratified", "controls"), conf_level = 0.95,
                        ci = c("normal", "percentile"), control = trop_control(), seed = NULL,
                        X = NULL, covariate_method = c("joint", "residualize"),
                        target = NULL, strata = NULL, boot_rows = NULL, demean = FALSE,
                        event_time = NULL) {
  t0 <- Sys.time()
  plan <- trop_prepare(Y, W, method = method, lambda = lambda, grid = grid, cv_cells = cv_cells,
                       n_placebo = n_placebo, placebo_size = placebo_size, max_cycles = max_cycles,
                       se = se, n_boot = n_boot, resample = resample, conf_level = conf_level,
                       ci = ci, control = control, seed = seed, X = X,
                       covariate_method = covariate_method, target = target, strata = strata,
                       boot_rows = boot_rows, demean = demean, event_time = event_time)
  fit <- trop_run(plan)
  fit$call <- match.call()
  fit$timing$total_seconds <- as.double(difftime(Sys.time(), t0, units = "secs"))
  fit
}

#' Examine the inputs and make all random draws for a TROP fit
#'
#' The function returns a plan. The plan contains the panel matrices, the
#' design, and all sampled inputs. The sampled inputs are the cross-validation
#' cells, the placebo sets, and the bootstrap rows, as integer matrices with
#' indices that start at 0. Then `trop_run()` makes only deterministic calls to
#' the C code. Thus, the benchmark can write the plan to a file for the C
#' program.
#' @noRd
trop_prepare <- function(Y, W, method = c("local", "global"), lambda = NULL, grid = trop_grid(),
                         cv_cells = NULL, n_placebo = 200L, placebo_size = 1L, max_cycles = 10L,
                         se = c("bootstrap", "none"), n_boot = 200L,
                         resample = c("stratified", "controls"), conf_level = 0.95,
                         ci = c("normal", "percentile"), control = trop_control(), seed = NULL,
                         X = NULL, covariate_method = c("joint", "residualize"),
                         target = NULL, strata = NULL, boot_rows = NULL, demean = FALSE,
                         event_time = NULL) {
  method <- match.arg(method)
  covariate_method <- match.arg(covariate_method)
  se <- match.arg(se)
  resample <- match.arg(resample)
  ci <- match.arg(ci)
  grid_auto <- is.character(grid) && length(grid) == 1L && grid == "auto"
  if (is.character(grid) && !grid_auto) {
    stop("`grid` must be a `trop_grid()` object, a list of values, or \"auto\".", call. = FALSE)
  }
  if (!grid_auto && !inherits(grid, "trop_grid")) grid <- do.call(trop_grid, as.list(grid))
  if (!inherits(control, "trop_control")) control <- do.call(trop_control, as.list(control))
  check_count(max_cycles, "max_cycles")
  check_flag(demean, "demean")
  if (!is.numeric(conf_level) || length(conf_level) != 1L || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be more than 0 and less than 1.", call. = FALSE)
  }
  if (!is.null(strata) && !is.null(boot_rows)) {
    stop("Use `strata` or `boot_rows`, not both. `boot_rows` has the draws, and `strata` makes them.",
         call. = FALSE)
  }
  if (!is.null(target) && method == "global") {
    stop("`target` is only for `method = \"local\"`. The global method fits all treated cells together.",
         call. = FALSE)
  }
  if (!is.null(lambda)) lambda <- as_lambda(lambda)
  p <- as_panel_inputs(Y, W)
  event <- as_event_time(event_time, p$W, p$obs, rownames(Y))
  n_excluded <- 0L
  excluded <- NULL
  if (!is.null(target)) {
    target <- as_target(target, p$W)
    drop <- p$W == 1L & !target
    if (!any(p$W == 1L & target & p$obs == 1L)) {
      stop("`target` has no treated cell with an outcome.", call. = FALSE)
    }
    excluded <- drop & p$obs == 1L
    n_excluded <- sum(excluded)
    p$Y[drop] <- NA_real_
    p$obs[drop] <- 0L
  }
  if (demean) p$Y <- demean_units(p$Y, p$W, p$obs)
  if (grid_auto) grid <- if (is.null(lambda)) trop_grid_scaled(p$Y, p$W) else trop_grid()
  design <- trop_design(p$W, p$obs)
  cov <- NULL
  Xa <- as_covariate_array(X, p$Y, p$obs)
  if (!is.null(Xa)) {
    cov <- prepare_covariates(Xa, p$W, p$obs)
    cov$method <- covariate_method
  }
  N <- nrow(p$Y)
  never <- setdiff(seq_len(N), design$treated)
  if (method == "global") {
    if (design$type != "block") {
      stop("`method = \"global\"` is only for a block design. In a block design, all treated units start in the same period and stay treated. For other designs, use `method = \"local\"`.",
           call. = FALSE)
    }
    if (length(never) < 2L) {
      stop("`method = \"global\"` must have at least 2 units that are never treated.", call. = FALSE)
    }
  }
  if (design$type != "block") {
    warning("Theorem 5.1 of the paper (triple robustness) has a proof only for a block design. This design is ",
            design$type, ".", call. = FALSE)
  }
  if (design$type == "general") {
    warning("The treatment stops for some units. The estimator assumes no carryover effects and no dynamic effects.",
            call. = FALSE)
  }
  if (se == "bootstrap") {
    if (is.null(boot_rows)) check_count(n_boot, "n_boot")
    if (design$n_treated_units < 5L) {
      warning("The bootstrap of Algorithm 3 is correct only with many treated units. This panel has only ",
              design$n_treated_units, " treated unit(s).", call. = FALSE)
    }
  }
  if (!is.null(strata)) strata <- check_strata(strata, N)
  if (se == "bootstrap" && !is.null(boot_rows)) boot_rows <- check_boot_rows(boot_rows, N)

  rng <- seed_scope(seed)
  on.exit(rng$restore(), add = TRUE)

  # With "residualize", the C code calculates phi once. Then all later steps
  # use Y - X phi and no covariates, but each bootstrap draw calculates phi
  # again. With "joint", all fits estimate phi together with the other terms.
  Xp <- NULL
  p$Yb <- p$Y
  Xb <- NULL
  resid <- FALSE
  if (!is.null(cov$X)) {
    Xb <- cov$X
    if (covariate_method == "residualize") {
      cov$phi_scaled <- .Call(C_residualize, p$Y, p$W, p$obs, cov$X, unclass(control))
      p$Y <- subtract_xphi(p$Y, cov$X, cov$phi_scaled)
      resid <- TRUE
    } else {
      Xp <- cov$X
    }
  }

  cells <- NULL
  sets <- NULL
  Yc <- oc <- Wz <- Xc <- NULL
  init <- NULL
  size <- NA_integer_
  if (is.null(lambda)) {
    if (method == "local") {
      if (!is.null(cv_cells)) {
        check_count(cv_cells, "cv_cells")
        pool <- which(p$W == 0L & p$obs == 1L) - 1L
        if (cv_cells < length(pool)) cells <- as.integer(sort(pool[sample.int(length(pool), cv_cells)]))
      }
    } else {
      size <- as.integer(placebo_size)
      check_count(size, "placebo_size")
      check_count(n_placebo, "n_placebo")
      if (size >= length(never)) {
        stop("`placebo_size` must be less than the number of units that are never treated.",
             call. = FALSE)
      }
      Yc <- p$Y[never, , drop = FALSE]
      oc <- p$obs[never, , drop = FALSE]
      Wz <- matrix(0L, nrow(Yc), ncol(Yc))
      if (!is.null(Xp)) Xc <- Xp[never, , , drop = FALSE]
      sets <- matrix(0L, n_placebo, size)
      for (m in seq_len(n_placebo)) sets[m, ] <- sample.int(length(never), size)
      sets <- sets - 1L
      init <- c(mean(grid$time), mean(grid$unit), mean(grid$nn))
    }
  }
  rows <- NULL
  if (se == "bootstrap") {
    if (is.null(boot_rows)) {
      rows <- boot_rows_matrix(design$treated, never, as.integer(n_boot), resample, strata)
    } else {
      rows <- boot_rows
      resample <- "custom"
    }
  }
  list(p = p, design = design, never = never, method = method, lambda = lambda, grid = grid,
       max_cycles = as.integer(max_cycles), cells = cells, Yc = Yc, oc = oc, Wz = Wz, Xc = Xc,
       X = Xp, Xb = Xb, resid = resid, cov = cov,
       sets = sets, set_size = size, init = init, se = se,
       rows = if (is.null(rows)) NULL else rows - 1L, rows1 = rows, resample = resample,
       conf_level = conf_level, ci = ci, control = control, ctrl = unclass(control),
       n_excluded = n_excluded, excluded = excluded, event = event, demean = demean,
       stratified = !is.null(strata), dimnames = dimnames(Y))
}

#' Run the C calculations for a prepared plan
#' @noRd
trop_run <- function(plan) {
  t0 <- Sys.time()
  p <- plan$p
  design <- plan$design
  ctrl <- plan$ctrl
  kernel <- 0
  phase <- c(cv = 0, estimate = 0, bootstrap = 0)
  cv <- NULL
  lambda <- plan$lambda
  if (is.null(lambda)) {
    if (plan$method == "local") {
      res <- .Call(C_local_search, p$Y, p$W, p$obs, plan$X, plan$grid$time, plan$grid$unit,
                   plan$grid$nn, plan$cells, plan$max_cycles, ctrl)
      cv <- list(criterion = "loocv", score = res$score, n_evals = res$n_evals,
                 trace = trace_frame(res$trace),
                 n_cells = if (is.null(plan$cells)) sum(p$W == 0L & p$obs == 1L) else length(plan$cells),
                 fixed_point = NA)
    } else {
      res <- .Call(C_global_search, plan$Yc, plan$oc, plan$Wz, plan$Xc, as.integer(design$k), plan$sets,
                   plan$grid$time, plan$grid$unit, plan$grid$nn, as.double(plan$init),
                   plan$max_cycles, ctrl)
      if (!isTRUE(res$fixed_point)) {
        warning("The placebo cross-validation did not find a fixed point in `max_cycles` cycles.",
                call. = FALSE)
      }
      # The placebo cells of all sets: the cells with an outcome of the units
      # of each set in the last k periods.
      last <- seq.int(ncol(plan$oc) - design$k + 1L, ncol(plan$oc))
      cv <- list(criterion = "placebo_rmse", score = res$score, n_evals = res$n_evals,
                 trace = trace_frame(res$trace), n_sets = nrow(plan$sets),
                 set_size = plan$set_size,
                 n_cells = sum(plan$oc[c(plan$sets) + 1L, last, drop = FALSE]),
                 fixed_point = isTRUE(res$fixed_point))
    }
    kernel <- kernel + res$kernel_seconds
    phase[["cv"]] <- res$kernel_seconds
    lambda <- res$lambda
    cv <- c(cv, cv_counts(cv$trace, lambda))
    cv_check(cv, plan$method)
    edge <- grid_boundary(lambda, plan$grid)
    cv$boundary <- edge
    if (length(edge) > 0L) {
      warning("The cross-validation chose the ", paste(edge, collapse = " and the "),
              " of the grid. A value outside of the grid can fit better. Use a wider grid or `grid = \"auto\"`.",
              call. = FALSE)
    }
  }

  if (plan$method == "local") {
    est <- .Call(C_local_estimate, p$Y, p$W, p$obs, plan$X, lambda, ctrl)
    n_skipped <- est$n_skipped
    n_nonconv <- est$n_nonconverged
  } else {
    est <- .Call(C_global_estimate, p$Y, p$W, p$obs, plan$X, as.integer(design$treated - 1L),
                 as.integer(design$k), lambda, ctrl)
    n_skipped <- est$n_skipped
    n_nonconv <- as.integer(!isTRUE(est$converged) && est$n_failed == 0L &&
                              any(is.finite(est$tau)))
  }
  n_failed <- est$n_failed
  kernel <- kernel + est$kernel_seconds
  phase[["estimate"]] <- est$kernel_seconds
  if (n_failed > 0L) {
    warning(n_failed, " fit(s) failed: the singular value decomposition failed or the result is not finite. The estimate does not use the cells of these fits.",
            call. = FALSE)
  }
  if (n_skipped > 0L) {
    warning("The estimator did not use ", n_skipped, " treated cell(s). The weighted untreated cells do not identify the unit effect and the period effect of these cells.",
            call. = FALSE)
  }
  if (n_nonconv > 0L) {
    warning(n_nonconv, " fit(s) did not converge. Use a larger `max_iter` in trop_control().",
            call. = FALSE)
  }
  tau <- est$tau
  dimnames(tau) <- plan$dimnames

  covariates <- NULL
  if (!is.null(plan$cov)) {
    cv_ <- plan$cov
    phi_s <- if (is.null(cv_$X)) double(0) else if (cv_$method == "residualize") cv_$phi_scaled else est$phi
    coef_x <- stats::setNames(phi_s / cv_$scale, cv_$used)
    covariates <- list(coef = coef_x, method = cv_$method, names = cv_$names, used = cv_$used,
                       dropped = cv_$dropped)
  }
  boot <- NULL
  std_error <- NA_real_
  ci <- c(NA_real_, NA_real_)
  if (plan$se == "bootstrap") {
    eligible <- matrix(as.integer(is.finite(est$tau)), nrow(p$W), ncol(p$W))
    bs <- .Call(C_bootstrap, p$Yb, p$W, p$obs, plan$Xb, plan$resid, if (plan$method == "local") 0L else 1L,
                as.integer(if (plan$method == "global") design$k else 0L), lambda, plan$rows, ctrl,
                eligible, NULL, 0L, NULL, FALSE)
    kernel <- kernel + bs$kernel_seconds
    phase[["bootstrap"]] <- bs$kernel_seconds
    diagnostics <- boot_diagnostics(bs$diagnostics)
    boot_warnings(diagnostics)
    draws <- bs$draws
    ok <- draws[is.finite(draws)]
    if (length(ok) >= 2L) {
      std_error <- stats::sd(ok)
      a <- (1 - plan$conf_level) / 2
      ci <- if (plan$ci == "normal") {
        est$att + c(-1, 1) * stats::qnorm(1 - a) * std_error
      } else {
        unname(stats::quantile(ok, c(a, 1 - a), type = 7, names = FALSE))
      }
    }
    if (length(ok) < length(draws)) {
      warning(sprintf("In %d of %d bootstrap draws, no treated cell has an effect. The standard error and the interval do not use these draws.",
                      length(draws) - length(ok), length(draws)), call. = FALSE)
    }
    if (length(ok) < 2L) {
      warning("Fewer than 2 bootstrap draws have an estimate. The fit has no standard error and no interval.",
              call. = FALSE)
    }
    boot <- list(draws = draws, n_failed = length(draws) - length(ok), resample = plan$resample,
                 strata = plan$stratified, rows = plan$rows1, diagnostics = diagnostics)
  }

  structure(
    list(
      estimate = est$att, std.error = std_error, conf.low = ci[1L], conf.high = ci[2L],
      conf.level = plan$conf_level, ci_type = plan$ci, lambda = lambda, method = plan$method,
      design = design$type,
      n_units = nrow(p$Y), n_periods = ncol(p$Y), n_treated_units = design$n_treated_units,
      n_treated_cells = design$n_treated_cells, n_excluded_cells = plan$n_excluded,
      n_skipped = n_skipped,
      n_nonconverged = n_nonconv, n_failed = n_failed, tau = tau, cv = cv, boot = boot,
      covariates = covariates, control = plan$control, demean = plan$demean,
      panel = list(Y = p$Yb, W = p$W, X = plan$Xb, resid = plan$resid, excluded = plan$excluded,
                   event = plan$event),
      package = package_info(),
      timing = list(kernel_seconds = kernel,
                    total_seconds = as.double(difftime(Sys.time(), t0, units = "secs")),
                    cv_seconds = phase[["cv"]], estimate_seconds = phase[["estimate"]],
                    bootstrap_seconds = phase[["bootstrap"]]),
      units = plan$dimnames[[1L]], times = plan$dimnames[[2L]], call = NULL
    ),
    class = "trop"
  )
}

trace_frame <- function(m) {
  if (is.null(m) || length(m) == 0L) return(NULL)
  m <- matrix(m, ncol = 7L)
  data.frame(time = m[, 1], unit = m[, 2], nn = m[, 3], score = m[, 4],
             n_valid = as.integer(m[, 5]), n_nonconverged = as.integer(m[, 6]),
             n_failed = as.integer(m[, 7]))
}

# Totals of the fits of the cross-validation, and the counts at the chosen
# tuning parameters.
cv_counts <- function(trace, lambda) {
  if (is.null(trace)) return(list(n_nonconverged = 0L, n_failed = 0L))
  at <- which(trace$time == lambda[["time"]] & trace$unit == lambda[["unit"]] &
                trace$nn == lambda[["nn"]])[1L]
  list(n_nonconverged = sum(trace$n_nonconverged), n_failed = sum(trace$n_failed),
       n_nonconverged_chosen = if (is.na(at)) NA_integer_ else trace$n_nonconverged[at],
       n_failed_chosen = if (is.na(at)) NA_integer_ else trace$n_failed[at],
       n_valid_chosen = if (is.na(at)) NA_integer_ else trace$n_valid[at],
       n_valid_max = max(c(0L, trace$n_valid[is.finite(trace$score)])))
}

# Stop if no tuning value has a usable score. Give a warning when the chosen
# values have a score for fewer cells than the criterion has, and when fits
# failed or did not converge. For the global method, the cells are the placebo
# cells of all sets.
cv_check <- function(cv, method) {
  unit <- if (method == "local") "cells" else "placebo cells"
  total <- cv$n_cells
  if (!is.finite(cv$score)) {
    stop(sprintf("The cross-validation found no tuning values with a usable score. For each value in the grid, the fits of the %s failed, or the untreated cells do not identify them. Examine the grid and the panel, or set `lambda`.",
                 unit), call. = FALSE)
  }
  if (isTRUE(cv$n_valid_chosen < cv$n_valid_max)) {
    warning(sprintf("The chosen tuning values have a score for %d of %d %s. Other values in the grid have a score for %d. Thus, the scores are not comparable.",
                    cv$n_valid_chosen, total, unit, cv$n_valid_max), call. = FALSE)
  } else if (isTRUE(cv$n_valid_chosen < total)) {
    warning(sprintf("The chosen tuning values have a score for %d of %d %s. The untreated cells do not identify the other %s, or their fits failed.",
                    cv$n_valid_chosen, total, unit, unit), call. = FALSE)
  }
  cv_warnings(cv, method)
}

cv_warnings <- function(cv, method) {
  what <- if (method == "local") "leave-one-out fit(s)" else "placebo fit(s)"
  if (cv$n_failed > 0L) {
    warning(sprintf("In the cross-validation, %d %s failed (%s at the chosen values). The criterion does not use these fits.",
                    cv$n_failed, what, format(cv$n_failed_chosen)), call. = FALSE)
  }
  if (cv$n_nonconverged > 0L) {
    warning(sprintf("In the cross-validation, %d %s did not converge (%s at the chosen values). Use a larger `max_iter` in trop_control().",
                    cv$n_nonconverged, what, format(cv$n_nonconverged_chosen)), call. = FALSE)
  }
  invisible(NULL)
}

# Set a temporary seed. The function restore() gives the random number
# generator of the caller its earlier state again.
seed_scope <- function(seed) {
  if (is.null(seed)) return(list(restore = function() invisible(NULL)))
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  set.seed(seed)
  list(restore = function() {
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
}
