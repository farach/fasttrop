#' TROP estimator for outcome and treatment matrices
#'
#' Matrix interface to [trop()]. Rows are units and columns are periods.
#'
#' @param Y Numeric matrix of outcomes (units x periods). `NA` marks
#'   unobserved cells.
#' @param W Matrix of 0/1 treatment indicators with the same dimensions.
#' @inheritParams trop
#' @return An object of class `trop`.
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
                        ci = c("normal", "percentile"), control = trop_control(), seed = NULL) {
  t0 <- Sys.time()
  plan <- trop_prepare(Y, W, method = method, lambda = lambda, grid = grid, cv_cells = cv_cells,
                       n_placebo = n_placebo, placebo_size = placebo_size, max_cycles = max_cycles,
                       se = se, n_boot = n_boot, resample = resample, conf_level = conf_level,
                       ci = ci, control = control, seed = seed)
  fit <- trop_run(plan)
  fit$call <- match.call()
  fit$timing$total_seconds <- as.double(difftime(Sys.time(), t0, units = "secs"))
  fit
}

#' Validate inputs and draw every random quantity for a TROP fit
#'
#' Returns a plan with the panel matrices, design, and all sampled inputs
#' (cross-validation cells, placebo sets, bootstrap rows) as 0-based integer
#' matrices. [trop_run()] then makes only deterministic native calls, so the
#' plan can be serialized for the standalone C benchmark driver.
#' @noRd
trop_prepare <- function(Y, W, method = c("local", "global"), lambda = NULL, grid = trop_grid(),
                         cv_cells = NULL, n_placebo = 200L, placebo_size = 1L, max_cycles = 10L,
                         se = c("bootstrap", "none"), n_boot = 200L,
                         resample = c("stratified", "controls"), conf_level = 0.95,
                         ci = c("normal", "percentile"), control = trop_control(), seed = NULL) {
  method <- match.arg(method)
  se <- match.arg(se)
  resample <- match.arg(resample)
  ci <- match.arg(ci)
  if (!inherits(grid, "trop_grid")) grid <- do.call(trop_grid, as.list(grid))
  if (!inherits(control, "trop_control")) control <- do.call(trop_control, as.list(control))
  check_count(max_cycles, "max_cycles")
  if (!is.numeric(conf_level) || length(conf_level) != 1L || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be in (0, 1).", call. = FALSE)
  }
  if (!is.null(lambda)) lambda <- as_lambda(lambda)
  p <- as_panel_inputs(Y, W)
  design <- trop_design(p$W, p$obs)
  N <- nrow(p$Y)
  never <- setdiff(seq_len(N), design$treated)
  if (method == "global") {
    if (design$type != "block") {
      stop("`method = \"global\"` requires a block design (all treated units start together and stay treated). Use `method = \"local\"`.",
           call. = FALSE)
    }
    if (length(never) < 2L) stop("`method = \"global\"` needs at least 2 never-treated units.", call. = FALSE)
  }
  if (design$type != "block") {
    warning("TROP's triple-robustness result (Theorem 5.1) is proven for block assignment only; this design is ",
            design$type, ".", call. = FALSE)
  }
  if (design$type == "general") {
    warning("Treatment switches off for some units; the estimator assumes no carryover or dynamic effects.",
            call. = FALSE)
  }
  if (se == "bootstrap") {
    check_count(n_boot, "n_boot")
    if (design$n_treated_units < 5L) {
      warning("Bootstrap validity requires a growing number of treated units (paper Algorithm 3); only ",
              design$n_treated_units, " unit(s) are treated.", call. = FALSE)
    }
  }

  rng <- seed_scope(seed)
  on.exit(rng$restore(), add = TRUE)

  cells <- NULL
  sets <- NULL
  Yc <- oc <- Wz <- NULL
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
        stop("`placebo_size` must be smaller than the number of never-treated units.", call. = FALSE)
      }
      Yc <- p$Y[never, , drop = FALSE]
      oc <- p$obs[never, , drop = FALSE]
      Wz <- matrix(0L, nrow(Yc), ncol(Yc))
      sets <- matrix(0L, n_placebo, size)
      for (m in seq_len(n_placebo)) sets[m, ] <- sample.int(length(never), size)
      sets <- sets - 1L
      init <- c(mean(grid$time), mean(grid$unit), mean(grid$nn))
    }
  }
  rows <- NULL
  if (se == "bootstrap") {
    tr <- design$treated
    rows <- matrix(0L, n_boot, N)
    for (b in seq_len(n_boot)) {
      a <- if (resample == "stratified") tr[sample.int(length(tr), length(tr), replace = TRUE)] else tr
      rows[b, ] <- c(a, never[sample.int(length(never), length(never), replace = TRUE)])
    }
    rows <- rows - 1L
  }
  list(p = p, design = design, never = never, method = method, lambda = lambda, grid = grid,
       max_cycles = as.integer(max_cycles), cells = cells, Yc = Yc, oc = oc, Wz = Wz,
       sets = sets, set_size = size, init = init, se = se, rows = rows, resample = resample,
       conf_level = conf_level, ci = ci, control = control, ctrl = unclass(control),
       dimnames = dimnames(Y))
}

#' Run the native kernels for a prepared plan
#' @noRd
trop_run <- function(plan) {
  t0 <- Sys.time()
  p <- plan$p
  design <- plan$design
  ctrl <- plan$ctrl
  kernel <- 0
  cv <- NULL
  lambda <- plan$lambda
  if (is.null(lambda)) {
    if (plan$method == "local") {
      res <- .Call(C_local_search, p$Y, p$W, p$obs, plan$grid$time, plan$grid$unit,
                   plan$grid$nn, plan$cells, plan$max_cycles, ctrl)
      cv <- list(criterion = "loocv", score = res$score, n_evals = res$n_evals,
                 trace = trace_frame(res$trace),
                 n_cells = if (is.null(plan$cells)) sum(p$W == 0L & p$obs == 1L) else length(plan$cells),
                 fixed_point = NA)
    } else {
      res <- .Call(C_global_search, plan$Yc, plan$oc, plan$Wz, as.integer(design$k), plan$sets,
                   plan$grid$time, plan$grid$unit, plan$grid$nn, as.double(plan$init),
                   plan$max_cycles, ctrl)
      if (!isTRUE(res$fixed_point)) {
        warning("Placebo cross-validation did not reach a fixed point within `max_cycles`.",
                call. = FALSE)
      }
      cv <- list(criterion = "placebo_rmse", score = res$score, n_evals = res$n_evals,
                 trace = trace_frame(res$trace), n_sets = nrow(plan$sets),
                 set_size = plan$set_size, fixed_point = isTRUE(res$fixed_point))
    }
    kernel <- kernel + res$kernel_seconds
    lambda <- res$lambda
  }

  if (plan$method == "local") {
    est <- .Call(C_local_estimate, p$Y, p$W, p$obs, lambda, ctrl)
    n_skipped <- est$n_skipped
    n_nonconv <- est$n_nonconverged
  } else {
    est <- .Call(C_global_estimate, p$Y, p$W, p$obs, as.integer(design$treated - 1L),
                 as.integer(design$k), lambda, ctrl)
    n_skipped <- 0L
    n_nonconv <- as.integer(!isTRUE(est$converged))
  }
  kernel <- kernel + est$kernel_seconds
  if (n_skipped > 0L) {
    warning(n_skipped, " treated cell(s) were skipped because their unit and period effects are not identified by the weighted untreated cells.",
            call. = FALSE)
  }
  if (n_nonconv > 0L) {
    warning(n_nonconv, " fit(s) did not converge; consider a larger `max_iter` in trop_control().",
            call. = FALSE)
  }
  tau <- est$tau
  dimnames(tau) <- plan$dimnames

  boot <- NULL
  std_error <- NA_real_
  ci <- c(NA_real_, NA_real_)
  if (plan$se == "bootstrap") {
    bs <- .Call(C_bootstrap, p$Y, p$W, p$obs, if (plan$method == "local") 0L else 1L,
                as.integer(if (plan$method == "global") design$k else 0L), lambda, plan$rows, ctrl)
    kernel <- kernel + bs$kernel_seconds
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
      warning(length(draws) - length(ok), " of ", length(draws),
              " bootstrap draws failed and were excluded.", call. = FALSE)
    }
    boot <- list(draws = draws, n_failed = length(draws) - length(ok), resample = plan$resample)
  }

  structure(
    list(
      estimate = est$att, std.error = std_error, conf.low = ci[1L], conf.high = ci[2L],
      conf.level = plan$conf_level, ci_type = plan$ci, lambda = lambda, method = plan$method,
      design = design$type,
      n_units = nrow(p$Y), n_periods = ncol(p$Y), n_treated_units = design$n_treated_units,
      n_treated_cells = design$n_treated_cells, n_skipped = n_skipped,
      n_nonconverged = n_nonconv, tau = tau, cv = cv, boot = boot, control = plan$control,
      timing = list(kernel_seconds = kernel,
                    total_seconds = as.double(difftime(Sys.time(), t0, units = "secs"))),
      units = plan$dimnames[[1L]], times = plan$dimnames[[2L]], call = NULL
    ),
    class = "trop"
  )
}

trace_frame <- function(m) {
  if (is.null(m) || length(m) == 0L) return(NULL)
  m <- matrix(m, ncol = 4L)
  data.frame(time = m[, 1], unit = m[, 2], nn = m[, 3], score = m[, 4])
}

# Run code with a temporary seed and restore the caller's RNG state.
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
