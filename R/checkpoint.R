# Checkpoint files for the bootstrap of trop_effects().

check_checkpoint_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`checkpoint` must be one file path.", call. = FALSE)
  }
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    stop(sprintf("The folder of the checkpoint file does not exist: %s", dir), call. = FALSE)
  }
  if (dir.exists(path)) stop("`checkpoint` must be a file, not a folder.", call. = FALSE)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  # Make sure that the folder is writable before the first expensive block.
  # The test file has a new random name. Thus, the test does not change or
  # remove a file that is already in the folder.
  probe <- own_temp_file(path, ".probe")
  why <- character()
  ok <- tryCatch(withCallingHandlers(file.create(probe), warning = function(w) {
    why <<- c(why, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    why <<- c(why, conditionMessage(e))
    FALSE
  })
  if (!isTRUE(ok)) {
    stop(sprintf("fasttrop cannot write to the folder of the checkpoint file: %s. %s",
                 dirname(path), paste(why, collapse = " ")), call. = FALSE)
  }
  unlink(probe)
  path
}

# A path for a new file in the folder of `path`. No file has this path. The
# name is short, so that it is not longer than the name of the checkpoint
# file in most cases.
own_temp_file <- function(path, ext) {
  repeat {
    f <- tempfile(pattern = ".ft", tmpdir = dirname(path), fileext = ext)
    if (!file.exists(f)) return(f)
  }
}

# A hash of an R object from the C code. It is for checkpoint files only.
fingerprint <- function(x) .Call(C_fingerprint, x)

# The fingerprint of a bootstrap of trop_effects(). The number of threads does
# not change the results, so it is not in the fingerprint.
effects_fingerprint <- function(fit, rows, gid, G, wt) {
  P <- fit$panel
  ctrl <- fit$control
  pkg <- package_info()
  fingerprint(list(
    format = 1L,
    package = c(pkg$version, pkg$sha),
    Y = P$Y, W = P$W, X = P$X, resid = isTRUE(P$resid),
    lambda = unname(as.double(fit$lambda)),
    control = as.double(c(ctrl$max_iter, ctrl$max_inner, ctrl$tol, ctrl$warm, ctrl$solver,
                          if (is.null(ctrl$fault)) 0 else ctrl$fault)),
    rows = rows, gid = gid, n_groups = as.integer(G),
    unit_weights = if (is.null(wt$u)) NULL else as.double(wt$u),
    per_unit = isTRUE(wt$per_unit)
  ))
}

new_boot_state <- function(key, B, G) {
  list(format = 1L, key = key, n_boot = as.integer(B), n_groups = as.integer(G),
       done = logical(B), att = rep(NA_real_, B),
       group_draws = matrix(NA_real_, B, G), diag = matrix(0L, B, 5L),
       package = package_info())
}

# The hash of the results in a state. A change of the results in the file
# after the save changes the hash.
payload_hash <- function(st) {
  fingerprint(list(st$done, st$att, st$group_draws, st$diag))
}

# A state that has the same key can still be damaged. Each part must have the
# correct type and size, and the hash of the results must agree. A draw that
# is not done has NA values and counts of 0. A draw that is done has an ATT
# only if it used a treated cell. A value that is not a number (NaN) shows a
# draw without an estimate. NA shows a draw that is not done.
valid_boot_state <- function(st, B, G) {
  shape <- is.logical(st$done) && length(st$done) == B && !anyNA(st$done) &&
    is.double(st$att) && length(st$att) == B &&
    is.matrix(st$group_draws) && is.double(st$group_draws) &&
    identical(dim(st$group_draws), c(as.integer(B), as.integer(G))) &&
    is.matrix(st$diag) && is.integer(st$diag) && identical(dim(st$diag), c(as.integer(B), 5L)) &&
    !anyNA(st$diag) && all(st$diag >= 0L) &&
    is.character(st$check) && length(st$check) == 1L
  if (!shape || !identical(st$check, payload_hash(st))) return(FALSE)
  d <- st$done
  gd <- st$group_draws
  blank <- function(v) is.na(v) & !is.nan(v)
  value <- function(v) is.finite(v) | is.nan(v)
  all(blank(st$att[!d])) && all(blank(gd[!d, ])) && all(st$diag[!d, ] == 0L) &&
    all(value(st$att[d])) && all(value(gd[d, ])) &&
    all(is.finite(st$att[d]) == (st$diag[d, 1L] > 0L)) &&
    all(is.nan(gd[d & st$diag[, 1L] == 0L, ]))
}

checkpoint_open <- function(path, key, B, G) {
  if (!file.exists(path)) return(new_boot_state(key, B, G))
  st <- tryCatch(readRDS(path), error = function(e) {
    stop(sprintf("fasttrop cannot read the checkpoint file %s: %s", path, conditionMessage(e)),
         call. = FALSE)
  })
  ok <- is.list(st) && identical(st$format, 1L) && is.character(st$key) &&
    is.logical(st$done) && is.matrix(st$group_draws) && is.matrix(st$diag)
  if (!ok) stop(sprintf("The file %s is not a fasttrop checkpoint file.", path), call. = FALSE)
  if (!identical(st$key, key) || !identical(st$n_boot, as.integer(B)) ||
      !identical(st$n_groups, as.integer(G))) {
    stop(sprintf("The checkpoint file %s is for a different calculation. Use a new file, or delete this file to start again.",
                 path), call. = FALSE)
  }
  if (!valid_boot_state(st, B, G)) {
    stop(sprintf("The checkpoint file %s is damaged: its parts do not have the correct types, sizes, or values. Delete the file to start again.",
                 path), call. = FALSE)
  }
  st
}

# Write the file in one step: first to a new temporary file in the same
# folder, then rename it. Thus, a stop during the write does not damage the
# file. The function removes only its own temporary file.
checkpoint_save <- function(state, path) {
  state$check <- payload_hash(state)
  tmp <- own_temp_file(path, ".tmp")
  why <- character()
  res <- tryCatch(withCallingHandlers(saveRDS(state, tmp), warning = function(w) {
    why <<- c(why, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(res, "error")) {
    unlink(tmp)
    stop(sprintf("fasttrop cannot write the checkpoint file %s: %s", path,
                 paste(c(why, conditionMessage(res)), collapse = " ")), call. = FALSE)
  }
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop(sprintf("fasttrop cannot replace the checkpoint file %s.", path), call. = FALSE)
  }
  invisible(path)
}
