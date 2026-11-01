# Covariates: Section 6.2 and Eq. 13 of the paper.
#
# The C core gets X as a double array with N x T x p entries. Covariate k of
# cell (i, t) is X[i, t, k]. The functions in this file examine X, remove
# columns that the fixed effects absorb, and scale the other columns. The
# scaling makes the pivot tolerance of the C code independent of the units of
# the covariates. trop_run() divides the coefficients by the scale again.

# Make an N x T x p array from long data. A row with an outcome must have
# finite values for all covariates.
covariate_array <- function(data, covariates, i, t, N, TT, units, times, has_y) {
  if (!is.character(covariates) || length(covariates) < 1L || anyNA(covariates)) {
    stop("`covariates` must be a character vector of column names.", call. = FALSE)
  }
  miss <- setdiff(covariates, names(data))
  if (length(miss)) {
    stop("These covariate columns are not in `data`: ", paste(miss, collapse = ", "), ".",
         call. = FALSE)
  }
  X <- array(0, c(N, TT, length(covariates)),
             dimnames = list(as.character(units), as.character(times), covariates))
  idx <- cbind(i, t)
  for (k in seq_along(covariates)) {
    v <- data[[covariates[k]]]
    if (is.logical(v)) v <- as.double(v)
    if (!is.numeric(v)) {
      stop(sprintf("Covariate `%s` must be numeric. For a factor, use model.matrix() to make indicator columns.",
                   covariates[k]), call. = FALSE)
    }
    if (any(!is.finite(v) & has_y)) {
      stop(sprintf("Covariate `%s` has `NA` or infinite values in rows with an outcome.",
                   covariates[k]), call. = FALSE)
    }
    v[!is.finite(v)] <- 0
    X[cbind(idx, k)] <- as.double(v)
  }
  X
}

# Examine X for trop_matrix(). A list of N x T matrices is also correct.
as_covariate_array <- function(X, Y, obs) {
  if (is.null(X)) return(NULL)
  if (is.list(X) && !is.data.frame(X)) {
    nms <- names(X)
    X <- lapply(X, function(m) {
      if (!is.matrix(m) || !all(dim(m) == dim(Y))) {
        stop("Each element of `X` must be a matrix with the same dimensions as `Y`.", call. = FALSE)
      }
      m
    })
    X <- array(unlist(X, use.names = FALSE), c(dim(Y), length(X)),
               dimnames = list(rownames(Y), colnames(Y), nms))
  }
  if (is.matrix(X)) X <- array(X, c(dim(X), 1L), dimnames = c(dimnames(X), list(NULL)))
  if (!is.array(X) || length(dim(X)) != 3L || !all(dim(X)[1:2] == dim(Y))) {
    stop("`X` must be an N x T x p array, an N x T matrix, or a list of N x T matrices.",
         call. = FALSE)
  }
  if (is.logical(X)) storage.mode(X) <- "double"
  if (!is.numeric(X)) stop("`X` must be numeric.", call. = FALSE)
  if (dim(X)[3L] < 1L) return(NULL)
  nm <- dimnames(X)[[3L]]
  if (is.null(nm)) nm <- paste0("x", seq_len(dim(X)[3L]))
  nm[is.na(nm) | !nzchar(nm)] <- paste0("x", which(is.na(nm) | !nzchar(nm)))
  has_y <- as.vector(obs == 1L)
  out <- array(0, dim(X))
  for (k in seq_len(dim(X)[3L])) {
    v <- as.double(X[, , k])
    if (any(!is.finite(v) & has_y)) {
      stop(sprintf("Covariate `%s` has `NA` or infinite values in cells with an outcome.", nm[k]),
           call. = FALSE)
    }
    v[!is.finite(v)] <- 0
    out[, , k] <- v
  }
  dimnames(out) <- list(NULL, NULL, nm)
  out
}

# Remove the unit and period means from x on the cells in `mask`. The function
# uses alternating projections. Thus, it is correct for an unbalanced mask.
# The attributes "a" and "b" give the unit and period effects.
within_fe <- function(x, mask, tol = 1e-12, max_iter = 1000L) {
  r <- ifelse(mask, x, 0)
  tol <- tol * max(1, max(abs(r)))
  nu <- pmax(rowSums(mask), 1)
  nt <- pmax(colSums(mask), 1)
  A <- double(nrow(r))
  B <- double(ncol(r))
  for (it in seq_len(max_iter)) {
    a <- rowSums(r) / nu
    r <- (r - a) * mask
    b <- colSums(r) / nt
    r <- sweep(r, 2L, b) * mask
    A <- A + a
    B <- B + b
    if (max(abs(a)) < tol && max(abs(b)) < tol) break
  }
  structure(r, a = A, b = B)
}

# Remove the covariates that the fixed effects absorb on the untreated cells,
# and the covariates that are linear combinations of other covariates. Then
# scale the other covariates.
prepare_covariates <- function(X, W, obs) {
  nm <- dimnames(X)[[3L]]
  mask <- W == 0L & obs == 1L
  n <- sum(mask)
  fe <- lapply(seq_len(dim(X)[3L]), function(k) within_fe(X[, , k], mask))
  R <- matrix(vapply(fe, function(r) r[mask], double(n)), nrow = n)
  raw <- vapply(seq_len(dim(X)[3L]), function(k) {
    v <- X[, , k][mask]
    sqrt(sum((v - mean(v))^2))
  }, double(1))
  ss <- sqrt(colSums(R^2))
  keep <- ss > 1e-8 * pmax(raw, 1e-300) & ss > 0
  dropped <- nm[!keep]
  if (any(keep)) {
    qr_r <- qr(sweep(R[, keep, drop = FALSE], 2L, ss[keep], "/"), tol = 1e-8)
    piv <- which(keep)
    ok <- piv[sort(qr_r$pivot[seq_len(qr_r$rank)])]
    dropped <- c(dropped, nm[setdiff(piv, ok)])
    keep <- seq_along(nm) %in% ok
  }
  if (length(dropped)) {
    warning("The estimator does not use these covariates. On the untreated cells, the unit effects, the period effects, or other covariates explain them completely: ",
            paste(dropped, collapse = ", "), ".", call. = FALSE)
  }
  if (!any(keep)) {
    return(list(X = NULL, names = nm, used = character(0), dropped = dropped,
                scale = double(0)))
  }
  # The fixed effects absorb a term a_i + b_t in X in each fit. Thus, the
  # estimates do not change when the function removes the unit and period
  # effects of the untreated cells from all cells. After this step, the rank
  # test of the C code measures only the collinearity between the covariates.
  sc <- ss[keep] / sqrt(n)
  Xs <- X[, , keep, drop = FALSE]
  for (j in seq_len(sum(keep))) {
    k <- which(keep)[j]
    v <- Xs[, , j] - attr(fe[[k]], "a") - rep(attr(fe[[k]], "b"), each = nrow(mask))
    Xs[, , j] <- v / sc[j]
  }
  Xs[!is.finite(Xs)] <- 0
  dim(Xs) <- c(dim(X)[1:2], sum(keep))
  storage.mode(Xs) <- "double"
  list(X = Xs, names = nm, used = nm[keep], dropped = dropped, scale = sc)
}

# Subtract X %*% phi from Y. NA outcomes stay NA.
subtract_xphi <- function(Y, X, phi) {
  for (k in seq_along(phi)) Y <- Y - phi[k] * X[, , k]
  Y
}
