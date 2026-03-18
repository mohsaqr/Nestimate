#' Graphical VAR Estimation
#'
#' Estimate a graphical vector autoregressive (GVAR) model from time series or
#' panel data. Jointly estimates a sparse temporal network (L1-penalized VAR
#' coefficients) and a sparse contemporaneous network (graphical lasso on
#' residuals) using EBIC model selection over a lambda grid.
#'
#' Follows the two-step approach of Epskamp et al. (2018):
#' \enumerate{
#'   \item Estimate sparse temporal coefficients (beta) via L1-penalized
#'     regression across a lambda grid.
#'   \item For each temporal model, estimate the contemporaneous precision
#'     matrix (kappa) via graphical lasso on the residual covariance.
#'   \item Select the (lambda_beta, lambda_kappa) pair minimizing joint EBIC.
#' }
#'
#' @param data A data.frame or matrix with columns for variables, and optionally
#'   id, day, beep columns for panel/ESM data.
#' @param vars Character vector of variable names.
#' @param id Character. Name of the person-ID column. If NULL, assumes single
#'   subject.
#' @param day Character. Name of the day/session column. Default: NULL.
#' @param beep Character. Name of the beep/measurement column. Default: NULL.
#' @param n_lambda Integer. Number of lambda values per penalty dimension.
#'   Default: 30.
#' @param gamma Numeric. EBIC hyperparameter (0 = BIC, higher = sparser).
#'   Default: 0.5.
#' @param scale Logical. Whether to standardize variables. Default: TRUE.
#' @param center_within Logical. Whether to center within person when id is
#'   provided (removes between-person variance). Default: TRUE.
#' @param lambda_min_ratio Numeric. Ratio of min/max lambda for both beta and
#'   kappa grids. Default: 0.01.
#' @param penalize_diagonal Logical. Penalize autoregressive diagonal in beta.
#'   Default: TRUE.
#'
#' @return A list of class \code{gvar_result} containing:
#' \describe{
#'   \item{beta}{Temporal coefficient matrix (p x p). Each entry represents the
#'     effect of a variable at t-1 on another at t. Sparse (L1-penalized).}
#'   \item{kappa}{Precision matrix (p x p, symmetric). Inverse of residual
#'     covariance. Sparse (graphical lasso).}
#'   \item{PCC}{Partial contemporaneous correlations (p x p, symmetric).
#'     Derived from kappa: \code{-cov2cor(kappa)}, diagonal zeroed.}
#'   \item{PDC}{Partial directed correlations (p x p). Standardized temporal
#'     coefficients scaled by residual covariance.}
#'   \item{temporal}{Alias for \code{beta}.}
#'   \item{contemporaneous}{Alias for \code{PCC}.}
#'   \item{labels}{Variable names.}
#'   \item{n_obs}{Number of valid lag-pair observations.}
#'   \item{lambda_beta}{Selected lambda for temporal penalty.}
#'   \item{lambda_kappa}{Selected lambda for contemporaneous penalty.}
#'   \item{gamma}{EBIC gamma used.}
#'   \item{EBIC}{Best EBIC value.}
#' }
#'
#' @references
#' Epskamp, S., Waldorp, L. J., Mottus, R., & Borsboom, D. (2018).
#' The Gaussian Graphical Model in Cross-Sectional and Time-Series Data.
#' \emph{Multivariate Behavioral Research}, 53(4), 453-480.
#'
#' @examples
#' \dontrun{
#' # Single subject
#' set.seed(1)
#' data <- data.frame(matrix(rnorm(200 * 4), 200, 4))
#' names(data) <- c("V1", "V2", "V3", "V4")
#' res <- graphical_var(data, vars = names(data))
#' res$beta    # temporal
#' res$PCC     # contemporaneous
#' }
#'
#' @importFrom glasso glasso
#' @importFrom stats cov cov2cor
#' @export
graphical_var <- function(data,
                          vars,
                          id = NULL,
                          day = NULL,
                          beep = NULL,
                          n_lambda = 30L,
                          gamma = 0.5,
                          scale = TRUE,
                          center_within = TRUE,
                          lambda_min_ratio = 0.01,
                          penalize_diagonal = TRUE) {

  stopifnot(is.data.frame(data) || is.matrix(data))
  stopifnot(is.character(vars), length(vars) >= 2L)
  stopifnot(is.numeric(gamma), gamma >= 0)
  stopifnot(is.numeric(n_lambda), n_lambda >= 2L)

  data <- as.data.frame(data)
  d <- length(vars)

  # ---- 1. Build lag pairs ----
  lag_data <- .gvar_build_lag_pairs(data, vars, id, day, beep)
  Y <- lag_data$Y
  X <- lag_data$X
  n <- nrow(Y)

  if (n < d + 1L) {
    stop("Too few lag pairs (", n, ") for ", d, " variables.", call. = FALSE)
  }

  # ---- 2. Center / scale ----
  if (!is.null(id) && center_within) {
    ids <- lag_data$id_vec
    uid <- unique(ids)
    Y_means <- matrix(0, n, d)
    X_means <- matrix(0, n, d)
    for (u in uid) {
      idx <- which(ids == u)
      Y_means[idx, ] <- rep(colMeans(Y[idx, , drop = FALSE]), each = length(idx))
      X_means[idx, ] <- rep(colMeans(X[idx, , drop = FALSE]), each = length(idx))
    }
    Y <- Y - Y_means
    X <- X - X_means
  }

  if (scale) {
    Y_sd <- apply(Y, 2, sd)
    X_sd <- apply(X, 2, sd)
    Y_sd[Y_sd == 0] <- 1
    X_sd[X_sd == 0] <- 1
    Y <- t(t(Y) / Y_sd)
    X <- t(t(X) / X_sd)
  }

  # ---- 3-4. Fit temporal models across lambda grid ----
  cross_cor <- abs(crossprod(X, Y)) / n
  if (!penalize_diagonal) diag(cross_cor) <- 0
  lambda_max_beta <- max(cross_cor)
  if (lambda_max_beta == 0) lambda_max_beta <- 1
  lambda_beta_seq <- exp(seq(log(lambda_max_beta),
                              log(lambda_max_beta * lambda_min_ratio),
                              length.out = n_lambda))

  # Variable-by-variable lasso via coordinate descent
  beta_fits <- lapply(lambda_beta_seq, function(lam) {
    B <- .gvar_multivariate_lasso(X, Y, lam, penalize_diagonal)
    resid <- Y - X %*% B
    list(beta = B, residuals = resid)
  })

  # ---- 5. Lambda grid for kappa + joint EBIC search ----
  best_ebic <- Inf
  best_result <- NULL

  for (bi in seq_along(beta_fits)) {
    fit <- beta_fits[[bi]]
    S <- crossprod(fit$residuals) / n
    S <- (S + t(S)) / 2
    eig <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    if (min(eig) <= 0) S <- S + (abs(min(eig)) + 1e-6) * diag(d)

    # Lambda grid for kappa
    lambda_max_kappa <- max(abs(S[upper.tri(S)]))
    if (lambda_max_kappa == 0) lambda_max_kappa <- 1
    lam_kappa_seq <- exp(seq(log(lambda_max_kappa),
                              log(lambda_max_kappa * lambda_min_ratio),
                              length.out = n_lambda))

    prev_wi <- prev_w <- NULL
    for (lk in lam_kappa_seq) {
      gl <- tryCatch(
        glasso::glasso(S, rho = lk, penalize.diagonal = FALSE,
                       start = if (!is.null(prev_wi)) "warm" else "cold",
                       wi.init = prev_wi, w.init = prev_w),
        error = function(e) NULL
      )
      if (is.null(gl)) next
      prev_wi <- gl$wi
      prev_w <- gl$w

      kappa_fit <- gl$wi
      n_beta <- sum(fit$beta != 0)
      n_kappa <- sum(kappa_fit[upper.tri(kappa_fit)] != 0)
      n_params <- n_beta + n_kappa

      # Log-likelihood: multivariate normal
      loglik <- (n / 2) * (determinant(kappa_fit, logarithm = TRUE)$modulus -
                             sum(diag(S %*% kappa_fit)) - d * log(2 * pi))

      # EBIC (Foygel & Drton, 2010; graphicalVAR convention)
      ebic <- -2 * as.numeric(loglik) +
        n_params * log(n) +
        4 * gamma * (n_beta * log(d * d) + n_kappa * log(d))

      if (ebic < best_ebic) {
        best_ebic <- ebic
        best_result <- list(
          beta = fit$beta,
          kappa = kappa_fit,
          lambda_beta = bi,
          lambda_kappa = lk,
          ebic = ebic,
          n_obs = n
        )
      }
    }
  }

  if (is.null(best_result)) {
    stop("Graphical VAR estimation failed: no valid model found.", call. = FALSE)
  }

  # ---- 6. Extract PCC, PDC ----
  # Transpose beta to standard convention: beta[outcome, predictor]
  beta <- t(best_result$beta)
  kappa <- best_result$kappa
  colnames(beta) <- rownames(beta) <- vars
  colnames(kappa) <- rownames(kappa) <- vars

  pcc <- .gvar_compute_pcc(kappa)
  pdc <- .gvar_compute_pdc(beta, kappa)
  colnames(pcc) <- rownames(pcc) <- vars
  colnames(pdc) <- rownames(pdc) <- vars

  result <- list(
    beta            = beta,
    kappa           = kappa,
    PCC             = pcc,
    PDC             = pdc,
    temporal        = beta,
    contemporaneous = pcc,
    labels          = vars,
    n_obs           = best_result$n_obs,
    lambda_beta     = best_result$lambda_beta,
    lambda_kappa    = best_result$lambda_kappa,
    gamma           = gamma,
    EBIC            = best_result$ebic
  )
  class(result) <- "gvar_result"
  result
}


# ============================================================
# Internal: Lag pair construction
# ============================================================

#' Build lag-1 pairs for GVAR
#'
#' Handles single-subject and panel data with optional day/beep boundaries.
#' @param data data.frame with vars and optional id/day/beep columns.
#' @param vars character vector of variable names.
#' @param id character or NULL.
#' @param day character or NULL.
#' @param beep character or NULL.
#' @return List with Y (n x d outcome), X (n x d predictor), id_vec.
#' @noRd
.gvar_build_lag_pairs <- function(data, vars, id, day, beep) {
  d <- length(vars)

  # Coerce vars to numeric, drop NA rows
  check_cols <- c(vars, id, day, beep)
  check_cols <- check_cols[!is.null(check_cols)]
  complete <- complete.cases(data[, check_cols, drop = FALSE])
  data <- data[complete, , drop = FALSE]

  for (v in vars) data[[v]] <- as.numeric(data[[v]])

  n <- nrow(data)
  if (n < 3L) stop("Fewer than 3 complete rows.", call. = FALSE)

  # Sort
  if (!is.null(id)) {
    order_cols <- id
    if (!is.null(day)) order_cols <- c(order_cols, day)
    if (!is.null(beep)) order_cols <- c(order_cols, beep)
    data <- data[do.call(order, data[, order_cols, drop = FALSE]), ]
  }

  # Lag-1 pairing
  idx_t <- 2:n
  idx_lag <- 1:(n - 1L)

  # Same subject
  valid <- rep(TRUE, n - 1L)
  if (!is.null(id)) valid <- valid & (data[[id]][idx_t] == data[[id]][idx_lag])
  if (!is.null(day)) valid <- valid & (data[[day]][idx_t] == data[[day]][idx_lag])
  if (!is.null(beep)) {
    valid <- valid & (data[[beep]][idx_t] - data[[beep]][idx_lag] == 1L)
  }

  Y <- as.matrix(data[idx_t[valid], vars, drop = FALSE])
  X <- as.matrix(data[idx_lag[valid], vars, drop = FALSE])
  id_vec <- if (!is.null(id)) data[[id]][idx_t[valid]] else rep(1L, sum(valid))

  if (nrow(Y) < 3L) stop("Fewer than 3 valid lag pairs.", call. = FALSE)

  list(Y = Y, X = X, id_vec = id_vec)
}


# ============================================================
# Internal: Multivariate lasso via coordinate descent
# ============================================================

#' Fit multivariate lasso (variable-by-variable)
#'
#' Estimates beta matrix where Y_j = X * beta_j + epsilon_j for each
#' outcome variable j, using coordinate descent with L1 penalty.
#'
#' @param X n x p predictor matrix.
#' @param Y n x d outcome matrix.
#' @param lambda L1 penalty.
#' @param penalize_diagonal Logical. If FALSE, diagonal of beta is unpenalized.
#' @return d x d beta matrix.
#' @noRd
.gvar_multivariate_lasso <- function(X, Y, lambda, penalize_diagonal = TRUE) {
  n <- nrow(X)
  d <- ncol(X)
  beta <- matrix(0, d, d)

  # Precompute X'X / n and X'Y / n
  XtX <- crossprod(X) / n
  XtY <- crossprod(X, Y) / n
  diag_XtX <- diag(XtX)

  # Coordinate descent for each outcome variable (sequential by nature)
  for (j in seq_len(d)) {
    b <- rep(0, d)
    pen <- rep(lambda, d)
    if (!penalize_diagonal) pen[j] <- 0

    for (iter in seq_len(200L)) {
      b_old <- b
      for (k in seq_len(d)) {
        r_k <- XtY[k, j] - sum(XtX[k, ] * b) + XtX[k, k] * b[k]
        b[k] <- sign(r_k) * max(0, abs(r_k) - pen[k]) / diag_XtX[k]
      }
      if (max(abs(b - b_old)) < 1e-8) break
    }
    beta[, j] <- b
  }

  beta
}


# ============================================================
# Internal: PCC and PDC
# ============================================================

#' Partial Contemporaneous Correlations from precision matrix
#' @noRd
.gvar_compute_pcc <- function(kappa) {
  pcc <- -stats::cov2cor(kappa)
  diag(pcc) <- 0
  pcc <- (pcc + t(pcc)) / 2
  pcc
}

#' Partial Directed Correlations from beta and kappa
#' @noRd
.gvar_compute_pdc <- function(beta, kappa) {
  sigma <- tryCatch(solve(kappa), error = function(e) diag(nrow(kappa)))
  sigma_diag <- diag(sigma)
  kappa_diag <- diag(kappa)
  denom <- sqrt(sigma_diag %o% kappa_diag + beta^2)
  denom[denom == 0] <- 1
  pdc <- t(beta / denom)
  pdc
}


# ============================================================
# S3 Methods
# ============================================================

#' @export
print.gvar_result <- function(x, ...) {
  d <- length(x$labels)
  n_temp <- sum(x$beta != 0)
  n_contemp <- sum(x$PCC[upper.tri(x$PCC)] != 0)

  cat("Graphical VAR Result\n")
  cat(sprintf("  Variables:      %d (%s)\n", d, paste(x$labels, collapse = ", ")))
  cat(sprintf("  Observations:   %d\n", x$n_obs))
  cat(sprintf("  Temporal edges: %d / %d\n", n_temp, d * d))
  cat(sprintf("  Contemp edges:  %d / %d\n", n_contemp, d * (d - 1) / 2))
  cat(sprintf("  EBIC:           %.2f (gamma=%.2f)\n", x$EBIC, x$gamma))
  cat(sprintf("  Lambda:         beta=%.4f, kappa=%.4f\n",
              x$lambda_beta, x$lambda_kappa))
  invisible(x)
}

#' @export
summary.gvar_result <- function(object, ...) {
  cat("=== Temporal Network (beta) ===\n")
  print(round(object$beta, 4))
  cat("\n=== Contemporaneous Network (PCC) ===\n")
  print(round(object$PCC, 4))
  cat("\n=== Partial Directed Correlations (PDC) ===\n")
  print(round(object$PDC, 4))
  invisible(object)
}
