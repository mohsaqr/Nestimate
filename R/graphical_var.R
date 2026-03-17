#' Graphical VAR Estimation
#'
#' @description
#' Estimate a graphical vector autoregressive (GVAR) model from time series or
#' panel data. Jointly estimates a sparse temporal network (L1-penalized VAR
#' coefficients) and a sparse contemporaneous network (GLASSO on residuals)
#' using EBIC model selection over a lambda grid.
#'
#' Follows the algorithm of Epskamp et al. (2018): alternating optimization
#' of temporal (beta) and contemporaneous (kappa) parameters, with joint EBIC
#' for model selection.
#'
#' @param data A data.frame with columns for variables, and optionally id, day,
#'   beep columns for panel/ESM data.
#' @param vars Character vector of variable names.
#' @param id Character. Name of the person-ID column. If NULL, assumes single
#'   subject.
#' @param day Character. Name of the day/session column. Default: NULL.
#' @param beep Character. Name of the beep/measurement column. Default: NULL.
#' @param n_lambda Integer. Number of lambda values per penalty dimension.
#'   Creates an n_lambda x n_lambda grid. Default: 30.
#' @param gamma Numeric. EBIC hyperparameter (0 = BIC, higher = sparser).
#'   Default: 0.5.
#' @param scale Logical. Whether to standardize variables. Default: TRUE.
#' @param center_within Logical. Whether to center within person (removes
#'   between-person variance). Default: TRUE.
#' @param maxit_in Integer. Max inner iterations for beta update. Default: 100.
#' @param maxit_out Integer. Max outer alternating iterations. Default: 100.
#' @param lambda_min_kappa Numeric. Ratio of min/max lambda for kappa.
#'   Default: 0.01.
#' @param lambda_min_beta Numeric. Ratio of min/max lambda for beta.
#'   Default: 0.01.
#' @param penalize_diagonal Logical. Penalize autoregressive diagonal in beta.
#'   Default: TRUE.
#'
#' @return A list of class \code{gvar_result} containing:
#' \describe{
#'   \item{beta}{Temporal coefficient matrix (p x p). \code{beta[i,j]} =
#'     effect of variable j at t-1 on variable i at t. Sparse (L1-penalized).}
#'   \item{kappa}{Precision matrix (p x p, symmetric). Inverse of residual
#'     covariance. Sparse (GLASSO-penalized).}
#'   \item{PCC}{Partial contemporaneous correlations (p x p, symmetric).
#'     Derived from kappa: \code{-cov2cor(kappa)}, diagonal zeroed.}
#'   \item{PDC}{Partial directed correlations (p x p). Standardized temporal
#'     coefficients scaled by residual covariance.}
#'   \item{labels}{Variable names.}
#'   \item{n_obs}{Number of valid observations used.}
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
#' data <- matrix(rnorm(500 * 4), 500, 4)
#' colnames(data) <- c("V1", "V2", "V3", "V4")
#' res <- graphical_var(as.data.frame(data), vars = colnames(data))
#' res$beta    # temporal
#' res$PCC     # contemporaneous
#' }
#'
#' @importFrom glasso glasso
#' @importFrom stats cov2cor lm.fit
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
                          maxit_in = 100L,
                          maxit_out = 100L,
                          lambda_min_kappa = 0.05,
                          lambda_min_beta = lambda_min_kappa,
                          penalize_diagonal = TRUE) {

  stopifnot(is.data.frame(data) || is.matrix(data))
  stopifnot(is.character(vars), length(vars) >= 2L)
  stopifnot(is.numeric(gamma), gamma >= 0)
  stopifnot(is.numeric(n_lambda), n_lambda >= 2L)

  if (!requireNamespace("graphicalVAR", quietly = TRUE)) {
    stop("Package 'graphicalVAR' is required. Install with: ",
         "install.packages('graphicalVAR')", call. = FALSE)
  }

  data <- as.data.frame(data)
  d <- length(vars)

  # ---- Build arguments for graphicalVAR::graphicalVAR ----
  gvar_args <- list(
    data       = data,
    vars       = vars,
    nLambda    = as.integer(n_lambda),
    gamma      = gamma,
    scale      = scale,
    verbose    = FALSE,
    maxit.in   = as.integer(maxit_in),
    maxit.out  = as.integer(maxit_out),
    lambda_min_kappa    = lambda_min_kappa,
    lambda_min_beta     = lambda_min_beta,
    penalize.diagonal   = penalize_diagonal,
    centerWithin        = center_within
  )
  if (!is.null(id))   gvar_args$idvar   <- id
  if (!is.null(day))  gvar_args$dayvar  <- day
  if (!is.null(beep)) gvar_args$beepvar <- beep

  # ---- Run graphicalVAR ----
  fit <- tryCatch(
    suppressWarnings(do.call(graphicalVAR::graphicalVAR, gvar_args)),
    error = function(e) {
      stop("graphicalVAR estimation failed: ", e$message, call. = FALSE)
    }
  )

  # ---- Extract and repackage ----
  beta_full <- fit$beta  # d x (d+1), col 1 = intercept
  beta_net <- beta_full[, -1, drop = FALSE]
  colnames(beta_net) <- rownames(beta_net) <- vars

  kappa <- fit$kappa
  colnames(kappa) <- rownames(kappa) <- vars

  pcc <- fit$PCC
  colnames(pcc) <- rownames(pcc) <- vars

  pdc <- fit$PDC
  colnames(pdc) <- rownames(pdc) <- vars

  # Extract selected lambdas from the best result
  all_ebics <- vapply(fit$allResults, `[[`, numeric(1), "EBIC")
  best_idx <- which.min(all_ebics)

  result <- list(
    beta            = beta_net,
    kappa           = kappa,
    PCC             = pcc,
    PDC             = pdc,
    labels          = vars,
    n_obs           = fit$N,
    lambda_beta     = fit$path$beta[best_idx],
    lambda_kappa    = fit$path$kappa[best_idx],
    gamma           = gamma,
    EBIC            = min(all_ebics, na.rm = TRUE)
  )
  class(result) <- "gvar_result"
  result
}




#' Multilevel Graphical VAR (Panel GVAR)
#'
#' @description
#' Estimate a multilevel graphical VAR model from panel/ESM data with multiple
#' subjects. Produces three network layers:
#' \itemize{
#'   \item \strong{Temporal}: Group-level VAR coefficients (fixed effects)
#'   \item \strong{Contemporaneous}: Group-level partial correlations among
#'     within-person residuals
#'   \item \strong{Between-subjects}: Partial correlations of person means
#' }
#' Optionally estimates individual (subject-level) networks.
#'
#' Wraps \code{graphicalVAR::mlGraphicalVAR()} for exact numerical equivalence.
#'
#' @param data A data.frame with columns for variables, id, and optionally
#'   day/beep.
#' @param vars Character vector of variable names.
#' @param id Character. Name of the person-ID column.
#' @param day Character. Day/session column. Default: NULL.
#' @param beep Character. Beep/measurement column. Default: NULL.
#' @param gamma Numeric. EBIC hyperparameter. Default: 0.5.
#' @param scale Logical. Standardize variables. Default: TRUE.
#' @param center_within Logical. Center within person. Default: TRUE.
#' @param subject_networks Logical. Estimate per-subject networks. Default: TRUE.
#' @param verbose Logical. Show progress. Default: FALSE.
#' @param ... Additional arguments passed to
#'   \code{graphicalVAR::mlGraphicalVAR()}.
#'
#' @return A list of class \code{ml_graphical_var_result} containing:
#' \describe{
#'   \item{temporal}{Group-level temporal coefficient matrix (p x p).}
#'   \item{PCC}{Group-level partial contemporaneous correlations (p x p).}
#'   \item{PDC}{Group-level partial directed correlations (p x p).}
#'   \item{between}{Between-subjects partial correlation network (p x p).}
#'   \item{subject_PCC}{List of per-subject contemporaneous networks (if
#'     \code{subject_networks = TRUE}).}
#'   \item{subject_PDC}{List of per-subject temporal networks (if
#'     \code{subject_networks = TRUE}).}
#'   \item{labels}{Variable names.}
#'   \item{ids}{Subject identifiers.}
#'   \item{n_subjects}{Number of subjects.}
#'   \item{gamma}{EBIC gamma used.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- ml_graphical_var(esm_data, vars = c("happy", "sad", "anxious"),
#'                   id = "subject", day = "day", beep = "beep")
#' res$temporal  # group-level temporal
#' res$PCC       # group-level contemporaneous
#' res$between   # between-subjects
#' }
#'
#' @export
ml_graphical_var <- function(data,
                       vars,
                       id,
                       day = NULL,
                       beep = NULL,
                       gamma = 0.5,
                       scale = TRUE,
                       center_within = TRUE,
                       subject_networks = TRUE,
                       verbose = FALSE,
                       ...) {

  if (!requireNamespace("graphicalVAR", quietly = TRUE)) {
    stop("Package 'graphicalVAR' is required. Install with: ",
         "install.packages('graphicalVAR')", call. = FALSE)
  }

  stopifnot(is.data.frame(data))
  stopifnot(is.character(vars), length(vars) >= 2L)
  stopifnot(is.character(id), length(id) == 1L, id %in% names(data))

  data <- as.data.frame(data)
  d <- length(vars)

  # Build arguments
  ml_args <- list(
    data            = data,
    vars            = vars,
    idvar           = id,
    gamma           = gamma,
    scale           = scale,
    centerWithin    = center_within,
    subjectNetworks = subject_networks,
    verbose         = verbose,
    ...
  )
  if (!is.null(day))  ml_args$dayvar  <- day
  if (!is.null(beep)) ml_args$beepvar <- beep

  fit <- tryCatch(
    suppressWarnings(do.call(graphicalVAR::mlGraphicalVAR, ml_args)),
    error = function(e) {
      stop("mlGraphicalVAR estimation failed: ", e$message, call. = FALSE)
    }
  )

  # Extract temporal beta (remove intercept column)
  beta_full <- fit$fixedResults$beta
  beta_net <- beta_full[, -1, drop = FALSE]
  colnames(beta_net) <- rownames(beta_net) <- vars

  # Fixed PCC / PDC
  pcc <- fit$fixedPCC
  colnames(pcc) <- rownames(pcc) <- vars
  pdc <- fit$fixedPDC
  colnames(pdc) <- rownames(pdc) <- vars

  # Between-subjects
  between <- fit$betweenNet
  colnames(between) <- rownames(between) <- vars

  # Subject-level networks
  subj_pcc <- NULL
  subj_pdc <- NULL
  if (subject_networks && !is.null(fit$subjectPCC)) {
    subj_pcc <- fit$subjectPCC
    subj_pdc <- fit$subjectPDC
    names(subj_pcc) <- names(subj_pdc) <- fit$ids
  }

  result <- list(
    temporal         = beta_net,
    PCC              = pcc,
    PDC              = pdc,
    between          = between,
    subject_PCC      = subj_pcc,
    subject_PDC      = subj_pdc,
    labels           = vars,
    ids              = fit$ids,
    n_subjects       = length(fit$ids),
    gamma            = gamma
  )
  class(result) <- "ml_graphical_var_result"
  result
}


#' @export
print.ml_graphical_var_result <- function(x, ...) {
  d <- length(x$labels)
  n_temp <- sum(x$temporal != 0)
  n_contemp <- sum(x$PCC[upper.tri(x$PCC)] != 0)
  n_between <- sum(x$between[upper.tri(x$between)] != 0)

  cat("Panel Graphical VAR Result\n")
  cat(sprintf("  Variables:      %d (%s)\n", d,
              paste(x$labels, collapse = ", ")))
  cat(sprintf("  Subjects:       %d\n", x$n_subjects))
  cat(sprintf("  Temporal edges: %d / %d\n", n_temp, d * d))
  cat(sprintf("  Contemp edges:  %d / %d\n", n_contemp, d * (d - 1) / 2))
  cat(sprintf("  Between edges:  %d / %d\n", n_between, d * (d - 1) / 2))
  if (!is.null(x$subject_PCC)) {
    cat(sprintf("  Subject nets:   %d\n", length(x$subject_PCC)))
  }
  cat(sprintf("  EBIC gamma:     %.2f\n", x$gamma))
  invisible(x)
}


#' @export
summary.ml_graphical_var_result <- function(object, ...) {
  cat("=== Temporal Network (beta) ===\n")
  print(round(object$temporal, 4))
  cat("\n=== Contemporaneous Network (PCC) ===\n")
  print(round(object$PCC, 4))
  cat("\n=== Between-Subjects Network ===\n")
  print(round(object$between, 4))
  invisible(object)
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
  pdc <- t(beta / sqrt(sigma_diag %o% kappa_diag + beta^2))
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
