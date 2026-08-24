# ---- Group Iterative Multiple Model Estimation (GIMME) ----

#' GIMME: Group Iterative Multiple Model Estimation
#'
#' @description
#' Estimates person-specific directed networks from intensive longitudinal data
#' using the unified Structural Equation Modeling (uSEM) framework. Implements
#' a data-driven search that identifies:
#' \enumerate{
#'   \item \strong{Group-level paths}: Directed edges present for a majority
#'     (default 75\%) of individuals.
#'   \item \strong{Individual-level paths}: Additional edges specific to each
#'     person, found after group paths are established.
#' }
#'
#' Estimation is delegated to [idiographic::fit_gimme()], the clean-room home
#' of the temporal idiographic estimators, whose search reproduces the
#' upstream \pkg{gimme} package (>= 10.0) exactly (verified at tolerance 0 on
#' path counts and per-person coefficient matrices in idiographic's own
#' parity suite). Uses \code{lavaan} for SEM estimation and modification
#' indices. Accepts a single data frame with an ID column (not CSV
#' directories).
#'
#' @param data A \code{data.frame} in long format with columns for person ID,
#'   time-varying variables, and optionally a time/beep column.
#' @param vars Character vector of variable names to model.
#' @param id Character string naming the person-ID column.
#' @param time Character string naming the time/order column, or \code{NULL}.
#'   When provided, data is sorted by \code{id} then \code{time} before lagging.
#' @param ar Logical. If \code{TRUE} (default), autoregressive paths
#'   (each variable predicting itself at lag 1) are included as fixed paths.
#' @param standardize Logical. If \code{TRUE} (default \code{FALSE}),
#'   variables are standardized per person before estimation.
#' @param groupcutoff Numeric between 0 and 1. Proportion of individuals for
#'   whom a path must be significant to be added at group level.
#'   Default \code{0.75}.
#' @param subcutoff Numeric. Not used (reserved for future subgrouping).
#' @param paths Character vector of lavaan-syntax paths to force into the model
#'   (e.g., \code{"V2~V1lag"}). Default \code{NULL}.
#' @param exogenous Character vector of variable names to treat as exogenous.
#'   Default \code{NULL}.
#' @param hybrid Logical. If \code{TRUE}, also searches residual covariances.
#'   Default \code{FALSE}.
#' @param rmsea_cutoff Numeric. RMSEA threshold for excellent fit (default 0.05).
#' @param srmr_cutoff Numeric. SRMR threshold for excellent fit (default 0.05).
#' @param nnfi_cutoff Numeric. NNFI/TLI threshold for excellent fit (default 0.95).
#' @param cfi_cutoff Numeric. CFI threshold for excellent fit (default 0.95).
#' @param n_excellent Integer. Number of fit indices that must be excellent to
#'   stop individual search. Default \code{2}.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#'
#' @return An S3 object of class \code{"net_gimme"} containing:
#' \describe{
#'   \item{\code{temporal}}{p x p matrix of group-level temporal (lagged)
#'     path counts -- entry \code{[i,j]} = number of individuals with path j(t-1)->i(t).}
#'   \item{\code{contemporaneous}}{p x p matrix of group-level contemporaneous
#'     path counts -- entry \code{[i,j]} = number of individuals with path j(t)->i(t).}
#'   \item{\code{temporal_avg}, \code{contemporaneous_avg}}{p x p group-average
#'     coefficient matrices.}
#'   \item{\code{coefs}}{List of per-person p x 2p coefficient matrices
#'     (rows = endogenous, cols = \code{[lagged, contemporaneous]}).}
#'   \item{\code{psi}}{List of per-person residual covariance matrices.}
#'   \item{\code{fit}}{Data frame of per-person fit indices (chisq, df, pvalue,
#'     rmsea, srmr, nnfi, cfi, bic, aic, logl, status).}
#'   \item{\code{path_counts}}{p x 2p matrix: how many individuals have each path.}
#'   \item{\code{paths}}{List of per-person character vectors of lavaan path syntax.}
#'   \item{\code{group_paths}}{Character vector of group-level paths found.}
#'   \item{\code{individual_paths}}{List of per-person character vectors of
#'     individual-level paths (beyond group).}
#'   \item{\code{syntax}}{List of per-person full lavaan syntax strings.}
#'   \item{\code{labels}}{Character vector of variable names.}
#'   \item{\code{n_subjects}}{Integer. Number of individuals.}
#'   \item{\code{n_obs}}{Integer vector. Time points per individual.}
#'   \item{\code{config}}{List of configuration parameters.}
#' }
#' The object additionally carries idiographic's netobject fields
#' (\code{weights}, \code{nodes}, \code{edges}, ...) so it renders directly
#' with cograph, and dispatches to idiographic's \code{print}, \code{summary},
#' and \code{plot} methods.
#'
#' @examplesIf requireNamespace("lavaan", quietly = TRUE)
#' \donttest{
#' # Create simple panel data (3 subjects, 4 variables, 50 time points).
#' set.seed(42)
#' n_sub <- 3; n_t <- 50; vars <- paste0("V", 1:4)
#' rows <- lapply(seq_len(n_sub), function(i) {
#'   d <- as.data.frame(matrix(rnorm(n_t * 4), ncol = 4))
#'   names(d) <- vars; d$id <- i; d
#' })
#' panel <- do.call(rbind, rows)
#' res <- build_gimme(panel, vars = vars, id = "id")
#' print(res)
#' }
#'
#' @seealso \code{\link{build_network}}
#'
#' @export
build_gimme <- function(data,
                        vars,
                        id,
                        time = NULL,
                        ar = TRUE,
                        standardize = FALSE,
                        groupcutoff = 0.75,
                        subcutoff = 0.50,
                        paths = NULL,
                        exogenous = NULL,
                        hybrid = FALSE,
                        rmsea_cutoff = 0.05,
                        srmr_cutoff = 0.05,
                        nnfi_cutoff = 0.95,
                        cfi_cutoff = 0.95,
                        n_excellent = 2L,
                        seed = NULL) {

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for build_gimme(). ", # nocov
         "Install it with install.packages('lavaan').", call. = FALSE) # nocov
  }

  if (!is.null(seed)) set.seed(seed)

  # --- Input validation (Nestimate's error surface, unchanged) ---
  stopifnot(is.data.frame(data))
  stopifnot(is.character(vars), length(vars) >= 2L)
  stopifnot(is.character(id), length(id) == 1L, id %in% names(data))
  if (!all(vars %in% names(data))) {
    missing_v <- setdiff(vars, names(data))
    stop("Variables not found in data: ", paste(missing_v, collapse = ", "),
         call. = FALSE)
  }
  stopifnot(is.numeric(groupcutoff), groupcutoff > 0, groupcutoff <= 1)
  if (!is.null(exogenous)) {
    stopifnot(is.character(exogenous))
    bad_exog <- setdiff(exogenous, vars)
    if (length(bad_exog) > 0L) {
      stop("`exogenous` names must be among `vars`: ",
           paste(bad_exog, collapse = ", "), call. = FALSE)
    }
    if (length(setdiff(vars, exogenous)) < 1L) {
      stop("`exogenous` cannot include every variable -- ",
           "at least one endogenous variable is required.", call. = FALSE)
    }
  }
  if (length(unique(data[[id]])) < 2L) {
    stop("build_gimme() requires at least 2 individuals.", call. = FALSE)
  }

  # Estimation is owned by idiographic; its search is upstream-gimme-exact
  # (>= 10.0), which the former in-package search was not -- path selection
  # can therefore differ from Nestimate <= 0.8.5 on the same data.
  # `subcutoff` is not forwarded: it is a no-op reserved formal here and in
  # idiographic (subgrouping is not implemented in either).
  idiographic::fit_gimme(
    data, vars = vars, id = id, time = time,
    ar = ar, standardize = standardize,
    groupcutoff = groupcutoff,
    paths = paths, exogenous = exogenous, hybrid = hybrid,
    rmsea_cutoff = rmsea_cutoff, srmr_cutoff = srmr_cutoff,
    nnfi_cutoff = nnfi_cutoff, cfi_cutoff = cfi_cutoff,
    n_excellent = as.integer(n_excellent), seed = seed
  )
}
