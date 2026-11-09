# ---- Group Iterative Multiple Model Estimation (GIMME) ----

#' GIMME: Group Iterative Multiple Model Estimation
#'
#' @description
#' `build_gimme()` estimates a directed network for each person from
#' intensive longitudinal data with GIMME (group iterative multiple model
#' estimation; Gates & Molenaar, 2012), in the unified structural equation
#' model (uSEM): each variable at time \eqn{t} is regressed on the variables
#' at \eqn{t - 1} (temporal, lagged paths) and on the other variables at
#' \eqn{t} (contemporaneous paths). The search runs in two stages:
#' \enumerate{
#'   \item \strong{Group-level paths}: directed paths added for everyone
#'     because they would improve the model for a large share of the people
#'     (more than the share `groupcutoff`, by default `0.75`, of the people
#'     whose models converged).
#'   \item \strong{Individual-level paths}: further paths specific to each
#'     person, added after the group paths until that person's model fits
#'     well (see `n_excellent`).
#' }
#' It takes one long data frame with a person-ID column, not a directory of
#' CSV files.
#'
#' The temporal and contemporaneous count matrices say how many people have
#' each path, which shows how common a path is in the sample, not how strong
#' it is; the coefficients are estimated per person. Group paths are shared
#' by construction, while individual paths are what the search adds for one
#' person's series: before reading them, check each person's number of time
#' points and fit statistics, which the result carries.
#'
#' @details
#' Estimation is delegated to [idiographic::fit_gimme()], the clean-room home
#' of the temporal idiographic estimators, whose search reproduces the
#' upstream \pkg{gimme} package (>= 10.0) exactly (verified at tolerance 0 on
#' path counts and per-person coefficient matrices in idiographic's own
#' parity suite). The models are fitted with \pkg{lavaan}, which also gives
#' the modification indices the search uses; `build_gimme()` raises an error
#' when \pkg{lavaan} is not installed. The settings `build_gimme()` does not
#' expose (such as the multiple-comparison corrections and `alpha`) keep
#' [idiographic::fit_gimme()]'s defaults.
#'
#' @section Results changed in 0.9.0:
#' Before 0.9.0 the search ran in an in-package implementation that was
#' \emph{not} upstream-gimme-exact. Delegating to
#' [idiographic::fit_gimme()] changed which paths the search selects on the
#' same data -- individual-level paths in particular -- so numeric results
#' are not comparable with Nestimate <= 0.8.5. The returned object also
#' gained fields (see \strong{Value}); nothing was removed.
#'
#' @param data A data frame in long format: one row per person and time
#'   point, with the person-ID column, the variables and, optionally, a time
#'   or beep column.
#' @param vars The names of the variables to model, at least two, all
#'   columns of `data`; a missing name raises an error.
#' @param id The name of the person-ID column. The data must hold at least
#'   two people; a person with fewer than three time points is left out of
#'   the analysis without a warning.
#' @param time The name of the time or order column. When given, the rows
#'   are sorted by `id` and then `time` before the lagged variables are
#'   formed; with `NULL` (default), the row order is used.
#' @param ar Whether each variable's autoregressive path (itself at lag 1)
#'   is fixed into every model. Default `TRUE`.
#' @param standardize Whether the variables are standardised within each
#'   person before estimation. Default `FALSE`. The coefficients come from
#'   lavaan's standardised solution either way, so this changes only the
#'   scale-dependent fit statistics (the log-likelihood, AIC and BIC); the
#'   chi-square, RMSEA, SRMR, NNFI and CFI do not change.
#' @param groupcutoff A path is added at the group level when it is
#'   significant for more than this share of the people whose models
#'   converged, a number in (0, 1]. No group path is added when half or fewer
#'   of the models converge. Default `0.75`.
#' @param subcutoff Not used: reserved for subgrouping, which neither
#'   Nestimate nor [idiographic::fit_gimme()] implements, and accepted for
#'   compatibility with \pkg{gimme}; it is not passed on. Default `0.50`.
#' @param paths Paths, in lavaan syntax, to fix into every model, such as
#'   `"V2~V1lag"`. Default `NULL`.
#' @param exogenous The variables treated as exogenous (predictors only),
#'   all among `vars`; at least one variable must remain endogenous, and
#'   either rule broken raises an error. Default `NULL`.
#' @param hybrid Whether the search also considers residual covariances
#'   between the variables at the same time point. Default `FALSE`.
#' @param rmsea_cutoff The RMSEA at or below which a person's fit counts as
#'   excellent. Default `0.05`.
#' @param srmr_cutoff The SRMR at or below which a person's fit counts as
#'   excellent. Default `0.05`.
#' @param nnfi_cutoff The NNFI (TLI) at or above which a person's fit counts
#'   as excellent. Default `0.95`.
#' @param cfi_cutoff The CFI at or above which a person's fit counts as
#'   excellent. Default `0.95`.
#' @param n_excellent How many of the four fit indices must be excellent for
#'   a person's individual search to stop, a whole number from 1 to 4.
#'   Default `2`.
#' @param seed An integer that makes the estimation reproducible; the
#'   session's own random stream is left as it was. Default `NULL`.
#'
#' @return The `net_gimme` object (class
#'   `c("net_gimme", "cograph_network", "list")`) that
#'   [idiographic::fit_gimme()] returns; its `print()`, `summary()` and
#'   `plot()` methods are \pkg{idiographic}'s, not Nestimate's. `print()`
#'   shows the group-level paths, the number of individual paths and the
#'   share of people with each path; [summary()] returns a data frame with
#'   one row per network layer (temporal, contemporaneous) and its size,
#'   density and weights; `plot()` draws the network. [idiographic::edges()]
#'   gives one row per path with its level (group or individual), and
#'   [idiographic::coefs()] one row per person and path with its
#'   coefficient. The object is a \emph{superset} of the pre-0.9.0
#'   in-package field contract: every field below is present, alongside
#'   idiographic's additions (\code{contemp_cov}, \code{contemp_cov_avg},
#'   \code{contemp_is_cov}). Its fields:
#' \describe{
#'   \item{\code{temporal}}{p x p matrix of temporal (lagged) path counts:
#'     entry \code{[i, j]} is the number of people with the path
#'     j(t-1) -> i(t), group and individual paths alike.}
#'   \item{\code{contemporaneous}}{p x p matrix of contemporaneous path
#'     counts: entry \code{[i, j]} is the number of people with the path
#'     j(t) -> i(t).}
#'   \item{\code{temporal_avg}, \code{contemporaneous_avg}}{p x p
#'     group-average coefficient matrices.}
#'   \item{\code{coefs}}{List of per-person p x 2p coefficient matrices
#'     (rows the endogenous variables, columns the lagged and then the
#'     contemporaneous variables).}
#'   \item{\code{psi}}{List of per-person residual covariance matrices.}
#'   \item{\code{fit}}{Data frame of per-person fit indices (\code{file},
#'     the person, then chisq, df, pvalue, rmsea, srmr, nnfi, cfi, bic, aic,
#'     logl and status).}
#'   \item{\code{path_counts}}{p x 2p matrix: how many people have each
#'     path.}
#'   \item{\code{paths}}{List of per-person character vectors of lavaan path
#'     syntax.}
#'   \item{\code{group_paths}}{Character vector of the group-level paths
#'     found.}
#'   \item{\code{individual_paths}}{List of per-person character vectors of
#'     the individual-level paths (beyond the group paths).}
#'   \item{\code{syntax}}{List of per-person full lavaan syntax strings.}
#'   \item{\code{labels}}{Character vector of variable names.}
#'   \item{\code{n_subjects}}{Number of people.}
#'   \item{\code{n_obs}}{Integer vector of the lagged time points per
#'     person (its rows minus one).}
#'   \item{\code{config}}{List of the settings used, including the
#'     corrections and cutoffs that are not arguments here.}
#' }
#' The object also carries idiographic's network fields (\code{nodes},
#' \code{edges}, \code{directed}, \code{data}, \code{meta},
#' \code{node_groups}, and a \code{weights} entry that is NULL) for
#' rendering with cograph.
#'
#' @examplesIf requireNamespace("lavaan", quietly = TRUE)
#' \donttest{
#' # Create simple panel data (3 subjects, 4 variables, 30 time points).
#' set.seed(42)
#' n_sub <- 3; n_t <- 30; vars <- paste0("V", 1:4)
#' rows <- lapply(seq_len(n_sub), function(i) {
#'   d <- as.data.frame(matrix(rnorm(n_t * 4), ncol = 4))
#'   names(d) <- vars; d$id <- i; d
#' })
#' panel <- do.call(rbind, rows)
#' res <- build_gimme(panel, vars = vars, id = "id")
#' print(res)
#' }
#'
#' @references
#' Gates, K. M., & Molenaar, P. C. M. (2012). Group search algorithm
#' recovers effective connectivity maps for individuals in homogeneous and
#' heterogeneous samples. \emph{NeuroImage}, 63(1), 310--319.
#' \doi{10.1016/j.neuroimage.2012.06.026}
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

  if (!is.null(seed)) .local_seed(seed)

  # --- Input validation (Nestimate's error surface, unchanged) ---
  .assert(is.data.frame(data))
  .assert(is.character(vars), length(vars) >= 2L)
  .assert(is.character(id), length(id) == 1L, id %in% names(data))
  if (!all(vars %in% names(data))) {
    missing_v <- setdiff(vars, names(data))
    stop(errorCondition(
      paste0(
        "Variables not found in data: ",
        paste(missing_v, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  .assert(is.numeric(groupcutoff), groupcutoff > 0, groupcutoff <= 1)
  if (!is.null(exogenous)) {
    .assert(is.character(exogenous))
    bad_exog <- setdiff(exogenous, vars)
    if (length(bad_exog) > 0L) {
      stop(errorCondition(paste0("`exogenous` names must be among `vars`: ",
           paste(bad_exog, collapse = ", ")),
                          class = "nestimate_bad_input", call = NULL))
    }
    if (length(setdiff(vars, exogenous)) < 1L) {
      stop(errorCondition(
        paste0(
          "`exogenous` cannot include every variable -- ",
          "at least one endogenous variable is required."),
        class = "nestimate_bad_input", call = NULL))
    }
  }
  if (length(unique(data[[id]])) < 2L) {
    stop(errorCondition("build_gimme() requires at least 2 individuals.",
                        class = "nestimate_bad_input", call = NULL))
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
