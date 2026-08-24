# ---- Multilevel Vector Autoregression (mlVAR) ----

#' Build a Multilevel Vector Autoregression (mlVAR) network
#'
#' @description Estimates three networks from ESM/EMA panel data, matching
#'   `mlVAR::mlVAR()` with `estimator = "lmer"`, `temporal = "fixed"`,
#'   `contemporaneous = "fixed"` at machine precision: (1) a directed
#'   temporal network of fixed-effect lagged regression coefficients, (2)
#'   an undirected contemporaneous network of partial correlations among
#'   residuals, and (3) an undirected between-subjects network of partial
#'   correlations derived from the person-mean fixed effects.
#'
#' @details Estimation is delegated to [idiographic::fit_mlvar()] (the
#' clean-room home of the temporal idiographic estimators), called with
#' `estimator = "lmer"`, `temporal = "fixed"`,
#' `contemporaneous = "fixed"`. The pipeline follows mlVAR's lmer path
#' exactly:
#' \enumerate{
#'   \item Drop rows with NA in id/day/beep and optionally grand-mean
#'         standardize each variable.
#'   \item Expand the per-(id, day) beep grid and right-join original
#'         values, producing the augmented panel (`augData`).
#'   \item Add within-person lagged predictors (`L1_*`) and person-mean
#'         predictors (`PM_*`).
#'   \item For each outcome variable fit
#'         `lmer(y ~ within + between-except-own-PM + (1 | id))` with
#'         `REML = FALSE`. Collect the fixed-effect temporal matrix `B`,
#'         between-effect matrix `Gamma`, random-intercept SDs (`mu_SD`),
#'         and lmer residual SDs.
#'   \item Contemporaneous network:
#'         `cor2pcor(D %*% cov2cor(cor(resid)) %*% D)`.
#'   \item Between-subjects network:
#'         `cor2pcor(pseudoinverse(forcePositive(D (I - Gamma))))`.
#' }
#'
#' Validated to machine precision (max_diff < 1e-10) against
#' `mlVAR::mlVAR()` on 25 real ESM datasets from `openesm` and 20 simulated
#' configurations, and to exact equality (max_diff == 0) against the
#' pre-delegation Nestimate implementation on all layers, coefficients,
#' and observation counts across lag/standardize/day/beep configurations.
#'
#' @param data A `data.frame` containing the panel data.
#' @param vars Character vector of variable column names to model.
#' @param id Character string naming the person-ID column.
#' @param day Character string naming the day/session column, or `NULL`.
#'   When provided, lag pairs are only formed within the same day.
#' @param beep Character string naming the measurement-occasion column, or
#'   `NULL`. When `NULL`, row position within each (id, day) is used.
#' @param lag Integer. The lag order (default 1).
#' @param standardize Logical. If `TRUE`, each variable is grand-mean
#'   centered and divided by its pooled SD *before* augmentation. Default
#'   `FALSE`, matching `mlVAR::mlVAR(scale = FALSE)` - the only setting for
#'   which numerical equivalence has been validated.
#'
#' @return A dual-class `c("net_mlvar", "netobject_group")` object - a
#'   named list of three full netobjects, one per network, plus
#'   model-level metadata stored as attributes. Each element is a
#'   standard `c("netobject", "cograph_network")` weight-matrix wrapper
#'   (no raw `$data`), so `print()`, `summary()`, [coefs()], and
#'   `cograph::splot(fit$temporal)` work directly. The three constituents
#'   are matrix-wrapped and carry no underlying panel data, so
#'   data-resampling verbs such as [bootstrap_network()] (and
#'   reliability/stability) cannot iterate over them - extract a single
#'   constituent and rebuild via [build_network()] if you need those.
#'   Structure:
#'   \describe{
#'     \item{`fit$temporal`}{Directed netobject for the `d x d` matrix of
#'       fixed-effect lagged coefficients. `$weights[i, j]` is the effect
#'       of variable j at t-lag on variable i at t. `method =
#'       "mlvar_temporal"`, `directed = TRUE`.}
#'     \item{`fit$contemporaneous`}{Undirected netobject for the `d x d`
#'       partial-correlation network of within-person lmer residuals.
#'       `method = "mlvar_contemporaneous"`, `directed = FALSE`.}
#'     \item{`fit$between`}{Undirected netobject for the `d x d`
#'       partial-correlation network of person means, derived from
#'       `D (I - Gamma)`. `method = "mlvar_between"`, `directed = FALSE`.}
#'     \item{`attr(fit, "coefs")` / [coefs()]}{Tidy `data.frame` with one
#'       row per `(outcome, predictor)` pair and columns `outcome`,
#'       `predictor`, `beta`, `se`, `t`, `p`, `ci_lower`, `ci_upper`,
#'       `significant`. Filter, sort, or plot with base R or the tidyverse.
#'       Retrieve with `coefs(fit)`.}
#'     \item{`attr(fit, "n_obs")`}{Number of rows in the augmented panel
#'       after na.omit.}
#'     \item{`attr(fit, "n_subjects")`}{Number of unique subjects remaining.}
#'     \item{`attr(fit, "lag")`}{Lag order used.}
#'     \item{`attr(fit, "standardize")`}{Logical; whether pre-augmentation
#'       standardization was applied.}
#'   }
#'
#' @examples
#' \dontrun{
#' d <- simulate_data("mlvar", seed = 1)
#' fit <- build_mlvar(d, vars = attr(d, "vars"),
#'                    id = "id", day = "day", beep = "beep")
#' print(fit)
#' summary(fit)
#' }
#'
#' @seealso [build_network()]
#' @export
build_mlvar <- function(data, vars, id,
                        day = NULL, beep = NULL,
                        lag = 1L,
                        standardize = FALSE) {
  # ---- Input validation ----
  stopifnot(
    is.data.frame(data),
    is.character(vars), length(vars) >= 2L,
    is.character(id), length(id) == 1L,
    is.numeric(lag), length(lag) == 1L, lag >= 1L,
    is.logical(standardize), length(standardize) == 1L
  )

  required <- c(vars, id)
  if (!is.null(day))  required <- c(required, day)
  if (!is.null(beep)) required <- c(required, beep)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop("Columns not found in data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  lag <- as.integer(lag)

  # Estimation is owned by idiographic; verified exactly equal (max abs
  # delta == 0 on all three layers, coefs, and n_obs) to the former
  # in-package lmer pipeline across lag/standardize/day/beep configs.
  fit <- idiographic::fit_mlvar(
    data, vars = vars, id = id, day = day, beep = beep,
    lags = lag,
    estimator = "lmer", temporal = "fixed", contemporaneous = "fixed",
    scale = standardize
  )

  stopifnot(
    "idiographic::fit_mlvar() did not return the three network layers" =
      all(c("temporal", "contemporaneous", "between") %in% names(fit)),
    "idiographic::fit_mlvar() result carries no coefficient table" =
      is.data.frame(attr(fit, "coefs"))
  )

  # Re-wrap each layer through Nestimate's own constructor so the return
  # contract (fields, meta, class) is byte-identical to the pre-delegation
  # object regardless of idiographic's internal netobject flavour.
  nets <- list(
    temporal        = .wrap_netobject(fit$temporal$weights,
                                      method   = "mlvar_temporal",
                                      directed = TRUE),
    contemporaneous = .wrap_netobject(fit$contemporaneous$weights,
                                      method   = "mlvar_contemporaneous",
                                      directed = FALSE),
    between         = .wrap_netobject(fit$between$weights,
                                      method   = "mlvar_between",
                                      directed = FALSE)
  )

  # Model-level metadata lives in attributes so the list stays a pure
  # netobject_group (each element is a netobject). Use coefs(fit) to
  # retrieve the tidy coefs data.frame.
  attr(nets, "coefs")       <- attr(fit, "coefs")
  attr(nets, "n_obs")       <- attr(fit, "n_obs")
  attr(nets, "n_subjects")  <- attr(fit, "n_subjects")
  attr(nets, "lag")         <- lag
  attr(nets, "standardize") <- standardize
  attr(nets, "group_col")   <- "network_type"

  class(nets) <- c("net_mlvar", "netobject_group")
  nets
}

#' Tidy coefficients from a fitted mlvar model
#'
#' Generic accessor for the tidy coefficient table stored on a
#' [build_mlvar()] result. Returns a `data.frame` with one row per
#' `(outcome, predictor)` pair and columns `outcome`, `predictor`,
#' `beta`, `se`, `t`, `p`, `ci_lower`, `ci_upper`, `significant`.
#'
#' Only the within-person (temporal) coefficients are tabulated -
#' these are the lagged fixed effects that populate `fit$temporal`.
#' The between-subjects effects that go into `fit$between` are handled
#' via the `D (I - Gamma)` transformation and are not exposed as a
#' separate tidy table.
#'
#' @param x A fitted model object - currently only `net_mlvar` is supported.
#' @param ... Unused.
#' @return A tidy `data.frame` of coefficient estimates.
#' @inherit build_mlvar examples
#' @export
coefs <- function(x, ...) {
  UseMethod("coefs")
}

#' @rdname coefs
#' @export
coefs.net_mlvar <- function(x, ...) {
  attr(x, "coefs")
}

#' @rdname coefs
#' @export
coefs.default <- function(x, ...) {
  stop("No coefs() method for object of class '",
       class(x)[1], "'", call. = FALSE)
}

# ---- S3 methods --------------------------------------------------------

#' Print method for net_mlvar
#'
#' @param x A `net_mlvar` object returned by [build_mlvar()].
#' @param ... Unused; present for S3 consistency.
#' @return Invisibly returns `x`.
#' @inherit build_mlvar examples
#' @export
print.net_mlvar <- function(x, ...) {
  coef_df <- attr(x, "coefs")
  d <- nrow(x$temporal$weights)
  n_sig <- sum(coef_df$significant, na.rm = TRUE)
  n_tot <- nrow(coef_df)

  cat(sprintf("mlVAR result: %d subjects, %d observations, %d variables (lag %d)\n",
              attr(x, "n_subjects"),
              attr(x, "n_obs"),
              d,
              attr(x, "lag")))
  cat(sprintf("  $temporal        %d x %d directed    (%d/%d edges significant at p<0.05)\n",
              d, d, n_sig, n_tot))
  cat(sprintf("  $contemporaneous %d x %d undirected\n", d, d))
  cat(sprintf("  $between         %d x %d undirected\n", d, d))
  cat("Access: fit$temporal, fit$contemporaneous, fit$between (each a netobject);",
      "coefs(fit) for tidy coefs.\n")
  invisible(x)
}

#' Summary method for net_mlvar
#'
#' @param object A `net_mlvar` object returned by [build_mlvar()].
#' @param ... Unused; present for S3 consistency.
#' @return Invisibly returns `object`.
#' @inherit build_mlvar examples
#' @export
summary.net_mlvar <- function(object, ...) {
  coef_df <- attr(object, "coefs")
  vars    <- rownames(object$temporal$weights)
  d       <- length(vars)

  cat("=== mlVAR Summary ===\n")
  cat("Subjects:",      attr(object, "n_subjects"),
      " | Observations:", attr(object, "n_obs"),
      " | Variables:",  d,
      " | Lag:",        attr(object, "lag"), "\n")
  cat("Standardized:",  attr(object, "standardize"), "\n\n")

  cat("--- Temporal Network (B matrix) ---\n")
  print(round(object$temporal$weights, 4))
  cat("\n")

  sig_rows <- coef_df[coef_df$significant, , drop = FALSE]
  if (nrow(sig_rows) > 0L) {
    cat("Significant temporal edges (p < 0.05):\n")
    sig_print <- sig_rows[, c("outcome", "predictor",
                              "beta", "se", "t", "p"), drop = FALSE]
    sig_print$beta <- round(sig_print$beta, 4)
    sig_print$se   <- round(sig_print$se, 4)
    sig_print$t    <- round(sig_print$t, 3)
    sig_print$p    <- round(sig_print$p, 4)
    rownames(sig_print) <- NULL
    print(sig_print)
  } else {
    cat("No significant temporal edges at p < 0.05.\n")
  }
  cat("\n")

  cat("--- Contemporaneous Network ---\n")
  print(round(object$contemporaneous$weights, 4))
  cat("\n")

  cat("--- Between-Subjects Network ---\n")
  print(round(object$between$weights, 4))
  cat("\n")

  coef_df
}
