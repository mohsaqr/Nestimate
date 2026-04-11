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
#'   `mlvar()` is a short-cut alias for `build_mlvar()`.
#'
#' @details The algorithm follows mlVAR's lmer pipeline exactly:
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
#' configurations (seeds 201-220). See `tmp/mlvar_equivalence_real20.R` and
#' `tmp/mlvar_equivalence_20seeds.R`.
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
#'   `FALSE`, matching `mlVAR::mlVAR(scale = FALSE)` — the only setting for
#'   which numerical equivalence has been validated.
#'
#' @return A dual-class `c("net_mlvar", "cograph_network")` object — a list
#'   with both the three mlvar networks and the `cograph_network` fields
#'   (`$weights`, `$nodes`, `$edges`, `$directed`, `$meta`). The default
#'   cograph view is the temporal network, so calling `cograph::splot(fit)`
#'   (if cograph is installed) plots the temporal matrix directly. cograph
#'   handles all plot dispatch; Nestimate does not provide a plot method.
#'   Components:
#'   \describe{
#'     \item{`temporal`}{`d x d` matrix of fixed-effect temporal
#'       coefficients. Entry `[i, j]` is the effect of variable j at
#'       t-lag on variable i at t.}
#'     \item{`contemporaneous`}{`d x d` symmetric matrix of partial
#'       correlations among within-person lmer residuals.}
#'     \item{`between`}{`d x d` symmetric matrix of partial correlations
#'       among person means, derived from `D (I - Gamma)`.}
#'     \item{`coefs`}{Tidy `data.frame` with one row per `(outcome,
#'       predictor)` pair and columns `outcome`, `predictor`, `beta`,
#'       `se`, `t`, `p`, `ci_lower`, `ci_upper`, `significant`. Filter,
#'       sort, or plot with base R or the tidyverse.}
#'     \item{`n_obs`}{Number of rows in the augmented panel after na.omit.}
#'     \item{`n_subjects`}{Number of unique subjects remaining.}
#'     \item{`lag`}{Lag order used.}
#'     \item{`standardize`}{Logical; whether pre-augmentation standardization
#'       was applied.}
#'     \item{`weights`, `nodes`, `edges`, `directed`, `meta`}{`cograph_network`
#'       fields, populated with the default temporal view.}
#'   }
#'
#' @examples
#' \dontrun{
#' d <- simulate_data("mlvar", seed = 1)
#' fit <- build_mlvar(d, vars = attr(d, "vars"),
#'                    id = "id", day = "day", beep = "beep")
#' print(fit)
#' summary(fit)
#'
#' # Short-cut alias
#' fit2 <- mlvar(d, vars = attr(d, "vars"),
#'               id = "id", day = "day", beep = "beep")
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

  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for build_mlvar().", call. = FALSE)
  }
  if (!requireNamespace("corpcor", quietly = TRUE)) {
    stop("Package 'corpcor' is required for build_mlvar().", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for build_mlvar().", call. = FALSE)
  }

  lag <- as.integer(lag)

  prepared <- .mlvar_prepare_data(data, vars, id, day, beep, standardize)
  aug      <- .mlvar_augment_data(prepared, vars, id, day, beep, lag)
  Res      <- .mlvar_estimate_lmer(aug$data, aug$predModel, vars, id)

  # Attach the TEMPORAL network as the default cograph_network view so
  # that `cograph::splot(net)` plots it directly. Nestimate never calls
  # cograph; plot dispatch is cograph's responsibility. See cograph's
  # CLAUDE.md "Nestimate net_mlvar" section for the spec of how cograph
  # should handle the three networks and per-type styling.
  cg <- .mlvar_cograph_fields(
    mat = Res$temporal$B, vars = vars,
    directed = TRUE, type = "temporal"
  )

  result <- list(
    temporal        = Res$temporal$B,
    contemporaneous = Res$contemporaneous,
    between         = Res$between,
    coefs           = Res$temporal$coefs,
    n_obs           = nrow(aug$data),
    n_subjects      = length(unique(aug$data[[id]])),
    lag             = lag,
    standardize     = standardize,
    # cograph_network fields (default view = temporal)
    weights         = cg$weights,
    nodes           = cg$nodes,
    edges           = cg$edges,
    directed        = cg$directed,
    n_nodes         = cg$n_nodes,
    n_edges         = cg$n_edges,
    meta            = cg$meta,
    node_groups     = NULL
  )
  class(result) <- c("net_mlvar", "cograph_network")
  result
}

#' @rdname build_mlvar
#' @export
mlvar <- build_mlvar

#' Build cograph_network fields from a single mlvar matrix
#'
#' Lightweight helper that packages one of the three mlvar matrices
#' (temporal / contemporaneous / between) as a cograph-compatible view
#' (`$weights`, `$nodes`, `$edges`, `$directed`, `$meta`). Used both to
#' attach the default temporal view at build time and to swap views inside
#' `plot.net_mlvar()`.
#' @noRd
.mlvar_cograph_fields <- function(mat, vars, directed, type) {
  nodes_df <- data.frame(
    id    = seq_along(vars),
    label = vars,
    name  = vars,
    stringsAsFactors = FALSE
  )
  edges <- .extract_edges_from_matrix(mat, directed = directed)
  list(
    weights  = mat,
    nodes    = nodes_df,
    edges    = edges,
    directed = directed,
    n_nodes  = length(vars),
    n_edges  = nrow(edges),
    meta = list(
      source = "nestimate",
      layout = NULL,
      tna    = list(method = paste0("mlvar_", type))
    )
  )
}


# ---- Internal helpers --------------------------------------------------

#' Drop rows with NA metadata and optionally grand-mean standardize
#' @noRd
.mlvar_prepare_data <- function(data, vars, id, day, beep, scale) {
  df <- as.data.frame(data)

  md_cols <- c(id,
               if (!is.null(day))  day,
               if (!is.null(beep)) beep)
  df <- df[stats::complete.cases(df[, md_cols, drop = FALSE]), , drop = FALSE]

  if (isTRUE(scale)) {
    for (v in vars) {
      x <- as.numeric(df[[v]])
      sd_val <- stats::sd(x, na.rm = TRUE)
      if (is.na(sd_val) || sd_val == 0) {
        df[[v]] <- 0
      } else {
        df[[v]] <- (x - mean(x, na.rm = TRUE)) / sd_val
      }
    }
  }
  df
}

#' Beep-grid augmentation + within/between predictor construction
#'
#' Hybrid implementation: `data.table::CJ` + keyed join for the expensive
#' grid construction, base R `ave()` for the within-group lag/center/mean
#' arithmetic. The split matters because `data.table`'s optimized group
#' aggregators (`gmean` et al.) accumulate sums in a different order than
#' base R, which drifts at 1e-16 and amplifies through lmer into 1e-10
#' coefficient diffs against mlVAR.
#' @noRd
.mlvar_augment_data <- function(data, vars, id, day, beep, lag) {
  # Silence R CMD check NOTE for data.table NSE symbols
  .first <- .last <- NULL

  id_col   <- id
  day_col  <- if (is.null(day))  ".day"  else day
  beep_col <- if (is.null(beep)) ".beep" else beep

  dt <- data.table::as.data.table(data)
  if (is.null(day)) {
    dt[, (day_col) := 1L]
  }
  if (is.null(beep)) {
    dt[, (beep_col) := seq_len(.N), by = c(id_col, day_col)]
  }

  # Per-(id, day) beep range
  first_last <- dt[, .(.first = min(get(beep_col), na.rm = TRUE),
                       .last  = max(get(beep_col), na.rm = TRUE)),
                   by = c(id_col, day_col)]

  # Global beep grid, filtered down to each (id, day)'s actual range
  gb_min <- min(dt[[beep_col]], na.rm = TRUE)
  gb_max <- max(dt[[beep_col]], na.rm = TRUE)
  allBeeps <- data.table::CJ(
    V1 = unique(dt[[id_col]]),
    V2 = unique(dt[[day_col]]),
    V3 = seq(gb_min, gb_max),
    sorted = FALSE
  )
  data.table::setnames(allBeeps, c("V1", "V2", "V3"),
                       c(id_col, day_col, beep_col))
  allBeeps <- allBeeps[first_last, on = c(id_col, day_col), nomatch = NULL]
  allBeeps <- allBeeps[get(beep_col) >= .first & get(beep_col) <= .last]
  allBeeps[, c(".first", ".last") := NULL]

  # Right-join original data onto the grid (keeps all grid rows)
  data.table::setkeyv(dt, c(id_col, day_col, beep_col))
  data.table::setkeyv(allBeeps, c(id_col, day_col, beep_col))
  augData <- dt[allBeeps, on = c(id_col, day_col, beep_col),
                allow.cartesian = TRUE]
  data.table::setkeyv(augData, c(id_col, day_col, beep_col))

  # Drop back to a plain data.frame so subsequent `ave()` calls match
  # mlVAR's accumulation order bit-for-bit.
  augData <- as.data.frame(augData)
  rownames(augData) <- NULL

  predModel <- list()

  # Within (lagged, person-centered) predictors
  for (v in vars) {
    p_id <- paste0("L", lag, "_", v)
    augData[[p_id]] <- stats::ave(
      augData[[v]], augData[[id_col]], augData[[day_col]],
      FUN = function(x) .mlvar_aveLag(x, lag)
    )
    augData[[p_id]] <- stats::ave(
      augData[[p_id]], augData[[id_col]],
      FUN = function(x) x - mean(x, na.rm = TRUE)
    )
    predModel[[length(predModel) + 1L]] <- list(
      dep = vars, pred = v, id = p_id, type = "within"
    )
  }

  # Between (person-mean) predictors
  for (v in vars) {
    p_id <- paste0("PM_", v)
    augData[[p_id]] <- stats::ave(
      augData[[v]], augData[[id_col]],
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    predModel[[length(predModel) + 1L]] <- list(
      dep = vars, pred = v, id = p_id, type = "between"
    )
  }

  involved <- unique(c(vars, vapply(predModel, `[[`, character(1), "id")))
  augData <- stats::na.omit(
    augData[, c(involved, id_col, day_col, beep_col), drop = FALSE]
  )
  rownames(augData) <- NULL

  list(data = augData, predModel = predModel)
}

#' Fit d outcome-specific lmer models and assemble the three networks
#'
#' Matches `mlVAR:::lmer_mlVAR` with `temporal = "fixed"`,
#' `contemporaneous = "fixed"`. For each outcome k fits
#' `outcome_k ~ L1_v1 + ... + L1_vd + PM_v_{-k} + (1 | id)` with
#' `REML = FALSE`, then assembles Beta, Gamma, mu_SD, residuals.
#' @noRd
.mlvar_estimate_lmer <- function(augData, predModel, vars, id) {
  d <- length(vars)
  n_obs <- nrow(augData)

  B             <- matrix(0, d, d, dimnames = list(vars, vars))
  Gamma         <- matrix(0, d, d, dimnames = list(vars, vars))
  mu_SD         <- stats::setNames(numeric(d), vars)
  sigma_vec     <- stats::setNames(numeric(d), vars)
  residuals_mat <- matrix(NA_real_, n_obs, d, dimnames = list(NULL, vars))

  within_ids <- vapply(Filter(function(m) m$type == "within",  predModel),
                       `[[`, character(1), "id")
  between_ids <- vapply(Filter(function(m) m$type == "between", predModel),
                        `[[`, character(1), "id")
  var_to_within  <- stats::setNames(within_ids,  vars)
  var_to_between <- stats::setNames(between_ids, vars)

  z975 <- stats::qnorm(0.975)

  # Tidy coefs: one row per (outcome, predictor) pair — fills d * d rows.
  # Faster and cleaner than growing a list of per-outcome data.frames and
  # `do.call(rbind, ...)` at the end.
  n_coef_rows <- d * d
  coefs_tidy <- data.frame(
    outcome     = rep(vars, each = d),
    predictor   = rep(vars, times = d),
    beta        = numeric(n_coef_rows),
    se          = numeric(n_coef_rows),
    t           = numeric(n_coef_rows),
    p           = numeric(n_coef_rows),
    ci_lower    = numeric(n_coef_rows),
    ci_upper    = numeric(n_coef_rows),
    significant = logical(n_coef_rows),
    stringsAsFactors = FALSE
  )

  for (k in seq_len(d)) {
    outcome <- vars[k]
    # Own PM excluded — matches mlVAR's `getModel` filter on `dep == outcome`
    fixed_preds <- c(within_ids, var_to_between[-k])

    # Random intercept only — matches mlVAR temporal="fixed".
    fm_str <- paste0(outcome, " ~ ",
                     paste(fixed_preds, collapse = " + "),
                     " + (1 | ", id, ")")
    fit <- suppressMessages(suppressWarnings(
      lme4::lmer(stats::as.formula(fm_str), data = augData, REML = FALSE)
    ))

    fe <- lme4::fixef(fit)
    B[k, ]        <- fe[var_to_within[vars]]
    Gamma[k, -k]  <- fe[var_to_between[vars[-k]]]

    vc <- lme4::VarCorr(fit)
    ri_var <- as.numeric(vc[[id]][1, 1])
    mu_SD[k] <- if (!is.na(ri_var) && ri_var > 0) sqrt(ri_var) else 0

    sigma_vec[k] <- stats::sigma(fit)

    # Align residuals to augData row order (lmer drops any NA rows)
    res <- stats::residuals(fit)
    row_names <- rownames(augData)
    if (!is.null(row_names) && !is.null(names(res))) {
      residuals_mat[, k] <- res[match(row_names, names(res))]
    } else {
      residuals_mat[, k] <- res
    }

    sfe    <- summary(fit)$coefficients
    beta_k <- as.numeric(fe[var_to_within[vars]])
    se_k   <- as.numeric(sfe[var_to_within[vars], "Std. Error"])
    t_k    <- as.numeric(sfe[var_to_within[vars], "t value"])
    p_k    <- 2 * (1 - stats::pnorm(abs(t_k)))

    rows <- ((k - 1L) * d + 1L):(k * d)
    coefs_tidy$beta[rows]        <- beta_k
    coefs_tidy$se[rows]          <- se_k
    coefs_tidy$t[rows]           <- t_k
    coefs_tidy$p[rows]           <- p_k
    coefs_tidy$ci_lower[rows]    <- beta_k - z975 * se_k
    coefs_tidy$ci_upper[rows]    <- beta_k + z975 * se_k
    coefs_tidy$significant[rows] <- p_k < 0.05
  }

  contemporaneous <- .mlvar_contemporaneous_fixed(residuals_mat, sigma_vec, vars)
  between         <- .mlvar_compute_between_from_gamma(Gamma, mu_SD, vars)

  list(temporal = list(B = B, coefs = coefs_tidy, residuals = residuals_mat),
       contemporaneous = contemporaneous,
       between = between)
}

#' Within-group lag (matches mlVAR:::aveLag)
#'
#' Uses a logical `NA` (not `NA_real_`) for the prepended entries so that
#' integer input columns retain integer type. Preserving the integer type
#' is critical because base R's `mean()` uses a two-pass summation
#' correction for numeric input but a simple sum/n for integer input — the
#' two paths drift by ~1.4e-14, which then amplifies through lmer into
#' ~1e-10 coefficient diffs against mlVAR's integer-typed pipeline.
#' @noRd
.mlvar_aveLag <- function(x, lag = 1L) {
  n <- length(x)
  if (lag >= n) return(rep(NA, n))
  c(rep(NA, lag), x[seq_len(n - lag)])
}

#' Force a symmetric matrix to be positive-definite — byte-for-byte replica
#' of `mlVAR:::forcePositive`.
#'
#' Note the scalar-recycling quirk in the upstream implementation. In
#' `x - (diag(n) * min_ev - 0.001)`, the `0.001` scalar is subtracted from
#' every element of the diagonal matrix — so the final operation adds
#' `|min_ev|` to the diagonal *and* `+0.001` to every off-diagonal element.
#' This looks unintentional upstream but has to be replicated exactly for
#' equivalence with `mlVAR::mlVAR()`.
#' @noRd
.mlvar_force_positive <- function(x) {
  x <- (x + t(x)) / 2
  ev <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev < 0)) {
    x - (diag(nrow(x)) * min(ev) - 0.001)
  } else {
    x
  }
}

#' Between-subjects partial correlation from Gamma + mu_SD
#'
#' Matches mlVAR's Omega_mu branch:
#'   `D = diag(1 / mu_SD^2)`
#'   `inv = forcePositive(D (I - Gamma))`
#'   `cov = corpcor::pseudoinverse(inv)`
#'   `pcor = corpcor::cor2pcor(cov)`
#' @noRd
.mlvar_compute_between_from_gamma <- function(Gamma, mu_SD, vars) {
  d <- length(vars)
  if (any(mu_SD == 0)) {
    return(matrix(0, d, d, dimnames = list(vars, vars)))
  }

  D <- diag(1 / mu_SD^2)
  inv <- D %*% (diag(d) - Gamma)
  inv <- (inv + t(inv)) / 2
  inv <- .mlvar_force_positive(inv)

  mu_cov <- corpcor::pseudoinverse(inv)
  pcor <- corpcor::cor2pcor(mu_cov)
  diag(pcor) <- 0
  rownames(pcor) <- colnames(pcor) <- vars
  pcor
}

#' Contemporaneous partial correlation via mlVAR's "fixed" path
#'
#' Replicates `mlVAR:::lmer_mlVAR` Theta assembly for
#' `contemporaneous = "fixed"`: rescale the residual correlation by the
#' per-outcome lmer residual SDs and take `cor2pcor` directly. No
#' EBIC-GLASSO regularization. Note `cor2pcor` is scale-invariant, so the
#' `D %*% . %*% D` rescaling does not affect the pcor output — it is kept
#' only for parity with mlVAR's `cov`/`prec` slots.
#' @noRd
.mlvar_contemporaneous_fixed <- function(residuals_mat, sigma_vec, vars) {
  d <- length(vars)
  R <- stats::cor(residuals_mat, use = "pairwise.complete.obs")
  if (any(is.na(R))) {
    return(matrix(0, d, d, dimnames = list(vars, vars)))
  }
  D <- diag(sigma_vec)
  Theta_cov <- D %*% stats::cov2cor(R) %*% D
  pcor <- corpcor::cor2pcor(Theta_cov)
  diag(pcor) <- 0
  rownames(pcor) <- colnames(pcor) <- vars
  pcor
}


# ---- S3 methods --------------------------------------------------------

#' Print method for net_mlvar
#'
#' @param x A `net_mlvar` object returned by [build_mlvar()].
#' @param ... Unused; present for S3 consistency.
#' @return Invisibly returns `x`.
#' @export
print.net_mlvar <- function(x, ...) {
  d <- nrow(x$temporal)
  n_sig <- sum(x$coefs$significant, na.rm = TRUE)
  n_tot <- nrow(x$coefs)

  cat(sprintf("mlVAR result: %d subjects, %d observations, %d variables (lag %d)\n",
              x$n_subjects, x$n_obs, d, x$lag))
  cat(sprintf("  Temporal network:        %d x %d directed  (%d/%d edges significant at p<0.05)\n",
              d, d, n_sig, n_tot))
  cat(sprintf("  Contemporaneous network: %d x %d undirected\n", d, d))
  cat(sprintf("  Between network:         %d x %d undirected\n", d, d))
  invisible(x)
}

#' Summary method for net_mlvar
#'
#' @param object A `net_mlvar` object returned by [build_mlvar()].
#' @param ... Unused; present for S3 consistency.
#' @return Invisibly returns `object`.
#' @export
summary.net_mlvar <- function(object, ...) {
  vars <- rownames(object$temporal)
  d <- length(vars)

  cat("=== mlVAR Summary ===\n")
  cat("Subjects:", object$n_subjects,
      " | Observations:", object$n_obs,
      " | Variables:", d,
      " | Lag:", object$lag, "\n")
  cat("Standardized:", object$standardize, "\n\n")

  cat("--- Temporal Network (B matrix) ---\n")
  print(round(object$temporal, 4))
  cat("\n")

  sig_rows <- object$coefs[object$coefs$significant, , drop = FALSE]
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
  print(round(object$contemporaneous, 4))
  cat("\n")

  cat("--- Between-Subjects Network ---\n")
  print(round(object$between, 4))
  cat("\n")

  invisible(object)
}

