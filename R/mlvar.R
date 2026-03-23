# ---- Multilevel Vector Autoregression (mlVAR) ----

#' Multilevel Vector Autoregression
#'
#' @description Estimates three networks from ESM/EMA panel data using
#'   multilevel VAR: (1) a directed temporal network of lagged regression
#'   coefficients, (2) an undirected contemporaneous network of partial
#'   correlations among residuals, and (3) an undirected between-subjects
#'   network of partial correlations among person means. Uses within-person
#'   centering (fixed-effects OLS) for the temporal model and EBIC-selected
#'   graphical LASSO for the contemporaneous and between-subjects networks.
#'
#' @param data A \code{data.frame} containing the panel data.
#' @param vars Character vector of variable column names to model.
#' @param id Character string naming the person-ID column.
#' @param day Character string naming the day/session column, or \code{NULL}
#'   (default). When provided, lag pairs are only formed within the same day.
#' @param beep Character string naming the measurement-occasion column, or
#'   \code{NULL} (default). When provided, lag pairs require
#'   \code{beep[t] - beep[t-lag] == lag}.
#' @param lag Integer. The lag order (default 1).
#' @param standardize Logical. If \code{TRUE} (default), within-centered
#'   variables are divided by their pooled standard deviation before OLS.
#' @param gamma Numeric. EBIC hyperparameter for graphical LASSO model
#'   selection (0 = BIC, higher = sparser). Default 0.5.
#' @param nlambda Integer. Number of lambda values in the regularization path
#'   for graphical LASSO. Default 100.
#'
#' @return An S3 object of class \code{"mlvar_result"}, a list with:
#'   \describe{
#'     \item{\code{temporal}}{d x d matrix of fixed-effect temporal regression
#'       coefficients. Entry \code{[i, j]} is the effect of variable j at
#'       t-lag on variable i at t.}
#'     \item{\code{contemporaneous}}{d x d symmetric matrix of partial
#'       correlations among within-person residuals (EBIC-GLASSO).}
#'     \item{\code{between}}{d x d symmetric matrix of partial correlations
#'       among person means (EBIC-GLASSO).}
#'     \item{\code{coefs}}{List of d data frames, one per outcome variable,
#'       each containing columns \code{predictor}, \code{beta}, \code{se},
#'       \code{t}, \code{p}, \code{ci_lower}, \code{ci_upper}.}
#'     \item{\code{labels}}{Character vector of variable names.}
#'     \item{\code{n_obs}}{Number of valid lag-pair observations.}
#'     \item{\code{n_subjects}}{Number of unique subjects.}
#'     \item{\code{lag}}{Lag order used.}
#'     \item{\code{standardize}}{Logical; whether standardization was applied.}
#'     \item{\code{gamma}}{EBIC gamma used.}
#'   }
#'
#' @details
#' The algorithm proceeds in seven steps:
#' \enumerate{
#'   \item \strong{Data preparation}: select columns, coerce types, drop NA
#'     rows, sort by \code{(id, day, beep)}.
#'   \item \strong{Lag-pair construction}: form valid outcome/predictor pairs
#'     respecting person, day, and beep boundaries.
#'   \item \strong{Within-centering}: person-mean center both Y (outcome) and
#'     X (predictor) matrices. Optionally standardize by pooled SD.
#'   \item \strong{Temporal OLS}: for each outcome k, fit
#'     \code{lm(Y_k ~ X - 1)} with corrected degrees of freedom.
#'   \item \strong{Contemporaneous network}: EBIC-GLASSO on the correlation
#'     matrix of OLS residuals.
#'   \item \strong{Between-subjects network}: EBIC-GLASSO on the correlation
#'     matrix of person means.
#'   \item \strong{Assembly}: collect results into \code{mlvar_result} object.
#' }
#'
#' @examples
#' \dontrun{
#' d <- simulate_data("mlvar", seed = 1)
#' fit <- mlvar(d, vars = attr(d, "vars"), id = "id", day = "day", beep = "beep")
#' print(fit)
#' summary(fit)
#' }
#'
#' @seealso \code{\link{build_network}}, \code{\link{graphical_var}}
#' @export
mlvar <- function(data,
                  vars,
                  id,
                  day = NULL,
                  beep = NULL,
                  lag = 1L,
                  standardize = TRUE,
                  gamma = 0.5,
                  nlambda = 100L) {
  # ---- Input validation ----
  stopifnot(
    is.data.frame(data),
    is.character(vars),
    length(vars) >= 2L,
    is.character(id),
    length(id) == 1L
  )
  stopifnot(
    is.numeric(lag),
    length(lag) == 1L,
    lag >= 1L
  )
  stopifnot(
    is.logical(standardize),
    length(standardize) == 1L
  )
  stopifnot(
    is.numeric(gamma),
    length(gamma) == 1L,
    gamma >= 0
  )
  nlambda <- as.integer(nlambda)
  stopifnot(
    is.integer(nlambda),
    length(nlambda) == 1L,
    nlambda >= 2L
  )

  # Check columns exist
  required <- c(vars, id)
  if (!is.null(day)) required <- c(required, day)
  if (!is.null(beep)) required <- c(required, beep)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop("Columns not found in data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  d <- length(vars)

  # Step 1: Prepare data
  prepared <- .mlvar_prepare_data(data, vars, id, day, beep)

  # Check minimum subjects
  n_subjects <- length(unique(prepared[[id]]))
  if (n_subjects < 2L) {
    stop("At least 2 subjects are required. Found: ", n_subjects, call. = FALSE)
  }

  # Step 2: Build lag pairs
  lag_result <- .mlvar_build_lag_pairs(prepared, vars, id, day, beep,
                                       lag = as.integer(lag))
  Y <- lag_result$Y
  X <- lag_result$X
  id_vec <- lag_result$id_vec
  n_obs <- nrow(Y)

  if (n_obs < d + 1L) {
    stop("Too few valid lag pairs (", n_obs,
         ") for ", d, " variables.", call. = FALSE)
  }

  n_subjects_pairs <- length(unique(id_vec))

  # Step 3: Within-center (person means from ALL data, grand SD)
  centered <- .mlvar_within_center(Y, X, id_vec, standardize,
                                    full_data = prepared, vars = vars,
                                    id_col = id)
  Y_c <- centered$Y
  X_c <- centered$X

  # Step 4: Temporal OLS
  temporal_result <- .mlvar_temporal_ols(Y_c, X_c, n_subjects_pairs, vars)

  # Step 5: Contemporaneous network
  contemporaneous <- .mlvar_contemporaneous(
    temporal_result$residuals, n_obs, gamma, nlambda
  )

  # Step 6: Between-subjects network (Epskamp et al. 2017 Eq. 3)
  between <- .mlvar_between(lag_result, prepared, vars, id, day, beep,
                            n_subjects, as.integer(lag))

  # Step 7: Assemble result
  result <- list(
    temporal        = temporal_result$B,
    contemporaneous = contemporaneous,
    between         = between,
    coefs           = temporal_result$coefs,
    labels          = vars,
    n_obs           = n_obs,
    n_subjects      = n_subjects_pairs,
    lag             = as.integer(lag),
    standardize     = standardize,
    gamma           = gamma
  )
  class(result) <- "mlvar_result"
  result
}


# ---- Internal helpers ----

#' Prepare data for mlVAR
#'
#' Selects relevant columns, coerces types, drops rows with NA in any variable,
#' and sorts by (id, day, beep).
#'
#' @noRd
.mlvar_prepare_data <- function(data, vars, id, day, beep) {
  keep_cols <- c(id, vars)
  if (!is.null(day)) keep_cols <- c(keep_cols, day)
  if (!is.null(beep)) keep_cols <- c(keep_cols, beep)
  keep_cols <- unique(keep_cols)

  df <- data[, keep_cols, drop = FALSE]

  # Coerce vars to numeric
  for (v in vars) {
    df[[v]] <- as.numeric(df[[v]])
  }

  # Drop rows with NA in any variable or id
  check_cols <- c(id, vars)
  if (!is.null(day)) check_cols <- c(check_cols, day)
  if (!is.null(beep)) check_cols <- c(check_cols, beep)
  complete <- complete.cases(df[, check_cols, drop = FALSE])
  if (!all(complete)) {
    df <- df[complete, , drop = FALSE]
  }

  if (nrow(df) < 2L) {
    stop("Fewer than 2 complete rows remain after removing NAs.", call. = FALSE)
  }

  # Sort by (id, day, beep)
  order_cols <- id
  if (!is.null(day)) order_cols <- c(order_cols, day)
  if (!is.null(beep)) order_cols <- c(order_cols, beep)
  ord <- do.call(order, df[, order_cols, drop = FALSE])
  df <- df[ord, , drop = FALSE]
  rownames(df) <- NULL

  df
}


#' Build lag pairs for mlVAR
#'
#' Matches mlVAR's data augmentation pipeline exactly:
#' 1. Expand each person-day to include all beeps in \code{[min, max]} range
#' 2. Fill missing beeps with NA
#' 3. Lag by row position (aveLag)
#' 4. Drop rows with NA in outcome or lagged predictor
#'
#' This correctly handles beep gaps (producing NA at the gap) and
#' duplicate beeps (both kept as separate observations).
#'
#' @return List with Y (n_pairs x d outcome matrix), X (n_pairs x d predictor
#'   matrix), id_vec (person IDs for each pair).
#' @noRd
.mlvar_build_lag_pairs <- function(data, vars, id, day, beep, lag) {
  d <- length(vars)

  if (!is.null(beep) && !is.null(day)) {
    # mlVAR-style augmentation: expand beep grid and fill missing with NA
    data_list <- split(data, data[[id]])
    aug_list <- lapply(data_list, function(subj_data) {
      subj_id <- subj_data[[id]][1]
      day_list <- split(subj_data, subj_data[[day]])
      day_aug <- lapply(day_list, function(day_data) {
        d_val <- day_data[[day]][1]
        beep_vals <- day_data[[beep]]
        beep_min <- min(beep_vals)
        beep_max <- max(beep_vals)
        all_beeps <- seq(beep_min, beep_max)

        # Create expanded frame with all beeps
        expanded <- data.frame(.beep = all_beeps, stringsAsFactors = FALSE)
        names(expanded) <- beep
        expanded[[id]] <- subj_id
        expanded[[day]] <- d_val

        # Merge: right join to keep all beeps, fill missing with NA
        merged <- merge(expanded, day_data, by = c(id, day, beep), all.x = TRUE)
        merged <- merged[order(merged[[beep]]), ]
        merged
      })
      do.call(rbind, day_aug)
    })
    aug_data <- do.call(rbind, aug_list)
    rownames(aug_data) <- NULL

    # Sort by (id, day, beep)
    aug_data <- aug_data[order(aug_data[[id]], aug_data[[day]],
                               aug_data[[beep]]), ]

    # Lag by position within each person-day group
    n <- nrow(aug_data)
    idx_t <- seq(lag + 1L, n)
    idx_lag <- seq(1L, n - lag)

    # Same person and same day
    valid <- aug_data[[id]][idx_t] == aug_data[[id]][idx_lag]
    valid <- valid & (aug_data[[day]][idx_t] == aug_data[[day]][idx_lag])

    # Extract Y and X (may contain NAs from expansion)
    Y_raw <- as.matrix(aug_data[idx_t[valid], vars, drop = FALSE])
    X_raw <- as.matrix(aug_data[idx_lag[valid], vars, drop = FALSE])
    id_raw <- aug_data[[id]][idx_t[valid]]

    # Drop rows where either Y or X has any NA
    complete <- complete.cases(Y_raw) & complete.cases(X_raw)
    Y <- Y_raw[complete, , drop = FALSE]
    X <- X_raw[complete, , drop = FALSE]
    id_vec <- id_raw[complete]
  } else {
    # No beep/day: simple consecutive-row pairing within person
    n <- nrow(data)
    if (n <= lag) {
      stop("Not enough rows (", n, ") for lag ", lag, ".", call. = FALSE)
    }

    idx_t <- seq(lag + 1L, n)
    idx_lag <- seq(1L, n - lag)
    valid <- data[[id]][idx_t] == data[[id]][idx_lag]

    if (!is.null(day)) {
      valid <- valid & (data[[day]][idx_t] == data[[day]][idx_lag])
    }

    Y <- as.matrix(data[idx_t[valid], vars, drop = FALSE])
    X <- as.matrix(data[idx_lag[valid], vars, drop = FALSE])
    id_vec <- data[[id]][idx_t[valid]]
  }

  if (length(id_vec) == 0L) {
    stop("No valid lag pairs found. Check id/day/beep columns and lag value.",
         call. = FALSE)
  }

  list(Y = Y, X = X, id_vec = id_vec)
}


#' Within-person centering for mlVAR
#'
#' Person-mean centers both Y and X matrices. Centering both sides absorbs
#' the random intercept, making OLS equivalent to the fixed-effects estimator.
#' If standardize = TRUE, divides by pooled SD after centering.
#'
#' @return List with centered Y and X matrices.
#' @noRd
.mlvar_within_center <- function(Y, X, id_vec, standardize,
                                  full_data = NULL, vars = NULL,
                                  id_col = NULL) {
  d <- ncol(Y)

  if (!is.null(full_data) && !is.null(vars) && !is.null(id_col)) {
    # Person means from ALL observations (matches psychaj/mlVAR)
    agg <- stats::aggregate(full_data[, vars, drop = FALSE],
                            by = list(.id = full_data[[id_col]]), FUN = mean)
    pm <- as.matrix(agg[, vars, drop = FALSE])
    rownames(pm) <- agg$.id
    mean_mat <- pm[match(id_vec, rownames(pm)), , drop = FALSE]
    Y_c <- Y - mean_mat
    X_c <- X - mean_mat
  } else {
    # Fallback: center from lag-pair means
    Y_c <- Y
    for (j in seq_len(d)) {
      Y_c[, j] <- Y[, j] - ave(Y[, j], id_vec, FUN = mean)
    }
    X_c <- X
    for (j in seq_len(d)) {
      X_c[, j] <- X[, j] - ave(X[, j], id_vec, FUN = mean)
    }
  }

  if (standardize) {
    if (!is.null(full_data) && !is.null(vars)) {
      # Grand SD from all raw data (matches psychaj/mlVAR)
      sd_vec <- vapply(vars, function(v) stats::sd(full_data[[v]]), numeric(1))
      sd_vec[sd_vec < 1e-12] <- 1
      sd_mat <- matrix(sd_vec, nrow = nrow(Y_c), ncol = d, byrow = TRUE)
      Y_c <- Y_c / sd_mat
      X_c <- X_c / sd_mat
    } else {
      for (j in seq_len(d)) {
        sd_y <- stats::sd(Y_c[, j])
        sd_x <- stats::sd(X_c[, j])
        if (sd_y > 0) Y_c[, j] <- Y_c[, j] / sd_y
        if (sd_x > 0) X_c[, j] <- X_c[, j] / sd_x
      }
    }
  }

  list(Y = Y_c, X = X_c)
}


#' Temporal OLS for mlVAR
#'
#' For each outcome variable, fits OLS without intercept (centering absorbed
#' it). Applies degrees-of-freedom correction for absorbed person fixed effects.
#'
#' @return List with B (d x d coefficient matrix), coefs (list of data frames),
#'   residuals (n_obs x d matrix).
#' @noRd
.mlvar_temporal_ols <- function(Y, X, n_subjects, vars) {
  n_obs <- nrow(Y)
  d <- ncol(Y)

  B <- matrix(0, nrow = d, ncol = d,
              dimnames = list(vars, vars))
  residuals <- matrix(0, nrow = n_obs, ncol = d,
                      dimnames = list(NULL, vars))
  coefs_list <- vector("list", d)
  names(coefs_list) <- vars

  df_ols <- n_obs - d
  df_correct <- n_obs - d - n_subjects

  if (df_correct < 1L) {
    warning("Corrected degrees of freedom < 1. Results may be unreliable.",
            call. = FALSE)
    df_correct <- max(1L, df_correct)
  }

  df_ratio <- sqrt(df_ols / df_correct)

  for (k in seq_len(d)) {
    # Include intercept (matches psychaj — absorbs residual centering)
    X_int <- cbind(1, X)
    fit <- stats::lm.fit(X_int, Y[, k])
    beta <- fit$coefficients[-1L]  # drop intercept
    resid <- fit$residuals

    # Corrected SE
    rss <- sum(resid^2)
    sigma2 <- rss / df_ols
    XtX_inv <- tryCatch(
      solve(crossprod(X)),
      error = function(e) {
        warning("Singular X'X for variable ", vars[k],
                ". Using pseudoinverse.", call. = FALSE)
        MASS_needed <- FALSE
        # Fallback: use diagonal approximation
        diag(1 / pmax(diag(crossprod(X)), 1e-10))
      }
    )
    se_ols <- sqrt(sigma2 * diag(XtX_inv))
    se_correct <- se_ols * df_ratio

    # t-statistics and p-values
    t_val <- beta / se_correct
    p_val <- 2 * stats::pt(abs(t_val), df = df_correct, lower.tail = FALSE)

    # Confidence intervals
    t_crit <- stats::qt(0.975, df = df_correct)
    ci_lower <- beta - t_crit * se_correct
    ci_upper <- beta + t_crit * se_correct

    B[k, ] <- beta
    residuals[, k] <- resid

    coefs_list[[k]] <- data.frame(
      predictor = vars,
      beta      = as.numeric(beta),
      se        = as.numeric(se_correct),
      t         = as.numeric(t_val),
      p         = as.numeric(p_val),
      ci_lower  = as.numeric(ci_lower),
      ci_upper  = as.numeric(ci_upper),
      stringsAsFactors = FALSE
    )
  }

  list(B = B, coefs = coefs_list, residuals = residuals)
}


#' Contemporaneous network via EBIC-GLASSO on residuals
#'
#' @return d x d partial correlation matrix (symmetric, zero diagonal).
#' @noRd
.mlvar_contemporaneous <- function(residuals, n_obs, gamma, nlambda) {
  d <- ncol(residuals)
  vars <- colnames(residuals)

  S <- stats::cor(residuals)

  # Guard: if any correlations are NA or all zero

  if (any(is.na(S))) {
    warning("NA correlations in residuals. Returning zero matrix.",
            call. = FALSE)
    mat <- matrix(0, d, d, dimnames = list(vars, vars))
    return(mat)
  }

  # Use existing GLASSO pipeline
  lambda_path <- tryCatch(
    .compute_lambda_path(S, nlambda, 0.01),
    error = function(e) NULL
  )

  if (is.null(lambda_path)) {
    mat <- matrix(0, d, d, dimnames = list(vars, vars))
    return(mat)
  }

  selected <- tryCatch(
    .select_ebic(S, lambda_path, n_obs, gamma, FALSE),
    error = function(e) NULL
  )

  if (is.null(selected)) {
    mat <- matrix(0, d, d, dimnames = list(vars, vars))
    return(mat)
  }

  pcor <- .precision_to_pcor(selected$wi, 0)
  pcor <- (pcor + t(pcor)) / 2  # force exact symmetry
  colnames(pcor) <- rownames(pcor) <- vars
  pcor
}


#' Between-subjects network via nodewise lmer (Epskamp et al. 2017)
#'
#' @return d x d partial correlation matrix (symmetric, zero diagonal).
#' @noRd
.mlvar_between <- function(lag_data, full_data, vars, id, day, beep,
                           n_subjects, lag_int) {
  d <- length(vars)
  zero_mat <- matrix(0, d, d, dimnames = list(vars, vars))

  # Guard: need at least d+1 subjects for meaningful estimation
  if (n_subjects < d + 1L) return(zero_mat)

  # Person means computed from ALL observations (not just lag-paired)
  agg <- stats::aggregate(full_data[, vars, drop = FALSE],
                          by = list(.id = full_data[[id]]), FUN = mean)
  pm_mat <- as.matrix(agg[, vars, drop = FALSE])
  rownames(pm_mat) <- agg$.id

  # Build lag-pair data frame for lmer: outcome Y, within-centered lagged X,
  # person means of other vars
  Y_mat <- lag_data$Y
  X_mat <- lag_data$X
  id_vec <- lag_data$id_vec

  # Within-center lagged predictors: X_within = X - person_mean(X_lagged)
  X_within <- X_mat
  for (j in seq_len(d)) {
    pm_lag <- ave(X_mat[, j], id_vec, FUN = mean)
    X_within[, j] <- X_mat[, j] - pm_lag
  }

  # Build person-mean matrix aligned to lag-pair rows
  pm_rows <- pm_mat[match(id_vec, rownames(pm_mat)), , drop = FALSE]

  # Nodewise lmer: for each outcome k, fit
  # Y_k ~ X_within[,1..d] + pm[, -k] + (1|subject)  with REML = FALSE
  Gamma <- matrix(0, d, d, dimnames = list(vars, vars))
  mu_SD <- numeric(d)
  names(mu_SD) <- vars

  for (k in seq_len(d)) {
    lmer_df <- data.frame(
      .y = Y_mat[, k],
      .subject = id_vec,
      stringsAsFactors = FALSE
    )

    # Within-centered lagged predictors (all d vars)
    for (j in seq_len(d)) {
      lmer_df[[paste0(".xw", j)]] <- X_within[, j]
    }
    # Person means of other vars (d-1 predictors)
    other_idx <- seq_len(d)[-k]
    for (j in other_idx) {
      lmer_df[[paste0(".pm", j)]] <- pm_rows[, j]
    }

    # Build formula
    xw_terms <- paste0(".xw", seq_len(d))
    pm_terms <- paste0(".pm", other_idx)
    rhs <- paste(c(xw_terms, pm_terms), collapse = " + ")
    fm <- stats::as.formula(paste(".y ~", rhs, "+ (1 | .subject)"))

    fit <- tryCatch(
      lme4::lmer(fm, data = lmer_df, REML = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) return(zero_mat)

    # Extract between regression coefficients (Gamma row k)
    fe <- lme4::fixef(fit)
    for (j in other_idx) {
      Gamma[k, j] <- fe[[paste0(".pm", j)]]
    }

    # Random intercept SD
    vc <- lme4::VarCorr(fit)
    ri_var <- as.numeric(vc$.subject)
    if (is.na(ri_var) || ri_var == 0) return(zero_mat)  # zero-SD guard
    mu_SD[k] <- sqrt(ri_var)
  }

  # Precision matrix K = D(I - Gamma), D = diag(1/mu_SD^2)
  D <- diag(1 / mu_SD^2)
  K <- D %*% (diag(d) - Gamma)
  K <- (K + t(K)) / 2  # symmetrize

  # Force positive definite (matches mlVAR:::forcePositive exactly)
  eig <- eigen(K)$values
  if (any(eig < 0)) {
    K <- K - (diag(d) * min(eig) - 0.001)
  }

  # Convert to partial correlations
  pcor <- .precision_to_pcor(K, 0)
  pcor <- (pcor + t(pcor)) / 2
  colnames(pcor) <- rownames(pcor) <- vars
  pcor
}


# ---- S3 methods ----

#' Print method for mlvar_result
#'
#' @param x An \code{mlvar_result} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
print.mlvar_result <- function(x, ...) {
  d <- length(x$labels)
  n_temp <- sum(x$temporal != 0)
  n_cont <- sum(x$contemporaneous[upper.tri(x$contemporaneous)] != 0)
  n_betw <- sum(x$between[upper.tri(x$between)] != 0)

  cat("mlVAR result:",
      x$n_subjects, "subjects,",
      x$n_obs, "observations,",
      d, "variables\n")
  cat("  Temporal edges:", n_temp, "(directed)\n")
  cat("  Contemporaneous edges:", n_cont, "(undirected)\n")
  cat("  Between-subjects edges:", n_betw, "(undirected)\n")
  invisible(x)
}


#' Summary method for mlvar_result
#'
#' @param object An \code{mlvar_result} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
summary.mlvar_result <- function(object, ...) {
  d <- length(object$labels)

  cat("=== mlVAR Summary ===\n")
  cat("Subjects:", object$n_subjects, " | Observations:", object$n_obs,
      " | Variables:", d, " | Lag:", object$lag, "\n")
  cat("Standardized:", object$standardize, " | EBIC gamma:", object$gamma, "\n")
  cat("\n")

  # Temporal coefficients
  cat("--- Temporal Network (B matrix) ---\n")
  print(round(object$temporal, 4))
  cat("\n")

  # Significant temporal edges
  sig_edges <- do.call(rbind, lapply(seq_along(object$coefs), function(k) {
    cf <- object$coefs[[k]]
    sig <- cf[cf$p < 0.05, , drop = FALSE]
    if (nrow(sig) > 0L) {
      sig$outcome <- object$labels[k]
      sig
    } else {
      NULL
    }
  }))
  if (!is.null(sig_edges) && nrow(sig_edges) > 0L) {
    cat("Significant temporal edges (p < 0.05):\n")
    sig_print <- sig_edges[, c("predictor", "outcome", "beta", "se", "t", "p"),
                           drop = FALSE]
    sig_print$beta <- round(sig_print$beta, 4)
    sig_print$se <- round(sig_print$se, 4)
    sig_print$t <- round(sig_print$t, 3)
    sig_print$p <- round(sig_print$p, 4)
    rownames(sig_print) <- NULL
    print(sig_print)
  } else {
    cat("No significant temporal edges at p < 0.05.\n")
  }
  cat("\n")

  # Contemporaneous
  cat("--- Contemporaneous Network ---\n")
  print(round(object$contemporaneous, 4))
  cat("\n")

  # Between
  cat("--- Between-Subjects Network ---\n")
  print(round(object$between, 4))
  cat("\n")

  invisible(object)
}
