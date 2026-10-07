# ---- Mixed Graphical Models (MGM) ----

#' MGM estimator hook for the registry
#'
#' Registered under \code{method = "mgm"} (aliases \code{"mixed"},
#' \code{"mixed_graphical"}). Type detection and input validation live here;
#' the nodewise L1 regressions, EBIC selection, LW thresholding and
#' symmetrization are \code{psychnets::mgm_fit()} with its \pkg{glmnet} engine
#' (\code{native = FALSE}), which reproduces Nestimate's former in-package
#' implementation exactly. The \pkg{mgm} package is never loaded.
#'
#' @param data Data frame or matrix, one column per variable.
#' @param type Character vector of \code{"g"}/\code{"c"} per column, or
#'   \code{NULL} (default) to auto-detect: factor or character becomes
#'   \code{"c"}; a numeric column whose non-NA values are whole numbers with
#'   at most 10 distinct values becomes \code{"c"} (and a warning names it,
#'   because that silently switches the node from Gaussian to multinomial);
#'   everything else becomes \code{"g"}. The rule follows \code{mgm::mgm()}.
#' @param level Integer vector of category counts per column, or \code{NULL}
#'   (default) to derive it from \code{type}. Kept only for
#'   \code{mgm::mgm()} API parity -- it is length-validated and echoed back,
#'   never used in the fit. \strong{Name
#'   collision:} \code{build_network()} has its own \code{level} argument
#'   for multilevel decomposition, so \code{level} passed through
#'   \code{build_network()} does not reach this estimator. Call
#'   \code{.estimator_mgm()} directly with explicit \code{type} and
#'   \code{level}, or pre-coerce the categorical columns with
#'   \code{as.factor()}, when the auto-detection is wrong.
#' @param lambdaGam EBIC gamma. Default 0.25.
#' @param ruleReg Symmetrization rule, \code{"AND"} (default) or \code{"OR"}.
#' @param threshold Thresholding rule, \code{"LW"} (default) or
#'   \code{"none"}.
#' @param scale Must be \code{TRUE} (the default): continuous columns,
#'   including each node's own outcome, are standardized before fitting.
#'   \code{FALSE} is an error, because the estimator does not support it.
#' @param ... Ignored.
#' @return The standard estimator list: \code{matrix} (the symmetric
#'   weighted adjacency), \code{nodes}, \code{directed} (always
#'   \code{FALSE}), \code{cleaned_data}, plus the resolved \code{type} and
#'   \code{level} vectors.
#' @noRd
.estimator_mgm <- function(data, type = NULL, level = NULL,
                            lambdaGam = 0.25, ruleReg = "AND",
                            threshold = "LW", scale = TRUE, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Method 'mgm' requires the 'glmnet' package.", call. = FALSE)
  }
  ruleReg <- match.arg(ruleReg, c("AND", "OR"))
  threshold <- match.arg(threshold, c("LW", "none"))
  if (!isTRUE(scale)) {
    stop(errorCondition(paste0("`scale = FALSE` is not supported by the mgm estimator: ",
                               "continuous columns are always standardized."),
                        class = "nestimate_mgm_unscaled", call = NULL))
  }
  data <- as.data.frame(data)
  p <- ncol(data)

  # Auto-detect types if not provided: factor / character / integer with
  # few unique values -> categorical; otherwise gaussian
  if (is.null(type)) {
    type <- vapply(data, function(col) {
      if (is.factor(col) || is.character(col)) "c"
      else if (is.numeric(col) &&
               length(unique(col[!is.na(col)])) <= 10 &&
               all(col == round(col), na.rm = TRUE)) "c"
      else "g"
    }, character(1))
    # Detection mirrors mgm::mgm()'s default, but the surprising case -- a
    # numeric column auto-classified CATEGORICAL via the <=10-integer rule
    # (e.g. a Likert/count item) -- must never be silent: it changes the
    # node's model from Gaussian to multinomial. Announce it; the user
    # overrides with explicit `type=`/`level=`.
    num_cat <- vapply(data, is.numeric, logical(1)) & type == "c"
    if (any(num_cat)) {
      warning(
        "mgm: numeric column(s) auto-detected as CATEGORICAL ",
        "(integer with <=10 distinct values): ",
        paste(names(data)[num_cat], collapse = ", "),
        ". If these are continuous, pass explicit `type=`/`level=`.",
        call. = FALSE
      )
    }
  }
  stopifnot(
    "`type` must have one \"g\" or \"c\" per column" =
      is.character(type) && length(type) == p && all(type %in% c("g", "c"))
  )
  if (is.null(level)) {
    level <- vapply(seq_along(type), function(i) {
      if (type[i] == "c") length(unique(data[[i]])) else 1L
    }, integer(1))
  }
  stopifnot("`level` must be numeric with one value per column" =
              is.numeric(level) && length(level) == p)

  # Categorical columns go to the fit as factors, as mgm::mgm() treats them.
  fit_data <- data
  fit_data[type == "c"] <- lapply(fit_data[type == "c"], as.factor)
  # Missing values: each nodewise regression uses the rows complete for its
  # variables (pairwise). The pre-0.9.23 in-package code failed on any NA.
  fit <- psychnets::mgm_fit(fit_data, gamma = lambdaGam, types = unname(type),
                            threshold = threshold, rule = ruleReg,
                            na_method = "pairwise", native = FALSE)
  wadj <- fit$weights
  dimnames(wadj) <- list(colnames(data), colnames(data))

  list(
    matrix       = wadj,
    nodes        = colnames(data),
    directed     = FALSE,
    cleaned_data = data,
    type         = type,
    level        = level
  )
}
