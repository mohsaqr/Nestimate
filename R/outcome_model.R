# ==============================================================================
# outcome_model() -- unit-level outcome regression with honest inference.
#
# Written to fix four defects in the usual "patterns -> outcome" routine:
#   * it is logistic-only, so a continuous outcome must be dichotomised;
#   * it ranks predictors by their association with the outcome and then
#     reports p-values from a model fitted to the SAME rows, which is
#     selective inference -- those p-values are not valid;
#   * it returns a raw model object rather than a tidy effect table;
#   * it corrects for nothing when many predictors are tested at once.
#
# Here: family is detected from the outcome, selection (when used) happens on
# a disjoint half of the data from the fit, every effect is reported with a
# confidence interval, and p-values carry a multiplicity correction.
# ==============================================================================

utils::globalVariables(c("term", "estimate", "ci_lower", "ci_upper"))

#' Model Unit-Level Outcomes from Sequence or Network Predictors
#'
#' @description
#' `outcome_model()` fits a regression of a unit-level outcome on predictor
#' columns (pattern indicators, topological features from
#' [simplicial_features()], or any other numeric covariates) and returns a
#' tidy effect table with confidence intervals and multiplicity-corrected
#' p-values. A two-valued outcome gets a logistic regression and a numeric
#' one a linear regression; a `group` column adds a random intercept per
#' group, and `select = "split"` chooses the predictors on one half of the
#' data and fits them on the other.
#'
#' Each effect is conditional on the other terms in the model; for a
#' two-valued outcome, the odds ratio is the factor by which the odds of the
#' outcome value coded 1 change per unit of the predictor. The p-values are
#' corrected across every term except the intercept, the `adjust` covariates
#' included. A predictor that counts events also measures how long the unit
#' was observed; enter that exposure through `adjust`.
#'
#' @section Honest inference:
#' Choosing predictors by their association with the outcome and then
#' testing them on the same rows invalidates the p-values. With
#' `select = "split"` the data is halved: predictors are ranked on one half
#' and the reported model is fitted on the other, so the returned inference
#' is valid for the selected set. `select = "none"` (default) fits every
#' supplied predictor and needs no split.
#'
#' @param data A data frame with one row per unit of analysis. Rows with a
#'   missing value in the outcome, a predictor, an `adjust` covariate or the
#'   `group` column are left out; when none is left, an error of class
#'   `nestimate_empty_model_frame` is raised. A column that is not in `data`
#'   raises an error.
#' @param outcome The outcome column. With `family = "auto"`, an outcome
#'   with two distinct values is modelled with the binomial family (the
#'   second value in sort order, level order for a factor, or `TRUE`,
#'   coded 1) and any other with the gaussian
#'   family, which needs a numeric outcome.
#' @param predictors The predictor columns. A predictor that is not numeric
#'   (a logical column included; convert it with `as.integer()`) or has no
#'   variance in the complete rows is dropped and listed in the result; when
#'   every predictor is dropped, an error of class `nestimate_no_variance` is
#'   raised.
#' @param group A grouping column, such as the actor the units are nested
#'   in. When \pkg{lme4} is installed, a random intercept per group is added,
#'   the right treatment for units nested in actors. \pkg{lme4} is a
#'   suggested package: without it the random intercept is dropped, a
#'   warning of class `nestimate_no_lme4` is raised, and a plain
#'   [stats::glm()] is fitted. Default `NULL`.
#' @param adjust Covariates entered before the predictors. Use it for
#'   exposure: a unit observed longer contains more of every pattern, so an
#'   unadjusted effect can be volume in disguise. Default `NULL`.
#' @param family The model family: `"auto"` (default) reads it from the
#'   outcome, `"binomial"` and `"gaussian"` force one. `"binomial"` codes the
#'   second-smallest outcome value (or `TRUE`) as 1 and every other value as
#'   0, so give it a two-valued outcome.
#' @param select How the predictors are chosen: `"none"` (default) fits all
#'   of them; `"split"` ranks them on a random half of the rows by the
#'   absolute z or t statistic of a one-predictor `glm()` (without `adjust`
#'   or `group`), keeps the best `n_select`, and fits the reported model on
#'   the other half.
#' @param n_select The number of predictors kept when `select = "split"`.
#'   Default `10`. Ignored otherwise.
#' @param correction The multiplicity correction, one of
#'   `stats::p.adjust.methods`, passed to [stats::p.adjust()]. Default
#'   `"BH"`.
#' @param ci_level The confidence level of the intervals, between 0 and 1.
#'   Default `0.95`.
#' @param seed An integer that makes the split reproducible; used only with
#'   `select = "split"`. A random-number state that existed before the call
#'   is restored afterwards. Default `NULL`.
#'
#' @return A `net_outcome_model` object. `print()` shows the outcome and
#'   family, the number of units (and groups), the number of selected terms,
#'   the adjustments, the dropped predictors and the correction, followed by
#'   the effect table, and returns its input invisibly. `summary()` and
#'   [effects_table()] return the effect table without the intercept,
#'   rounded to three decimals except the p-values; [effects_table()] also
#'   keeps the intercept, or only the significant terms, on request.
#'   `plot()` returns a `ggplot` forest plot of the effects (odds ratios for
#'   a binomial fit, against a dashed line at no effect), with filled markers
#'   for terms whose corrected p-value is below 0.05. The object's fields:
#' \describe{
#'   \item{effects}{The effect table: a data frame with one row per model
#'     term, the intercept and the `adjust` covariates included, and columns
#'     `term`, `estimate`, `std_error`, `statistic`, `ci_lower`, `ci_upper`,
#'     `p_value` and `p_adj` (`NA` on the intercept row, which is left out
#'     of the correction), plus `odds_ratio`, `or_lower` and `or_upper` for a
#'     binomial fit.}
#'   \item{model}{The fitted model: a `glm`, or an \pkg{lme4} fit when a
#'     random intercept was added.}
#'   \item{family}{`"binomial"` or `"gaussian"`.}
#'   \item{n}{The number of rows the reported model was fitted on.}
#'   \item{n_groups}{The number of groups; `NA` unless a random intercept
#'     was fitted.}
#'   \item{selected}{The predictors in the reported model.}
#'   \item{dropped}{The predictors dropped as non-numeric or constant.}
#'   \item{adjust, select, correction, ci_level, outcome}{The settings
#'     used.}
#'   \item{mixed}{Whether a random intercept was fitted.}
#' }
#'
#' @details
#' Ordinary gaussian fits use Student t confidence limits with the model's
#' residual degrees of freedom; binomial and mixed fits use normal (Wald)
#' limits. A linear mixed model (\code{lme4::lmer()}) reports no p-values,
#' so they are computed from the normal distribution of its t statistics.
#'
#' @seealso \code{\link{simplicial_features}}, \code{\link{effects_table}}
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   hint    = rbinom(300, 1, 0.4),
#'   think   = rbinom(300, 1, 0.3),
#'   n_events = rpois(300, 20),
#'   actor   = rep(letters[1:10], each = 30)
#' )
#' d$success <- rbinom(300, 1, plogis(-0.5 + 0.8 * d$hint))
#' fit <- outcome_model(d, outcome = "success",
#'                      predictors = c("hint", "think"),
#'                      adjust = "n_events")
#' effects_table(fit)
#' @export
outcome_model <- function(data, outcome, predictors, group = NULL,
                          adjust = NULL,
                          family = c("auto", "binomial", "gaussian"),
                          select = c("none", "split"), n_select = 10L,
                          correction = "BH", ci_level = 0.95, seed = NULL) {
  family <- match.arg(family)
  select <- match.arg(select)
  .assert(
    "`data` must be a data.frame" = is.data.frame(data),
    "`outcome` must be a single column name" =
      is.character(outcome) && length(outcome) == 1L,
    "`predictors` must be a character vector" =
      is.character(predictors) && length(predictors) >= 1L,
    "`ci_level` must be between 0 and 1" =
      is.numeric(ci_level) && length(ci_level) == 1L &&
      ci_level > 0 && ci_level < 1,
    "`n_select` must be a single positive integer" =
      length(n_select) == 1L && is.finite(n_select) && n_select >= 1
  )
  needed <- c(outcome, predictors, adjust, group)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop(errorCondition(paste0("Column(s) not found in `data`: ",
         paste(utils::head(missing_cols, 5L), collapse = ", ")),
                        class = "nestimate_bad_input", call = NULL))
  }
  d <- data[stats::complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  if (nrow(d) == 0L) {
    stop(errorCondition(
      "No complete rows across the outcome, predictors and covariates.",
      class = "nestimate_empty_model_frame", call = NULL))
  }

  y <- d[[outcome]]
  fam <- .om_family(y, family)
  if (identical(fam, "binomial")) d[[outcome]] <- .om_binary(y)

  # Drop predictors with no variance -- they cannot carry an effect, and in a
  # binomial fit they produce separation rather than an error.
  keep <- vapply(predictors, function(p) {
    v <- d[[p]]
    is.numeric(v) && stats::var(v) > 0
  }, logical(1L))
  dropped <- predictors[!keep]
  predictors <- predictors[keep]
  if (length(predictors) == 0L) {
    stop(errorCondition(
      "Every predictor is constant after dropping incomplete rows.",
      class = "nestimate_no_variance", call = NULL))
  }

  selected <- predictors
  fit_rows <- seq_len(nrow(d))
  if (identical(select, "split")) {
    sp <- .om_split(d, outcome, predictors, fam, n_select, seed)
    selected <- sp$selected
    fit_rows <- sp$fit_rows
  }

  fit <- .om_fit(d[fit_rows, , drop = FALSE], outcome, selected, adjust,
                 group, fam)
  eff <- .om_effects(fit$model, fam, ci_level, correction, fit$mixed)

  structure(
    list(effects = eff, model = fit$model, family = fam,
         n = length(fit_rows), n_groups = fit$n_groups,
         selected = selected, dropped = dropped, adjust = adjust,
         select = select, correction = correction, ci_level = ci_level,
         mixed = fit$mixed, outcome = outcome),
    class = "net_outcome_model")
}

# ---- internals --------------------------------------------------------------

.om_family <- function(y, family) {
  if (!identical(family, "auto")) return(family)
  u <- unique(y[!is.na(y)])
  if (length(u) == 2L) "binomial" else "gaussian"
}

# Map a two-valued outcome to 0/1, keeping the lower level as the reference.
.om_binary <- function(y) {
  if (is.logical(y)) return(as.integer(y))
  lv <- sort(unique(y[!is.na(y)]))
  as.integer(y == lv[2L])
}

# Rank predictors on one half, fit on the other. Ranking uses the univariate
# association of each predictor with the outcome; only the ranking sees the
# selection half, so the fitted half stays untouched by the choice.
.om_split <- function(d, outcome, predictors, fam, n_select, seed) {
  if (!is.null(seed)) .local_seed(seed)
  n <- nrow(d)
  sel_rows <- sample.int(n, size = floor(n / 2))
  fit_rows <- setdiff(seq_len(n), sel_rows)
  sel <- d[sel_rows, , drop = FALSE]

  strength <- vapply(predictors, function(p) {
    f <- stats::as.formula(paste(outcome, "~", p))
    m <- stats::glm(f, data = sel,
                    family = if (fam == "binomial") stats::binomial() else stats::gaussian())
    cf <- stats::coef(summary(m))
    if (nrow(cf) < 2L) return(0)
    abs(cf[2L, 3L])                       # |z| or |t|
  }, numeric(1L))

  ord <- order(strength, decreasing = TRUE)
  list(selected = predictors[utils::head(ord, as.integer(n_select))],
       fit_rows = fit_rows)
}

.om_fit <- function(d, outcome, predictors, adjust, group, fam) {
  rhs <- paste(c(adjust, predictors), collapse = " + ")
  mixed <- !is.null(group) && requireNamespace("lme4", quietly = TRUE)
  if (!is.null(group) && !mixed) {
    warning(warningCondition(
      paste0("lme4 is not installed, so the random intercept was dropped ",
             "and a plain glm was fitted instead. Install lme4 to model ",
             "the grouping."),
      class = "nestimate_no_lme4"))
  }
  if (mixed) {
    f <- stats::as.formula(sprintf("%s ~ %s + (1 | %s)", outcome, rhs, group))
    m <- if (fam == "binomial") {
      lme4::glmer(f, data = d, family = stats::binomial())
    } else {
      lme4::lmer(f, data = d)
    }
    return(list(model = m, mixed = TRUE, n_groups = nlevels(factor(d[[group]]))))
  }
  f <- stats::as.formula(sprintf("%s ~ %s", outcome, rhs))
  m <- stats::glm(f, data = d,
                  family = if (fam == "binomial") stats::binomial() else stats::gaussian())
  list(model = m, mixed = FALSE, n_groups = NA_integer_)
}

# Tidy coefficient table: estimate, Wald interval, corrected p-value, and for
# a binomial fit the odds ratio on its own scale.
.om_effects <- function(model, fam, ci_level, correction, mixed) {
  cf <- if (mixed) {
    summary(model)$coefficients
  } else {
    stats::coef(summary(model))
  }
  est <- cf[, 1L]
  se  <- cf[, 2L]
  stat <- cf[, 3L]
  p <- if (ncol(cf) >= 4L) {
    cf[, 4L]
  } else {
    2 * stats::pnorm(abs(stat), lower.tail = FALSE)   # lmer gives no p column
  }
  z <- if (!mixed && identical(fam, "gaussian")) {
    stats::qt(1 - (1 - ci_level) / 2, df = stats::df.residual(model))
  } else {
    stats::qnorm(1 - (1 - ci_level) / 2)
  }

  out <- data.frame(
    term      = rownames(cf),
    estimate  = as.numeric(est),
    std_error = as.numeric(se),
    statistic = as.numeric(stat),
    ci_lower  = as.numeric(est - z * se),
    ci_upper  = as.numeric(est + z * se),
    p_value   = as.numeric(p),
    stringsAsFactors = FALSE
  )
  is_term <- out$term != "(Intercept)"
  out$p_adj <- NA_real_
  out$p_adj[is_term] <- stats::p.adjust(out$p_value[is_term], method = correction)

  if (identical(fam, "binomial")) {
    out$odds_ratio <- exp(out$estimate)
    out$or_lower   <- exp(out$ci_lower)
    out$or_upper   <- exp(out$ci_upper)
  }
  rownames(out) <- NULL
  out
}

# ---- accessor + methods -----------------------------------------------------

#' Effect Table of a Fitted Outcome Model
#'
#' @description
#' `effects_table()` returns the effect table of an [outcome_model()] fit:
#' one row per term, with its estimate, confidence interval and corrected
#' p-value, rounded for reading. Arguments keep the intercept or only the
#' significant terms.
#'
#' @param x An [outcome_model()] result; anything else raises an error.
#' @param intercept Whether to keep the intercept row. Default `FALSE`.
#' @param significant Whether to keep only the terms whose corrected p-value
#'   (`p_adj`) is below `alpha`. Default `FALSE`.
#' @param alpha The threshold for `significant`, between 0 and 1. Default
#'   `0.05`.
#' @param digits The decimal places of the numeric columns. Default `3`;
#'   `p_value` and `p_adj` are never rounded.
#' @return A data frame with the columns of the model's effect table (`term`,
#'   `estimate`, `std_error`, `statistic`, `ci_lower`, `ci_upper`, `p_value`,
#'   `p_adj`, and for a binomial fit `odds_ratio`, `or_lower` and
#'   `or_upper`), one row per retained term. The intercept row is dropped
#'   unless `intercept = TRUE`, and every term is kept unless
#'   `significant = TRUE` restricts them to `p_adj < alpha`.
#' @examples
#' set.seed(1)
#' d <- data.frame(hint = rbinom(200, 1, 0.5))
#' d$success <- rbinom(200, 1, plogis(-0.3 + 0.9 * d$hint))
#' effects_table(outcome_model(d, outcome = "success", predictors = "hint"))
#' @export
effects_table <- function(x, intercept = FALSE, significant = FALSE,
                          alpha = 0.05, digits = 3) {
  .assert(
    "`x` must be a net_outcome_model" = inherits(x, "net_outcome_model"),
    "`alpha` must be between 0 and 1" =
      is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1
  )
  out <- x$effects
  if (!isTRUE(intercept)) out <- out[out$term != "(Intercept)", , drop = FALSE]
  if (isTRUE(significant)) {
    out <- out[!is.na(out$p_adj) & out$p_adj < alpha, , drop = FALSE]
  }
  num <- setdiff(names(out)[vapply(out, is.numeric, logical(1L))],
                 c("p_value", "p_adj"))
  out[num] <- lapply(out[num], round, digits = digits)
  rownames(out) <- NULL
  out
}

#' @rdname outcome_model
#' @param x For the `print()` and `plot()` methods: an `outcome_model()`
#'   result.
#' @param ... For the `print()`, `summary()` and `plot()` methods: not used.
#' @param object For `summary()`: an `outcome_model()` result.
#' @export
print.net_outcome_model <- function(x, ...) {
  cat("Outcome Model\n=============\n")
  cat("Outcome  :", x$outcome, "|", x$family, "\n")
  cat("Units    :", format(x$n, big.mark = ","),
      if (isTRUE(x$mixed)) sprintf("in %s groups", format(x$n_groups, big.mark = ",")) else "",
      "\n")
  cat("Terms    :", length(x$selected),
      if (identical(x$select, "split")) "(selected on a held-out half)" else "",
      "\n")
  if (length(x$adjust)) cat("Adjusted :", paste(x$adjust, collapse = ", "), "\n")
  if (length(x$dropped)) {
    cat("Dropped  :", paste(x$dropped, collapse = ", "), "(no variance)\n")
  }
  cat("p-values :", x$correction, "corrected\n\n")
  print(effects_table(x))
  invisible(x)
}

#' @rdname outcome_model
#' @export
summary.net_outcome_model <- function(object, ...) effects_table(object)

#' @rdname outcome_model
#' @export
plot.net_outcome_model <- function(x, ...) {
  d <- effects_table(x)
  ratio <- identical(x$family, "binomial")
  if (ratio) {
    d$estimate <- d$odds_ratio
    d$ci_lower <- d$or_lower
    d$ci_upper <- d$or_upper
  }
  d$term <- factor(d$term, levels = d$term[order(d$estimate)])
  d$sig  <- !is.na(d$p_adj) & d$p_adj < 0.05
  ref <- if (ratio) 1 else 0

  ggplot2::ggplot(d, ggplot2::aes(x = estimate, y = term)) +
    ggplot2::geom_vline(xintercept = ref, linetype = "dashed",
                        colour = "grey50") +
    # geom_errorbarh() is superseded and rewrites `height` to `width` with a
    # message; the horizontal orientation of geom_errorbar() is the current form.
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper), orientation = "y",
      width = 0.2, colour = "#0072B2") +
    ggplot2::geom_point(ggplot2::aes(shape = sig), size = 2.6,
                        colour = "#0072B2", fill = "#0072B2") +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 19),
                                name = sprintf("%s < 0.05", x$correction)) +
    ggplot2::labs(
      x = if (ratio) "Odds ratio (95% CI)" else "Estimate (95% CI)",
      y = NULL,
      title = sprintf("%s: %s", x$outcome, x$family)) +
    ggplot2::theme_minimal(base_size = 12)
}
