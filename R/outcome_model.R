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
#' Fits a regression of an outcome on predictor columns -- pattern
#' indicators, topological features from
#' \code{\link{simplicial_features}}, or any other numeric covariates --
#' and returns a tidy effect table with confidence intervals and
#' multiplicity-corrected p-values.
#'
#' @section Honest inference:
#' Choosing predictors by their association with the outcome and then
#' testing them on the same rows invalidates the p-values. With
#' \code{select = "split"} the data is halved: predictors are ranked on one
#' half and the reported model is fitted on the other, so the returned
#' inference is valid for the selected set. \code{select = "none"}
#' (default) fits every supplied predictor and needs no split.
#'
#' @param data A \code{data.frame} with one row per unit of analysis.
#' @param outcome Name of the outcome column. A two-valued outcome is
#'   modelled with a binomial family, a numeric one with gaussian; override
#'   with \code{family}.
#' @param predictors Character vector of predictor column names.
#' @param group Optional column name giving a grouping factor. When
#'   supplied and \pkg{lme4} is installed, a random intercept per group is
#'   added -- the right treatment for units nested in actors. \pkg{lme4} is
#'   a suggested package: when it is not installed the random intercept is
#'   dropped, a \code{"nestimate_no_lme4"} warning is raised, and a plain
#'   \code{\link[stats]{glm}} is fitted instead.
#' @param adjust Optional character vector of covariates entered before the
#'   predictors. Use it for exposure: a unit observed longer contains more
#'   of every pattern, so an unadjusted effect can be volume in disguise.
#' @param family \code{"auto"} (default), \code{"binomial"} or
#'   \code{"gaussian"}.
#' @param select \code{"none"} (default) fits all predictors;
#'   \code{"split"} ranks them on half the data and fits on the other half.
#' @param n_select Number of predictors kept when \code{select = "split"}.
#'   Default \code{10}.
#' @param correction Multiplicity correction passed to
#'   \code{\link[stats]{p.adjust}}. Default \code{"BH"}.
#' @param ci_level Confidence level for the intervals. Default \code{0.95}.
#' @param seed Optional integer seed for the split, so the result is
#'   reproducible. The RNG state is restored on exit.
#'
#' @return An object of class \code{net_outcome_model}: a list whose
#'   \code{$effects} element is the tidy \code{data.frame}, one row per
#'   model term (the intercept included), with columns \code{term},
#'   \code{estimate}, \code{std_error}, \code{statistic}, \code{ci_lower},
#'   \code{ci_upper}, \code{p_value} and \code{p_adj} (\code{NA} on the
#'   intercept row, which is excluded from the correction), plus
#'   \code{odds_ratio}, \code{or_lower} and \code{or_upper} for a binomial
#'   fit. The remaining elements are the fitted \code{$model} (a
#'   \code{glm}, or an \pkg{lme4} fit when a random intercept was added),
#'   \code{$family}, \code{$n} (rows the reported model was fitted on),
#'   \code{$n_groups} (\code{NA} unless mixed), \code{$selected},
#'   \code{$dropped} (zero-variance predictors), \code{$adjust},
#'   \code{$select}, \code{$correction}, \code{$ci_level}, \code{$mixed}
#'   and \code{$outcome}. Retrieve the table with
#'   \code{\link{effects_table}}.
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
  stopifnot(
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
    stop("Column(s) not found in `data`: ",
         paste(utils::head(missing_cols, 5L), collapse = ", "), call. = FALSE)
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
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv())) {
      old <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old, envir = globalenv()),
              add = TRUE, after = FALSE)
    }
    set.seed(seed)
  }
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
  z <- stats::qnorm(1 - (1 - ci_level) / 2)

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
#' The tidy one-row-per-term table of estimates, confidence intervals and
#' corrected p-values.
#'
#' @param x A \code{net_outcome_model} from \code{\link{outcome_model}}.
#' @param intercept Keep the intercept row? Default \code{FALSE}.
#' @param significant Keep only terms whose corrected p-value is below
#'   \code{alpha}? Default \code{FALSE}.
#' @param alpha Threshold used by \code{significant}. Default \code{0.05}.
#' @param digits Rounding for the numeric columns. Default \code{3};
#'   \code{p_value} and \code{p_adj} are never rounded.
#' @return A \code{data.frame} with the same columns as the model's effect
#'   table, one row per retained term. The intercept row is dropped unless
#'   \code{intercept = TRUE}, and every term is kept unless
#'   \code{significant = TRUE} restricts them to \code{p_adj < alpha}.
#' @examples
#' set.seed(1)
#' d <- data.frame(hint = rbinom(200, 1, 0.5))
#' d$success <- rbinom(200, 1, plogis(-0.3 + 0.9 * d$hint))
#' effects_table(outcome_model(d, outcome = "success", predictors = "hint"))
#' @export
effects_table <- function(x, intercept = FALSE, significant = FALSE,
                          alpha = 0.05, digits = 3) {
  stopifnot(
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

#' @param x A \code{net_outcome_model}, for the \code{print} and \code{plot}
#'   methods.
#' @param ... Unused.
#' @rdname outcome_model
#' @return \code{print} returns its input invisibly.
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
#' @param object A \code{net_outcome_model}.
#' @return \code{summary} returns the tidy effect table.
#' @export
summary.net_outcome_model <- function(object, ...) effects_table(object)

#' @rdname outcome_model
#' @return \code{plot} returns a \code{ggplot} forest of the effects.
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
