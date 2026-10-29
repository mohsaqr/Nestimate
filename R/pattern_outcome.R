#' Model a Binary Outcome From Sequential Patterns
#'
#' Selects the \code{n} patterns most associated with a two-valued outcome
#' and fits one logistic regression of the outcome on their presence. This
#' reproduces \code{codyna::analyze_outcome()}: the same pattern discovery
#' and filters, the same ranking, and the same model, so the coefficients
#' agree. The result is a tidy table, one row per selected pattern.
#'
#' @details
#' Steps:
#' \enumerate{
#'   \item Patterns are found as in \code{\link{sequence_patterns}()} and
#'     filtered by \code{min_freq}, \code{min_support}, \code{start},
#'     \code{end} and \code{contain}. Patterns absent from all sequences of
#'     one outcome value are dropped.
#'   \item Each pattern gets the chi-square statistic of the number of
#'     sequences containing it under each outcome value, against an equal
#'     split across the values (codyna's statistic; it does not weight by
#'     the size of each outcome group).
#'   \item The top \code{n} patterns by \code{priority} enter one model,
#'     \code{outcome ~ adjust + pattern_1 + ... + pattern_n}, as presence
#'     indicators (or occurrence counts with \code{occurrences = TRUE}).
#'     With \code{group}, a random intercept per group is added
#'     (\code{lme4::glmer}).
#' }
#' The patterns are chosen by their association with the outcome in the
#' same data the model is fitted to, so the p-values do not account for
#' that selection and overstate the evidence. Treat them as a ranking, or
#' confirm the selected patterns on new data.
#'
#' @param x Sequence data, in any form accepted by
#'   \code{\link{extract_sequences}()}: prepared data, a long event log, a
#'   wide data frame, a \code{netobject}, \code{netobject_group} or
#'   \code{tna} model.
#' @param outcome The outcome, with exactly two distinct values: the name of
#'   a column stored with the sequences (a covariate of a wide data frame or
#'   long event log, or of a network's metadata), or a vector with one value
#'   per sequence. Use a vector for a model that stores no covariates, such
#'   as a \code{tna} model. \code{outcome = "last_obs"} uses each sequence's
#'   last observed state as the outcome and removes the states that end some
#'   sequence from all sequences (the default, as in codyna; for data whose
#'   final state codes the outcome, such as an appended \code{"pass"} or
#'   \code{"fail"}).
#' @param type,len,gap,min_freq,min_support,start,end,contain Pattern
#'   definition and filters; see \code{\link{sequence_patterns}()}. The
#'   defaults are codyna's: \code{len = 1:2}, \code{gap = 1},
#'   \code{min_freq = 5}, \code{min_support = 0.01}.
#' @param n Number of patterns entered into the model. Default \code{10}.
#' @param priority Column(s) ranking the patterns for selection, applied in
#'   order: \code{"chisq"} (default), \code{"p_value"}, \code{"support"},
#'   \code{"n_sequences"}, \code{"occurrences"}, \code{"proportion"},
#'   \code{"lift"}, \code{"length"}.
#' @param desc Logical, recycled to \code{priority}: rank in descending
#'   order. Default \code{TRUE}.
#' @param occurrences Logical. Enter each pattern as its number of
#'   occurrences in the sequence instead of presence (0/1). Default
#'   \code{FALSE}.
#' @param adjust Optional character vector of covariate columns added to
#'   the model.
#' @param adjust_length Logical. Also add the number of observed states per
#'   sequence as a covariate. Default \code{FALSE} (codyna does not).
#' @param group Optional column identifying the actor (or other cluster)
#'   each sequence belongs to; adds a random intercept. Requires \pkg{lme4}.
#' @param reference The reference outcome value; the odds ratios are for
#'   the other value. Default: the first value in sorted order.
#' @param ci_level Confidence level of the Wald intervals. Default
#'   \code{0.95}.
#' @param cols State columns of a wide data frame; see
#'   \code{\link{extract_sequences}()}.
#' @param ... Passed to \code{\link{extract_sequences}()}, e.g.
#'   \code{actor}, \code{action} and \code{time} for a long event log.
#'
#' @return A \code{data.frame} with one row per selected pattern, in the
#'   model's term order (pattern name): \code{pattern}, \code{length},
#'   \code{occurrences}, \code{n_sequences}, \code{support}, \code{lift}
#'   (as in \code{\link{sequence_patterns}()}), the support within each
#'   outcome value (\code{support_<value>}), \code{chisq} (the selection
#'   statistic),
#'   \code{estimate} (log odds ratio, adjusted for the other terms),
#'   \code{std_error}, \code{ci_lower}, \code{ci_upper}, \code{p_value}
#'   (Wald z), \code{odds_ratio}, \code{or_lower}, \code{or_upper}.
#'
#'   Conditions: \code{nestimate_bad_input} for a missing column, an outcome
#'   without exactly two values, or an unknown \code{reference};
#'   \code{nestimate_no_patterns} when no pattern passes the filters;
#'   \code{nestimate_no_lme4} when \code{group} is set and \pkg{lme4} is not
#'   installed; warning \code{nestimate_missing_outcome} when sequences with
#'   a missing outcome, covariate or group are dropped; warning
#'   \code{nestimate_aliased_pattern} when a selected pattern is collinear
#'   with the other terms (its row is \code{NA}). Warnings from the model
#'   fit (e.g. separation) are passed on unchanged.
#'
#' @references
#' Saqr, M., Tikka, S., & López-Pernas, S. (2026). codyna: Complex Dynamic
#' Systems. R package version 1.0.0.
#'
#' @examples
#' pattern_outcome(group_regulation_long, outcome = "Achiever",
#'                 actor = "Actor", action = "Action", time = "Time")
#'
#' @seealso \code{\link{sequence_patterns}}
#' @export
pattern_outcome <- function(x, outcome = "last_obs",
                            type = c("ngram", "gapped", "repeated"),
                            len = 1:2, gap = 1, min_freq = 5,
                            min_support = 0.01, start = NULL, end = NULL,
                            contain = NULL, n = 10, priority = "chisq",
                            desc = TRUE, occurrences = FALSE, adjust = NULL,
                            adjust_length = FALSE, group = NULL,
                            reference = NULL, ci_level = 0.95, cols = NULL,
                            ...) {
  type <- match.arg(type)
  priority <- match.arg(priority, c("chisq", "p_value", "support",
                                    "n_sequences", "occurrences",
                                    "proportion", "lift", "length"),
                        several.ok = TRUE)
  stopifnot(
    "`outcome` must be a column name or a vector of outcome values" =
      is.atomic(outcome) && length(outcome) >= 1L,
    "`n` must be a single positive whole number" = is.numeric(n) &&
      length(n) == 1L && n >= 1 && n == round(n),
    "`desc` must be logical" = is.logical(desc) && !anyNA(desc),
    "`occurrences` must be TRUE or FALSE" =
      isTRUE(occurrences) || isFALSE(occurrences),
    "`adjust` must be NULL or character" = is.null(adjust) ||
      is.character(adjust),
    "`adjust_length` must be TRUE or FALSE" =
      isTRUE(adjust_length) || isFALSE(adjust_length),
    "`group` must be NULL or a single column name" = is.null(group) ||
      (is.character(group) && length(group) == 1L),
    "`ci_level` must be a single number in (0, 1)" = is.numeric(ci_level) &&
      length(ci_level) == 1L && ci_level > 0 && ci_level < 1
  )
  .sp_check_args(len, gap, min_freq, min_support, start, end, contain)
  if (!is.null(group) && !requireNamespace("lme4", quietly = TRUE)) {
    stop(errorCondition("`group` needs the lme4 package; install it.",
                        class = "nestimate_no_lme4", call = NULL))
  }

  last_obs <- identical(outcome, "last_obs")
  outcome_is_column <- !last_obs && is.character(outcome) &&
    length(outcome) == 1L
  seqs <- .as_sequences(x, cols = cols,
                        exclude = c(if (outcome_is_column) outcome, adjust,
                                    group), ...)
  units <- seqs$units
  if (last_obs) {
    split <- .sp_last_obs(seqs$states)
    seqs$states <- split$states
    units$.outcome <- split$labels
    outcome <- ".outcome"
  } else if (!outcome_is_column) {
    if (length(outcome) != nrow(seqs$states)) {
      stop(errorCondition(
        sprintf("`outcome` has %d values; the data has %d sequences.",
                length(outcome), nrow(seqs$states)),
        class = "nestimate_bad_input", call = NULL))
    }
    units$.outcome <- outcome
    outcome <- ".outcome"
  }
  needed <- c(outcome, adjust, group)
  .po_check_columns(x, setdiff(needed, names(units)), names(units))
  complete <- stats::complete.cases(units[needed])
  if (!all(complete)) {
    warning(warningCondition(
      sprintf(paste("%d of %d sequences have a missing outcome, covariate",
                    "or group and are dropped."),
              sum(!complete), length(complete)),
      class = "nestimate_missing_outcome"))
  }
  states <- seqs$states[complete, , drop = FALSE]
  units <- units[complete, , drop = FALSE]
  y <- .po_outcome(units[[outcome]], reference)

  # ---- Discover, filter, rank and select (codyna's order of operations) ----
  found <- .sp_find(states, type = type, len = len, gap = gap)
  tab <- .sp_filter(found$table, min_freq = min_freq,
                    min_support = min_support, start = start, end = end,
                    contain = contain, sort = FALSE)
  occ <- found$occurrences
  # Sequences containing each pattern, per outcome value (patterns x values)
  by_level <- .sp_group_counts(found, tab$pattern, y$label, y$levels)
  in_all <- rowSums(by_level > 0) == ncol(by_level)
  tab <- tab[in_all, , drop = FALSE]
  by_level <- by_level[in_all, , drop = FALSE]
  if (nrow(tab) == 0L) {
    stop(errorCondition(
      "no pattern passes the filters under both outcome values.",
      class = "nestimate_no_patterns", call = NULL))
  }
  expected <- tab$n_sequences / ncol(by_level)
  tab$chisq <- rowSums((by_level - expected)^2 / expected)
  tab$p_value <- stats::pchisq(tab$chisq, ncol(by_level) - 1L,
                               lower.tail = FALSE)

  # Ties keep codyna's order: occurrences (desc), then first appearance
  base_order <- order(-tab$occurrences, method = "radix")
  tab <- tab[base_order, , drop = FALSE]
  by_level <- by_level[base_order, , drop = FALSE]
  desc <- rep_len(desc, length(priority))
  keys <- Map(\(col, d) if (d) -xtfrm(tab[[col]]) else xtfrm(tab[[col]]),
              priority, desc)
  ranked <- do.call(order, c(unname(keys), list(method = "radix")))
  chosen <- utils::head(ranked, n)
  chosen <- chosen[order(tab$pattern[chosen], method = "radix")]
  tab <- tab[chosen, , drop = FALSE]
  by_level <- by_level[chosen, , drop = FALSE]

  # ---- One model with every selected pattern ----
  sel_ids <- match(tab$pattern, found$table$pattern)
  # Internal names for every model column: a data column can then neither
  # collide with a pattern term nor break the formula (spaces, symbols)
  terms <- sprintf(".p%d", seq_along(sel_ids))
  adjust_terms <- sprintf(".a%d", seq_along(adjust))
  model_data <- units[c(adjust, group)]
  names(model_data) <- c(adjust_terms, if (!is.null(group)) ".g")
  model_data$.y <- y$value
  if (adjust_length) model_data$.length <- rowSums(!is.na(states))
  model_data[terms] <- lapply(sel_ids, \(pid) {
    rows <- occ$id == pid
    value <- integer(nrow(states))
    value[occ$seq[rows]] <- if (occurrences) occ$times[rows] else 1L
    value
  })
  rhs <- c(adjust_terms, if (adjust_length) ".length", terms,
           if (!is.null(group)) "(1 | .g)")
  form <- stats::reformulate(rhs, response = ".y")
  fit <- if (is.null(group)) {
    stats::glm(form, data = model_data, family = stats::binomial())
  } else {
    lme4::glmer(form, data = model_data, family = stats::binomial())
  }
  # An aliased pattern (e.g. present in every sequence) has no estimate:
  # glm returns NA, glmer drops the column; both become NA here
  all_coef <- if (is.null(group)) stats::coef(fit) else lme4::fixef(fit)
  coef_tab <- stats::coef(summary(fit))
  estimate <- unname(all_coef[terms])
  std_error <- unname(coef_tab[match(terms, rownames(coef_tab)), "Std. Error"])
  if (anyNA(estimate)) {
    warning(warningCondition(
      sprintf(paste("pattern(s) collinear with the other terms have no",
                    "estimate: %s."),
              paste(tab$pattern[is.na(estimate)], collapse = ", ")),
      class = "nestimate_aliased_pattern"))
  }

  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  support_by <- as.data.frame(
    by_level / matrix(table(factor(y$label, levels = y$levels)),
                      nrow(by_level), ncol(by_level), byrow = TRUE),
    check.names = FALSE)
  names(support_by) <- paste0("support_", y$levels)

  out <- data.frame(
    tab[c("pattern", "length", "occurrences", "n_sequences", "support",
          "lift")],
    support_by,
    chisq = tab$chisq,
    estimate = estimate,
    std_error = std_error,
    ci_lower = estimate - z * std_error,
    ci_upper = estimate + z * std_error,
    p_value = 2 * stats::pnorm(abs(estimate / std_error), lower.tail = FALSE),
    check.names = FALSE
  )
  out$odds_ratio <- exp(out$estimate)
  out$or_lower <- exp(out$ci_lower)
  out$or_upper <- exp(out$ci_upper)
  rownames(out) <- NULL
  out
}


#' Refuse missing covariate columns, saying why they are missing
#' @noRd
.po_check_columns <- function(x, absent, available) {
  if (length(absent) == 0L) return(invisible(TRUE))
  as_states <- inherits(x, "netobject") && is.data.frame(x$data) &&
    any(absent %in% names(x$data))
  hint <- if (as_states) {
    sprintf(paste("%s was read as a state column when the network was",
                  "built; rebuild it with `metadata_cols = %s`."),
            paste(intersect(absent, names(x$data)), collapse = ", "),
            deparse(intersect(absent, names(x$data))))
  } else if (length(available) == 0L) {
    paste("This input stores no covariates; pass `outcome` as a vector",
          "with one value per sequence.")
  } else {
    sprintf("Available: %s.", paste(available, collapse = ", "))
  }
  stop(errorCondition(
    sprintf("column(s) not found with the sequences: %s. %s",
            paste(absent, collapse = ", "), hint),
    class = "nestimate_bad_input", call = NULL))
}


#' Encode a two-valued outcome (1 = the non-reference value)
#' @noRd
.po_outcome <- function(v, reference = NULL) {
  values <- sort(unique(as.character(v)))
  if (length(values) != 2L) {
    stop(errorCondition(
      sprintf("`outcome` must have exactly two distinct values; it has %d.",
              length(values)),
      class = "nestimate_bad_input", call = NULL))
  }
  reference <- reference %||% values[1L]
  if (!reference %in% values) {
    stop(errorCondition(
      sprintf("`reference` must be one of: %s.",
              paste(values, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  label <- as.character(v)
  list(value = as.integer(label != reference), label = label,
       levels = c(reference, setdiff(values, reference)))
}
