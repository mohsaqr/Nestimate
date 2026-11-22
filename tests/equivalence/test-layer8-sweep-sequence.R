# Layer 8, the argument sweep of the sequence verbs: markov_order_test() at every dof and max_order against
# tna::markov_order() (deterministic columns), with the permutation p-value, the selected order at several alpha
# and the information-criterion orders checked against their documented rules; build_hypa() at every min_count,
# alpha and p_adjust against the documented labelling rule, several orders against one fit per order, and the
# deprecated `k`; sequence_patterns() at every type, length, gap and filter against patterns counted in base R.
# test-markov_order_test.R, test-build_hypa.R and test-sequence_patterns.R check the defaults; this file adds the
# options.

sq_wide <- tn_simulate(n_actors = 40L, n_states = 3L, seq_length = 14L, seed = 841L, ragged = TRUE)
# The same sequences with interior gaps: markov_order_test() splits a sequence at a missing value.
sq_gaps <- sq_wide
set.seed(842)
sq_gaps[cbind(sample(nrow(sq_gaps), 12L), sample(3:6, 12L, replace = TRUE))] <- NA

test_that("every dof and max_order of markov_order_test() equals tna::markov_order() on the deterministic columns", {
  skip_without_reference("tna")
  sw_skip_old_tna()
  cases <- expand.grid(input = c("ragged", "gaps"), dof = c("paths", "observed"), max_order = 1:3, stringsAsFactors = FALSE)
  inputs <- list(ragged = sq_wide, gaps = sq_gaps)
  invisible(Map(\(input, dof, max_order) {
    found <- summary(markov_order_test(inputs[[input]], max_order = max_order, n_perm = 19L, seed = 1L, dof = dof))
    reference <- tna::markov_order(tna::tna(inputs[[input]]), max_order = max_order, iter = 20L, dof = dof)$stats
    expected <- data.frame(order = reference$order, loglik = reference$loglik, AIC = reference$aic, BIC = reference$bic,
                           g2 = reference$g2, df = reference$df, p_asymptotic = reference$p_asymptotic)
    expect_reference(found, expected, reference = "tna", covers = c("markov_order_test()", "summary.net_markov_order()"),
                     quantity = sprintf("%s, dof %s, max_order %d", input, dof, max_order), tolerance = 1e-10,
                     measure = "relative", key = "order")
  }, cases$input, cases$dof, cases$max_order))
})

test_that("the permutation p-value, the selected order at every alpha and the criterion orders follow their rules", {
  fit_at <- \(alpha) markov_order_test(sq_wide, max_order = 3L, n_perm = 99L, seed = 3L, alpha = alpha)
  invisible(lapply(c(0.01, 0.05, 0.2, 0.5), \(alpha) {
    fit <- fit_at(alpha)
    table <- summary(fit)
    tested <- table[table$order >= 1L, , drop = FALSE]
    # (1 + B) / (1 + n_perm), B the permutations whose G2 is at least the observed one (?markov_order_test).
    nulls <- fit$permutation_null[sprintf("order_%d", tested$order)]
    expected_p <- vapply(seq_len(nrow(tested)), \(i) (1 + sum(nulls[[i]] >= tested$g2[[i]])) / (1 + 99), numeric(1))
    expect_reference(data.frame(order = tested$order, p_permutation = tested$p_permutation),
                     data.frame(order = tested$order, p_permutation = expected_p),
                     reference = "independent: (1 + B) / (1 + n_perm) from the stored permutation null",
                     covers = "markov_order_test()", quantity = sprintf("permutation p-values, alpha %s", alpha),
                     tolerance = 1e-12, key = "order")
    rejected <- !is.na(tested$p_permutation) & tested$p_permutation < alpha
    selected <- if (all(rejected)) max(tested$order) else tested$order[which(!rejected)[1L]] - 1L
    expect_reference(c(optimal = fit$optimal_order, aic = fit$aic_order, bic = fit$bic_order),
                     c(optimal = selected, aic = table$order[which.min(table$AIC)], bic = table$order[which.min(table$BIC)]),
                     reference = "independent: smallest order whose next test is not rejected; lowest AIC and BIC",
                     covers = "markov_order_test()", quantity = sprintf("selected orders, alpha %s", alpha), tolerance = 0)
  }))
  long <- tn_simulate(n_actors = 10L, n_states = 3L, seq_length = 4L, seed = 843L)
  expect_message(capped <- markov_order_test(long, max_order = 6L, n_perm = 5L, seed = 1L))
  expect_reference(capped$max_order, 3, reference = "independent: ?markov_order_test, max_order lowered to the longest length minus one",
                   covers = "markov_order_test()", quantity = "max_order capped", tolerance = 0)
  expect_reference(summary(markov_order_test(sq_wide, n_perm = 19L, seed = 2L)),
                   summary(markov_order_test(sq_wide, max_order = 3L, n_perm = 19L, alpha = 0.05, parallel = FALSE,
                                             n_cores = 2L, seed = 2L, dof = "paths")),
                   reference = "independent: the documented defaults written out", covers = "markov_order_test()",
                   quantity = "defaults", tolerance = 0, key = "order")
})

test_that("every min_count, alpha and p_adjust of build_hypa() labels the paths by the documented rule", {
  sequences <- tn_sequences(tn_simulate(n_actors = 60L, n_states = 4L, seq_length = 10L, seed = 844L, ragged = TRUE))
  cases <- expand.grid(min_count = c(1L, 3L, 8L), alpha = c(0.01, 0.05, 0.2), p_adjust = c("none", "BH", "holm", "bonferroni"),
                       order = 1:2, stringsAsFactors = FALSE)
  unadjusted <- lapply(1:2, \(o) build_hypa(sequences, order = o, min_count = 1L)$scores)
  invisible(Map(\(min_count, alpha, p_adjust, order) {
    found <- build_hypa(sequences, order = order, alpha = alpha, min_count = min_count, p_adjust = p_adjust)$scores
    base <- unadjusted[[order]]
    under <- if (identical(p_adjust, "none")) base$p_under else stats::p.adjust(base$p_under, p_adjust)
    over <- if (identical(p_adjust, "none")) base$p_over else stats::p.adjust(base$p_over, p_adjust)
    label <- ifelse(base$expected >= min_count & under < alpha, "under",
                    ifelse(base$observed >= min_count & over <= alpha, "over", "normal"))
    expect_reference(data.frame(path = found$path, anomaly = found$anomaly),
                     data.frame(path = base$path, anomaly = label),
                     reference = "independent: ?build_hypa labelling rule on stats::p.adjust() of the raw tails",
                     covers = "build_hypa()",
                     quantity = sprintf("labels, order %d, min_count %d, alpha %s, p_adjust %s", order, min_count, alpha,
                                        p_adjust),
                     tolerance = 0, key = "path")
  }, cases$min_count, cases$alpha, cases$p_adjust, cases$order))
  several <- build_hypa(sequences, order = c(3L, 1L, 2L), min_count = 2L)$scores
  invisible(lapply(1:3, \(o) {
    alone <- build_hypa(sequences, order = o, min_count = 2L)$scores
    expect_reference(several[several$order == o, c("path", "observed", "expected", "p_under", "p_over")],
                     alone[c("path", "observed", "expected", "p_under", "p_over")],
                     reference = "independent: the fit of that order alone", covers = "build_hypa()",
                     quantity = sprintf("order %d of order = c(3, 1, 2)", o), tolerance = 1e-12, key = "path")
  }))
  expect_warning(deprecated <- build_hypa(sequences, k = 2L, min_count = 2L), class = "deprecatedWarning")
  expect_reference(deprecated$scores[c("path", "observed", "expected")],
                   build_hypa(sequences, order = 2L, min_count = 2L)$scores[c("path", "observed", "expected")],
                   reference = "independent: ?build_hypa, k is the former name of order", covers = "build_hypa()",
                   quantity = "deprecated k", tolerance = 0, key = "path")
  expect_reference(build_hypa(sequences)$scores[c("path", "expected", "p_under")],
                   build_hypa(sequences, order = 2L, alpha = 0.05, min_count = 5L, p_adjust = "none",
                              k = NULL)$scores[c("path", "expected", "p_under")],
                   reference = "independent: the documented defaults written out", covers = "build_hypa()",
                   quantity = "defaults", tolerance = 0, key = "path")
})

# Patterns of one sequence set counted from the definition: contiguous runs of `len` states (ngram), two states
# `gap` positions apart with wildcards between them (gapped), runs of one repeated state (repeated); a window
# with a missing cell is skipped.
sq_patterns <- function(wide, type, size) {
  rows <- lapply(seq_len(nrow(wide)), \(i) {
    states <- unlist(wide[i, ], use.names = FALSE)
    span <- if (identical(type, "gapped")) size + 2L else size
    starts <- seq_len(max(0L, length(states) - span + 1L))
    grams <- vapply(starts, \(s) {
      window <- states[s:(s + span - 1L)]
      if (anyNA(window)) return(NA_character_)
      if (identical(type, "gapped")) {
        return(paste(c(window[1L], strrep("*", size), window[span]), collapse = "->"))
      }
      if (identical(type, "repeated") && length(unique(window)) > 1L) return(NA_character_)
      paste(window, collapse = "->")
    }, character(1))
    grams <- grams[!is.na(grams)]
    data.frame(sequence = rep(i, length(grams)), pattern = grams, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

sq_expected <- function(wide, type, sizes, min_freq, min_support, start = NULL, end = NULL, contain = NULL) {
  n <- nrow(wide)
  state_support <- vapply(sort(unique(stats::na.omit(unlist(wide)))), \(s) mean(apply(wide, 1L, \(r) s %in% r)),
                          numeric(1))
  tables <- lapply(sizes, \(size) {
    grams <- sq_patterns(wide, type, size)
    if (!nrow(grams)) return(NULL)
    occurrences <- table(grams$pattern)
    containing <- table(unique(grams)$pattern)
    states <- strsplit(names(occurrences), "->", fixed = TRUE)
    fixed <- lapply(states, \(s) s[!grepl("^\\*+$", s)])
    data.frame(pattern = names(occurrences), occurrences = as.vector(occurrences),
               proportion = as.vector(occurrences) / nrow(grams),
               n_sequences = as.vector(containing[names(occurrences)]),
               support = as.vector(containing[names(occurrences)]) / n,
               lift = as.vector(containing[names(occurrences)]) / n /
                 vapply(fixed, \(s) prod(state_support[s]), numeric(1)),
               first = vapply(fixed, \(s) s[1L], character(1)), last = vapply(fixed, \(s) s[length(s)], character(1)),
               stringsAsFactors = FALSE)
  })
  table <- do.call(rbind, tables)
  keep <- table$occurrences >= min_freq & table$support >= min_support &
    (is.null(start) | table$first %in% start) & (is.null(end) | table$last %in% end) &
    (is.null(contain) | vapply(strsplit(table$pattern, "->", fixed = TRUE), \(s) any(s %in% contain), logical(1)))
  table[keep, c("pattern", "occurrences", "proportion", "n_sequences", "support", "lift"), drop = FALSE]
}

test_that("every type, length, gap and filter of sequence_patterns() equals the patterns counted in base R", {
  wide <- tn_simulate(n_actors = 30L, n_states = 3L, seq_length = 10L, seed = 845L, ragged = TRUE)
  settings <- list(
    list(type = "ngram", sizes = 2L, min_freq = 1, min_support = 0),
    list(type = "ngram", sizes = 2:4, min_freq = 2, min_support = 0.01),
    list(type = "ngram", sizes = 3L, min_freq = 5, min_support = 0.2),
    list(type = "ngram", sizes = 2:3, min_freq = 1, min_support = 0, start = "A"),
    list(type = "ngram", sizes = 2:3, min_freq = 1, min_support = 0, end = c("B", "C")),
    list(type = "ngram", sizes = 3L, min_freq = 1, min_support = 0, contain = "C"),
    list(type = "gapped", sizes = 1L, min_freq = 1, min_support = 0),
    list(type = "gapped", sizes = 1:3, min_freq = 2, min_support = 0.05),
    list(type = "gapped", sizes = 2L, min_freq = 1, min_support = 0, start = "B", end = "A"),
    list(type = "repeated", sizes = 2:3, min_freq = 1, min_support = 0),
    list(type = "repeated", sizes = 2L, min_freq = 3, min_support = 0.1)
  )
  invisible(lapply(settings, \(s) {
    args <- c(list(wide, type = s$type, min_freq = s$min_freq, min_support = s$min_support),
              if (identical(s$type, "gapped")) list(gap = s$sizes) else list(len = s$sizes),
              s[intersect(names(s), c("start", "end", "contain"))])
    found <- do.call(sequence_patterns, args)
    label <- sprintf("%s, sizes %s, min_freq %s, min_support %s%s", s$type, paste(s$sizes, collapse = ":"), s$min_freq,
                     s$min_support, paste(sprintf(", %s %s", names(s)[names(s) %in% c("start", "end", "contain")],
                                                  vapply(s[intersect(names(s), c("start", "end", "contain"))],
                                                         paste, character(1), collapse = "/")), collapse = ""))
    expect_reference(found, sq_expected(wide, s$type, s$sizes, s$min_freq, s$min_support, s$start, s$end, s$contain),
                     reference = "independent: patterns counted in base R, then the documented filters",
                     covers = "sequence_patterns()", quantity = label, tolerance = 1e-12, key = "pattern")
  }))
  all_rows <- sequence_patterns(wide, len = 2:3, min_freq = 1, min_support = 0)
  ordered <- all_rows[order(-all_rows$support, -all_rows$occurrences, all_rows$pattern), , drop = FALSE]
  expect_reference(data.frame(rank = seq_len(nrow(all_rows)), pattern = all_rows$pattern),
                   data.frame(rank = seq_len(nrow(ordered)), pattern = ordered$pattern),
                   reference = "independent: ?sequence_patterns, ordered by support, then occurrences, then pattern",
                   covers = "sequence_patterns()", quantity = "row order", tolerance = 0, key = "rank")
  expect_reference(sequence_patterns(wide, len = 2:3, min_freq = 1, min_support = 0, top = 5L),
                   utils::head(all_rows, 5L), reference = "independent: ?sequence_patterns, top keeps the first rows",
                   covers = "sequence_patterns()", quantity = "top 5", tolerance = 0, key = "pattern")
  expect_reference(sequence_patterns(wide),
                   sequence_patterns(wide, type = "ngram", len = 2:5, gap = 1:3, min_freq = 2, min_support = 0.01,
                                     start = NULL, end = NULL, contain = NULL, pattern = NULL, group = NULL,
                                     reference = NULL, adjust = "BH", top = NULL),
                   reference = "independent: the documented defaults written out", covers = "sequence_patterns()",
                   quantity = "defaults", tolerance = 0, key = "pattern")
  expect_error(sequence_patterns(wide, pattern = "A->->B"), class = "nestimate_bad_input")
})
