# Layer 6, adversarial data: the sequence verbs (sequence_metrics(), sequence_patterns()), the network comparison
# compare_model() and the multi-cluster network build_mcml() on the adversarial sequence sets. Every column of
# every table is compared.
# References: TraMineR 2.2 for sequence_metrics() (the reference of test-sequence_metrics.R, copied here; where a
# row has an internal gap, the columns TraMineR computes on joined runs are a known difference and every OTHER
# column of those rows is compared, which test-sequence_metrics.R leaves out); n-grams counted in base R for
# sequence_patterns() (occurrences, sequences, support, share of the n-grams of that length, lift = support over the
# product of the states' supports, the documented definition); tna::compare() for compare_model() on networks with
# zero cells; the sequences recoded to clusters (base R and tna) for build_mcml(), whose summary() is compared in
# all its columns (test-build_mcml.R leaves out within_total). build_mcml()'s initial probabilities follow its
# documented rule, the shares of the FIRST-COLUMN states with a missing cell left out (tna's rule), which differs
# from build_network()'s first OBSERVED state on rows whose first position is missing.
# codyna 1.0 (sequence_indices()) is the reference for the six indices TraMineR has no counterpart for; no other
# test compares them.
# Tolerances, fixed before comparing: counts exact; shares and probabilities 1e-12; TraMineR 1e-12 (as
# test-sequence_metrics.R); compare_model() 1e-10 (as test-compare_model.R).

sets <- adv_sequence_sets()

# ---- sequence_metrics() --------------------------------------------------------------------------------------

adv_traminer <- function(wide, favorable, omega) {
  alphabet <- sort(unique(stats::na.omit(unlist(wide, use.names = FALSE))))
  sequences <- suppressMessages(TraMineR::seqdef(wide, left = NA, gaps = NA, right = "DEL", alphabet = alphabet,
                                                 cpal = rep("#999999", length(alphabet))))
  indices <- TraMineR::seqindic(sequences, indic = c("lgth", "nonm", "visited", "transp", "meand", "entr"))
  shares <- TraMineR::seqistatd(sequences, prop = TRUE)
  spells <- as.matrix(TraMineR::seqdss(sequences))
  spells[spells %in% c("%", "*")] <- NA
  durations <- as.matrix(TraMineR::seqdur(sequences))
  rows <- seq_len(nrow(spells))
  dominant <- alphabet[max.col(shares, ties.method = "first")]
  variation <- apply(durations, 1L, \(d) {
    d <- d[!is.na(d)]
    if (length(d) == 1L) NA_real_ else min(stats::sd(d) / mean(d), 1)
  })
  data.frame(
    sequence = rows, valid_n = indices$NonM, valid_proportion = indices$NonM / indices$Lgth,
    unique_states = indices$Visited, mean_spell_duration = indices$MeanD,
    max_spell_duration = apply(durations, 1L, max, na.rm = TRUE), longitudinal_entropy = indices$Entr,
    simpson_diversity = 1 - rowSums(shares^2), self_loop_tendency = 1 - indices$Transp, transition_rate = indices$Transp,
    initial_state_persistence = durations[, 1L] / indices$Lgth,
    initial_state_proportion = shares[cbind(rows, match(spells[, 1L], alphabet))],
    first_state = spells[, 1L], last_state = spells[cbind(rows, rowSums(!is.na(spells)))],
    dominant_state = dominant, dominant_proportion = apply(shares, 1L, max),
    dominant_max_spell = vapply(rows, \(i) max(durations[i, which(spells[i, ] == dominant[i])]), numeric(1)),
    integrative_potential = as.vector(TraMineR::seqintegr(sequences, state = favorable, pow = omega)),
    complexity_index = ifelse(is.na(variation), 0, 0.4 * indices$Entr + 0.4 * indices$Transp + 0.2 * variation),
    stringsAsFactors = FALSE)
}

# Rows with a missing position between two observed ones.
adv_has_gap <- function(wide) {
  observed <- !is.na(as.matrix(wide))
  rowSums(observed) < max.col(observed, "last") - max.col(observed, "first") + 1L
}

gap_sensitive <- c("mean_spell_duration", "max_spell_duration", "self_loop_tendency", "transition_rate",
                   "initial_state_persistence", "dominant_max_spell", "complexity_index")

# Which rows of a wide table differ from TraMineR by definition, and in which columns.
adv_row_classes <- function(wide) {
  observed <- !is.na(as.matrix(wide))
  gap <- adv_has_gap(wide)
  list(gap = gap, leading = !gap & !observed[, 1L], single = !gap & observed[, 1L] & rowSums(observed) == 1L)
}

test_that("sequence_metrics() equals TraMineR on every adversarial set, every column", {
  skip_without_reference("TraMineR")
  reasons <- list(
    gap = paste("TraMineR (gaps = NA) joins the states before and after a missing position into one spell and one",
                "transition; sequence_metrics() ends a spell at a gap (documented)."),
    leading = paste("TraMineR (left = NA) counts the missing leading positions in the length, so its transition",
                    "rate is transitions / (length - 1) over all positions; sequence_metrics() divides by the observed",
                    "positions minus one (a missing position is not an observation, documented)."),
    single = paste("A sequence of one observed state has no transition: TraMineR gives a transition rate of 0 (and",
                   "self-loop tendency 1), sequence_metrics() NA (0 / 0)."))
  columns <- list(gap = gap_sensitive, leading = c("self_loop_tendency", "transition_rate", "complexity_index"),
                  single = c("self_loop_tendency", "transition_rate"))
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    favorable <- adv_states(wide)[[1L]]
    # a row with no observed state is no sequence for TraMineR: it is compared with the documented rule below
    nonempty <- rowSums(!is.na(wide)) > 0L
    found <- cbind(sequence = seq_len(nrow(wide)), sequence_metrics(wide, favorable = favorable))[nonempty, ]
    expected <- adv_traminer(wide[nonempty, , drop = FALSE], favorable, 1)
    expected$sequence <- which(nonempty)
    classes <- lapply(adv_row_classes(wide), \(rows) rows[nonempty])
    ordinary <- !Reduce(`|`, classes)
    expect_reference(found[ordinary, , drop = FALSE], expected[ordinary, , drop = FALSE], reference = "TraMineR",
                     covers = "sequence_metrics()", quantity = sprintf("ordinary rows, every column, %s", name),
                     tolerance = 1e-12, key = "sequence")
    invisible(lapply(names(classes), \(class) {
      rows <- classes[[class]]
      if (!any(rows)) return(NULL)
      kept <- setdiff(names(expected), columns[[class]])
      expect_reference(found[rows, kept, drop = FALSE], expected[rows, kept, drop = FALSE], reference = "TraMineR",
                       covers = "sequence_metrics()",
                       quantity = sprintf("%s rows, the columns the definitions share, %s", class, name),
                       tolerance = 1e-12, key = "sequence")
      expect_known_difference(found[rows, c("sequence", columns[[class]]), drop = FALSE],
                              expected[rows, c("sequence", columns[[class]]), drop = FALSE], reference = "TraMineR",
                              covers = "sequence_metrics()",
                              quantity = sprintf("%s rows, the columns the definitions differ on, %s", class, name),
                              tolerance = 1e-12, key = "sequence", reason = reasons[[class]])
    }))
  }))
})

test_that("a row with no observed state has valid_n and unique_states 0 and every other index NA", {
  found <- sequence_metrics(sets$empty_row, favorable = "A")
  empty <- which(rowSums(!is.na(sets$empty_row)) == 0L)
  row <- found[empty, , drop = FALSE]
  numeric_columns <- names(row)[vapply(row, is.numeric, logical(1))]
  expected <- stats::setNames(ifelse(numeric_columns %in% c("valid_n", "unique_states"), 0, NA_real_), numeric_columns)
  expect_reference(unlist(row[numeric_columns]), expected,
                   reference = "independent: the documented rule for a sequence with no observed state",
                   covers = "sequence_metrics()", quantity = "numeric indices of an empty row", tolerance = 0)
  expect_true(all(is.na(unlist(row[setdiff(names(row), numeric_columns)]))), info = "state columns of an empty row are NA")
})

# The indices TraMineR has no counterpart for: codyna 1.0's sequence_indices(), the implementation sequence_metrics()
# ports (codyna reads position 1 as the first state and NA as a value, so rows with a missing first position, an
# internal gap or no observed state are compared separately).
codyna_only <- c("transition_complexity", "initial_state_influence_decay", "cyclic_feedback_strength",
                 "emergent_state", "emergent_state_persistence", "emergent_state_proportion")

test_that("the six indices without a TraMineR counterpart equal codyna on rows without a leading or internal gap", {
  skip_without_reference("codyna")
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    observed <- !is.na(as.matrix(wide))
    plain <- observed[, 1L] & !adv_has_gap(wide)
    # codyna stops on a table holding a row with a missing first position, so both sides read the plain rows only
    rows <- wide[plain, , drop = FALSE]
    found <- cbind(sequence = which(plain), sequence_metrics(rows))[c("sequence", codyna_only)]
    reference <- cbind(sequence = which(plain), as.data.frame(codyna::sequence_indices(rows)))
    expect_reference(found, reference[c("sequence", codyna_only)], reference = "codyna",
                     authority = "open-source", covers = "sequence_metrics()",
                     quantity = sprintf("codyna-only indices, rows from position 1 without a gap, %s", name),
                     tolerance = 1e-12, key = "sequence")
  }))
})

# ---- sequence_patterns() -------------------------------------------------------------------------------------

# Contiguous n-grams of observed states (a gap ends a run), with every column sequence_patterns() returns.
adv_patterns <- function(wide, lens) {
  m <- as.matrix(wide)
  occurrences <- do.call(rbind, lapply(lens, \(len) {
    starts <- seq_len(max(ncol(m) - len + 1L, 0L))
    do.call(rbind, c(list(data.frame(sequence = integer(), pattern = character(), length = integer())),
                     lapply(starts, \(s) {
                       window <- m[, s:(s + len - 1L), drop = FALSE]
                       keep <- rowSums(is.na(window)) == 0L
                       data.frame(sequence = which(keep), pattern = apply(window[keep, , drop = FALSE], 1L, paste,
                                                                          collapse = "->"),
                                  length = rep(len, sum(keep)), stringsAsFactors = FALSE)
                     })))
  }))
  n <- nrow(m)
  state_support <- vapply(adv_states(wide), \(s) mean(rowSums(m == s, na.rm = TRUE) > 0), numeric(1))
  per_pattern <- split(occurrences, occurrences$pattern)
  totals <- table(occurrences$length)
  rows <- lapply(per_pattern, \(x) {
    len <- x$length[[1L]]
    support <- length(unique(x$sequence)) / n
    parts <- strsplit(x$pattern[[1L]], "->", fixed = TRUE)[[1L]]
    data.frame(pattern = x$pattern[[1L]], length = len, occurrences = nrow(x),
               proportion = nrow(x) / as.numeric(totals[[as.character(len)]]),
               n_sequences = length(unique(x$sequence)), support = support,
               lift = support / prod(state_support[parts]), stringsAsFactors = FALSE)
  })
  do.call(rbind, unname(rows))
}

test_that("sequence_patterns() n-grams of lengths 2 and 3 equal a count of the observed runs, every column", {
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    expect_reference(sequence_patterns(wide, type = "ngram", len = 2:3, min_freq = 1, min_support = 0),
                     adv_patterns(wide, 2:3), reference = "independent: contiguous n-grams of observed states, base R",
                     covers = "sequence_patterns()", quantity = sprintf("n-grams of lengths 2 and 3, %s", name),
                     tolerance = 1e-12, key = "pattern")
  }))
})

# ---- compare_model() -------------------------------------------------------------------------------------------

test_that("compare_model() of networks with zero cells equals tna::compare(), every table", {
  skip_without_reference("tna")
  # tna::compare() before 1.5 computed the distance correlation differently; the comparison is with the current tna
  skip_if(utils::packageVersion("tna") < "1.5.7", "needs tna >= 1.5.7 (the definition compared here)")
  pairs <- list(c("leading_na", "trailing_na"), c("identical", "ties"), c("one_actor", "short"),
                c("internal_na", "identical"))
  invisible(lapply(pairs, \(pair) {
    label <- paste(pair, collapse = " vs ")
    found <- compare_model(build_network(sets[[pair[[1L]]]], method = "relative"),
                           build_network(sets[[pair[[2L]]]], method = "relative"))
    reference <- tna::compare(tna::tna(sets[[pair[[1L]]]]), tna::tna(sets[[pair[[2L]]]]))
    expect_reference(found$summary_metrics, as.data.frame(reference$summary_metrics), reference = "tna",
                     covers = c("compare_model()", "compare_model.netobject()"),
                     quantity = sprintf("summary metrics, %s", label), tolerance = 1e-10, key = c("category", "metric"))
    # Three edge ratios are not defined where their denominator is 0 (tna Inf or NaN, Nestimate NA): those
    # cells are recorded as not comparable, every other cell of every column is compared.
    ratios <- c("difference_index", "similarity_strength_index", "standardized_score_inflation")
    edges_found <- found$edge_metrics
    edges_tna <- as.data.frame(reference$edge_metrics)
    expect_reference(edges_found, edges_tna[setdiff(names(edges_tna), ratios)], reference = "tna",
                     covers = c("compare_model()", "compare_model.netobject()"),
                     quantity = sprintf("edge metrics, %s", label), tolerance = 1e-10, key = c("source", "target"))
    long <- \(edges) do.call(rbind, lapply(ratios, \(ratio) {
      data.frame(edges[c("source", "target")], ratio = ratio, value = edges[[ratio]])
    }))
    ratios_found <- long(edges_found)
    ratios_tna <- long(edges_tna)
    ratio_key <- c("source", "target", "ratio")
    defined <- is.finite(ratios_tna$value)
    found_defined <- do.call(paste, ratios_found[ratio_key]) %in% do.call(paste, ratios_tna[defined, ratio_key])
    expect_reference(ratios_found[found_defined, ], ratios_tna[defined, ], reference = "tna",
                     covers = c("compare_model()", "compare_model.netobject()"),
                     quantity = sprintf("edge ratios where defined, %s", label), tolerance = 1e-10, key = ratio_key)
    if (any(!defined)) {
      expect_true(all(is.na(ratios_found$value[!found_defined])),
                  info = sprintf("an undefined ratio is NA, not Inf or NaN: %s", label))
      expect_not_comparable(ratios_found[!found_defined, ], ratios_tna[!defined, ], reference = "tna",
                            covers = "compare_model()", quantity = sprintf("edge ratios with a zero denominator, %s", label),
                            key = ratio_key,
                            reason = "A ratio with a zero denominator is not defined: tna returns Inf or NaN, Nestimate NA.")
    }
    expect_reference(found$network_metrics, as.data.frame(reference$network_metrics), reference = "tna",
                     covers = c("compare_model()", "compare_model.netobject()"),
                     quantity = sprintf("network metrics, %s", label), tolerance = 1e-10, key = "metric")
    expect_reference(found$difference_matrix, unclass(reference$difference_matrix), reference = "tna",
                     covers = c("compare_model()", "compare_model.netobject()"),
                     quantity = sprintf("difference matrix, %s", label), tolerance = 1e-10)
  }))
})

test_that("compare_model() of two constant matrices equals tna::compare(); undefined correlations on both sides", {
  skip_without_reference("tna")
  # the known difference below is against tna >= 1.5.7; older tna agrees with the NA by accident of its definition
  skip_if(utils::packageVersion("tna") < "1.5.7", "needs tna >= 1.5.7 (the definition compared here)")
  constant <- matrix(2, 3L, 3L, dimnames = list(letters[1:3], letters[1:3]))
  found <- compare_model(constant, 2 * constant)
  # tna's cor() warns "the standard deviation is zero" (diagnosed: the reference side only, a constant input)
  reference <- withCallingHandlers(tna::compare(constant, 2 * constant, scaling = "none"),
                                   warning = \(w) if (grepl("standard deviation is zero", conditionMessage(w)))
                                     invokeRestart("muffleWarning"))
  summary_tna <- as.data.frame(reference$summary_metrics)
  is_distance <- \(table) table$metric == "Distance"
  # a constant matrix has no variance: its correlations, CV ratio and RV coefficient are undefined; a value NA or
  # NaN on both sides is recorded as not comparable, every other value is compared
  expect_reference(found$summary_metrics[!is_distance(found$summary_metrics), ], summary_tna[!is_distance(summary_tna), ],
                   reference = "tna", covers = c("compare_model()", "compare_model.matrix()"),
                   quantity = "summary metrics but the distance correlation, constant matrices", tolerance = 1e-10,
                   key = c("category", "metric"))
  expect_known_difference(found$summary_metrics[is_distance(found$summary_metrics), ], summary_tna[is_distance(summary_tna), ],
                          reference = "tna", covers = "compare_model()",
                          quantity = "distance correlation, constant matrices", tolerance = 1e-10,
                          key = c("category", "metric"),
                          reason = paste("A zero distance variance: Nestimate returns NA on purpose",
                                         "(R/compare_network.R:671, 'like the Pearson branch'); tna returns 0, which is",
                                         "the value Szekely, Rizzo and Bakirov (2007, Definition 3) assign."))
  expect_reference(found$edge_metrics, as.data.frame(reference$edge_metrics), reference = "tna",
                   covers = c("compare_model()", "compare_model.matrix()"), quantity = "edge metrics, constant matrices",
                   tolerance = 1e-10, key = c("source", "target"))
})

# ---- build_mcml() ---------------------------------------------------------------------------------------------

# Shares of each level among `values`; all zero when there is none.
adv_first_shares <- function(values, levels) {
  counts <- as.numeric(table(factor(values, levels)))
  stats::setNames(if (sum(counts) > 0) counts / sum(counts) else counts, levels)
}

adv_partitions <- function(states) {
  half <- ceiling(length(states) / 2)
  list(two = list(g1 = states[seq_len(half)], g2 = states[-seq_len(half)]),
       singleton = list(alone = states[[1L]], rest = states[-1L]))
}

test_that("the macro and within-cluster networks and the summary of build_mcml() equal the recoded tables", {
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    states <- adv_states(wide)
    invisible(lapply(names(adv_partitions(states)), \(part) {
      clusters <- adv_partitions(states)[[part]]
      label <- sprintf("%s, %s partition", name, part)
      fit <- build_mcml(wide, clusters = clusters, type = "tna")
      recoded <- sq_recode(wide, clusters)
      macro <- adv_relative(adv_counts(recoded, names(clusters)))
      expect_reference(fit$macro$weights, macro,
                       reference = "independent: adjacent observed pairs of the recoded sequences, base R",
                       covers = "build_mcml()", quantity = sprintf("macro weights, %s", label), tolerance = 1e-12)
      # the documented rule (tna's): the shares of the FIRST-COLUMN states, a missing first cell left out
      first <- wide[[1L]]
      first <- first[!is.na(first)]
      expect_reference(fit$macro$inits, adv_first_shares(sq_recode(data.frame(f = first), clusters)$f, names(clusters)),
                       reference = "independent: shares of the first-column clusters, missing cells left out, base R",
                       covers = "build_mcml()", quantity = sprintf("macro initial probabilities, %s", label),
                       tolerance = 1e-12)
      within <- lapply(clusters, \(members) adv_relative(adv_counts(wide, states)[members, members, drop = FALSE]))
      invisible(lapply(names(clusters), \(cluster) {
        expect_reference(fit$clusters[[cluster]]$weights, within[[cluster]],
                         reference = "independent: adjacent observed pairs within the cluster, base R",
                         covers = "build_mcml()", quantity = sprintf("within weights of %s, %s", cluster, label),
                         tolerance = 1e-12)
        expect_reference(fit$clusters[[cluster]]$inits,
                         adv_first_shares(first[first %in% clusters[[cluster]]], clusters[[cluster]]),
                         reference = "independent: shares of the cluster's first-column states, base R",
                         covers = "build_mcml()", quantity = sprintf("within initial probabilities of %s, %s", cluster, label),
                         tolerance = 1e-12)
      }))
      expect_reference(summary(fit),
                       data.frame(cluster = names(clusters), size = lengths(clusters),
                                  within_total = vapply(within, sum, numeric(1)),
                                  between_out = rowSums(macro) - diag(macro), between_in = colSums(macro) - diag(macro)),
                       reference = "independent: adjacent observed pairs of the recoded sequences, base R",
                       covers = c("build_mcml()", "summary.mcml()"), quantity = sprintf("cluster summary, %s", label),
                       tolerance = 1e-12, key = "cluster")
    }))
  }))
})

test_that("the macro network equals tna::tna() on the recoded sequences", {
  skip_without_reference("tna")
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    clusters <- adv_partitions(adv_states(wide))$two
    fit <- build_mcml(wide, clusters = clusters, type = "tna")
    reference <- tna::tna(sq_recode(wide, clusters))
    expect_reference(fit$macro$weights, unclass(reference$weights), reference = "tna", covers = "build_mcml()",
                     quantity = sprintf("macro weights, %s", name), tolerance = 1e-12)
    expect_reference(fit$macro$inits, reference$inits, reference = "tna", covers = "build_mcml()",
                     quantity = sprintf("macro initial probabilities, %s", name), tolerance = 1e-12)
  }))
})
