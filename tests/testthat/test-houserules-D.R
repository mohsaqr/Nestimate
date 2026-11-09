# House-rule fixes in the plotting, sequence-comparison, co-occurrence and
# association-rule files: classed conditions, warnings for ignored arguments,
# corrected labels and printed text, and loops replaced by vectorised code
# (checked against values computed independently here).

.hrd_net <- function() {
  build_network(group_regulation_long, method = "relative",
                actor = "Actor", action = "Action", time = "Time")
}

.hrd_null_device <- function(env = parent.frame()) {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  do.call(on.exit, list(bquote(grDevices::dev.off(.(dev))), add = TRUE),
          envir = env)
}

# Text drawn by a base-graphics call, read from an uncompressed PDF.
.hrd_drawn_text <- function(code) {
  path <- tempfile(fileext = ".pdf")
  on.exit(unlink(path), add = TRUE)
  grDevices::pdf(path, compress = FALSE, useKerning = FALSE)
  dev <- grDevices::dev.cur()
  tryCatch(force(code), finally = grDevices::dev.off(dev))
  lines <- readLines(path, warn = FALSE)
  paste(lines[grepl("Tj$", lines, useBytes = TRUE)], collapse = "\n")
}

# ---- plot_state_frequencies.R ------------------------------------------------

test_that("classic mosaic_plot() names an argument it ignores", {
  tab <- table(group_regulation_long$Course, group_regulation_long$Achiever)
  plain <- mosaic_plot(tab, residuals = "asymptotic")
  expect_warning(
    typo <- mosaic_plot(tab, residuals = "asymptotic", valeus = TRUE),
    class = "nestimate_ignored_argument")
  expect_identical(typo$data, plain$data)
  expect_warning(mosaic_plot(tab, residuals = "asymptotic", valeus = TRUE),
                 "valeus")
  expect_no_warning(mosaic_plot(tab, residuals = "asymptotic", values = TRUE))
})

test_that("plot_state_frequencies() names an argument it ignores", {
  net <- .hrd_net()
  plain <- plot_state_frequencies(net)
  expect_warning(typo <- plot_state_frequencies(net, lable = "freq"),
                 class = "nestimate_ignored_argument")
  expect_identical(state_distribution(net), plain$table)
  expect_identical(typo$table, plain$table)
})

test_that("mosaic_plot() refuses a `fill` that names no column", {
  counts <- as.data.frame(
    table(group = group_regulation_long$Course,
          state = group_regulation_long$Achiever),
    responseName = "count", stringsAsFactors = FALSE)
  expect_error(
    mosaic_plot(counts, x_var = "group", y_var = "state", weight = "count",
                fill = "nope"),
    class = "nestimate_bad_input")
  expect_s3_class(
    mosaic_plot(counts, x_var = "group", y_var = "state", weight = "count"),
    "ggplot")
})

test_that("a one-state cluster without self-transitions is skipped, not fatal", {
  # No bundled data has a state that never follows itself, so this one frame
  # is written out: D never follows D.
  sequences <- data.frame(T1 = c("A", "B", "C", "A", "D"),
                          T2 = c("B", "C", "A", "D", "A"),
                          T3 = c("C", "A", "B", "B", "B"))
  fit <- build_mcml(sequences, clusters = list(G1 = c("A", "B", "C"), G2 = "D"))
  expect_warning(
    p <- mosaic_plot(fit, level = "clusters", residuals = "asymptotic"),
    class = "nestimate_skipped_cluster")
  expect_s3_class(p, "ggplot")
})

test_that("mosaic_plot(seed =) leaves the session's random stream alone", {
  tab <- table(group_regulation_long$Course, group_regulation_long$Achiever)
  set.seed(1)
  expected <- runif(1)
  set.seed(1)
  mosaic_plot(tab, seed = 5, n_perm = 20L)
  expect_identical(runif(1), expected)
})

# ---- sequence_plot.R / distribution_plot.R ------------------------------------

test_that("sequence_plot() argument checks are classed and worded", {
  .hrd_null_device()
  halves <- rep(c("a", "b"), length.out = nrow(trajectories))
  expect_error(sequence_plot(trajectories, k = 1),
               class = "nestimate_bad_input")
  expect_error(sequence_plot(trajectories, k = 1), "`k` must be")
  expect_error(sequence_plot(trajectories, type = "index", row_gap = 1),
               "`row_gap` must be", class = "nestimate_bad_input")
  expect_error(sequence_plot(trajectories, type = "index", group = 1:3),
               "`group` must have one value per sequence",
               class = "nestimate_bad_input")
  expect_error(sequence_plot(trajectories, type = "index", group = halves,
                             ncol = 1, nrow = 1),
               "must be at least the number of groups",
               class = "nestimate_bad_input")
})

test_that("distribution_plot() argument checks are classed and worded", {
  .hrd_null_device()
  halves <- rep(c("a", "b"), length.out = nrow(trajectories))
  expect_error(distribution_plot(trajectories, group = 1:3),
               "`group` must have one value per sequence",
               class = "nestimate_bad_input")
  expect_error(distribution_plot(trajectories, group = halves,
                                 ncol = 1, nrow = 1),
               "must be at least the number of groups",
               class = "nestimate_bad_input")
  expect_error(distribution_plot(trajectories, na = NA),
               "`na` must be TRUE or FALSE", class = "nestimate_bad_input")
})

test_that("panels of a clustering keep the \"Cluster\" title", {
  clustering <- build_clusters(trajectories, k = 2, seed = 1)
  index_text <- .hrd_drawn_text(sequence_plot(clustering, type = "index"))
  expect_match(index_text, "Cluster 1", fixed = TRUE)
  dist_text <- .hrd_drawn_text(distribution_plot(clustering))
  expect_match(dist_text, "Cluster 2", fixed = TRUE)
})

test_that("sequence_plot() names the arguments an mcml ignores", {
  expect_identical(.seqplot_mcml_ignored(c("x", "type", "trim", "main")),
                   character(0))
  expect_warning(.seqplot_mcml_ignored(c("x", "sort", "cols")),
                 "`sort`, `cols`", class = "nestimate_ignored_argument")
  fit <- build_mcml(trajectories,
                    clusters = list(On = c("Active", "Average"),
                                    Off = "Disengaged"))
  expect_warning(sequence_plot(fit, sort = "start"),
                 class = "nestimate_ignored_argument")
  expect_no_warning(sequence_plot(fit, type = "distribution"))
})

# ---- sequence_compare.R -------------------------------------------------------

.hrd_halves <- function() {
  rep(c("x", "y"), length.out = nrow(trajectories))
}

test_that("sequence_compare() argument checks are classed and worded", {
  expect_error(sequence_compare(trajectories, group = .hrd_halves(),
                                min_freq = 0),
               "`min_freq` must be", class = "nestimate_bad_input")
  expect_error(sequence_compare(trajectories, group = .hrd_halves(), sub = 0),
               "`sub` must be", class = "nestimate_bad_input")
  expect_error(sequence_compare(trajectories, group = .hrd_halves(), iter = 0),
               "`iter` must be", class = "nestimate_bad_input")
})

test_that("sequence_compare() warns when one permutation leaves no effect size", {
  expect_warning(
    res <- sequence_compare(trajectories, group = .hrd_halves(), sub = 2,
                            iter = 1),
    class = "nestimate_undefined_effect_size")
  expect_true(all(is.na(summary(res)$effect_size)))
  expect_no_warning(sequence_compare(trajectories, group = .hrd_halves(),
                                     sub = 2, iter = 2))
})

test_that("sequence_compare() says how many sequences have no group", {
  labels <- .hrd_halves()
  labels[1:10] <- NA
  expect_warning(
    res <- sequence_compare(trajectories, group = labels, sub = 2,
                            test = "none"),
    "10 of 138", class = "nestimate_missing_group")
  complete <- sequence_compare(trajectories[-(1:10), ], group = labels[-(1:10)],
                               sub = 2, test = "none")
  expect_identical(summary(res), summary(complete))
})

test_that("sequence_compare() residuals are the standardised residuals of each 2 x G table", {
  res <- summary(sequence_compare(trajectories, group = .hrd_halves(), sub = 2,
                                  min_freq = 1, test = "chisq",
                                  adjust = "none"))
  counts <- cbind(res$freq_x, res$freq_y)
  totals <- colSums(counts)
  reference <- t(vapply(seq_len(nrow(counts)), function(i) {
    two_by_g <- rbind(counts[i, ], totals - counts[i, ])
    fit <- suppressWarnings(stats::chisq.test(two_by_g))
    c(fit$stdres[1L, ], unname(fit$statistic), fit$p.value)
  }, numeric(4L)))
  expect_equal(res$resid_x, reference[, 1L], tolerance = 1e-12)
  expect_equal(res$resid_y, reference[, 2L], tolerance = 1e-12)
  expect_equal(res$statistic, reference[, 3L], tolerance = 1e-12)
  expect_equal(res$p_value, reference[, 4L], tolerance = 1e-12)
})

# ---- association_rules.R ------------------------------------------------------

test_that("association_rules() argument checks are classed and worded", {
  net <- build_network(trajectories, method = "relative")
  expect_error(association_rules(net, max_length = 1),
               "`max_length` must be", class = "nestimate_bad_input")
  expect_error(association_rules(net, min_support = 2),
               "`min_support` must be", class = "nestimate_bad_input")
  expect_error(association_rules(net, min_confidence = -1),
               "`min_confidence` must be", class = "nestimate_bad_input")
})

test_that("association_rules() warns that format = \"binary\" is ignored for a network", {
  net <- build_network(trajectories, method = "relative")
  plain <- association_rules(net, min_support = 0.2, min_confidence = 0.2,
                             min_lift = 0)
  expect_warning(
    declared <- association_rules(net, min_support = 0.2, min_confidence = 0.2,
                                  min_lift = 0, format = "binary"),
    class = "nestimate_ignored_argument")
  expect_identical(summary(declared), summary(plain))
})

test_that("association_rules() support, confidence and lift match direct counts", {
  net <- build_network(trajectories, method = "relative")
  rules <- summary(association_rules(net, min_support = 0.05,
                                     min_confidence = 0, min_lift = 0))
  expect_gt(nrow(rules), 0L)
  states <- c("Active", "Average", "Disengaged")
  present <- vapply(states, function(s) rowSums(trajectories == s,
                                                na.rm = TRUE) > 0,
                    logical(nrow(trajectories)))
  has_all <- function(items) {
    rowSums(present[, strsplit(items, ", ", fixed = TRUE)[[1L]],
                    drop = FALSE]) == lengths(strsplit(items, ", ",
                                                       fixed = TRUE))
  }
  reference <- t(vapply(seq_len(nrow(rules)), function(i) {
    lhs <- has_all(rules$antecedent[i])
    rhs <- has_all(rules$consequent[i])
    c(mean(lhs & rhs), mean(lhs & rhs) / mean(lhs),
      mean(lhs & rhs) / (mean(lhs) * mean(rhs)))
  }, numeric(3L)))
  expect_equal(rules$support, reference[, 1L], tolerance = 1e-12)
  expect_equal(rules$confidence, reference[, 2L], tolerance = 1e-12)
  expect_equal(rules$lift, reference[, 3L], tolerance = 1e-12)
  # Every pair and the triple, in both directions: 3 * 2 + 6 rules at most.
  expect_false(anyDuplicated(rules[c("antecedent", "consequent")]) > 0L)
})

test_that("print.net_association_rules() lists one numbered line per rule", {
  net <- build_network(trajectories, method = "relative")
  rules <- association_rules(net, min_support = 0.05, min_confidence = 0,
                             min_lift = 0)
  shown <- capture.output(print(rules))
  numbered <- grep("^    [0-9]+\\. ", shown, value = TRUE)
  expect_length(numbered, min(10L, rules$n_rules))
  top <- summary(rules)
  expect_identical(
    numbered[1L],
    sprintf("    1. %s -> %s  (sup=%.3f conf=%.3f lift=%.2f)",
            top$antecedent[1L], top$consequent[1L], top$support[1L],
            top$confidence[1L], top$lift[1L]))
})

# ---- co_occurrence.R ----------------------------------------------------------

.hrd_log <- function() {
  utils::head(group_regulation_long, 600L)
}

test_that("cooccurrence() names an argument it ignores", {
  plain <- cooccurrence(.hrd_log(), field = "Action", by = "Actor")
  expect_warning(
    typo <- cooccurrence(.hrd_log(), field = "Action", by = "Actor",
                         similarty = "jaccard"),
    "similarty", class = "nestimate_ignored_argument")
  expect_identical(typo$weights, plain$weights)
})

test_that("cooccurrence() warns that `by` is ignored beside `sep`", {
  rows <- cooccurrence(.hrd_log(), field = "Action", sep = ";")
  expect_warning(
    both <- cooccurrence(.hrd_log(), field = "Action", by = "Actor", sep = ";"),
    "`by` is ignored", class = "nestimate_ignored_argument")
  expect_identical(both$weights, rows$weights)
})

test_that("cooccurrence() refuses several item columns in a long table", {
  expect_error(
    cooccurrence(.hrd_log(), field = c("Action", "Achiever"), by = "Actor"),
    class = "nestimate_bad_input")
})

test_that("cooccurrence() argument checks are classed and worded", {
  expect_error(cooccurrence(.hrd_log(), field = "Action", by = "Actor",
                            threshold = -1),
               "`threshold` must be", class = "nestimate_bad_input")
  expect_error(cooccurrence(.hrd_log(), field = "Action", by = "Actor",
                            top_n = "a"),
               "`top_n` must be", class = "nestimate_bad_input")
  expect_error(cooccurrence(.hrd_log(), field = "Action", by = "Actor",
                            diagonal = "yes"),
               "`diagonal` must be", class = "nestimate_bad_input")
})

test_that("cooccurrence() counts equal the cross-product of the presence table", {
  net <- cooccurrence(.hrd_log(), field = "Action", by = "Actor")
  presence <- table(.hrd_log()$Actor, .hrd_log()$Action) > 0
  reference <- crossprod(presence * 1)
  weights <- net$weights
  expect_equal(unname(weights[rownames(reference), colnames(reference)]),
               unname(unclass(reference)), tolerance = 0)
  expect_true(isSymmetric(unname(weights)))
})

# ---- sequence_terminal.R ------------------------------------------------------

test_that("mark_*_state() conditions are classed", {
  expect_error(mark_terminal_state(trajectories, state = c("a", "b")),
               "`state` must be", class = "nestimate_bad_input")
  expect_error(mark_first_state(trajectories, state = ""),
               class = "nestimate_bad_input")
  expect_error(mark_terminal_state(1:3), "`data` must be",
               class = "nestimate_bad_input")
  expect_error(actor_endpoints(1:3), class = "nestimate_bad_input")
  expect_warning(marked <- mark_terminal_state(trajectories, state = "Active"),
                 class = "nestimate_state_collision")
  expect_identical(attr(marked, "terminal_state"), "Active_1")
  expect_warning(mark_first_state(trajectories, state = "Average"),
                 class = "nestimate_state_collision")
})

test_that("mark_terminal_state() fills exactly the cells after the last state", {
  original <- as.matrix(trajectories)
  marked <- as.matrix(mark_terminal_state(trajectories))
  n_steps <- ncol(original)
  last_seen <- apply(original, 1L, function(r) {
    seen <- which(!is.na(r))
    if (length(seen)) max(seen) else NA_integer_
  })
  expect_true(any(last_seen < n_steps, na.rm = TRUE))
  after <- !is.na(last_seen)[row(original)] & col(original) > last_seen[row(original)]
  expect_true(all(marked[after] == "End"))
  expect_identical(marked[!after], original[!after])
  expect_identical(sum(marked == "End", na.rm = TRUE), sum(after))
})

test_that("mark_first_state() fills exactly the cells before the first state", {
  # Reversing the columns turns the drop-out tails of `trajectories` into
  # leading gaps, which the bundled data does not have on its own.
  reversed <- trajectories[, rev(seq_len(ncol(trajectories)))]
  colnames(reversed) <- colnames(trajectories)
  original <- as.matrix(reversed)
  marked <- as.matrix(mark_first_state(reversed))
  first_seen <- apply(original, 1L, function(r) {
    seen <- which(!is.na(r))
    if (length(seen)) min(seen) else NA_integer_
  })
  before <- !is.na(first_seen)[row(original)] &
    col(original) < first_seen[row(original)]
  expect_gt(sum(before), 0L)
  expect_true(all(marked[before] == "Start"))
  expect_identical(marked[!before], original[!before])
  # The two markers mirror each other.
  mirrored <- as.matrix(mark_terminal_state(trajectories, state = "Start"))
  expect_identical(unname(marked),
                   unname(mirrored[, rev(seq_len(ncol(mirrored)))]))
})

# ---- state_colors.R / plot-utils.R ---------------------------------------------

test_that("set_state_colors() says an unmatched name stays in the palette", {
  net <- .hrd_net()
  expect_message(
    painted <- set_state_colors(net, c(plan = "#0072B2", Approve = "#D55E00")),
    "They stay in its palette")
  expect_message(set_state_colors(net, c(plan = "#0072B2", Approve = "#D55E00")),
                 "1 of 2 names are not drawn")
  expect_no_message(set_state_colors(net, c(plan = "#0072B2")))
  expect_message(set_state_colors(net, c(Approve = "#D55E00")),
                 "no name matches a state of this object")
  shown <- state_colors(painted)
  expect_identical(subset(shown, state == "plan")$color, "#0072B2")
  expect_false("Approve" %in% shown$state)
})

test_that(".encode_states() refuses anything but a frame or matrix by class", {
  expect_error(.encode_states(1:3), class = "nestimate_bad_input")
  encoded <- .encode_states(trajectories)
  expect_identical(encoded$levels, c("Active", "Average", "Disengaged"))
})

# ---- mosaic_analysis.R / magnitude_difference.R --------------------------------

test_that("print.mosaic_analysis() does not report an underflowed p-value as 0", {
  res <- mosaic_analysis(group_regulation_long, "Course", "Achiever")
  shown <- capture.output(print(res))
  expect_true(any(grepl("p < 2.2e-16", shown, fixed = TRUE)))
  expect_false(any(grepl("p = 0$", shown)))
  weak <- mosaic_analysis(group_regulation_long, "Achiever", "Action")
  expect_true(any(grepl("p = ", capture.output(print(weak)), fixed = TRUE)) ||
                any(grepl("p < 2.2e-16", capture.output(print(weak)),
                          fixed = TRUE)))
})

test_that("magnitude_difference() refuses non-frame data by class", {
  expect_error(magnitude_difference(1:3), "`data` must be a data frame",
               class = "nestimate_bad_input")
})
