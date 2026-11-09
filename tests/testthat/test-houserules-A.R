# House-rule fixes in transition_entropy.R, markov.R, markov_order.R,
# chain_structure.R, path_dependence.R and stability_landscape.R:
# classed conditions, clear errors in place of obscure ones, corrected labels.

.hra_net <- function() {
  build_network(as.data.frame(trajectories), method = "relative")
}

.hra_groups <- function() {
  build_network(group_regulation_long, action = "Action", actor = "Actor",
                time = "Time", group = "Achiever", method = "relative")
}

test_that("row normalisation warns with a class and the same result", {
  freq <- build_network(as.data.frame(trajectories), method = "frequency")
  expect_warning(cs <- chain_structure(freq),
                 class = "nestimate_rows_normalized")
  expect_warning(pt <- passage_time(freq), class = "nestimate_rows_normalized")
  expect_warning(te <- transition_entropy(freq),
                 class = "nestimate_rows_normalized")
  expect_equal(pt$matrix, passage_time(.hra_net())$matrix)
  expect_equal(rowSums(cs$P), c(Active = 1, Average = 1, Disengaged = 1))
})

test_that("passage-time plots refuse a bad log_scale or digits by class", {
  pt <- passage_time(.hra_net())
  expect_error(plot(pt, log_scale = "a"), class = "nestimate_bad_input")
  expect_error(plot(pt, digits = "a"), class = "nestimate_bad_input")
  expect_error(plot(passage_time(.hra_groups()), log_scale = NA),
               class = "nestimate_bad_input")
  expect_s3_class(plot(pt, log_scale = FALSE), "ggplot")
})

test_that("plot() of a grouped markov_stability warns that combined is ignored", {
  ms <- markov_stability(.hra_groups())
  expect_warning(p <- plot(ms, combined = FALSE),
                 class = "nestimate_ignored_argument")
  expect_s3_class(p, "ggplot")
  expect_no_warning(plot(ms))
})

test_that("group print methods print every group once and return the input", {
  groups <- .hra_groups()
  pt <- passage_time(groups)
  out <- capture.output(res <- print(pt))
  expect_identical(res, pt)
  expect_identical(sum(grepl("^--- (High|Low) ---$", out)), 2L)
  cs <- chain_structure(groups)
  out <- capture.output(res <- print(cs))
  expect_identical(res, cs)
  expect_identical(sum(grepl("^--- (High|Low) ---$", out)), 2L)
})

test_that("chain_structure hitting probabilities and periods are unchanged", {
  # Hand-computed: a absorbing; b <-> c closed with period 2; d transient,
  # leaving for b with probability 1 after a geometric stay.
  P <- matrix(0, 4, 4, dimnames = list(letters[1:4], letters[1:4]))
  P["a", "a"] <- 1
  P["b", "c"] <- 1
  P["c", "b"] <- 1
  P["d", "b"] <- 0.5
  P["d", "d"] <- 0.5
  cs <- chain_structure(P)
  expect_identical(unname(cs$period), c(1L, 2L, 2L, NA_integer_))
  expect_equal(unname(diag(cs$hitting_probabilities)), c(1, 1, 1, 0.5))
  expect_equal(cs$hitting_probabilities["d", "b"], 1)
  expect_equal(cs$hitting_probabilities["d", "a"], 0)
})

test_that("summary() of a chain_structure_group fills absent columns with NA", {
  irreducible <- chain_structure(.hra_net())
  P <- rbind(c(0.5, 0.5, 0), c(0.2, 0.3, 0.5), c(0, 0, 1))
  dimnames(P) <- list(c("a", "b", "c"), c("a", "b", "c"))
  mixed <- structure(list(one = irreducible, two = chain_structure(P)),
                     class = c("chain_structure_group", "list"))
  tab <- summary(mixed)
  expect_s3_class(tab, "data.frame")
  expect_identical(nrow(tab), 6L)
  expect_true(all(is.na(subset(tab, group == "one")$mean_absorption_time)))
  expect_true(all(is.na(subset(tab, group == "two")$stationary_probability)))
})

test_that("plot.chain_structure refuses a non-logical show_values by class", {
  cs <- chain_structure(.hra_net())
  expect_error(plot(cs, show_values = "a"), class = "nestimate_bad_input")
  expect_error(plot(cs, digits = "a"), class = "nestimate_bad_input")
  expect_s3_class(plot(cs, show_values = FALSE), "ggplot")
})

test_that("markov_order_test names a bad n_cores when it runs in parallel", {
  skip_on_os("windows")
  expect_error(
    markov_order_test(as.data.frame(trajectories), max_order = 1, n_perm = 9,
                      parallel = TRUE, n_cores = "x"),
    class = "nestimate_bad_input")
})

test_that("plot.net_markov_order says when there is no permutation null", {
  res <- markov_order_test(list(rep("A", 6), rep("A", 5)), max_order = 2,
                           n_perm = 9)
  expect_error(plot(res, panel = "permutation"), class = "nestimate_bad_input")
  expect_error(plot(res), class = "nestimate_bad_input")
  expect_s3_class(plot(res, panel = "ic"), "ggplot")
})

test_that("the permutation panel labels a non-rejection as not rejected", {
  res <- markov_order_test(as.data.frame(trajectories), max_order = 1,
                           n_perm = 9, seed = 1)
  p <- plot(res, panel = "permutation")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_identical(unique(as.character(p$layers[[2]]$data$decision)),
                   "do not reject H0")
})

test_that("path_dependence says when no context reaches min_count", {
  expect_error(path_dependence(as.data.frame(trajectories), min_count = 10000),
               class = "nestimate_bad_input")
  expect_error(path_dependence(as.data.frame(trajectories), min_count = 10000),
               "min_count")
})

test_that("a stream shorter than the window warns with a class", {
  short <- head(group_regulation_long, 300)
  expect_warning(
    tr <- entropy_trajectory(short, action = "Action", actor = "Actor",
                             window = 5000),
    class = "nestimate_short_stream")
  expect_identical(nrow(summary(tr)), 1L)
})

test_that("landscape_trend refuses single-state sequences by class", {
  one_state <- matrix("A", 10, 6, dimnames = list(NULL, paste0("T", 1:6)))
  expect_error(landscape_trend(one_state, iter = 0),
               class = "nestimate_bad_input")
})

test_that("landscape_trend labels the phase axis in one spelling", {
  lt <- landscape_trend(.hra_net(), iter = 0)
  expect_identical(plot(lt)$labels$x,
                   "Programme phase (centre of a 5-step window)")
  expect_match(capture.output(print(lt))[1], "programme phase", fixed = TRUE)
})
