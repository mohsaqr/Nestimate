# ---- stability_landscape() / landscape_surface() / landscape_trend() ----

# A A B A C B C A -> AA 1, AB 1, AC 1, BA 1, BC 1, CA 1, CB 1.
# Stay probabilities: A 1/3, B 0, C 0. Occupancy: A 4/8, B 2/8, C 2/8.
toy <- data.frame(T1 = "A", T2 = "A", T3 = "B", T4 = "A", T5 = "C",
                  T6 = "B", T7 = "C", T8 = "A")

test_that("cup depth is the log of the hand-computed expected stay", {
  cups <- summary(stability_landscape(toy, iter = 0))
  expect_equal(cups$state, c("A", "B", "C"))
  expect_equal(cups$stay, c(1 / 3, 0, 0), tolerance = 1e-12)
  expect_equal(cups$dwell, c(1.5, 1, 1), tolerance = 1e-12)
  expect_equal(cups$depth, c(log(1.5), 0, 0), tolerance = 1e-12)
  expect_equal(cups$share, c(4, 2, 2) / 8, tolerance = 1e-12)
  expect_true(all(is.na(cups$ci_lower)))          # no resamples, no interval
  exits <- summary(stability_landscape(toy, iter = 0), table = "exits")
  expect_equal(exits$probability[exits$from == "A"], c(1, 1) / 3,
               tolerance = 1e-12)
})

test_that("depth follows the diagonal of build_network() and exp(depth) = dwell", {
  net <- build_network(trajectories, method = "relative")
  cups <- summary(stability_landscape(net, iter = 0))
  stay <- unname(diag(net$weights)[cups$state])
  expect_equal(cups$depth, -log(1 - stay), tolerance = 1e-12)
  expect_equal(exp(cups$depth), cups$dwell, tolerance = 1e-12)
})

test_that("three states: the circulation is the net flux on every edge", {
  set.seed(11)
  df <- as.data.frame(matrix(sample(c("A", "B", "C"), 30 * 8, replace = TRUE,
                                    prob = c(0.5, 0.3, 0.2)), nrow = 30))
  circ <- summary(stability_landscape(df, iter = 0), table = "circulation")
  P <- build_network(df, method = "relative")$weights[c("A", "B", "C"),
                                                     c("A", "B", "C")]
  pi <- passage_time(P)$stationary
  edge <- \(i, j) pi[i] * P[i, j] - pi[j] * P[j, i]
  expect_equal(circ$current, abs(unname(edge(1, 2))), tolerance = 1e-12)
  expect_equal(unname(edge(1, 2)), unname(edge(2, 3)), tolerance = 1e-12)
  expect_equal(unname(edge(1, 2)), unname(edge(3, 1)), tolerance = 1e-12)
  expect_error(summary(stability_landscape(
    data.frame(T1 = c("A", "B"), T2 = c("B", "A")), iter = 0),
    table = "circulation"), class = "nestimate_bad_input")
})

test_that("states sets the order and is validated", {
  land <- stability_landscape(toy, states = c("C", "A", "B"), iter = 0)
  expect_equal(summary(land)$state, c("C", "A", "B"))
  expect_error(stability_landscape(toy, states = c("A", "B"), iter = 0),
               class = "nestimate_bad_input")
  expect_error(stability_landscape(toy, states = c("A", "A", "B", "C"),
                                   iter = 0),
               class = "nestimate_bad_input")
})

test_that("non-transition networks and never-left states are refused", {
  set.seed(3)
  num <- as.data.frame(matrix(rnorm(60), ncol = 3))
  expect_error(stability_landscape(build_network(num, method = "cor")),
               class = "nestimate_not_markov")
  absorbing <- data.frame(T1 = c("A", "B"), T2 = c("B", "C"),
                          T3 = c("A", "C"), T4 = c("B", "C"))
  expect_error(stability_landscape(absorbing, iter = 0),
               class = "nestimate_landscape_unestimable")
})

test_that("grouped networks give one landscape per group", {
  df <- data.frame(T1 = c("A", "B", "A", "C", "A", "B"),
                   T2 = c("B", "A", "A", "A", "C", "C"),
                   T3 = c("A", "C", "B", "B", "B", "A"),
                   T4 = c("C", "B", "C", "A", "B", "C"),
                   grp = c("x", "x", "x", "y", "y", "y"))
  g <- build_network(df, method = "relative", group = "grp")
  cups <- summary(stability_landscape(g, iter = 0))
  expect_equal(unique(cups$group), c("x", "y"))
  expect_equal(nrow(cups), 6L)
})

test_that("seed makes the resampling reproducible and restores the RNG", {
  set.seed(99)
  before <- .Random.seed
  a <- summary(stability_landscape(trajectories, iter = 20, seed = 5))
  expect_identical(.Random.seed, before)
  b <- summary(stability_landscape(trajectories, iter = 20, seed = 5))
  expect_identical(a, b)
  expect_true(all(a$ci_lower <= a$depth & a$depth <= a$ci_upper))
})

test_that("a state left once can have an unbounded upper interval", {
  # C is left once: resamples without that sequence never leave C
  df <- data.frame(T1 = c("A", "A", "B", "C", "C"),
                   T2 = c("B", "A", "A", "C", "A"),
                   T3 = c("A", "B", "B", "C", "B"),
                   T4 = c("B", "A", "A", "C", "A"))
  cups <- summary(stability_landscape(df, iter = 200, seed = 1))
  expect_true(is.infinite(cups$ci_upper[cups$state == "C"]))
  expect_true(is.finite(cups$depth[cups$state == "C"]))
})

test_that("landscape_surface slices match an independent recomputation", {
  surf <- landscape_surface(trajectories, measure = "self_loop_tendency",
                            k = 20, n_grid = 5)
  full <- summary(surf, ends = FALSE)
  expect_equal(nrow(full), 5L * 3L)
  expect_true(all(full$effective_n <= 20 + 1e-9))
  met <- sequence_metrics(trajectories)
  keep <- which(met$valid_n >= 2)
  m <- met$self_loop_tendency[keep]
  y <- sort(unique(full$value))[3]
  d <- abs(m - y)
  h <- max(sort(d)[20], diff(range(m)) / 50) * 1.01
  w <- ifelse(d < h, (1 - (d / h)^3)^3, 0)
  st <- c("Active", "Average", "Disengaged")
  tabs <- lapply(keep, \(i) {
    v <- trajectories[i, ]
    v <- v[!is.na(v)]
    table(factor(utils::head(v, -1), st), factor(v[-1], st))
  })
  n <- Reduce(`+`, Map(`*`, tabs, w))
  slice <- full[abs(full$value - y) < 1e-12, ]
  expect_equal(slice$depth, unname(-log(1 - diag(n) / rowSums(n))),
               tolerance = 1e-10)
})

test_that("landscape_surface summary and errors", {
  surf <- landscape_surface(trajectories,
                            measure = c("self_loop_tendency",
                                        "longitudinal_entropy"),
                            n_grid = 4)
  ends <- summary(surf)
  expect_equal(nrow(ends), 2L * 3L)
  expect_true(all(ends$low_value < ends$high_value))
  expect_error(landscape_surface(trajectories, measure = "dominant_state"),
               class = "nestimate_bad_input")
  expect_error(landscape_surface(trajectories, measure = "no_such_measure"),
               class = "nestimate_bad_input")
})

test_that("landscape_trend over phase: one full-length window is the landscape", {
  net <- build_network(trajectories, method = "relative")
  one <- summary(landscape_trend(net, window = 15, iter = 0))
  expect_equal(one$depth, summary(stability_landscape(net, iter = 0))$depth,
               tolerance = 1e-12)
  tr <- summary(landscape_trend(net, window = 5, iter = 0))
  expect_equal(nrow(tr), 11L * 3L)
  expect_equal(unique(tr$label)[1], "1-5")
  # window 2 at position 1 is the first transition of every sequence
  w2 <- summary(landscape_trend(net, window = 2, iter = 0))
  first <- table(factor(trajectories[, 1], c("Active", "Average", "Disengaged")),
                 factor(trajectories[, 2], c("Active", "Average", "Disengaged")))
  expect_equal(w2$depth[w2$value == 1.5],
               unname(-log(1 - diag(first) / rowSums(first))), tolerance = 1e-12)
})

test_that("landscape_trend over a covariate matches an independent slice", {
  set.seed(4)
  df <- data.frame(trajectories, grade = round(stats::rnorm(138, 70, 10), 1))
  net <- build_network(df, method = "relative")
  tr <- summary(landscape_trend(net, along = "grade", k = 25, n_grid = 4,
                                iter = 0))
  expect_equal(nrow(tr), 4L * 3L)
  seqs <- as.matrix(df[, 1:15])
  keep <- which(rowSums(!is.na(seqs)) >= 2)
  m <- df$grade[keep]
  y <- sort(unique(tr$value))[2]
  d <- abs(m - y)
  h <- max(sort(d)[25], diff(range(m)) / 50) * 1.01
  w <- ifelse(d < h, (1 - (d / h)^3)^3, 0)
  st <- c("Active", "Average", "Disengaged")
  tabs <- lapply(keep, \(i) {
    v <- seqs[i, ]
    v <- v[!is.na(v)]
    table(factor(utils::head(v, -1), st), factor(v[-1], st))
  })
  n <- Reduce(`+`, Map(`*`, tabs, w))
  expect_equal(tr$depth[abs(tr$value - y) < 1e-12],
               unname(-log(1 - diag(n) / rowSums(n))), tolerance = 1e-10)
})

test_that("landscape_trend occupancy: cups -log pi and flux ridges", {
  net <- build_network(trajectories, method = "relative")
  st <- c("Disengaged", "Average", "Active")
  one <- summary(landscape_trend(net, depth = "occupancy", window = 15,
                                 states = st, iter = 0))
  P <- net$weights[st, st]
  pi <- passage_time(net, states = st)$stationary
  flux <- \(i, j) (pi[i] * P[i, j] + pi[j] * P[j, i]) / 2
  expect_equal(one$potential[one$feature == "well"], unname(-log(pi)),
               tolerance = 1e-10)
  expect_equal(one$potential[one$feature == "ridge"],
               unname(-log(c(flux(1, 2), flux(2, 3), flux(3, 1)))),
               tolerance = 1e-10)
  expect_equal(one$state[one$feature == "ridge"],
               c("Disengaged | Average", "Average | Active",
                 "Active | Disengaged"))
})

test_that("landscape_trend validates along and window, and scopes the seed", {
  net <- build_network(trajectories, method = "relative")
  expect_error(landscape_trend(net, along = "grade", iter = 0),
               class = "nestimate_bad_input")
  expect_error(landscape_trend(net, window = 16, iter = 0),
               class = "nestimate_bad_input")
  set.seed(7)
  before <- .Random.seed
  a <- summary(landscape_trend(net, iter = 10, seed = 3))
  expect_identical(.Random.seed, before)
  expect_identical(a, summary(landscape_trend(net, iter = 10, seed = 3)))
})

test_that("print and plot methods run", {
  skip_on_cran()
  land <- stability_landscape(trajectories,
                              states = c("Disengaged", "Average", "Active"),
                              iter = 20, seed = 1)
  expect_output(print(land), "Net circulation")
  p <- plot(land)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  surf <- landscape_surface(trajectories, measure = "self_loop_tendency",
                            n_grid = 4)
  expect_output(print(surf), "Cup depth")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(surf), surf)
  tr <- landscape_trend(trajectories, iter = 10, seed = 1)
  expect_output(print(tr), "Cup depth")
  expect_silent(ggplot2::ggplot_build(plot(tr)))
  occ <- landscape_trend(trajectories, depth = "occupancy", window = 6,
                         iter = 5, seed = 1)
  expect_identical(plot(occ, type = "surface"), occ)
})
