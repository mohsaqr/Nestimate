# bayes_compare(): the Bayesian comparison of two transition networks under a row-wise Dirichlet posterior. The
# posterior means of each edge and their difference are compared with the closed form (counts + prior) / row
# total, the counts with a plain table of consecutive pairs, and the posterior probability of the more likely
# direction, p_difference = max(P(X > Y), P(X < Y)) of two independent Beta marginals, with its value by numerical
# integration of dbeta * pbeta. Re-expressed from local_testing_and_equivalence/test-equiv-bayes-compare.R, layers
# 1 and 3; its layers 2 and 4 replay the package's own random stream and stay there as a standing study.

# Counts of consecutive pairs of a wide table of sequences, rows = from, columns = to.
rs_pair_counts <- function(sequences) {
  sequences <- as.matrix(sequences)
  states <- sort(unique(as.vector(sequences)))
  unclass(table(factor(as.vector(sequences[, -ncol(sequences)]), states),
                factor(as.vector(sequences[, -1L]), states)))
}

# The Dirichlet posterior means of two count matrices, as one edge table.
rs_posterior <- function(counts_x, counts_y, prior) {
  alpha_x <- counts_x + prior
  alpha_y <- counts_y + prior
  mean_x <- alpha_x / rowSums(alpha_x)
  mean_y <- alpha_y / rowSums(alpha_y)
  rs_edge_table(count_x = counts_x, count_y = counts_y, weight_x = mean_x, weight_y = mean_y, diff = mean_x - mean_y)
}

configs <- list(list(states = 3L, n = 30L, length = 5L, prior = 0.5, seed = 11L),
                list(states = 5L, n = 60L, length = 7L, prior = 1, seed = 12L),
                list(states = 7L, n = 120L, length = 4L, prior = 0.5, seed = 13L))

test_that("posterior means and counts equal the closed-form Dirichlet mean", {
  invisible(lapply(seq_along(configs), \(i) {
    cfg <- configs[[i]]
    x <- rs_random_sequences(cfg$n, cfg$states, cfg$length, cfg$seed)
    y <- rs_random_sequences(cfg$n, cfg$states, cfg$length, cfg$seed + 100L)
    found <- bayes_compare(build_network(x, method = "relative"), build_network(y, method = "relative"),
                           prior = cfg$prior, draws = 200L, seed = 1)
    expected <- rs_posterior(rs_pair_counts(x), rs_pair_counts(y), cfg$prior)
    expect_reference(summary(found), expected, reference = "independent: Dirichlet posterior mean, closed form",
                     covers = c("bayes_compare()", "summary.net_bayes()"),
                     quantity = sprintf("config %d, posterior means", i), tolerance = 1e-12, key = c("from", "to"))
    expect_reference(extract_transition_matrix(as_netdifference(found, significant_only = FALSE)),
                     with(expected, matrix(diff, cfg$states, dimnames = list(sort(unique(from)), sort(unique(to))))),
                     reference = "independent: Dirichlet posterior mean, closed form",
                     covers = c("as_netdifference()", "as_netdifference.net_bayes()"),
                     quantity = sprintf("config %d, posterior mean difference as a netdifference", i),
                     tolerance = 1e-12)
  }))
})

test_that("p_difference equals P(X > Y) of the Beta marginals by numerical integration", {
  x <- rs_random_sequences(80L, 4L, 6L, 21L)
  y <- rs_random_sequences(80L, 4L, 6L, 22L)
  draws <- 20000L
  found <- bayes_compare(build_network(x, method = "relative"), build_network(y, method = "relative"),
                         draws = draws, seed = 7)
  alpha_x <- rs_pair_counts(x) + 0.5
  alpha_y <- rs_pair_counts(y) + 0.5
  beta_x <- rowSums(alpha_x) - alpha_x
  beta_y <- rowSums(alpha_y) - alpha_y
  # P(X > Y) = integral of f_X(t) F_Y(t) dt over [0, 1]: no code shared with the package.
  greater <- vapply(seq_along(alpha_x), \(e) stats::integrate(
    \(t) stats::dbeta(t, alpha_x[e], beta_x[e]) * stats::pbeta(t, alpha_y[e], beta_y[e]), 0, 1, rel.tol = 1e-9)$value,
    numeric(1))
  expected <- rs_edge_table(p_difference = matrix(pmax(greater, 1 - greater), nrow(alpha_x), dimnames = dimnames(alpha_x)))
  # A Monte Carlo estimate from 20000 draws: its standard error is at most sqrt(0.25 / 20000) = 0.0035, and 0.025
  # is seven standard errors, so a cell outside it is a wrong formula, not sampling noise.
  expect_reference(summary(found), expected, reference = "independent: Beta marginals, numerical integration",
                   covers = c("bayes_compare()", "summary.net_bayes()"), quantity = "p_difference, 20000 draws",
                   tolerance = 0.025, key = c("from", "to"))
})

test_that("a grouped network compares its two groups with the closed-form posterior means", {
  x <- rs_random_sequences(60L, 4L, 6L, 31L)
  group <- rep(c("g1", "g2"), length.out = nrow(x))
  found <- bayes_compare(build_network(cbind(x, grp = group), method = "relative", group = "grp"),
                         draws = 100L, seed = 1)
  expected <- rs_posterior(rs_pair_counts(x[group == "g1", ]), rs_pair_counts(x[group == "g2", ]), 0.5)
  expect_reference(summary(found), cbind(comparison = "g1 vs g2", expected),
                   reference = "independent: Dirichlet posterior mean, closed form",
                   covers = c("bayes_compare()", "summary.net_bayes_group()"), quantity = "two groups, posterior means",
                   tolerance = 1e-12, key = c("comparison", "from", "to"))
})

# Outside reference: the same posterior simulated with gtools::rdirichlet() (seeded Monte Carlo, three seeds). The
# tolerances are fixed before the run at 4 standard errors of the difference of two Monte Carlo estimates with
# `draws` draws each: a probability sqrt(2 * 0.25 / draws), its largest value; a mean 4 sd sqrt(2 / draws) with the
# largest posterior SD of the reference's draws.
test_that("posterior means and p_difference equal a gtools::rdirichlet() simulation of the posterior", {
  skip_without_reference("gtools")
  draws <- 20000L
  x <- rs_random_sequences(80L, 4L, 6L, 21L)
  y <- rs_random_sequences(80L, 4L, 6L, 22L)
  invisible(lapply(c(101L, 102L, 103L), \(seed) {
    sample_x <- out_dirichlet_draws(rs_pair_counts(x) + 0.5, draws, seed)
    sample_y <- out_dirichlet_draws(rs_pair_counts(y) + 0.5, draws, seed + 1000L)
    greater <- apply(sample_x > sample_y, c(2L, 3L), mean)
    dimnames(greater) <- dimnames(rs_pair_counts(x))
    mean_x <- apply(sample_x, c(2L, 3L), mean)
    mean_y <- apply(sample_y, c(2L, 3L), mean)
    dimnames(mean_x) <- dimnames(mean_y) <- dimnames(greater)
    sd_max <- max(apply(sample_x, c(2L, 3L), stats::sd), apply(sample_y, c(2L, 3L), stats::sd))
    found <- bayes_compare(build_network(x, method = "relative"), build_network(y, method = "relative"),
                           draws = draws, seed = seed)
    expect_reference(summary(found), rs_edge_table(p_difference = pmax(greater, 1 - greater)), reference = "gtools",
                     covers = c("bayes_compare()", "summary.net_bayes()"),
                     quantity = sprintf("p_difference, %d draws, reference seed %d, gtools", draws, seed),
                     tolerance = 4 * sqrt(2 * 0.25 / draws), key = c("from", "to"))
    expect_reference(summary(found), rs_edge_table(weight_x = mean_x, weight_y = mean_y, diff = mean_x - mean_y),
                     reference = "gtools", covers = c("bayes_compare()", "summary.net_bayes()"),
                     quantity = sprintf("posterior means, reference seed %d, gtools", seed),
                     tolerance = 2 * 4 * sd_max * sqrt(2 / draws), key = c("from", "to"),
                     note = "The difference of two means has twice the tolerance of one.")
    expect_reference(extract_transition_matrix(as_netdifference(found, significant_only = FALSE)), mean_x - mean_y,
                     reference = "gtools", covers = c("as_netdifference()", "as_netdifference.net_bayes()"),
                     quantity = sprintf("posterior mean difference as a netdifference, reference seed %d, gtools", seed),
                     tolerance = 2 * 4 * sd_max * sqrt(2 / draws))
  }))
  group <- rep(c("g1", "g2"), length.out = nrow(x))
  grouped <- bayes_compare(build_network(cbind(x, grp = group), method = "relative", group = "grp"), draws = draws,
                           seed = 5)
  sample_1 <- out_dirichlet_draws(rs_pair_counts(x[group == "g1", ]) + 0.5, draws, 201L)
  sample_2 <- out_dirichlet_draws(rs_pair_counts(x[group == "g2", ]) + 0.5, draws, 202L)
  greater <- apply(sample_1 > sample_2, c(2L, 3L), mean)
  dimnames(greater) <- dimnames(rs_pair_counts(x))
  expect_reference(summary(grouped), cbind(comparison = "g1 vs g2", rs_edge_table(p_difference = pmax(greater, 1 - greater))),
                   reference = "gtools", covers = c("bayes_compare()", "summary.net_bayes_group()"),
                   quantity = "two groups, p_difference, gtools", tolerance = 4 * sqrt(2 * 0.25 / draws),
                   key = c("comparison", "from", "to"))
})
