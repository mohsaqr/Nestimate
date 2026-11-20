# entropy_bayes(): the posterior of the entropy of a Markov chain under independent Dirichlet(counts + prior) rows.
# The draws are taken inside the function with only a seed, so the comparison is statistical (Monte Carlo
# summaries), declared as such:
# - the posterior mean and SD of each row's entropy against their closed forms (Wolpert & Wolf 1995; Hutter 2002),
#   exact for any number of draws;
# - the chain-level quantities (entropy rate, stationary entropy, redundancy) and the per-edge contributions
#   against an independent simulation of the same posterior written here in base R (rgamma rows, the stationary
#   distribution from eigen()), with its own draws.
# Tolerances are fixed before the run from the Monte Carlo error: 4 standard errors of the difference of the two
# estimates (a mean: sd / sqrt(draws) on each side; an SD: sd * sqrt((kurtosis - 1) / (4 draws)), the kurtosis read
# from the reference's draws). A seed sweep of 20 seeds before fixing these passed 20 of 20 for every summary.

eb_draws <- 20000L
eb_states <- c("A", "B", "C", "D")
eb_counts <- matrix(c(30, 12, 5, 3,
                      8, 25, 10, 2,
                      4, 9, 18, 11,
                      6, 3, 14, 20), 4L, 4L, byrow = TRUE, dimnames = list(eb_states, eb_states))

# Closed-form posterior mean and SD of the entropy (in bits) of p ~ Dirichlet(a): E[H] = psi(A + 1) -
# sum(a / A psi(a + 1)) (Wolpert & Wolf 1995, eq. 15) and E[H^2] (Hutter 2002, Theorem 4, for the entropy).
eb_dirichlet_entropy <- function(a) {
  total <- sum(a)
  mean_nats <- digamma(total + 1) - sum(a / total * digamma(a + 1))
  centred <- digamma(a + 1) - digamma(total + 2)
  same <- sum(a * (a + 1) / (total * (total + 1)) *
                ((digamma(a + 2) - digamma(total + 2))^2 + trigamma(a + 2) - trigamma(total + 2)))
  cross <- (sum(a * centred)^2 - sum((a * centred)^2)) / (total * (total + 1)) -
    (total^2 - sum(a^2)) / (total * (total + 1)) * trigamma(total + 2)
  c(mean = mean_nats / log(2), sd = sqrt(same + cross - mean_nats^2) / log(2))
}

# An independent simulation of the posterior: chain-level quantities and per-edge contributions, in bits.
eb_reference <- function(counts, prior, draws, seed) {
  set.seed(seed)
  n <- nrow(counts)
  one <- function() {
    gamma <- matrix(stats::rgamma(n * n, shape = counts + prior), n, n)
    p <- gamma / rowSums(gamma)
    leading <- eigen(t(p))$vectors[, 1L]
    stationary <- Re(leading) / sum(Re(leading))
    row_entropy <- -rowSums(p * log2(p))
    rate <- sum(stationary * row_entropy)
    stationary_entropy <- -sum(stationary * log2(stationary))
    c(entropy_rate = rate, stationary_entropy = stationary_entropy, redundancy = stationary_entropy - rate,
      as.vector(-stationary * p * log2(p)))
  }
  vapply(seq_len(draws), \(d) one(), numeric(3L + n * n))
}

eb_kurtosis <- function(v) mean((v - mean(v))^4) / stats::var(v)^2

test_that("posterior mean and SD of each row's entropy equal the Dirichlet closed forms (Monte Carlo summary)", {
  invisible(lapply(c(11L, 12L), \(seed) {
    fit <- entropy_bayes(eb_counts, prior = 0.5, draws = eb_draws, seed = seed)
    closed <- t(vapply(eb_states, \(s) eb_dirichlet_entropy(eb_counts[s, ] + 0.5), numeric(2)))
    # API gap: the per-state table is the field `states` (summary() returns the edge table).
    found <- subset(fit$states, select = c(state, row_entropy, row_entropy_sd))
    expected <- data.frame(state = eb_states, row_entropy = closed[, "mean"], row_entropy_sd = closed[, "sd"])
    # The SD's standard error uses a kurtosis of 4, above that of every row's entropy in a simulation of 1e5 draws
    # (2.98 to 3.48).
    expect_reference(subset(found, select = c(state, row_entropy)), subset(expected, select = c(state, row_entropy)),
                     reference = "published: Wolpert & Wolf (1995), posterior mean of the entropy under a Dirichlet",
                     covers = "entropy_bayes()", quantity = sprintf("Monte Carlo summary: row entropy mean, seed %d", seed),
                     tolerance = 4 * max(closed[, "sd"]) / sqrt(eb_draws), key = "state")
    expect_reference(subset(found, select = c(state, row_entropy_sd)), subset(expected, select = c(state, row_entropy_sd)),
                     reference = "published: Hutter (2002), posterior variance of the entropy under a Dirichlet",
                     covers = "entropy_bayes()", quantity = sprintf("Monte Carlo summary: row entropy SD, seed %d", seed),
                     tolerance = 4 * max(closed[, "sd"]) * sqrt((4 - 1) / (4 * eb_draws)), key = "state")
  }))
})

test_that("chain-level posterior means and per-edge contributions equal an independent simulation (Monte Carlo summary)", {
  reference <- eb_reference(eb_counts, prior = 0.5, draws = eb_draws, seed = 99L)
  means <- rowMeans(reference)
  sds <- apply(reference, 1L, stats::sd)
  # 4 standard errors of the difference of two independent means of eb_draws draws each.
  tolerance <- 4 * sds * sqrt(2 / eb_draws)
  fit <- entropy_bayes(eb_counts, prior = 0.5, draws = eb_draws, seed = 21L)
  chain <- c("entropy_rate", "stationary_entropy", "redundancy")
  # API gap: the chain-level table is the field `summary` (summary() prints it and returns the edge table).
  invisible(lapply(seq_along(chain), \(i) {
    expect_reference(subset(fit$summary, quantity == chain[[i]], select = c(quantity, mean)),
                     data.frame(quantity = chain[[i]], mean = means[[i]]),
                     reference = "independent: simulated Dirichlet posterior of the chain, base R", covers = "entropy_bayes()",
                     quantity = sprintf("Monte Carlo summary: %s mean", chain[[i]]), tolerance = tolerance[[i]],
                     key = "quantity")
    expect_reference(subset(fit$summary, quantity == chain[[i]], select = c(quantity, sd)),
                     data.frame(quantity = chain[[i]], sd = sds[[i]]),
                     reference = "independent: simulated Dirichlet posterior of the chain, base R", covers = "entropy_bayes()",
                     quantity = sprintf("Monte Carlo summary: %s SD", chain[[i]]),
                     tolerance = 4 * sds[[i]] * sqrt(2 * (eb_kurtosis(reference[i, ]) - 1) / (4 * eb_draws)),
                     key = "quantity")
  }))
  cells <- expand.grid(from = eb_states, to = eb_states, stringsAsFactors = FALSE)
  edge_rows <- 3L + seq_len(nrow(cells))
  expected <- data.frame(from = cells$from, to = cells$to, contribution = means[edge_rows])
  edges <- summary(fit)
  # One tolerance per edge: each edge is judged against its own Monte Carlo error.
  invisible(lapply(seq_len(nrow(cells)), \(i) {
    expect_reference(subset(edges, from == cells$from[[i]] & to == cells$to[[i]], select = c(from, to, contribution)),
                     expected[i, ], reference = "independent: simulated Dirichlet posterior of the chain, base R",
                     covers = c("entropy_bayes()", "summary.net_entropy_bayes()"),
                     quantity = sprintf("Monte Carlo summary: contribution mean %s -> %s", cells$from[[i]], cells$to[[i]]),
                     tolerance = tolerance[[edge_rows[[i]]]], key = c("from", "to"))
  }))
})

test_that("a grouped network gives each group's posterior (Monte Carlo summary)", {
  wide <- tn_simulate(n_actors = 60L, n_states = 3L, seq_length = 15L, seed = 431L)
  group <- rep(c("g1", "g2"), length.out = nrow(wide))
  network <- build_network(cbind(wide, grp = group), method = "frequency", group = "grp")
  found <- summary(entropy_bayes(network, prior = 0.5, draws = eb_draws, seed = 3L))
  invisible(lapply(c("g1", "g2"), \(g) {
    counts <- tn_pair_counts(tn_sequences(wide[group == g, ]))
    reference <- eb_reference(counts, prior = 0.5, draws = eb_draws, seed = 98L)
    expect_reference(subset(found, group == g & quantity == "entropy_rate", select = c(group, quantity, mean)),
                     data.frame(group = g, quantity = "entropy_rate", mean = mean(reference[1L, ])),
                     reference = "independent: simulated Dirichlet posterior of the chain, base R",
                     covers = c("entropy_bayes()", "summary.net_entropy_bayes_group()"),
                     quantity = sprintf("Monte Carlo summary: entropy rate mean, group %s", g),
                     tolerance = 4 * stats::sd(reference[1L, ]) * sqrt(2 / eb_draws), key = c("group", "quantity"))
  }))
})
