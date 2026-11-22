# entropy_network(), path_dependence(), entropy_trajectory() and the group summaries of transition_entropy() and
# chain_structure(): information-theoretic summaries of a Markov chain (Cover and Thomas 2006, chapters 2 and 4;
# Krejtz et al. 2015, 2025). References: the stationary distribution from markovchain::steadyStates(), entropies
# from entropy::entropy.plugin(), and the per-edge weights, conditional entropies, Kullback-Leibler divergences and
# windowed empirical entropies written from their definitions in base R. New: these functions had no comparison
# with a reference. Tables a function only returns as a field or a list element are read with `$` (API gap).

matrix_p <- ho_stochastic(5L, 1201L)

# Consecutive-pair counts of a list of sequences, as a states x states matrix.
ho_pair_counts <- function(sequences, states) {
  pairs <- do.call(rbind, lapply(sequences, \(s) data.frame(from = utils::head(s, -1L), to = utils::tail(s, -1L))))
  unclass(table(factor(pairs$from, states), factor(pairs$to, states)))
}

test_that("contribution, surprisal and production weights equal their definitions", {
  skip_without_reference("markovchain")
  stationary <- ho_stationary(matrix_p)[rownames(matrix_p)]
  flow <- stationary * matrix_p
  expected <- list(contribution = flow * log2(1 / matrix_p), surprisal = log2(1 / matrix_p),
                   production = flow * log2(flow / t(flow)))
  diag(expected$production) <- 0
  invisible(lapply(names(expected), \(weight) {
    expect_reference(extract_transition_matrix(entropy_network(matrix_p, weight = weight)), expected[[weight]],
                     reference = "markovchain", covers = "entropy_network()", quantity = sprintf("%s weights", weight),
                     tolerance = 1e-10, note = "stationary distribution from markovchain::steadyStates(); weights from their definitions")
  }))
})

test_that("context entropies and Kullback-Leibler divergences of path_dependence() equal their definitions", {
  skip_without_reference("entropy")
  wide <- ho_sequences(40L, 25L, ho_stochastic(3L, 1202L), 1203L)
  sequences <- lapply(seq_len(nrow(wide)), \(i) unlist(wide[i, ], use.names = FALSE))
  states <- sort(unique(unlist(sequences)))
  first <- ho_pair_counts(sequences, states)
  first <- first / rowSums(first)
  triples <- do.call(rbind, lapply(sequences, \(s) {
    positions <- seq_len(length(s) - 2L)
    data.frame(context = paste(s[positions], s[positions + 1L], sep = " -> "), last = s[positions + 1L],
               next_state = s[positions + 2L])
  }))
  rows <- lapply(split(triples, triples$context), \(rows) {
    p_k <- as.vector(table(factor(rows$next_state, states))) / nrow(rows)
    p_1 <- first[rows$last[[1L]], ]
    used <- p_k > 0
    data.frame(context = rows$context[[1L]], n = nrow(rows), H_order1 = entropy::entropy.plugin(p_1, unit = "log2"),
               H_orderk = entropy::entropy.plugin(p_k, unit = "log2"),
               KL = sum(p_k[used] * log2(p_k[used] / p_1[used])))
  })
  expected <- subset(do.call(rbind, rows), n >= 5L)
  expected$H_drop <- expected$H_order1 - expected$H_orderk
  found <- summary(path_dependence(wide, order = 2L, min_count = 5L))
  expect_reference(subset(found$contexts, select = c(context, n, H_order1, H_orderk, H_drop, KL)), expected,
                   reference = "entropy", covers = c("path_dependence()", "summary.net_path_dependence()"),
                   quantity = "order-2 contexts", tolerance = 1e-10, key = "context",
                   note = "entropies by entropy::entropy.plugin(); the divergence written from its definition")
  expect_reference(c(KL_weighted = found$chain$KL_weighted, H_drop_weighted = found$chain$H_drop_weighted),
                   c(KL_weighted = stats::weighted.mean(expected$KL, expected$n),
                     H_drop_weighted = stats::weighted.mean(expected$H_drop, expected$n)),
                   reference = "entropy", covers = c("path_dependence()", "summary.net_path_dependence()"),
                   quantity = "count-weighted chain summaries", tolerance = 1e-10)
})

test_that("windowed entropies of entropy_trajectory() equal the empirical conditional entropy of each window", {
  sequence <- unlist(ho_sequences(1L, 601L, ho_stochastic(4L, 1204L), 1205L), use.names = FALSE)
  log <- data.frame(actor = "a", action = sequence)
  found <- entropy_trajectory(log, action = "action", actor = "actor", window = 100L, step = 50L)
  from <- utils::head(sequence, -1L)
  to <- utils::tail(sequence, -1L)
  starts <- seq(1L, length(from) - 99L, by = 50L)
  entropy <- vapply(starts, \(start) {
    window <- start + 0:99
    counts <- table(from[window], to[window])
    shares <- counts / sum(counts)
    conditional <- counts / rowSums(counts)
    -sum(shares[counts > 0] * log2(conditional[counts > 0]))
  }, numeric(1))
  expect_reference(subset(found$trajectory, select = c(window, entropy)),
                   data.frame(window = seq_along(starts), entropy = entropy),
                   reference = "independent: empirical conditional entropy of each window of transitions, base R",
                   covers = "entropy_trajectory()", quantity = "window entropies", tolerance = 1e-10, key = "window")
  expect_reference(summary(found),
                   data.frame(group = "all", windows = length(entropy), mean = mean(entropy), sd = stats::sd(entropy), min = min(entropy),
                              max = max(entropy), first = entropy[[1L]], last = utils::tail(entropy, 1L),
                              change = utils::tail(entropy, 1L) - entropy[[1L]]),
                   reference = "independent: empirical conditional entropy of each window of transitions, base R",
                   covers = c("entropy_trajectory()", "summary.net_entropy_trajectory()"), quantity = "trajectory summary",
                   tolerance = 1e-10, key = "group")
})

test_that("group summaries of transition_entropy() and chain_structure() equal markovchain on each group's matrix", {
  skip_without_reference("markovchain")
  skip_without_reference("entropy")
  wide <- rbind(ho_sequences(30L, 20L, ho_stochastic(4L, 1206L), 1207L), ho_sequences(30L, 20L, ho_stochastic(4L, 1208L), 1209L))
  wide$grp <- rep(c("a", "b"), each = 30L)
  networks <- build_network(wide, method = "relative", group = "grp")
  per_group <- lapply(c(a = "a", b = "b"), \(label) {
    rows <- subset(wide, grp == label, select = -grp)
    counts <- ho_pair_counts(lapply(seq_len(nrow(rows)), \(i) unlist(rows[i, ], use.names = FALSE)), paste0("s", 1:4))
    counts / rowSums(counts)
  })
  chain <- do.call(rbind, lapply(names(per_group), \(label) {
    p <- per_group[[label]]
    stationary <- ho_stationary(p)[rownames(p)]
    rows <- vapply(rownames(p), \(state) entropy::entropy.plugin(p[state, ], unit = "log2"), numeric(1))
    rate <- sum(stationary * rows)
    stationary_entropy <- entropy::entropy.plugin(stationary, unit = "log2")
    data.frame(group = label, entropy_rate = rate, stationary_entropy = stationary_entropy,
               redundancy = stationary_entropy - rate)
  }))
  expect_reference(summary(transition_entropy(networks)), chain, reference = "markovchain",
                   covers = c("transition_entropy()", "summary.net_transition_entropy_group()"),
                   quantity = "chain entropies per group", tolerance = 1e-8, key = "group",
                   note = "stationary distribution from markovchain::steadyStates(); entropies from entropy::entropy.plugin()")
  states <- do.call(rbind, lapply(names(per_group), \(label) {
    p <- per_group[[label]]
    data.frame(group = label, state = rownames(p), persistence = diag(p),
               stationary_probability = ho_stationary(p)[rownames(p)])
  }))
  # summary() rounds persistence and the stationary probability to four decimals.
  expect_reference(subset(summary(chain_structure(networks)), select = c(group, state, persistence, stationary_probability)),
                   states, reference = "markovchain", covers = c("chain_structure()", "summary.chain_structure_group()"),
                   quantity = "persistence and stationary probability per group", tolerance = 5e-5, key = c("group", "state"))
})

test_that("windowed entropies equal H(X, Y) - H(X) by entropy::entropy.plugin()", {
  skip_without_reference("entropy")
  invisible(lapply(c(1205L, 1206L), \(seed) {
    sequence <- unlist(ho_sequences(1L, 601L, ho_stochastic(4L, 1204L), seed), use.names = FALSE)
    found <- entropy_trajectory(data.frame(actor = "a", action = sequence), action = "action", actor = "actor",
                                window = 100L, step = 50L)
    from <- utils::head(sequence, -1L)
    to <- utils::tail(sequence, -1L)
    starts <- seq(1L, length(from) - 99L, by = 50L)
    entropy <- vapply(starts, \(start) {
      window <- start + 0:99
      entropy::entropy.plugin(table(from[window], to[window]), unit = "log2") -
        entropy::entropy.plugin(table(from[window]), unit = "log2")
    }, numeric(1))
    expect_reference(subset(found$trajectory, select = c(window, entropy)),
                     data.frame(window = seq_along(starts), entropy = entropy), reference = "entropy",
                     covers = "entropy_trajectory()", quantity = sprintf("window entropies, seed %d, entropy", seed),
                     tolerance = 1e-10, key = "window")
    expect_reference(summary(found),
                     data.frame(group = "all", windows = length(entropy), mean = mean(entropy), sd = stats::sd(entropy),
                                min = min(entropy), max = max(entropy), first = entropy[[1L]],
                                last = utils::tail(entropy, 1L), change = utils::tail(entropy, 1L) - entropy[[1L]]),
                     reference = "entropy", covers = "summary.net_entropy_trajectory()",
                     quantity = sprintf("trajectory summary, seed %d, entropy", seed), tolerance = 1e-10, key = "group")
  }))
})
