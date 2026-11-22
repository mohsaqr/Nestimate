# Layer 6, adversarial data: build_hypa() and markov_order_test() on the adversarial sequence sets, compared in
# full. HYPA (LaRock et al. 2020) is rebuilt here from its definition, independently of the package's fitted
# fields: the order-k De Bruijn graph of the sequences split at missing positions, its possible edges (a source with
# out-strength, a target with in-strength, overlapping in k - 1 states, whose last step is an observed first-order
# pair), the propensities Xi fitted by iterative proportional scaling to the observed strengths (Xi_ij in units of
# s_out s_in, summing to m^2), and the hypergeometric tails of every possible path, unobserved ones included. Every
# row of `scores` and every column is compared, with no filter on the observed count (test-build_hypa.R compares
# `subset(scores, observed > 0)` for the counts, against the package's own Xi for the tails).
# markov_order_test(): every deterministic column of its table at every order (log-likelihood of the multi-order
# chain, AIC and BIC with path or observed degrees of freedom, G2, df, asymptotic p-value), the selected orders and
# the counts (test-markov_order_test.R compares g2 and df only, of orders >= 1). The permutation p-values are drawn
# inside the function: only their range is asserted.
# Tolerances, fixed before comparing: counts exact; Xi and the expected counts 1e-8 (two scaling iterations, each
# run to its own stopping rule, 1e-13 here); hypergeometric tails 1e-8; log-likelihoods 1e-10 relative.

sets <- adv_sequence_sets()

# The observed runs of each row: split at a missing position, runs shorter than two states dropped.
adv_pieces <- function(wide) {
  pieces <- unlist(lapply(seq_len(nrow(wide)), \(i) {
    r <- unlist(wide[i, ], use.names = FALSE)
    observed <- !is.na(r)
    unname(split(r[observed], cumsum(!observed)[observed]))
  }), recursive = FALSE)
  Filter(\(p) length(p) >= 2L, pieces)
}

# Every window of `width` consecutive states of the pieces, as a character matrix (one row per window).
adv_windows <- function(pieces, width) {
  rows <- lapply(pieces, \(p) {
    n <- length(p)
    if (n < width) return(NULL)
    matrix(unlist(lapply(seq_len(n - width + 1L), \(s) p[s:(s + width - 1L)])), ncol = width, byrow = TRUE)
  })
  do.call(rbind, c(list(matrix(character(), 0L, width)), rows))
}

# Iterative proportional scaling of `start` to the row sums `r` and column sums `s`, until the largest margin error
# is below tol * sum(r) or after max_iter sweeps (a fixed point; no for loop: a bounded recursion).
adv_scale <- function(start, r, s, tol = 1e-13, max_iter = 100000L) {
  safe <- \(target, current) ifelse(current > 0, target / current, 0)
  sweep_once <- \(e) {
    e <- e * safe(r, rowSums(e))
    t(t(e) * safe(s, colSums(e)))
  }
  error <- \(e) max(abs(rowSums(e) - r), abs(colSums(e) - s))
  # sweeps in blocks of 500, so that the recursion stays shallow
  block <- \(e) Reduce(\(x, i) if (error(x) < tol * sum(r)) x else sweep_once(x), seq_len(500L), e)
  run <- \(e, done) if (error(e) < tol * sum(r) || done >= max_iter) e else run(block(e), done + 500L)
  run(start, 0L)
}

# HYPA scores of every possible path of order k, from the definition.
adv_hypa <- function(wide, k, alpha = 0.05, min_count = 1L) {
  pieces <- adv_pieces(wide)
  key <- \(m) apply(m, 1L, paste, collapse = "\r")
  grams <- adv_windows(pieces, k + 1L)
  nodes <- sort(unique(key(adv_windows(pieces, k))))
  source <- key(grams[, seq_len(k), drop = FALSE])
  target <- key(grams[, 1L + seq_len(k), drop = FALSE])
  adjacency <- matrix(as.numeric(table(factor(source, nodes), factor(target, nodes))), length(nodes),
                      dimnames = list(nodes, nodes))
  parts <- strsplit(nodes, "\r", fixed = TRUE)
  first_order <- unique(key(adv_windows(pieces, 2L)))
  s_out <- rowSums(adjacency)
  s_in <- colSums(adjacency)
  cells <- expand.grid(i = seq_along(nodes), j = seq_along(nodes))
  overlap <- vapply(seq_len(nrow(cells)), \(r) {
    a <- parts[[cells$i[r]]]
    b <- parts[[cells$j[r]]]
    identical(a[-1L], b[-k]) && paste(a[k], b[k], sep = "\r") %in% first_order
  }, logical(1))
  possible <- if (k == 1L) adjacency > 0 else
    matrix(overlap, length(nodes)) & outer(s_out > 0, s_in > 0)
  m <- sum(adjacency)
  start <- outer(s_out, s_in) * possible
  fitted <- adv_scale(start * m / sum(start), s_out, s_in)
  xi <- fitted * m
  index <- which(xi > 0, arr.ind = TRUE)
  pool <- round(sum(xi))
  draws <- min(m, pool)
  observed <- adjacency[index]
  marked <- pmin(round(xi[index]), pool)
  expected <- m * xi[index] / sum(xi)
  p_under <- stats::phyper(observed, marked, pool - marked, draws)
  p_over <- stats::phyper(observed, marked, pool - marked, draws, lower.tail = FALSE)
  context <- vapply(parts[index[, 1L]], paste, character(1), collapse = " -> ")
  following <- vapply(parts[index[, 2L]], \(p) p[k], character(1))
  data.frame(path = paste(context, following, sep = " -> "), from = context, to = following,
             observed = as.integer(observed), expected = expected,
             ratio = ifelse(expected > 0, observed / expected, Inf), p_value = p_under, p_under = p_under,
             p_over = p_over,
             anomaly = ifelse(expected >= min_count & p_under < alpha, "under",
                              ifelse(observed >= min_count & p_over <= alpha, "over", "normal")),
             p_adjusted_under = p_under, p_adjusted_over = p_over, order = k, stringsAsFactors = FALSE)
}

test_that("HYPA scores every possible path, observed or not, as the definition does, orders 1 to 3", {
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    longest <- max(lengths(adv_pieces(wide)))
    orders <- seq_len(min(3L, longest - 1L))
    expected <- do.call(rbind, lapply(orders, \(k) adv_hypa(wide, k)))
    found <- build_hypa(wide, order = orders, min_count = 1L)
    expect_reference(found$scores, expected, reference = "independent: HYPA from LaRock et al. (2020), base R",
                     covers = "build_hypa()", quantity = sprintf("every column of the scores, %s", name),
                     tolerance = 1e-8, key = "path")
    expect_reference(c(n_edges = found$n_edges, n_over = found$n_over, n_under = found$n_under,
                       n_anomalous = found$n_anomalous),
                     c(n_edges = nrow(expected), n_over = sum(expected$anomaly == "over"),
                       n_under = sum(expected$anomaly == "under"), n_anomalous = sum(expected$anomaly != "normal")),
                     reference = "independent: HYPA from LaRock et al. (2020), base R", covers = "build_hypa()",
                     quantity = sprintf("path and anomaly counts, %s", name), tolerance = 0)
  }))
})

test_that("unobserved possible paths are scored: an observed count of 0 with its lower tail", {
  # many_states: most possible paths of order 2 never occur
  expected <- adv_hypa(sets$many_states, 2L)
  found <- build_hypa(sets$many_states, order = 2L, min_count = 1L)
  unobserved <- \(scores) scores[scores$observed == 0L, c("path", "observed", "expected", "p_under", "p_over")]
  expect_reference(unobserved(found$scores), unobserved(expected),
                   reference = "independent: HYPA from LaRock et al. (2020), base R", covers = "build_hypa()",
                   quantity = "unobserved possible paths, many_states", tolerance = 1e-8, key = "path")
})

# ---- markov_order_test() ---------------------------------------------------------------------------------------

# G2, df and asymptotic p-value of order k against k - 1 (conditional independence of the first and the next state
# given the k - 1 between), summed over the contexts.
adv_order_g2 <- function(pieces, k) {
  windows <- adv_windows(pieces, k + 1L)
  context <- if (k == 1L) rep("", nrow(windows)) else apply(windows[, 1L + seq_len(k - 1L), drop = FALSE], 1L, paste,
                                                            collapse = "\r")
  parts <- lapply(split(seq_len(nrow(windows)), context), \(rows) {
    counts <- table(windows[rows, 1L], windows[rows, k + 1L])
    if (nrow(counts) < 2L || ncol(counts) < 2L) return(c(g2 = 0, df = 0))
    expected <- outer(rowSums(counts), colSums(counts)) / sum(counts)
    seen <- counts > 0
    c(g2 = 2 * sum(counts[seen] * log(counts[seen] / expected[seen])), df = (nrow(counts) - 1) * (ncol(counts) - 1))
  })
  total <- Reduce(`+`, parts)
  c(total, p = if (total[["df"]] > 0) stats::pchisq(total[["g2"]], total[["df"]], lower.tail = FALSE) else NA_real_)
}

# Log-likelihood of the multi-order chain of order k (Gote et al. 2020): the first state of a piece from the
# marginal distribution, the state at step t from the order min(t - 1, k) chain, all estimated by counting.
adv_multi_order_loglik <- function(pieces, k) {
  states <- unlist(pieces)
  marginal <- table(states) / length(states)
  steps <- do.call(rbind, lapply(pieces, \(p) {
    t <- seq_len(length(p) - 1L) + 1L
    used <- pmin(t - 1L, k)
    data.frame(order = used,
               context = vapply(seq_along(t), \(s) paste(p[t[s] - rev(seq_len(used[s]))], collapse = "\r"), character(1)),
               following = p[t], stringsAsFactors = FALSE)
  }))
  probability <- vapply(seq_len(nrow(steps)), \(r) {
    if (steps$order[r] == 0L) return(as.numeric(marginal[[steps$following[r]]]))
    windows <- adv_windows(pieces, steps$order[r] + 1L)
    contexts <- apply(windows[, seq_len(steps$order[r]), drop = FALSE], 1L, paste, collapse = "\r")
    in_context <- contexts == steps$context[r]
    sum(in_context & windows[, steps$order[r] + 1L] == steps$following[r]) / sum(in_context)
  }, numeric(1))
  firsts <- vapply(pieces, `[`, character(1), 1L)
  sum(log(marginal[firsts])) + sum(log(probability))
}

# Free parameters of each layer: order 0 has (states - 1); the order-j layer has, with dof = "paths", one per
# walk of j - 1 steps in the observed first-order graph ending in state i, times (out-degree of i - 1); with
# dof = "observed", the sum over the order-j contexts of (distinct next states - 1).
adv_layer_dofs <- function(pieces, max_order, dof) {
  states <- sort(unique(unlist(pieces)))
  pairs <- unique(adv_windows(pieces, 2L))
  out_degree <- table(factor(pairs[, 1L], states))
  extra <- pmax(as.numeric(out_degree) - 1, 0)
  names(extra) <- states
  edge <- paste(pairs[, 1L], pairs[, 2L], sep = "\r")
  walks_ending <- \(steps) {
    if (steps == 0L) return(stats::setNames(rep(1, length(states)), states))
    walks <- expand.grid(rep(list(states), steps + 1L), stringsAsFactors = FALSE)
    valid <- Reduce(`&`, lapply(seq_len(steps), \(s) paste(walks[[s]], walks[[s + 1L]], sep = "\r") %in% edge))
    as.numeric(table(factor(walks[[steps + 1L]][valid], states)))
  }
  layers <- vapply(seq_len(max_order), \(j) {
    if (dof == "paths") return(sum(walks_ending(j - 1L) * extra))
    windows <- adv_windows(pieces, j + 1L)
    contexts <- apply(windows[, seq_len(j), drop = FALSE], 1L, paste, collapse = "\r")
    distinct <- tapply(windows[, j + 1L], contexts, \(x) length(unique(x)))
    sum(pmax(distinct - 1, 0))
  }, numeric(1))
  c(length(states) - 1, layers)
}

adv_order_table <- function(wide, max_order, dof) {
  pieces <- adv_pieces(wide)
  loglik <- vapply(0:max_order, \(k) adv_multi_order_loglik(pieces, k), numeric(1))
  cumulative <- cumsum(adv_layer_dofs(pieces, max_order, dof))
  n_obs <- length(unlist(pieces))
  tests <- lapply(seq_len(max_order), \(k) adv_order_g2(pieces, k))
  data.frame(order = 0:max_order, loglik = loglik, AIC = 2 * cumulative - 2 * loglik,
             BIC = log(n_obs) * cumulative - 2 * loglik,
             df = c(NA, vapply(tests, \(x) x[["df"]], numeric(1))),
             g2 = c(NA, vapply(tests, \(x) x[["g2"]], numeric(1))),
             p_asymptotic = c(NA, vapply(tests, \(x) x[["p"]], numeric(1))))
}

test_that("every deterministic column of markov_order_test() equals the definitions, both dof rules", {
  invisible(lapply(names(sets), \(name) {
    wide <- sets[[name]]
    longest <- max(lengths(adv_pieces(wide)))
    max_order <- min(3L, longest - 1L)
    invisible(lapply(c("paths", "observed"), \(dof) {
      found <- markov_order_test(wide, max_order = max_order, n_perm = 19L, seed = 1L, dof = dof)
      expected <- adv_order_table(wide, max_order, dof)
      table <- summary(found)
      expect_reference(table, expected, reference = "independent: multi-order chain and stratified G2, base R",
                       covers = c("markov_order_test()", "summary.net_markov_order()"),
                       quantity = sprintf("order table, %s, dof %s", name, dof), tolerance = 1e-10,
                       measure = "relative", key = "order")
      expect_reference(c(bic_order = found$bic_order, aic_order = found$aic_order, n_sequences = found$n_sequences,
                         n_observations = found$n_observations),
                       c(bic_order = expected$order[which.min(expected$BIC)],
                         aic_order = expected$order[which.min(expected$AIC)],
                         n_sequences = length(adv_pieces(wide)), n_observations = length(unlist(adv_pieces(wide)))),
                       reference = "independent: multi-order chain and stratified G2, base R",
                       covers = "markov_order_test()", quantity = sprintf("selected orders and counts, %s, dof %s", name, dof),
                       tolerance = 0)
      p <- table$p_permutation[-1L]
      tested <- !is.na(expected$p_asymptotic[-1L])
      expect_true(all(p[tested] >= 1 / 20 & p[tested] <= 1) && all(is.na(p[!tested])),
                  info = sprintf("permutation p-values in [1/(n_perm + 1), 1], NA where df = 0: %s, %s", name, dof))
      expect_identical(table$significant[-1L], ifelse(is.na(p), FALSE, p < 0.05),
                       info = sprintf("significant is p_permutation < alpha: %s, %s", name, dof))
    }))
  }))
})
