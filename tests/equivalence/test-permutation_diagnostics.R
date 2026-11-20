# permutation_diagnostics(): a sequence-level and an actor-level permutation test of two networks side by side, with
# the intraclass correlation that explains their difference. The ICC is deterministic and is written here from its
# documented definition (?permutation, "ICC and design effect"): the one-way ANOVA ICC of each sequence's transition
# shares within each network, (MSB - MSW) / (MSB + (n0 - 1) MSW) with the unbalanced cluster size n0 = (N - sum
# n_k^2 / N) / (k - 1), the mean squares taken from stats::anova() of a one-way model; the two networks' ICCs of an
# edge pooled by their numbers of actors; edges averaged weighted by their counts; jackknife bias-corrected over
# actors with the interval estimate +/- t(k - 1) x jackknife SE (Efron & Tibshirani 1993, ch. 11). The observed
# differences are those of the two networks' relative weights. The p-values come from permutations drawn inside the
# function with only a seed: they are checked by a property of any permutation p-value (Phipson & Smyth 2010), and
# the smallest p-value the actor test can reach by its documented value.

pd_iter <- 200L
pd_states <- c("A", "B", "C")

# 16 students with four sequences each, eight per group; each student leans towards one state.
pd_data <- function(seed) {
  set.seed(seed)
  student <- rep(sprintf("s%02d", seq_len(16L)), each = 4L)
  lean <- as.integer(substr(student, 2L, 3L)) %% 3L
  wide <- as.data.frame(do.call(rbind, lapply(lean, \(l) {
    sample(pd_states, 12L, replace = TRUE, prob = 1 + 1.5 * (l == 0:2))
  })))
  names(wide) <- paste0("T", seq_len(12L))
  list(wide = wide, student = student, grp = ifelse(as.integer(substr(student, 2L, 3L)) <= 8L, "g1", "g2"))
}

# One-way ANOVA ICC of y across clusters, unbalanced n0; NA when y does not vary.
pd_anova_icc <- function(y, cluster) {
  if (stats::var(y) == 0) return(NA_real_)
  mean_squares <- stats::anova(stats::lm(y ~ factor(cluster)))[["Mean Sq"]]
  sizes <- as.vector(table(cluster))
  n0 <- (sum(sizes) - sum(sizes^2) / sum(sizes)) / (length(sizes) - 1L)
  (mean_squares[[1L]] - mean_squares[[2L]]) / (mean_squares[[1L]] + (n0 - 1) * mean_squares[[2L]])
}

pd_icc <- function(data) {
  sequences <- tn_sequences(data$wide)
  cells <- expand.grid(to = pd_states, from = pd_states, stringsAsFactors = FALSE)
  counts <- t(vapply(sequences, \(s) as.vector(t(tn_pair_counts(list(s), pd_states))), numeric(nrow(cells))))
  shares <- counts / rowSums(counts)
  edge_icc <- function(rows) vapply(seq_len(nrow(cells)), \(e) {
    per <- lapply(split(rows, data$grp[rows]), \(r) {
      c(icc = pd_anova_icc(shares[r, e], data$student[r]), n = length(unique(data$student[r])))
    })
    icc <- vapply(per, \(p) p[["icc"]], numeric(1))
    n <- vapply(per, \(p) p[["n"]], numeric(1))
    if (all(is.na(icc))) NA_real_ else sum((icc * n)[!is.na(icc)]) / sum(n[!is.na(icc)])
  }, numeric(1))
  pooled <- function(rows) {
    icc <- edge_icc(rows)
    weight <- colSums(counts[rows, , drop = FALSE])
    ok <- !is.na(icc) & weight > 0
    sum(icc[ok] * weight[ok]) / sum(weight[ok])
  }
  rows <- seq_along(sequences)
  blocks <- split(rows, data$student)
  k <- length(blocks)
  jackknife <- vapply(seq_len(k), \(j) pooled(unlist(blocks[-j], use.names = FALSE)), numeric(1))
  corrected <- k * pooled(rows) - (k - 1) * mean(jackknife)
  se <- sqrt((k - 1) / k * sum((jackknife - mean(jackknife))^2))
  list(overall = data.frame(pair = "g1 vs g2", icc = corrected,
                            icc_ci_lower = corrected - stats::qt(0.975, k - 1) * se,
                            icc_ci_upper = corrected + stats::qt(0.975, k - 1) * se),
       edges = data.frame(pair = "g1 vs g2", from = cells$from, to = cells$to, icc = edge_icc(rows)))
}

test_that("the ICC, its interval and the per-edge ICCs equal the ANOVA ICC from stats::anova() with a jackknife", {
  invisible(lapply(c(461L, 462L), \(seed) {
    data <- pd_data(seed)
    network <- build_network(cbind(data$wide, student = data$student, grp = data$grp), method = "relative",
                             group = "grp")
    expected <- pd_icc(data)
    expect_reference(permutation_diagnostics(network, actor = "student", iter = pd_iter, seed = 1), expected$overall,
                     reference = "stats", covers = "permutation_diagnostics()",
                     quantity = sprintf("jackknife ICC and interval, seed %d", seed), tolerance = 1e-10, key = "pair",
                     note = "mean squares from stats::anova(); n0, pooling and jackknife written from the definition")
    edges <- permutation_diagnostics(network, actor = "student", iter = pd_iter, seed = 1, level = "edges")
    expect_reference(subset(edges, select = c(pair, from, to, icc)), expected$edges, reference = "stats",
                     covers = "permutation_diagnostics()", quantity = sprintf("per-edge ANOVA ICC, seed %d", seed),
                     tolerance = 1e-10, key = c("pair", "from", "to"),
                     note = "mean squares from stats::anova(); n0 and pooling written from the definition")
  }))
})

test_that("observed differences are those of the relative weights; p-values lie where any permutation p-value lies", {
  data <- pd_data(461L)
  network <- build_network(cbind(data$wide, student = data$student, grp = data$grp), method = "relative", group = "grp")
  weights <- lapply(c(g1 = "g1", g2 = "g2"), \(g) tn_row_normalise(tn_pair_counts(tn_sequences(data$wide[data$grp == g, ]),
                                                                                    pd_states)))
  difference <- weights$g1 - weights$g2
  edges <- permutation_diagnostics(network, actor = "student", iter = pd_iter, seed = 1, level = "edges")
  expect_reference(subset(edges, select = c(from, to, diff)),
                   data.frame(from = rep(pd_states, each = 3L), to = rep(pd_states, 3L), diff = as.vector(t(difference))),
                   reference = "independent: difference of the row-normalised pair tables, base R",
                   covers = "permutation_diagnostics()", quantity = "observed edge difference", tolerance = 1e-12,
                   key = c("from", "to"))
  # A permutation p-value (count + 1) / (B + 1) lies in [1 / (B + 1), 1] for any draw (Phipson & Smyth 2010).
  p_values <- subset(edges, select = c(from, to, p_sequence, p_actor))
  bounded <- within(p_values, {
    p_sequence <- pmin(pmax(p_sequence, 1 / (pd_iter + 1)), 1)
    p_actor <- pmin(pmax(p_actor, 1 / (pd_iter + 1)), 1)
  })
  expect_reference(p_values, bounded,
                   reference = "published: Phipson & Smyth (2010), a permutation p-value lies in [1/(B+1), 1]",
                   covers = "permutation_diagnostics()", quantity = "Monte Carlo summary: p-values within their bounds",
                   tolerance = 0, key = c("from", "to"))
  # 16 actors split 8 / 8 have choose(16, 8) = 12870 arrangements, so the smallest p-value is 1 / (iter + 1).
  overall <- permutation_diagnostics(network, actor = "student", iter = pd_iter, seed = 1)
  expect_reference(subset(overall, select = c(pair, min_p_actor)),
                   data.frame(pair = "g1 vs g2", min_p_actor = max(1 / choose(16, 8), 1 / (pd_iter + 1))),
                   reference = "independent: documented smallest attainable p-value", covers = "permutation_diagnostics()",
                   quantity = "smallest attainable actor-level p-value", tolerance = 1e-15, key = "pair")
})
