# casedrop_reliability() and its summaries: edge-weight case-dropping stability (the CS-coefficient of Epskamp,
# Borsboom & Fried 2018, applied to the edge weights). No established package computes it for edge weights
# (bootnet's case-drop is for centralities), and the subsets are drawn inside the function with only a seed, so the
# comparison is statistical (Monte Carlo summaries), declared as such: the documented procedure (?casedrop_reliability)
# is written here in base R with its own draws: for each drop proportion p, subsets of n - floor(n p) sequences
# drawn without replacement, the relative transition network of the subset, and its existing off-diagonal edges
# set against the original's (Spearman correlation, mean, median and largest absolute change). The mean and SD of
# every metric at every proportion are compared with a tolerance fixed before the run: 4 standard errors of the
# difference of the two Monte Carlo estimates (a mean: sd * sqrt(2 / B); an SD: sd * sqrt(2 (kurtosis - 1) / (4 B)),
# the kurtosis read from the reference's draws). A seed sweep of 20 seeds before fixing these passed 20 of 20.

cd_iter <- 1000L
cd_drop <- c(0.2, 0.5, 0.7)
cd_metrics <- c("mean_abs_dev", "median_abs_dev", "correlation", "max_abs_dev")

# Per-iteration metrics of the documented procedure, one matrix per metric (iterations x proportions).
cd_reference <- function(sequences, states, seed) {
  set.seed(seed)
  per_sequence <- lapply(sequences, \(s) tn_pair_counts(list(s), states))
  original <- tn_row_normalise(Reduce(`+`, per_sequence))
  cells <- original != 0 & row(original) != col(original)
  n <- length(sequences)
  one <- function(p) {
    kept <- sample.int(n, n - floor(n * p))
    changes <- tn_row_normalise(Reduce(`+`, per_sequence[kept]))[cells] - original[cells]
    c(mean_abs_dev = mean(abs(changes)), median_abs_dev = stats::median(abs(changes)),
      correlation = stats::cor(original[cells], original[cells] + changes, method = "spearman"),
      max_abs_dev = max(abs(changes)))
  }
  draws <- lapply(cd_drop, \(p) vapply(seq_len(cd_iter), \(i) one(p), numeric(4)))
  list(draws = draws, n_edges = sum(cells))
}

cd_kurtosis <- function(v) mean((v - mean(v))^4) / stats::var(v)^2

# Compare one result's summary() with the reference draws: mean and SD per metric and proportion.
cd_compare <- function(found, reference, covers, label) {
  invisible(lapply(seq_along(cd_drop), \(j) lapply(cd_metrics, \(m) {
    v <- reference$draws[[j]][m, ]
    row <- subset(found, metric == m & abs(drop_prop - cd_drop[[j]]) < 1e-9, select = c(metric, drop_prop, mean, sd))
    expect_reference(subset(row, select = c(metric, mean)), data.frame(metric = m, mean = mean(v)),
                     reference = "independent: case-dropping written from the definition, base R", covers = covers,
                     quantity = sprintf("Monte Carlo summary: %s, %s mean, drop %.1f", label, m, cd_drop[[j]]),
                     tolerance = 4 * stats::sd(v) * sqrt(2 / cd_iter), key = "metric")
    expect_reference(subset(row, select = c(metric, sd)), data.frame(metric = m, sd = stats::sd(v)),
                     reference = "independent: case-dropping written from the definition, base R", covers = covers,
                     quantity = sprintf("Monte Carlo summary: %s, %s SD, drop %.1f", label, m, cd_drop[[j]]),
                     tolerance = 4 * stats::sd(v) * sqrt(2 * (cd_kurtosis(v) - 1) / (4 * cd_iter)), key = "metric")
  })))
}

test_that("case-dropping metrics of a relative network equal the documented procedure (Monte Carlo summary)", {
  wide <- tn_simulate(n_actors = 30L, n_states = 4L, seq_length = 20L, seed = 441L)
  sequences <- tn_sequences(wide)
  states <- sort(unique(unlist(sequences)))
  reference <- cd_reference(sequences, states, seed = 9L)
  found <- casedrop_reliability(build_network(wide, method = "relative"), iter = cd_iter, drop_prop = cd_drop, seed = 4)
  cd_compare(summary(found), reference, c("casedrop_reliability()", "summary.net_casedrop_reliability()"),
             "30 sequences")
})

test_that("a grouped network gives each group's case-dropping metrics (Monte Carlo summary)", {
  wide <- tn_simulate(n_actors = 50L, n_states = 3L, seq_length = 20L, seed = 442L)
  group <- rep(c("g1", "g2"), length.out = nrow(wide))
  found <- casedrop_reliability(build_network(cbind(wide, grp = group), method = "relative", group = "grp"),
                                iter = cd_iter, drop_prop = cd_drop, seed = 4)
  references <- lapply(c(g1 = "g1", g2 = "g2"), \(g) {
    sequences <- tn_sequences(wide[group == g, ])
    cd_reference(sequences, sort(unique(unlist(sequences))), seed = 10L)
  })
  invisible(lapply(c("g1", "g2"), \(g) {
    cd_compare(summary(found[[g]]), references[[g]], c("casedrop_reliability()", "summary.net_casedrop_reliability()"),
               sprintf("group %s", g))
  }))
  # The group summary formats each mean and SD as text ("0.912 +/- 0.040", API gap); its number of edges compared is
  # a number: the existing off-diagonal edges of each group's network.
  expect_reference(data.frame(group = c("g1", "g2"), n_edges = summary(found, drop_prop = 0.5)$n_edges),
                   data.frame(group = c("g1", "g2"), n_edges = vapply(references, \(r) r$n_edges, numeric(1))),
                   reference = "independent: case-dropping written from the definition, base R",
                   covers = c("casedrop_reliability()", "summary.net_casedrop_reliability_group()"),
                   quantity = "edges compared per group", tolerance = 0, key = "group")
})
