# summary() of the bootstrap and of the permutation test of a mixed window network (wtna(method = "both")): the
# transition and co-occurrence networks resampled together. The resamples are drawn inside bootstrap_network() and
# permutation() with only a seed, so the resampled quantities are compared statistically (Monte Carlo summaries),
# declared as such. The reference counts come from tna::import_onehot() (tumbling windows of two rows) per actor,
# with tna::ftna() for transitions and tna::build_model(type = "co-occurrence") for co-occurrences (the owner's
# co-authored package, a separate implementation; open-source, to be confirmed by the owner); their sum over the
# actors equals wtna()'s networks (compared in test-wtna.R).
# - bootstrap: the actors are the cases (?bootstrap_network); boot::boot() resamples the actors with its own draws
#   and sums their counts. The bootstrap mean and SD of every edge are compared with a tolerance fixed before the
#   run: 4 standard errors of the difference of the two Monte Carlo estimates (a mean: sd * sqrt(2 / B); an SD:
#   sd * sqrt(2 (kurtosis - 1) / (4 B)), the kurtosis read from boot's replicates). A seed sweep of 20 seeds before
#   fixing these passed 20 of 20.
# - permutation: the observed weights of both networks and their difference are no draws and are compared exactly;
#   the p-values by the bounds that hold for any permutation p-value (Phipson & Smyth 2010).

wr_iter <- 2000L
wr_codes <- c("A", "B", "C", "D")

wr_data <- function(seed, n_actors = 20L, rate = 0.4) {
  set.seed(seed)
  data <- as.data.frame(matrix(stats::rbinom(4L * 10L * n_actors, 1L, rate), ncol = 4L,
                               dimnames = list(NULL, wr_codes)))
  data$actor <- rep(seq_len(n_actors), each = 10L)
  data
}

# Each actor's transition and co-occurrence counts as tna builds them.
wr_per_actor <- function(data) {
  lapply(split(data, data$actor), \(rows) {
    imported <- suppressWarnings(tna::import_onehot(rows, cols = wr_codes, actor = "actor", window_size = 2,
                                                    window_type = "tumbling"))
    list(transition = wr_full(unclass(tna::ftna(imported)$weights)),
         cooccurrence = wr_full(unclass(tna::build_model(imported, type = "co-occurrence", scaling = character(0))$weights)))
  })
}

# An actor's matrix on all four codes: tna leaves out a code the actor never shows, which then has no count.
wr_full <- function(m) {
  full <- matrix(0, length(wr_codes), length(wr_codes), dimnames = list(wr_codes, wr_codes))
  full[rownames(m), colnames(m)] <- m
  full
}

wr_kurtosis <- function(v) mean((v - mean(v))^4) / stats::var(v)^2

test_that("bootstrap means and SDs of both networks equal boot::boot() over the actors (Monte Carlo summary)", {
  skip_without_reference("tna")
  skip_without_reference("boot")
  data <- wr_data(481L)
  per_actor <- wr_per_actor(data)
  found <- summary(bootstrap_network(wtna(data, method = "both", codes = wr_codes, actor = "actor", window_size = 2,
                                          mode = "non-overlapping"), iter = wr_iter, seed = 2))
  invisible(lapply(c("transition", "cooccurrence"), \(part) {
    set.seed(8)
    replicates <- boot::boot(seq_along(per_actor), \(actors, rows) {
      as.vector(Reduce(`+`, lapply(per_actor[actors[rows]], \(a) a[[part]])))
    }, R = wr_iter)$t
    cells <- expand.grid(from = wr_codes, to = wr_codes, stringsAsFactors = FALSE)
    # The co-occurrence summary lists each unordered pair once (from <= to).
    keep <- if (identical(part, "cooccurrence")) cells$from <= cells$to else rep(TRUE, nrow(cells))
    v <- replicates[, keep, drop = FALSE]
    expected <- data.frame(cells[keep, ], mean = colMeans(v), sd = apply(v, 2L, stats::sd))
    sds <- expected$sd
    kurtosis <- apply(v, 2L, wr_kurtosis)
    # API gap: summary() returns a list of two tables, read here by name.
    table <- found[[part]]
    invisible(lapply(seq_len(nrow(expected)), \(i) {
      edge <- subset(table, from == expected$from[[i]] & to == expected$to[[i]])
      expect_reference(subset(edge, select = c(from, to, mean)), subset(expected[i, ], select = c(from, to, mean)),
                       reference = "boot", covers = "summary.wtna_boot_mixed()",
                       quantity = sprintf("Monte Carlo summary: %s bootstrap mean %s -> %s", part, expected$from[[i]],
                                          expected$to[[i]]),
                       tolerance = 4 * sds[[i]] * sqrt(2 / wr_iter), key = c("from", "to"),
                       note = "boot::boot() over the actors, summing tna's per-actor counts")
      expect_reference(subset(edge, select = c(from, to, sd)), subset(expected[i, ], select = c(from, to, sd)),
                       reference = "boot", covers = "summary.wtna_boot_mixed()",
                       quantity = sprintf("Monte Carlo summary: %s bootstrap SD %s -> %s", part, expected$from[[i]],
                                          expected$to[[i]]),
                       tolerance = 4 * sds[[i]] * sqrt(2 * (kurtosis[[i]] - 1) / (4 * wr_iter)), key = c("from", "to"),
                       note = "boot::boot() over the actors, summing tna's per-actor counts")
    }))
  }))
})

test_that("permutation of two mixed networks: observed weights equal tna's, p-values lie within their bounds", {
  skip_without_reference("tna")
  first <- wr_data(482L)
  second <- wr_data(483L, rate = 0.3)
  build <- \(data) wtna(data, method = "both", codes = wr_codes, actor = "actor", window_size = 2,
                        mode = "non-overlapping")
  found <- summary(permutation(build(first), build(second), iter = 500L, seed = 2))
  totals <- lapply(list(x = first, y = second), \(data) {
    per_actor <- wr_per_actor(data)
    lapply(c(transition = "transition", cooccurrence = "cooccurrence"),
           \(part) Reduce(`+`, lapply(per_actor, \(a) a[[part]])))
  })
  invisible(lapply(c("transition", "cooccurrence"), \(part) {
    cells <- expand.grid(from = wr_codes, to = wr_codes, stringsAsFactors = FALSE)
    keep <- if (identical(part, "cooccurrence")) cells$from <= cells$to else rep(TRUE, nrow(cells))
    weight_x <- as.vector(totals$x[[part]])[keep]
    weight_y <- as.vector(totals$y[[part]])[keep]
    expect_reference(subset(found[[part]], select = c(from, to, weight_x, weight_y, diff)),
                     data.frame(cells[keep, ], weight_x = weight_x, weight_y = weight_y, diff = weight_x - weight_y),
                     reference = "tna", covers = "summary.wtna_perm_mixed()",
                     quantity = sprintf("%s: observed weights and difference", part), tolerance = 1e-10,
                     key = c("from", "to"), note = "tna::import_onehot() per actor, summed")
    # A permutation p-value (count + 1) / (B + 1) lies in [1 / (B + 1), 1] for any draw (Phipson & Smyth 2010).
    p_values <- subset(found[[part]], select = c(from, to, p_value))
    expect_reference(p_values, within(p_values, p_value <- pmin(pmax(p_value, 1 / 501), 1)),
                     reference = "published: Phipson & Smyth (2010), a permutation p-value lies in [1/(B+1), 1]",
                     covers = "summary.wtna_perm_mixed()",
                     quantity = sprintf("Monte Carlo summary: %s p-values within their bounds", part), tolerance = 0,
                     key = c("from", "to"))
  }))
})
