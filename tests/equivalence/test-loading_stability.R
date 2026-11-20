# loading_stability() and its summary(): the case bootstrap of the signed item weights of a build_mcml_pc() fit.
# The resamples are drawn inside the function with only a seed, so the comparison is statistical (Monte Carlo
# summaries), declared as such: boot::boot() resamples the same rows with its own draws and recomputes the
# documented weight (?build_mcml_pc, aggregation = "loadings": the item's mean absolute correlation with the other
# members of its cluster, normalised to sum to one in the cluster, signed by the leading eigenvector of the
# within-cluster block, oriented to a non-negative sum). The bootstrap mean and SD of each weight are compared with
# a tolerance fixed before the run: 4 standard errors of the difference of the two Monte Carlo estimates (a mean:
# sd * sqrt(2 / B); an SD: sd * sqrt(2 (kurtosis - 1) / (4 B)), the kurtosis read from boot's replicates). A seed
# sweep of 20 seeds before fixing these passed 20 of 20.

ls_iter <- 2000L
ls_clusters <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3", "b4"))

ls_items <- function(n, seed) {
  set.seed(seed)
  f <- stats::rnorm(n)
  g <- stats::rnorm(n)
  data.frame(a1 = f + stats::rnorm(n), a2 = 0.6 * f + stats::rnorm(n), a3 = -f + stats::rnorm(n),
             b1 = g + stats::rnorm(n), b2 = 0.5 * g + stats::rnorm(n), b3 = g + stats::rnorm(n),
             b4 = 0.3 * g + stats::rnorm(n))
}

# The signed strength weight of every item, from the correlation network of the rows given.
ls_signed_weights <- function(data, rows) {
  network <- stats::cor(data[rows, ])
  diag(network) <- 0
  unlist(lapply(ls_clusters, \(members) {
    block <- network[members, members]
    leading <- eigen(block, symmetric = TRUE)$vectors[, 1L]
    if (sum(leading) < 0) leading <- -leading
    own <- rowSums(abs(block)) / (length(members) - 1L)
    stats::setNames(own / sum(own) * ifelse(leading < 0, -1, 1), members)
  }), use.names = FALSE)
}

ls_kurtosis <- function(v) mean((v - mean(v))^4) / stats::var(v)^2

test_that("bootstrap mean and SD of the signed item weights equal boot::boot() of the same weights (Monte Carlo summary)", {
  skip_without_reference("boot")
  data <- ls_items(150L, 51L)
  items <- unlist(ls_clusters, use.names = FALSE)
  set.seed(7)
  replicates <- boot::boot(data, \(d, rows) ls_signed_weights(d, rows), R = ls_iter)$t
  means <- colMeans(replicates)
  sds <- apply(replicates, 2L, stats::sd)
  kurtosis <- apply(replicates, 2L, ls_kurtosis)
  fit <- build_mcml_pc(data, ls_clusters, aggregation = "loadings", method = "cor")
  found <- summary(loading_stability(fit, iter = ls_iter, seed = 3))
  invisible(lapply(seq_along(items), \(i) {
    expect_reference(subset(found, node == items[[i]], select = c(node, boot_mean)),
                     data.frame(node = items[[i]], boot_mean = means[[i]]), reference = "boot",
                     covers = c("loading_stability()", "summary.pc_loading_stability()"),
                     quantity = sprintf("Monte Carlo summary: bootstrap mean of the weight of %s", items[[i]]),
                     tolerance = 4 * sds[[i]] * sqrt(2 / ls_iter), key = "node",
                     note = "boot::boot() with the weight written from the documented definition")
    expect_reference(subset(found, node == items[[i]], select = c(node, boot_sd)),
                     data.frame(node = items[[i]], boot_sd = sds[[i]]), reference = "boot",
                     covers = c("loading_stability()", "summary.pc_loading_stability()"),
                     quantity = sprintf("Monte Carlo summary: bootstrap SD of the weight of %s", items[[i]]),
                     tolerance = 4 * sds[[i]] * sqrt(2 * (kurtosis[[i]] - 1) / (4 * ls_iter)), key = "node",
                     note = "boot::boot() with the weight written from the documented definition")
  }))
  # The observed weight is no draw: it is the weight of all rows.
  expect_reference(subset(found, select = c(node, weight)),
                   data.frame(node = items, weight = ls_signed_weights(data, seq_len(nrow(data)))),
                   reference = "independent: signed strength weight of the correlation network, base R",
                   covers = c("loading_stability()", "summary.pc_loading_stability()"),
                   quantity = "observed signed weight", tolerance = 1e-12, key = "node")
})
