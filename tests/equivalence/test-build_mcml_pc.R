# build_mcml_pc(), summary.mcml_pc(), item_loadings(): the cluster-level (macro) network of a psychometric network
# and its item diagnostics. References: the first canonical correlation of two item blocks (stats::cancor()),
# Escoufier's (1973) RV coefficient and the block mean of the partial correlation network written from their
# definitions, and the item weights and cross-loadings from the corpcor partial correlations. New: these
# functions had no comparison with a reference.

clusters <- list(A = c("V1", "V2", "V3"), B = c("V4", "V5", "V6"), C = c("V7", "V8"))
data <- psy_continuous(300L, 8L, 1.5, 601L)
blocks <- lapply(clusters, \(members) as.matrix(data[members]))
cluster_pairs <- utils::combn(names(clusters), 2L, simplify = FALSE)

# One row per pair of clusters with the value `between(a, b)`.
psy_macro <- function(between) {
  data.frame(from = vapply(cluster_pairs, \(pair) pair[[1L]], character(1)),
             to = vapply(cluster_pairs, \(pair) pair[[2L]], character(1)),
             weight = vapply(cluster_pairs, \(pair) between(pair[[1L]], pair[[2L]]), numeric(1)))
}

# Escoufier's RV coefficient of two column-centred blocks: ||A'B||^2 / sqrt(||A'A||^2 ||B'B||^2).
psy_rv <- function(a, b) {
  a <- scale(a, scale = FALSE)
  b <- scale(b, scale = FALSE)
  sum(crossprod(a, b)^2) / sqrt(sum(crossprod(a)^2) * sum(crossprod(b)^2))
}

test_that("the canonical-correlation macro network equals stats::cancor()", {
  expected <- psy_macro(\(a, b) stats::cancor(blocks[[a]], blocks[[b]])$cor[[1L]])
  expect_reference(summary(build_mcml_pc(data, clusters = clusters, aggregation = "cancor")), expected,
                   reference = "stats", covers = c("build_mcml_pc()", "summary.mcml_pc()"),
                   quantity = "first canonical correlation", tolerance = 1e-10, key = c("from", "to"))
})

test_that("the Escoufier macro network equals the RV coefficient", {
  expect_reference(summary(build_mcml_pc(data, clusters = clusters, aggregation = "escoufier")),
                   psy_macro(\(a, b) psy_rv(blocks[[a]], blocks[[b]])),
                   reference = "independent: Escoufier's RV coefficient, base R",
                   covers = c("build_mcml_pc()", "summary.mcml_pc()"), quantity = "RV coefficient", tolerance = 1e-10,
                   key = c("from", "to"))
})

test_that("the averaged macro network and the item diagnostics follow from the corpcor partial correlations", {
  skip_without_reference("corpcor")
  pcor <- corpcor::cor2pcor(stats::cor(data))
  dimnames(pcor) <- list(names(data), names(data))
  expect_reference(summary(build_mcml_pc(data, clusters = clusters, aggregation = "average", method = "pcor")),
                   psy_macro(\(a, b) mean(pcor[clusters[[a]], clusters[[b]]])), reference = "corpcor",
                   covers = c("build_mcml_pc()", "summary.mcml_pc()"), quantity = "mean partial correlation between clusters",
                   tolerance = 1e-10, key = c("from", "to"), note = "the block mean written in the test on cor2pcor()")
  # Each item's mean absolute partial correlation with the other members of each cluster; its cross-loading is the
  # largest value for another cluster. (The item weights follow `weighting`, "equal" by default, so they are 1 / size.)
  strengths <- do.call(rbind, lapply(names(clusters), \(cluster) {
    do.call(rbind, lapply(clusters[[cluster]], \(node) {
      to <- vapply(names(clusters), \(target) mean(abs(pcor[node, setdiff(clusters[[target]], node)])), numeric(1))
      data.frame(node = node, cluster = cluster, own = to[[cluster]], max_cross = max(to[names(to) != cluster]))
    }))
  }))
  expected <- data.frame(node = strengths$node, max_cross = strengths$max_cross)
  expect_reference(item_loadings(build_mcml_pc(data, clusters = clusters, aggregation = "composite", method = "pcor")),
                   expected, reference = "corpcor", covers = c("item_loadings()", "item_loadings.mcml_pc()"),
                   quantity = "cross-loading", tolerance = 1e-10, key = "node")
})
