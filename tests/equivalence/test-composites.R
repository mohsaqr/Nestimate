# composites() and as_networks(): the cluster scores of a re-estimated build_mcml_pc() fit and the group of its
# macro and within-cluster networks. The scores are written from the documented definition (?build_mcml_pc,
# ?composites): each item standardised ("scaled", "loadings") or raw ("mean", "median"), signed by the leading
# eigenvector of the within-cluster block of the node-level network (oriented so the vector sums to a non-negative
# number), weighted equally ("scaled", "mean") or by its mean absolute connection to the other members of its
# cluster ("loadings": the strength weights, normalised to sum to one), and averaged over the items observed in
# the row with the absolute weights renormalised; "median" is the row median of the signed raw items. The
# networks are compared with stats::cor() and corpcor::cor2pcor(). New: these functions had no comparison.

# Three item blocks: a reverse-keyed item in A, missing values in A only, one cluster of two items.
pc_items <- function(n, seed) {
  set.seed(seed)
  f <- stats::rnorm(n)
  g <- stats::rnorm(n)
  h <- stats::rnorm(n)
  data <- data.frame(a1 = f + stats::rnorm(n), a2 = f + stats::rnorm(n), a3 = -f + stats::rnorm(n),
                     b1 = g + stats::rnorm(n), b2 = g + stats::rnorm(n), b3 = g + stats::rnorm(n),
                     c1 = h + stats::rnorm(n), c2 = h + stats::rnorm(n))
  data$a2[sample.int(n, 4L)] <- NA
  data$a1[sample.int(n, 2L)] <- NA
  data
}
pc_clusters <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"), C = c("c1", "c2"))

# Item signs and weights of one cluster from the node-level correlation network (zero diagonal, complete rows).
pc_signs_weights <- function(data, members, strength) {
  network <- stats::cor(data[stats::complete.cases(data), ])
  diag(network) <- 0
  block <- network[members, members]
  leading <- eigen(block, symmetric = TRUE)$vectors[, 1L]
  if (sum(leading) < 0) leading <- -leading
  own <- rowSums(abs(block)) / (length(members) - 1L)
  list(sign = ifelse(leading < 0, -1, 1),
       weight = if (strength) own / sum(own) else rep(1 / length(members), length(members)))
}

# The documented score of every row on every cluster.
pc_scores <- function(data, clusters, aggregation) {
  scores <- lapply(clusters, \(members) {
    items <- as.matrix(data[members])
    if (aggregation %in% c("scaled", "loadings")) items <- base::scale(items)
    sw <- pc_signs_weights(data, members, strength = identical(aggregation, "loadings"))
    signed <- sweep(items, 2L, sw$sign, "*")
    if (identical(aggregation, "median")) return(apply(signed, 1L, stats::median, na.rm = TRUE))
    weights <- matrix(sw$weight, nrow(items), ncol(items), byrow = TRUE)
    weights[is.na(items)] <- 0
    rowSums(ifelse(is.na(signed), 0, signed) * weights) / rowSums(weights)
  })
  data.frame(row = seq_len(nrow(data)), scores)
}

# A matrix's upper triangle as a from/to/weight table.
pc_upper <- function(m) {
  cells <- which(upper.tri(m), arr.ind = TRUE)
  data.frame(from = rownames(m)[cells[, 1L]], to = colnames(m)[cells[, 2L]], weight = m[cells])
}

# An undirected netobject's weights as the same table (the zero diagonal left out).
pc_network_upper <- function(network) pc_upper(extract_transition_matrix(network))

test_that("cluster scores equal the documented sign-corrected weighted means, four aggregations, two data sets", {
  invisible(lapply(list(list(n = 200L, seed = 31L), list(n = 120L, seed = 32L)), \(input) {
    data <- pc_items(input$n, input$seed)
    invisible(lapply(c("scaled", "mean", "median", "loadings"), \(aggregation) {
      fit <- suppressWarnings(build_mcml_pc(data, pc_clusters, aggregation = aggregation, method = "cor"))
      expect_reference(data.frame(row = seq_len(nrow(data)), composites(fit)), pc_scores(data, pc_clusters, aggregation),
                       reference = "independent: sign-corrected weighted item mean, base R",
                       covers = c("composites()", "composites.mcml_pc()"),
                       quantity = sprintf("cluster score, %s, n = %d", aggregation, input$n), tolerance = 1e-10,
                       key = "row")
    }))
  }))
})

test_that("the item signs found agree with the sign of the reverse-keyed item", {
  data <- pc_items(200L, 31L)
  fit <- suppressWarnings(build_mcml_pc(data, pc_clusters, aggregation = "scaled", method = "cor"))
  signs <- do.call(rbind, lapply(names(pc_clusters), \(cluster) {
    data.frame(node = pc_clusters[[cluster]], sign = pc_signs_weights(data, pc_clusters[[cluster]], FALSE)$sign)
  }))
  expect_reference(item_loadings(fit), signs, reference = "independent: sign of the leading eigenvector, base R",
                   covers = c("item_loadings()", "item_loadings.mcml_pc()"), quantity = "item sign", tolerance = 0,
                   key = "node")
})

test_that("as_networks() holds the macro network of the scores and the within-cluster networks of the items", {
  data <- pc_items(200L, 31L)
  fit <- suppressWarnings(build_mcml_pc(data, pc_clusters, aggregation = "scaled", method = "cor"))
  networks <- as_networks(fit)
  scores <- composites(fit)
  expect_reference(pc_network_upper(networks$macro), pc_upper(stats::cor(scores, use = "complete.obs")),
                   reference = "stats", covers = c("as_networks()", "as_networks.mcml_pc()"),
                   quantity = "macro network: correlation of the cluster scores", tolerance = 1e-12, key = c("from", "to"))
  invisible(lapply(names(pc_clusters), \(cluster) {
    items <- data[pc_clusters[[cluster]]]
    expect_reference(pc_network_upper(networks[[cluster]]),
                     pc_upper(stats::cor(items[stats::complete.cases(items), ])),
                     reference = "stats", covers = c("as_networks()", "as_networks.mcml_pc()"),
                     quantity = sprintf("within-cluster correlation network, %s", cluster), tolerance = 1e-12,
                     key = c("from", "to"))
  }))
})

test_that("as_networks() of a partial-correlation fit holds corpcor::cor2pcor() of each cluster's items", {
  skip_without_reference("corpcor")
  data <- pc_items(200L, 33L)
  fit <- suppressWarnings(build_mcml_pc(data, pc_clusters, aggregation = "mean", method = "pcor"))
  networks <- as_networks(fit)
  invisible(lapply(names(pc_clusters)[1:2], \(cluster) {
    items <- data[pc_clusters[[cluster]]]
    partial <- corpcor::cor2pcor(stats::cor(items[stats::complete.cases(items), ]))
    dimnames(partial) <- list(names(items), names(items))
    expect_reference(pc_network_upper(networks[[cluster]]), pc_upper(partial), reference = "corpcor",
                     covers = c("as_networks()", "as_networks.mcml_pc()"),
                     quantity = sprintf("within-cluster partial correlation network, %s", cluster), tolerance = 1e-10,
                     key = c("from", "to"))
  }))
})

test_that("as_networks() returns a netobject_group as it is: its networks equal stats::cor() of each group", {
  data <- pc_items(200L, 34L)
  complete <- data[stats::complete.cases(data), ]
  group <- rep(c("x", "y"), length.out = nrow(complete))
  networks <- as_networks(build_network(cbind(complete, grp = group), method = "cor", group = "grp"))
  invisible(lapply(c("x", "y"), \(g) {
    expect_reference(pc_network_upper(networks[[g]]), pc_upper(stats::cor(complete[group == g, ])),
                     reference = "stats", covers = c("as_networks()", "as_networks.default()"),
                     quantity = sprintf("group %s correlation network", g), tolerance = 1e-12, key = c("from", "to"))
  }))
})

# Outside reference: psych::scoreItems() with the keys of the design (a3 reverse-keyed), the mean of the items
# answered in each row (impute = "none"). psych reverses an item as min + max - x; with min = -max that is -x, the
# sign correction of composites().
test_that("mean and scaled cluster scores equal psych::scoreItems() with the design's keys", {
  skip_without_reference("psych")
  keys <- list(A = c("a1", "a2", "-a3"), B = c("b1", "b2", "b3"), C = c("c1", "c2"))
  invisible(lapply(list(list(n = 200L, seed = 31L), list(n = 120L, seed = 32L)), \(input) {
    data <- pc_items(input$n, input$seed)
    invisible(lapply(c("mean", "scaled"), \(aggregation) {
      items <- if (identical(aggregation, "scaled")) as.data.frame(base::scale(data)) else data
      bound <- 10 * max(abs(items), na.rm = TRUE)
      scored <- psych::scoreItems(keys, items, totals = FALSE, missing = TRUE, impute = "none",
                                  min = -bound, max = bound)
      fit <- suppressWarnings(build_mcml_pc(data, pc_clusters, aggregation = aggregation, method = "cor"))
      expect_reference(data.frame(row = seq_len(nrow(data)), composites(fit)),
                       data.frame(row = seq_len(nrow(data)), as.data.frame(unclass(scored$scores))),
                       reference = "psych", covers = c("composites()", "composites.mcml_pc()"),
                       quantity = sprintf("cluster score, %s, n = %d, psych scoreItems()", aggregation, input$n),
                       tolerance = 1e-10, key = "row")
    }))
  }))
})
