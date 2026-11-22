# build_clusters() (and its old name cluster_data()), cluster_network(), cluster_choice(), cluster_diagnostics(),
# summary.net_clustering(): dissimilarity-based clustering of sequences. References: stringdist's string
# distances and tna::cluster_data() for the dissimilarities; cluster::pam() and stats::hclust() with
# cluster::silhouette() for the partition and its quality (build_clusters() itself calls cluster::pam(), so pam
# checks the hand-over); tna::tna() for the per-cluster networks; nnet::multinom() and stats::glm() for the
# covariate model. Re-expressed from local_testing_and_equivalence/test-build_clusters-equiv-extracted.R,
# test-equiv-moved-build_clusters.R and test-equiv-moved3-build_clusters.R. The distance, the assignments and the
# covariate table are read from the result's fields (API gap: no accessor returns them).

regulation <- tna::group_regulation[seq_len(50L), ]
random <- sq_random(30L, 10L, 4L, 42L)
strings <- apply(random, 1L, paste, collapse = "")

test_that("string dissimilarities equal stringdist::stringdistmatrix()", {
  skip_without_reference("stringdist")
  invisible(lapply(c("hamming", "osa", "lv", "dl", "lcs", "qgram", "cosine", "jaccard", "jw"), \(metric) {
    found <- build_clusters(random, k = 2, dissimilarity = metric, q = 2L, p = 0.1)
    expected <- stringdist::stringdistmatrix(strings, method = metric, q = 2, p = 0.1)
    expect_reference(sq_pairs(found$distance), sq_pairs(expected), reference = "stringdist", covers = "build_clusters()",
                     quantity = sprintf("%s dissimilarity", metric), tolerance = 1e-10, key = c("i", "j"))
  }))
})

test_that("dissimilarities, PAM and hierarchical partitions equal tna::cluster_data()", {
  skip_without_reference("tna")
  invisible(lapply(c("hamming", "qgram", "cosine", "jaccard"), \(metric) {
    expect_reference(sq_pairs(build_clusters(regulation, k = 2, dissimilarity = metric, q = 2L)$distance),
                     sq_pairs(tna::cluster_data(regulation, k = 2, dissimilarity = metric, q = 2)$distance),
                     reference = "tna", covers = "build_clusters()", quantity = sprintf("%s dissimilarity, group_regulation", metric),
                     tolerance = 1e-10, key = c("i", "j"))
  }))
  invisible(lapply(c(0.5, 2), \(lambda) {
    expect_reference(sq_pairs(build_clusters(regulation, k = 2, dissimilarity = "hamming", weighted = TRUE, lambda = lambda)$distance),
                     sq_pairs(tna::cluster_data(regulation, k = 2, dissimilarity = "hamming", weighted = TRUE, lambda = lambda)$distance),
                     reference = "tna", covers = "build_clusters()",
                     quantity = sprintf("weighted hamming dissimilarity, lambda %g", lambda), tolerance = 1e-10, key = c("i", "j"))
  }))
  invisible(lapply(c("pam", "complete", "average"), \(method) {
    found <- build_clusters(regulation, k = 3, dissimilarity = "hamming", method = method)
    reference <- tna::cluster_data(regulation, k = 3, dissimilarity = "hamming", method = method)
    expect_reference(data.frame(sequence = seq_along(found$assignments), cluster = found$assignments),
                     data.frame(sequence = seq_along(reference$assignments), cluster = reference$assignments),
                     reference = "tna", covers = "build_clusters()", quantity = sprintf("assignment, %s", method),
                     tolerance = 0, key = "sequence")
  }))
  old_name <- suppressWarnings(cluster_data(regulation, k = 3, dissimilarity = "hamming"))
  expect_reference(data.frame(sequence = seq_along(old_name$assignments), cluster = old_name$assignments),
                   data.frame(sequence = seq_len(nrow(regulation)),
                              cluster = tna::cluster_data(regulation, k = 3, dissimilarity = "hamming")$assignments),
                   reference = "tna", covers = "cluster_data()", quantity = "assignment, cluster_data()", tolerance = 0,
                   key = "sequence", note = "cluster_data() is the deprecated name of build_clusters() and warns")
})

test_that("PAM partition, silhouettes and within-cluster distances equal cluster::pam() and cluster::silhouette()", {
  skip_without_reference("cluster")
  found <- build_clusters(regulation, k = 3, dissimilarity = "hamming")
  distance <- as.matrix(found$distance)
  medoids <- cluster::pam(found$distance, k = 3, diss = TRUE)
  widths <- cluster::silhouette(medoids$clustering, found$distance)
  expect_reference(data.frame(sequence = seq_along(found$assignments), cluster = found$assignments),
                   data.frame(sequence = seq_along(medoids$clustering), cluster = unname(medoids$clustering)),
                   reference = "cluster", covers = "build_clusters()", quantity = "PAM assignment", tolerance = 0,
                   key = "sequence")
  expect_reference(found$silhouette, mean(widths[, "sil_width"]), reference = "cluster", covers = "build_clusters()",
                   quantity = "average silhouette width", tolerance = 1e-12)
  members <- split(seq_len(nrow(distance)), medoids$clustering)
  # mean distance between distinct members of a cluster
  within <- vapply(members, \(m) mean(distance[m, m][upper.tri(diag(length(m)))]), numeric(1))
  expected <- data.frame(cluster = as.integer(names(members)), size = lengths(members), mean_within_dist = within)
  expect_reference(summary(found), expected, reference = "cluster", covers = c("build_clusters()", "summary.net_clustering()"),
                   quantity = "cluster size and mean within-cluster distance", tolerance = 1e-12, key = "cluster",
                   note = "the mean of the pairwise distances written in the test on cluster::pam()'s partition")
  silhouettes <- data.frame(cluster = as.integer(names(members)),
                            sil_mean = as.vector(tapply(widths[, "sil_width"], widths[, "cluster"], mean)),
                            mean_within_dist = within)
  expect_reference(as.data.frame(cluster_diagnostics(found)), silhouettes, reference = "cluster",
                   covers = c("cluster_diagnostics()", "cluster_diagnostics.net_clustering()",
                              "as.data.frame.net_cluster_diagnostics()"),
                   quantity = "per-cluster silhouette and within distance", tolerance = 1e-12, key = "cluster")
  expect_reference(summary(cluster_diagnostics(found)), silhouettes, reference = "cluster",
                   covers = c("cluster_diagnostics()", "summary.net_cluster_diagnostics()"),
                   quantity = "per-cluster silhouette, summary()", tolerance = 1e-12, key = "cluster")
})

test_that("the cluster_choice() sweep equals stats::hclust(ward.D2) with cluster::silhouette()", {
  skip_without_reference("cluster")
  distance <- build_clusters(regulation, k = 2, dissimilarity = "hamming")$distance
  tree <- stats::hclust(distance, method = "ward.D2")
  expected <- data.frame(k = 2:4, silhouette = vapply(2:4, \(k) mean(cluster::silhouette(stats::cutree(tree, k), distance)[, "sil_width"]),
                                                      numeric(1)),
                         min_size = vapply(2:4, \(k) min(table(stats::cutree(tree, k))), numeric(1)),
                         max_size = vapply(2:4, \(k) max(table(stats::cutree(tree, k))), numeric(1)))
  found <- cluster_choice(regulation, k = 2:4)
  expect_reference(found, expected, reference = "cluster", covers = "cluster_choice()",
                   quantity = "silhouette and cluster sizes per k", tolerance = 1e-12, key = "k")
  expect_reference(summary(found), expected, reference = "cluster", covers = c("cluster_choice()", "summary.cluster_choice()"),
                   quantity = "silhouette and cluster sizes per k, summary()", tolerance = 1e-12, key = "k")
})

test_that("per-cluster transition networks equal tna::tna() on each cluster's sequences", {
  skip_without_reference("tna")
  data <- tna::group_regulation[seq_len(60L), ]
  networks <- cluster_network(data, k = 2)
  assignments <- build_clusters(data, k = 2)$assignments
  invisible(lapply(1:2, \(cluster) {
    expect_reference(extract_transition_matrix(networks[[sprintf("Cluster %d", cluster)]]),
                     tna::tna(data[assignments == cluster, ])$weights, reference = "tna", covers = "cluster_network()",
                     quantity = sprintf("transition probabilities, cluster %d", cluster), tolerance = 1e-12)
  }))
  expect_reference(as.data.frame(cluster_diagnostics(networks)),
                   data.frame(cluster = 1:2, size = as.vector(table(tna::cluster_data(data, k = 2)$assignments))),
                   reference = "tna", covers = c("cluster_diagnostics()", "cluster_diagnostics.netobject_group()"),
                   quantity = "cluster sizes of a cluster_network() result", tolerance = 0, key = "cluster")
})

test_that("covariate coefficients equal nnet::multinom() (k = 3) and stats::glm() (k = 2)", {
  skip_without_reference("nnet")
  set.seed(77)
  data <- data.frame(sq_random(50L, 3L, 3L, 77L), X1 = stats::rnorm(50), X2 = stats::rnorm(50))
  found <- suppressWarnings(build_clusters(data, k = 3, covariates = c("X1", "X2"), estimator = "multinom"))
  fit <- nnet::multinom(cluster ~ X1 + X2, data = data.frame(cluster = factor(found$assignments), data[c("X1", "X2")]),
                        trace = FALSE)
  fitted <- summary(fit)
  estimate <- sq_cells(fitted$coefficients, "estimate")
  se <- sq_cells(fitted$standard.errors, "std_error")
  z <- estimate$estimate / se$std_error
  expected <- data.frame(cluster = estimate$from, variable = estimate$to, estimate = estimate$estimate,
                         std_error = se$std_error, odds_ratio = exp(estimate$estimate), z = z,
                         p = 2 * stats::pnorm(-abs(z)))
  # The same nnet fit on both sides; 1e-5 relative is the original test's tolerance (BFGS stopping rule).
  expect_reference(found$covariates$coefficients, expected, reference = "nnet", covers = "build_clusters()",
                   quantity = "multinomial covariate model, k = 3", tolerance = 1e-5, measure = "relative",
                   key = c("cluster", "variable"))
  two <- suppressWarnings(build_clusters(data, k = 2, covariates = c("X1", "X2"), estimator = "multinom"))
  table_glm <- summary(stats::glm(cluster ~ X1 + X2, family = stats::binomial,
                                  data = data.frame(cluster = factor(two$assignments), data[c("X1", "X2")])))$coefficients
  expected_two <- data.frame(variable = rownames(table_glm), estimate = table_glm[, "Estimate"],
                             std_error = table_glm[, "Std. Error"], z = table_glm[, "z value"], p = table_glm[, "Pr(>|z|)"])
  # Maximum likelihood by nnet's BFGS (the estimator asked for) against glm's IRLS: nnet stops at its relative
  # change 1e-8, so its estimates and Hessian-based standard errors differ from glm's in the fifth digit;
  # nnet::multinom() run directly on the same data shows the same difference from glm (established on the Mac).
  expect_known_difference(two$covariates$coefficients, expected_two, reference = "stats", covers = "build_clusters()",
                          reason = paste("estimator = 'multinom' fits by nnet's BFGS, which stops at a relative change",
                                         "of 1e-8; glm's IRLS converges further, so standard errors and z differ",
                                         "by about 1e-5 relative. nnet::multinom() itself differs from glm the same way."),
                          agrees_with = "nnet::multinom() (the comparison below)",
                          quantity = "binomial covariate model, k = 2", tolerance = 1e-5, measure = "relative", key = "variable")
  binomial_nnet <- summary(nnet::multinom(cluster ~ X1 + X2, trace = FALSE,
                                          data = data.frame(cluster = factor(two$assignments), data[c("X1", "X2")])))
  expect_known_difference(two$covariates$coefficients,
                          data.frame(variable = names(binomial_nnet$coefficients), estimate = unname(binomial_nnet$coefficients),
                                     std_error = unname(binomial_nnet$standard.errors)),
                          reference = "nnet", covers = "build_clusters()", quantity = "binomial covariate model, k = 2, nnet",
                          tolerance = 1e-10, measure = "relative", key = "variable",
                          reason = paste("build_clusters() fits nnet::multinom() on standardised covariates and transforms",
                                         "back; BFGS from that start stops (relative change 1e-8) at a different point",
                                         "than nnet on the raw covariates: estimates differ by about 4e-6 relative.",
                                         "The k = 3 comparison above agrees within 1e-5."))
})
