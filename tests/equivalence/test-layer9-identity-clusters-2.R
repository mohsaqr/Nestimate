# Part 2 of test-layer9-identity-clusters.R, split so that one Prova case holds at most 7 references;
# the shared set-up is repeated.

# Layer 9, composite identity: every route that builds networks for subsets of the sequences — the groups of
# build_network(group = ), the clusters of cluster_network(), of build_network(<net_clustering>) and of
# build_network(<net_mmm>) — must return, for each subset, exactly the network build_network() fits directly on
# that subset's rows with the same settings (start/end, threshold, scaling, windows, estimator arguments, one-hot
# codes). The cluster networks of a network are built afresh with its method; with inherit = TRUE they take its
# threshold, scaling and estimator settings.

wide <- tn_simulate(n_actors = 60L, n_states = 4L, seq_length = 10L, seed = 921L)
long <- tn_long(wide)
codes <- c("Reading", "Video", "Forum", "Quiz", "Coding", "Review")
onehot <- subset(learning_activities, student <= 40)
students <- sort(unique(onehot$student))

# The wide rows of one cluster: the direct fit's input.
members <- function(assignments, cluster) wide[assignments == cluster, , drop = FALSE]

test_that("a clustered long event log equals the direct fit on each cluster's actors", {
  source <- build_network(long, method = "relative", actor = "actor", action = "action", order = "time")
  lapply(c("relative", "frequency"), \(method) {
    networks <- if (method == "relative") cluster_network(source, k = 2L) else cluster_network(source, k = 2L, method = method)
    assignments <- attr(networks, "clustering")$assignments
    lapply(1:2, \(cluster) {
      direct <- build_network(subset(long, actor %in% which(assignments == cluster)), method = method,
                              actor = "actor", action = "action", order = "time")
      li_expect_same_network(networks[[cluster]], direct, route = "build_network() on the cluster's actors",
                             covers = "cluster_network()", quantity = sprintf("cluster %d, %s", cluster, method),
                             parts = c("weights", "edges", "initial"))
    })
  })
  # a documented limitation, a design decision (owner, 2026-11-22), not a defect: the call's outcome (the class of
  # the condition it raises) is recorded beside the class of the result the relative method returns
  relative_class <- class(cluster_network(source, k = 2L, cluster_by = "mmm", seed = 1L, n_starts = 3L))[[1L]]
  frequency_outcome <- rf_outcome(cluster_network(source, k = 2L, cluster_by = "mmm", method = "frequency", seed = 1L,
                                                  n_starts = 3L), expected = relative_class)
  expect_known_difference(frequency_outcome, relative_class,
                          reference = li_ref("mixture clusters of the long event log with the relative method"),
                          covers = "cluster_network()", quantity = "frequency networks of a mixture of a long event log",
                          tolerance = 0, authority = "internal",
                          reason = "cluster_network(cluster_by = 'mmm') of a long-format network with a non-relative method passes the long column names to the wide per-cluster data (documented limitation)",
                          agrees_with = "method = 'relative' (the default) on a long-format network")
})

test_that("the clusters of a one-hot network equal the direct one-hot fit on each cluster's actors", {
  configurations <- expand.grid(method = c("relative", "frequency", "co_occurrence"), window = c(1L, 3L),
                                stringsAsFactors = FALSE)
  # a design decision (owner, 2026-11-22), not a defect
  order_reason <- "cluster networks of a one-hot network list the codes alphabetically, the one-hot network in the order of `codes`"
  Map(\(method, window) {
    # window_size 1 without sessions warns that the network is sparse; that warning is the subject of other tests
    fit_onehot <- function(rows) {
      withCallingHandlers(build_network(rows, method = method, codes = codes, actor = "student", window_size = window),
                          nestimate_sparse_onehot = \(condition) invokeRestart("muffleWarning"))
    }
    source <- fit_onehot(onehot)
    routes <- list(
      cluster_network = \() cluster_network(source, k = 2L, seed = 1L),
      build_clusters = \() build_network(build_clusters(source, k = 2L, seed = 1L)),
      build_mmm = \() build_network(build_mmm(source, k = 2L, seed = 3L, n_starts = 3L),
                                    method = if (method == "relative") "frequency" else method)
    )
    covers <- list(cluster_network = "cluster_network()", build_clusters = c("build_network()", "build_clusters()"),
                   build_mmm = c("build_network()", "build_mmm()"))
    Map(\(route, name) {
      networks <- route()
      assignments <- attr(networks, "clustering")$assignments
      rebuilt <- if (name == "build_mmm" && method == "relative") "frequency" else method
      direct <- withCallingHandlers(
        build_network(subset(onehot, student %in% students[assignments == 1L]), method = rebuilt, codes = codes,
                      actor = "student", window_size = window),
        nestimate_sparse_onehot = \(condition) invokeRestart("muffleWarning"))
      quantity <- sprintf("%s cluster 1, %s, window %d", name, method, window)
      li_expect_same_network(networks[[1L]], direct, route = "one-hot build_network() on the cluster's actors",
                             covers = covers[[name]], quantity = quantity, parts = c("weights", "edges"))
      expect_known_difference(rownames(networks[[1L]]$weights), rownames(direct$weights),
                              reference = li_ref("one-hot build_network() on the cluster's actors"), covers = covers[[name]],
                              quantity = paste(quantity, "node order"), tolerance = 0, authority = "internal",
                              reason = order_reason, agrees_with = "the codes sorted alphabetically")
    }, routes, names(routes))
  }, configurations$method, configurations$window)
})

test_that("the relative networks of build_network(<net_mmm>) are the fitted chains unless reshaped", {
  source <- build_network(wide, method = "relative", threshold = 0.3, scaling = "minmax")
  fit <- build_mmm(source, k = 2L, seed = 5L, n_starts = 5L)
  # minmax over every cell, then the threshold: given in the call, or the source network's with inherit = TRUE
  expected <- lapply(fit$models, \(model) {
    weights <- model$weights
    scaled <- (weights - min(weights)) / (max(weights) - min(weights))
    scaled[scaled < 0.3] <- 0
    scaled
  })
  routes <- list(
    list(route = "build_network()", label = "scaling and threshold given in the call",
         networks = build_network(fit, threshold = 0.3, scaling = "minmax"), expected = expected),
    list(route = "build_network()", label = "inherit = TRUE",
         networks = build_network(fit, inherit = TRUE), expected = expected),
    list(route = "cluster_network()", label = "inherit = TRUE",
         networks = cluster_network(source, k = 2L, cluster_by = "mmm", seed = 5L, n_starts = 5L, inherit = TRUE),
         expected = expected),
    list(route = "build_network()", label = "the default, the fitted chains",
         networks = build_network(fit), expected = lapply(fit$models, \(model) model$weights)),
    list(route = "cluster_network()", label = "the default, the fitted chains",
         networks = cluster_network(source, k = 2L, cluster_by = "mmm", seed = 5L, n_starts = 5L),
         expected = lapply(fit$models, \(model) model$weights))
  )
  lapply(routes, \(case) {
    Map(\(net, weights, cluster) {
      internal <- grepl("default", case$label)
      reference <- if (internal) li_ref("the fitted chain of the mixture") else
        "independent: minmax over every cell, then threshold, base R"
      expect_reference(net$weights, weights, reference = reference, authority = if (internal) "internal",
                       covers = c(case$route, "build_mmm()"),
                       quantity = sprintf("%s, %s, mixture cluster %d of a thresholded, minmax-scaled network",
                                          case$route, case$label, cluster), tolerance = 0)
    }, case$networks, case$expected, 1:2)
  })
})
