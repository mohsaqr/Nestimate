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

test_that("each group of build_network(group = ) equals the direct fit on that group's rows", {
  grouped <- transform(wide, g = rep(c("x", "y", "z"), 20L))
  settings <- list(
    list(quantity = "relative", args = list(method = "relative")),
    list(quantity = "frequency with start, threshold 1 and windows of 2",
         args = list(method = "frequency", start = TRUE, threshold = 1, params = list(window_size = 2L))),
    list(quantity = "relative with end and minmax scaling", args = list(method = "relative", end = TRUE, scaling = "minmax")),
    list(quantity = "attention, lambda 2", args = list(method = "attention", lambda = 2))
  )
  lapply(settings, \(setting) {
    fit <- do.call(build_network, c(list(grouped, group = "g"), setting$args))
    lapply(c("x", "y", "z"), \(level) {
      rows <- subset(grouped, g == level, select = -g)
      direct <- do.call(build_network, c(list(rows), setting$args))
      # the grouped route records the positions and the shared alphabet it resolved once for all groups
      li_expect_same_network(fit[[level]], direct, route = "build_network() on the group's rows",
                             covers = "build_network()", quantity = sprintf("group %s, %s", level, setting$quantity),
                             drop_params = c("cols", "alphabet"))
    })
  })
})

test_that("each group of a long event log equals the direct fit on that group's events", {
  fit <- build_network(group_regulation_long, method = "relative", actor = "Actor", action = "Action",
                       time = "Time", group = "Course")
  lapply(c("A", "B", "C"), \(course) {
    direct <- build_network(subset(group_regulation_long, Course == course), method = "relative",
                            actor = "Actor", action = "Action", time = "Time")
    li_expect_same_network(fit[[course]], direct, route = "build_network() on the course's events",
                           covers = "build_network()", quantity = sprintf("course %s", course),
                           drop_params = c("cols", "alphabet"))
  })
  glasso_data <- transform(psy_continuous(n = 300L, p = 5L, rho = 0.8, seed = 922L), g = rep(c("a", "b"), 150L))
  grouped_glasso <- build_network(glasso_data, method = "glasso", group = "g")
  lapply(c("a", "b"), \(level) {
    li_expect_same_network(grouped_glasso[[level]], build_network(subset(glasso_data, g == level, select = -g), method = "glasso"),
                           route = "build_network() on the group's rows", covers = "build_network()",
                           quantity = sprintf("glasso group %s", level))
  })
})

test_that("cluster_network() on sequences equals the direct fit on each cluster's rows", {
  cases <- list(
    list(quantity = "PAM, relative", call = list(cluster_by = "pam"), direct = list(method = "relative")),
    list(quantity = "ward.D2, frequency with threshold 2 and start in the call",
         call = list(cluster_by = "ward.D2", method = "frequency", threshold = 2, start = TRUE),
         direct = list(method = "frequency", threshold = 2, start = TRUE)),
    list(quantity = "PAM, attention lambda 3 in the call", call = list(cluster_by = "pam", method = "attention", params = list(lambda = 3)),
         direct = list(method = "attention", params = list(lambda = 3)))
  )
  lapply(cases, \(case) {
    networks <- do.call(cluster_network, c(list(wide, k = 2L), case$call))
    assignments <- attr(networks, "clustering")$assignments
    lapply(1:2, \(cluster) {
      direct <- do.call(build_network, c(list(members(assignments, cluster)), case$direct))
      li_expect_same_network(networks[[cluster]], direct, route = "build_network() on the cluster's rows",
                             covers = "cluster_network()", quantity = sprintf("cluster %d, %s", cluster, case$quantity))
    })
  })
})

# A source network and what a fresh rebuild (inherit = FALSE, the default) and an inheriting one (inherit = TRUE) of
# its clusters must equal: build_network() on the cluster's rows with the network's method only, or with all of the
# network's settings. Windows on sequence data travel as attributes of the data, so both rebuilds count in them, but
# a fresh rebuild does not record them in its $params (`fresh_drop`; an open question for the owner).
sources <- list(
  list(quantity = "frequency with start, end and windows of 2",
       args = list(method = "frequency", start = TRUE, end = TRUE, params = list(window_size = 2L)),
       fresh = list(method = "frequency", params = list(window_size = 2L)), fresh_drop = "window_size"),
  list(quantity = "attention, lambda 3", args = list(method = "attention", lambda = 3),
       fresh = list(method = "attention")),
  list(quantity = "relative with start and end", args = list(method = "relative", start = TRUE, end = TRUE),
       fresh = list(method = "relative")),
  list(quantity = "relative, threshold 0.2 and minmax scaling",
       args = list(method = "relative", threshold = 0.2, scaling = "minmax"), fresh = list(method = "relative"))
)
rebuilds <- list(list(inherit = FALSE, expected = "fresh", label = "the network's method only"),
                 list(inherit = TRUE, expected = "args", label = "the source network's settings"))

test_that("cluster_network() and build_network(<net_clustering>) rebuild each cluster fresh, or with inherit = TRUE", {
  lapply(sources, \(setting) {
    source <- do.call(build_network, c(list(wide), setting$args))
    lapply(rebuilds, \(rebuild) {
      clustering <- build_clusters(source, k = 2L)
      routes <- list(cluster_network = cluster_network(source, k = 2L, inherit = rebuild$inherit),
                     build_network = build_network(clustering, inherit = rebuild$inherit))
      Map(\(networks, route) {
        assignments <- attr(networks, "clustering")$assignments
        lapply(1:2, \(cluster) {
          direct <- do.call(build_network, c(list(members(assignments, cluster)), setting[[rebuild$expected]]))
          li_expect_same_network(networks[[cluster]], direct,
                                 route = sprintf("build_network() on the cluster's rows, %s", rebuild$label),
                                 drop_params = if (rebuild$inherit) character() else setting$fresh_drop %||% character(),
                                 covers = if (route == "cluster_network") "cluster_network()" else c("build_network()", "build_clusters()"),
                                 quantity = sprintf("%s, inherit = %s, cluster %d of a %s network", route, rebuild$inherit,
                                                    cluster, setting$quantity))
        })
      }, routes, names(routes))
    })
  })
})

test_that("build_network(<net_mmm>) rebuilds each cluster fresh, or with inherit = TRUE", {
  # A rebuilt frequency or relative network takes the mixture's own initial probabilities (documented), so the
  # initial probabilities are not part of this identity.
  parts <- c("weights", "edges", "settings")
  windowed <- build_network(wide, method = "relative", params = list(window_size = 2L))
  fit <- build_mmm(windowed, k = 2L, seed = 5L, n_starts = 5L)
  lapply(1:2, \(cluster) {
    li_expect_same_network(build_network(fit, method = "frequency")[[cluster]],
                           build_network(members(fit$assignments, cluster), method = "frequency", params = list(window_size = 2L)),
                           route = "build_network() on the cluster's rows, windows in the data", covers = c("build_network()", "build_mmm()"),
                           quantity = sprintf("mixture cluster %d, frequency of a windowed network", cluster), parts = parts,
                           drop_params = "window_size")
  })
  mixtures <- list(
    list(quantity = "attention network, lambda 3", args = list(method = "attention", lambda = 3),
         fresh = list(method = "attention")),
    list(quantity = "frequency network, threshold 22", args = list(method = "frequency", threshold = 22),
         fresh = list(method = "frequency"))
  )
  lapply(mixtures, \(setting) {
    source <- do.call(build_network, c(list(wide), setting$args))
    mixture <- build_mmm(source, k = 2L, seed = 5L, n_starts = 5L)
    lapply(rebuilds, \(rebuild) {
      networks <- build_network(mixture, inherit = rebuild$inherit)
      lapply(1:2, \(cluster) {
        li_expect_same_network(networks[[cluster]],
                               do.call(build_network, c(list(members(mixture$assignments, cluster)), setting[[rebuild$expected]])),
                               route = sprintf("build_network() on the cluster's rows, %s", rebuild$label),
                               covers = c("build_network()", "build_mmm()"),
                               quantity = sprintf("inherit = %s, mixture cluster %d of a %s", rebuild$inherit, cluster, setting$quantity),
                               parts = parts)
      })
    })
  })
})

test_that("cluster_network(cluster_by = 'mmm') equals build_network() of the same mixture fit", {
  source <- build_network(wide, method = "frequency")
  one_step <- cluster_network(source, k = 2L, cluster_by = "mmm", seed = 5L, n_starts = 5L)
  two_step <- build_network(build_mmm(source, k = 2L, seed = 5L, n_starts = 5L))
  Map(\(a, e, cluster) li_expect_same_network(a, e, route = "build_network(build_mmm())", covers = c("cluster_network()", "build_mmm()"),
                                            quantity = sprintf("mixture cluster %d", cluster)),
      one_step, two_step, 1:2)
})
