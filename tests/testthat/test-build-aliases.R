# build_network() is a hub: data, method, the shared arguments, and each
# method's arguments by name in `...`. Each method has an alias that names
# the arguments of its kind of network and passes on the ones its caller set.
testthat::skip_on_cran()

grl <- group_regulation_long
srl_bin <- as.data.frame(lapply(srl_strategies, function(v) {
  as.integer(v > stats::median(v))
}))
oh_codes <- c("Reading", "Video", "Forum", "Quiz", "Coding", "Review")
la <- subset(learning_activities, student <= 20)

test_that("build_network() names only what every network shares", {
  expect_identical(names(formals(build_network)),
                   c("data", "method", "...", "group", "scaling", "threshold",
                     "labels", "inherit"))
})

test_that("build_network() binds positional arguments after method as before", {
  expect_identical(
    build_network(grl, "relative", "Actor", "Action", "Time"),
    build_network(grl, method = "relative", actor = "Actor", action = "Action",
                  time = "Time"))
})

test_that("build_network() keeps a missing argument missing", {
  clustering <- build_clusters(build_network(trajectories, method = "attention",
                                             lambda = 2), k = 2, seed = 1)
  # no method given: the clustering's own
  networks <- build_network(clustering)
  expect_identical(networks[[1]]$method, "attention")
  wrapped <- function(d, cutoff) build_network(d, method = "frequency",
                                               threshold = cutoff)
  expect_error(wrapped(trajectories), "cutoff")
})

test_that("each sequence alias equals build_network() with its method", {
  expect_identical(build_ngram(trajectories, n_gram = 3L),
                   build_network(trajectories, method = "ngram",
                                 params = list(n_gram = 3L)))
  expect_identical(build_gap(trajectories, max_gap = 2L),
                   build_network(trajectories, method = "gap", max_gap = 2L))
  expect_identical(build_reverse(trajectories, weighted = TRUE),
                   build_network(trajectories, method = "reverse",
                                 weighted = TRUE))
  expect_identical(build_atna(trajectories, lambda = 2, direction = "both"),
                   build_network(trajectories, method = "attention",
                                 params = list(lambda = 2,
                                               direction = "both")))
  expect_identical(build_cna(trajectories, type = "relative"),
                   build_network(trajectories, method = "co_occurrence",
                                 params = list(type = "relative")))
})

test_that("a method setting named on the alias is stored in params", {
  expect_identical(build_atna(trajectories, lambda = 3)$params$lambda, 3)
  expect_identical(build_ngram(trajectories, n_gram = 3L)$params$n_gram, 3L)
  expect_identical(build_glasso(srl_strategies, gamma = 0.25)$params$gamma,
                   0.25)
})

test_that("build_latency() passes floor, cap, aggregate and scope on", {
  course <- subset(grl, Course == "A")
  expect_identical(
    build_latency(course, action = "Action", time = "Time", actor = "Actor",
                  floor = 0.3, cap = "60s", aggregate = "mean"),
    build_network(course, method = "latency", action = "Action",
                  time = "Time", actor = "Actor", latency_floor = 0.3,
                  latency_cap = "60s", latency_aggregate = "mean"))
  net <- build_latency(course, "Action", "Time", "Actor", floor = 0.2)
  expect_identical(net$params$latency_floor, 0.2)
})

test_that("build_wtna() equals build_network(method = 'wtna')", {
  expect_identical(
    build_wtna(la, codes = oh_codes, actor = "student", type = "relative"),
    build_network(la, method = "wtna", codes = oh_codes, actor = "student",
                  params = list(type = "relative")))
  la$part <- paste(la$student, (stats::ave(seq_len(nrow(la)), la$student,
                                           FUN = seq_along) - 1L) %/% 10L)
  by_part <- build_wtna(la, codes = oh_codes, actor = "student",
                        session = "part")
  expect_identical(unname(by_part$weights),
                   unname(build_wtna(la, codes = oh_codes,
                                     actor = c("student", "part"))$weights))
})

test_that("each association alias equals build_network() with its method", {
  expect_identical(build_cor(srl_strategies, cor_method = "spearman"),
                   build_network(srl_strategies, method = "cor",
                                 params = list(cor_method = "spearman")))
  expect_identical(build_pcor(srl_strategies, threshold = 0.05),
                   build_network(srl_strategies, method = "pcor",
                                 threshold = 0.05))
  expect_identical(build_glasso(srl_strategies, gamma = 0.25, nlambda = 50L),
                   build_network(srl_strategies, method = "glasso",
                                 params = list(gamma = 0.25, nlambda = 50L)))
  skip_if_not_installed("glmnet")
  expect_identical(build_ising(srl_bin, rule = "OR"),
                   build_network(srl_bin, method = "ising",
                                 params = list(rule = "OR")))
  expect_identical(build_mgm(srl_strategies[, 1:5], ruleReg = "OR"),
                   build_network(srl_strategies[, 1:5], method = "mgm",
                                 params = list(ruleReg = "OR")))
})

test_that("every alias names its kind's arguments and the shared ones", {
  shared <- c("group", "scaling", "threshold", "labels", "params", "...")
  sequence <- c("start", "end", "actor", "action", "time", "session", "order",
                "time_threshold", "timezone", "state_cols", "metadata_cols")
  onehot <- c("codes", "window_size", "mode")
  has <- function(fn, args) all(args %in% names(formals(fn)))
  invisible(lapply(list(build_tna, build_ftna, build_cna), function(fn) {
    expect_true(has(fn, c(sequence, onehot, shared)))
  }))
  invisible(lapply(list(build_atna, build_ngram, build_gap, build_reverse),
                   function(fn) expect_true(has(fn, c(sequence, shared)))))
  expect_true(has(build_latency, c("action", "time", "actor", "session",
                                   "floor", "cap", "aggregate", "scope")))
  expect_true(has(build_wtna, c(onehot, "type", "actor", "session")))
  invisible(lapply(list(build_cor, build_pcor, build_glasso), function(fn) {
    expect_true(has(fn, c("actor", "level", "predictability",
                          "metadata_cols", shared)))
  }))
  expect_true(has(build_ising, c("actor", "gamma", "rule", shared)))
  expect_true(has(build_mgm, c("actor", "level", "lambdaGam", "ruleReg",
                               shared)))
})

test_that("every alias refuses `method` with a classed error", {
  aliases <- list(build_ngram = trajectories, build_gap = trajectories,
                  build_reverse = trajectories, build_wtna = la,
                  build_mgm = srl_strategies)
  invisible(Map(function(name, d) {
    expect_error(getExportedValue("Nestimate", name)(d, method = "relative"),
                 class = "nestimate_bad_input")
  }, names(aliases), aliases))
  expect_error(build_latency(grl, action = "Action", time = "Time",
                             method = "reverse"),
               class = "nestimate_bad_input")
})

test_that("a latency network asks splot() for widths by rank, other methods ask nothing", {
  course <- subset(grl, Course == "A")
  net <- build_latency(course, action = "Action", time = "Time", actor = "Actor")
  expect_identical(net$meta$splot$defaults, list(edge_scale_mode = "rank"))
  expect_null(build_tna(trajectories)$meta$splot)
  grouped <- build_network(grl, method = "latency", action = "Action", time = "Time",
                           actor = "Actor", group = "Course")
  invisible(lapply(grouped, function(n) {
    expect_identical(n$meta$splot$defaults$edge_scale_mode, "rank")
  }))
})

test_that("a latency network averages its replies by default", {
  course <- subset(grl, Course == "A")
  net <- build_latency(course, action = "Action", time = "Time", actor = "Actor")
  expect_identical(net$params$latency_aggregate, "mean")
  expect_true(all(net$weights >= 0 & net$weights <= 1))
  summed <- build_latency(course, action = "Action", time = "Time", actor = "Actor",
                          aggregate = "sum")
  replied <- net$n_replies > 0
  expect_equal(summed$weights[replied] / net$n_replies[replied], net$weights[replied])
})

# Before 1.0.21 the aliases ended their positional arguments at `session`
# (`level` for cor, pcor and glasso; `actor` for ising) with `...`, so later
# positional arguments filled build_network()'s next ones: order, codes, group.
test_that("positional arguments after session bind as before the hub", {
  onehot <- data.frame(id = rep(1:3, each = 6),
                       A = rep(c(1, 0, 1, 0, 1, 0), 3),
                       B = rep(c(0, 1, 1, 0, 0, 1), 3))
  by_position <- build_tna(onehot, FALSE, FALSE, "id", NULL, NULL, NULL, NULL,
                           c("A", "B"))
  by_name <- build_tna(onehot, actor = "id", codes = c("A", "B"))
  expect_equal(by_position$weights, by_name$weights)
  expect_equal(unname(by_position$weights), matrix(0.5, 2, 2))

  wide <- data.frame(T1 = rep(c("A", "B"), each = 5),
                     T2 = rep(c("B", "A"), each = 5),
                     T3 = rep(c("C", "B"), each = 5),
                     g = rep(c("up", "down"), each = 5))
  grouped <- build_tna(wide, FALSE, FALSE, NULL, NULL, NULL, NULL, NULL, NULL,
                       "g")
  expect_s3_class(grouped, "netobject_group")
  expect_identical(names(grouped), c("up", "down"))

  dots_after <- vapply(c("build_tna", "build_ftna", "build_cna", "build_atna",
                         "build_cor", "build_pcor", "build_glasso",
                         "build_ising"), function(f) {
    args <- names(formals(getExportedValue("Nestimate", f)))
    args[match("...", args) - 1L]
  }, character(1))
  expect_identical(unname(dots_after),
                   c(rep("session", 4L), rep("level", 3L), "actor"))
})
