# bootstrap_network() on latency networks: sequences are resampled, every
# edge's weight is recomputed from the drawn replies, a percentile cap the
# network resolved from its own replies is resolved again, and a mean network
# is tested for stability on its speed score.
testthat::skip_on_cran()

lb_log <- function() {
  .local_seed(11)
  data.frame(chat = rep(sprintf("c%02d", 1:24), each = 15),
             sender = sample(c("A", "B", "C"), 360, replace = TRUE),
             at = cumsum(stats::rexp(360, 1 / 40)) + rep(seq(0, by = 1e5, length.out = 24), each = 15))
}
lb_net <- function(...) {
  build_network(lb_log(), method = "latency", actor = "chat", action = "sender", time = "at", ...)
}

test_that("a resample equals the latency estimator on the drawn sequences", {
  net <- lb_net(latency_floor = 0.2)
  .local_seed(7)
  boot <- .bootstrap_latency(net, 1L)
  .local_seed(7)
  drawn <- sample.int(nrow(net$data), nrow(net$data), replace = TRUE)
  rows <- cbind(net$data[drawn, , drop = FALSE], net$latency_time[drawn, , drop = FALSE])
  direct <- .estimator_latency(rows, alphabet = rownames(net$weights),
                               latency_floor = 0.2, latency_cap = 0.95)$matrix
  expect_equal(boot$weights[1, ], as.vector(direct), tolerance = 1e-12)
  expect_identical(boot$cap, "resolved")
})

test_that("a cap in seconds and a cap pooled over groups are held", {
  expect_identical(.bootstrap_latency(lb_net(latency_cap = "60s"), 2L)$cap, "held")
  log <- transform(lb_log(), team = rep(c("x", "y"), each = 180))
  grouped <- build_network(log, method = "latency", actor = "chat", action = "sender",
                           time = "at", group = "team")
  expect_identical(.bootstrap_latency(grouped$x, 2L)$cap, "held")
  boot <- bootstrap_network(grouped, iter = 20, seed = 1)
  expect_s3_class(boot, "net_bootstrap_group")
})

test_that("the summary adds replies, sequences and the CV of the speed score", {
  boot <- bootstrap_network(lb_net(), iter = 50, seed = 1)
  s <- summary(boot)
  expect_true(all(c("n_replies", "n_sequences", "cv") %in% names(s)))
  expect_true(all(s$n_sequences <= 24 & s$n_sequences >= 1))
  expect_identical(boot$latency_cap, "resolved")
  # a summed network keeps the test on its weights
  summed <- summary(bootstrap_network(lb_net(latency_aggregate = "sum"), iter = 20, seed = 1))
  expect_false("cv" %in% names(summed))
})

test_that("the stability test runs on the speed score, not the weight", {
  boot <- bootstrap_network(lb_net(), iter = 100, seed = 1)
  s <- summary(boot)
  floor <- 0.5
  score <- (s$weight - floor) / (1 - floor)
  # the band is 75-125% of the score, mapped back to the weight scale
  expect_equal(s$cr_lower, floor + (1 - floor) * 0.75 * score, tolerance = 1e-12)
  expect_equal(s$cr_upper, floor + (1 - floor) * 1.25 * score, tolerance = 1e-12)
})

test_that("an edge on fewer than 5 sequences is never stable, nor one stuck at the floor", {
  log <- lb_log()
  # D replies to A only in two sequences
  log$sender[c(2, 17)] <- "D"
  net <- build_network(log, method = "latency", actor = "chat", action = "sender",
                       time = "at", latency_cap = "0.01s")
  s <- summary(bootstrap_network(net, iter = 50, seed = 1))
  rare <- subset(s, from == "D" | to == "D")
  expect_true(all(rare$p_value == 1))
  # a cap of 0.01 s puts every reply at the floor: score 0, never stable
  expect_false(any(s$sig))
})

test_that("per-unit caps and actor are refused with a class", {
  expect_error(bootstrap_network(suppressMessages(lb_net(latency_scope = "actor")), iter = 10),
               class = "nestimate_latency_resampling")
  expect_error(bootstrap_network(lb_net(), iter = 10, actor = "chat"),
               class = "nestimate_actor_unsupported")
})

test_that("a latency bootstrap asks splot() for widths by rank", {
  boot <- bootstrap_network(lb_net(), iter = 20, seed = 1)
  expect_identical(boot$meta$splot$defaults$edge_scale_mode, "rank")
  expect_identical(boot$model$meta$splot$defaults$edge_scale_mode, "rank")
  # beside the width mode, the network's colours
  expect_identical(boot$meta$splot$defaults$node_fill,
                   lb_net()$meta$splot$defaults$node_fill)
})

test_that("a latency bootstrap says it is experimental, once per call", {
  expect_message(bootstrap_network(lb_net(), iter = 10, seed = 1),
                 class = "nestimate_experimental")
  log <- transform(lb_log(), team = rep(c("x", "y"), each = 180))
  grouped <- build_network(log, method = "latency", actor = "chat", action = "sender",
                           time = "at", group = "team")
  said <- 0L
  withCallingHandlers(bootstrap_network(grouped, iter = 10, seed = 1),
                      nestimate_experimental = function(m) {
                        said <<- said + 1L
                        invokeRestart("muffleMessage")
                      })
  expect_identical(said, 1L)
  expect_no_message(bootstrap_network(build_tna(trajectories), iter = 10, seed = 1),
                    class = "nestimate_experimental")
})

test_that("a grouped bootstrap keeps every group's edges in one table", {
  log <- transform(lb_log(), team = rep(c("x", "y"), each = 180))
  grouped <- build_network(log, method = "latency", actor = "chat", action = "sender",
                           time = "at", group = "team")
  boot <- suppressMessages(bootstrap_network(grouped, iter = 10, seed = 1))
  table <- attr(boot, "summary")
  expect_s3_class(table, "data.frame")
  expect_identical(names(table)[1], "group")
  expect_setequal(unique(table$group), c("x", "y"))
  expect_identical(summary(boot), table)
  # the object still iterates as its two groups
  expect_length(boot, 2L)
})

# labels = renames the weights and the matrices beside them, which the
# bootstrap summary reads by the displayed names.
test_that("a relabelled latency network bootstraps under its labels", {
  replies <- data.frame(id = rep(1:6, each = 2), state = rep(c("A", "B"), 6),
                        t = rep(c(0, 1), 6))
  net <- suppressMessages(build_latency(replies, action = "state", time = "t",
                                        actor = "id", cap = "2s",
                                        labels = c(A = "Alpha", B = "Beta")))
  expect_identical(rownames(net$n_replies), c("Alpha", "Beta"))
  expect_identical(rownames(net$median_latency), c("Alpha", "Beta"))
  edges <- suppressMessages(summary(bootstrap_network(net, iter = 20,
                                                      seed = 1)))
  expect_identical(c(edges$from, edges$to), c("Beta", "Alpha"))
  expect_equal(edges$weight, 0.75)
  expect_equal(edges$n_replies, 6)
  expect_equal(edges$n_sequences, 6)
})

# A resample of mostly instant replies resolves a percentile cap of 0, which
# build_network() refuses too; the bootstrap stops instead of dropping the
# draws as NaN.
test_that("a resample with a zero percentile cap is refused", {
  replies <- data.frame(id = rep(1:2, each = 2), state = rep(c("A", "B"), 2),
                        t = c(0, 0, 0, 8))
  net <- suppressMessages(build_latency(replies, action = "state", time = "t",
                                        actor = "id"))
  expect_error(suppressMessages(bootstrap_network(net, iter = 20, seed = 1)),
               class = "nestimate_latency_zero_cap")
  held <- suppressMessages(build_latency(replies, action = "state", time = "t",
                                         actor = "id", cap = "8s"))
  expect_s3_class(suppressMessages(bootstrap_network(held, iter = 20,
                                                     seed = 1)),
                  "net_bootstrap")
})
