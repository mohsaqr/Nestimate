# Skipped on CRAN to keep the check under its time limit; runs locally and in CI.
testthat::skip_on_cran()

# ---- bootstrap_network(block = ): cluster bootstrap for nested data ----

# Sessions nested in persons, long format; person-specific transition
# matrices with log-scale SD `tau`.
.sim_boot_nested <- function(n_person = 20L, n_session = 6L, len = 10L,
                             tau = 2, seed = 1) {
  set.seed(seed)
  states <- c("A", "B", "C")
  base <- matrix(c(.5, .3, .2, .2, .5, .3, .3, .2, .5), 3, byrow = TRUE)
  sessions <- lapply(seq_len(n_person), function(p) {
    P <- exp(log(base) + matrix(rnorm(9, sd = tau), 3))
    P <- P / rowSums(P)
    lapply(seq_len(n_session), function(s) {
      seq_int <- Reduce(function(prev, i) sample.int(3L, 1L, prob = P[prev, ]),
                        seq_len(len - 1L), accumulate = TRUE,
                        sample.int(3L, 1L))
      data.frame(person = p, session = s, t = seq_len(len),
                 Action = states[seq_int])
    })
  })
  build_network(do.call(rbind, unlist(sessions, recursive = FALSE)),
                method = "relative", actor = "person", session = "session",
                order = "t", action = "Action")
}

test_that("one sequence per block reproduces the ordinary bootstrap exactly", {
  net <- .sim_boot_nested(tau = 1, seed = 2)
  plain <- bootstrap_network(net, iter = 100L, seed = 3)
  single <- bootstrap_network(net, iter = 100L, seed = 3,
                              actor = ".session_label")
  expect_identical(single$mean, plain$mean)
  expect_identical(single$sd, plain$sd)
  expect_identical(single$p_values, plain$p_values)
  expect_identical(single$ci_lower, plain$ci_lower)
  expect_identical(single$n_actors, nrow(net$data))
})

test_that("nesting widens the bootstrap and is reported", {
  net <- .sim_boot_nested(tau = 2, seed = 4)
  plain <- bootstrap_network(net, iter = 300L, seed = 1)
  blocked <- bootstrap_network(net, iter = 300L, seed = 1, actor = "person")
  expect_identical(blocked$actor, "person")
  expect_identical(blocked$n_actors, 20L)
  expect_s3_class(blocked$clustering, "data.frame")
  expect_identical(names(blocked$clustering),
                   c("n_sequences", "n_actors", "icc", "icc_ci_lower",
                     "icc_ci_upper", "deff_edges"))
  expect_gt(blocked$clustering$icc, 0.1)
  expect_gt(blocked$clustering$deff_edges, 1.5)
  expect_gt(median(blocked$sd[blocked$sd > 0] / plain$sd[blocked$sd > 0]), 1.2)
  expect_identical(blocked$clustering_edges$from, blocked$summary$from)
  expect_equal(blocked$clustering_edges$deff,
               (blocked$clustering_edges$sd_actor /
                  blocked$clustering_edges$sd_sequence)^2)
  # observed network is untouched by the resampling unit
  expect_identical(blocked$original$weights, plain$original$weights)
})

test_that("blocked result prints the nesting lines; plain does not", {
  net <- .sim_boot_nested(seed = 5)
  out <- capture.output(print(bootstrap_network(net, iter = 50L, seed = 1,
                                                actor = "person")))
  expect_true(any(grepl("Actor      : person \\(20 actors\\)", out)))
  expect_true(any(grepl("Nesting    : ICC = ", out)))
  plain_out <- capture.output(print(bootstrap_network(net, iter = 50L,
                                                      seed = 1)))
  expect_false(any(grepl("Nesting", plain_out)))
  expect_null(bootstrap_network(net, iter = 20L, seed = 1)$clustering)
})

test_that("block is forwarded through grouped dispatch", {
  skip_on_cran()  # full 2000-sequence bundled data
  teams <- build_network(group_regulation_long, method = "relative",
                         actor = "Actor", action = "Action", time = "Time",
                         group = "Achiever")
  res <- bootstrap_network(teams, iter = 20L, seed = 1, actor = "Group")
  expect_s3_class(res, "net_bootstrap_group")
  expect_identical(res[[1]]$n_actors, 100L)
  expect_true(any(grepl("100 actors", capture.output(print(res)))))
})

test_that("bad block inputs raise classed errors", {
  net <- .sim_boot_nested(seed = 6)
  expect_error(bootstrap_network(net, iter = 10L, actor = "nope"),
               class = "nestimate_actor_missing")
  expect_error(bootstrap_network(net, iter = 10L, actor = c("a", "b")),
               class = "nestimate_bad_actor")
  net_one <- net
  net_one$metadata$one <- 1L
  expect_error(bootstrap_network(net_one, iter = 10L, actor = "one"),
               class = "nestimate_bad_actor")
  set.seed(8)
  cor_net <- build_network(data.frame(a = rnorm(40), b = rnorm(40),
                                      c = rnorm(40)), method = "cor")
  expect_error(bootstrap_network(cor_net, iter = 10L, actor = "id"),
               class = "nestimate_actor_unsupported")
})
