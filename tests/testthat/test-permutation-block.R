# ---- permutation(block = ): cluster-level permutation for nested data ----

# Sessions nested in persons, long format. Each person has an idiosyncratic
# transition matrix (log-scale deviations with SD `tau`), so sessions from one
# person are correlated. `design = "between"`: condition is a person label;
# `"within"`: condition alternates over a person's sessions. Groups never
# differ, so H0 is true.
.sim_nested_long <- function(n_person = 20L, n_session = 6L, len = 10L,
                             tau = 2, design = "between", seed = 1) {
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
      cond <- if (design == "between") p %% 2L else s %% 2L
      data.frame(person = p, session = s, t = seq_len(len),
                 Action = states[seq_int], cond = paste0("g", cond))
    })
  })
  do.call(rbind, unlist(sessions, recursive = FALSE))
}

.nested_net <- function(...) {
  build_network(.sim_nested_long(...), method = "relative", actor = "person",
                session = "session", order = "t", action = "Action",
                group = "cond")
}


# ---- the permutation scheme itself (property tests) ----

test_that("block scheme keeps pure blocks whole and crossed counts fixed", {
  ids <- c("p1", "p1", "p2", "p3", "p3", "p4", "p4", "p4", "p5", "p5")
  is_x <- c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
  # 3 pure arrangements x 3 x 2 crossed = 18 < 1 / alpha
  expect_warning(design <- .block_design(ids, is_x, block = "b", alpha = 0.05),
                 class = "nestimate_few_blocks")
  set.seed(1)
  draws <- replicate(2000L, .block_permute(design))

  # pure blocks (p1, p2, p3) move as units
  expect_true(all(draws[1, ] == draws[2, ]))
  expect_true(all(draws[4, ] == draws[5, ]))
  # two of the three pure blocks are always in x
  expect_true(all(draws[1, ] + draws[3, ] + draws[4, ] == 2L))
  # crossed blocks keep their per-block x counts (p4: 2 of 3, p5: 1 of 2)
  expect_true(all(colSums(draws[6:8, ]) == 2L))
  expect_true(all(colSums(draws[9:10, ]) == 1L))
  # ... and are actually shuffled, uniformly (2/3 and 1/2 marginally)
  expect_equal(mean(draws[8, ]), 2 / 3, tolerance = 0.05)
  expect_equal(mean(draws[10, ]), 1 / 2, tolerance = 0.05)
  # every pure arrangement occurs (choose(3, 2) = 3)
  expect_length(unique(apply(draws[c(1, 3, 4), ], 2, paste, collapse = "")), 3L)
})

test_that("singleton blocks reduce to an ordinary label shuffle", {
  is_x <- rep(c(TRUE, FALSE), c(6L, 4L))
  design <- .block_design(as.character(seq_along(is_x)), is_x,
                          block = "b", alpha = 0.05)
  set.seed(2)
  draws <- replicate(500L, .block_permute(design))
  expect_true(all(colSums(draws) == 6L))
  expect_equal(rowMeans(draws), rep(0.6, 10L), tolerance = 0.1)
})


# ---- permutation() surface ----

test_that("block changes the null, not the observed differences", {
  net <- .nested_net(seed = 3)
  plain <- permutation(net, iter = 99L, seed = 1)[[1]]
  blocked <- permutation(net, iter = 99L, seed = 1, block = "person")[[1]]
  expect_identical(blocked$diff, plain$diff)
  expect_identical(blocked$block, "person")
  expect_identical(blocked$n_blocks, 20L)
  expect_null(plain$block)
  expect_null(plain$n_blocks)
  expect_output(print(blocked), "Blocked by: person \\(20 blocks\\)")
})

test_that("blocked result carries tidy clustering tables and prints them", {
  net <- .nested_net(tau = 2, seed = 10)
  res <- permutation(net, iter = 99L, seed = 1, block = "person")[[1]]
  expect_s3_class(res$clustering, "data.frame")
  expect_identical(nrow(res$clustering), 1L)
  expect_identical(names(res$clustering), c(
    "n_sequences", "n_blocks", "design", "icc", "icc_ci_lower",
    "icc_ci_upper", "deff_edges", "deff_global", "min_p"))
  expect_identical(res$clustering$design, "between")
  expect_gt(res$clustering$icc, 0.1)
  expect_gt(res$clustering$deff_global, 1.5)
  expect_identical(res$clustering_edges$from, res$summary$from)
  expect_equal(res$clustering_edges$deff,
               (res$clustering_edges$null_sd_blocked /
                  res$clustering_edges$null_sd_unblocked)^2)
  out <- capture.output(print(res))
  expect_true(any(grepl("Global test", out)))
  expect_true(any(grepl("Clustering by person: ICC = ", out)))
  expect_true(any(grepl("Design effect", out)))
  # unblocked results print the global test but no clustering lines
  plain_out <- capture.output(print(permutation(net, iter = 20L, seed = 1)[[1]]))
  expect_true(any(grepl("Global test", plain_out)))
  expect_false(any(grepl("Clustering", plain_out)))
  expect_null(permutation(net, iter = 20L, seed = 1)[[1]]$clustering)
})

test_that("grouped print shows every pair in full", {
  net <- .nested_net(seed = 11)
  out <- capture.output(print(permutation(net, iter = 20L, seed = 1,
                                          block = "person")))
  expect_true(any(grepl("Grouped Permutation", out)))
  expect_true(any(grepl("^-- .* vs .* --$", out)))
  expect_true(any(grepl("Blocked by: person", out)))
})

test_that("blocking removes the false positives of clustered between-person data", {
  net <- .nested_net(tau = 2, design = "between", seed = 4)
  plain <- permutation(net, iter = 199L, seed = 1)[[1]]
  blocked <- permutation(net, iter = 199L, seed = 1, block = "person")[[1]]
  # H0 is true: the plain test is fooled by the person effect, the blocked
  # test is not
  expect_lt(plain$global$p_value[1], 0.05)
  expect_gt(blocked$global$p_value[1], 0.05)
  expect_gt(mean(blocked$p_values), mean(plain$p_values))
})

test_that("within-person designs use within-block shuffling", {
  net <- .nested_net(tau = 2, design = "within", seed = 5)
  blocked <- permutation(net, iter = 99L, seed = 1, block = "person")[[1]]
  expect_identical(blocked$n_blocks, 20L)
  expect_true(all(blocked$p_values > 0 & blocked$p_values <= 1))
})

test_that("block is forwarded through grouped dispatch", {
  net <- .nested_net(seed = 6)
  res <- permutation(net, iter = 20L, seed = 1, block = "person")
  expect_s3_class(res, "net_permutation_group")
  expect_identical(res[[1]]$block, "person")
  both <- permutation(net, net, iter = 20L, seed = 1, block = "person")
  expect_identical(both[[1]]$block, "person")
})

test_that("block works with a team label on bundled data", {
  net <- build_network(group_regulation_long, method = "relative",
                       actor = "Actor", action = "Action", time = "Time",
                       group = "Achiever")
  res <- permutation(net, iter = 20L, seed = 1, block = "Group")[[1]]
  expect_identical(res$n_blocks, 200L)
})


# ---- error and warning paths, by class ----

test_that("bad block arguments raise classed errors", {
  net <- .nested_net(seed = 7)
  expect_error(permutation(net, iter = 10L, block = "nope"),
               class = "nestimate_block_missing")
  expect_error(permutation(net, iter = 10L, block = c("a", "b")),
               class = "nestimate_bad_block")
  expect_error(permutation(net, iter = 10L, block = 1),
               class = "nestimate_bad_block")
  expect_error(permutation(net[[1]], net[[1]], iter = 10L, paired = TRUE,
                           block = "person"),
               class = "nestimate_bad_block")
})

test_that("block is refused for association networks", {
  set.seed(8)
  d <- data.frame(a = rnorm(60), b = rnorm(60), c = rnorm(60))
  n1 <- build_network(d, method = "cor")
  n2 <- build_network(d, method = "cor")
  expect_error(permutation(n1, n2, iter = 10L, block = "id"),
               class = "nestimate_block_unsupported")
})

test_that("too few blocks to reach alpha warns", {
  net <- .nested_net(n_person = 4L, seed = 9)
  expect_warning(permutation(net, iter = 10L, block = "person"),
                 class = "nestimate_few_blocks")
})
