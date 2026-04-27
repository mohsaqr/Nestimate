# Smoke + identity tests for chain_structure().
# Numerical equivalence vs markovchain lives in
# local_testing_and_equivalence/test-equiv-chain_structure.R.

test_that("chain_structure on a 3-cycle has period 3 and is not aperiodic", {
  P <- matrix(0, 3, 3, dimnames = list(c("A","B","C"), c("A","B","C")))
  P["A","B"] <- 1; P["B","C"] <- 1; P["C","A"] <- 1
  cs <- chain_structure(P)
  expect_true(cs$is_irreducible)
  expect_false(cs$is_aperiodic)
  expect_false(cs$is_regular)
  expect_equal(unname(cs$period), c(3L, 3L, 3L))
  expect_length(cs$communicating_classes, 1L)
})

test_that("chain_structure on an absorbing chain identifies the absorbing state", {
  P <- matrix(c(0.5, 0.3, 0.2,
                0.4, 0.4, 0.2,
                0.0, 0.0, 1.0),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("X","Y","Done"), c("X","Y","Done")))
  cs <- chain_structure(P)
  expect_equal(cs$absorbing_states, "Done")
  expect_equal(unname(cs$classification),
               c("transient", "transient", "absorbing"))
  expect_false(cs$is_irreducible)
  expect_false(cs$is_regular)
  expect_true(is.na(cs$is_reversible))
  # Both transient states absorb into Done with probability 1.
  expect_equal(unname(cs$absorption_probabilities[, "Done"]), c(1, 1))
  # Mean absorption time positive and finite.
  expect_true(all(cs$mean_absorption_time > 0 & is.finite(cs$mean_absorption_time)))
})

test_that("chain_structure on two disjoint closed classes has zero cross-hitting", {
  P <- matrix(0, 4, 4,
              dimnames = list(c("a1","a2","b1","b2"), c("a1","a2","b1","b2")))
  P["a1","a2"] <- 1; P["a2","a1"] <- 1
  P["b1","b2"] <- 1; P["b2","b1"] <- 1
  cs <- chain_structure(P)
  expect_false(cs$is_irreducible)
  expect_length(cs$communicating_classes, 2L)
  # Cross-class hitting probability is 0
  expect_equal(cs$hitting_probabilities["a1", "b1"], 0)
  expect_equal(cs$hitting_probabilities["b1", "a1"], 0)
  # Within-class hitting probability is 1
  expect_equal(cs$hitting_probabilities["a1", "a2"], 1)
})

test_that("chain_structure accepts a netobject and a sequence data.frame", {
  net <- build_network(as.data.frame(trajectories), method = "relative")
  cs_net <- chain_structure(net)
  expect_s3_class(cs_net, "chain_structure")
  expect_true(cs_net$is_irreducible)

  cs_df <- chain_structure(as.data.frame(trajectories))
  expect_s3_class(cs_df, "chain_structure")
  expect_equal(cs_net$states, cs_df$states)
  expect_equal(cs_net$P, cs_df$P, tolerance = 1e-12)
})

test_that("print and summary dispatch without error", {
  net <- build_network(as.data.frame(trajectories), method = "relative")
  cs <- chain_structure(net)
  expect_output(print(cs), "Chain structure")
  s <- summary(cs)
  expect_s3_class(s, "data.frame")
  expect_true("classification" %in% colnames(s))
})

test_that("chain_structure dispatches over netobject_group", {
  seqs <- data.frame(
    V1 = rep(c("A","B","C"), 6L),
    V2 = rep(c("B","C","A"), 6L),
    V3 = rep(c("C","A","B"), 6L),
    grp = rep(c("g1", "g2"), each = 9L),
    stringsAsFactors = FALSE
  )
  grp_net <- build_network(seqs, method = "relative", group = "grp")
  cs_grp <- chain_structure(grp_net)
  expect_s3_class(cs_grp, "chain_structure_group")
  expect_named(cs_grp, names(grp_net))
  expect_true(all(vapply(cs_grp, inherits, logical(1), "chain_structure")))

  s <- summary(cs_grp)
  expect_s3_class(s, "data.frame")
  expect_true(all(c("group", "state", "classification") %in% colnames(s)))
  expect_equal(nrow(s),
               sum(vapply(cs_grp, function(x) length(x$states), integer(1))))
  expect_invisible(print(cs_grp))
})

test_that("plot returns a ggplot for regular, absorbing, and multi-class chains", {
  skip_if_not_installed("ggplot2")
  # Regular
  net <- build_network(as.data.frame(trajectories), method = "relative")
  expect_s3_class(suppressWarnings(plot(chain_structure(net))), "ggplot")
  # Absorbing
  P <- matrix(c(0.5, 0.3, 0.2, 0.4, 0.4, 0.2, 0, 0, 1), 3, 3, byrow = TRUE,
              dimnames = list(c("X","Y","Done"), c("X","Y","Done")))
  expect_s3_class(suppressWarnings(plot(chain_structure(P))), "ggplot")
  # Multi-class
  P2 <- matrix(0, 4, 4,
               dimnames = list(c("a1","a2","b1","b2"), c("a1","a2","b1","b2")))
  P2["a1","a2"] <- 1; P2["a2","a1"] <- 1
  P2["b1","b2"] <- 1; P2["b2","b1"] <- 1
  expect_s3_class(suppressWarnings(plot(chain_structure(P2))), "ggplot")
})
