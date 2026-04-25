# Helper: simulate an order-k Markov sequence
.sim_order_k <- function(k, n_seqs = 20L, len = 80L, states = letters[1:4],
                          seed = 1L) {
  set.seed(seed)
  S <- length(states)
  if (k == 1L) {
    tm <- matrix(runif(S * S), S, S, dimnames = list(states, states))
    tm <- tm / rowSums(tm)
    lapply(seq_len(n_seqs), function(.) {
      s <- character(len); s[1L] <- sample(states, 1L)
      for (i in 2L:len) s[i] <- sample(states, 1L, prob = tm[s[i - 1L], ])
      s
    })
  } else if (k == 2L) {
    keys <- as.vector(outer(states, states, paste, sep = "|"))
    tm <- matrix(runif(length(keys) * S), length(keys), S,
                 dimnames = list(keys, states))
    tm <- tm / rowSums(tm)
    lapply(seq_len(n_seqs), function(.) {
      s <- character(len); s[1:2] <- sample(states, 2L, replace = TRUE)
      for (i in 3L:len) {
        key <- paste(s[i - 2L], s[i - 1L], sep = "|")
        s[i] <- sample(states, 1L, prob = tm[key, ])
      }
      s
    })
  } else stop("only k = 1 or 2 supported in helper")
}

test_that("markov_order_test returns proper structure", {
  seqs <- .sim_order_k(1L, n_seqs = 10L, len = 30L, seed = 42L)
  res <- markov_order_test(seqs, max_order = 2L, n_perm = 50L, seed = 1L)

  expect_s3_class(res, "net_markov_order")
  expect_true(all(c("optimal_order", "test_table", "permutation_null",
                    "logliks", "transition_matrices", "states",
                    "n_perm", "alpha", "max_order") %in% names(res)))
  expect_s3_class(res$test_table, "data.frame")
  expect_true(all(c("order", "loglik", "df", "g2",
                    "p_permutation", "p_asymptotic", "significant") %in%
                  names(res$test_table)))
  expect_equal(nrow(res$test_table), 3L)
  expect_length(res$permutation_null, 2L)
})

test_that("summary returns tidy df with selected-order attribute", {
  seqs <- .sim_order_k(1L, n_seqs = 8L, len = 25L, seed = 7L)
  res <- markov_order_test(seqs, max_order = 2L, n_perm = 50L, seed = 1L)
  s <- summary(res)
  expect_s3_class(s, "data.frame")
  expect_equal(attr(s, "optimal_order"), res$optimal_order)
})

test_that("first-order data selects order 1", {
  seqs <- .sim_order_k(1L, n_seqs = 30L, len = 80L, seed = 123L)
  res <- markov_order_test(seqs, max_order = 3L, n_perm = 200L, seed = 1L)
  expect_equal(res$optimal_order, 1L)
  expect_lt(res$test_table$p_permutation[2L], 0.05)
  expect_gt(res$test_table$p_permutation[3L], 0.05)
})

test_that("second-order data selects order >= 2", {
  seqs <- .sim_order_k(2L, n_seqs = 40L, len = 100L, seed = 321L)
  res <- markov_order_test(seqs, max_order = 3L, n_perm = 200L, seed = 1L)
  expect_gte(res$optimal_order, 2L)
})

test_that("plot returns a ggplot (single panel)", {
  seqs <- .sim_order_k(1L, n_seqs = 8L, len = 25L, seed = 7L)
  res <- markov_order_test(seqs, max_order = 2L, n_perm = 50L, seed = 1L)
  g <- plot(res, panel = "ic")
  expect_s3_class(g, "ggplot")
  g2 <- plot(res, panel = "permutation")
  expect_s3_class(g2, "ggplot")
})

test_that("print runs without error", {
  seqs <- .sim_order_k(1L, n_seqs = 6L, len = 20L, seed = 7L)
  res <- markov_order_test(seqs, max_order = 2L, n_perm = 30L, seed = 1L)
  expect_output(print(res), "Markov Order Test")
  expect_output(print(res), "Selected order")
})
