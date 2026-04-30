# ---- Tests for path_dependence() ----

# Synthetic example: 4 sequences, 5 columns
.sim_seqs <- function() {
  data.frame(
    t1 = c("A", "A", "B", "C"),
    t2 = c("B", "B", "A", "B"),
    t3 = c("A", "C", "B", "A"),
    t4 = c("B", "A", "C", "B"),
    t5 = c("C", "B", "A", "A"),
    stringsAsFactors = FALSE
  )
}

# ---- structure ----

test_that("path_dependence returns net_path_dependence with all fields", {
  d  <- .sim_seqs()
  pd <- path_dependence(d, order = 2L, min_count = 1L)
  expect_s3_class(pd, "net_path_dependence")
  expect_named(pd, c("contexts", "chain", "order", "base", "min_count",
                     "states"))
  expect_equal(pd$order, 2L)
  expect_equal(pd$base, 2)
  expect_equal(pd$min_count, 1L)
  expect_true(is.data.frame(pd$contexts))
  expect_named(pd$contexts, c("context", "n", "H_order1", "H_orderk",
                              "H_drop", "KL", "top_o1", "top_ok", "flips"))
})

test_that("contexts are sorted by KL descending", {
  data(trajectories, package = "Nestimate")
  pd <- path_dependence(as.data.frame(trajectories), order = 2L,
                        min_count = 5L)
  expect_equal(pd$contexts$KL,
               sort(pd$contexts$KL, decreasing = TRUE))
})

# ---- numerical identities ----

test_that("KL = 0 when order-k matches order-1 (i.i.d. data)", {
  ## Construct sequences where the next state depends ONLY on current,
  ## not on the prior — order-1 and order-2 should agree exactly.
  set.seed(42)
  states <- c("A", "B", "C")
  P <- matrix(c(0.6, 0.2, 0.2,
                0.1, 0.7, 0.2,
                0.3, 0.3, 0.4),
              nrow = 3, byrow = TRUE,
              dimnames = list(states, states))
  n_seq <- 200; n_t <- 50
  m <- matrix(NA_character_, n_seq, n_t)
  m[, 1L] <- sample(states, n_seq, replace = TRUE)
  for (t in 2:n_t) {
    for (i in seq_len(n_seq)) {
      m[i, t] <- sample(states, 1L, prob = P[m[i, t - 1L], ])
    }
  }
  pd <- path_dependence(as.data.frame(m, stringsAsFactors = FALSE),
                         order = 2L, min_count = 30L)
  ## Random sampling means small sampling-error KL, but should be << 0.05 bits
  expect_lt(pd$chain$KL_weighted, 0.05)
})

test_that("H_drop equals H_order1 - H_orderk per row", {
  data(trajectories, package = "Nestimate")
  pd <- path_dependence(as.data.frame(trajectories), order = 2L,
                        min_count = 5L)
  expect_equal(pd$contexts$H_drop,
               pd$contexts$H_order1 - pd$contexts$H_orderk,
               tolerance = 1e-12)
})

test_that("KL is non-negative for finite values", {
  data(trajectories, package = "Nestimate")
  pd <- path_dependence(as.data.frame(trajectories), order = 2L,
                        min_count = 5L)
  expect_true(all(pd$contexts$KL[is.finite(pd$contexts$KL)] >= 0 - 1e-10))
})

test_that("base 2 vs base e differ by ln(2)", {
  data(trajectories, package = "Nestimate")
  pd2 <- path_dependence(as.data.frame(trajectories), order = 2L,
                         min_count = 5L, base = 2)
  pde <- path_dependence(as.data.frame(trajectories), order = 2L,
                         min_count = 5L, base = exp(1))
  expect_equal(pd2$contexts$KL * log(2),
               pde$contexts$KL,
               tolerance = 1e-12)
  expect_equal(pd2$chain$KL_weighted * log(2),
               pde$chain$KL_weighted,
               tolerance = 1e-12)
})

# ---- input validation ----

test_that("path_dependence validates order and min_count", {
  d <- .sim_seqs()
  expect_error(path_dependence(d, order = 1L), "order' must be >= 2")
  expect_error(path_dependence(d, min_count = 0L), "min_count' must be >= 1")
  expect_error(path_dependence(d, base = 1), "base != 1|stopifnot")
})

test_that("order higher than sequence length errors informatively", {
  d <- .sim_seqs()
  expect_error(path_dependence(d, order = 10L),
               regexp = "Sequences too short")
})

# ---- order-3 works ----

test_that("path_dependence works at order 3", {
  data(trajectories, package = "Nestimate")
  pd3 <- path_dependence(as.data.frame(trajectories), order = 3L,
                         min_count = 5L)
  expect_equal(pd3$order, 3L)
  expect_true(nrow(pd3$contexts) > 0L)
  ## Order-3 contexts have format "X -> Y -> Z" (two arrows)
  arrows <- vapply(pd3$contexts$context,
                   function(s) length(gregexpr(" -> ", s, fixed = TRUE)[[1]]),
                   integer(1))
  expect_true(all(arrows == 2L))
})

# ---- S3 dispatch ----

test_that("print/summary/plot dispatch correctly", {
  data(trajectories, package = "Nestimate")
  pd <- path_dependence(as.data.frame(trajectories), order = 2L,
                        min_count = 5L)
  expect_output(print(pd), "Path Dependence")
  s <- summary(pd)
  expect_s3_class(s, "summary.net_path_dependence")
  expect_output(print(s), "Path Dependence Summary")
  g <- plot(pd)
  expect_true(inherits(g, "gg"))
})
