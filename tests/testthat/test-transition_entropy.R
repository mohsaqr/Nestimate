# ---- Tests for transition_entropy() ----

.P3 <- matrix(
  c(0.7, 0.2, 0.1,
    0.3, 0.5, 0.2,
    0.2, 0.3, 0.5),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
)

# ---- structure ----

test_that("transition_entropy returns net_transition_entropy with all fields", {
  te <- transition_entropy(.P3, base = 2)
  expect_s3_class(te, "net_transition_entropy")
  expect_named(te, c("row_entropy", "row_entropy_norm",
                     "stationary",
                     "stationary_entropy", "stationary_entropy_norm",
                     "entropy_rate", "entropy_rate_norm",
                     "redundancy", "redundancy_norm",
                     "max_entropy", "base", "states"))
  expect_length(te$row_entropy, 3L)
  expect_length(te$stationary,  3L)
  expect_equal(te$states, c("A", "B", "C"))
  expect_equal(te$base, 2)
})

test_that("normalised fields equal raw / log_b(n) and lie in [0, 1]", {
  te <- transition_entropy(.P3, base = 2)
  H_max <- log2(3)
  expect_equal(te$max_entropy,             H_max,                       tolerance = 1e-12)
  expect_equal(te$entropy_rate_norm,       te$entropy_rate       / H_max, tolerance = 1e-12)
  expect_equal(te$stationary_entropy_norm, te$stationary_entropy / H_max, tolerance = 1e-12)
  expect_equal(unname(te$row_entropy_norm),
               unname(te$row_entropy) / H_max, tolerance = 1e-12)
  expect_true(te$entropy_rate_norm       >= 0 && te$entropy_rate_norm       <= 1 + 1e-10)
  expect_true(te$stationary_entropy_norm >= 0 && te$stationary_entropy_norm <= 1 + 1e-10)
  expect_true(all(te$row_entropy_norm >= 0 - 1e-10) &&
              all(te$row_entropy_norm <= 1 + 1e-10))
})

test_that("redundancy_norm is fraction of marginal uncertainty captured by memory", {
  te <- transition_entropy(.P3, base = 2)
  expect_equal(te$redundancy_norm,
               te$redundancy / te$stationary_entropy,
               tolerance = 1e-12)
  expect_true(te$redundancy_norm >= 0 - 1e-10 &&
              te$redundancy_norm <= 1 + 1e-10)
})

# ---- numerical identities ----

test_that("entropy_rate equals sum(pi * row_entropy)", {
  te <- transition_entropy(.P3, base = 2)
  expect_equal(te$entropy_rate,
               sum(te$stationary * te$row_entropy),
               tolerance = 1e-12)
})

test_that("redundancy equals stationary_entropy minus entropy_rate", {
  te <- transition_entropy(.P3, base = 2)
  expect_equal(te$redundancy,
               te$stationary_entropy - te$entropy_rate,
               tolerance = 1e-12)
})

test_that("entropy_rate <= stationary_entropy (information inequality)", {
  te <- transition_entropy(.P3, base = 2)
  expect_lte(te$entropy_rate, te$stationary_entropy + 1e-10)
  expect_gte(te$redundancy, -1e-10)
})

test_that("uniform rows give zero redundancy and h = log_b(n)", {
  P <- matrix(1 / 3, 3, 3,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  te <- transition_entropy(P, base = 2)
  expect_equal(te$redundancy, 0, tolerance = 1e-10)
  expect_equal(te$entropy_rate,       log2(3), tolerance = 1e-10)
  expect_equal(te$stationary_entropy, log2(3), tolerance = 1e-10)
})

test_that("deterministic permutation chain has zero entropy rate", {
  P <- matrix(c(0, 1, 0,
                0, 0, 1,
                1, 0, 0),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  te <- transition_entropy(P, base = 2)
  expect_equal(te$entropy_rate, 0, tolerance = 1e-12)
  expect_equal(te$row_entropy, c(A = 0, B = 0, C = 0), tolerance = 1e-12)
})

test_that("base 2 vs base e differ by ln(2) factor", {
  te2 <- transition_entropy(.P3, base = 2)
  tee <- transition_entropy(.P3, base = exp(1))
  expect_equal(te2$entropy_rate * log(2), tee$entropy_rate,
               tolerance = 1e-12)
  expect_equal(te2$row_entropy * log(2), tee$row_entropy,
               tolerance = 1e-12)
})

test_that("row entropies bounded by log_b(n_nonzero_columns)", {
  te <- transition_entropy(.P3, base = 2)
  expect_true(all(te$row_entropy <= log2(3) + 1e-10))
  expect_true(all(te$row_entropy >= 0 - 1e-10))
})

# ---- stationary distribution checks ----

test_that("stationary distribution sums to 1 and satisfies pi P = pi", {
  te <- transition_entropy(.P3, base = 2)
  expect_equal(sum(te$stationary), 1, tolerance = 1e-10)
  expect_equal(as.vector(te$stationary %*% .P3),
               unname(te$stationary), tolerance = 1e-8)
})

# ---- input dispatch ----

test_that("transition_entropy accepts a netobject", {
  data(trajectories, package = "Nestimate")
  net <- build_network(as.data.frame(trajectories), method = "relative")
  te  <- transition_entropy(net, base = 2)
  expect_s3_class(te, "net_transition_entropy")
  expect_length(te$states, ncol(net$weights))
})

test_that("transition_entropy accepts a wide data.frame", {
  data(trajectories, package = "Nestimate")
  te <- transition_entropy(as.data.frame(trajectories), base = 2)
  expect_s3_class(te, "net_transition_entropy")
})

test_that("transition_entropy normalises rows that don't sum to 1", {
  P <- .P3 * 2
  expect_warning(te <- transition_entropy(P, base = 2),
                 regexp = "normalizing")
  te_ref <- transition_entropy(.P3, base = 2)
  expect_equal(te$entropy_rate, te_ref$entropy_rate, tolerance = 1e-10)
})

test_that("base argument is validated", {
  expect_error(transition_entropy(.P3, base = 1),
               regexp = "base != 1|stopifnot")
  expect_error(transition_entropy(.P3, base = -2))
  expect_error(transition_entropy(.P3, base = c(2, 10)))
})

# ---- group dispatch ----

test_that("transition_entropy dispatches on netobject_group", {
  data(trajectories, package = "Nestimate")
  df <- as.data.frame(trajectories)
  df$grp <- rep(c("g1", "g2"), length.out = nrow(df))
  grp <- build_network(df, method = "relative", group = "grp")
  te_grp <- transition_entropy(grp, base = 2)
  expect_s3_class(te_grp, "net_transition_entropy_group")
  expect_length(te_grp, 2L)
  expect_true(all(vapply(te_grp, inherits, logical(1),
                         "net_transition_entropy")))
})

# ---- S3 methods dispatch ----

test_that("print/summary/plot dispatch correctly", {
  te <- transition_entropy(.P3, base = 2)
  expect_output(print(te), "Transition Entropy")
  s <- summary(te)
  expect_s3_class(s, "summary.net_transition_entropy")
  expect_true(is.data.frame(s$table))
  expect_true(is.data.frame(s$chain))
  ## Tidy table has 6 columns including normalised + percent
  expect_named(s$table, c("state", "stationary", "row_entropy",
                          "row_entropy_norm",
                          "contribution", "contribution_pct"))
  ## Sorted by contribution_pct descending
  expect_equal(s$table$contribution_pct,
               sort(s$table$contribution_pct, decreasing = TRUE))
  ## Contributions sum to entropy_rate; pct sums to 100
  expect_equal(sum(s$table$contribution), te$entropy_rate, tolerance = 1e-10)
  expect_equal(sum(s$table$contribution_pct), 100, tolerance = 1e-10)
  ## Chain table has 4 rows: h(P), H(pi), redundancy, ceiling
  expect_equal(nrow(s$chain), 4L)
  g <- plot(te)
  expect_true(inherits(g, "gg"))
})
