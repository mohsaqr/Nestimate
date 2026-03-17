# ---- wtna() tests ----

test_that("wtna transition counts match manual crossprod", {
  # 3-code, 5-row one-hot data
  df <- data.frame(
    A = c(1, 0, 1, 0, 1),
    B = c(0, 1, 0, 1, 0),
    C = c(0, 0, 0, 0, 0)
  )

  net <- wtna(df, method = "transition", type = "frequency")
  expect_s3_class(net, "netobject")
  expect_true(net$directed)

  # Manual: crossprod(X[-5,], X[-1,])
  X <- as.matrix(df)
  expected <- crossprod(X[-5, ], X[-1, ])
  expect_equal(unname(net$matrix), unname(expected))
})

test_that("wtna cooccurrence counts match crossprod", {
  df <- data.frame(
    A = c(1, 0, 1, 0, 1),
    B = c(0, 1, 0, 1, 0),
    C = c(0, 0, 1, 0, 0)
  )

  net <- wtna(df, method = "cooccurrence", type = "frequency")
  expect_s3_class(net, "netobject")
  expect_false(net$directed)

  # Manual: crossprod(X)
  X <- as.matrix(df)
  expected <- crossprod(X)
  expect_equal(unname(net$matrix), unname(expected))
})

test_that("wtna method='both' returns list of two networks", {
  df <- data.frame(
    A = c(1, 0, 1, 0),
    B = c(0, 1, 0, 1)
  )

  result <- wtna(df, method = "both")
  expect_true(is.list(result))
  expect_true("transition" %in% names(result))
  expect_true("cooccurrence" %in% names(result))
  expect_s3_class(result$transition, "netobject")
  expect_s3_class(result$cooccurrence, "netobject")
  expect_true(result$transition$directed)
  expect_false(result$cooccurrence$directed)
})

test_that("wtna non-overlapping window aggregation", {
  df <- data.frame(
    A = c(1, 0, 0, 1),
    B = c(0, 1, 1, 0)
  )

  net <- wtna(df, method = "transition", type = "frequency",
              window_size = 2, mode = "non-overlapping")
  expect_s3_class(net, "netobject")

  # Non-overlapping: rows 1-2 -> window 0, rows 3-4 -> window 1
  # Window 0: A=1,B=1 (both active); Window 1: A=1,B=1
  # Transitions: (1,1) -> (1,1) = crossprod should be 1 on all pairs
  X <- as.matrix(df)
  wid <- c(0L, 0L, 1L, 1L)
  agg <- rowsum(X, wid)
  agg[agg > 0] <- 1
  expected <- crossprod(agg[-nrow(agg), , drop = FALSE], agg[-1, , drop = FALSE])
  expect_equal(unname(net$matrix), unname(expected))
})

test_that("wtna overlapping window aggregation", {
  df <- data.frame(
    A = c(1, 0, 0, 1),
    B = c(0, 1, 1, 0)
  )

  net <- wtna(df, method = "transition", type = "frequency",
              window_size = 2, mode = "overlapping")
  expect_s3_class(net, "netobject")

  # Overlapping with ws=2: windows at rows (1,2), (2,3), (3,4)
  X <- as.matrix(df)
  n <- nrow(X)
  combined <- X[1:3, ] | X[2:4, ]
  storage.mode(combined) <- "integer"
  expected <- crossprod(combined[-3, ], combined[-1, ])
  expect_equal(unname(net$matrix), unname(expected))
})

test_that("wtna per-actor grouping", {
  df <- data.frame(
    actor = c(1, 1, 1, 2, 2, 2),
    A = c(1, 0, 1, 0, 1, 0),
    B = c(0, 1, 0, 1, 0, 1)
  )

  net <- wtna(df, method = "transition", type = "frequency",
              codes = c("A", "B"), actor = "actor")
  expect_s3_class(net, "netobject")

  # Manual: compute per actor and sum
  g1 <- as.matrix(df[df$actor == 1, c("A", "B")])
  g2 <- as.matrix(df[df$actor == 2, c("A", "B")])
  t1 <- crossprod(g1[-3, ], g1[-1, ])
  t2 <- crossprod(g2[-3, ], g2[-1, ])
  expected <- t1 + t2
  expect_equal(unname(net$matrix), unname(expected))
})

test_that("wtna type='relative' row-normalizes", {
  df <- data.frame(
    A = c(1, 0, 1, 0, 1),
    B = c(0, 1, 0, 1, 0)
  )

  net <- wtna(df, method = "transition", type = "relative")
  rs <- rowSums(net$matrix)
  # Non-zero rows should sum to 1
  non_zero <- rs > 0
  if (any(non_zero)) {
    expect_true(all(abs(rs[non_zero] - 1) < 1e-10))
  }
})

test_that("wtna auto-detects one-hot columns", {
  df <- data.frame(
    id = c(1, 2, 3, 4),
    A = c(1, 0, 1, 0),
    B = c(0, 1, 0, 1),
    score = c(3.5, 2.1, 4.0, 1.5)
  )

  # id and score are not binary, should be excluded
  # Actually id column has values in 0/1 range... let's use different ids
  df2 <- data.frame(
    name = c("x", "y", "z", "w"),
    A = c(1, 0, 1, 0),
    B = c(0, 1, 0, 1),
    score = c(3.5, 2.1, 4.0, 1.5)
  )

  net <- wtna(df2, method = "transition")
  # Should only use A and B (binary columns)
  expect_equal(sort(net$nodes), c("A", "B"))
})

test_that("wtna single row edge case", {
  df <- data.frame(A = 1, B = 0)
  net <- wtna(df, method = "transition", type = "frequency")
  expect_s3_class(net, "netobject")
  expect_true(all(net$matrix == 0))
})

test_that("wtna all-zeros edge case", {
  df <- data.frame(
    A = c(0, 0, 0),
    B = c(0, 0, 0)
  )
  net <- wtna(df, method = "transition", type = "frequency")
  expect_s3_class(net, "netobject")
  expect_true(all(net$matrix == 0))
})

test_that("wtna single actor with actor param", {
  df <- data.frame(
    actor = c(1, 1, 1),
    A = c(1, 0, 1),
    B = c(0, 1, 0)
  )
  net <- wtna(df, method = "transition", codes = c("A", "B"), actor = "actor")
  expect_s3_class(net, "netobject")
})

test_that("wtna cooccurrence is symmetric", {
  df <- data.frame(
    A = c(1, 0, 1, 1),
    B = c(0, 1, 1, 0),
    C = c(1, 0, 0, 1)
  )

  net <- wtna(df, method = "cooccurrence", type = "frequency")
  expect_true(isSymmetric(unname(net$matrix)))
})

test_that("wtna validates inputs", {
  df <- data.frame(A = c(1, 0), B = c(0, 1))
  expect_error(wtna(df, codes = c("X", "Y")))
  expect_error(wtna(df, codes = c("A"), method = "transition"))
})

# ---- Registry integration tests ----

test_that("wtna works via build_network registry", {
  df <- data.frame(
    A = c(1, 0, 1, 0, 1),
    B = c(0, 1, 0, 1, 0),
    C = c(0, 0, 1, 0, 0)
  )

  net <- build_network(df, method = "wtna")
  expect_s3_class(net, "netobject")
  expect_true(net$directed)
  expect_equal(net$method, "wtna")

  net2 <- build_network(df, method = "cna")
  expect_s3_class(net2, "netobject")
  expect_false(net2$directed)
  expect_equal(net2$method, "co_occurrence")
})

test_that("wtna registry with window_size param", {
  df <- data.frame(
    A = c(1, 0, 0, 1, 1, 0),
    B = c(0, 1, 1, 0, 0, 1)
  )

  net <- build_network(df, method = "wtna",
                       params = list(window_size = 2,
                                     mode = "non-overlapping"))
  expect_s3_class(net, "netobject")
})

test_that("wtna registry with actor param", {
  df <- data.frame(
    actor = c(1, 1, 1, 2, 2, 2),
    A = c(1, 0, 1, 0, 1, 0),
    B = c(0, 1, 0, 1, 0, 1)
  )

  net <- build_network(df, method = "wtna",
                       params = list(codes = c("A", "B"), actor = "actor"))
  expect_s3_class(net, "netobject")
})

test_that("wtna registered, wcna is alias for co_occurrence", {
  estimators <- list_estimators()
  expect_true("wtna" %in% estimators$name)
  expect_true(estimators$directed[estimators$name == "wtna"])

  # wcna resolves to co_occurrence
  expect_equal(.resolve_method_alias("wcna"), "co_occurrence")
  expect_equal(.resolve_method_alias("cna"), "co_occurrence")
})

test_that("wtna bootstrap works via wtna()", {
  set.seed(42)
  df <- data.frame(
    A = sample(0:1, 30, replace = TRUE),
    B = sample(0:1, 30, replace = TRUE),
    C = sample(0:1, 30, replace = TRUE)
  )

  net <- wtna(df, method = "transition")
  boot <- bootstrap_network(net, iter = 20, seed = 1)
  expect_s3_class(boot, "net_bootstrap")
  expect_true(boot$original$directed)
  expect_true(nrow(boot$summary) > 0)
})

test_that("wtna bootstrap via build_network falls back to association path", {
  set.seed(42)
  df <- data.frame(
    A = sample(0:1, 30, replace = TRUE),
    B = sample(0:1, 30, replace = TRUE),
    C = sample(0:1, 30, replace = TRUE)
  )

  # build_network(method="wtna") stores method as "wtna" (not co_occurrence),

  # so bootstrap dispatches via association path — no error but prefer wtna()
  boot <- bootstrap_network(build_network(df, method = "wtna"),
                            iter = 5, seed = 1)
  expect_s3_class(boot, "net_bootstrap")
})

test_that("wcna bootstrap works via wtna()", {
  set.seed(42)
  df <- data.frame(
    A = sample(0:1, 30, replace = TRUE),
    B = sample(0:1, 30, replace = TRUE),
    C = sample(0:1, 30, replace = TRUE)
  )

  net <- wtna(df, method = "cooccurrence")
  boot <- bootstrap_network(net, iter = 20, seed = 1)
  expect_s3_class(boot, "net_bootstrap")
  expect_equal(boot$method, "co_occurrence")
  expect_false(boot$original$directed)
})

test_that("wcna bootstrap errors with build_network (no stored data)", {
  set.seed(42)
  df <- data.frame(
    A = sample(0:1, 30, replace = TRUE),
    B = sample(0:1, 30, replace = TRUE),
    C = sample(0:1, 30, replace = TRUE)
  )

  expect_error(
    bootstrap_network(build_network(df, method = "cna"), iter = 20, seed = 1),
    "requires the original data"
  )
})

test_that("wtna permutation test works via wtna()", {
  set.seed(42)
  df1 <- data.frame(
    A = sample(0:1, 20, replace = TRUE),
    B = sample(0:1, 20, replace = TRUE)
  )
  df2 <- data.frame(
    A = sample(0:1, 20, replace = TRUE),
    B = sample(0:1, 20, replace = TRUE)
  )

  net1 <- wtna(df1, method = "transition")
  net2 <- wtna(df2, method = "transition")
  perm <- permutation_test(net1, net2, iter = 20, seed = 1)

  expect_s3_class(perm, "net_permutation")
  expect_true(is.matrix(perm$p_values))
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))
})

test_that("wtna permutation via build_network falls back to association path", {
  set.seed(42)
  df1 <- data.frame(A = sample(0:1, 20, replace = TRUE), B = sample(0:1, 20, replace = TRUE))
  df2 <- data.frame(A = sample(0:1, 20, replace = TRUE), B = sample(0:1, 20, replace = TRUE))

  # build_network(method="wtna") stores method as "wtna" (not co_occurrence),
  # so permutation dispatches via association path
  perm <- permutation_test(build_network(df1, method = "wtna"),
                           build_network(df2, method = "wtna"),
                           iter = 5, seed = 1)
  expect_s3_class(perm, "net_permutation")
})

test_that("wcna permutation test works via wtna()", {
  set.seed(42)
  df1 <- data.frame(
    A = sample(0:1, 20, replace = TRUE),
    B = sample(0:1, 20, replace = TRUE)
  )
  df2 <- data.frame(
    A = sample(0:1, 20, replace = TRUE),
    B = sample(0:1, 20, replace = TRUE)
  )

  net1 <- wtna(df1, method = "cooccurrence")
  net2 <- wtna(df2, method = "cooccurrence")
  perm <- permutation_test(net1, net2, iter = 20, seed = 1)

  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "co_occurrence")
})

test_that("wcna permutation errors with build_network (no stored data)", {
  set.seed(42)
  df1 <- data.frame(
    A = sample(0:1, 20, replace = TRUE),
    B = sample(0:1, 20, replace = TRUE)
  )
  df2 <- data.frame(
    A = sample(0:1, 20, replace = TRUE),
    B = sample(0:1, 20, replace = TRUE)
  )

  expect_error(
    permutation_test(build_network(df1, method = "cna"),
                     build_network(df2, method = "cna"),
                     iter = 20, seed = 1),
    "requires the original data"
  )
})
