testthat::skip_on_cran()

# ---- Tests for centrality_stability() ----

test_that("centrality_stability works for cor (association) method", {
  set.seed(5)
  df <- as.data.frame(matrix(rpois(100 * 5, 10), nrow = 100))
  colnames(df) <- paste0("V", 1:5)
  net <- build_network(df, method = "cor")
  cs <- centrality_stability(net, iter = 20L,
                             measures = c("InStrength", "OutStrength"),
                             drop_prop = c(0.2, 0.4),
                             seed = 5)

  expect_s3_class(cs, "net_stability")
  expect_true(all(cs$cs %in% c(0, cs$drop_prop)))
})


# ---- association build_matrix function (L181-197) ----

test_that("centrality_stability association path tolerates estimator errors gracefully", {
  set.seed(9)
  df <- as.data.frame(matrix(rpois(60 * 4, 10), nrow = 60))
  colnames(df) <- paste0("V", 1:4)
  net <- build_network(df, method = "pcor")
  cs <- centrality_stability(net, iter = 20L,
                             measures = c("InStrength"),
                             drop_prop = c(0.3, 0.6),
                             seed = 9)

  expect_s3_class(cs, "net_stability")
  # CS should be 0 or a valid drop_prop value
  expect_true(cs$cs["InStrength"] %in% c(0, cs$drop_prop))
})


# ---- single-measure storage path (L216 and L222) ----

test_that(".calculate_cs returns 0 when no prop_above meets certainty", {
  # Build a correlation matrix where certainty is never met
  iter <- 10L
  n_prop <- 3L
  corr_mat <- matrix(0, nrow = iter, ncol = n_prop)  # all zeros < threshold
  result <- Nestimate:::.calculate_cs(corr_mat, threshold = 0.7, certainty = 0.95,
                                      drop_prop = c(0.1, 0.3, 0.5))
  expect_equal(result, 0)
})

test_that(".calculate_cs returns max valid drop_prop when certainty is met", {
  iter <- 20L
  n_prop <- 3L
  corr_mat <- matrix(1, nrow = iter, ncol = n_prop)  # all ones >= threshold
  result <- Nestimate:::.calculate_cs(corr_mat, threshold = 0.7, certainty = 0.95,
                                      drop_prop = c(0.1, 0.3, 0.5))
  expect_equal(result, 0.5)
})


# ---- cograph_network input (L76) ----

test_that("centrality.netobject returns correct directed defaults (L163-177)", {
  seqs <- data.frame(
    V1 = c("A","B","A","C","B","A"),
    V2 = c("B","C","B","A","C","B"),
    V3 = c("C","A","C","B","A","C")
  )
  net <- build_network(seqs, method = "relative")
  c1 <- net_centrality(net)
  expect_true(is.data.frame(c1))
  expect_equal(nrow(c1), 3)
  expect_true(all(c("InStrength", "OutStrength", "Betweenness") %in% names(c1)))
})

test_that("centrality.netobject returns correct undirected defaults (L163-177)", {
  set.seed(42)
  panel <- data.frame(V1 = rnorm(50), V2 = rnorm(50), V3 = rnorm(50))
  net_ud <- build_network(panel, method = "cor")
  c2 <- net_centrality(net_ud)
  expect_true(is.data.frame(c2))
  expect_true(all(c("Closeness", "Betweenness") %in% names(c2)))
})

test_that("centrality.netobject_group returns list of data frames (L185-188)", {
  seqs <- data.frame(
    V1 = c("A","B","A","C","B","A"),
    V2 = c("B","C","B","A","C","B"),
    V3 = c("C","A","C","B","A","C"),
    grp = c("X","X","X","Y","Y","Y")
  )
  nets <- build_network(seqs, method = "relative", group = "grp")
  c3 <- net_centrality(nets)
  expect_true(is.list(c3))
  expect_equal(length(c3), 2)
  expect_true(all(vapply(c3, is.data.frame, logical(1))))
})

test_that(".betweenness returns zeros for n < 3 (L53)", {
  W <- matrix(c(0, 1, 1, 0), nrow = 2, dimnames = list(c("A","B"), c("A","B")))
  btw <- Nestimate:::.betweenness(W, directed = TRUE)
  expect_equal(unname(btw), c(0, 0))
  expect_equal(names(btw), c("A", "B"))
})

test_that(".compute_centralities handles external centrality_fn (L325-340)", {
  seqs <- data.frame(
    V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
    V3 = c("C","A","C","B")
  )
  net <- build_network(seqs, method = "relative")
  custom_fn <- function(mat) {
    list(MyMeasure = setNames(rowSums(abs(mat)), rownames(mat)))
  }
  c4 <- net_centrality(net, measures = c("InStrength", "MyMeasure"),
                    centrality_fn = custom_fn)
  expect_true("MyMeasure" %in% names(c4))
  expect_true(is.data.frame(c4))
})

test_that(".compute_centralities errors when external measure lacks centrality_fn (L325-329)", {
  seqs <- data.frame(
    V1 = c("A","B","A"), V2 = c("B","C","B"), V3 = c("C","A","C")
  )
  net <- build_network(seqs, method = "relative")
  expect_error(
    net_centrality(net, measures = c("InStrength", "BadMeasure")),
    "centrality_fn is required"
  )
})

# ---- Regression: non-square matrix $data (pre-2026-04-21 bug) ----

test_that("centrality_stability works when $data is a raw numeric matrix", {
  # For association methods (glasso/pcor/cor), build_network() stores $data
  # as a numeric matrix (not data.frame). When centrality_stability resamples
  # rows and re-invokes the estimator, .prepare_association_input was
  # failing the nrow==ncol check and silently producing NULL, which showed
  # up as all-NaN correlations or a zero-variance warning. Fixed by making
  # the matrix branch coerce non-square inputs to data.frame.
  skip_if_not_installed("glasso")
  set.seed(1)
  n <- 150; p <- 6
  Sigma <- diag(p)
  for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
    Sigma[i, j] <- Sigma[j, i] <- 0.4 ^ (j - i)
  }
  L   <- chol(Sigma)
  df  <- as.data.frame(matrix(rnorm(n * p), n, p) %*% L)
  colnames(df) <- paste0("V", seq_len(p))
  net <- build_network(df, method = "glasso", params = list(nlambda = 20))
  # Critical invariant: $data IS a matrix for association methods
  expect_true(is.matrix(net$data))
  # Fix should allow the re-estimation loop to run without all-NaN fallout
  cs <- suppressWarnings(suppressMessages(
    centrality_stability(net, iter = 5L)
  ))
  expect_true(inherits(cs, "net_stability"))
  # At least one correlation must be finite — if the bug regressed, every
  # re-estimation returns NULL and all correlations are NaN.
  corrs <- cs$correlations
  if (is.data.frame(corrs)) corrs <- corrs$correlation
  expect_true(any(is.finite(unlist(corrs))))
})
