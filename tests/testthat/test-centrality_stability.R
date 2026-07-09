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
  expect_s3_class(c1, "net_centrality")
  expect_equal(nrow(c1), 3)
  # default is the compact trio
  expect_identical(setdiff(names(c1), "state"),
                   c("InStrength", "Betweenness", "Diffusion"))
  # measures = "all" expands to every built-in measure
  c_all <- net_centrality(net, measures = "all")
  expect_true(all(c("InStrength", "OutStrength", "Betweenness",
                    "BetweennessRSP", "Diffusion", "Clustering") %in%
                    names(c_all)))
})

test_that("centrality.netobject returns correct undirected defaults (L163-177)", {
  set.seed(42)
  panel <- data.frame(V1 = rnorm(50), V2 = rnorm(50), V3 = rnorm(50))
  net_ud <- build_network(panel, method = "cor")
  c2 <- net_centrality(net_ud)
  expect_true(is.data.frame(c2))
  expect_s3_class(c2, "net_centrality")
  # Betweenness is part of the default trio
  expect_true("Betweenness" %in% names(c2))
  # the full closeness family appears under measures = "all"
  expect_true(all(c("Closeness", "Betweenness") %in%
                    names(net_centrality(net_ud, measures = "all"))))
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
  expect_s3_class(c3, "net_centrality_group")
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

test_that("net_centrality matches tna centralities for all shared measures", {
  skip_if_not_installed("tna")
  measures <- c("OutStrength", "InStrength", "ClosenessIn", "ClosenessOut",
                "Closeness", "Betweenness", "BetweennessRSP", "Diffusion",
                "Clustering")
  mat <- matrix(c(
    0,   .20, 0,   .35,
    .90, 0,   .40, .10,
    .10, .30, 0,   .60,
    .20, .15, .50, 0
  ), nrow = 4L, byrow = TRUE,
  dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  net <- .wrap_netobject(mat, method = "relative", directed = TRUE,
                         data = NULL)

  ours <- suppressMessages(net_centrality(
    net, measures = measures, normalize_diffusion = FALSE
  ))
  ref <- tna::centralities(mat, measures = measures)

  expect_equal(as.data.frame(ours)[, measures],
               as.data.frame(ref)[, measures],
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("Diffusion is normalized by default but can match raw tna output", {
  skip_if_not_installed("tna")
  mat <- matrix(c(
    0,   .20, 0,   .35,
    .90, 0,   .40, .10,
    .10, .30, 0,   .60,
    .20, .15, .50, 0
  ), nrow = 4L, byrow = TRUE,
  dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  net <- .wrap_netobject(mat, method = "relative", directed = TRUE,
                         data = NULL)

  def <- suppressMessages(net_centrality(net, measures = "Diffusion"))
  raw <- suppressMessages(net_centrality(
    net, measures = "Diffusion", normalize_diffusion = FALSE
  ))
  ref_raw <- tna::centralities(mat, measures = "Diffusion")
  ref_norm <- tna::centralities(mat, measures = "Diffusion",
                                normalize = TRUE)

  expect_equal(def$Diffusion, ref_norm$Diffusion,
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(raw$Diffusion, ref_raw$Diffusion,
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("net_centrality plots work for single and grouped outputs", {
  seqs <- data.frame(
    V1 = c("A","B","A","C","B","A"),
    V2 = c("B","C","B","A","C","B"),
    V3 = c("C","A","C","B","A","C"),
    grp = c("X","X","X","Y","Y","Y")
  )
  net <- build_network(seqs[1:3], method = "relative")
  cent <- suppressMessages(net_centrality(
    net, measures = c("InStrength", "OutStrength", "Diffusion")
  ))
  expect_s3_class(plot(cent), "ggplot")
  expect_s3_class(plot(cent, type = "line"), "ggplot")
  expect_s3_class(plot(cent, type = "profile"), "ggplot")  # back-compat alias
  expect_s3_class(plot(cent, type = "heatmap"), "ggplot")
  expect_s3_class(plot(cent, type = "bar", drop_zero = TRUE), "ggplot")

  nets <- build_network(seqs, method = "relative", group = "grp")
  cents <- suppressMessages(net_centrality(
    nets, measures = c("InStrength", "OutStrength", "Diffusion")
  ))
  expect_s3_class(plot(cents), "ggplot")
  expect_s3_class(plot(cents, type = "line"), "ggplot")
  expect_s3_class(plot(cents, type = "profile"), "ggplot")  # back-compat alias
  expect_s3_class(plot(cents, type = "delta"), "ggplot")
})

test_that("centrality heatmap and delta views behave correctly", {
  seqs <- data.frame(
    V1 = c("A","B","A","C","B","A"),
    V2 = c("B","C","B","A","C","B"),
    V3 = c("C","A","C","B","A","C"),
    grp = c("X","X","X","Y","Y","Y")
  )
  # drop_zero removes an all-zero measure panel
  net <- build_network(seqs[1:3], method = "relative")
  cent <- suppressMessages(net_centrality(
    net, measures = c("OutStrength", "Betweenness")
  ))
  # Betweenness is all zero on this tiny ring -> dropped
  p_keep <- plot(cent, type = "bar", drop_zero = FALSE)
  p_drop <- plot(cent, type = "bar", drop_zero = TRUE)
  expect_true(length(levels(p_keep$data$measure)) >=
              length(levels(droplevels(p_drop$data$measure))))

  # delta supports 3+ groups via deviation-from-mean
  three <- within(seqs, grp <- c("X","X","Y","Y","Z","Z"))
  nets3 <- build_network(three, method = "relative", group = "grp")
  cents3 <- suppressMessages(net_centrality(nets3, measures = "OutStrength"))
  expect_s3_class(plot(cents3, type = "delta", labels = TRUE), "ggplot")
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

test_that("centrality_stability default is the trio and accepts measures = 'all'", {
  set.seed(1)
  seqs <- as.data.frame(matrix(sample(c("A","B","C","D"), 160, TRUE), ncol = 4))
  net <- build_network(seqs, method = "relative")

  # Default measures. Kept at the 0.6.0 trio (OutStrength, not Diffusion):
  # htna's CRAN release compares its own explicit trio against this default,
  # so changing it breaks the reverse dependency. See CLAUDE.md, "Reverse
  # Dependency: htna".
  expect_identical(eval(formals(centrality_stability)$measures),
                   c("InStrength", "OutStrength", "Betweenness"))
  s_def <- suppressWarnings(suppressMessages(
    centrality_stability(net, iter = 40, seed = 1)))
  expect_identical(s_def$measures, c("InStrength", "OutStrength", "Betweenness"))

  # "all" expands to every built-in measure (and no longer errors on
  # degenerate resamples where Diffusion/RSP can be NA)
  s_all <- suppressWarnings(suppressMessages(
    centrality_stability(net, measures = "all", iter = 40, seed = 1)))
  expect_true(all(c("OutStrength", "InStrength", "Closeness", "Betweenness",
                    "BetweennessRSP", "Diffusion", "Clustering") %in%
                    s_all$measures))
})
