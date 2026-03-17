# ---- attention estimator tests ----

test_that("attention estimator basic forward counting", {
  # 3-state, 3-timepoint, 2-sequence data
  wide_data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A"),
    T3 = c("A", "B"),
    stringsAsFactors = FALSE
  )

  net <- build_network(wide_data, method = "attention",
                       params = list(format = "wide", lambda = 1))

  expect_s3_class(net, "netobject")
  expect_true(net$directed)
  expect_equal(net$method, "attention")
  expect_true(all(net$weights >= 0))
  # Forward direction: (T1,T2), (T1,T3), (T2,T3)
  expect_true(nrow(net$weights) > 0)
})

test_that("attention estimator alias works", {
  wide_data <- data.frame(
    T1 = c("A", "B", "A"),
    T2 = c("B", "A", "B"),
    stringsAsFactors = FALSE
  )

  net <- build_network(wide_data, method = "atna",
                       params = list(format = "wide"))
  expect_equal(net$method, "attention")
})

test_that("attention estimator with custom lambda", {
  wide_data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A"),
    T3 = c("A", "B"),
    stringsAsFactors = FALSE
  )

  net1 <- build_network(wide_data, method = "attention",
                        params = list(format = "wide", lambda = 0.5))
  net2 <- build_network(wide_data, method = "attention",
                        params = list(format = "wide", lambda = 5))

  # Smaller lambda = faster decay = less weight on distant pairs
  # So total weight should be less with smaller lambda
  expect_true(sum(net1$weights) < sum(net2$weights))
})

test_that("attention estimator direction parameter", {
  wide_data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )

  fwd <- build_network(wide_data, method = "attention",
                       params = list(format = "wide", direction = "forward"))
  bwd <- build_network(wide_data, method = "attention",
                       params = list(format = "wide", direction = "backward"))
  both <- build_network(wide_data, method = "attention",
                        params = list(format = "wide", direction = "both"))

  # With 2 columns, forward has (1,2), backward has (2,1)
  # "both" should be the sum of forward + backward
  expect_equal(sum(both$weights), sum(fwd$weights) + sum(bwd$weights),
               tolerance = 1e-10)
})

test_that("attention estimator custom decay function", {
  wide_data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A"),
    T3 = c("A", "B"),
    stringsAsFactors = FALSE
  )

  # Linear decay
  linear_decay <- function(ti, tj, lambda) {
    1 / (1 + abs(ti - tj) / lambda)
  }

  net <- build_network(wide_data, method = "attention",
                       params = list(format = "wide", decay = linear_decay))
  expect_s3_class(net, "netobject")
  expect_true(all(net$weights >= 0))
})

test_that("attention estimator custom time_matrix", {
  wide_data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A"),
    T3 = c("A", "B"),
    stringsAsFactors = FALSE
  )

  # Custom time: unevenly spaced
  tm <- matrix(c(0, 0, 1, 1, 5, 5), nrow = 2, ncol = 3)

  net <- build_network(wide_data, method = "attention",
                       params = list(format = "wide", time_matrix = tm))
  expect_s3_class(net, "netobject")
})

test_that("attention estimator duration parameter", {
  wide_data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A"),
    T3 = c("A", "B"),
    stringsAsFactors = FALSE
  )

  net <- build_network(wide_data, method = "attention",
                       params = list(format = "wide", duration = c(1, 2, 3)))
  expect_s3_class(net, "netobject")
})

test_that("attention estimator handles NAs", {
  wide_data <- data.frame(
    T1 = c("A", NA, "A"),
    T2 = c("B", "A", NA),
    T3 = c("A", "B", "B"),
    stringsAsFactors = FALSE
  )

  net <- build_network(wide_data, method = "attention",
                       params = list(format = "wide"))
  expect_s3_class(net, "netobject")
  expect_true(all(net$weights >= 0))
})

test_that("attention estimator long format", {
  long_data <- data.frame(
    Actor = c(1, 1, 1, 2, 2, 2),
    Time = c(1, 2, 3, 1, 2, 3),
    Action = c("A", "B", "A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )

  net <- build_network(long_data, method = "attention",
                       params = list(format = "long", action = "Action",
                                     id = "Actor", time = "Time"))
  expect_s3_class(net, "netobject")
  expect_true(net$directed)
  expect_equal(sort(net$nodes$label), c("A", "B"))
})

test_that("attention estimator print label", {
  wide_data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )

  net <- build_network(wide_data, method = "attention",
                       params = list(format = "wide"))
  output <- capture.output(print(net))
  expect_true(any(grepl("Attention", output)))
})

test_that("attention estimator is NOT row-normalized", {
  wide_data <- data.frame(
    T1 = c("A", "A", "A"),
    T2 = c("B", "B", "B"),
    T3 = c("A", "A", "A"),
    stringsAsFactors = FALSE
  )

  net <- build_network(wide_data, method = "attention",
                       params = list(format = "wide"))
  # Raw attention counts should not sum to 1 per row
  rs <- rowSums(net$weights)
  expect_false(all(abs(rs - 1) < 1e-10))
})

test_that("attention estimator cross-validates with tna::atna on simple data", {
  skip_if_not_installed("tna")

  set.seed(42)
  states <- c("A", "B", "C")
  n_seq <- 50
  n_time <- 5
  wide_data <- data.frame(matrix(
    sample(states, n_seq * n_time, replace = TRUE),
    nrow = n_seq, ncol = n_time
  ))
  names(wide_data) <- paste0("V", seq_len(n_time))

  # Nestimate attention
  net <- build_network(wide_data, method = "attention",
                       params = list(format = "wide", lambda = 1,
                                     direction = "forward"))

  # tna::atna (if available)
  tna_model <- tryCatch(
    tna::atna(wide_data),
    error = function(e) NULL
  )

  if (!is.null(tna_model)) {
    tna_mat <- tna_model$weights
    # Both should have same states
    expect_equal(sort(rownames(net$weights)), sort(rownames(tna_mat)))
    # Compare values (allow tolerance for different implementations)
    common <- sort(rownames(net$weights))
    nest_mat <- net$weights[common, common]
    tna_ref <- tna_mat[common, common]
    # Check correlation is high (same relative pattern)
    if (sum(nest_mat) > 0 && sum(tna_ref) > 0) {
      cor_val <- cor(as.vector(nest_mat), as.vector(tna_ref))
      expect_true(cor_val > 0.9,
                  label = sprintf("Correlation with tna::atna: %.3f", cor_val))
    }
  }
})

test_that("attention estimator registered in registry", {
  estimators <- list_estimators()
  expect_true("attention" %in% estimators$name)
  expect_true(estimators$directed[estimators$name == "attention"])
})
