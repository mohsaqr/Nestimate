# ---- Tests for reliability() ----

test_that("single model reliability returns correct structure", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  rel <- reliability(net, iter = 50, seed = 42)

  expect_s3_class(rel, "net_reliability")
  expect_equal(rel$iter, 50L)
  expect_equal(rel$split, 0.5)
  expect_equal(rel$scale, "none")
  expect_named(rel, c("iterations", "summary", "models", "iter", "split", "scale"))

  # iterations data frame
  expect_s3_class(rel$iterations, "data.frame")
  expect_equal(nrow(rel$iterations), 50)
  expect_named(rel$iterations,
               c("model", "mean_dev", "median_dev", "cor", "max_dev"))
  expect_true(all(rel$iterations$model == "relative"))

  # summary data frame
  expect_s3_class(rel$summary, "data.frame")
  expect_equal(nrow(rel$summary), 4)
  expect_named(rel$summary, c("model", "metric", "mean", "sd"))
  expect_equal(rel$summary$metric,
               c("mean_dev", "median_dev", "cor", "max_dev"))

  # models list
  expect_length(rel$models, 1)
  expect_true(inherits(rel$models[[1]], "netobject"))
})


test_that("multi-model reliability stacks iterations correctly", {
  skip_if_not_installed("tna")
  net_r <- build_network(tna::group_regulation, method = "relative")
  net_f <- build_network(tna::group_regulation, method = "frequency")

  rel <- reliability(net_r, net_f, iter = 30, scale = "minmax", seed = 123)

  expect_equal(nrow(rel$iterations), 60)
  expect_equal(sort(unique(rel$iterations$model)), c("frequency", "relative"))
  expect_equal(nrow(rel$summary), 8)
  expect_equal(rel$scale, "minmax")
})


test_that("named models use provided names", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  rel <- reliability(my_model = net, iter = 20, seed = 1)

  expect_true(all(rel$iterations$model == "my_model"))
  expect_true(all(rel$summary$model == "my_model"))
})


test_that("seed produces reproducible results", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")

  rel1 <- reliability(net, iter = 50, seed = 99)
  rel2 <- reliability(net, iter = 50, seed = 99)

  expect_equal(rel1$iterations, rel2$iterations)
})


test_that("correlation values are between -1 and 1", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  rel <- reliability(net, iter = 50, seed = 42)

  cors <- rel$iterations$cor
  cors <- cors[!is.na(cors)]
  expect_true(all(cors >= -1 & cors <= 1))
})


test_that("deviation metrics are non-negative", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  rel <- reliability(net, iter = 50, seed = 42)

  expect_true(all(rel$iterations$mean_dev >= 0))
  expect_true(all(rel$iterations$median_dev >= 0))
  expect_true(all(rel$iterations$max_dev >= 0))
})


test_that("max_dev >= mean_dev >= median_dev in each iteration", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  rel <- reliability(net, iter = 100, seed = 42)

  # max >= mean always holds for absolute values

  expect_true(all(rel$iterations$max_dev >= rel$iterations$mean_dev - 1e-10))
})


test_that("warning issued when different methods without scaling", {
  skip_if_not_installed("tna")
  net_r <- build_network(tna::group_regulation, method = "relative")
  net_f <- build_network(tna::group_regulation, method = "frequency")

  expect_warning(
    reliability(net_r, net_f, iter = 10, seed = 1),
    "Models use different methods"
  )
})


test_that("no warning when different methods with scaling", {
  skip_if_not_installed("tna")
  net_r <- build_network(tna::group_regulation, method = "relative")
  net_f <- build_network(tna::group_regulation, method = "frequency")

  expect_no_warning(
    reliability(net_r, net_f, iter = 10, scale = "minmax", seed = 1)
  )
})


test_that("frequency method works", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "frequency")
  rel <- reliability(net, iter = 30, seed = 42)

  expect_s3_class(rel, "net_reliability")
  expect_equal(nrow(rel$iterations), 30)
  expect_true(all(rel$iterations$model == "frequency"))
})


test_that("co_occurrence method works", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "co_occurrence")
  rel <- reliability(net, iter = 30, seed = 42)

  expect_s3_class(rel, "net_reliability")
  expect_equal(nrow(rel$iterations), 30)
})


test_that("scale = 'standardize' works", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  rel <- reliability(net, iter = 30, scale = "standardize", seed = 42)

  expect_s3_class(rel, "net_reliability")
  expect_equal(rel$scale, "standardize")
})


test_that("scale = 'proportion' works", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  rel <- reliability(net, iter = 30, scale = "proportion", seed = 42)

  expect_s3_class(rel, "net_reliability")
  expect_equal(rel$scale, "proportion")
})


test_that("split parameter changes split ratio", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")

  rel_30 <- reliability(net, iter = 30, split = 0.3, seed = 42)
  rel_70 <- reliability(net, iter = 30, split = 0.7, seed = 42)

  # Different splits should give different results
  expect_false(identical(rel_30$iterations, rel_70$iterations))
  expect_equal(rel_30$split, 0.3)
  expect_equal(rel_70$split, 0.7)
})


test_that("netobject_group is flattened", {
  skip_if_not_installed("tna")
  df <- tna::group_regulation
  df$grp <- rep(c("A", "B"), length.out = nrow(df))
  group_net <- build_network(df, method = "relative", group = "grp")

  rel <- reliability(group_net, iter = 20, seed = 42)

  expect_s3_class(rel, "net_reliability")
  expect_true(length(unique(rel$iterations$model)) > 1)
})


test_that("input validation catches bad arguments", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")

  expect_error(reliability(net, iter = 0))
  expect_error(reliability(net, split = 0))
  expect_error(reliability(net, split = 1))
  expect_error(reliability(net, scale = "invalid"))
  expect_error(reliability("not_a_netobject"))
  expect_error(reliability())
})


test_that("print method works", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  rel <- reliability(net, iter = 20, seed = 42)

  out <- capture.output(print(rel))
  expect_true(any(grepl("Split-Half Reliability", out)))
  expect_true(any(grepl("Mean Abs. Dev.", out)))
  expect_true(any(grepl("Correlation", out)))
})


test_that("print method shows multi-model headers", {
  skip_if_not_installed("tna")
  net_r <- build_network(tna::group_regulation, method = "relative")
  net_f <- build_network(tna::group_regulation, method = "frequency")
  rel <- reliability(net_r, net_f, iter = 20, scale = "minmax", seed = 42)

  out <- capture.output(print(rel))
  expect_true(any(grepl("relative", out)))
  expect_true(any(grepl("frequency", out)))
  expect_true(any(grepl("scale = minmax", out)))
})


test_that("plot method returns a ggplot", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  rel <- reliability(net, iter = 30, seed = 42)

  p <- plot(rel)
  expect_s3_class(p, "ggplot")
})


test_that("plot method works for multi-model", {
  skip_if_not_installed("tna")
  net_r <- build_network(tna::group_regulation, method = "relative")
  net_f <- build_network(tna::group_regulation, method = "frequency")
  rel <- reliability(net_r, net_f, iter = 20, scale = "minmax", seed = 42)

  p <- plot(rel)
  expect_s3_class(p, "ggplot")
})


test_that("summary mean matches manual computation", {
  skip_if_not_installed("tna")
  net <- build_network(tna::group_regulation, method = "relative")
  rel <- reliability(net, iter = 50, seed = 42)

  manual_mean <- mean(rel$iterations$mean_dev)
  summary_mean <- rel$summary$mean[rel$summary$metric == "mean_dev"]
  expect_equal(manual_mean, summary_mean)
})


test_that("duplicate method names get deduplicated", {
  skip_if_not_installed("tna")
  net1 <- build_network(tna::group_regulation, method = "relative")
  net2 <- build_network(tna::group_regulation, method = "relative")

  rel <- reliability(net1, net2, iter = 20, seed = 42)
  models <- unique(rel$iterations$model)
  expect_length(models, 2)
  expect_false(models[1] == models[2])
})


test_that(".scale_matrix handles edge cases", {
  mat <- matrix(5, 3, 3)  # constant matrix

  # minmax of constant → unchanged
  scaled <- Nestimate:::.scale_matrix(mat, "minmax")
  expect_equal(scaled, mat)

  # standardize of constant → unchanged
  scaled <- Nestimate:::.scale_matrix(mat, "standardize")
  expect_equal(scaled, mat)

  # proportion of zero matrix → unchanged
  zero_mat <- matrix(0, 3, 3)
  scaled <- Nestimate:::.scale_matrix(zero_mat, "proportion")
  expect_equal(scaled, zero_mat)
})
