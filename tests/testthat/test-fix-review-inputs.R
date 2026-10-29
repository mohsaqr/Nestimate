testthat::skip_on_cran()

test_that("grouped wide input resolves state columns before its alphabet", {
  d <- data.frame(T1 = factor(c("A", "B", "A", "B")),
                  T2 = c("B", "A", "B", "A"),
                  age = c(20, 30, 40, 50), g = c("x", "x", "y", "y"))
  grouped <- build_network(d, "frequency", group = "g")
  invisible(lapply(names(grouped), function(g) {
    child <- grouped[[g]]
    expected <- build_network(d[d$g == g, c("T1", "T2")], "frequency")
    expect_equal(child$weights, expected$weights)
    expect_identical(names(child$data), c("T1", "T2"))
    expect_identical(child$metadata$age, d$age[d$g == g])
  }))
  # Explicit numeric positions are states, and an explicit alphabet is kept.
  numeric_states <- data.frame(T1 = c(1, 3, 1, 3), T2 = c(3, 1, 3, 1), g = d$g)
  expect_setequal(rownames(build_network(numeric_states, "frequency",
                                         group = "g")[[1]]$weights), c("1", "3"))
  alphabet <- build_network(d, "frequency", group = "g",
                            params = list(alphabet = c("A", "B", "C")))
  expect_identical(rownames(alphabet[[1]]$weights), c("A", "B", "C"))
  overlap <- d
  overlap$T1 <- as.character(overlap$T1)
  overlap$T1[1] <- "20"
  overlap_net <- build_network(overlap, "frequency", group = "g")
  expect_identical(names(overlap_net[[1]]$data), c("T1", "T2"))
  expect_identical(overlap_net[[1]]$metadata$age, c(20, 30))
})

test_that("MGM level is checked against final variables, excluding metadata and ids", {
  skip_if_not_installed("glmnet")
  set.seed(2701)
  d <- data.frame(id = seq_len(100), x = rnorm(100), y = rnorm(100),
                  z = rnorm(100), group = rep(c("a", "b"), each = 50))
  expected <- build_network(d[c("x", "y", "z")], "mgm")
  by_metadata <- build_network(d[setdiff(names(d), "group")], "mgm",
                                metadata_cols = "id", level = rep(1L, 3))
  by_actor <- build_network(d[setdiff(names(d), "group")], "mgm",
                             actor = "id", level = rep(1L, 3))
  expect_equal(by_metadata$weights, expected$weights)
  expect_equal(by_actor$weights, expected$weights)
  grouped <- build_network(d, "mgm", group = "group", actor = "id",
                             level = rep(1L, 3))
  expect_true(all(vapply(grouped, function(n) identical(rownames(n$weights),
                                                        c("x", "y", "z")), logical(1))))
  expect_error(build_network(d, "mgm", group = "group", metadata_cols = "id",
                              level = rep(1L, 4)), "length 3")
})

test_that("NCT refuses infinity before pairwise or listwise estimation", {
  skip_if_not_installed("Matrix")
  set.seed(2702)
  d <- matrix(rnorm(120), 40, 3, dimnames = list(NULL, c("x", "y", "z")))
  invisible(lapply(c(Inf, -Inf), function(value) {
    bad <- d
    bad[1, 1] <- value
    expect_error(nct(bad, d, iter = 3), class = "nestimate_bad_input")
    expect_error(nct(d, bad, iter = 3, na_method = "listwise"),
                  class = "nestimate_bad_input")
  }))
  # Supported missing observations are still handled by the selected policy.
  missing <- d
  missing[1, 1] <- NA_real_
  set.seed(2702)
  fit <- nct(missing, d, iter = 3)
  expect_true(all(is.finite(fit$nw1)))
})

test_that("mosaic_plot network calls honor labels and reject data-frame styling", {
  d <- data.frame(T1 = c("A", "A", "B", "B", "C", "C"),
                  T2 = c("B", "C", "A", "C", "A", "B"))
  net <- build_network(d, "frequency")
  set.seed(2703)
  p <- mosaic_plot(net, xlab = "Destination", ylab = "Source")
  expect_identical(p$labels$x, "Destination")
  expect_identical(p$labels$y, "Source")
  expect_error(mosaic_plot(net, colors = c(A = "red")),
                class = "nestimate_bad_input")
  expect_error(mosaic_plot(net, show_labels = FALSE),
                class = "nestimate_bad_input")
  expect_error(mosaic_plot(as.table(net$weights), x_var = "a"),
                class = "nestimate_bad_input")
  # The data-frame method keeps its own option contract.
  tidy <- data.frame(group = c("g1", "g1", "g2", "g2"),
                     state = c("A", "B", "A", "B"), count = c(3, 1, 1, 3))
  p_df <- mosaic_plot(tidy, x_var = "group", y_var = "state",
                      weight = "count", colors = c(A = "red", B = "blue"),
                      xlab = "Groups", show_labels = FALSE)
  expect_identical(p_df$labels$x, "Groups")
})

test_that("comparison global labels describe their actual statistic", {
  d <- data.frame(T1 = c("A", "B", "A", "B"),
                  T2 = c("B", "A", "A", "B"))
  net <- build_network(d, "frequency")
  cmp <- compare_networks(first = net, second = net, test = "permutation",
                           iter = 5, seed = 2704)
  row <- cmp$global[cmp$global$key == "perm_M", ]
  expect_identical(row$metric, "M (total absolute difference)")
  expect_equal(row$value, 0)
})

test_that("custom variable estimators retain the numeric observations for permutation", {
  estimator <- function(data, ...) {
    weights <- stats::cor(as.matrix(data))
    diag(weights) <- 0
    list(matrix = weights, nodes = colnames(weights), directed = FALSE,
          cleaned_data = data)
  }
  name <- "review_variable_correlation"
  register_estimator(name, estimator, "Numeric variable regression control", FALSE)
  on.exit(remove_estimator(name), add = TRUE)
  set.seed(2705)
  d <- data.frame(x = rnorm(40), y = rnorm(40), z = rnorm(40))
  net <- build_network(d, name)
  expect_identical(net$data, d)
  comparison <- permutation(net, net, iter = 5, seed = 2705)
  expect_true(all(comparison$p_values == 1))
  expect_true(all(comparison$global$p_value == 1))
})

test_that("MGM's documented imputation matches explicitly prepared observations", {
  skip_if_not_installed("glmnet")
  set.seed(2706)
  d <- data.frame(x = rnorm(100), y = rnorm(100), z = rnorm(100))
  d$y <- d$x + d$y / 2
  d$x[seq_len(20)] <- NA_real_
  imputed <- d
  imputed$x[is.na(imputed$x)] <- mean(imputed$x, na.rm = TRUE)
  expect_equal(build_network(d, "mgm", type = rep("g", 3))$weights,
                build_network(imputed, "mgm", type = rep("g", 3))$weights)
})

test_that("custom sequence states can share names with positional columns", {
  name <- "review_custom_sequence"
  frequency <- get_estimator("frequency")$fn
  register_estimator(name, frequency, "Custom sequence regression control", TRUE)
  on.exit(remove_estimator(name), add = TRUE)
  d <- data.frame(V1 = rep("V1", 4), V2 = rep("V2", 4), V3 = rep("V1", 4))
  net <- build_network(d, name)
  expect_identical(net$data, d)
  expect_equal(bootstrap_network(net, iter = 5, seed = 2707)$mean, net$weights)
  numeric <- data.frame(`1` = rep(1, 4), `2` = rep(2, 4), `3` = rep(1, 4),
                        check.names = FALSE)
  numeric_net <- build_network(numeric, name)
  expect_identical(numeric_net$data, numeric)
  expect_equal(bootstrap_network(numeric_net, iter = 5, seed = 2707)$mean,
                numeric_net$weights)
})

test_that("wide session identifiers stay out of sequence positions", {
  d <- data.frame(T1 = c("A", "B", "A", "B"),
                  T2 = c("B", "A", "B", "A"), session = seq_len(4))
  expected <- build_network(d[c("T1", "T2")], "frequency")
  auto <- build_network(d, "frequency")
  expect_equal(auto$weights, expected$weights)
  expect_identical(names(auto$data), c("T1", "T2"))
  expect_identical(auto$metadata$session, d$session)
  d$id <- c("p1", "p1", "p2", "p2")
  names(d)[names(d) == "session"] <- "day"
  explicit <- build_network(d, "frequency", actor = "id", session = "day")
  expect_equal(explicit$weights, expected$weights)
  expect_identical(names(explicit$data), c("T1", "T2"))
})

test_that("MGM retains identifiers until multilevel decomposition uses them", {
  skip_if_not_installed("glmnet")
  set.seed(2708)
  d <- data.frame(id = rep(seq_len(20), each = 3),
                  x = rnorm(60), y = rnorm(60), z = rnorm(60))
  between_data <- stats::aggregate(d[c("x", "y", "z")], list(id = d$id), mean)
  between_data <- between_data[c("x", "y", "z")]
  within_data <- as.data.frame(lapply(d[c("x", "y", "z")], function(col) {
    col - ave(col, d$id, FUN = mean)
  }))
  both <- build_network(d, "mgm", actor = "id", level = "both",
                         params = list(level = rep(1L, 3)))
  expect_identical(rownames(both$between$weights), c("x", "y", "z"))
  expect_equal(both$between$weights, build_network(between_data, "mgm")$weights)
  expect_equal(both$within$weights, build_network(within_data, "mgm")$weights)
})
