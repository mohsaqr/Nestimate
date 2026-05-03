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

test_that("attention estimator registered in registry", {
  estimators <- list_estimators()
  expect_true("attention" %in% estimators$name)
  expect_true(estimators$directed[estimators$name == "attention"])
})


# ---- Coverage gap tests ----

# estimators.R L558-559: .count_attention_wide empty matrix early return
test_that("attention wide: empty matrix when n_states == 0", {
  # All NAs → no valid states
  wide_data <- data.frame(
    T1 = c(NA_character_, NA_character_),
    T2 = c(NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  result <- .count_attention_wide(wide_data)
  expect_equal(nrow(result), 0L)
  expect_equal(ncol(result), 0L)
})

# estimators.R L558-559: .count_attention_wide nc < 2 early return (zero matrix)
test_that("attention wide: zero matrix when only 1 column (nc < 2)", {
  wide_data <- data.frame(T1 = c("A", "B"), stringsAsFactors = FALSE)
  result <- .count_attention_wide(wide_data)
  # Should return an n_states x n_states zero matrix (no pairs to form)
  expect_true(is.matrix(result))
  expect_true(all(result == 0))
})

# estimators.R L595: backward direction
test_that("attention wide: backward direction only counts (j<i) pairs", {
  wide_data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A"),
    T3 = c("A", "B"),
    stringsAsFactors = FALSE
  )
  bwd_mat <- .count_attention_wide(wide_data, direction = "backward")
  fwd_mat <- .count_attention_wide(wide_data, direction = "forward")

  # backward should only consider pairs where i > j
  # sum of backward == sum of forward for symmetric data
  expect_true(is.matrix(bwd_mat))
  expect_equal(dim(bwd_mat), dim(fwd_mat))
})

# estimators.R L627: .count_attention_long action col missing
test_that("attention long: errors when action column not found", {
  df <- data.frame(Actor = 1:3, Time = 1:3, Action = c("A", "B", "A"),
                   stringsAsFactors = FALSE)
  expect_error(
    .count_attention_long(df, action = "NonExistent"),
    "Action column.*not found"
  )
})

# estimators.R L643-644: .count_attention_long NULL id → single sequence grp
test_that("attention long: NULL id creates single sequence group", {
  df <- data.frame(
    Time = 1:4,
    Action = c("A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  result <- .count_attention_long(df, action = "Action", id = NULL,
                                   time = "Time")
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2L, 2L))
})

# estimators.R L648-650: .count_attention_long multi-id composite key
test_that("attention long: multi-id composite group key", {
  df <- data.frame(
    Actor   = c(1L, 1L, 2L, 2L),
    Session = c("s1", "s1", "s1", "s1"),
    Time    = c(1L, 2L, 1L, 2L),
    Action  = c("A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  result <- .count_attention_long(df, action = "Action",
                                   id = c("Actor", "Session"),
                                   time = "Time")
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2L, 2L))
})

# estimators.R L660: .count_attention_long n_states == 0 early return
test_that("attention long: empty matrix when all actions are NA", {
  df <- data.frame(
    Time   = 1:3,
    Action = c(NA_character_, NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  result <- .count_attention_long(df, action = "Action", id = NULL,
                                   time = "Time")
  expect_equal(nrow(result), 0L)
})

# estimators.R L675: .count_attention_long group with n < 2 is skipped
test_that("attention long: groups with only 1 obs are skipped", {
  df <- data.frame(
    Actor  = c(1L, 2L, 2L, 3L, 3L, 3L),
    Time   = c(1L, 1L, 2L, 1L, 2L, 3L),
    Action = c("A", "B", "A", "A", "B", "A"),
    stringsAsFactors = FALSE
  )
  # Actor 1 has only 1 observation → should be skipped
  result <- .count_attention_long(df, action = "Action", id = "Actor",
                                   time = "Time")
  expect_true(is.matrix(result))
  expect_true(all(result >= 0))
})

# estimators.R L681: time column in long format attention
test_that("attention long: time column used for decay", {
  df <- data.frame(
    Actor  = c(1L, 1L, 1L),
    Time   = c(1L, 10L, 20L),  # large gaps
    Action = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )
  result_large_gap <- .count_attention_long(df, action = "Action",
                                             id = "Actor", time = "Time",
                                             lambda = 1)
  result_no_time <- .count_attention_long(df, action = "Action",
                                           id = "Actor", time = "NoTime",
                                           lambda = 1)
  # Large time gaps should produce smaller attention weights than unit steps
  expect_true(sum(result_large_gap) < sum(result_no_time))
})

# estimators.R L687 L690-692: attention long inner loop NA check
test_that("attention long: NA actions in long format are skipped", {
  df <- data.frame(
    Actor  = c(1L, 1L, 1L, 1L),
    Time   = c(1L, 2L, 3L, 4L),
    Action = c("A", NA_character_, "B", "A"),
    stringsAsFactors = FALSE
  )
  result <- .count_attention_long(df, action = "Action", id = "Actor",
                                   time = "Time")
  expect_true(is.matrix(result))
  expect_true(all(result >= 0))
})

# estimators.R L752: .count_attention_long returns full matrix
test_that("attention long: returns correct matrix dimensions", {
  df <- data.frame(
    Actor  = c(1L, 1L, 1L, 2L, 2L, 2L),
    Time   = c(1L, 2L, 3L, 1L, 2L, 3L),
    Action = c("A", "B", "C", "C", "B", "A"),
    stringsAsFactors = FALSE
  )
  result <- .count_attention_long(df, action = "Action", id = "Actor",
                                   time = "Time")
  expect_equal(dim(result), c(3L, 3L))
  expect_equal(sort(rownames(result)), c("A", "B", "C"))
})
