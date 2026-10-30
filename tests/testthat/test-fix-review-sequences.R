testthat::skip_on_cran()

test_that("logical missing positions preserve numeric sequence roles", {
  d <- data.frame(id = c("s1", "s2"), T1 = c(1, 2), T2 = NA,
                  T3 = c(2, 1), T4 = c(1, 2))
  positions <- c("T1", "T2", "T3", "T4")
  extracted <- extract_sequences(d, meta = TRUE)
  expect_identical(names(extracted), c("id", positions))
  expect_identical(extracted$id, d$id)
  expect_identical(extracted$T1, c("1", "2"))
  expect_true(all(is.na(extracted$T2)))
  net <- build_network(d, method = "relative")
  reference <- build_network(d, method = "relative", cols = positions)
  expect_identical(net$weights, reference$weights)
  expect_setequal(net$nodes$label, c("1", "2"))
  expect_identical(names(net$metadata), "id")
  expect_equal(sequence_patterns(d, len = 1, min_freq = 1),
               sequence_patterns(d, cols = positions, len = 1, min_freq = 1))
})

test_that("generated outcomes cannot overwrite adjustment covariates", {
  d <- expand.grid(hit = 0:1, score = 0:1, y = 0:1)
  d <- d[rep(seq_len(8), c(9, 4, 3, 2, 2, 3, 7, 9)), ]
  d$T1 <- ifelse(d$hit == 1, "A", "B")
  d$T2 <- "B"
  fit <- function(data, adjustment, last = FALSE) {
    pattern_outcome(data, outcome = if (last) "last_obs" else data$y,
                    cols = if (last) c("T1", "T2", "T3") else c("T1", "T2"),
                    len = 2, n = 1, min_freq = 1, start = "A",
                    adjust = adjustment)
  }
  reference <- fit(d, "score")
  renamed <- d
  names(renamed)[names(renamed) == "score"] <- ".outcome"
  expect_equal(fit(renamed, ".outcome"), reference)
  expect_gt(abs(reference$estimate), 0.9)
  d$T3 <- ifelse(d$y == 1, "pass", "fail")
  renamed$T3 <- d$T3
  expect_equal(fit(renamed, ".outcome", last = TRUE),
               fit(d, "score", last = TRUE))
  # Missing adjustments still drop their original rows, rather than being
  # replaced by the fully observed outcome.
  d$score[1] <- NA
  renamed$.outcome[1] <- NA
  expect_warning(ref_missing <- fit(d, "score"),
                 class = "nestimate_missing_outcome")
  expect_warning(actual_missing <- fit(renamed, ".outcome"),
                 class = "nestimate_missing_outcome")
  expect_equal(actual_missing, ref_missing)
})

test_that("stars-only state labels are refused before pattern parsing", {
  d <- data.frame(T1 = c("**", "A"), T2 = c("B", "B"))
  expect_error(sequence_patterns(d, len = 1, min_freq = 1),
               class = "nestimate_bad_input")
  expect_error(sequence_patterns(d, len = 2, min_freq = 1),
               class = "nestimate_bad_input")
  expect_error(sequence_patterns(d, pattern = "*->B", min_freq = 1),
               class = "nestimate_bad_input")
  expect_error(pattern_outcome(d, outcome = c(0, 1), min_freq = 1),
               class = "nestimate_bad_input")
  d$T1[1] <- "A*"
  expect_no_error(sequence_patterns(d, len = 1, min_freq = 1))
  expect_no_error(sequence_patterns(d, pattern = "*->B", min_freq = 1))
})
