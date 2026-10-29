# State identity and label replay across model and resampling interfaces.
testthat::skip_on_cran()

.review_state_sequences <- function() {
  data.frame(V1 = c(1, 3, 5, 1), V2 = c(3, 5, 1, 5),
              V3 = c(5, 1, 3, 3))
}

test_that("numeric netobjects keep literal states in higher-order input", {
  raw <- .review_state_sequences()
  net <- build_network(raw, method = "relative")
  expect_identical(.coerce_sequence_input(net), raw)
  expect_false(anyNA(.coerce_sequence_input(net)))

  hypa_raw <- build_hypa(raw, order = 2L, min_count = 1L)
  hypa_net <- build_hypa(net, order = 2L, min_count = 1L)
  expect_gt(nrow(hypa_raw$scores), 0L)
  expect_equal(hypa_net$scores, hypa_raw$scores)
  expect_equal(build_hon(net, max_order = 2L)$weights,
                build_hon(raw, max_order = 2L)$weights)
  expect_equal(build_mogen(net, max_order = 2L)$weights,
                build_mogen(raw, max_order = 2L)$weights)

  grouped <- structure(list(first = net, second = net),
                       class = "netobject_group")
  parsed <- .cs_parse_input(grouped, NULL)
  expect_equal(unname(parsed$seq_mat), unname(rbind(as.matrix(raw), as.matrix(raw))))
})

test_that("tna state codes are decoded even when codes also name states", {
  skip_if_not_installed("tna")
  raw <- as.data.frame(matrix(rep(as.character(1:12), 8L),
                              nrow = 8L, byrow = TRUE))
  group <- rep(c("first", "second"), each = 4L)
  model <- tna::tna(raw)
  expect_false(identical(model$labels, as.character(1:12)))
  expect_identical(unname(as.matrix(.coerce_sequence_input(model))),
                   unname(as.matrix(raw)))
  actual <- sequence_compare(model, group = group, sub = 2L,
                              min_freq = 1L, test = "none")$patterns
  expected <- sequence_compare(raw, group = group, sub = 2L,
                                min_freq = 1L, test = "none")$patterns
  expect_equal(actual, expected)

  numeric_raw <- as.data.frame(lapply(raw, as.numeric))
  numeric_model <- tna::tna(numeric_raw)
  expect_identical(unname(as.matrix(.coerce_sequence_input(numeric_model))),
                   unname(as.matrix(raw)))
  expect_equal(sequence_compare(numeric_model, group = group, sub = 2L,
                                 min_freq = 1L, test = "none")$patterns,
                expected)
})

test_that("generic cograph numeric states require explicit encoding provenance", {
  raw <- .review_state_sequences()
  net <- build_network(raw, method = "relative")
  generic <- net
  class(generic) <- "cograph_network"
  expect_identical(.as_netobject(generic)$data, raw)
  expect_identical(.coerce_sequence_input(generic), raw)

  encoded <- generic
  encoded$data[] <- lapply(raw, function(col) match(as.character(col), c("1", "3", "5")))
  attr(encoded$data, "labels") <- c("1", "3", "5")
  expected <- as.data.frame(lapply(raw, as.character))
  expect_equal(.coerce_sequence_input(encoded), expected)

  skip_if_not_installed("tna")
  skip_if_not_installed("cograph")
  converted <- cograph::as_cograph(tna::tna(expected))
  expect_equal(.coerce_sequence_input(converted), expected)
})

test_that("numeric state relabeling preserves observed and bootstrap estimates", {
  raw <- data.frame(A = rep(1, 4L), B = rep(2, 4L), C = rep(1, 4L))
  labels <- c("1" = "Alpha", "2" = "Beta")
  net <- build_network(raw, method = "frequency", state_cols = names(raw),
                       labels = labels)
  expect_identical(names(net$data), names(raw))
  expect_identical(net$params$cols, names(raw))
  expect_identical(net$data$A, rep("Alpha", 4L))
  expect_identical(net$data$B, rep("Beta", 4L))
  expect_equal(bootstrap_network(net, iter = 6L, seed = 1L)$mean, net$weights)
  replay <- do.call(get_estimator(net$method)$fn,
                    c(list(data = net$data), net$params))$matrix
  expect_equal(replay, net$weights)
})

test_that("wide state labels do not rename positional sequence columns", {
  raw <- data.frame(A = rep("A", 4L), B = rep("B", 4L), C = rep("A", 4L))
  net <- build_network(raw, method = "frequency", state_cols = names(raw),
                       labels = c(A = "Alpha", B = "Beta"))
  expect_identical(names(net$data), names(raw))
  expect_identical(net$params$cols, names(raw))
  expect_equal(bootstrap_network(net, iter = 6L, seed = 1L)$mean, net$weights)
})

test_that("renamed marker and alphabet parameters replay their observed counts", {
  raw <- data.frame(V1 = rep("A", 4L), V2 = rep("B", 4L), V3 = rep("A", 4L))
  net <- build_network(raw, method = "frequency", start = TRUE, end = TRUE,
                       params = list(alphabet = c("A", "B")),
                       labels = c(A = "Alpha", Start = "Begin", End = "Finish"))
  expect_identical(net$params$begin_state, "Begin")
  expect_identical(net$params$end_state, "Finish")
  expect_identical(net$params$alphabet, c("Alpha", "B"))
  expect_equal(bootstrap_network(net, iter = 6L, seed = 1L)$mean, net$weights)
  replay <- do.call(get_estimator(net$method)$fn,
                    c(list(data = net$data), net$params))$matrix
  expect_equal(replay, net$weights)
})

test_that("onehot label remapping updates codes without changing indicators", {
  raw <- data.frame(id = rep(seq_len(4L), each = 6L),
                     A = rep(c(1, 0), 12L), B = rep(c(0, 1), 12L))
  plain <- build_network(raw, method = "frequency", format = "onehot",
                         actor = "id", codes = c("A", "B"), window_size = 2L)
  named <- build_network(raw, method = "frequency", format = "onehot",
                         actor = "id", codes = c("A", "B"), window_size = 2L,
                         labels = c(A = "Alpha", B = "Beta"))
  expect_equal(unname(as.matrix(named$data)), unname(as.matrix(plain$data)))
  expect_identical(names(named$data), c("Alpha", "Beta"))
  expect_identical(named$params$codes, c("Alpha", "Beta"))
  expect_identical(named$build_args$codes, c("Alpha", "Beta"))
  expect_null(named$build_args$action)
  expect_identical(named$params$actor, "id")
  expect_equal(unname(bootstrap_network(named, iter = 6L, seed = 1L)$mean),
                unname(bootstrap_network(plain, iter = 6L, seed = 1L)$mean))
})

test_that("association relabeling renames variables and preserves measurements", {
  set.seed(307L)
  raw <- data.frame(A = stats::rnorm(50L), B = stats::rnorm(50L),
                     C = stats::rnorm(50L))
  plain <- build_network(raw, method = "cor")
  named <- build_network(raw, method = "cor", labels = c(A = "Alpha", B = "Beta"))
  expect_equal(unname(named$data), unname(plain$data))
  expect_identical(colnames(named$data), c("Alpha", "Beta", "C"))
  expect_null(named$build_args$actor)
  expect_identical(names(named$predictability), c("Alpha", "Beta", "C"))
  expect_equal(unname(bootstrap_network(named, iter = 6L, seed = 1L)$mean),
                unname(bootstrap_network(plain, iter = 6L, seed = 1L)$mean))
})
