# ---- extract_sequences() ----

.es_wide <- function() {
  data.frame(T1 = c("A", "B", "A"), T2 = c("B", "B", NA),
             T3 = c("A", NA, NA), score = c(1, 2, 3),
             stringsAsFactors = FALSE)
}

test_that("wide data frame: states only, numeric covariates set aside", {
  out <- extract_sequences(.es_wide())
  expect_identical(names(out), c("T1", "T2", "T3"))
  expect_identical(out$T2, c("B", "B", NA))
  with_meta <- extract_sequences(.es_wide(), meta = TRUE)
  expect_identical(names(with_meta), c("score", "T1", "T2", "T3"))
})

test_that("explicit cols select the states", {
  out <- extract_sequences(.es_wide(), cols = c("T1", "T2"))
  expect_identical(names(out), c("T1", "T2"))
  expect_error(extract_sequences(.es_wide(), cols = "T9"),
               class = "nestimate_bad_input")
})

test_that("character matrix keeps its states and drop-out NAs", {
  mat <- matrix(c("A", "B", NA, "C", "A", "B"), nrow = 2, byrow = TRUE)
  out <- extract_sequences(mat)
  expect_identical(dim(out), c(2L, 3L))
  expect_true(is.na(out[1, 3]))
  expect_identical(unname(unlist(out[2, ])), c("C", "A", "B"))
})

test_that("TraMineR void and missing markers read as NA, as in build_network", {
  # an stslist pads with "%" (void) and marks missing with "*"
  d <- data.frame(T1 = c("A", "B"), T2 = c("*", "A"), T3 = c("B", "%"),
                  stringsAsFactors = FALSE)
  out <- extract_sequences(d)
  expect_identical(sum(is.na(out)), 2L)
  expect_identical(sort(unique(stats::na.omit(unlist(out)))), c("A", "B"))
  expect_identical(sequence_metrics(d)$valid_n, c(2L, 2L))
})

test_that("long log and its prepared data give identical sequences", {
  long <- data.frame(id = c(1, 1, 1, 2, 2),
                     act = c("a", "b", "a", "b", "c"),
                     t = c(1, 2, 3, 1, 2))
  from_long <- extract_sequences(long, actor = "id", action = "act",
                                 time = "t")
  from_prepared <- extract_sequences(prepare(long, actor = "id",
                                             action = "act", time = "t"))
  expect_identical(from_long, from_prepared)
  expect_identical(unname(unlist(from_long[1, ])), c("a", "b", "a"))
})

test_that("a network returns the sequences it was built from", {
  seqs <- data.frame(T1 = c("A", "B", "A", "B"), T2 = c("B", "A", "A", "B"),
                     T3 = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
  net <- build_network(seqs, method = "relative")
  expect_identical(extract_sequences(net), seqs)
})

test_that("a grouped network stacks its groups and names the group", {
  seqs <- data.frame(T1 = c("A", "B", "A", "B"), T2 = c("B", "A", "A", "B"),
                     grp = c("x", "x", "y", "y"), stringsAsFactors = FALSE)
  nets <- build_network(seqs, method = "relative", group = "grp")
  out <- extract_sequences(nets, meta = TRUE)
  expect_identical(nrow(out), 4L)
  expect_identical(sort(out$grp), c("x", "x", "y", "y"))
})

test_that("unsupported input is refused by class", {
  expect_error(extract_sequences(42), class = "nestimate_bad_input")
  expect_error(extract_sequences(data.frame(T1 = c(NA, NA))),
               class = "nestimate_bad_input")
})
