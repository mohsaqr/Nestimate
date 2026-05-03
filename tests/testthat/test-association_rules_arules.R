testthat::skip_on_cran()

test_that("multi-consequent support is correct by manual count", {
  trans <- list(c("A","B","C"), c("A","B","C"), c("A","B"), c("B","C"), c("A","C"))
  rules <- association_rules(trans, min_support = 0.1, min_confidence = 0.0, min_lift = 0)
  r <- rules$rules
  n_cons <- lengths(strsplit(r$consequent, ", ", fixed = TRUE))
  multi <- r[n_cons > 1L, , drop = FALSE]
  if (nrow(multi) == 0) skip("No multi-consequent rules")
  trans_mat <- matrix(FALSE, 5, 3, dimnames = list(NULL, c("A","B","C")))
  for (i in seq_along(trans)) trans_mat[i, trans[[i]]] <- TRUE
  for (i in seq_len(nrow(multi))) {
    items <- unique(c(strsplit(multi$antecedent[i], ", ", fixed = TRUE)[[1]],
                      strsplit(multi$consequent[i], ", ", fixed = TRUE)[[1]]))
    mask <- Reduce(`&`, lapply(items, function(it) trans_mat[, it]))
    expect_equal(multi$support[i], sum(mask) / 5)
  }
})

test_that("multi-consequent confidence = support / ante_support", {
  trans <- list(c("A","B","C"), c("A","B","C"), c("A","B"), c("B","C"), c("A","C"))
  rules <- association_rules(trans, min_support = 0.1, min_confidence = 0.0, min_lift = 0)
  r <- rules$rules
  n_cons <- lengths(strsplit(r$consequent, ", ", fixed = TRUE))
  multi <- r[n_cons > 1L, , drop = FALSE]
  if (nrow(multi) == 0) skip("No multi-consequent rules")
  trans_mat <- matrix(FALSE, 5, 3, dimnames = list(NULL, c("A","B","C")))
  for (i in seq_along(trans)) trans_mat[i, trans[[i]]] <- TRUE
  for (i in seq_len(nrow(multi))) {
    ante_items <- strsplit(multi$antecedent[i], ", ", fixed = TRUE)[[1]]
    ante_mask <- Reduce(`&`, lapply(ante_items, function(it) trans_mat[, it]))
    expect_equal(multi$confidence[i], multi$support[i] / (sum(ante_mask) / 5), tolerance = 1e-12)
  }
})

test_that("multi-consequent lift is correct", {
  trans <- list(c("A","B","C"), c("A","B","C"), c("A","B"), c("B","C"), c("A","C"))
  rules <- association_rules(trans, min_support = 0.1, min_confidence = 0.0, min_lift = 0)
  r <- rules$rules
  n_cons <- lengths(strsplit(r$consequent, ", ", fixed = TRUE))
  multi <- r[n_cons > 1L, , drop = FALSE]
  if (nrow(multi) == 0) skip("No multi-consequent rules")
  trans_mat <- matrix(FALSE, 5, 3, dimnames = list(NULL, c("A","B","C")))
  for (i in seq_along(trans)) trans_mat[i, trans[[i]]] <- TRUE
  for (i in seq_len(nrow(multi))) {
    ante_items <- strsplit(multi$antecedent[i], ", ", fixed = TRUE)[[1]]
    cons_items <- strsplit(multi$consequent[i], ", ", fixed = TRUE)[[1]]
    ante_sup <- sum(Reduce(`&`, lapply(ante_items, function(it) trans_mat[, it]))) / 5
    cons_sup <- sum(Reduce(`&`, lapply(cons_items, function(it) trans_mat[, it]))) / 5
    expect_equal(multi$lift[i], multi$support[i] / (ante_sup * cons_sup), tolerance = 1e-12)
  }
})

test_that("subscript out of bounds fixed: high min_support with no 2-itemsets", {
  trans <- lapply(seq_len(20), function(i) sample(LETTERS[1:15], sample(2:4, 1)))
  expect_no_error(association_rules(trans, min_support = 0.8, min_confidence = 0.0, min_lift = 0))
})

test_that("$frequent is a tidy data frame", {
  trans <- list(c("A","B","C"), c("A","B"), c("B","C","D"), c("A","C","D"))
  rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.0, min_lift = 0)
  expect_true(is.data.frame(rules$frequent))
  expect_equal(names(rules$frequent), c("itemset", "size", "support", "count"))
  expect_type(rules$frequent$itemset, "character")
  expect_type(rules$frequent$size, "integer")
  expect_true(all(rules$frequent$support > 0))
  expect_true(all(rules$frequent$size[rules$frequent$size == 1] == 1))
})

test_that("$rules has character antecedent/consequent columns", {
  trans <- list(c("A","B","C"), c("A","B"), c("B","C","D"), c("A","C","D"))
  rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.0, min_lift = 0)
  expect_type(rules$rules$antecedent, "character")
  expect_type(rules$rules$consequent, "character")
  expect_true(all(nchar(rules$rules$antecedent) > 0))
})
