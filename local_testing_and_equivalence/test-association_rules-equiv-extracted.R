# Equivalence test_that() blocks extracted from
# tests/testthat/test-association_rules.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

testthat::skip_on_cran()

# ---- Association Rules Tests ----

# Helper: standard test transactions
.make_ar_trans <- function() {
  list(
    c("bread", "milk", "eggs"),
    c("bread", "butter"),
    c("milk", "eggs", "butter"),
    c("bread", "milk", "eggs", "butter"),
    c("bread", "milk")
  )
}

# ---- 1. Basic functionality ----

test_that("association_rules works on tna::group_regulation via netobject", {
  skip_if_not_installed("tna")
  data(group_regulation, package = "tna")
  net <- build_network(group_regulation, method = "relative")
  rules <- association_rules(net, min_support = 0.3,
                             min_confidence = 0.5, min_lift = 1)
  expect_s3_class(rules, "net_association_rules")
  expect_true(rules$n_rules > 0)
})


# ---- 18. print and summary work ----

