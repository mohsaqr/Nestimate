# print()/summary() of build_mmm() name how covariates entered the fit

.mmm_label_data <- function() {
  set.seed(11)
  states <- c("a", "b", "c")
  seqs <- as.data.frame(matrix(sample(states, 60 * 6, replace = TRUE), 60, 6),
                        stringsAsFactors = FALSE)
  seqs$grp <- rep(c("g1", "g2"), each = 30)
  seqs
}

test_that("post-hoc covariates are labelled post-hoc, EM ones integrated", {
  d <- .mmm_label_data()
  posthoc <- build_mmm(d, k = 2, n_starts = 2, seed = 1, covariates = "grp",
                       covariate_effect = "posthoc")
  em <- build_mmm(d, k = 2, n_starts = 2, seed = 1, covariates = "grp",
                  covariate_effect = "em")
  expect_output(print(posthoc), "post-hoc")
  expect_output(print(em), "integrated")
  expect_output(summary(posthoc), "does not influence cluster membership")
  expect_output(summary(em), "influences cluster membership")
})
