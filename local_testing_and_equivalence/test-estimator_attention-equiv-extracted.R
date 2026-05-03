# Equivalence test_that() blocks extracted from
# tests/testthat/test-estimator_attention.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

# ---- attention estimator tests ----

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

