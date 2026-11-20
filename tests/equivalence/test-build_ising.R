# build_ising(): the eLasso Ising model (van Borkulo et al. 2014), L1-penalised nodewise logistic regressions
# with EBIC selection and the AND rule. Nestimate's arithmetic is psychnets' (an internal delegation, calling
# glmnet); the reference is IsingFit::IsingFit(). Re-expressed from local_testing_and_equivalence/
# test-equiv-estimators.R ("ising", 30 random inputs) and test-estimator_ising-equiv-extracted.R.
# The thresholds are read from the result's field `thresholds`: no accessor returns them (API gap).

inputs <- list(
  list(n = 100L, p = 4L, rho = 0.3, seed = 201L),
  list(n = 150L, p = 5L, rho = 0.45, seed = 202L),
  list(n = 250L, p = 6L, rho = 0.4, seed = 203L)
)

test_that("Ising weights and thresholds equal IsingFit::IsingFit()", {
  skip_without_reference("IsingFit")
  invisible(lapply(inputs, \(input) {
    data <- psy_binary(input$n, input$p, input$rho, input$seed)
    network <- build_ising(data)
    reference <- IsingFit::IsingFit(data, plot = FALSE, progressbar = FALSE)
    expect_reference(extract_transition_matrix(network), reference$weiadj, reference = "IsingFit", covers = "build_ising()",
                     quantity = sprintf("edge weight, n = %d, p = %d", input$n, input$p), tolerance = 1e-10)
    expect_reference(network$thresholds, reference$thresholds, reference = "IsingFit", covers = "build_ising()",
                     quantity = sprintf("threshold, n = %d, p = %d", input$n, input$p), tolerance = 1e-10)
  }))
})
