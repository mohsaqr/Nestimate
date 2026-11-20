# build_cor(), build_pcor(): Pearson correlation and partial correlation networks (the arithmetic is psychnets',
# an internal delegation). References: stats::cor() and corpcor::cor2pcor(), and the partial correlations of the
# inverse correlation matrix written from the definition. Re-expressed from
# local_testing_and_equivalence/test-equiv-estimators.R ("cor", "pcor": 100 random inputs, cell by cell);
# here four fixed inputs.

inputs <- list(
  list(n = 50L, p = 4L, rho = 0.5, seed = 101L),
  list(n = 100L, p = 5L, rho = 1, seed = 102L),
  list(n = 200L, p = 6L, rho = 1.5, seed = 103L),
  list(n = 120L, p = 8L, rho = 0.8, seed = 104L)
)

test_that("correlation networks equal stats::cor() with a zero diagonal", {
  invisible(lapply(inputs, \(input) {
    data <- psy_continuous(input$n, input$p, input$rho, input$seed)
    expected <- stats::cor(data)
    diag(expected) <- 0
    expect_reference(extract_transition_matrix(build_cor(data)), expected, reference = "stats", covers = "build_cor()",
                     quantity = sprintf("correlation, n = %d, p = %d", input$n, input$p), tolerance = 1e-10)
  }))
})

test_that("partial correlation networks equal corpcor::cor2pcor()", {
  skip_without_reference("corpcor")
  invisible(lapply(inputs, \(input) {
    data <- psy_continuous(input$n, input$p, input$rho, input$seed)
    expected <- corpcor::cor2pcor(stats::cor(data))
    dimnames(expected) <- list(names(data), names(data))
    diag(expected) <- 0
    expect_reference(extract_transition_matrix(build_pcor(data)), expected, reference = "corpcor", covers = "build_pcor()",
                     quantity = sprintf("partial correlation, n = %d, p = %d", input$n, input$p), tolerance = 1e-10)
  }))
})

test_that("partial correlation networks equal -K_ij / sqrt(K_ii K_jj) of the inverse correlation matrix", {
  invisible(lapply(inputs, \(input) {
    data <- psy_continuous(input$n, input$p, input$rho, input$seed)
    expected <- psy_partial(solve(stats::cor(data)))
    expect_reference(extract_transition_matrix(build_pcor(data)), expected,
                     reference = "independent: partial correlations of the inverse correlation matrix, base R",
                     covers = "build_pcor()", quantity = sprintf("partial correlation, n = %d, p = %d", input$n, input$p),
                     tolerance = 1e-10)
  }))
})
