# build_glasso(): the EBIC-selected graphical lasso (Friedman, Hastie and Tibshirani 2008; Foygel and Drton
# 2010). Nestimate solves it with psychnets' own solver (an internal delegation), so the references are the
# glasso package (the Fortran solver qgraph wraps, run here to thr = 1e-12) and qgraph::EBICglasso().
# Re-expressed from local_testing_and_equivalence/test-equiv-estimators.R ("glasso"), which compared 50 random
# inputs with qgraph through branches; here four fixed inputs, every number declared.
# The lambda path, EBIC path and precision matrix are read from the result's fields: no accessor returns them
# (API gap).

inputs <- list(
  list(n = 200L, p = 8L, rho = 0.8, seed = 32L),
  list(n = 250L, p = 7L, rho = 1, seed = 34L),
  list(n = 300L, p = 8L, rho = 1, seed = 31L)
)

test_that("partial correlations at the selected lambda equal glasso::glasso()", {
  skip_without_reference("glasso")
  invisible(lapply(inputs, \(input) {
    data <- psy_continuous(input$n, input$p, input$rho, input$seed)
    network <- build_glasso(data)
    fit <- glasso::glasso(stats::cor(data), rho = network$lambda_selected, penalize.diagonal = FALSE, thr = 1e-12)
    expected <- psy_partial(fit$wi)
    dimnames(expected) <- list(names(data), names(data))
    # Two coordinate-descent solvers, each run to its own threshold: they agree to about 1e-10.
    expect_reference(extract_transition_matrix(network), expected, reference = "glasso", covers = "build_glasso()",
                     quantity = sprintf("partial correlation, n = %d, p = %d", input$n, input$p), tolerance = 1e-8)
  }))
})

test_that("the EBIC path equals the extended BIC of glasso::glasso() fits", {
  skip_without_reference("glasso")
  invisible(lapply(inputs, \(input) {
    data <- psy_continuous(input$n, input$p, input$rho, input$seed)
    network <- build_glasso(data)
    expected <- psy_ebic_path(stats::cor(data), input$n, gamma = 0.5, lambdas = network$lambda_path)
    # EBIC values are of order 1e3; relative 1e-8 is the solvers' agreement carried through log det.
    expect_known_difference(network$ebic_path, expected, reference = "glasso", covers = "build_glasso()",
                            quantity = sprintf("EBIC path, n = %d, p = %d", input$n, input$p), tolerance = 1e-8,
                            measure = "relative",
                            reason = paste("The path is scanned with each lambda solved to 1e-4 (psychnets' scan_tol);",
                                           "only the selected lambda is refitted to 1e-8, so the EBIC of the other",
                                           "lambdas is that of an approximate fit (relative differences about 1e-6)."),
                            note = "EBIC written from Foygel and Drton (2010) on glasso fits at thr = 1e-12")
    expect_reference(network$lambda_selected, network$lambda_path[which.min(expected)], reference = "glasso",
                     covers = "build_glasso()", quantity = sprintf("selected lambda, n = %d, p = %d", input$n, input$p),
                     tolerance = 1e-12, measure = "relative")
  }))
})

test_that("the lambda path, the selected lambda and the network equal qgraph::EBICglasso()", {
  skip_without_reference("qgraph")
  invisible(lapply(inputs, \(input) {
    data <- psy_continuous(input$n, input$p, input$rho, input$seed)
    network <- build_glasso(data)
    reference <- qgraph::EBICglasso(stats::cor(data), n = input$n, gamma = 0.5, returnAllResults = TRUE)
    expect_reference(sort(network$lambda_path), sort(reference$lambda), reference = "qgraph", covers = "build_glasso()",
                     quantity = sprintf("lambda path, n = %d, p = %d", input$n, input$p), tolerance = 1e-12,
                     measure = "relative")
    expect_reference(network$lambda_selected, reference$lambda[which.min(reference$ebic)], reference = "qgraph",
                     covers = "build_glasso()", quantity = sprintf("selected lambda, n = %d, p = %d", input$n, input$p),
                     tolerance = 1e-12, measure = "relative")
    expected <- reference$optnet
    diag(expected) <- 0
    # qgraph runs glasso at its default thr = 1e-4, so its partial correlations are within about 1e-4 of the
    # optimum (the tolerance of the original test); glasso at thr = 1e-12 above is the tight comparison.
    expect_reference(extract_transition_matrix(network), expected, reference = "qgraph", covers = "build_glasso()",
                     quantity = sprintf("partial correlation, n = %d, p = %d", input$n, input$p), tolerance = 1e-4)
  }))
})

test_that("where qgraph's thr = 1e-4 leaves an edge at lambda_max, its selection differs", {
  skip_without_reference("qgraph")
  data <- psy_continuous(200L, 6L, 0.3, 11L)
  network <- build_glasso(data)
  reference <- qgraph::EBICglasso(stats::cor(data), n = 200L, gamma = 0.5, returnAllResults = TRUE)
  expect_known_difference(network$lambda_selected, reference$lambda[which.min(reference$ebic)], reference = "qgraph",
                          covers = "build_glasso()", quantity = "selected lambda, weak structure", tolerance = 1e-12,
                          measure = "relative",
                          reason = paste("At lambda_max = max|S_ij| the exact solution is the empty network; qgraph's glasso",
                                         "(thr = 1e-4) leaves one edge there, adds its EBIC penalty and selects a smaller",
                                         "lambda. Nestimate's EBIC at lambda_max is the exact one (n p)."))
  expect_reference(network$lambda_selected,
                   network$lambda_path[which.min(psy_ebic_path(stats::cor(data), 200L, 0.5, network$lambda_path))],
                   reference = "glasso", covers = "build_glasso()", quantity = "selected lambda, weak structure",
                   tolerance = 1e-12, measure = "relative")
})
