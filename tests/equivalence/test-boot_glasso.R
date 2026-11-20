# boot_glasso(): the EBIC graphical lasso with nonparametric bootstrap intervals, edge inclusion and case-dropping
# stability (Epskamp, Borsboom and Fried 2018). The original network is deterministic and is compared with the
# glasso package; the bootstrap summaries are compared with the percentile interval and inclusion share written
# from the definition on the bootstrap replicates the function returns (the resampling itself is drawn inside the
# function, with only a seed to control it, so no outside program can replay it).
# Re-expressed from local_testing_and_equivalence/test-equiv-boot_glasso.R, test-boot_glasso-equiv-extracted.R and
# test-equiv-moved-boot_glasso.R (qgraph agreement as a share of configurations; quantiles by their maxima).
# The replicates, the selected lambda and the original predictability are read from the result's fields (API gap).

inputs <- list(
  list(n = 150L, p = 5L, rho = 1.2, seed = 401L, gamma = 0.5),
  list(n = 300L, p = 6L, rho = 1, seed = 402L, gamma = 0.25)
)

test_that("the original partial correlations and the selected lambda equal glasso::glasso() with EBIC", {
  skip_without_reference("glasso")
  invisible(lapply(inputs, \(input) {
    data <- psy_continuous(input$n, input$p, input$rho, input$seed)
    found <- boot_glasso(data, iter = 10L, cs_iter = 2L, gamma = input$gamma, centrality = "strength", seed = 1L)
    cor_matrix <- stats::cor(data)
    ebic <- psy_ebic_path(cor_matrix, input$n, input$gamma, found$lambda_path)
    expect_reference(found$lambda_selected, found$lambda_path[which.min(ebic)], reference = "glasso",
                     covers = "boot_glasso()", quantity = sprintf("selected lambda, n = %d", input$n),
                     tolerance = 1e-12, measure = "relative")
    precision <- glasso::glasso(cor_matrix, rho = found$lambda_selected, penalize.diagonal = FALSE, thr = 1e-12)$wi
    weights <- psy_partial(precision)
    dimnames(weights) <- dimnames(cor_matrix)
    edges <- psy_pairs(weight = weights)
    expected <- data.frame(edge = paste(edges$from, "--", edges$to), weight = edges$weight)
    # Two glasso solvers, each run to its own threshold (agreement about 1e-10).
    expect_reference(summary(found), expected, reference = "glasso", covers = c("boot_glasso()", "summary.boot_glasso()"),
                     quantity = sprintf("original partial correlation, n = %d", input$n), tolerance = 1e-8, key = "edge")
    expect_reference(found$original_predictability, 1 - 1 / diag(precision), reference = "glasso",
                     covers = "boot_glasso()", quantity = sprintf("original predictability R2, n = %d", input$n),
                     tolerance = 1e-8)
  }))
})

test_that("bootstrap intervals and edge inclusion equal the percentile interval and share of non-zero replicates", {
  invisible(lapply(inputs, \(input) {
    data <- psy_continuous(input$n, input$p, input$rho, input$seed)
    found <- boot_glasso(data, iter = 50L, cs_iter = 2L, gamma = input$gamma, centrality = "strength", seed = 2L)
    replicates <- found$boot_edges
    expected <- data.frame(edge = colnames(replicates),
                           ci_lower = apply(replicates, 2L, stats::quantile, probs = 0.025, names = FALSE),
                           ci_upper = apply(replicates, 2L, stats::quantile, probs = 0.975, names = FALSE),
                           inclusion = colMeans(replicates != 0))
    expect_reference(summary(found), expected,
                     reference = "independent: percentile interval (type 7) and inclusion share of the bootstrap replicates",
                     covers = c("boot_glasso()", "summary.boot_glasso()"),
                     quantity = sprintf("bootstrap summary, n = %d", input$n), tolerance = 1e-10, key = "edge")
  }))
})
