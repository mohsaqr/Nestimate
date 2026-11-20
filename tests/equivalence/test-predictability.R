# predictability() (Haslbeck and Waldorp 2018): the share of each node's variance explained by its neighbours
# (R2) and the root mean squared error of that prediction on the data's scale. References: ordinary least
# squares with stats::lm() on the node's neighbours (correlation and partial correlation networks), and the
# glasso precision matrix from glasso::glasso() (graphical lasso networks: R2 = 1 - 1 / K_jj, the regression
# coefficients -K_ij / K_jj). Re-expressed from local_testing_and_equivalence/test-equiv-predictability.R (50 + 50
# random inputs compared by their largest difference); here fixed inputs, every node declared.
# A grouped or multilevel result is a list of tables, one per group or level, read with `$` (API gap: no
# tidy accessor stacks them).

# R2 and RMSE of the least-squares regression of each column on all the others (an intercept included).
psy_ols <- function(data) {
  rows <- lapply(names(data), \(node) {
    fit <- stats::lm(stats::reformulate(setdiff(names(data), node), response = node), data = data)
    data.frame(node = node, R2 = summary(fit)$r.squared, RMSE = sqrt(mean(stats::residuals(fit)^2)))
  })
  do.call(rbind, rows)
}

# R2 and RMSE from a precision matrix K on the data's scale.
psy_precision <- function(data, precision) {
  matrix_x <- as.matrix(data)
  z <- scale(matrix_x)
  rows <- lapply(seq_len(ncol(matrix_x)), \(j) {
    beta <- -precision[-j, j] / precision[j, j]
    fitted <- mean(matrix_x[, j]) + stats::sd(matrix_x[, j]) * (z[, -j, drop = FALSE] %*% beta)
    data.frame(node = colnames(matrix_x)[j], R2 = 1 - 1 / precision[j, j], RMSE = sqrt(mean((matrix_x[, j] - fitted)^2)))
  })
  do.call(rbind, rows)
}

inputs <- list(
  list(n = 80L, p = 5L, rho = 1, seed = 301L),
  list(n = 150L, p = 6L, rho = 1.5, seed = 302L),
  list(n = 300L, p = 7L, rho = 0.8, seed = 303L)
)

test_that("R2 and RMSE of correlation and partial correlation networks equal stats::lm()", {
  invisible(lapply(inputs, \(input) {
    data <- psy_continuous(input$n, input$p, input$rho, input$seed)
    expected <- psy_ols(data)
    expect_reference(predictability(build_cor(data), data = data), expected, reference = "stats",
                     covers = c("predictability()", "predictability.netobject()"),
                     quantity = sprintf("correlation network, n = %d, p = %d", input$n, input$p), tolerance = 1e-8,
                     key = "node")
    expect_reference(predictability(build_pcor(data), data = data), expected, reference = "stats",
                     covers = c("predictability()", "predictability.netobject()"),
                     quantity = sprintf("partial correlation network, n = %d, p = %d", input$n, input$p),
                     tolerance = 1e-8, key = "node")
  }))
})

test_that("R2 and RMSE of graphical lasso networks equal those of the glasso::glasso() precision matrix", {
  skip_without_reference("glasso")
  invisible(lapply(inputs, \(input) {
    data <- psy_continuous(input$n, input$p, input$rho, input$seed)
    network <- build_glasso(data)
    precision <- glasso::glasso(stats::cor(data), rho = network$lambda_selected, penalize.diagonal = FALSE,
                                thr = 1e-12)$wi
    # Two glasso solvers, each run to its own threshold (agreement about 1e-10 on the partial correlations).
    expect_reference(predictability(network, data = data), psy_precision(data, precision), reference = "glasso",
                     covers = c("predictability()", "predictability.netobject()"),
                     quantity = sprintf("graphical lasso network, n = %d, p = %d", input$n, input$p), tolerance = 1e-8,
                     key = "node")
  }))
})

test_that("predictability of grouped networks equals stats::lm() per group", {
  data <- psy_continuous(240L, 5L, 1, 304L)
  data$group <- rep(c("a", "b"), each = 120L)
  found <- predictability(build_network(data, method = "pcor", group = "group"))
  invisible(lapply(c("a", "b"), \(label) {
    rows <- subset(data, group == label, select = -group)
    expect_reference(found[[label]], psy_ols(rows), reference = "stats",
                     covers = c("predictability()", "predictability.netobject_group()"),
                     quantity = sprintf("partial correlation network, group %s", label), tolerance = 1e-8, key = "node")
  }))
})

test_that("R2 of a multilevel network equals stats::lm() on person means and on person-centred scores", {
  data <- psy_continuous(300L, 4L, 1, 305L)
  data$id <- rep(seq_len(30L), each = 10L)
  found <- predictability(build_network(data, method = "pcor", actor = "id", level = "both"))
  variables <- paste0("V", seq_len(4L))
  means <- stats::aggregate(data[variables], by = list(id = data$id), FUN = mean)
  centred <- as.data.frame(lapply(data[variables], \(column) column - stats::ave(column, data$id)))
  expect_reference(found$between, subset(psy_ols(means[variables]), select = c(node, R2)), reference = "stats",
                   covers = c("predictability()", "predictability.netobject_ml()"),
                   quantity = "between persons (person means)", tolerance = 1e-8, key = "node")
  expect_reference(found$within, subset(psy_ols(centred), select = c(node, R2)), reference = "stats",
                   covers = c("predictability()", "predictability.netobject_ml()"),
                   quantity = "within persons (person-centred)", tolerance = 1e-8, key = "node")
})
