# Layer 6, adversarial data: the association networks (cor, pcor, glasso, ising) on a constant column, two
# near-collinear columns (correlation 1 - 1e-8), more variables than rows (8 x 12), binary columns, and binary
# columns with a state seen in 3 of 200 rows. Every field of the result that carries numbers is compared: the
# weights, the correlation matrix, the precision matrix, n, p, the selected lambda and the predictability table.
# A verb that refuses an input is asserted by the class of its condition.
# References: stats::cor(), corpcor::cor2pcor(), the inverse correlation matrix, glasso::glasso() at the selected
# lambda (thr = 1e-12), the extended BIC written from Foygel and Drton (2010) (helper-psychometric.R),
# IsingFit::IsingFit(), stats::lm() for predictability.
# Tolerances, fixed before comparing: correlations 1e-10; partial correlations of a well-conditioned matrix 1e-10;
# of the near-collinear matrix (condition number about 1e8) 1e-6, the rounding of a solve() at that condition
# number; two glasso solvers each to its own threshold 1e-8 (as test-build_glasso.R); IsingFit 1e-10.

sets <- adv_association_sets()

# The columns with a non-zero variance: a constant column is dropped by the association estimators.
adv_varying <- function(data) data[vapply(data, \(column) stats::var(column) > 0, logical(1))]

adv_zero_diag <- function(m) {
  diag(m) <- 0
  m
}

# R2 and RMSE of the least-squares regression of each column on all the others (an intercept included).
adv_ols <- function(data) {
  rows <- lapply(names(data), \(node) {
    fit <- stats::lm(stats::reformulate(setdiff(names(data), node), response = node), data = data)
    data.frame(node = node, R2 = summary(fit)$r.squared, RMSE = sqrt(mean(stats::residuals(fit)^2)))
  })
  do.call(rbind, rows)
}

test_that("a constant column is dropped with a message and the rest equal stats::cor() and the inverse", {
  data <- sets$constant
  kept <- adv_varying(data)
  expect_message(net <- build_network(data, method = "cor"), "zero-variance")
  expect_reference(extract_transition_matrix(net), adv_zero_diag(stats::cor(kept)), reference = "stats",
                   covers = "build_network()", quantity = "correlation weights, constant column", tolerance = 1e-10)
  expect_reference(net$cor_matrix, stats::cor(kept), reference = "stats", covers = "build_network()",
                   quantity = "correlation matrix, constant column", tolerance = 1e-10)
  expect_reference(c(n = net$n, p = net$p), c(n = nrow(kept), p = ncol(kept)),
                   reference = "independent: rows and non-constant columns of the input", covers = "build_network()",
                   quantity = "n and p, constant column", tolerance = 0)
  expect_message(partial <- build_network(data, method = "pcor"), "zero-variance")
  expect_reference(extract_transition_matrix(partial), psy_partial(solve(stats::cor(kept))),
                   reference = "independent: -K_ij / sqrt(K_ii K_jj) of solve(cor), base R", covers = "build_network()",
                   quantity = "partial correlation weights, constant column", tolerance = 1e-10)
  expect_reference(partial$precision_matrix, solve(stats::cor(kept)), reference = "stats", covers = "build_network()",
                   quantity = "precision matrix, constant column", tolerance = 1e-10)
  expect_reference(predictability(partial, data = kept), adv_ols(kept), reference = "stats",
                   covers = c("predictability()", "predictability.netobject()"),
                   quantity = "predictability, constant column", tolerance = 1e-8, key = "node")
})

test_that("near-collinear columns: correlations and partial correlations equal stats and corpcor", {
  skip_without_reference("corpcor")
  data <- sets$near_collinear
  net <- build_network(data, method = "cor")
  expect_reference(extract_transition_matrix(net), adv_zero_diag(stats::cor(data)), reference = "stats",
                   covers = "build_network()", quantity = "correlation weights, near-collinear", tolerance = 1e-10)
  partial <- build_network(data, method = "pcor")
  expected <- corpcor::cor2pcor(stats::cor(data))
  dimnames(expected) <- list(names(data), names(data))
  expect_reference(extract_transition_matrix(partial), adv_zero_diag(expected), reference = "corpcor",
                   covers = "build_network()", quantity = "partial correlation weights, near-collinear",
                   tolerance = 1e-6)
  expect_reference(extract_transition_matrix(partial), psy_partial(solve(stats::cor(data))),
                   reference = "independent: -K_ij / sqrt(K_ii K_jj) of solve(cor), base R", covers = "build_network()",
                   quantity = "partial correlation weights, near-collinear, solve()", tolerance = 1e-6)
  expect_reference(predictability(partial, data = data), adv_ols(data), reference = "stats",
                   covers = c("predictability()", "predictability.netobject()"),
                   quantity = "predictability, near-collinear", tolerance = 1e-6, key = "node")
})

test_that("more variables than rows: cor is computed and warns on predictability, pcor is refused", {
  data <- sets$n_less_p
  expect_warning(net <- build_network(data, method = "cor"), class = "nestimate_singular_neighbours")
  expect_reference(extract_transition_matrix(net), adv_zero_diag(stats::cor(data)), reference = "stats",
                   covers = "build_network()", quantity = "correlation weights, n < p", tolerance = 1e-10)
  expect_error(build_network(data, method = "pcor"), class = "nestimate_bad_input")
})

test_that("the graphical lasso at the selected lambda equals glasso::glasso() and the EBIC selection, every set", {
  skip_without_reference("glasso")
  invisible(lapply(c("constant", "near_collinear", "n_less_p", "binary", "binary_rare"), \(name) {
    data <- adv_varying(sets[[name]])
    net <- if (name == "constant") {
      expect_message(fit <- build_network(sets[[name]], method = "glasso"), "zero-variance")
      fit
    } else build_network(data, method = "glasso")
    correlation <- stats::cor(data)
    fit <- glasso::glasso(correlation, rho = net$lambda_selected, penalize.diagonal = FALSE, thr = 1e-12)
    expected <- psy_partial(fit$wi)
    dimnames(expected) <- dimnames(correlation)
    expect_reference(extract_transition_matrix(net), expected, reference = "glasso", covers = "build_network()",
                     quantity = sprintf("glasso partial correlations, %s", name), tolerance = 1e-8)
    expect_reference(net$cor_matrix, correlation, reference = "stats", covers = "build_network()",
                     quantity = sprintf("glasso correlation matrix, %s", name), tolerance = 1e-10)
    ebic <- psy_ebic_path(correlation, nrow(data), gamma = 0.5, lambdas = net$lambda_path)
    expect_reference(net$lambda_selected, net$lambda_path[which.min(ebic)], reference = "glasso",
                     covers = "build_network()", quantity = sprintf("selected lambda, %s", name), tolerance = 1e-12,
                     measure = "relative")
    expect_reference(c(n = net$n, p = net$p), c(n = nrow(data), p = ncol(data)),
                     reference = "independent: rows and non-constant columns of the input", covers = "build_network()",
                     quantity = sprintf("glasso n and p, %s", name), tolerance = 0)
  }))
})

test_that("the Ising model equals IsingFit::IsingFit() on binary columns, a rare state included", {
  skip_without_reference("IsingFit")
  invisible(lapply(c("binary", "binary_rare"), \(name) {
    data <- sets[[name]]
    # glmnet warns, for the rare column, that a class has fewer than 8 observations; both sides call glmnet
    net <- withCallingHandlers(build_network(data, method = "ising"),
                               warning = \(w) if (grepl("fewer than 8", conditionMessage(w))) invokeRestart("muffleWarning"))
    reference <- withCallingHandlers(IsingFit::IsingFit(data, plot = FALSE, progressbar = FALSE),
                                     warning = \(w) if (grepl("fewer than 8", conditionMessage(w))) invokeRestart("muffleWarning"))
    expect_reference(extract_transition_matrix(net), reference$weiadj, reference = "IsingFit", covers = "build_network()",
                     quantity = sprintf("Ising weights, %s", name), tolerance = 1e-10)
    expect_reference(net$thresholds, reference$thresholds, reference = "IsingFit", covers = "build_network()",
                     quantity = sprintf("Ising thresholds, %s", name), tolerance = 1e-10)
  }))
})

test_that("the Ising model refuses non-binary data with a classed error", {
  invisible(lapply(c("near_collinear", "n_less_p"), \(name) {
    expect_error(build_network(sets[[name]], method = "ising"), class = "nestimate_bad_input")
  }))
})
