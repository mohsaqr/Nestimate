# build_network(method = "mgm"): the mixed graphical model (Haslbeck and Waldorp 2020), nodewise L1-penalised
# regressions with EBIC selection, AND rule and the LW threshold. Nestimate's arithmetic is psychnets'
# (an internal delegation, calling glmnet); the reference is mgm::mgm(). Re-expressed from
# local_testing_and_equivalence/test-equiv-mgm.R, which called the internal estimator on 30 inputs; here the
# exported build_network() on four fixed inputs, Gaussian only and mixed Gaussian/categorical.

# Gaussian columns with a random linear structure; the last `p_c` cut into `k` ordered categories (integers).
psy_mixed <- function(seed, n, p_g, p_c, k = 3L) {
  set.seed(seed)
  p <- p_g + p_c
  draws <- matrix(stats::rnorm(n * p), n, p) %*% (diag(p) + 0.4 * matrix(stats::rnorm(p * p, sd = 0.3), p, p))
  data <- as.data.frame(draws)
  categorical <- p_g + seq_len(p_c)
  data[categorical] <- lapply(data[categorical], \(column) {
    as.integer(cut(column, breaks = stats::quantile(column, probs = seq(0, 1, length.out = k + 1L)),
                   include.lowest = TRUE, labels = FALSE))
  })
  names(data) <- paste0("V", seq_len(p))
  data
}

inputs <- list(
  list(seed = 1L, n = 150L, p_g = 4L, p_c = 0L),
  list(seed = 2L, n = 300L, p_g = 2L, p_c = 2L),
  list(seed = 3L, n = 150L, p_g = 4L, p_c = 3L),
  list(seed = 4L, n = 300L, p_g = 3L, p_c = 2L)
)

test_that("mixed graphical model weights equal mgm::mgm()", {
  skip_without_reference("mgm")
  invisible(lapply(inputs, \(input) {
    data <- psy_mixed(input$seed, input$n, input$p_g, input$p_c)
    type <- c(rep("g", input$p_g), rep("c", input$p_c))
    level <- c(rep(1L, input$p_g), rep(3L, input$p_c))
    network <- build_network(data, method = "mgm",
                             params = list(type = type, level = level, lambdaGam = 0.25, ruleReg = "AND", threshold = "LW"))
    reference <- mgm::mgm(data = as.matrix(data), type = type, level = level, lambdaSel = "EBIC", lambdaGam = 0.25,
                          ruleReg = "AND", threshold = "LW", overparameterize = FALSE, scale = TRUE, pbar = FALSE,
                          signInfo = FALSE, warnings = FALSE)
    expected <- reference$pairwise$wadj
    dimnames(expected) <- list(names(data), names(data))
    expect_reference(extract_transition_matrix(network), expected, reference = "mgm", covers = "build_network()",
                     quantity = sprintf("edge weight, seed %d, %d Gaussian, %d categorical", input$seed, input$p_g, input$p_c),
                     tolerance = 1e-10, note = "build_network(method = \"mgm\")")
    alias <- build_mgm(data, lambdaGam = 0.25, ruleReg = "AND", params = list(type = type, level = level, threshold = "LW"))
    expect_reference(extract_transition_matrix(alias), expected, reference = "mgm", covers = "build_mgm()",
                     quantity = sprintf("edge weight, seed %d, %d Gaussian, %d categorical", input$seed, input$p_g, input$p_c),
                     tolerance = 1e-10, note = "build_mgm()")
  }))
})
