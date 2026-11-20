# build_mlvar(), coefs(), summary.net_mlvar(): the multilevel vector autoregression (Epskamp, Waldorp, Mottus and
# Borsboom 2018) with fixed temporal and contemporaneous effects, estimated by lme4. Nestimate delegates the fit to
# idiographic::fit_mlvar() (internal); the reference is mlVAR::mlVAR() with the same settings.
# Re-expressed from local_testing_and_equivalence/test-mlvar-equiv-extracted.R and test-equiv-mlvar.R (five
# simulated panels, compared by their largest difference); here three simulated panels, every cell declared.
# The three networks of a net_mlvar are read as its elements `temporal`, `contemporaneous`, `between` (API gap:
# no accessor returns one layer).

inputs <- list(list(n_id = 20L, days = 3L, beeps = 12L, p = 4L, seed = 201L),
               list(n_id = 30L, days = 2L, beeps = 15L, p = 5L, seed = 205L),
               list(n_id = 25L, days = 4L, beeps = 10L, p = 3L, seed = 210L))

fit_pair <- function(input) {
  data <- psy_panel(input$n_id, input$days, input$beeps, input$p, input$seed)
  vars <- paste0("V", seq_len(input$p))
  list(found = build_mlvar(data, vars = vars, id = "id", day = "day", beep = "beep"),
       reference = mlVAR::mlVAR(data, vars = vars, idvar = "id", dayvar = "day", beepvar = "beep", lags = 1,
                                estimator = "lmer", temporal = "fixed", contemporaneous = "fixed", scale = FALSE,
                                verbose = FALSE))
}

test_that("fixed temporal effects, standard errors and p-values equal mlVAR::mlVAR()", {
  skip_without_reference("mlVAR")
  invisible(lapply(inputs, \(input) {
    pair <- fit_pair(input)
    expected <- psy_mlvar_coefs(pair$reference)
    label <- sprintf("%d people, %d variables", input$n_id, input$p)
    # The same lme4 fits on both sides; relative because the estimates and standard errors span orders of magnitude.
    expect_reference(coefs(pair$found), subset(expected, select = -p), reference = "mlVAR",
                     covers = c("build_mlvar()", "coefs()", "coefs.net_mlvar()"),
                     quantity = sprintf("temporal coefficients, %s", label), tolerance = 1e-10, measure = "relative",
                     key = c("outcome", "predictor"))
    expect_reference(summary(pair$found), subset(expected, select = -p), reference = "mlVAR",
                     covers = c("build_mlvar()", "summary.net_mlvar()"), quantity = sprintf("summary, %s", label),
                     tolerance = 1e-10, measure = "relative", key = c("outcome", "predictor"))
    # mlVAR computes the p-value as 2 * (1 - pnorm(|t|)), whose absolute precision is about 1e-16 (it returns 0
    # where the p-value is 1e-28): its p-values are compared on the absolute scale. A first run compared them
    # relative to 1e-10 and failed only on p-values below 1e-8 (differences of 1e-16 to 1e-19), which is how the
    # cause was found. The relative comparison is made against the normal tail written in the test below.
    expect_reference(coefs(pair$found), subset(expected, select = c(outcome, predictor, p)), reference = "mlVAR",
                     covers = c("build_mlvar()", "coefs()", "coefs.net_mlvar()"),
                     quantity = sprintf("temporal p-values, %s", label), tolerance = 1e-14,
                     key = c("outcome", "predictor"))
    tail <- data.frame(outcome = expected$outcome, predictor = expected$predictor,
                       p = 2 * stats::pnorm(abs(expected$beta / expected$se), lower.tail = FALSE))
    expect_reference(coefs(pair$found), tail, reference = "independent: two-sided normal tail of mlVAR's beta / se, base R",
                     covers = c("build_mlvar()", "coefs()", "coefs.net_mlvar()"),
                     quantity = sprintf("temporal p-values from the normal tail, %s", label), tolerance = 1e-8,
                     measure = "relative", key = c("outcome", "predictor"))
  }))
})

test_that("the temporal, contemporaneous and between-person networks equal mlVAR::mlVAR()", {
  skip_without_reference("mlVAR")
  invisible(lapply(inputs, \(input) {
    pair <- fit_pair(input)
    results <- pair$reference$results
    label <- sprintf("%d people, %d variables", input$n_id, input$p)
    contemporaneous <- results$Theta$pcor$mean
    diag(contemporaneous) <- 0
    between <- results$Omega_mu$pcor$mean
    diag(between) <- 0
    # build_mlvar() stores the temporal network predictor (row) to outcome (column): mlVAR's Beta transposed.
    expect_reference(extract_transition_matrix(pair$found$temporal), t(results$Beta$mean[, , 1L]), reference = "mlVAR",
                     covers = "build_mlvar()", quantity = sprintf("temporal network, %s", label), tolerance = 1e-12)
    expect_reference(extract_transition_matrix(pair$found$contemporaneous), contemporaneous, reference = "mlVAR",
                     covers = "build_mlvar()", quantity = sprintf("contemporaneous network, %s", label), tolerance = 1e-12)
    # The between-person partial correlations pass through one more matrix inversion: 1e-10, the original test's.
    expect_reference(extract_transition_matrix(pair$found$between), between, reference = "mlVAR",
                     covers = "build_mlvar()", quantity = sprintf("between-person network, %s", label), tolerance = 1e-10)
  }))
})
