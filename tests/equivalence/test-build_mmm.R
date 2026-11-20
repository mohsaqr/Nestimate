# build_mmm(), cluster_mmm(), compare_mmm(), summary.net_mmm(), summary.mmm_compare(), cluster_diagnostics.net_mmm():
# the mixture of first-order Markov chains fitted by EM. EM from random starts is stochastic, so the fit is
# checked by what does not depend on the starts: seqHMM 2.2 evaluates the log-likelihood and the posterior cluster
# probabilities at Nestimate's fitted parameters, seqHMM's own EM started from those parameters stays there (the
# fit is a fixed point of EM), and the information criteria and classification statistics follow from that
# log-likelihood and posterior by their definitions. Re-expressed from local_testing_and_equivalence/
# test-equiv-build_mmm.R (identities by maxima over 100 inputs and a seqHMM parity share) and
# test-mmm-equiv-extracted.R. The parameters, posterior and criteria are read from the result's fields (API gap).

data_sets <- list(two_groups = sq_two_groups(30L, 10L, 42L), regulation = tna::group_regulation[seq_len(120L), ])

test_that("log-likelihood and posterior cluster probabilities equal seqHMM at the fitted parameters", {
  skip_without_reference("seqHMM")
  invisible(lapply(names(data_sets), \(name) {
    data <- data_sets[[name]]
    fit <- build_mmm(data, k = 2, n_starts = 3, seed = 1)
    model <- sq_seqhmm(fit, data)
    expect_reference(fit$log_likelihood, as.numeric(stats::logLik(model)), reference = "seqHMM", covers = "build_mmm()",
                     quantity = sprintf("log-likelihood, %s", name), tolerance = 1e-8, measure = "relative")
    expect_reference(fit$posterior, sq_posterior(model), reference = "seqHMM", covers = "build_mmm()",
                     quantity = sprintf("posterior cluster probability, %s", name), tolerance = 1e-8)
  }))
})

test_that("seqHMM's EM started at Nestimate's unsmoothed fit stays there", {
  skip_without_reference("seqHMM")
  data <- data_sets$two_groups
  fit <- build_mmm(data, k = 2, n_starts = 3, seed = 1, smooth = 0, tol = 1e-10, max_iter = 2000L)
  refit <- seqHMM::fit_model(sq_seqhmm(fit, data), em_step = TRUE, global_step = FALSE, local_step = FALSE,
                             control_em = list(reltol = 1e-12, maxeval = 10000))
  states <- fit$states
  # Both EMs stop by their own relative-change rules (1e-10, 1e-12): the parameters agree to about 1e-9.
  invisible(lapply(1:2, \(cluster) {
    expect_reference(fit$models[[cluster]]$weights[states, states], refit$model$transition_probs[[cluster]],
                     reference = "seqHMM", covers = "build_mmm()",
                     quantity = sprintf("transition probability at the EM fixed point, cluster %d", cluster), tolerance = 1e-6)
    expect_reference(fit$models[[cluster]]$initial[states], refit$model$initial_probs[[cluster]], reference = "seqHMM",
                     covers = "build_mmm()", quantity = sprintf("initial probability at the EM fixed point, cluster %d", cluster),
                     tolerance = 1e-6)
  }))
  expect_reference(fit$log_likelihood, refit$logLik, reference = "seqHMM", covers = "build_mmm()",
                   quantity = "log-likelihood at the EM fixed point", tolerance = 1e-8, measure = "relative")
})

test_that("information criteria, classification statistics and the summaries follow from seqHMM's likelihood", {
  skip_without_reference("seqHMM")
  data <- data_sets$regulation
  fit <- build_mmm(data, k = 2, n_starts = 3, seed = 1)
  model <- sq_seqhmm(fit, data)
  quality <- sq_mixture_quality(as.numeric(stats::logLik(model)), sq_posterior(model), length(fit$states))
  criteria <- data.frame(criterion = c("n_params", "AIC", "BIC", "ICL"),
                         value = c(quality$n_params, quality$AIC, quality$BIC, quality$ICL))
  expect_reference(data.frame(criterion = c("n_params", "AIC", "BIC", "ICL"), value = c(fit$n_params, fit$AIC, fit$BIC, fit$ICL)),
                   criteria, reference = "seqHMM", covers = "build_mmm()", quantity = "information criteria",
                   tolerance = 1e-8, measure = "relative", key = "criterion",
                   note = "AIC, BIC and ICL written from their definitions on seqHMM's log-likelihood and posterior")
  sizes <- as.vector(table(factor(quality$assignments, 1:2)))
  expect_reference(summary(fit), data.frame(component = 1:2, n_assigned = sizes, avepp = quality$avepp),
                   reference = "seqHMM", covers = c("build_mmm()", "summary.net_mmm()"), quantity = "components",
                   tolerance = 1e-8, key = "component")
  expect_reference(as.data.frame(cluster_diagnostics(fit)),
                   data.frame(cluster = 1:2, size = sizes, avepp = quality$avepp, class_err_pct = 100 * (1 - quality$avepp)),
                   reference = "seqHMM", covers = c("cluster_diagnostics()", "cluster_diagnostics.net_mmm()"),
                   quantity = "per-cluster diagnostics", tolerance = 1e-8, key = "cluster")
  # The clustering a mixture's networks carry (attribute "clustering": API gap, no accessor returns it).
  expect_reference(as.data.frame(cluster_diagnostics(attr(build_network(fit), "clustering"))),
                   data.frame(cluster = 1:2, size = sizes, avepp = quality$avepp), reference = "seqHMM",
                   covers = c("cluster_diagnostics()", "cluster_diagnostics.net_mmm_clustering()"),
                   quantity = "per-cluster diagnostics of the networks' clustering", tolerance = 1e-8, key = "cluster")
})

test_that("cluster_mmm() and the compare_mmm() table follow from seqHMM's likelihood", {
  skip_without_reference("seqHMM")
  data <- data_sets$two_groups
  clustered <- cluster_mmm(data, k = 2, n_starts = 3, seed = 1)
  expect_reference(clustered$log_likelihood, as.numeric(stats::logLik(sq_seqhmm(clustered, data))), reference = "seqHMM",
                   covers = "cluster_mmm()", quantity = "log-likelihood, cluster_mmm()", tolerance = 1e-8, measure = "relative")
  compared <- compare_mmm(data, k = 2:3, n_starts = 3, seed = 1, return_fits = TRUE)
  rows <- lapply(attr(compared, "fits"), \(fit) {
    model <- sq_seqhmm(fit, data)
    loglik <- as.numeric(stats::logLik(model))
    quality <- sq_mixture_quality(loglik, sq_posterior(model), length(fit$states))
    data.frame(k = fit$k, log_likelihood = loglik, AIC = quality$AIC, BIC = quality$BIC, ICL = quality$ICL,
               AvePP = quality$AvePP, Entropy = quality$Entropy)
  })
  expected <- do.call(rbind, rows)
  expect_reference(compared, expected, reference = "seqHMM", covers = "compare_mmm()", quantity = "model comparison table",
                   tolerance = 1e-8, measure = "relative", key = "k")
  expect_reference(summary(compared), expected, reference = "seqHMM", covers = c("compare_mmm()", "summary.mmm_compare()"),
                   quantity = "model comparison table, summary()", tolerance = 1e-8, measure = "relative", key = "k")
})
