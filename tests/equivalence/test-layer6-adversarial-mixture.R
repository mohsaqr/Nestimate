# Layer 6, adversarial data: build_mmm(), the EM mixture of first-order Markov chains, on the adversarial sequence
# sets. EM from random starts is stochastic, so the fit is checked by what does not depend on the starts: at the
# fitted parameters (mixing, initial and transition probabilities of each cluster) the log-likelihood and the
# posterior are rewritten from the documented model, in base R: a sequence's initial term is its FIRST POSITION
# (none when that position is missing), its transitions are the adjacent observed pairs (a gap breaks a pair).
# Everything that follows from the likelihood and the posterior is compared in full: the information criteria, the
# assignments, every column of summary() and of cluster_diagnostics(), and the quality statistics. Two seeds per set.
# seqHMM 2.2 is the second reference where its model is the same one: it marginalises a missing position (a hidden
# state at a gap), so on sets with a leading or internal gap the two differ by design (a known difference).
# Tolerances, fixed before comparing: log-likelihood 1e-10 relative (a sum of about 100 logs); posterior,
# AvePP and the shares 1e-10; seqHMM 1e-8 relative (its forward algorithm, as test-build_mmm.R).

sets <- adv_sequence_sets()
fitted_sets <- setdiff(names(sets), "one_actor")

# Fit and keep every warning, so that a test can assert which warnings were raised instead of hiding them.
adv_fit_mmm <- function(wide, ...) {
  seen <- character()
  fit <- withCallingHandlers(build_mmm(wide, k = 2L, n_starts = 3L, ...),
                             warning = \(w) {
                               seen <<- c(seen, conditionMessage(w))
                               invokeRestart("muffleWarning")
                             })
  list(fit = fit, warnings = seen)
}

# Log-likelihood and posterior of a mixture of Markov chains at given parameters, from the definition.
adv_mmm_likelihood <- function(fit, wide) {
  states <- fit$states
  m <- as.matrix(wide)
  first <- m[, 1L]
  pairs <- adv_pairs(wide)
  rows <- factor(pairs$row, levels = seq_len(nrow(m)))
  terms <- vapply(seq_len(fit$k), \(k) {
    transition <- fit$models[[k]]$weights[states, states]
    initial <- fit$models[[k]]$initial[states]
    log_transitions <- tapply(log(transition[cbind(pairs$from, pairs$to)]), rows, sum)
    log_transitions[is.na(log_transitions)] <- 0
    log_initial <- ifelse(is.na(first), 0, log(initial[ifelse(is.na(first), states[[1L]], first)]))
    log(fit$mixing[[k]]) + as.vector(log_transitions) + unname(log_initial)
  }, numeric(nrow(m)))
  top <- apply(terms, 1L, max)
  per_sequence <- top + log(rowSums(exp(terms - top)))
  list(loglik = sum(per_sequence), posterior = exp(terms - per_sequence))
}

# Everything that follows from the likelihood and the posterior, by its definition.
adv_mmm_expected <- function(fit, wide) {
  likelihood <- adv_mmm_likelihood(fit, wide)
  posterior <- likelihood$posterior
  quality <- sq_mixture_quality(likelihood$loglik, posterior, length(fit$states))
  k <- fit$k
  n <- nrow(posterior)
  sizes <- tabulate(quality$assignments, nbins = k)
  modal <- posterior[cbind(seq_len(n), quality$assignments)]
  list(
    likelihood = likelihood,
    criteria = data.frame(criterion = c("log_likelihood", "n_params", "AIC", "BIC", "ICL"),
                          value = c(likelihood$loglik, quality$n_params, quality$AIC, quality$BIC, quality$ICL)),
    assignments = quality$assignments,
    summary = data.frame(component = seq_len(k), prior = fit$mixing, n_assigned = sizes,
                         mean_posterior = quality$avepp, avepp = quality$avepp),
    diagnostics = data.frame(cluster = seq_len(k), size = sizes, pct = 100 * sizes / n, mix_pct = 100 * fit$mixing,
                             avepp = quality$avepp, class_err_pct = 100 * (1 - quality$avepp)),
    quality = c(avepp_overall = quality$AvePP, entropy = quality$Entropy, relative_entropy = 1 - quality$Entropy,
                classification_error = mean(1 - modal), class_entropy = -sum(log(modal))))
}

test_that("at the fitted parameters, the likelihood, posterior and everything derived equal the definition", {
  invisible(lapply(fitted_sets, \(name) {
    wide <- sets[[name]]
    invisible(lapply(c(1L, 2L), \(seed) {
      run <- adv_fit_mmm(wide, seed = seed)
      fit <- run$fit
      # the only warning build_mmm() may raise here is the empty-cluster one (identical rows, a deterministic cycle)
      expect_true(all(grepl("empty hard-assignment cluster", run$warnings)), info = paste(name, toString(run$warnings)))
      expected <- adv_mmm_expected(fit, wide)
      label <- sprintf("%s, seed %d", name, seed)
      expect_reference(data.frame(criterion = c("log_likelihood", "n_params", "AIC", "BIC", "ICL"),
                                  value = c(fit$log_likelihood, fit$n_params, fit$AIC, fit$BIC, fit$ICL)),
                       expected$criteria, reference = "independent: mixture likelihood at the fitted parameters, base R",
                       covers = "build_mmm()", quantity = sprintf("log-likelihood and criteria, %s", label),
                       tolerance = 1e-10, measure = "relative", key = "criterion")
      expect_reference(fit$posterior, expected$likelihood$posterior,
                       reference = "independent: mixture likelihood at the fitted parameters, base R",
                       covers = "build_mmm()", quantity = sprintf("posterior, %s", label), tolerance = 1e-10)
      expect_reference(fit$assignments, expected$assignments,
                       reference = "independent: modal posterior cluster, base R", covers = "build_mmm()",
                       quantity = sprintf("assignments, %s", label), tolerance = 0)
      expect_reference(summary(fit), expected$summary,
                       reference = "independent: mixture likelihood at the fitted parameters, base R",
                       covers = c("build_mmm()", "summary.net_mmm()"), quantity = sprintf("summary(), %s", label),
                       tolerance = 1e-10, key = "component")
      expect_reference(as.data.frame(cluster_diagnostics(fit)), expected$diagnostics,
                       reference = "independent: mixture likelihood at the fitted parameters, base R",
                       covers = c("cluster_diagnostics()", "cluster_diagnostics.net_mmm()"),
                       quantity = sprintf("cluster diagnostics, %s", label), tolerance = 1e-10, key = "cluster")
      expect_reference(unlist(fit$quality[names(expected$quality)]), expected$quality,
                       reference = "independent: mixture likelihood at the fitted parameters, base R",
                       covers = "build_mmm()", quantity = sprintf("quality statistics, %s", label), tolerance = 1e-10)
    }))
  }))
})

test_that("the log-likelihood and posterior equal seqHMM at the fitted parameters where the two models agree", {
  skip_without_reference(c("seqHMM", "TraMineR"))
  # the row of empty_row with no observed state is no TraMineR sequence; the hand likelihood above covers it
  gapless <- setdiff(fitted_sets, c("leading_na", "internal_na", "empty_row"))
  invisible(lapply(gapless, \(name) {
    wide <- sets[[name]]
    fit <- adv_fit_mmm(wide, seed = 1L)$fit
    model <- sq_seqhmm(fit, wide)
    expect_reference(fit$log_likelihood, as.numeric(stats::logLik(model)), reference = "seqHMM", covers = "build_mmm()",
                     quantity = sprintf("log-likelihood, %s", name), tolerance = 1e-8, measure = "relative")
    expect_reference(fit$posterior, sq_posterior(model), reference = "seqHMM", covers = "build_mmm()",
                     quantity = sprintf("posterior, %s", name), tolerance = 1e-8)
  }))
  invisible(lapply(c("leading_na", "internal_na"), \(name) {
    wide <- sets[[name]]
    fit <- adv_fit_mmm(wide, seed = 1L)$fit
    expect_known_difference(fit$log_likelihood, as.numeric(stats::logLik(sq_seqhmm(fit, wide))), reference = "seqHMM",
                            covers = "build_mmm()", quantity = sprintf("log-likelihood, %s", name), tolerance = 1e-8,
                            measure = "relative",
                            reason = paste("seqHMM sums over the hidden state at a missing position (the pair across a",
                                           "gap becomes a two-step transition; a missing first position takes the",
                                           "initial distribution and a transition into the first observed state);",
                                           "build_mmm() drops the pair across a gap and the initial term of a row",
                                           "whose first position is missing, as its documentation states."),
                            agrees_with = "the likelihood of adjacent observed pairs, compared above")
  }))
})

test_that("seqHMM's EM started at an unsmoothed fit stays there", {
  skip_without_reference(c("seqHMM", "TraMineR"))
  invisible(lapply(c("two_states", "numeric_labels"), \(name) {
    wide <- sets[[name]]
    fit <- adv_fit_mmm(wide, seed = 1L, smooth = 0, tol = 1e-10, max_iter = 2000L)$fit
    refit <- seqHMM::fit_model(sq_seqhmm(fit, wide), em_step = TRUE, global_step = FALSE, local_step = FALSE,
                               control_em = list(reltol = 1e-12, maxeval = 10000))
    expect_reference(fit$log_likelihood, refit$logLik, reference = "seqHMM", covers = "build_mmm()",
                     quantity = sprintf("log-likelihood at the EM fixed point, %s", name), tolerance = 1e-8,
                     measure = "relative")
  }))
})

test_that("one actor is refused with a classed error", {
  expect_error(build_mmm(sets$one_actor, k = 2L, n_starts = 3L, seed = 1L), class = "nestimate_bad_input")
})
