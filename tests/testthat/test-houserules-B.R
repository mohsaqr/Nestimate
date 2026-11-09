# House-rule fixes in cluster_data.R, cluster_choice.R, cluster_diagnostics.R,
# mmm.R and mcml_pc.R: classed conditions, corrected printed text, seeds that
# leave the session's random stream alone, and warnings for what used to be
# dropped or left unconverged without a word.
testthat::skip_on_cran()

study_group <- data.frame(
  study_group = rep(c("campus", "online"), length.out = nrow(trajectories)),
  stringsAsFactors = FALSE
)
strategy_clusters <- list(
  metacognitive = c("Planning", "Monitoring", "Evaluating"),
  cognitive     = c("Elaboration", "Organization", "Rehearsal"),
  resource      = c("Help_Seeking", "Time_Mgmt", "Effort_Reg")
)
stream_unmoved <- function(call_with_seed) {
  set.seed(1)
  expected <- runif(1)
  set.seed(1)
  call_with_seed()
  identical(runif(1), expected)
}

# ---- cluster_data.R ---------------------------------------------------------

test_that("the covariate line names the covariates, after chisq and a logit", {
  chisq <- build_clusters(trajectories, k = 2, covariates = study_group,
                          estimator = "chisq")
  logit <- build_clusters(trajectories, k = 2, covariates = study_group,
                          estimator = "multinom") |>
    suppressWarnings(classes = "simpleWarning")
  expect_output(print(chisq), "Covariates: study_group (post-hoc, 1 covariate)",
                fixed = TRUE)
  expect_output(print(logit), "Covariates: study_group (post-hoc, 1 covariate)",
                fixed = TRUE)
})

test_that("the odds-ratio plot of a chisq analysis points at summary()", {
  chisq <- build_clusters(trajectories, k = 2, covariates = study_group,
                          estimator = "chisq")
  err <- expect_error(plot(chisq, type = "predictors"),
                      class = "nestimate_bad_input")
  expect_match(conditionMessage(err), "summary()", fixed = TRUE)
  expect_no_match(conditionMessage(err), "$", fixed = TRUE)
  expect_match(conditionMessage(err), "firth", fixed = TRUE)
})

test_that("cluster_network() without k raises a classed error", {
  expect_error(cluster_network(trajectories), class = "nestimate_bad_input")
})

test_that("build_clusters(seed =) leaves the session's random stream alone", {
  expect_true(stream_unmoved(function() {
    build_clusters(trajectories, k = 2, seed = 5)
  }))
})

test_that("covariate profiles hold one row per cluster, variable and level", {
  fit <- build_clusters(trajectories, k = 2, covariates = study_group,
                        estimator = "chisq")
  profile <- fit$covariates$profiles$categorical
  expect_identical(profile$cluster, c(1L, 1L, 2L, 2L))
  expect_identical(profile$level, rep(c("campus", "online"), 2L))
  campus_1 <- sum(fit$assignments == 1L & study_group$study_group == "campus")
  expect_identical(profile$count[1L], campus_1)
  expect_identical(sum(profile$count), nrow(trajectories))
})

# ---- cluster_choice.R -------------------------------------------------------

test_that("plotting a cluster_choice with nothing swept is a classed error", {
  one <- cluster_choice(trajectories, k = 3)
  expect_error(plot(one), class = "nestimate_bad_input")
  expect_error(plot(one, type = "tradeoff"), class = "nestimate_bad_input")
})

test_that("the sweep keeps one row per configuration, k innermost", {
  sweep <- cluster_choice(trajectories, k = 2:3,
                          dissimilarity = c("hamming", "lcs"))
  expect_identical(sweep$k, c(2L, 3L, 2L, 3L))
  expect_identical(sweep$dissimilarity, rep(c("hamming", "lcs"), each = 2L))
  direct <- build_clusters(trajectories, k = 3, dissimilarity = "lcs",
                           method = "ward.D2")
  expect_equal(sweep$silhouette[4L], direct$silhouette)
  expect_output(print(cluster_choice(trajectories, k = 2:3)), "sweep: k")
})

# ---- cluster_diagnostics.R --------------------------------------------------

test_that("diagnostics label mixture covariates by how they entered the fit", {
  posthoc <- build_mmm(trajectories, k = 2, n_starts = 1, max_iter = 20,
                       seed = 1, covariates = study_group,
                       covariate_effect = "posthoc", estimator = "chisq")
  expect_output(print(cluster_diagnostics(posthoc)),
                "Covariates: study_group (post-hoc, 1 covariate)", fixed = TRUE)
  distance <- build_clusters(trajectories, k = 2, covariates = study_group,
                             estimator = "chisq")
  expect_output(print(cluster_diagnostics(distance)),
                "Covariates: study_group (post-hoc, 1 covariate)", fixed = TRUE)
})

test_that("per-cluster silhouettes are the means of the silhouette widths", {
  fit <- build_clusters(trajectories, k = 3)
  widths <- cluster::silhouette(fit$assignments, dist = fit$distance)
  by_hand <- as.numeric(tapply(widths[, 3L], widths[, 1L], mean))
  expect_equal(cluster_diagnostics(fit)$per_cluster$sil_mean, by_hand)
})

test_that("a group without a clustering names the routes that attach one", {
  plain <- structure(list(), class = "netobject_group")
  err <- expect_error(cluster_diagnostics(plain), class = "nestimate_bad_input")
  expect_match(conditionMessage(err), "build_network() on a build_clusters()",
               fixed = TRUE)
})

# ---- mmm.R ------------------------------------------------------------------

test_that("a chisq covariate analysis prints its covariate and refuses the forest", {
  fit <- build_mmm(trajectories, k = 2, n_starts = 1, max_iter = 20, seed = 1,
                   covariates = study_group, covariate_effect = "posthoc",
                   estimator = "chisq")
  expect_output(print(fit), "Covariates: study_group (post-hoc, 1 covariate)",
                fixed = TRUE)
  expect_error(plot(fit, type = "covariates"), class = "nestimate_bad_input")
})

test_that("build_mmm(seed =) is reproducible and leaves the stream alone", {
  first <- build_mmm(trajectories, k = 2, n_starts = 3, max_iter = 20, seed = 4)
  again <- build_mmm(trajectories, k = 2, n_starts = 3, max_iter = 20, seed = 4)
  expect_identical(first$assignments, again$assignments)
  expect_equal(first$log_likelihood, again$log_likelihood)
  expect_true(stream_unmoved(function() {
    build_mmm(trajectories, k = 2, n_starts = 2, max_iter = 5, seed = 5)
  }))
})

test_that("the BIC and ICL markers of an mmm_compare start in one column", {
  comparison <- structure(
    data.frame(k = 2:3, BIC = c(10, 20), ICL = c(20, 10)),
    class = c("mmm_compare", "data.frame")
  )
  shown <- capture.output(print(comparison))
  expect_identical(regexpr("<-- BIC", shown[grepl("<-- BIC", shown)])[[1L]],
                   regexpr("<-- ICL", shown[grepl("<-- ICL", shown)])[[1L]])
})

# ---- mcml_pc.R --------------------------------------------------------------

test_that("arguments a factor backend does not read are named in a warning", {
  expect_warning(
    build_mcml_pc(srl_strategies, strategy_clusters, method = "cor",
                  weighting = "factor", fa_method = "minres", rotation = 1),
    class = "nestimate_ignored_argument"
  )
})

test_that("unconverged principal axis factoring warns and keeps its loadings", {
  unconverged <- 0L
  warned <- withCallingHandlers(
    build_mcml_pc(srl_strategies, strategy_clusters, method = "cor",
                  weighting = "factor", fa_method = "paf", max_iter = 1),
    nestimate_no_converge = function(w) {
      unconverged <<- unconverged + 1L
      invokeRestart("muffleWarning")
    }
  )
  # One warning per cluster whose extraction stopped at the bound.
  expect_identical(unconverged, length(strategy_clusters))
  quiet <- suppressWarnings(
    build_mcml_pc(srl_strategies, strategy_clusters, method = "cor",
                  weighting = "factor", fa_method = "paf", max_iter = 1),
    classes = "nestimate_no_converge"
  )
  expect_identical(warned$loadings, quiet$loadings)
  expect_false(isTRUE(all.equal(
    warned$loadings$weight,
    suppressWarnings(build_mcml_pc(
      srl_strategies, strategy_clusters, method = "cor",
      weighting = "factor", fa_method = "paf", max_iter = 50
    ), classes = "nestimate_no_converge")$loadings$weight
  )))
})

test_that("loading_stability() has a summary table, one warning, a local seed", {
  fit <- build_mcml_pc(srl_strategies, strategy_clusters, method = "cor",
                       weighting = "factor", fa_method = "paf", max_iter = 1) |>
    suppressWarnings(classes = "nestimate_no_converge")
  warnings_seen <- 0L
  stability <- withCallingHandlers(
    loading_stability(fit, iter = 4, seed = 2),
    nestimate_no_converge = function(w) {
      warnings_seen <<- warnings_seen + 1L
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(warnings_seen, 1L)
  summary_method <- get("summary.pc_loading_stability",
                        envir = asNamespace("Nestimate"))
  expect_identical(summary_method(stability), stability$summary)
  expect_identical(nrow(stability$summary), 9L)

  plain <- build_mcml_pc(srl_strategies, strategy_clusters, method = "cor")
  expect_true(stream_unmoved(function() {
    loading_stability(plain, iter = 3, seed = 5)
  }))
})

test_that("a descriptive fit prints its within estimator and no flipped items", {
  reversed <- srl_strategies
  reversed$Evaluating <- max(reversed$Evaluating) - reversed$Evaluating
  fit <- suppressWarnings(
    build_mcml_pc(reversed, strategy_clusters, aggregation = "average")
  )
  expect_gt(fit$meta$n_flipped, 0L)
  shown <- capture.output(print(fit))
  expect_false(any(grepl("flipped", shown)))
  expect_true(any(grepl("Within: reestimate (pcor)", shown, fixed = TRUE)))

  scored <- suppressWarnings(build_mcml_pc(reversed, strategy_clusters))
  expect_output(print(scored), "reverse-keyed item\\(s\\) flipped")
})

test_that("mcml_pc input errors carry a class", {
  expect_error(build_mcml_pc(srl_strategies, strategy_clusters["cognitive"]),
               class = "nestimate_bad_input")
  expect_error(build_mcml_pc(srl_strategies, strategy_clusters, id_col = "zz"),
               class = "nestimate_bad_input")
  expect_error(as_networks(1), class = "nestimate_bad_input")
  fit <- build_mcml_pc(srl_strategies, strategy_clusters, method = "cor")
  expect_error(loading_stability(fit, iter = 1), class = "nestimate_bad_input")
})
