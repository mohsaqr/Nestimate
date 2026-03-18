# Tests for mmm.R: Mixed Markov Model

# ============================================
# Synthetic data helper
# ============================================

.make_mmm_data <- function(n_per_group = 30, n_cols = 10, seed = 42) {
  set.seed(seed)
  g1 <- data.frame(matrix(sample(c("A", "B", "C"), n_per_group * n_cols,
                                  replace = TRUE, prob = c(0.7, 0.2, 0.1)),
                           ncol = n_cols))
  g2 <- data.frame(matrix(sample(c("A", "B", "C"), n_per_group * n_cols,
                                  replace = TRUE, prob = c(0.1, 0.2, 0.7)),
                           ncol = n_cols))
  rbind(g1, g2)
}

# ============================================
# build_mmm basic
# ============================================

test_that("build_mmm returns correct class and structure", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 3, seed = 1)

  expect_s3_class(mmm, "net_mmm")
  expect_equal(mmm$k, 2)
  expect_equal(mmm$n_sequences, 60)
  expect_equal(length(mmm$states), 3)
  expect_equal(length(mmm$models), 2)
  expect_equal(length(mmm$mixing), 2)
  expect_equal(nrow(mmm$posterior), 60)
  expect_equal(ncol(mmm$posterior), 2)
  expect_equal(length(mmm$assignments), 60)
})

test_that("build_mmm models are netobjects", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 2, seed = 1)

  expect_s3_class(mmm$models[[1]], "netobject")
  expect_s3_class(mmm$models[[2]], "netobject")
  expect_equal(nrow(mmm$models[[1]]$weights), 3)
  expect_equal(ncol(mmm$models[[1]]$weights), 3)
  expect_equal(mmm$models[[1]]$nodes$label, mmm$states)
})

test_that("build_mmm transition matrices are row-normalized", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 2, seed = 1)

  lapply(mmm$models, function(m) {
    rs <- rowSums(m$weights)
    expect_equal(unname(rs), rep(1, 3), tolerance = 1e-10)
  })
})

test_that("build_mmm mixing proportions sum to 1", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 2, seed = 1)

  expect_equal(sum(mmm$mixing), 1, tolerance = 1e-10)
})

test_that("build_mmm posteriors sum to 1 per row", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 2, seed = 1)

  rs <- rowSums(mmm$posterior)
  expect_equal(rs, rep(1, 60), tolerance = 1e-10)
})

test_that("build_mmm assignments are in valid range", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 3, n_starts = 2, seed = 1)

  expect_true(all(mmm$assignments %in% 1:3))
})

test_that("build_mmm recovers two distinct groups", {
  data <- .make_mmm_data(n_per_group = 50, n_cols = 15, seed = 123)
  mmm <- build_mmm(data, k = 2, n_starts = 5, seed = 1)

  diff <- sum(abs(mmm$models[[1]]$weights - mmm$models[[2]]$weights))
  expect_true(diff > 0.5)
  expect_true(all(mmm$mixing > 0.2))
})

test_that("build_mmm is reproducible with seed", {
  data <- .make_mmm_data()
  mmm1 <- build_mmm(data, k = 2, n_starts = 3, seed = 42)
  mmm2 <- build_mmm(data, k = 2, n_starts = 3, seed = 42)

  expect_equal(mmm1$log_likelihood, mmm2$log_likelihood)
  expect_equal(mmm1$mixing, mmm2$mixing)
  expect_equal(mmm1$assignments, mmm2$assignments)
})

test_that("build_mmm computes BIC, AIC, and ICL", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 2, seed = 1)

  expect_true(is.numeric(mmm$BIC))
  expect_true(is.numeric(mmm$AIC))
  expect_true(is.numeric(mmm$ICL))
  expect_true(mmm$BIC > mmm$AIC)
  expect_true(mmm$ICL >= mmm$BIC)  # ICL adds entropy penalty
  expect_true(is.finite(mmm$log_likelihood))
})

test_that("build_mmm converges", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 2, seed = 1)

  expect_true(mmm$converged)
  expect_true(mmm$iterations > 1)
  expect_true(mmm$iterations < 200)
})

# ============================================
# Quality metrics
# ============================================

test_that("quality metrics are computed", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 3, seed = 1)
  q <- mmm$quality

  expect_true(is.list(q))
  expect_equal(length(q$avepp), 2)
  expect_true(all(q$avepp > 0 & q$avepp <= 1))
  expect_true(q$avepp_overall > 0 & q$avepp_overall <= 1)
  expect_true(q$entropy >= 0 & q$entropy <= 1)
  expect_true(q$relative_entropy >= 0 & q$relative_entropy <= 1)
  expect_equal(q$entropy + q$relative_entropy, 1, tolerance = 1e-10)
  expect_true(q$classification_error >= 0 & q$classification_error <= 1)
})

test_that("well-separated data has high AvePP and low entropy", {
  data <- .make_mmm_data(n_per_group = 50, n_cols = 15, seed = 123)
  mmm <- build_mmm(data, k = 2, n_starts = 5, seed = 1)

  expect_true(mmm$quality$avepp_overall > 0.7)
  expect_true(mmm$quality$entropy < 0.5)
})

# ============================================
# Input types
# ============================================

test_that("build_mmm errors on bad k", {
  data <- .make_mmm_data()
  expect_error(build_mmm(data, k = 1), "k")
})

test_that("build_mmm errors on too few sequences", {
  data <- data.frame(X1 = "A", X2 = "B")
  expect_error(build_mmm(data, k = 2), "sequences")
})

test_that("build_mmm works with netobject input", {
  data <- .make_mmm_data()
  net <- build_network(data, method = "relative")
  mmm <- build_mmm(net, k = 2, n_starts = 2, seed = 1)

  expect_s3_class(mmm, "net_mmm")
  expect_equal(mmm$k, 2)
})

test_that("build_mmm works with tna input", {
  skip_if_not_installed("tna")
  model <- tna::tna(tna::group_regulation)
  mmm <- build_mmm(model, k = 2, n_starts = 3, seed = 1)

  expect_s3_class(mmm, "net_mmm")
  expect_equal(mmm$k, 2)
  expect_true(length(mmm$states) > 2)
})

# ============================================
# compare_mmm
# ============================================

test_that("compare_mmm returns comparison table", {
  data <- .make_mmm_data()
  comp <- compare_mmm(data, k = 2:3, n_starts = 2, seed = 1)

  expect_s3_class(comp, "mmm_compare")
  expect_equal(nrow(comp), 2)
  expect_true(all(c("k", "log_likelihood", "AIC", "BIC", "ICL",
                     "AvePP", "Entropy") %in% names(comp)))
})

test_that("compare_mmm BIC favors correct k", {
  data <- .make_mmm_data(n_per_group = 50)
  comp <- compare_mmm(data, k = 2:4, n_starts = 3, seed = 1)

  expect_equal(comp$k[which.min(comp$BIC)], 2)
})

# ============================================
# S3 methods
# ============================================

test_that("print.net_mmm shows quality metrics", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 2, seed = 1)

  expect_output(print(mmm), "Mixed Markov Model")
  expect_output(print(mmm), "AvePP")
  expect_output(print(mmm), "Entropy")
})

test_that("summary.net_mmm shows transition matrices", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 2, seed = 1)

  expect_output(summary(mmm), "Cluster 1")
  expect_output(summary(mmm), "Cluster 2")
})

test_that("print.mmm_compare produces output", {
  data <- .make_mmm_data()
  comp <- compare_mmm(data, k = 2:3, n_starts = 2, seed = 1)

  expect_output(print(comp), "MMM Model Comparison")
  expect_output(print(comp), "BIC")
})

# ============================================
# Covariate-integrated MMM
# ============================================

.make_mmm_cov_data <- function(n = 80, seq_len = 15, seed = 42) {
  set.seed(seed)
  states <- LETTERS[1:3]
  P1 <- matrix(c(0.1,0.7,0.2, 0.3,0.1,0.6, 0.5,0.2,0.3), 3, 3, byrow = TRUE)
  P2 <- matrix(c(0.6,0.2,0.2, 0.1,0.6,0.3, 0.2,0.3,0.5), 3, 3, byrow = TRUE)

  Age <- rnorm(n)
  true_p2 <- 1 / (1 + exp(-1.5 * Age))
  true_cl <- vapply(seq_len(n), function(i) {
    sample(2L, 1L, prob = c(1 - true_p2[i], true_p2[i]))
  }, integer(1L))

  seqs <- matrix(NA_character_, n, seq_len)
  for (i in seq_len(n)) {
    P <- if (true_cl[i] == 1L) P1 else P2
    s <- sample(3L, 1L)
    seqs[i, 1L] <- states[s]
    for (t in 2:seq_len) {
      s <- sample(3L, 1L, prob = P[s, ])
      seqs[i, t] <- states[s]
    }
  }
  df <- as.data.frame(seqs, stringsAsFactors = FALSE)
  colnames(df) <- paste0("T", seq_len(seq_len))
  df$Age <- Age
  list(data = df, true_cl = true_cl, true_beta = 1.5)
}

test_that("build_mmm with covariates returns $covariates", {
  sim <- .make_mmm_cov_data(n = 60)
  mmm <- build_mmm(sim$data, k = 2, n_starts = 3, seed = 42,
                    covariates = "Age")
  expect_s3_class(mmm, "net_mmm")
  expect_true(!is.null(mmm$covariates))
  expect_true(is.data.frame(mmm$covariates$coefficients))
  expect_true(!is.null(mmm$covariates$beta))
  expect_true(!is.null(mmm$covariates$profiles))
  expect_equal(length(mmm$states), 3L)
})

test_that("coefficient direction is recovered", {
  sim <- .make_mmm_cov_data(n = 100, seed = 99)
  mmm <- build_mmm(sim$data, k = 2, n_starts = 5, seed = 99,
                    covariates = "Age")
  # The Age coefficient should be non-zero and match true direction
  # Account for label switching: check that |beta| is positive somewhere
  age_coef <- mmm$covariates$beta[1L, "Age"]
  expect_true(abs(age_coef) > 0.3,
              info = sprintf("Age coef = %.3f, expected |beta| > 0.3", age_coef))
})

test_that("covariate model LL >= plain model LL", {
  sim <- .make_mmm_cov_data(n = 60, seed = 77)
  seq_only <- sim$data[, grep("^T", names(sim$data))]
  mmm_plain <- build_mmm(seq_only, k = 2, n_starts = 3, seed = 77)
  mmm_cov <- build_mmm(sim$data, k = 2, n_starts = 3, seed = 77,
                        covariates = "Age")
  # Nested model: covariate LL should be >= plain LL
  expect_true(mmm_cov$log_likelihood >= mmm_plain$log_likelihood - 0.1,
              info = sprintf("cov LL=%.2f, plain LL=%.2f",
                             mmm_cov$log_likelihood, mmm_plain$log_likelihood))
})

test_that("param count increases with covariates", {
  sim <- .make_mmm_cov_data(n = 60)
  seq_only <- sim$data[, grep("^T", names(sim$data))]
  mmm_plain <- build_mmm(seq_only, k = 2, n_starts = 3, seed = 42)
  mmm_cov <- build_mmm(sim$data, k = 2, n_starts = 3, seed = 42,
                        covariates = "Age")
  # Covariate model has (k-1)*(p+1) mixing params instead of (k-1)
  # With 1 covariate: (2-1)*(1+1) = 2 vs (2-1) = 1 → +1 param
  expect_true(mmm_cov$n_params > mmm_plain$n_params)
})

test_that("print and summary show covariates", {
  sim <- .make_mmm_cov_data(n = 60)
  mmm <- build_mmm(sim$data, k = 2, n_starts = 3, seed = 42,
                    covariates = "Age")
  expect_output(print(mmm), "Covariates:")
  expect_output(print(mmm), "integrated")
  expect_output(summary(mmm), "Covariate Analysis")
  expect_output(summary(mmm), "Predictors of Membership")
})

test_that("plot type='covariates' works", {
  sim <- .make_mmm_cov_data(n = 60)
  mmm <- build_mmm(sim$data, k = 2, n_starts = 3, seed = 42,
                    covariates = "Age")
  p <- plot(mmm, type = "covariates")
  expect_s3_class(p, "ggplot")
})

test_that("plot type='covariates' errors without covariates", {
  data <- .make_mmm_data(n_per_group = 20)
  mmm <- build_mmm(data, k = 2, n_starts = 2, seed = 1)
  expect_error(plot(mmm, type = "covariates"), "No covariate analysis")
})

test_that("covariates with data.frame input form works", {
  sim <- .make_mmm_cov_data(n = 60)
  cov_df <- data.frame(Age = sim$data$Age)
  seq_only <- sim$data[, grep("^T", names(sim$data))]
  mmm <- build_mmm(seq_only, k = 2, n_starts = 3, seed = 42,
                    covariates = cov_df)
  expect_true(!is.null(mmm$covariates))
})

test_that("covariate with NAs warns and works", {
  sim <- .make_mmm_cov_data(n = 60)
  sim$data$Age[c(1, 5, 10)] <- NA
  expect_warning(
    mmm <- build_mmm(sim$data, k = 2, n_starts = 3, seed = 42,
                      covariates = "Age"),
    "Dropped 3 rows"
  )
  expect_true(!is.null(mmm$covariates))
  expect_equal(mmm$n_sequences, 57L)
})

# ============================================
# Coverage gap: .mmm_em called with NULL init_posterior (L29-30)
# ============================================

test_that("build_mmm runs without init_posterior (random init path)", {
  data <- .make_mmm_data()
  # n_starts = 1 → only one init via .mmm_init_kmeans, rest random
  mmm <- build_mmm(data, k = 2, n_starts = 1, seed = 7)
  expect_s3_class(mmm, "net_mmm")
})

# ============================================
# Coverage gap: .mmm_em called with NULL from_ind (L39-41)
# ============================================

test_that("build_mmm handles NULL from_ind path via direct EM call", {
  # build_mmm always computes from_ind; test via a single EM run
  data <- .make_mmm_data(n_per_group = 10, n_cols = 5)
  states <- c("A", "B", "C")
  n_states <- 3L
  counts_raw <- matrix(0L, 20, n_states^2)
  counts_raw[1, 1] <- 1L  # minimal data
  init_st <- rep(1L, 20)
  # Call .mmm_em directly with from_ind = NULL to trigger the computation path
  result <- Nestimate:::.mmm_em(counts_raw, init_st, 2L, 10L, 1e-6, 0.1,
                                 n_states, init_posterior = NULL,
                                 from_ind = NULL, cov_df = NULL)
  expect_true(is.list(result))
  expect_true(!is.null(result$posterior))
})

# ============================================
# Coverage gap: .mmm_init_kmeans when N <= n_comp (L149-151)
# ============================================

test_that(".mmm_init_kmeans handles N <= n_comp edge case", {
  # N = 2, n_comp = 3 → cannot cluster, use round-robin assignment
  counts <- matrix(c(1, 0, 0, 0, 0, 0, 0, 1, 0), nrow = 2, ncol = 9)
  post <- Nestimate:::.mmm_init_kmeans(counts, 3L)
  expect_equal(nrow(post), 2L)
  expect_equal(ncol(post), 3L)
  # Each row should be a valid probability vector
  expect_equal(unname(rowSums(post)), c(1, 1), tolerance = 1e-10)
})

# ============================================
# Coverage gap: .mmm_init_kmeans when kmeans fails (L160-161)
# ============================================

test_that(".mmm_init_kmeans returns random post when kmeans fails", {
  # Pass degenerate all-zero data to force kmeans failure
  counts <- matrix(0, 5, 9)
  post <- Nestimate:::.mmm_init_kmeans(counts, 3L)
  expect_equal(nrow(post), 5L)
  expect_equal(ncol(post), 3L)
  # Should still be valid posterior (rows sum to 1)
  rs <- rowSums(post)
  expect_equal(unname(rs), rep(1, 5), tolerance = 1e-10)
})

# ============================================
# Coverage gap: build_mmm with k=2 covariates Newton step (L266, L270-282)
# ============================================

test_that("build_mmm with k=2 covariate Newton step runs without error", {
  sim <- .make_mmm_cov_data(n = 80, seed = 11)
  mmm <- build_mmm(sim$data, k = 2, n_starts = 2, seed = 11,
                    covariates = "Age")
  expect_s3_class(mmm, "net_mmm")
  expect_true(!is.null(mmm$covariates$beta))
  expect_equal(ncol(mmm$covariates$beta), 2L)  # intercept + Age
})

# ============================================
# Coverage gap: build_mmm with k=3 covariates general Newton step (L396-397)
# ============================================

test_that("build_mmm with k=3 runs general covariate Newton step", {
  sim <- .make_mmm_cov_data(n = 90, seed = 21)
  mmm <- build_mmm(sim$data, k = 3, n_starts = 2, seed = 21,
                    covariates = "Age")
  expect_s3_class(mmm, "net_mmm")
  expect_equal(mmm$k, 3L)
  expect_true(!is.null(mmm$covariates))
})

# ============================================
# Coverage gap: compare_mmm print with matching BIC and ICL (L499-500, L507, L512, L515)
# ============================================

test_that("print.mmm_compare handles case where best_bic == best_icl", {
  data <- .make_mmm_data()
  comp <- compare_mmm(data, k = 2:3, n_starts = 2, seed = 1)
  # Manually force best_bic == best_icl
  comp$BIC[1] <- -1000
  comp$ICL[1] <- -1000
  out <- capture.output(print(comp))
  expect_true(any(grepl("BIC", out)))
})

test_that("print.mmm_compare handles case where best_bic != best_icl", {
  data <- .make_mmm_data(n_per_group = 50)
  comp <- compare_mmm(data, k = 2:4, n_starts = 2, seed = 5)
  # Force mismatch
  comp$BIC[1] <- min(comp$BIC) - 1
  comp$ICL[2] <- min(comp$ICL) - 1
  out <- capture.output(print(comp))
  expect_true(any(grepl("<-- BIC", out)))
  expect_true(any(grepl("<-- ICL", out)))
})

# ============================================
# Coverage gap: plot.net_mmm type='posterior' (L727-728)
# ============================================

test_that("plot.net_mmm type='posterior' returns ggplot invisibly", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 2, seed = 1)
  p <- plot(mmm, type = "posterior")
  expect_s3_class(p, "ggplot")
})

# ============================================
# Coverage gap: plot.net_mmm type='networks' removed
# ============================================

test_that("plot.net_mmm rejects removed 'networks' type", {
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 2, seed = 1)
  expect_error(plot(mmm, type = "networks"), "should be one of")
})

# ============================================
# Coverage gap: plot.mmm_compare returns ggplot (L755-756, L759-773)
# ============================================

test_that("plot.mmm_compare returns ggplot", {
  data <- .make_mmm_data()
  comp <- compare_mmm(data, k = 2:3, n_starts = 2, seed = 1)
  p <- plot(comp)
  expect_s3_class(p, "ggplot")
})

test_that("plot.mmm_compare errors without ggplot2", {
  skip_if(requireNamespace("ggplot2", quietly = TRUE),
          "ggplot2 installed; skipping test")
  data <- .make_mmm_data()
  comp <- compare_mmm(data, k = 2:3, n_starts = 2, seed = 1)
  expect_error(plot(comp), "ggplot2")
})

# ============================================
# Coverage gap: .plot_mmm_posterior requires ggplot2 (L783)
# ============================================

test_that(".plot_mmm_posterior errors without ggplot2", {
  skip_if(requireNamespace("ggplot2", quietly = TRUE),
          "ggplot2 installed; skipping test")
  data <- .make_mmm_data()
  mmm <- build_mmm(data, k = 2, n_starts = 2, seed = 1)
  expect_error(plot(mmm, type = "posterior"), "ggplot2")
})

# ============================================
# Coverage gap: summary.net_mmm with covariates calls .print_covariate_profiles (L790-791)
# ============================================

test_that("summary.net_mmm with covariates prints profiles", {
  sim <- .make_mmm_cov_data(n = 60)
  mmm <- build_mmm(sim$data, k = 2, n_starts = 2, seed = 42,
                    covariates = "Age")
  out <- capture.output(summary(mmm))
  expect_true(any(grepl("Covariate Analysis", out)))
})

# ============================================
# Coverage gap: build_mmm cograph_network input (L794-811 region: decode path)
# ============================================

test_that("build_mmm works with cograph_network (decode path)", {
  data <- .make_mmm_data()
  net <- build_network(data, method = "relative")
  # net inherits both netobject and cograph_network
  mmm <- build_mmm(net, k = 2, n_starts = 2, seed = 1)
  expect_s3_class(mmm, "net_mmm")
  expect_equal(mmm$k, 2L)
})
