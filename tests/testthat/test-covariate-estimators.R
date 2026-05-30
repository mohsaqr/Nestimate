# Tests for build_clusters / build_mmm `estimator` argument
# (firth | multinom | chisq).

# ---------------------------------------------------------------------------
# Helpers (small worked examples; no external dependency on real datasets).
# ---------------------------------------------------------------------------

make_seq_data <- function(n_per_cluster = 30, seed = 1) {
  set.seed(seed)
  # Two distinct sequence prototypes -> hierarchical clustering picks them up
  # cleanly, so the partition is stable across runs and tests can rely on it.
  proto1 <- c("A", "A", "B", "B", "C")
  proto2 <- c("C", "C", "B", "A", "A")
  perturb <- function(seq, eps = 0.1) {
    states <- c("A", "B", "C")
    flip <- runif(length(seq)) < eps
    seq[flip] <- sample(states, sum(flip), replace = TRUE)
    seq
  }
  rows1 <- t(replicate(n_per_cluster, perturb(proto1)))
  rows2 <- t(replicate(n_per_cluster, perturb(proto2)))
  df <- as.data.frame(rbind(rows1, rows2), stringsAsFactors = FALSE)
  names(df) <- paste0("T", seq_len(ncol(df)))
  df
}


# ---------------------------------------------------------------------------
# Firth (default): finite ORs even with sparse covariates.
# ---------------------------------------------------------------------------

test_that("estimator = 'firth' is the default and uses brglm2", {
  skip_if_not_installed("brglm2")
  set.seed(42)
  df <- make_seq_data()
  cov <- data.frame(label = c(rep("L1", 59), "rare"))  # one rare cell

  cl <- build_clusters(df, k = 2, dissimilarity = "hamming",
                       covariates = cov)
  expect_equal(cl$covariates$fit$method, "firth")
  # Coefficients must be finite under separation.
  expect_true(all(is.finite(cl$covariates$coefficients$odds_ratio)))
  expect_true(all(is.finite(cl$covariates$coefficients$ci_upper)))
})


test_that("estimator = 'firth' produces tighter ORs than 'multinom' under separation", {
  skip_if_not_installed("brglm2")
  set.seed(7)
  df <- make_seq_data()
  # Construct deliberate separation: the "rare" label only appears in the
  # second half (cluster 2 territory).
  cov <- data.frame(label = c(rep("base", 30), rep("base", 28), rep("rare", 2)))

  expect_warning(
    cl_ml <- build_clusters(df, k = 2, dissimilarity = "hamming",
                             covariates = cov, estimator = "multinom"),
    "estimator = 'multinom'"
  )
  cl_fi <- build_clusters(df, k = 2, dissimilarity = "hamming",
                           covariates = cov, estimator = "firth")

  rare_ml <- cl_ml$covariates$coefficients[
    grepl("rare", cl_ml$covariates$coefficients$variable), "odds_ratio"]
  rare_fi <- cl_fi$covariates$coefficients[
    grepl("rare", cl_fi$covariates$coefficients$variable), "odds_ratio"]

  # ML's OR is artificially huge under separation; Firth's stays bounded.
  # Use 1e6 as a soft cap that any sensible Firth fit on n = 60 sits well below.
  expect_true(any(rare_ml > 1e3),
              info = "ML should blow up on rare-cell separation (sanity check)")
  expect_true(all(rare_fi < 1e3),
              info = "Firth must produce finite, sane ORs on rare-cell data")
})


# ---------------------------------------------------------------------------
# multinom: warns loudly.
# ---------------------------------------------------------------------------

test_that("estimator = 'multinom' emits a separation-risk warning", {
  set.seed(3)
  df <- make_seq_data()
  cov <- data.frame(score = rnorm(60))
  expect_warning(
    build_clusters(df, k = 2, dissimilarity = "hamming",
                   covariates = cov, estimator = "multinom"),
    "diverge under quasi-complete separation"
  )
})


# ---------------------------------------------------------------------------
# chisq: WeightedCluster-style descriptive output.
# ---------------------------------------------------------------------------

test_that("estimator = 'chisq' returns tests + residuals, no ORs", {
  set.seed(11)
  df <- make_seq_data()
  cov <- data.frame(
    factor_var  = factor(c(rep("A", 30), rep("B", 30))),
    numeric_var = rnorm(60)
  )
  cl <- build_clusters(df, k = 2, dissimilarity = "hamming",
                       covariates = cov, estimator = "chisq")
  expect_equal(cl$covariates$fit$method, "chisq")
  expect_null(cl$covariates$coefficients)         # no ORs
  expect_null(cl$covariates$model)                 # no fitted model
  expect_true(is.data.frame(cl$covariates$tests))
  expect_setequal(cl$covariates$tests$variable,
                  c("factor_var", "numeric_var"))
  # factor row carries Cramer's V; numeric row carries eta^2 (KW).
  ftype <- cl$covariates$tests[cl$covariates$tests$type == "factor", ]
  ntype <- cl$covariates$tests[cl$covariates$tests$type == "numeric", ]
  expect_equal(ftype$effect_label, "Cramer's V")
  expect_equal(ntype$effect_label, "eta^2 (KW)")
  # Residuals frame exists for the factor covariate.
  expect_true(is.data.frame(cl$covariates$residuals))
  expect_setequal(cl$covariates$residuals$variable, "factor_var")
})


test_that("estimator = 'chisq' rejects non-factor non-numeric covariates", {
  set.seed(17)
  df <- make_seq_data()
  cov <- data.frame(weird = as.complex(seq(60)))
  expect_error(
    build_clusters(df, k = 2, dissimilarity = "hamming",
                   covariates = cov, estimator = "chisq"),
    "unsupported type"
  )
})


# ---------------------------------------------------------------------------
# build_mmm threading
# ---------------------------------------------------------------------------

test_that("build_mmm threads estimator through the covariate fit", {
  skip_if_not_installed("brglm2")
  set.seed(23)
  df <- make_seq_data()
  cov <- data.frame(score = rnorm(60))
  # Explicit firth: the default estimator = "auto" would pick multinom for
  # a numeric covariate (no rare-cell risk on continuous data), bypassing
  # the firth threading we want to verify here.
  m <- build_mmm(df, k = 2, n_starts = 1, max_iter = 20, seed = 1,
                 covariates = cov, estimator = "firth")
  expect_equal(m$covariates$fit$method, "firth")
})
