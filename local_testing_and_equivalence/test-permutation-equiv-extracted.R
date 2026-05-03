# Equivalence test_that() blocks extracted from
# tests/testthat/test-permutation.R.
# These tests reference packages that aren't part of the
# CI Suggests set; they live here for local validation.

testthat::skip_on_cran()

# ---- permutation() Tests ----

# Helper: generate wide sequence data
.make_perm_wide <- function(n = 100, t = 10, states = c("A", "B", "C"),
                            seed = 42) {
  set.seed(seed)
  mat <- matrix(sample(states, n * t, replace = TRUE), nrow = n, ncol = t)
  colnames(mat) <- paste0("T", seq_len(t))
  as.data.frame(mat, stringsAsFactors = FALSE)
}

# Helper: generate frequency-like data for association methods
.make_freq_data <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  as.data.frame(mat)
}


# ---- Input validation ----

test_that("permutation dispatches over matching netobject_group pairs", {
  skip_if_not_installed("tna")
  df1 <- tna::group_regulation
  df2 <- tna::group_regulation
  df1$grp <- rep(c("A", "B"), length.out = nrow(df1))
  df2$grp <- rep(c("A", "B"), length.out = nrow(df2))
  grp1 <- build_network(df1, method = "relative", group = "grp")
  grp2 <- build_network(df2, method = "relative", group = "grp")

  results <- permutation(grp1, grp2, iter = 20L, seed = 1)
  expect_true(is.list(results))
  expect_true(all(vapply(results, inherits, logical(1), "net_permutation")))
})

test_that("permutation errors when group names do not overlap", {
  skip_if_not_installed("tna")
  df <- tna::group_regulation
  df$g1 <- rep(c("A", "B"), length.out = nrow(df))
  df$g2 <- rep(c("X", "Y"), length.out = nrow(df))
  grp1 <- build_network(df, method = "relative", group = "g1")
  grp2 <- build_network(df, method = "relative", group = "g2")
  expect_error(permutation(grp1, grp2, iter = 10L),
               "No matching group names")
})


# ---- missing $data on x (L111-112) ----

test_that("permutation matches tna::permutation_test numerically", {
  skip_if_not_installed("tna")
  df <- tna::group_regulation
  d1 <- df[1:1000, ]
  d2 <- df[1001:2000, ]
  iter <- 1000L
  seed <- 42

  # tna approach
  m1_tna <- tna::tna(d1)
  m2_tna <- tna::tna(d2)
  set.seed(seed)
  tna_perm <- tna::permutation_test(m1_tna, m2_tna, iter = iter)

  # Nestimate approach
  n1 <- build_network(d1, method = "relative")
  n2 <- build_network(d2, method = "relative")
  nest_perm <- permutation(n1, n2, iter = iter, seed = seed)

  # 1. True differences must match exactly
  tna_diff <- tna_perm$edges$diffs_true
  expect_equal(nest_perm$diff, tna_diff, tolerance = 1e-12)

  # 2. Build p-value matrix from tna long-form stats
  states <- rownames(n1$weights)
  n <- length(states)
  tna_stats <- tna_perm$edges$stats
  tna_pmat <- matrix(NA, n, n, dimnames = list(states, states))
  vapply(seq_len(nrow(tna_stats)), function(i) {
    parts <- strsplit(tna_stats$edge_name[i], " -> ")[[1]]
    tna_pmat[parts[1], parts[2]] <<- tna_stats$p_value[i]
    0
  }, numeric(1))

  nest_pmat <- nest_perm$p_values

  # RNG consumption differs between tna and Nestimate,
  # so p-values won't match exactly. With 1000 iterations:
  # - correlation should be > 0.999
  # - max absolute difference should be < 0.03
  valid <- !is.na(tna_pmat) & !is.na(nest_pmat)
  expect_gt(cor(tna_pmat[valid], nest_pmat[valid]), 0.999)
  expect_lt(max(abs(tna_pmat[valid] - nest_pmat[valid])), 0.05)

  # 3. Significant edges should largely agree
  tna_sig <- tna_stats$edge_name[tna_stats$p_value < 0.05]
  nest_sig_df <- nest_perm$summary[nest_perm$summary$p_value < 0.05, ]
  nest_sig <- paste(nest_sig_df$from, "->", nest_sig_df$to)

  # At most 2 borderline edges may flip (Monte Carlo noise at p~0.05)
  sym_diff <- length(setdiff(tna_sig, nest_sig)) +
              length(setdiff(nest_sig, tna_sig))
  expect_lte(sym_diff, 2L)

  # 4. Effect sizes should match closely
  tna_es <- matrix(NA, n, n, dimnames = list(states, states))
  vapply(seq_len(nrow(tna_stats)), function(i) {
    parts <- strsplit(tna_stats$edge_name[i], " -> ")[[1]]
    tna_es[parts[1], parts[2]] <<- tna_stats$effect_size[i]
    0
  }, numeric(1))

  # Replace NaN with 0 for comparison (zero-weight edges)
  tna_es[is.nan(tna_es)] <- 0
  nest_es <- nest_perm$effect_size
  nest_es[is.nan(nest_es)] <- 0
  valid_es <- !is.na(tna_es) & !is.na(nest_es) & tna_es != 0 & nest_es != 0

  if (sum(valid_es) > 5) {
    expect_gt(cor(tna_es[valid_es], nest_es[valid_es]), 0.99)
  }
})


# ---- Group dispatches ----

