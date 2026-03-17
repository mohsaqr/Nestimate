# ===========================================================================
# Tests for bootstrap_mcml() — MCML Bootstrap
# ===========================================================================

# --- Helpers ---
.make_boot_wide <- function(seed = 1) {
  set.seed(seed)
  states <- c("plan", "monitor", "adapt", "discuss", "evaluate", "reflect")
  as.data.frame(matrix(
    sample(states, 500, replace = TRUE), nrow = 100, ncol = 5,
    dimnames = list(NULL, paste0("T", 1:5))
  ))
}

.boot_clusters <- list(
  Regulation = c("plan", "adapt"),
  Cognition = c("monitor", "evaluate"),
  Social = c("discuss", "reflect")
)


# ===========================================================================
# Section 1: Input validation
# ===========================================================================
test_that("bootstrap_mcml rejects matrix input", {
  mat <- matrix(1, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_error(
    bootstrap_mcml(mat, clusters = list(G1 = c("A", "B"), G2 = "C")),
    "data.frame"
  )
})

test_that("bootstrap_mcml rejects single cluster", {
  wide <- .make_boot_wide()
  expect_error(
    bootstrap_mcml(wide, clusters = list(
      All = c("plan", "monitor", "adapt", "discuss", "evaluate", "reflect")
    ))
  )
})

test_that("bootstrap_mcml rejects unnamed clusters", {
  wide <- .make_boot_wide()
  expect_error(
    bootstrap_mcml(wide, clusters = list(
      c("plan", "adapt"), c("monitor", "evaluate")
    )),
    "named list"
  )
})


# ===========================================================================
# Section 2: Structure
# ===========================================================================
test_that("bootstrap_mcml returns correct class", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  expect_s3_class(boot, "mcml_bootstrap")
})

test_that("bootstrap_mcml between is saqr_bootstrap", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  expect_s3_class(boot$between, "saqr_bootstrap")
})

test_that("bootstrap_mcml within is named list of saqr_bootstrap", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  expect_true(is.list(boot$within))
  expect_true(length(boot$within) > 0)
  for (w in boot$within) {
    expect_s3_class(w, "saqr_bootstrap")
  }
})

test_that("bootstrap_mcml original is mcml_network", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  expect_s3_class(boot$original, "mcml_network")
})


# ===========================================================================
# Section 3: Between-cluster bootstrap
# ===========================================================================
test_that("between bootstrap has cluster-level nodes", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  between_nodes <- boot$between$original$nodes
  expect_true(setequal(between_nodes, names(.boot_clusters)))
})

test_that("between bootstrap p-values are in [0, 1]", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  pvals <- boot$between$p_values
  expect_true(all(pvals >= 0 & pvals <= 1))
})

test_that("between bootstrap has correct dimensions", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  n_cl <- length(.boot_clusters)
  expect_equal(nrow(boot$between$mean), n_cl)
  expect_equal(ncol(boot$between$mean), n_cl)
})


# ===========================================================================
# Section 4: Within-cluster bootstrap
# ===========================================================================
test_that("within bootstrap per-cluster has correct states", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  for (cl_name in names(boot$within)) {
    states <- boot$within[[cl_name]]$original$nodes
    expect_true(all(states %in% .boot_clusters[[cl_name]]))
  }
})

test_that("within bootstrap has correct matrix dimensions", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  for (cl_name in names(boot$within)) {
    n_states <- length(.boot_clusters[[cl_name]])
    expect_equal(nrow(boot$within[[cl_name]]$mean), n_states)
    expect_equal(ncol(boot$within[[cl_name]]$mean), n_states)
  }
})

test_that("within bootstrap p-values in [0, 1]", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  for (cl_name in names(boot$within)) {
    pvals <- boot$within[[cl_name]]$p_values
    expect_true(all(pvals >= 0 & pvals <= 1))
  }
})


# ===========================================================================
# Section 5: Reproducibility
# ===========================================================================
test_that("bootstrap_mcml is reproducible with same seed", {
  wide <- .make_boot_wide()
  boot1 <- bootstrap_mcml(wide, clusters = .boot_clusters,
                            iter = 20, seed = 123)
  boot2 <- bootstrap_mcml(wide, clusters = .boot_clusters,
                            iter = 20, seed = 123)
  expect_equal(boot1$between$mean, boot2$between$mean)
  expect_equal(boot1$between$p_values, boot2$between$p_values)
})


# ===========================================================================
# Section 6: S3 methods
# ===========================================================================
test_that("print.mcml_bootstrap produces output", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  out <- capture.output(print(boot))
  expect_true(any(grepl("MCML Bootstrap", out)))
  expect_true(any(grepl("Between", out)))
  expect_true(any(grepl("Within", out)))
})

test_that("summary.mcml_bootstrap returns between summary", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  s <- summary(boot, level = "between")
  expect_true(is.data.frame(s))
})

test_that("summary.mcml_bootstrap returns within summary", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  s <- summary(boot, level = "within")
  expect_true(is.list(s))
})

test_that("plot.mcml_bootstrap runs without error for between", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  expect_no_error(plot(boot, type = "between"))
})

test_that("plot.mcml_bootstrap runs without error for within", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 20, seed = 42)
  expect_no_error(suppressMessages(plot(boot, type = "within")))
})


# ===========================================================================
# Section 7: Single-node cluster handling
# ===========================================================================
test_that("single-node cluster is skipped in within bootstrap", {
  wide <- .make_boot_wide()
  clusters_with_single <- list(
    Regulation = c("plan", "adapt"),
    Cognition = c("monitor", "evaluate"),
    Solo = "discuss",
    Social = c("reflect")
  )
  boot <- bootstrap_mcml(wide, clusters = clusters_with_single,
                          iter = 20, seed = 42)
  # Single-node clusters should not appear in within
  expect_false("Solo" %in% names(boot$within))
  expect_false("Social" %in% names(boot$within))
})

test_that("single-node clusters still in between bootstrap", {
  wide <- .make_boot_wide()
  clusters_with_single <- list(
    Regulation = c("plan", "adapt", "monitor"),
    Cognition = c("evaluate", "reflect"),
    Solo = "discuss"
  )
  boot <- bootstrap_mcml(wide, clusters = clusters_with_single,
                          iter = 20, seed = 42)
  # Between should include all clusters
  expect_true("Solo" %in% boot$between$original$nodes)
})


# ===========================================================================
# Section 8: Integration with group_regulation
# ===========================================================================
test_that("bootstrap_mcml works with group_regulation data", {
  skip_if_not_installed("tna")
  gr <- tna::group_regulation

  clusters <- list(
    Social = c("discuss", "consensus", "cohesion"),
    Cognitive = c("plan", "adapt", "monitor"),
    Meta = c("emotion", "coregulate", "synthesis")
  )

  boot <- bootstrap_mcml(gr, clusters = clusters,
                          iter = 30, seed = 42)

  expect_s3_class(boot, "mcml_bootstrap")
  expect_equal(length(boot$within), 3)
  expect_s3_class(boot$between, "saqr_bootstrap")

  # Between has 3 cluster-level nodes
  expect_equal(boot$between$original$n_nodes, 3)
  expect_true(setequal(boot$between$original$nodes, names(clusters)))
})

test_that("bootstrap_mcml integration: within has correct nodes", {
  skip_if_not_installed("tna")
  gr <- tna::group_regulation

  clusters <- list(
    Social = c("discuss", "consensus", "cohesion"),
    Cognitive = c("plan", "adapt", "monitor"),
    Meta = c("emotion", "coregulate", "synthesis")
  )

  boot <- bootstrap_mcml(gr, clusters = clusters,
                          iter = 30, seed = 42)

  for (cl_name in names(boot$within)) {
    nodes <- boot$within[[cl_name]]$original$nodes
    expect_true(all(nodes %in% clusters[[cl_name]]))
  }
})


# ===========================================================================
# Section 9: Cluster format support in bootstrap
# ===========================================================================
test_that("bootstrap_mcml works with named vector clusters", {
  wide <- .make_boot_wide()
  named_vec <- c(
    plan = "Regulation", adapt = "Regulation",
    monitor = "Cognition", evaluate = "Cognition",
    discuss = "Social", reflect = "Social"
  )
  boot <- bootstrap_mcml(wide, clusters = named_vec,
                          iter = 20, seed = 42)
  expect_s3_class(boot, "mcml_bootstrap")
  expect_equal(boot$between$original$n_nodes, 3)
})

test_that("bootstrap_mcml works with data.frame clusters", {
  wide <- .make_boot_wide()
  cl_df <- data.frame(
    node = c("plan", "adapt", "monitor", "evaluate", "discuss", "reflect"),
    group = c("Regulation", "Regulation", "Cognition", "Cognition",
              "Social", "Social")
  )
  boot <- bootstrap_mcml(wide, clusters = cl_df,
                          iter = 20, seed = 42)
  expect_s3_class(boot, "mcml_bootstrap")
  expect_equal(boot$between$original$n_nodes, 3)
})


# ===========================================================================
# Section 10: Configuration fields
# ===========================================================================
test_that("bootstrap_mcml stores correct config", {
  wide <- .make_boot_wide()
  boot <- bootstrap_mcml(wide, clusters = .boot_clusters,
                          iter = 50, ci_level = 0.10,
                          inference = "stability",
                          consistency_range = c(0.5, 1.5),
                          seed = 42)
  expect_equal(boot$iter, 50L)
  expect_equal(boot$ci_level, 0.10)
  expect_equal(boot$inference, "stability")
  expect_equal(boot$consistency_range, c(0.5, 1.5))
  expect_equal(boot$method, "relative")
})
