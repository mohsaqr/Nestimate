# Skipped on CRAN to keep the check under its time limit; runs locally and in CI.
testthat::skip_on_cran()

# ---- compare_networks(actor = ): nesting-aware permutation backend ----

.cn_actor_nets <- function() {
  build_network(group_regulation_long, method = "relative", actor = "Actor",
                action = "Action", time = "Time", group = "Achiever")
}

test_that("actor is passed to permutation() and reported as Nesting rows", {
  nets <- .cn_actor_nets()
  cmp <- compare_networks(nets, test = "permutation", iter = 99L, seed = 1,
                          actor = "Group")
  ref <- permutation(nets[[1]], nets[[2]], iter = 99L, seed = 2,
                     actor = "Group")
  g <- global_differences(cmp, digits = 10)
  nesting <- g[g$category == "Nesting", ]
  expect_identical(nesting$key, c("icc", "deff_edges", "deff_global"))
  expect_equal(nesting$value,
               c(ref$clustering$icc, ref$clustering$deff_edges,
                 ref$clustering$deff_global), tolerance = 1e-8)
  m_row <- g[g$key == "perm_M", ]
  expect_equal(m_row$perm_p, ref$global$p_value[1])
  expect_true(any(grepl("actor = Group", capture.output(print(cmp)))))
})

test_that("without actor there are no Nesting rows", {
  cmp <- compare_networks(.cn_actor_nets(), test = "permutation", iter = 19L,
                          seed = 1)
  expect_false("Nesting" %in% global_differences(cmp)$category)
  expect_null(cmp$actor)
})

test_that("actor needs the permutation backend", {
  expect_error(compare_networks(.cn_actor_nets(), actor = "Group"),
               class = "nestimate_compare_actor_needs_permutation")
  expect_error(compare_networks(.cn_actor_nets(), test = "permutation",
                                actor = c("a", "b")),
               "single column name")
})

test_that("the global view and summary accept the Nesting rows", {
  cmp <- compare_networks(.cn_actor_nets(), test = "permutation", iter = 19L,
                          seed = 1, actor = "Group")
  expect_s3_class(summary(cmp), "data.frame")
  p <- plot(cmp, type = "global")
  expect_silent(ggplot2::ggplot_build(p))
})
