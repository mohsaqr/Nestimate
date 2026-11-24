testthat::skip_on_cran()

.sp_net <- function() {
  w <- matrix(c(0, .5, .3, .2,
                .4, 0, .4, .2,
                .3, .3, 0, .4,
                .2, .5, .3, 0), 4, byrow = TRUE,
              dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  build_network(w, method = "relative")
}

test_that("the persistence plot uses Okabe-Ito with a linetype and shape per dimension", {
  skip_if_not_installed("gridExtra")
  ph <- persistent_homology(.sp_net(), n_steps = 8)
  panels <- plot(ph, combined = FALSE)
  curve <- ggplot2::ggplot_build(panels[[1]])$data[[1]]
  expect_true(all(unique(curve$colour) %in% .okabe_ito))
  expect_equal(length(unique(curve$linetype)), length(unique(curve$colour)))
  diagram <- ggplot2::ggplot_build(panels[[2]])$data[[2]]
  expect_true(all(unique(diagram$colour) %in% .okabe_ito))
  expect_equal(length(unique(diagram$shape)), length(unique(diagram$colour)))
})

test_that("the simplicial complex and Q-analysis plots use Okabe-Ito only", {
  skip_if_not_installed("gridExtra")
  sc <- build_simplicial(.sp_net(), threshold = 0.25)
  panels <- plot(sc, combined = FALSE)
  bar_fills <- unlist(lapply(panels[c("f_vector", "betti", "degree")], \(p)
    ggplot2::ggplot_build(p)$data[[1]]$fill))
  expect_true(all(unique(bar_fills) %in% .okabe_ito))
  # the degree heatmap is a ramp from white to Okabe-Ito orange
  built <- ggplot2::ggplot_build(panels$degree_heatmap)$data
  top <- built[[1]]$fill[which.max(as.numeric(built[[2]]$label))]
  expect_equal(toupper(top), "#E69F00")
  q <- plot(q_analysis(sc), combined = FALSE)
  qcols <- c(ggplot2::ggplot_build(q[[1]])$data[[1]]$colour,
             ggplot2::ggplot_build(q[[1]])$data[[2]]$colour,
             ggplot2::ggplot_build(q[[2]])$data[[1]]$fill)
  expect_true(all(unique(qcols) %in% .okabe_ito))
})
