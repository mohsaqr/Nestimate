testthat::skip_on_cran()

test_that("pyramid uses the house colours, leaves thin bars unlabelled and does not clip", {
  # mostly A->A in both groups, with a rare B->A in each: the B->A bars are
  # thin, the A->A bars wide
  rows <- c(rep("A A A A", 19), "A B A A", rep("A A A A", 19), "B B A A")
  seqs <- as.data.frame(do.call(rbind, strsplit(rows, " ")),
                        stringsAsFactors = FALSE)
  names(seqs) <- paste0("T", 1:4)
  seqs$grp <- rep(c("x", "y"), each = 20)
  res <- sequence_compare(seqs, group = "grp", sub = 2, min_freq = 1L,
                          test = "none")
  p <- plot(res, style = "pyramid", show_residuals = TRUE)
  built <- ggplot2::ggplot_build(p)
  fill <- p$scales$get_scales("fill")
  expect_identical(fill$palette(0), "#D33F6A")   # low end of the gradient
  expect_identical(fill$palette(1), "#4A6FE3")   # high end
  expect_identical(p$coordinates$clip, "off")
  # the residual layer: a label exactly where the bar (width = the pattern's
  # proportion) is at least 20% of the widest bar
  rect <- built$data[[1]]
  resid_layer <- built$data[[3]]
  width <- rect$xmax - rect$xmin
  wide <- width >= max(width) * 0.2
  expect_true(any(!wide))               # the data do contain a thin bar
  expect_identical(nzchar(resid_layer$label), wide)
})

test_that("grouped sequence plots name numbered clusters and keep named groups", {
  expect_identical(.group_panel_label(c("1", "12")), c("Cluster 1", "Cluster 12"))
  expect_identical(.group_panel_label(c("early", "late")), c("early", "late"))
  expect_identical(.group_panel_label("Cluster 2"), "Cluster 2")
})

test_that("clustering plots use Okabe-Ito colours and distinguish clusters by shape", {
  seqs <- data.frame(T1 = rep(c("A", "B", "C"), 10), T2 = rep(c("B", "C", "A"), 10),
                     T3 = rep(c("C", "A", "B"), 10), T4 = rep(c("A", "A", "B"), 10),
                     stringsAsFactors = FALSE)
  cl <- build_clusters(seqs, k = 3)
  sil <- plot(cl, type = "silhouette")
  fills <- unique(ggplot2::ggplot_build(sil)$data[[1]]$fill)
  expect_true(all(fills %in% .okabe_ito))
  mds <- plot(cl, type = "mds")
  built <- ggplot2::ggplot_build(mds)$data[[1]]
  expect_true(all(unique(built$colour) %in% .okabe_ito))
  expect_equal(length(unique(built$shape)), length(unique(built$colour)))
})

test_that("the compare_mmm plot uses Okabe-Ito colours and one shape per criterion", {
  seqs <- data.frame(T1 = rep(c("A", "B", "C"), 10), T2 = rep(c("B", "C", "A"), 10),
                     T3 = rep(c("C", "A", "B"), 10), T4 = rep(c("A", "A", "B"), 10),
                     stringsAsFactors = FALSE)
  cmp <- compare_mmm(seqs, k = 2:3, seed = 1, n_starts = 2)
  pdf(NULL); p <- plot(cmp); grDevices::dev.off()
  built <- ggplot2::ggplot_build(p)$data[[2]]
  expect_true(all(unique(built$colour) %in% .okabe_ito))
  expect_equal(length(unique(built$shape)), 3L)
})
