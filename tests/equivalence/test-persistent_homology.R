# persistent_homology(type = "vr"), persistence_landscape(), bottleneck_distance(): Vietoris-Rips persistence of a
# point cloud's distance matrix and two summaries of a persistence diagram. References: TDAstats (ripser) for the
# diagram, TDA for the landscape (Bubenik 2015) and the bottleneck distance. Re-expressed from
# local_testing_and_equivalence/test-equiv-persistent_homology.R (30 clouds against TDAstats, by maxima); the
# landscape and bottleneck comparisons are new. TDAstats reports no essential class (death beyond the threshold);
# Nestimate reports them with death = Inf, so the finite pairs are compared. The diagram is read from the result's
# field `persistence` (API gap: no accessor).

clouds <- lapply(c(741L, 742L, 743L), \(seed) {
  set.seed(seed)
  matrix(stats::rnorm(24), 12L, 2L)
})

# A diagram as rows numbered in (dimension, birth, death) order, for a comparison by position.
hg_diagram <- function(dimension, birth, death) {
  ordered <- order(dimension, birth, death)
  data.frame(pair = seq_along(ordered), dimension = dimension[ordered], birth = birth[ordered], death = death[ordered])
}

test_that("finite Vietoris-Rips persistence pairs equal TDAstats::calculate_homology()", {
  skip_without_reference("TDAstats")
  invisible(lapply(seq_along(clouds), \(i) {
    distances <- as.matrix(stats::dist(clouds[[i]]))
    scale <- max(distances) * 1.05
    found <- persistent_homology(distances, n_steps = 5L, max_dim = 2L, type = "vr", max_scale = scale)
    finite <- subset(found$persistence, is.finite(death) & dimension <= 1L)
    reference <- as.data.frame(TDAstats::calculate_homology(clouds[[i]], dim = 1L, threshold = scale))
    expect_reference(hg_diagram(finite$dimension, finite$birth, finite$death),
                     hg_diagram(reference$dimension, reference$birth, reference$death), reference = "TDAstats",
                     covers = "persistent_homology()", quantity = sprintf("persistence pairs, cloud %d", i),
                     tolerance = 1e-10, key = "pair")
  }))
})

test_that("persistence landscapes and bottleneck distances equal TDA", {
  skip_without_reference("TDA")
  diagrams <- lapply(clouds, \(cloud) {
    distances <- as.matrix(stats::dist(cloud))
    persistent_homology(distances, n_steps = 5L, max_dim = 2L, type = "vr", max_scale = max(distances) * 1.05)
  })
  as_tda <- \(ph) {
    finite <- subset(ph$persistence, is.finite(death))
    cbind(dimension = finite$dimension, Birth = finite$birth, Death = finite$death)
  }
  grid <- seq(0, 2, length.out = 41L)
  invisible(lapply(seq_along(diagrams), \(i) {
    landscape <- persistence_landscape(diagrams[[i]], k_max = 2L, dimension = 0L, t_grid = grid)
    expected <- data.frame(k = rep(1:2, each = length(grid)), t = rep(grid, 2L),
                           value = as.vector(vapply(1:2, \(k) as.vector(TDA::landscape(as_tda(diagrams[[i]]), dimension = 0L,
                                                                                         KK = k, tseq = grid)),
                                                    numeric(length(grid)))))
    expect_reference(landscape$landscape, expected, reference = "TDA", covers = "persistence_landscape()",
                     quantity = sprintf("landscape, dimension 0, cloud %d", i), tolerance = 1e-12, key = c("k", "t"))
  }))
  # Clouds 1 and 2 have no 1-dimensional class; cloud 3 has three. A point on the diagonal is at bottleneck
  # distance 0 from every diagram, so adding one to both sides leaves the distance unchanged; TDA needs it when one
  # side has no point of the dimension (below).
  on_diagonal <- \(diagram) rbind(diagram, c(1, 0.5, 0.5))
  bottleneck <- \(a, b) bottleneck_distance(diagrams[[a]], diagrams[[b]], dimension = 1L)
  expect_reference(unname(bottleneck(1L, 2L)), TDA::bottleneck(as_tda(diagrams[[1L]]), as_tda(diagrams[[2L]]), dimension = 1L),
                   reference = "TDA", covers = "bottleneck_distance()",
                   quantity = "bottleneck distance, dimension 1, clouds 1 and 2 (both empty)", tolerance = 1e-10)
  invisible(lapply(1:2, \(i) {
    expect_reference(unname(bottleneck(i, 3L)),
                     TDA::bottleneck(on_diagonal(as_tda(diagrams[[i]])), on_diagonal(as_tda(diagrams[[3L]])), dimension = 1L),
                     reference = "TDA", covers = "bottleneck_distance()",
                     quantity = sprintf("bottleneck distance, dimension 1, clouds %d and 3", i), tolerance = 1e-10,
                     note = "a point on the diagonal added to both diagrams")
    expect_known_difference(unname(bottleneck(i, 3L)),
                            TDA::bottleneck(as_tda(diagrams[[i]]), as_tda(diagrams[[3L]]), dimension = 1L),
                            reference = "TDA", covers = "bottleneck_distance()",
                            quantity = sprintf("bottleneck distance, dimension 1, clouds %d and 3, one side empty", i),
                            tolerance = 1e-10,
                            reason = paste("TDA::bottleneck() returns 0.5 when one diagram has no point of the dimension;",
                                           "the distance to an empty diagram is the largest half-persistence, (d - b) / 2,",
                                           "which TDA returns once a point on the diagonal is added to both sides."),
                            agrees_with = "TDA with a diagonal point added to both diagrams")
    finite <- subset(diagrams[[3L]]$persistence, is.finite(death) & dimension == 1L)
    expect_reference(unname(bottleneck(i, 3L)), max((finite$death - finite$birth) / 2),
                     reference = "independent: largest half-persistence (distance to the empty diagram)",
                     covers = "bottleneck_distance()",
                     quantity = sprintf("bottleneck distance to an empty diagram, clouds %d and 3", i), tolerance = 1e-12)
  }))
})
