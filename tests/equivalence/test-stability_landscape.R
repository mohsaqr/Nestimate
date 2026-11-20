# stability_landscape(): the cup depth of each state, -log(1 - P_ii) (Kramers 1940; Scheffer et al. 2009), its
# expected dwell time 1 / (1 - P_ii), the stay probability P_ii, the number of transitions leaving the state and
# the state's share of the observed time steps, written from the weights of tna::tna() and tna::ftna(). The percentile intervals come from
# resamples drawn inside the function with no argument to supply them, and are not compared.

test_that("depth, dwell time, stay probability and outgoing transitions equal their values from tna's weights", {
  skip_without_reference("tna")
  sequences <- tna::group_regulation
  stay <- diag(tna::tna(sequences)$weights)
  outgoing <- rowSums(tna::ftna(sequences)$weights)
  # share: the state's share of all observed time steps (every non-missing position of the sequences).
  positions <- table(factor(unlist(sequences, use.names = FALSE), levels = names(stay)))
  expected <- data.frame(group = "network", state = names(stay), depth = -log(1 - stay), dwell = 1 / (1 - stay),
                         stay = unname(stay), share = as.vector(positions / sum(positions)), n_from = unname(outgoing))
  found <- stability_landscape(build_network(sequences, method = "relative"), iter = 20L, seed = 1)
  expect_reference(summary(found), expected, reference = "tna",
                   covers = c("stability_landscape()", "summary.net_landscape()"),
                   quantity = "group_regulation", tolerance = 1e-12, key = c("group", "state"),
                   note = "the reference is the documented definition applied to tna::tna() and tna::ftna() weights")
})

# 90 sequences of three states whose stay probability rises with the sequence, so stability varies along them.
sl_wide <- function(seed) {
  set.seed(seed)
  states <- c("A", "B", "C")
  wide <- as.data.frame(do.call(rbind, lapply(seq_len(90L), \(i) {
    stay <- 0.2 + 0.7 * i / 90
    Reduce(\(previous, t) if (stats::runif(1L) < stay) previous else sample(setdiff(states, previous), 1L),
           seq_len(19L), accumulate = TRUE, sample(states, 1L))
  })))
  names(wide) <- paste0("T", seq_len(20L))
  wide
}

# Each slice of the documented surface: tricube weights over the k nearest values of the measure, the half-width the
# distance to the k-th nearest (at least 1/50 of the range) times `widen`, and the depth -log(1 - P_ii) of the
# weighted pooled transitions; effective_n = (sum w)^2 / sum w^2.
sl_surface <- function(wide, measure, k, n_grid, widen) {
  sequences <- tn_sequences(wide)
  states <- sort(unique(unlist(sequences)))
  per_sequence <- lapply(sequences, \(s) tn_pair_counts(list(s), states))
  grid <- seq(stats::quantile(measure, 0.05), stats::quantile(measure, 0.95), length.out = n_grid)
  do.call(rbind, lapply(grid, \(y) {
    distance <- abs(measure - y)
    half_width <- max(sort(distance)[[k]], diff(range(measure)) / 50) * widen
    w <- ifelse(distance < half_width, (1 - (distance / half_width)^3)^3, 0)
    pooled <- Reduce(`+`, Map(`*`, per_sequence, w))
    data.frame(value = y, state = states, depth = -log(1 - diag(pooled) / rowSums(pooled)),
               effective_n = sum(w)^2 / sum(w^2))
  }))
}

test_that("the surface along longitudinal entropy follows the tricube-weighted pooled transitions (TraMineR's entropy)", {
  skip_without_reference("TraMineR")
  wide <- sl_wide(471L)
  sequences <- suppressMessages(TraMineR::seqdef(wide))
  entropy <- TraMineR::seqindic(sequences, indic = "entr")$Entr
  found <- summary(landscape_surface(wide, measure = "longitudinal_entropy", k = 20L, n_grid = 8L), ends = FALSE)
  found$slice <- match(found$value, unique(found$value))
  expected <- function(widen) {
    surface <- sl_surface(wide, entropy, k = 20L, n_grid = 8L, widen = widen)
    surface$slice <- match(surface$value, unique(surface$value))
    surface
  }
  # The package widens the half-width by 1 % so the k-th nearest sequence keeps a positive weight.
  expect_reference(subset(found, select = c(slice, state, value, depth, effective_n)), expected(1.01),
                   reference = "independent: tricube-weighted pooled transitions on TraMineR's entropy, base R", covers = c("landscape_surface()", "summary.net_landscape_surface()"),
                   quantity = "surface slices, half-width 1.01 x the k-th nearest distance", tolerance = 1e-10,
                   key = c("slice", "state"),
                   note = "the axis is TraMineR::seqindic() Entr; the kernel and the depths are written from the definition")
  expect_reference(subset(found, select = c(slice, state, depth, effective_n)),
                   subset(expected(1), select = c(slice, state, depth, effective_n)),
                   reference = "TraMineR", covers = c("landscape_surface()", "summary.net_landscape_surface()"),
                   quantity = "surface slices, half-width as documented", tolerance = 1e-10, key = c("slice", "state"),
                   defect = paste("the kernel half-width is 1.01 times the distance to the k-th nearest sequence, not that",
                                  "distance as ?landscape_surface states; FINDINGS.md"))
})

# Stationary distribution of a transition matrix: the left eigenvector of eigenvalue 1, summing to 1.
sl_stationary <- function(p) {
  leading <- Re(eigen(t(p))$vectors[, 1L])
  leading / sum(leading)
}

test_that("phase windows give the depth and the occupancy landscape of tna::tna() on each window's columns", {
  skip_without_reference("tna")
  wide <- sl_wide(472L)
  window <- 6L
  starts <- seq_len(ncol(wide) - window + 1L)
  models <- lapply(starts, \(t) tna::tna(wide[, t:(t + window - 1L)]))
  counts <- lapply(starts, \(t) tna::ftna(wide[, t:(t + window - 1L)]))
  states <- c("A", "B", "C")
  labels <- sprintf("%d-%d", starts, starts + window - 1L)
  stay <- do.call(rbind, Map(\(model, count, label) {
    p <- unclass(model$weights)
    data.frame(label = label, state = states, depth = -log(1 - diag(p)), n_from = rowSums(unclass(count$weights)))
  }, models, counts, labels))
  found <- summary(landscape_trend(wide, along = "phase", window = window, iter = 0L))
  expect_reference(subset(found, select = c(label, state, depth, n_from)), stay, reference = "tna",
                   covers = c("landscape_trend()", "summary.net_landscape_trend()"), quantity = "depth per window",
                   tolerance = 1e-10, key = c("label", "state"),
                   note = "-log(1 - P_ii) of tna::tna() and the outgoing counts of tna::ftna() on each window's columns")
  # Occupancy: wells -log(pi_i); ridges -log((pi_i P_ij + pi_j P_ji) / 2) between neighbours in state order and between
  # the last state and the first (Wang, Xu & Wang 2008).
  occupancy <- do.call(rbind, Map(\(model, label) {
    p <- unclass(model$weights)
    pi <- sl_stationary(p)
    following <- c(2L, 3L, 1L)
    flux <- (pi * p[cbind(1:3, following)] + pi[following] * p[cbind(following, 1:3)]) / 2
    data.frame(label = label, state = c(states, paste(states, states[following], sep = " | ")),
               potential = c(-log(pi), -log(flux)))
  }, models, labels))
  found <- summary(landscape_trend(wide, along = "phase", depth = "occupancy", window = window, iter = 0L))
  expect_reference(subset(found, select = c(label, state, potential)), occupancy, reference = "tna",
                   covers = c("landscape_trend()", "summary.net_landscape_trend()"),
                   quantity = "occupancy wells and ridges per window", tolerance = 1e-10, key = c("label", "state"),
                   note = "stationary distribution of tna::tna()'s weights from eigen(), written from the definition")
})
