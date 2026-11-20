# certainty(): the closed-form Dirichlet-Multinomial posterior of each transition probability. With a
# Dirichlet(c_i + prior) posterior for the outgoing transitions of state i, each edge is marginally
# Beta(c_ij + prior, n_i + K * prior - c_ij - prior); its mean, standard deviation, equal-tailed interval and the
# posterior mass outside the band consistency_range * weight are written here from that definition with
# stats::qbeta() and stats::pbeta(). The posterior table is read with summary() (the net_bootstrap method the
# result inherits).

wide <- tn_simulate(n_actors = 40L, n_states = 5L, seq_length = 10L, seed = 401L, ragged = TRUE)

# The posterior of every non-zero edge, one row per edge, from the transition counts.
tn_certainty <- function(counts, prior, ci_level, band) {
  k <- ncol(counts)
  totals <- rowSums(counts)
  weights <- tn_row_normalise(counts)
  cells <- which(weights > 0, arr.ind = TRUE)
  alpha <- counts[cells] + prior
  beta <- totals[cells[, "row"]] + k * prior - alpha
  edge_weight <- weights[cells]
  data.frame(from = rownames(counts)[cells[, "row"]], to = colnames(counts)[cells[, "col"]],
             weight = edge_weight, mean = alpha / (alpha + beta),
             sd = sqrt(alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))),
             ci_lower = stats::qbeta(ci_level / 2, alpha, beta), ci_upper = stats::qbeta(1 - ci_level / 2, alpha, beta),
             cr_lower = min(band) * edge_weight, cr_upper = max(band) * edge_weight,
             p_value = 1 - (stats::pbeta(max(band) * edge_weight, alpha, beta) -
                              stats::pbeta(min(band) * edge_weight, alpha, beta)))
}

test_that("posterior means, standard deviations, intervals and tail masses equal the Beta marginals", {
  counts <- tn_pair_counts(tn_sequences(wide))
  cases <- list(list(prior = 0.5, ci_level = 0.05, band = c(0.75, 1.25)),
                list(prior = 1, ci_level = 0.1, band = c(0.5, 1.5)))
  invisible(lapply(cases, \(case) {
    found <- summary(certainty(build_tna(wide), prior = case$prior, ci_level = case$ci_level,
                               consistency_range = case$band))
    # Closed-form Beta quantities: agreement to rounding, but qbeta() inverts numerically (about 1e-12).
    expect_reference(found, tn_certainty(counts, case$prior, case$ci_level, case$band),
                     reference = "independent: Beta marginals of the Dirichlet posterior, stats::qbeta() and pbeta()",
                     covers = "certainty()", quantity = sprintf("posterior, prior %s, level %s", case$prior, case$ci_level),
                     tolerance = 1e-10, key = c("from", "to"))
  }))
})
