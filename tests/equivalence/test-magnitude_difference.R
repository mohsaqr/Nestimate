# magnitude_difference(): the per-edge disagreement of the frequency (FTNA) and probability (TNA) views of one set
# of sequences, written from its definition on the weights of tna::ftna() and tna::tna(): the counts rescaled
# linearly into the range of the probabilities ("tna_range"), or both scaled to [0, 1] ("minmax"), then the
# absolute difference, or the absolute difference of the two vectors scaled to unit length ("chord_dist").
# The edge table is a field of the result (API gap).

rs_magnitude <- function(counts, probabilities, scale, metric) {
  rescale <- \(v, low, high) low + (v - min(v)) / (max(v) - min(v)) * (high - low)
  scaled <- switch(scale,
                   tna_range = list(ftna = rescale(counts, min(probabilities), max(probabilities)), tna = probabilities),
                   minmax = list(ftna = rescale(counts, 0, 1), tna = rescale(probabilities, 0, 1)))
  value <- switch(metric,
                  abs_diff = abs(scaled$tna - scaled$ftna),
                  chord_dist = abs(scaled$tna / sqrt(sum(scaled$tna^2)) - scaled$ftna / sqrt(sum(scaled$ftna^2))))
  rs_edge_table(ftna = scaled$ftna, tna = scaled$tna, signed = scaled$tna - scaled$ftna, value = value)
}

test_that("magnitude differences equal their definition on tna's frequency and probability weights", {
  skip_without_reference("tna")
  sequences <- tna::group_regulation
  counts <- tna::ftna(sequences)$weights
  probabilities <- tna::tna(sequences)$weights
  settings <- list(c("tna_range", "abs_diff"), c("minmax", "abs_diff"), c("tna_range", "chord_dist"))
  invisible(lapply(settings, \(setting) {
    found <- magnitude_difference(sequences, scale = setting[1L], metric = setting[2L])
    expect_reference(found$edges, rs_magnitude(counts, probabilities, setting[1L], setting[2L]),
                     reference = "tna", covers = "magnitude_difference()",
                     quantity = sprintf("group_regulation, scale %s, metric %s", setting[1L], setting[2L]),
                     tolerance = 1e-12, key = c("from", "to"),
                     note = "the reference is the documented definition applied to tna::ftna() and tna::tna() weights")
  }))
})
