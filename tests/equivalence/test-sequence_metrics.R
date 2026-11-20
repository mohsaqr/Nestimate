# sequence_metrics(): per-sequence indices (length, spells, entropy, transitions, first, last and dominant states,
# integrative potential, the composite complexity index). Reference: TraMineR 2.2 (seqdef with left = NA, gaps = NA,
# right = "DEL"; seqindic, seqistatd, seqdss, seqdur, seqintegr), and for complexity_index the composite its
# documentation states (0.4 entropy + 0.4 transition rate + 0.2 min(CV of spell durations, 1); not TraMineR's
# seqici) written on TraMineR's spells. Re-expressed from local_testing_and_equivalence/
# test-equiv-sequence-metrics-traminer.R (ten datasets compared by agreement shares); here three datasets.
# Where a sequence has an internal gap TraMineR joins the states on both sides of the gap into one run and
# Nestimate does not: those sequences are compared separately as a known difference.

# TraMineR's indices in sequence_metrics()' columns, one row per sequence (`sequence` = row number).
sq_traminer <- function(wide, favorable, omega) {
  alphabet <- sort(unique(stats::na.omit(unlist(wide, use.names = FALSE))))
  sequences <- suppressMessages(TraMineR::seqdef(wide, left = NA, gaps = NA, right = "DEL", alphabet = alphabet,
                                                 cpal = rep("#999999", length(alphabet))))
  indices <- TraMineR::seqindic(sequences, indic = c("lgth", "nonm", "visited", "transp", "meand", "entr"))
  shares <- TraMineR::seqistatd(sequences, prop = TRUE)
  spells <- as.matrix(TraMineR::seqdss(sequences))
  spells[spells %in% c("%", "*")] <- NA
  durations <- as.matrix(TraMineR::seqdur(sequences))
  rows <- seq_len(nrow(spells))
  dominant <- alphabet[max.col(shares, ties.method = "first")]
  variation <- apply(durations, 1L, \(d) {
    d <- d[!is.na(d)]
    if (length(d) == 1L) NA_real_ else min(stats::sd(d) / mean(d), 1)
  })
  data.frame(
    sequence = rows, valid_n = indices$NonM, valid_proportion = indices$NonM / indices$Lgth,
    unique_states = indices$Visited, mean_spell_duration = indices$MeanD,
    max_spell_duration = apply(durations, 1L, max, na.rm = TRUE), longitudinal_entropy = indices$Entr,
    simpson_diversity = 1 - rowSums(shares^2), self_loop_tendency = 1 - indices$Transp, transition_rate = indices$Transp,
    initial_state_persistence = durations[, 1L] / indices$Lgth,
    initial_state_proportion = shares[cbind(rows, match(spells[, 1L], alphabet))],
    first_state = spells[, 1L], last_state = spells[cbind(rows, rowSums(!is.na(spells)))],
    dominant_state = dominant, dominant_proportion = apply(shares, 1L, max),
    dominant_max_spell = vapply(rows, \(i) max(durations[i, which(spells[i, ] == dominant[i])]), numeric(1)),
    integrative_potential = as.vector(TraMineR::seqintegr(sequences, state = favorable, pow = omega)),
    complexity_index = ifelse(is.na(variation), 0, 0.4 * indices$Entr + 0.4 * indices$Transp + 0.2 * variation),
    stringsAsFactors = FALSE)
}

# Sequences with an internal gap: a missing position between two observed ones.
sq_has_gap <- function(wide) {
  observed <- !is.na(as.matrix(wide))
  rowSums(observed) < max.col(observed, "last") - max.col(observed, "first") + 1L
}

sq_numbered <- function(metrics) cbind(sequence = seq_len(nrow(metrics)), metrics)

env <- new.env()
utils::data(list = c("biofam", "ex1"), package = "TraMineR", envir = env)
data_sets <- list(
  biofam = as.data.frame(lapply(env$biofam[seq_len(300L), 10:25], as.character), stringsAsFactors = FALSE),
  # tna stores the void positions at the end of a sequence as "%"
  engagement = as.data.frame(lapply(tna::engagement, \(column) {
    column <- as.character(column)
    ifelse(column == "%", NA_character_, column)
  }), stringsAsFactors = FALSE)
)

test_that("sequence indices equal TraMineR on complete sequences", {
  skip_without_reference("TraMineR")
  invisible(lapply(names(data_sets), \(name) {
    wide <- data_sets[[name]]
    favorable <- names(which.max(table(unlist(wide))))
    invisible(lapply(c(1, 2), \(omega) {
      expect_reference(sq_numbered(sequence_metrics(wide, favorable = favorable, omega = omega)),
                       sq_traminer(wide, favorable, omega), reference = "TraMineR", covers = "sequence_metrics()",
                       quantity = sprintf("%s, omega %g", name, omega), tolerance = 1e-12, key = "sequence")
    }))
  }))
})

test_that("on sequences with leading missing positions and internal gaps, TraMineR joins the runs across a gap", {
  skip_without_reference("TraMineR")
  wide <- as.data.frame(lapply(env$ex1[seq_len(6L), 1:13], as.character), stringsAsFactors = FALSE)
  gap <- sq_has_gap(wide)
  found <- sq_numbered(sequence_metrics(wide, favorable = "A"))
  expected <- sq_traminer(wide, "A", 1)
  expect_reference(subset(found, !gap), subset(expected, !gap), reference = "TraMineR", covers = "sequence_metrics()",
                   quantity = "ex1, no internal gap", tolerance = 1e-12, key = "sequence")
  sensitive <- c("sequence", "mean_spell_duration", "max_spell_duration", "self_loop_tendency", "transition_rate",
                 "initial_state_persistence", "dominant_max_spell", "complexity_index")
  expect_known_difference(subset(found, gap, select = sensitive), subset(expected, gap, select = sensitive),
                          reference = "TraMineR", covers = "sequence_metrics()", quantity = "ex1, internal gap",
                          tolerance = 1e-12, key = "sequence",
                          reason = paste("TraMineR (gaps = NA) joins the states before and after a missing position into",
                                         "one spell and one transition; sequence_metrics() ends a spell at a gap, so",
                                         "spell lengths and transition shares differ on those sequences."))
})
