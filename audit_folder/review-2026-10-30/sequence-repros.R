pkgload::load_all('.', quiet = TRUE)
# S1: all-NA logical numeric position breaks role inference
numeric_gap <- data.frame(id = c('s1', 's2'), T1 = c(1, 2), T2 = NA,
                          T3 = c(2, 1))
print(extract_sequences(numeric_gap))
print(build_network(numeric_gap, method = 'relative')$nodes)
stopifnot(identical(names(extract_sequences(numeric_gap)), c('id', 'T2')))
# Explicit selectors are a valid control/workaround
stopifnot(identical(names(extract_sequences(numeric_gap,
                           cols = c('T1', 'T2', 'T3'))), c('T1','T2','T3')))
# S2: outcome vector overwrites existing adjustment named .outcome
adjustment <- expand.grid(hit = 0:1, score = 0:1, y = 0:1)
adjustment <- adjustment[rep(seq_len(8), c(9,4,3,2,2,3,7,9)), ]
adjustment$T1 <- ifelse(adjustment$hit == 1, 'A', 'B')
adjustment$T2 <- 'B'
fit_pattern <- function(d, a) pattern_outcome(d, outcome = d$y,
  cols = c('T1','T2'), len = 2, n = 1, min_freq = 1,
  start = 'A', adjust = a)$estimate
original <- fit_pattern(adjustment, 'score')
names(adjustment)[names(adjustment) == 'score'] <- '.outcome'
renamed <- fit_pattern(adjustment, '.outcome')
print(c(original = original, renamed = renamed))
stopifnot(abs(original - renamed) > 0.9)
# S3: stars-only literal states are silently interpreted as wildcards
star_states <- data.frame(T1 = c('**', 'A'), T2 = c('B','B'))
print(try(sequence_patterns(star_states, len = 1, min_freq = 1,
                            min_support = 0)))
star_bigrams <- sequence_patterns(star_states, len = 2, min_freq = 1,
                                  min_support = 0)
print(star_bigrams)
stopifnot(star_bigrams$lift[star_bigrams$pattern == '**->B'] == 0.5)
