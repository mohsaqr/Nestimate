pkgload::load_all(quiet = TRUE)

# F-INFER-1: symmetry requires 1/sqrt(n), but norm squaring underflows.
n <- 150L
hg <- bipartite_groups(data.frame(player = seq_len(n), group = 'G'),
                       player = 'player', group = 'group')
warnings <- character()
z <- withCallingHandlers(hypergraph_centrality(hg, type = 'Z')$Z,
                         warning = function(w) warnings <<- c(warnings, conditionMessage(w)))
cat('Z range:', range(z), '\nExpected:', 1 / sqrt(n),
    '\nL2 norm:', sqrt(sum(z^2)), '\nWarnings:', length(warnings), '\n')
stopifnot(max(abs(z - 1/sqrt(n))) > 0.08, length(warnings) == 0L)

# F-INFER-2: path DOFs deliberately use doubles, summary truncates to integer.
set.seed(4)
seqs <- replicate(100, sample(LETTERS[1:10], 14, TRUE), simplify = FALSE)
m <- build_mogen(seqs, max_order = 10L)
sm <- NULL
invisible(capture.output(sm <- summary(m)))
cat('Stored layer DOF:', tail(m$layer_dof, 2), '\nSummary layer DOF:', tail(sm$layer_dof, 2), '\n')
stopifnot(all(is.finite(tail(m$layer_dof, 2))), all(is.na(tail(sm$layer_dof, 2))))

# F-INFER-3: original-zero exclusion ignores new association edges.
set.seed(19)
d <- as.data.frame(matrix(rnorm(80 * 5), 80, 5))
d$V2 <- d$V1 + rnorm(80, sd = .7)
d$V3 <- .4 * d$V1 + rnorm(80)
n <- build_network(d, method = 'cor', threshold = .2)
set.seed(2)
keep <- sample(seq_len(nrow(d)), 40)
sub <- build_network(d[keep, ], method = 'cor', threshold = .2)
es <- casedrop_reliability(n, iter = 2, drop_prop = .5, seed = 2)
cat('New nonzero directed cells:', sum(n$weights == 0 & sub$weights != 0),
    '\nReported max deviation:', es$metrics$max_abs_dev[1, 1],
    '\nWhole-network max deviation:', max(abs(n$weights - sub$weights)), '\n')
stopifnot(sum(n$weights == 0 & sub$weights != 0) == 4L,
          es$metrics$max_abs_dev[1, 1] < max(abs(n$weights - sub$weights)))
