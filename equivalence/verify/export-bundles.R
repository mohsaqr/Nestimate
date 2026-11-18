# Exports Nestimate comparison cases as bundles for Prova Verify (format:
# plan/24-verify.md in the prova repository), with Prova's R writer.
#
# Every side of a case is an R script, stored in the bundle and run there in its
# own R process: its console is the side's raw output and the table it leaves in
# `result` holds the side's values. Quantities are matched across sides by name.
#
# Run from the package root:
#   PROVA_VERIFY_R=<prova>/reporters/verify.R \
#     Rscript equivalence/verify/export-bundles.R [case ...]
source(Sys.getenv("PROVA_VERIFY_R", file.path("..", "prova", "reporters", "verify.R")))

.verify_base <- file.path("equivalence", "verify")
.verify_generator <- "Nestimate equivalence/verify/export-bundles.R"

# Tolerances, fixed before any case was run.
#
# 1e-10 is the tolerance of the package's own equivalence tests against tna,
# stats and corpcor (`TOL` in local_testing_and_equivalence/test-equiv-estimators.R
# and test-equiv-bootstrap.R): the quantities are counts, ratios of counts, a
# matrix inverse of a well-conditioned 9 x 9 correlation matrix, shortest paths,
# or resampled versions of these under one shared random stream.
.tol_exact <- 1e-10
# 1e-4 is `TOL_GLASSO` of the same tests. qgraph::EBICglasso() runs the glasso
# solver at its default convergence threshold of 1e-4, so two correct solutions
# of the same penalised problem may differ by about that much in a weight.
.tol_glasso <- 1e-4

#' The data files the cases read, written to a temporary folder
#'
#' `sequences.csv` and `sequences_b.csv` are sequences 1 to 300 and 301 to 600
#' of `tna::group_regulation` (nine states, up to 26 steps, missing values where
#' a sequence is shorter), the data most of the package's equivalence tests use.
#' `srl_strategies.csv` is Nestimate's bundled `srl_strategies` (250 students,
#' nine continuous scores).
#' @return A named character vector of file paths.
.verify_data <- function() {
  folder <- file.path(tempdir(), "nestimate-verify-data")
  dir.create(folder, showWarnings = FALSE)
  regulation <- as.data.frame(tna::group_regulation)
  strategies <- new.env()
  load(file.path("data", "srl_strategies.rda"), envir = strategies)
  files <- c(sequences = file.path(folder, "sequences.csv"),
             sequences_b = file.path(folder, "sequences_b.csv"),
             strategies = file.path(folder, "srl_strategies.csv"))
  utils::write.csv(utils::head(regulation, 300L), files[["sequences"]], row.names = FALSE)
  utils::write.csv(utils::head(utils::tail(regulation, -300L), 300L), files[["sequences_b"]], row.names = FALSE)
  utils::write.csv(strategies$srl_strategies, files[["strategies"]], row.names = FALSE)
  files
}

.read_sequences <- "sequences <- read.csv(\"../data/sequences.csv\")"
.read_groups <- c("group_a <- read.csv(\"../data/sequences.csv\")",
                  "group_b <- read.csv(\"../data/sequences_b.csv\")")
.read_strategies <- "strategies <- read.csv(\"../data/srl_strategies.csv\")"

.nestimate <- function(script, key = NULL) {
  verify_side("nestimate", package = "Nestimate", primary = TRUE, load = ".", key = key,
              script = c("library(Nestimate)", script))
}
.reference <- function(package, script, key = NULL, name = package) {
  verify_side(tolower(package), name = name, package = package, authority = "open-source", key = key,
              script = c(sprintf("library(%s)", package), script))
}

#' Every case: its words, its data, its tolerance and the script of each side
#'
#' A script's paths are relative to its own folder inside the bundle. Every
#' script leaves a data frame `result`: long (`quantity`, `value`, `kind`), or
#' wide with one key column and one column per measure.
#' @return A named list of cases.
.verify_cases <- function() {
  list(
    transition_networks = list(
      title = "Transition networks from sequences",
      description = paste(
        "Relative (transition probability), frequency (transition count) and attention (decay-weighted)",
        "networks of 300 sequences over nine states, against tna::tna(), tna::ftna() and tna::atna().",
        "Every non-zero edge is compared, self-transitions included; an edge that only one side has would",
        "show as missing."),
      covers = c("build_network()", "extract_edges()"),
      needs = "tna", data = "sequences", tolerance = .tol_exact,
      sides = list(
        .nestimate(c(
          .read_sequences,
          "relative <- extract_edges(build_network(sequences, method = \"relative\"), include_self = TRUE)",
          "frequency <- extract_edges(build_network(sequences, method = \"frequency\"), include_self = TRUE)",
          "attention <- extract_edges(build_network(sequences, method = \"attention\"), include_self = TRUE)",
          "result <- rbind(transform(relative, kind = \"relative\"),",
          "                transform(frequency, kind = \"frequency\"),",
          "                transform(attention, kind = \"attention\"))",
          "result <- transform(result, quantity = paste0(kind, \": \", from, \" -> \", to), value = weight)",
          "result")),
        .reference("tna", c(
          .read_sequences,
          "edges_of <- function(model, kind) {",
          "  edges <- as.data.frame(as.table(model$weights), stringsAsFactors = FALSE)",
          "  names(edges) <- c(\"from\", \"to\", \"value\")",
          "  edges <- edges[edges$value != 0, ]",
          "  data.frame(quantity = paste0(kind, \": \", edges$from, \" -> \", edges$to), kind = kind,",
          "             value = edges$value)",
          "}",
          "result <- rbind(edges_of(tna(sequences), \"relative\"),",
          "                edges_of(ftna(sequences), \"frequency\"),",
          "                edges_of(atna(sequences), \"attention\"))",
          "result")))),

    correlation_networks = list(
      title = "Correlation and partial correlation networks",
      description = paste(
        "Pearson correlation and partial correlation networks of nine continuous scores (250 students).",
        "Base R computes the partial correlations from the inverse of the correlation matrix; corpcor",
        "computes both matrices with its own routines at zero shrinkage."),
      covers = c("build_network()", "build_cor()", "build_pcor()"),
      needs = "corpcor", data = "strategies", tolerance = .tol_exact,
      sides = list(
        .nestimate(c(
          .read_strategies,
          "correlation <- build_network(strategies, method = \"cor\")",
          "partial <- build_network(strategies, method = \"pcor\")",
          "result <- transform(rbind(as.data.frame(correlation), as.data.frame(partial)),",
          "                    quantity = paste0(network, \": \", from, \" -- \", to), kind = network,",
          "                    value = weight)",
          "result")),
        .reference("stats", name = "stats (base R)", c(
          .read_strategies,
          "correlation <- cor(strategies)",
          "partial <- -cov2cor(solve(correlation))",
          "pairs <- which(upper.tri(correlation), arr.ind = TRUE)",
          "labels <- paste(colnames(correlation)[pairs[, \"row\"]], \"--\", colnames(correlation)[pairs[, \"col\"]])",
          "result <- data.frame(quantity = c(paste0(\"cor: \", labels), paste0(\"pcor: \", labels)),",
          "                     kind = rep(c(\"cor\", \"pcor\"), each = nrow(pairs)),",
          "                     value = c(correlation[pairs], partial[pairs]))",
          "result")),
        .reference("corpcor", c(
          .read_strategies,
          "correlation <- cor.shrink(strategies, lambda = 0, verbose = FALSE)",
          "partial <- pcor.shrink(strategies, lambda = 0, verbose = FALSE)",
          "pairs <- which(upper.tri(correlation), arr.ind = TRUE)",
          "labels <- paste(colnames(correlation)[pairs[, \"row\"]], \"--\", colnames(correlation)[pairs[, \"col\"]])",
          "result <- data.frame(quantity = c(paste0(\"cor: \", labels), paste0(\"pcor: \", labels)),",
          "                     kind = rep(c(\"cor\", \"pcor\"), each = nrow(pairs)),",
          "                     value = c(correlation[pairs], partial[pairs]))",
          "result")))),

    glasso_ebic = list(
      title = "Regularised partial correlation network (EBICglasso)",
      description = paste(
        "Graphical lasso network of nine continuous scores with the penalty chosen by EBIC (gamma 0.5,",
        "100 penalties). Nestimate uses its own pure R solver. qgraph::EBICglasso() runs the glasso",
        "package at its default convergence threshold of 1e-4, which is why the tolerance is 1e-4. The",
        "third side repeats the EBIC selection by hand on glasso::glassopath() solved to 1e-10. Only",
        "non-zero edges are listed, so a difference in which edges are kept would show as missing."),
      covers = c("build_network()", "build_glasso()"),
      needs = c("qgraph", "glasso"), data = "strategies", tolerance = .tol_glasso,
      sides = list(
        .nestimate(c(
          .read_strategies,
          "network <- build_network(strategies, method = \"glasso\")",
          "network",
          "result <- transform(as.data.frame(network), quantity = paste(from, \"--\", to),",
          "                    kind = \"partial correlation\", value = weight)",
          "result")),
        .reference("qgraph", c(
          .read_strategies,
          "fit <- EBICglasso(cor(strategies), n = nrow(strategies), gamma = 0.5, returnAllResults = TRUE)",
          "fit$lambda[which.min(fit$ebic)]",
          "network <- fit$optnet",
          "pairs <- which(upper.tri(network) & network != 0, arr.ind = TRUE)",
          "result <- data.frame(quantity = paste(colnames(network)[pairs[, \"row\"]], \"--\",",
          "                                      colnames(network)[pairs[, \"col\"]]),",
          "                     kind = \"partial correlation\", value = network[pairs])",
          "result")),
        .reference("glasso", c(
          .read_strategies,
          "S <- cor(strategies)",
          "n <- nrow(strategies)",
          "p <- ncol(S)",
          "largest <- max(abs(S[upper.tri(S)]))",
          "penalties <- exp(seq(log(0.01 * largest), log(largest), length.out = 100))",
          "path <- glassopath(S, rholist = penalties, thr = 1e-10, maxit = 1e5, penalize.diagonal = FALSE,",
          "                   trace = 0)",
          "ebic <- apply(path$wi, 3, function(K) {",
          "  edges <- sum(K[lower.tri(K)] != 0)",
          "  loglik <- n / 2 * (log(det(K)) - sum(diag(K %*% S)))",
          "  -2 * loglik + edges * log(n) + 4 * edges * 0.5 * log(p)",
          "})",
          "penalties[which.min(ebic)]",
          "network <- -cov2cor(path$wi[, , which.min(ebic)])",
          "dimnames(network) <- dimnames(S)",
          "pairs <- which(upper.tri(network) & network != 0, arr.ind = TRUE)",
          "result <- data.frame(quantity = paste(colnames(network)[pairs[, \"row\"]], \"--\",",
          "                                      colnames(network)[pairs[, \"col\"]]),",
          "                     kind = \"partial correlation\", value = network[pairs])",
          "result")))),

    centralities = list(
      title = "Centralities of a transition network",
      description = paste(
        "Strength, closeness and betweenness of the relative transition network of 300 sequences,",
        "self-transitions excluded and a path's length the sum of inverse weights. tna estimates the",
        "network itself; the igraph side counts the transitions in base R and uses igraph's measures."),
      covers = c("net_centrality()", "build_network()"),
      needs = c("tna", "igraph"), data = "sequences", tolerance = .tol_exact, key = "state",
      sides = list(
        .nestimate(key = "state", c(
          .read_sequences,
          "network <- build_network(sequences, method = \"relative\")",
          "result <- net_centrality(network, measures = c(\"OutStrength\", \"InStrength\", \"ClosenessIn\",",
          "                                               \"ClosenessOut\", \"Closeness\", \"Betweenness\"))",
          "result")),
        .reference("tna", key = "state", c(
          .read_sequences,
          "model <- tna(sequences)",
          "result <- as.data.frame(centralities(model, measures = c(\"OutStrength\", \"InStrength\", \"ClosenessIn\",",
          "                                                         \"ClosenessOut\", \"Closeness\", \"Betweenness\")))",
          "result")),
        .reference("igraph", key = "state", c(
          .read_sequences,
          "states <- sort(unique(na.omit(unlist(sequences))))",
          "from <- factor(unlist(sequences[-ncol(sequences)]), levels = states)",
          "to <- factor(unlist(sequences[-1]), levels = states)",
          "weights <- unclass(prop.table(table(from, to), margin = 1))",
          "diag(weights) <- 0",
          "graph <- graph_from_adjacency_matrix(weights, mode = \"directed\", weighted = TRUE)",
          "lengths <- 1 / E(graph)$weight",
          "result <- data.frame(state = V(graph)$name,",
          "                     OutStrength = strength(graph, mode = \"out\"),",
          "                     InStrength = strength(graph, mode = \"in\"),",
          "                     ClosenessIn = closeness(graph, mode = \"in\", weights = lengths),",
          "                     ClosenessOut = closeness(graph, mode = \"out\", weights = lengths),",
          "                     Closeness = closeness(graph, mode = \"all\", weights = lengths),",
          "                     Betweenness = betweenness(graph, weights = lengths))",
          "result")))),

    centralities_tna = list(
      title = "Centralities only tna also computes",
      description = paste(
        "Randomised shortest path betweenness, diffusion and weighted clustering of the same network.",
        "The installed tna returns diffusion on its raw scale, so Nestimate is asked for the raw scale",
        "too (normalize_diffusion = FALSE; its default is the 0 to 1 scale)."),
      covers = "net_centrality()",
      variant = list(of = "centralities", option = "measures", value = "BetweennessRSP, Diffusion, Clustering"),
      needs = "tna", data = "sequences", tolerance = .tol_exact,
      sides = list(
        .nestimate(key = "state", c(
          .read_sequences,
          "network <- build_network(sequences, method = \"relative\")",
          "result <- net_centrality(network, measures = c(\"BetweennessRSP\", \"Diffusion\", \"Clustering\"),",
          "                         normalize_diffusion = FALSE)",
          "result")),
        .reference("tna", key = "state", c(
          .read_sequences,
          "model <- tna(sequences)",
          "result <- as.data.frame(centralities(model, measures = c(\"BetweennessRSP\", \"Diffusion\", \"Clustering\")))",
          "result")))),

    bootstrap_edges = list(
      title = "Bootstrap of the edges of a transition network",
      description = paste(
        "200 bootstrap resamples of the sequences, seed 1 on both sides. Both programs draw the same",
        "resamples from the same random stream, so the bootstrap mean, standard deviation, stability",
        "p-value, significance, percentile interval and consistency range of every edge are compared",
        "directly, not as Monte Carlo estimates."),
      covers = c("bootstrap_network()", "build_network()"),
      needs = "tna", data = "sequences", tolerance = .tol_exact,
      sides = list(
        .nestimate(key = "edge", c(
          .read_sequences,
          "network <- build_network(sequences, method = \"relative\")",
          "boot <- bootstrap_network(network, iter = 200, seed = 1)",
          "boot",
          "result <- subset(transform(summary(boot), edge = paste(from, \"->\", to)),",
          "                 select = c(edge, weight, mean, sd, p_value, sig, ci_lower, ci_upper, cr_lower, cr_upper))",
          "result")),
        .reference("tna", key = "edge", c(
          .read_sequences,
          "model <- tna(sequences)",
          "set.seed(1)",
          "boot <- bootstrap(model, iter = 200)",
          "edges <- boot$summary",
          "at <- cbind(match(edges$from, model$labels), match(edges$to, model$labels))",
          "result <- data.frame(edge = paste(edges$from, \"->\", edges$to), weight = edges$weight,",
          "                     mean = boot$weights_mean[at], sd = boot$weights_sd[at],",
          "                     p_value = edges$p_value, sig = edges$sig,",
          "                     ci_lower = edges$ci_lower, ci_upper = edges$ci_upper,",
          "                     cr_lower = edges$cr_lower, cr_upper = edges$cr_upper)",
          "result")))),

    permutation_differences = list(
      title = "Permutation test: observed edge differences",
      description = paste(
        "Two groups of 300 sequences compared edge by edge. Compared here: the two networks and their",
        "observed differences, which involve no randomness. The p-values and effect sizes are NOT",
        "compared: the two programs draw their permutations differently (tna permutes the sequences,",
        "Nestimate permutes precomputed per-sequence counts), so under one seed they are two Monte",
        "Carlo estimates of the same quantity, not the same numbers. Both are printed in the raw output."),
      covers = c("permutation()", "build_network()"),
      needs = "tna", data = c("sequences", "sequences_b"), tolerance = .tol_exact,
      sides = list(
        .nestimate(key = "edge", c(
          .read_groups,
          "test <- permutation(build_network(group_a, method = \"relative\"),",
          "                    build_network(group_b, method = \"relative\"), iter = 200, seed = 1)",
          "test",
          "summary(test)",
          "result <- subset(transform(summary(test), edge = paste(from, \"->\", to)),",
          "                 select = c(edge, weight_x, weight_y, diff))",
          "result")),
        .reference("tna", key = "edge", c(
          .read_groups,
          "model_a <- tna(group_a)",
          "model_b <- tna(group_b)",
          "set.seed(1)",
          "test <- permutation_test(model_a, model_b, iter = 200)",
          "as.data.frame(test$edges$stats)",
          "tested <- which(model_a$weights != 0 | model_b$weights != 0, arr.ind = TRUE)",
          "result <- data.frame(edge = paste(model_a$labels[tested[, \"row\"]], \"->\", model_a$labels[tested[, \"col\"]]),",
          "                     weight_x = model_a$weights[tested], weight_y = model_b$weights[tested],",
          "                     diff = test$edges$diffs_true[tested])",
          "result")))),

    centrality_stability = list(
      title = "Case-dropping stability of centralities",
      description = paste(
        "Sequences are dropped in proportions 0.1 to 0.9, 200 times each, seed 1 on both sides, and the",
        "centralities of each reduced network are correlated with the original ones. Compared: the mean",
        "and standard deviation of those correlations and the proportion at or above 0.7, per measure",
        "and proportion. The CS coefficient is a function of the last. Agreement at this tolerance",
        "requires both programs to drop the same sequences from the same random stream."),
      covers = "centrality_stability()",
      needs = "tna", data = "sequences", tolerance = .tol_exact,
      sides = list(
        .nestimate(key = "drop", c(
          .read_sequences,
          "network <- build_network(sequences, method = \"relative\")",
          "stability <- centrality_stability(network, measures = c(\"InStrength\", \"OutStrength\", \"Betweenness\"),",
          "                                  iter = 200, seed = 1)",
          "stability",
          "result <- subset(transform(summary(stability), drop = paste0(measure, \", dropping \", drop_prop)),",
          "                 select = c(drop, mean_cor, sd_cor, prop_above))",
          "result")),
        .reference("tna", key = "drop", c(
          .read_sequences,
          "model <- tna(sequences)",
          "set.seed(1)",
          "stability <- estimate_cs(model, measures = c(\"InStrength\", \"OutStrength\", \"Betweenness\"), iter = 200)",
          "stability",
          "result <- do.call(rbind, lapply(names(stability), function(measure) {",
          "  correlations <- stability[[measure]]$correlations",
          "  data.frame(drop = paste0(measure, \", dropping \", attr(stability, \"drop_prop\")),",
          "             mean_cor = colMeans(correlations), sd_cor = apply(correlations, 2, sd),",
          "             prop_above = colMeans(correlations >= 0.7))",
          "}))",
          "result"))))
  )
}

#' Export one case
#'
#' A case whose reference package is not installed is not exported: a classed
#' message (`nestimate_verify_skipped`) says which package is missing.
#' @param name The case's name.
#' @param data The data files, as [.verify_data()] returns them.
#' @return One row: the bundle, its sides and how its compared values came out,
#'   or `NULL` for a skipped case.
export_case <- function(name, data) {
  cases <- .verify_cases()
  stopifnot("`name` must be one exported case" = length(name) == 1L && name %in% names(cases))
  case <- cases[[name]]
  absent <- case$needs[!vapply(case$needs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent) > 0L) {
    message(messageCondition(
      sprintf("Case `%s` skipped: reference package not installed: %s\n", name, toString(absent)),
      class = "nestimate_verify_skipped"))
    return(NULL)
  }
  labels <- c(sequences = "Sequences (group A)", sequences_b = "Sequences (group B)",
              strategies = "Self-regulated learning scores")
  files <- stats::setNames(data[case$data], labels[case$data])
  verify_export(bundle = file.path(.verify_base, "bundles", name), title = case$title,
                description = case$description, sides = case$sides, data = files,
                tolerance = case$tolerance, measure = "absolute", covers = case$covers,
                variant = case$variant, generator = .verify_generator)
}

requested <- commandArgs(trailingOnly = TRUE)
if (length(requested) == 0L) requested <- names(.verify_cases())
stopifnot("tna is needed to write the sequence data files" = requireNamespace("tna", quietly = TRUE))
exported <- lapply(requested, export_case, data = .verify_data())
print(do.call(rbind, exported), row.names = FALSE)
