# ---- Bootstrap Network Estimation ----

#' Bootstrap a network estimate
#'
#' @description
#' `bootstrap_network()` resamples the cases a network was estimated from,
#' re-estimates the network on every resample with the same method and
#' settings, and reports for each edge its bootstrap mean, standard
#' deviation, confidence interval and a p-value, together with the network
#' pruned to the edges that pass. It reads any network from
#' [build_network()], built-in method or registered estimator, as well as a
#' `netobject_group` or `mcml` (every network is bootstrapped) and a
#' `wtna_mixed` (both of its components). The cases are the sequences of a
#' transition network (for one-hot data, its actors or sessions) and the
#' rows of the data otherwise; with `actor`, whole actors are resampled (see
#' *Nested data and `actor`*).
#'
#' The p-value comes from one of two tests. With `inference = "stability"`
#' it counts the replicates that fall outside a consistency range around
#' the original weight (by default 75% to 125% of it), so a significant edge
#' is one whose weight the resamples reproduce within that range. With
#' `inference = "threshold"` it counts the replicates whose absolute weight
#' falls below `edge_threshold`, so a significant edge stays above that
#' threshold. In both, the p-value is that count plus one, divided by the
#' number of replicates plus one, so it is never 0. Edges with a p-value
#' below `ci_level` are kept in the pruned network.
#'
#' The p-value counts replicates, so it cannot fall below
#' `1 / (iter + 1)`: with too few iterations no edge can pass a small
#' `ci_level`. The stability test asks whether a weight is reproduced
#' within a range, not whether it differs from zero, and because the range
#' is a proportion of the weight, a weak edge has to be reproduced within a
#' narrow band; read the interval beside the test. The bootstrap describes
#' how the estimate varies over resamples of one sample; to compare groups,
#' use [permutation()].
#'
#' @details
#' The two tests and the inclusive boundary rule follow `tna::bootstrap()`.
#'
#' A latency network (`method = "latency"`) is resampled by its sequences, and
#' every edge's weight is recomputed from the replies of the drawn sequences.
#' A percentile cap the network resolved from its own replies is resolved
#' again in every resample; a cap in seconds, or one pooled over the groups of
#' a grouped build, keeps its value. A resample whose resolved cap is 0
#' seconds (its replies mostly instant) stops the bootstrap with an error of
#' class `nestimate_latency_zero_cap`; give the cap in seconds. Per-unit caps (`latency_scope` other than
#' `"global"`) and `actor` raise an error of class
#' `nestimate_latency_resampling` or `nestimate_actor_unsupported`; to resample
#' within groups, build the network with `group` and bootstrap the group. With
#' the default `aggregate = "mean"`, the stability test is run on the speed
#' score `(weight - floor) / (1 - floor)`, which runs from 0 (every reply at or
#' past the cap) to 1 (every reply instant), because a band around a mean
#' weight between the floor and 1 passes almost every edge. An edge whose
#' replies come from fewer than 5 sequences cannot vary across resamples and
#' gets a p-value of 1. The summary adds each edge's `n_replies`,
#' `n_sequences` and `cv`, the coefficient of variation of its resampled
#' speed scores. These rules were checked by simulation against known reply
#' speeds: with a percentile cap the 95% intervals covered 0.91-0.95 of the
#' true weights from 30 sequences on, and an edge called stable was within
#' 25% of its true score in 97-100% of datasets from 300 sequences on.
#' Bootstrapping latency networks is experimental: every call says so in a
#' message of class `nestimate_experimental` (once for a group of networks).
#'
#' Transition networks (`"relative"`, `"frequency"`, `"co_occurrence"`, and
#' `wtna` networks, which hold the same one-hot counts) take a fast path:
#' the count matrix of every sequence is computed once, and each replicate
#' sums the counts of the sequences drawn, then applies the network's
#' normalisation, scaling and threshold. With `params = list(concat = k)`,
#' the consecutive groups of up to `k` rows merged by the estimator are the
#' sequences resampled; joins inside each merged sequence are retained.
#' With `actor`, a merged sequence must belong entirely to one actor. Every
#' other method (the association methods, `"attention"`, `"ngram"`, `"gap"`,
#' `"reverse"` and registered estimators) is re-estimated in full on
#' resampled rows. A replicate whose estimation fails, or
#' that loses a node, is left out of the p-values, with a warning of class
#' `nestimate_resample_failed` that counts them, and the bootstrap
#' statistics are computed over the remaining replicates; if every
#' replicate fails, every p-value is 1.
#'
#' A network built from an edge list (for example an `mcml` from
#' [build_mcml()] with `from` and `to`) has no sequences, so each
#' transition is resampled as an independent case, with a warning of class
#' `nestimate_edgelist_resampling` that points to the alternatives. An
#' `mcml` built with a `method` other than `"sum"` (such as `"mean"` or
#' `"max"`) combines the transitions of each resample by that method, so
#' every replicate is the `mcml` that [build_mcml()] gives on the resampled
#' rows; `actor` is then not supported and raises an error of class
#' `nestimate_actor_unsupported`.
#'
#' @section Nested data and \code{actor}:
#' The bootstrap resamples sequences as independent units. When sequences
#' are nested in actors (sessions in students, students in teams),
#' \code{actor} names the column identifying the actor, and whole actors
#' are resampled with replacement, keeping all their sequences together:
#' the cluster bootstrap that resamples at the top level only (Davison &
#' Hinkley, 1997, section 3.8; Field & Welsh, 2007). The number of
#' sequences per replicate then varies with the actors drawn.
#'
#' Several columns, outermost first (`actor = c("Group", "Student")`, students
#' within teams), name nested levels. Every column is checked, but only the
#' outermost level is resampled: whole teams, each with all its students
#' and sessions. Resampling inside the drawn clusters as well (students
#' within teams) adds their variation a second time and overstates the
#' variance (Davison & Hinkley, 1997, section 3.8). The ICC and design
#' effect are then those of the outermost level, and `print()` names it as
#' the resampled level.
#'
#' With \code{actor}, the result also reports the nesting effect. The ICC
#' is the proportion of the total variance that lies between actors (Shrout
#' & Fleiss, 1979); an ICC close to 0 indicates little evidence of a nesting
#' effect. It is computed as in \code{\link{permutation}}. The design effect
#' is the ratio of the variance under the nested design to the variance had
#' the sequences been sampled independently (Kish, 1965): here, the variance
#' of the edge weights over actor-level replicates divided by their variance
#' over sequence-level replicates drawn in the same run, reported as the
#' median over edges. \code{actor} is available for transition networks
#' (\code{"relative"}, \code{"frequency"}, \code{"co_occurrence"}).
#'
#' When to set \code{actor}: resampling sequences describes the variability
#' of these actors' sequences; to generalise to other actors (other
#' students, other teams), resample actors. A drawn actor brings all of
#' their sequences, so a replicate can contain an actor twice or not at
#' all. The design effect, not the ICC alone, decides how much this
#' matters: Kish's approximation `1 + (m - 1) * ICC` for clusters of m
#' equally informative units grows with the cluster size, so a small ICC
#' over large clusters still widens the intervals. In one data set of
#' sessions nested in students, the transitions of one student correlated
#' at about 0.01, yet with about 59 transitions per student (size-weighted
#' mean) the Rao-Scott design effect of the transition table was about 1.5
#' ([mosaic_test()]). The printed values, which take the sequence as the
#' unit, were an ICC of 0.05 and a design effect of 1.24; resampling
#' students gave intervals about a quarter wider than resampling sessions,
#' and fewer significant edges (18 against 22).
#'
#' @param x The network: a `netobject` that carries its data (any network from
#'   [build_network()]), or a `cograph_network`, which is converted first. A
#'   `netobject_group` or `mcml` bootstraps every network it holds (for an
#'   `mcml`, the macro network and each cluster's), and a `wtna_mixed`
#'   bootstraps both of its components (see *Value*). A `tna` model or a
#'   `group_tna` (package tna) is first rebuilt as a network from its own
#'   sequences (see [as_netobject()]). `actor` is read from the metadata the
#'   model keeps (`attr(model$data, "meta")`, kept when the model is built
#'   from `tna::prepare_data()`); a model without it raises
#'   `nestimate_actor_missing`. A network without its data raises an error,
#'   and one that rests on a single sequence, or one-hot data with no actor or
#'   session, has nothing to resample and raises an error of class
#'   `nestimate_single_sequence`. A sequence network with no observed
#'   transition raises an error of class `nestimate_no_transition`, and a
#'   network built with `level` one of class
#'   `nestimate_multilevel_resampling`. For the `print()` methods: the result of
#'   `bootstrap_network()`.
#' @param iter The number of bootstrap replicates, at least 2. Default
#'   `1000`.
#' @param ci_level The significance level, between 0 and 1: the confidence
#'   intervals cover `1 - ci_level`, and an edge is significant when its
#'   p-value is below `ci_level`. Default `0.05`.
#' @param inference The test behind the p-values. `"stability"` (default)
#'   counts the replicates outside the consistency range around the
#'   original weight; `"threshold"` counts the replicates whose absolute
#'   weight is below `edge_threshold`.
#' @param consistency_range For `inference = "stability"`: the two
#'   multipliers of the original weight that bound the consistency range.
#'   Default `c(0.75, 1.25)`. For a negative weight the bounds are swapped,
#'   so the range still runs from the lower to the upper value. Not used by
#'   `"threshold"`.
#' @param edge_threshold For `inference = "threshold"`: the absolute weight
#'   a replicate must reach. Default `NULL`, the 10th percentile of the
#'   non-zero absolute original weights. Not used by `"stability"`.
#' @param seed A number that makes the resamples reproducible. The session's
#'   own random stream is left as it was. For a `netobject_group` or `mcml`, each network
#'   is bootstrapped from the same seed; for a `wtna_mixed`, the seed is set
#'   once, before both components. Default `NULL`.
#' @param boundary For `inference = "stability"`: how a replicate that lands
#'   exactly on a bound of the consistency range counts. `"inclusive"`
#'   (default, as in tna) counts it as outside (\eqn{\le} / \eqn{\ge});
#'   `"strict"` counts only replicates strictly outside (\eqn{<} /
#'   \eqn{>}).
#' @param ci_method The confidence interval for the edge weights.
#'   `"percentile"` (default) uses the empirical bootstrap quantiles
#'   (Efron). `"basic"` reflects those quantiles around the observed weight,
#'   \eqn{(2\hat{\theta} - q_{1-\alpha/2}, 2\hat{\theta} - q_{\alpha/2})}
#'   (Davison & Hinkley 1997, eq. 5.6), which corrects first-order bootstrap
#'   bias but can give bounds outside the natural weight range near its
#'   limits (for example below 0 for a transition probability close to 0).
#' @param actor The column identifying the actor each sequence belongs to
#'   (for example `"student_id"` for sessions nested in students, `"Group"`
#'   for students nested in teams), looked up in the network's metadata or
#'   its wide sequence data. With it, whole actors are resampled; see
#'   *Nested data and `actor`*. Several column names, outermost first (such
#'   as `c("Group", "Student")`), are nested levels, of which the outermost
#'   is resampled; the reported number of actors, ICC and design effect are
#'   those of the outermost level, and levels that are not nested outermost
#'   first raise an error of class `nestimate_bad_actor`. Default `NULL`:
#'   sequences are resampled one
#'   by one. It is available for sequence-based transition networks
#'   (`"relative"`, `"frequency"`, `"co_occurrence"`); a method outside the
#'   fast path (see *Details*) raises an error of class
#'   `nestimate_actor_unsupported`. A value that is not one or more distinct
#'   column names, or whose outermost column names fewer than two actors, raises an error of class
#'   `nestimate_bad_actor`; a column that is absent or has missing values,
#'   one of class `nestimate_actor_missing`; and ids that do not give one
#'   actor per sequence, one of class `nestimate_actor_misaligned`.
#'
#' @return A `net_bootstrap` object. `print()` lists the five strongest
#'   significant edges (by absolute bootstrap mean, with their interval)
#'   and reports the iterations, the number of significant edges, the
#'   interval level and the test, and with `actor` the ICC and design
#'   effect; [summary()] returns one row per edge of the original network;
#'   `cograph::splot()` draws the bootstrapped network; [net_prune()] with
#'   `method = "bootstrap"` removes the non-significant edges from the
#'   original network. The fields:
#' \describe{
#'   \item{original}{The original `netobject`.}
#'   \item{mean, sd}{The bootstrap mean and standard deviation of each
#'     weight, as matrices.}
#'   \item{p_values}{The matrix of p-values.}
#'   \item{significant}{The original weights where the p-value is below
#'     `ci_level`, and 0 elsewhere, self-loops included.}
#'   \item{ci_lower, ci_upper}{The matrices of interval bounds.}
#'   \item{cr_lower, cr_upper}{The matrices of consistency-range bounds
#'     (used by `"stability"` only).}
#'   \item{summary}{A data frame with one row per non-zero original edge
#'     (one row per unordered pair in an undirected network) and the
#'     columns `from`, `to`, `weight`, `mean`, `sd`, `p_value`, `sig`,
#'     `ci_lower` and `ci_upper`, plus `cr_lower` and `cr_upper` with
#'     `inference = "stability"`.}
#'   \item{model}{The pruned `netobject`, whose weights are `significant`:
#'     its non-significant edges, self-loops included, are zero. (The
#'     network [net_prune()] returns keeps every self-loop.)}
#'   \item{method, params, iter, ci_level, inference, ci_method}{The
#'     bootstrap settings.}
#'   \item{consistency_range, edge_threshold}{The test settings.}
#'   \item{actor, n_actors}{The `actor` column (or columns, outermost first)
#'     and the number of actors resampled (at the outermost level); `NULL`
#'     without `actor`.}
#'   \item{clustering}{Only with `actor`: a one-row data frame of
#'     `n_sequences`, `n_actors`, `icc`, `icc_ci_lower`, `icc_ci_upper` and
#'     `deff_edges`.}
#'   \item{clustering_edges}{Only with `actor`: one row per edge of
#'     `summary`, with `from`, `to`, `icc`, `sd_actor`, `sd_sequence` and
#'     `deff`.}
#' }
#' A `netobject_group` or `mcml` gives a `net_bootstrap_group`, a named list
#' of `net_bootstrap` results: its `print()` reports each network's
#' significant and total edges and the significant edges all networks
#' share (when the networks have edges in common); its `summary()` is
#' described with the methods below. A `wtna_mixed` gives a
#' `wtna_boot_mixed` holding the `transition` and `cooccurrence` results.
#'
#' @examples
#' net <- build_network(data.frame(V1 = c("A","B","C"), V2 = c("B","C","A")),
#'   method = "relative")
#' boot <- bootstrap_network(net, iter = 10)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' boot <- bootstrap_network(net, iter = 100)
#' print(boot)
#' summary(boot)
#'
#' # Students nested in teams: resample whole teams
#' teams <- build_network(group_regulation_long, method = "relative",
#'                        actor = "Actor", action = "Action", time = "Time",
#'                        group = "Achiever")
#' bootstrap_network(teams, iter = 100, actor = "Group", seed = 1)
#' }
#'
#' @references
#' Davison, A. C., & Hinkley, D. V. (1997). \emph{Bootstrap Methods and
#' Their Application}. Cambridge University Press.
#'
#' Field, C. A., & Welsh, A. H. (2007). Bootstrapping clustered data.
#' \emph{Journal of the Royal Statistical Society: Series B}, 69(3),
#' 369-390.
#'
#' Kish, L. (1965). \emph{Survey Sampling}. Wiley.
#'
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
#' assessing rater reliability. \emph{Psychological Bulletin}, 86(2),
#' 420-428.
#'
#' @seealso \code{\link{certainty}} for the closed-form Bayesian counterpart
#'   (same result layout, no resampling);
#'   \code{\link{build_network}}
#'
#' @importFrom stats quantile sd
#' @export
bootstrap_network <- function(x,
                              iter = 1000L,
                              ci_level = 0.05,
                              inference = "stability",
                              consistency_range = c(0.75, 1.25),
                              edge_threshold = NULL,
                              seed = NULL,
                              boundary = c("inclusive", "strict"),
                              ci_method = c("percentile", "basic"),
                              actor = NULL) {
  # a tna model or group_tna is rebuilt as a Nestimate network (R/tna_input.R)
  x <- .tna_input(x, "bootstrap_network()", need_data = TRUE, actor = actor)
  boundary <- match.arg(boundary)
  ci_method <- match.arg(ci_method)

  # ---- wtna_mixed dispatch: bootstrap both components ----
  if (inherits(x, "wtna_mixed")) {
    if (!is.null(seed)) .local_seed(seed)
    result <- list(
      transition   = bootstrap_network(x$transition,   iter = iter,
                                       ci_level = ci_level, inference = inference,
                                       consistency_range = consistency_range,
                                       edge_threshold = edge_threshold,
                                       boundary = boundary,
                                       ci_method = ci_method, actor = actor),
      cooccurrence = bootstrap_network(x$cooccurrence, iter = iter,
                                       ci_level = ci_level, inference = inference,
                                       consistency_range = consistency_range,
                                       edge_threshold = edge_threshold,
                                       boundary = boundary,
                                       ci_method = ci_method, actor = actor)
    )
    class(result) <- "wtna_boot_mixed"
    return(result)
  }

  # ---- mcml dispatch: convert to netobject_group via as_tna ----
  if (inherits(x, "mcml")) {
    x <- as_tna(x)
  }

  # ---- netobject_group dispatch: bootstrap each element ----
  .stop_multilevel_resampling(list(x), "bootstrap_network()")
  if (inherits(x, "netobject_group")) {
    .stop_single_unit_group(x)
    # the experimental note of a latency bootstrap is said once for the group
    said <- FALSE
    results <- withCallingHandlers(lapply(x, function(net) {
      bootstrap_network(net, iter = iter, ci_level = ci_level,
                        inference = inference,
                        consistency_range = consistency_range,
                        edge_threshold = edge_threshold, seed = seed,
                        boundary = boundary, ci_method = ci_method,
                        actor = actor)
    }), nestimate_experimental = function(m) {
      if (said) invokeRestart("muffleMessage")
      said <<- TRUE
    })
    class(results) <- c("net_bootstrap_group", "list")
    # every group's edges in one table, kept beside the groups (an attribute,
    # so the object still iterates as its groups); summary() returns it
    attr(results, "summary") <- .bootstrap_group_table(results)
    return(results)
  }

  # Accept netobject or cograph_network
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)
  if (!inherits(x, "netobject")) {
    stop(errorCondition(
      "'x' must be a netobject from build_network().",
      class = "nestimate_bad_input", call = NULL)) # nocov
  }
  if (is.null(x$data)) {
    stop(errorCondition(
      "netobject does not contain $data. Rebuild with build_network().",
      class = "nestimate_bad_input", call = NULL))
  }

  # Note: edgelist-derived data has no actor / session grouping, so
  # row-resampling uses individual transitions as the unit -- which may
  # under-represent within-actor correlation. Case-dropping and actor-
  # level sequence reconstruction are both more robust alternatives.
  if (identical(attr(x$data, "source"), "edgelist")) {
    warning(warningCondition(paste0(
      "Bootstrapping a network built from edgelist input (no per-actor ",
      "sequence data). Each transition is treated as an independent case, ",
      "which has limitations when transitions are correlated within ",
      "actors or sessions.\n",
      "Better alternatives:\n",
      "  * Rebuild as a sequence-based network: build_network(data, ",
      "method = \"tna\", actor = <...>, action = <...>, time = <...>) ",
      "then feed the result to build_mcml(). This reconstructs per-actor ",
      "sequences (via prepare()) so bootstrap and case-dropping operate at ",
      "the actor level natively.\n",
      "  * Use case-dropping rather than bootstrap:\n",
      "      - casedrop_reliability(x)  - model-level edge-weight reliability\n",
      "      - centrality_stability(x)  - per-centrality CS coefficient\n",
      "  * permutation(x, y) for comparing two networks."),
      class = "nestimate_edgelist_resampling", call = NULL))
  }

  data <- x$data
  method <- x$method
  params <- x$params
  scaling <- x$scaling
  threshold <- x$threshold
  level <- x$level

  # ---- Input validation ----
  .assert(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(ci_level), length(ci_level) == 1, ci_level > 0, ci_level < 1,
    is.character(inference), length(inference) == 1,
    is.numeric(consistency_range), length(consistency_range) == 2
  )
  iter <- as.integer(iter)
  inference <- match.arg(inference, c("stability", "threshold"))
  .check_permutation_block(actor, paired = FALSE)

  if (!is.null(seed)) {
    .assert(is.numeric(seed), length(seed) == 1)
    .local_seed(seed)
  }

  # ---- Resolve method and get estimator ----
  method <- .resolve_method_alias(method)
  estimator <- get_estimator(method)
  directed <- estimator$directed

  # ---- Original network is x itself ----
  original <- x
  states <- original$nodes$label
  n_states <- length(states)

  .stop_single_unit(x)

  # ---- Dispatch bootstrap ----
  # Fast precomputed path only works when data has usable columns
  has_data <- is.data.frame(data) && ncol(data) > 0L && nrow(data) > 0L
  count_method <- .resampling_method(method, params)
  counts_relative <- .counts_are_relative(count_method, params)
  # An mcml layer whose pairs were combined by a mean, median, maximum,
  # minimum, geometric mean or density (as_tna() records it) is not a sum of
  # per-unit counts: each resample combines its own rows the same way.
  aggregation <- .param_get(params, "aggregation", "sum")
  latency_boot <- NULL
  if (identical(method, "latency")) {
    if (!is.null(actor)) {
      stop(errorCondition(paste0(
        "`actor` is not supported for a latency network: its sequences are ",
        "resampled. To resample within groups, build it with `group` and ",
        "bootstrap the grouped network."),
        class = "nestimate_actor_unsupported", call = NULL))
    }
    .latency_experimental_message()
    latency_boot <- .bootstrap_latency(x, iter)
    boot_matrices <- latency_boot$weights
  } else if (!identical(aggregation, "sum") && has_data &&
      count_method %in% c("relative", "frequency")) {
    if (!is.null(actor)) {
      stop(errorCondition(
        sprintf(paste0("`actor` is not supported for a network whose ",
                       "transitions were combined by method = \"%s\"."),
                aggregation),
        class = "nestimate_actor_unsupported", call = NULL))
    }
    boot_matrices <- .bootstrap_aggregated(
      .resampling_transition_data(data, x$metadata, params),
      method = count_method, params = params, states = states,
      scaling = scaling, threshold = threshold, iter = iter,
      is_relative = counts_relative,
      include_zeros = .scaling_includes_zeros(method)
    )
  } else if (count_method %in% c("relative", "frequency", "co_occurrence")) {
    if (!has_data) {
      stop(errorCondition(
        paste0(
          "Bootstrap requires the original data stored in the netobject. ",
          "For wtna/cna networks, use wtna() directly instead of ",
          "build_network(method='cna')."),
        class = "nestimate_bad_input", call = NULL))
    }
    transition_data <- .resampling_transition_data(data, x$metadata, params)
    trans_2d <- .precompute_per_sequence(transition_data, count_method, params,
                                         states)
    .stop_no_transition(trans_2d, "bootstrap_network()")
    # Several actor columns are nested levels: every level is checked, and
    # only the outermost is resampled (Davison & Hinkley, 1997, section 3.8:
    # resampling inside the clusters as well overstates the variance).
    if (length(actor) > 1L) {
      .check_actor_levels(lapply(actor, function(col) {
        .permutation_block_ids(x, col, nrow(trans_2d), "x")
      }), actor)
    }
    block_ids <- if (is.null(actor)) NULL else
      .permutation_block_ids(x, actor[[1L]], nrow(trans_2d), "x")
    # Cluster bootstrap: resample whole units (summed counts), top level only
    units <- if (is.null(block_ids)) trans_2d else
      rowsum(trans_2d, block_ids, reorder = FALSE)
    if (!is.null(block_ids) && nrow(units) < 2L) {
      stop(errorCondition(
        sprintf("`actor = \"%s\"` has %d actor; resampling needs at least 2.",
                actor[[1L]], nrow(units)),
        class = "nestimate_bad_actor", call = NULL))
    }
    # the network's own scaling domain: a wtna network resampled as counts
    # still scales its non-zero cells only, as build_network() did
    boot_matrices <- .bootstrap_counts(
      units, method = count_method, states = states, scaling = scaling,
      threshold = threshold, iter = iter, is_relative = counts_relative,
      include_zeros = .scaling_includes_zeros(method)
    )
    # sequence-level reference drawn in the same run, for the design effect
    reference <- if (is.null(block_ids)) NULL else .bootstrap_counts(
      trans_2d, method = count_method, states = states, scaling = scaling,
      threshold = threshold, iter = iter, is_relative = counts_relative,
      include_zeros = .scaling_includes_zeros(method)
    )
  } else {
    .stop_block_unsupported(actor, method)
    .stop_no_transition_refit(x, method, "bootstrap_network()")
    boot_matrices <- .bootstrap_association(
      data = data, estimator = estimator, params = params, states = states,
      scaling = scaling, threshold = threshold, iter = iter,
      level = level, id_col = params$id,
      include_zeros = .scaling_includes_zeros(method)
    )
  }

  # ---- Auto edge_threshold for threshold inference ----
  if (inference == "threshold" && is.null(edge_threshold)) {
    abs_weights <- abs(as.vector(original$weights))
    nz_weights <- abs_weights[abs_weights > 0]
    edge_threshold <- if (length(nz_weights) > 0) {
      quantile(nz_weights, probs = 0.10)
    } else {
      0
    }
  }

  # ---- Compute statistics ----
  stats <- .compute_bootstrap_stats(
    boot_matrices = boot_matrices,
    original_matrix = original$weights,
    states = states,
    directed = directed,
    iter = iter,
    ci_level = ci_level,
    inference = inference,
    consistency_range = consistency_range,
    edge_threshold = edge_threshold,
    boundary = boundary,
    ci_method = ci_method
  )
  # a mean latency network is tested on its speed score (the weight with the
  # floor taken off), and needs replies in several sequences
  if (!is.null(latency_boot$scores) && identical(inference, "stability")) {
    stats <- .latency_stability(stats, latency_boot, original$weights,
                                params$latency_floor, consistency_range,
                                boundary, ci_level)
  }

  # ---- Build summary data frame ----
  summary_df <- .build_bootstrap_summary(
    stats = stats,
    original_matrix = original$weights,
    states = states,
    directed = directed,
    ci_level = ci_level,
    inference = inference
  )
  if (!is.null(latency_boot)) {
    edge <- cbind(summary_df$from, summary_df$to)
    summary_df$n_replies <- original$n_replies[edge]
    summary_df$n_sequences <- matrix(latency_boot$n_sequences, n_states,
                                     n_states, dimnames = dimnames(original$weights))[edge]
    if (!is.null(stats$cv)) summary_df$cv <- stats$cv[edge]
  }

  # ---- Build pruned model ----
  pruned_matrix <- stats$significant
  pruned_edges <- .extract_edges_from_matrix(pruned_matrix, directed = directed)

  model <- list(
    weights = pruned_matrix,
    nodes = original$nodes,
    edges = pruned_edges,
    directed = directed,
    method = method,
    params = params,
    scaling = scaling,
    threshold = threshold,
    n_nodes = n_states,
    n_edges = nrow(pruned_edges),
    level = level,
    meta = c(list(source = "nestimate", layout = NULL,
                  tna = list(method = method)), .splot_contract(method)),
    node_groups = NULL
  )
  class(model) <- c("netobject", "cograph_network")

  # ---- Assemble result ----
  result <- list(
    original          = original,
    mean              = stats$mean,
    sd                = stats$sd,
    p_values          = stats$p_values,
    significant       = stats$significant,
    ci_lower          = stats$ci_lower,
    ci_upper          = stats$ci_upper,
    cr_lower          = stats$cr_lower,
    cr_upper          = stats$cr_upper,
    summary           = summary_df,
    model             = model,
    method            = method,
    params            = params,
    iter              = iter,
    ci_level          = ci_level,
    inference         = inference,
    consistency_range = consistency_range,
    edge_threshold    = edge_threshold,
    ci_method         = ci_method,
    actor             = actor,
    n_actors          = if (is.null(actor)) NULL else nrow(units)
  )
  if (!is.null(latency_boot)) {
    result$latency_cap <- latency_boot$cap
    # drawn as the network itself is (cograph's rendering contract)
    result$meta <- .splot_contract(method)
  }
  if (!is.null(actor)) {
    clus <- .build_bootstrap_clustering(stats$sd, reference, trans_2d,
                                        block_ids, summary_df, states)
    result$clustering <- clus$overall
    result$clustering_edges <- clus$edges
  }
  class(result) <- "net_bootstrap"
  result
}


#' Tidy nesting diagnostics of a blocked bootstrap
#'
#' Design effect = variance over unit-level replicates / variance over
#' sequence-level replicates (Kish, 1965); ICC as in permutation().
#' @noRd
.build_bootstrap_clustering <- function(sd_blocked, reference, trans_2d,
                                        block_ids, summary_df, states) {
  n_states <- length(states)
  sd_unblocked <- matrix(apply(reference, 2, stats::sd, na.rm = TRUE),
                         n_states, n_states, dimnames = list(states, states))
  icc <- .block_icc(trans_2d, block_ids, rep("x", length(block_ids)))
  icc_edges <- matrix(icc$per_edge, n_states, n_states, byrow = TRUE,
                      dimnames = list(states, states))
  at <- cbind(summary_df$from, summary_df$to)
  edges <- data.frame(
    from = summary_df$from, to = summary_df$to,
    icc = icc_edges[at],
    sd_actor = sd_blocked[at], sd_sequence = sd_unblocked[at],
    deff = (sd_blocked[at] / sd_unblocked[at])^2,
    stringsAsFactors = FALSE
  )
  overall <- data.frame(
    n_sequences = nrow(trans_2d),
    n_actors = length(unique(block_ids)),
    icc = icc$estimate,
    icc_ci_lower = icc$ci[1L], icc_ci_upper = icc$ci[2L],
    deff_edges = stats::median(edges$deff[is.finite(edges$deff)])
  )
  list(overall = overall, edges = edges)
}


#' Reattach grouping columns needed by transition resampling
#' @noRd
.resampling_transition_data <- function(data, metadata, params) {
  if (!is.data.frame(data) || !is.data.frame(metadata)) return(data)
  needed <- unique(c(params$actor, params$id, params$id_col))
  needed <- needed[!is.na(needed) & nzchar(needed)]
  add <- setdiff(intersect(needed, names(metadata)), names(data))
  if (length(add) == 0L) return(data)
  .keep_window_attrs(cbind(metadata[, add, drop = FALSE], data), data)
}

# A network estimated from binary indicator columns, by build_network() or wtna().
#' @noRd
.is_onehot_params <- function(params) {
  !is.null(.param_get(params, "codes")) ||
    identical(.param_get(params, "format"), "onehot")
}

# .prepare_onehot() records its windows as attributes, and selecting columns
# drops them. Without them the positions read as one unwindowed sequence.
#' @noRd
.keep_window_attrs <- function(subset, source) {
  kept <- intersect(names(attributes(source)),
                    c("windowed", "window_size", "window_span", "codes"))
  attributes(subset)[kept] <- attributes(source)[kept]
  subset
}


#' The count-based method a network is resampled by
#'
#' wtna() networks hold one-hot counts, the same as frequency, relative and
#' co_occurrence do, so their units are resampled the same way. Refitting on
#' shuffled rows would break the time order their transitions are counted on.
#' @param method The resolved method name.
#' @param params The network's stored parameters.
#' @return A method name.
#' @noRd
.resampling_method <- function(method, params) {
  if (!.is_onehot_params(params)) return(method)
  switch(method,
    wtna_cooccurrence = "co_occurrence",
    wtna = if (identical(.param_get(params, "type"), "relative")) "relative" else
      "frequency",
    method
  )
}

# One-hot co-occurrence can also be row-normalised, by `type = "relative"`.
#' @noRd
.counts_are_relative <- function(method, params) {
  method == "relative" ||
    (method == "co_occurrence" &&
       identical(.param_get(params, "type"), "relative"))
}

#' Refuse to resample a sequence network that holds a single unit
#'
#' Resampling draws whole sequences (or, for one-hot data, whole actors or
#' sessions). With one there is nothing to draw: every resample is the data,
#' which gives intervals of zero width and p-values of 1 that read as certainty.
#' @param x A netobject.
#' @param label Group name for the message, or NULL.
#' @noRd
.stop_single_unit <- function(x, label = NULL) {
  data <- x$data
  if (!is.data.frame(data) || is.null(x$method)) return(invisible(NULL))
  params <- x$params
  method <- .resampling_method(.resolve_method_alias(x$method), params)
  sequence_based <- method %in% c("relative", "frequency", "co_occurrence",
                                  "attention", "ngram", "gap", "reverse") ||
    .is_onehot_params(params) || isTRUE(attr(data, "windowed"))
  if (!sequence_based) return(invisible(NULL))
  n_units <- .transition_resampling_n_sequences(
    .resampling_transition_data(data, x$metadata, params), params
  )
  if (!is.na(n_units) && n_units == 1L) {
    stop(errorCondition(
      paste0(if (is.null(label)) "" else paste0(label, ": "),
             "This network rests on one sequence, so there is nothing to ",
             "resample. Rebuild it with `actor` or `session` set to the column ",
             "that identifies the units."),
      class = c("nestimate_single_sequence", "nestimate_bad_input"),
      call = NULL
    ))
  }
  invisible(NULL)
}

# A sequence network with no observed transition (no sequence has two
# adjacent observed states) gives all-zero resamples and zero-width
# intervals; it is refused, as transition_entropy() and the landscapes do.
#' @param counts Per-unit transition counts (.precompute_per_sequence()).
#' @param verb The verb, for the message.
#' @noRd
.stop_no_transition <- function(counts, verb) {
  if (is.null(counts) || sum(counts) > 0) return(invisible(NULL))
  stop(errorCondition(sprintf(paste0(
    "%s: no transition is observed (no sequence has two adjacent observed ",
    "states), so there is nothing to resample."), verb),
    class = c("nestimate_no_transition", "nestimate_bad_input"),
    call = NULL))
}

# attention, ngram, gap and reverse networks are resampled by refitting, so
# there are no per-unit counts to sum: the network's own estimator, run once
# on its stored data, says whether any transition is observed.
#' @param net A netobject.
#' @param method Its resolved method.
#' @param verb The verb, for the message.
#' @noRd
.stop_no_transition_refit <- function(net, method, verb) {
  if (!method %in% c("attention", "ngram", "gap", "reverse") ||
      !is.data.frame(net$data)) {
    return(invisible(NULL))
  }
  fit <- do.call(get_estimator(method)$fn,
                 c(list(data = net$data), net$params))
  .stop_no_transition(fit$matrix, verb)
}

#' @noRd
.stop_single_unit_group <- function(x) {
  invisible(Map(.stop_single_unit, x, label = names(x)))
}

#' Note that a transition network rests on one sequence
#'
#' Building one is legitimate; resampling it is refused by
#' \code{.stop_single_unit()}. Class \code{nestimate_single_sequence}.
#' @noRd
.single_sequence_note <- function() {
  # base R has no messageCondition(); add the class to a simpleMessage
  cond <- simpleMessage(paste0(
    "A network with one long sequence is not recommended and can't be ",
    "validated using bootstrap and other confirmatory testings.\n"))
  class(cond) <- c("nestimate_single_sequence", class(cond))
  message(cond)
}


#' @noRd
.transition_resampling_n_sequences <- function(data, params) {
  if (!is.data.frame(data)) {
    return(NA_integer_)
  }
  concat <- .wide_resampling_concat(params)
  if (concat > 1L) return(as.integer(ceiling(nrow(data) / concat)))
  # The units the resampler draws: for one-hot data, the actors (or
  # sessions) of `params$actor` (.precompute_per_sequence_onehot()), and
  # without one, the rows are time points of one sequence; otherwise every
  # row is a sequence, whatever id it carries (the estimators count each
  # row as its own sequence and the resampler draws rows).
  if (!.is_onehot_params(params)) return(nrow(data))
  actor <- .param_get(params, "actor")
  if (is.null(actor)) return(1L)
  if (!all(actor %in% names(data))) return(nrow(data))
  if (length(actor) == 1L) return(length(unique(data[[actor]])))
  length(unique(interaction(data[, actor, drop = FALSE], drop = TRUE)))
}


#' Number of wide rows merged into one resampling sequence
#'
#' Onehot and long estimators do not use the wide-row concat option.
#' @param params Stored estimator parameters.
#' @return A positive integer.
#' @noRd
.wide_resampling_concat <- function(params) {
  if (!is.null(.param_get(params, "codes")) ||
      .param_get(params, "format", "auto") %in% c("onehot", "long")) {
    return(1L)
  }
  concat <- as.integer(.param_get(params, "concat", 1L))
  .assert("'concat' must be a positive integer." =
              length(concat) == 1L && !is.na(concat) && concat >= 1L)
  concat
}


# ---- Transition fast path ----

#' Bootstrap transition networks via pre-computed per-sequence counts
#' @noRd
.bootstrap_transition <- function(data, method, params, states,
                                  scaling, threshold, iter) {
  counts <- .precompute_per_sequence(data, method, params, states)
  .stop_no_transition(counts, "bootstrap_network()")
  .bootstrap_counts(counts,
                    method = method, states = states, scaling = scaling,
                    threshold = threshold, iter = iter)
}


#' Bootstrap from a precomputed unit x (n_states^2) count matrix
#'
#' Rows are the resampling units: sequences, or whole blocks (summed
#' counts) for the cluster bootstrap.
#' @noRd
.bootstrap_counts <- function(trans_2d, method, states, scaling, threshold,
                              iter, is_relative = method == "relative",
                              include_zeros = .scaling_includes_zeros(method)) {
  n_states <- length(states)
  nbins <- n_states * n_states
  n_seq <- nrow(trans_2d)

  # vapply: each iteration resamples sequences + sums + post-processes
  boot_flat <- vapply(seq_len(iter), function(i) {
    idx <- sample.int(n_seq, n_seq, replace = TRUE)
    boot_counts <- colSums(trans_2d[idx, , drop = FALSE])
    # Inline post-processing (no dimnames -- we flatten immediately)
    mat <- matrix(boot_counts, n_states, n_states, byrow = TRUE)
    .finish_replicate(mat, is_relative, scaling, include_zeros, threshold)
  }, numeric(nbins))

  # With one state (nbins = 1) vapply returns a vector; keep it nbins x iter
  # so the result is always iter x nbins.
  t(matrix(boot_flat, nrow = nbins))
}


#' Row-normalise, scale and threshold one replicate, flattened column-major
#' @noRd
.finish_replicate <- function(mat, is_relative, scaling, include_zeros,
                              threshold) {
  if (is_relative) {
    rs <- rowSums(mat)
    nz <- rs > 0
    mat[nz, ] <- mat[nz, ] / rs[nz]
  }
  if (!is.null(scaling)) {
    mat <- .apply_scaling(mat, scaling, include_zeros = include_zeros)
  }
  if (threshold > 0) mat[abs(mat) < threshold] <- 0
  as.vector(mat)
}


#' Bootstrap an mcml layer whose transitions were combined other than by sum
#'
#' build_mcml() combines the weights of the transitions between two clusters
#' (or two states) with net_aggregate_weights(): a mean, median, maximum,
#' minimum, geometric mean or density. Those are not sums of per-unit
#' counts, so each resample combines the weights of the transitions its
#' drawn units hold, then applies what the fit applied after aggregating
#' (`symmetrize`, row normalisation for type = "tna"). Every transition of a
#' unit carries the unit's weight: an edge-list row holds one transition of
#' its row weight, a sequence holds transitions of weight 1.
#' @param data,method,params,states As for .precompute_per_sequence().
#' @param scaling,threshold,iter,is_relative,include_zeros As for
#'   .bootstrap_counts().
#' @return An iter x n_states^2 matrix of flattened (column-major) replicates.
#' @noRd
.bootstrap_aggregated <- function(data, method, params, states, scaling,
                                  threshold, iter, is_relative,
                                  include_zeros) {
  aggregation <- .param_get(params, "aggregation")
  n_states <- length(states)
  nbins <- n_states * n_states
  counts <- .precompute_per_sequence_counts(data, method, params, states)
  row_weight <- .param_get(params, "row_weight")
  unit_weight <- if (is.null(row_weight)) rep(1, nrow(counts)) else
    data[[row_weight]]
  n_possible <- .param_get(params, "n_possible")
  # cells are flattened row-major in `counts`
  cell_possible <- if (is.null(n_possible)) rep(1, nbins) else
    as.vector(t(n_possible[states, states, drop = FALSE]))
  # one record per (unit, cell) holding transitions of a non-zero weight;
  # net_aggregate_weights() leaves zero weights out
  counts[unit_weight == 0, ] <- 0
  at <- which(counts > 0, arr.ind = TRUE)
  rec_unit <- at[, 1L]
  rec_n <- counts[at]
  rec_w <- unit_weight[rec_unit]
  by_cell <- split(seq_along(rec_unit), at[, 2L])
  cells <- as.integer(names(by_cell))
  n_units <- nrow(counts)
  symmetrize <- isTRUE(.param_get(params, "symmetrize", FALSE))

  boot_flat <- vapply(seq_len(iter), function(i) {
    drawn <- tabulate(sample.int(n_units, n_units, replace = TRUE), n_units)
    mult <- drawn[rec_unit] * rec_n
    flat <- numeric(nbins)
    flat[cells] <- vapply(seq_along(cells), function(j) {
      r <- by_cell[[j]]
      .aggregate_repeated(rec_w[r], mult[r], aggregation,
                          cell_possible[cells[j]])
    }, numeric(1L))
    mat <- matrix(flat, n_states, n_states, byrow = TRUE)
    if (symmetrize) mat <- (mat + t(mat)) / 2
    .finish_replicate(mat, is_relative, scaling, include_zeros, threshold)
  }, numeric(nbins))
  t(matrix(boot_flat, nrow = nbins))
}


#' net_aggregate_weights() of values repeated `times` times
#'
#' Equals `net_aggregate_weights(rep(values, times), method, n_possible)`
#' for non-negative values, without building the repeated vector.
#' @noRd
.aggregate_repeated <- function(values, times, method, n_possible = 1) {
  keep <- times > 0 & values != 0
  values <- values[keep]
  times <- times[keep]
  if (length(values) == 0L) return(0)
  switch(method,
    mean = sum(values * times) / sum(times),
    density = sum(values * times) / n_possible,
    max = max(values),
    min = min(values),
    geomean = exp(sum(log(values) * times) / sum(times)),
    median = {
      ord <- order(values)
      values <- values[ord]
      ends <- cumsum(times[ord])
      total <- ends[length(ends)]
      at <- function(pos) values[which(ends >= pos)[1L]]
      (at(floor((total + 1) / 2)) + at(floor(total / 2) + 1)) / 2
    },
    stop(errorCondition(paste0("Unknown aggregation: ", method),
                        class = "nestimate_bad_input", call = NULL))
  )
}


#' Pre-compute per-sequence count matrices
#'
#' Returns an n_seq x (n_states^2) matrix where each row holds the
#' flattened transition (or co-occurrence) counts for one sequence.
#'
#' Only wide-format data is supported for the bootstrap fast path.
#' For long-format data, convert first using
#' \code{\link{convert_sequence_format}}.
#' @noRd
.precompute_per_sequence <- function(data, method, params, states) {
  # Per-unit counts add up to the network only when its pairs were summed. An
  # mcml layer combined by a mean, maximum, density, ... (as_tna() records it
  # as params$aggregation) is resampled by bootstrap_network() alone, which
  # combines each resample the same way; every other reader would treat it
  # as summed counts.
  aggregation <- .param_get(params, "aggregation", "sum")
  if (!identical(aggregation, "sum")) {
    stop(errorCondition(sprintf(paste0(
      "This network combines its transitions by method = \"%s\", so its ",
      "weights are not counts; of the resampling verbs only ",
      "bootstrap_network() reads it. Build the mcml with method = \"sum\" ",
      "for the others."), aggregation),
      class = "nestimate_bad_input", call = NULL))
  }
  .declared_count_transform(
    .precompute_per_sequence_counts(data, method, params, states),
    data, params, states)
}


#' Per-unit weights and symmetrisation a network declares in its params
#'
#' A weighted edge list (`params$row_weight`, set by as_tna() for an mcml
#' built with `weight =`) counts each row by its weight, and a network whose
#' weights are (W + t(W)) / 2 of the counts W (`params$symmetrize`, set by
#' as_tna() for an undirected mcml) has each unit's counts symmetrised the
#' same way. Both are linear, so the summed units equal the fit.
#' @noRd
.declared_count_transform <- function(counts, data, params, states) {
  row_weight <- .param_get(params, "row_weight")
  if (!is.null(row_weight)) counts <- counts * data[[row_weight]]
  if (isTRUE(.param_get(params, "symmetrize", FALSE))) {
    n <- length(states)
    transposed <- as.vector(matrix(seq_len(n * n), n, n, byrow = TRUE))
    counts <- (counts + counts[, transposed, drop = FALSE]) / 2
  }
  counts
}


#' Per-unit counts as the estimator counts them, before any declared
#' weighting or symmetrisation (see .precompute_per_sequence())
#' @noRd
.precompute_per_sequence_counts <- function(data, method, params, states) {
  # Determine format (mirrors estimator dispatch)
  format <- .param_get(params, "format", "auto")
  action <- .param_get(params, "action")
  codes <- .param_get(params, "codes")

  # Before the format: a one-hot column may be called "Action".
  if (identical(format, "onehot") || !is.null(codes)) {
    return(.precompute_per_sequence_onehot(data, method, params, states))
  }

  if (format == "auto") {
    format <- .format_from_action(action, data)
  }

  if (format == "long") {
    stop(errorCondition(paste0(
      "Bootstrap fast path requires wide-format data. ",
      "Convert with convert_sequence_format() first."),
      class = "nestimate_bad_input", call = NULL))
  }

  id_col <- .param_get(params, "id") %||% .param_get(params, "id_col")
  cols <- params$cols
  if (isTRUE(attr(data, "windowed"))) {
    return(.precompute_per_sequence_windowed(data, method, params, states))
  }
  .precompute_per_sequence_wide(data, method, cols, id_col, states,
                                begin_state = .param_get(params, "begin_state"),
                                end_state = .param_get(params, "end_state"),
                                weighted = .param_get(params, "weighted", FALSE),
                                concat = .param_get(params, "concat", 1L))
}


#' Pre-compute per-sequence counts from windowed sequences
#'
#' The estimator's own counting on one sequence at a time, so a window means
#' the same in a resample as in the network being resampled. With
#' `params$concat = k`, the estimator merges each run of `k` rows (the last
#' run may be shorter) into one sequence before windowing; such a merged
#' sequence is one resampling unit, counted with the same `concat`.
#' @noRd
.precompute_per_sequence_windowed <- function(data, method, params, states) {
  count <- if (method == "co_occurrence") .count_cooccurrence_wide else
    .count_transitions_wide
  concat <- .wide_resampling_concat(params)
  units <- split(seq_len(nrow(data)), ceiling(seq_len(nrow(data)) / concat))
  rows <- lapply(units, function(ii) {
    mat <- count(data[ii, , drop = FALSE], cols = params$cols,
                 alphabet = states, concat = concat,
                 weighted = .param_get(params, "weighted", FALSE),
                 begin_state = .param_get(params, "begin_state"),
                 end_state = .param_get(params, "end_state"))
    as.vector(t(mat[states, states, drop = FALSE]))
  })
  matrix(unlist(rows), nrow = length(units), ncol = length(states)^2L,
         byrow = TRUE)
}


#' Pre-compute per-sequence counts from one-hot format
#' @noRd
.precompute_per_sequence_onehot <- function(data, method, params, states) {
  codes <- .param_get(params, "codes", states)
  actor <- .param_get(params, "actor")
  window_size <- .param_get(params, "window_size", 1L)
  mode <- .param_get(params, "mode", "non-overlapping")
  wtna_method <- if (method == "co_occurrence") "cooccurrence" else "transition"
  n_states <- length(states)
  nbins <- n_states * n_states

  # One unit per actor, in the order the actors first appear in the rows (not
  # the order of a factor's levels), so the i-th unit of two networks is
  # their i-th actor, as a paired permutation pairs them.
  groups <- if (is.null(actor)) {
    list(data)
  } else {
    key <- if (length(actor) == 1L) as.character(data[[actor]]) else
      as.character(interaction(data[, actor, drop = FALSE], drop = TRUE))
    split(data, factor(key, levels = unique(key)), drop = TRUE)
  }

  rows <- lapply(groups, function(g) {
    mat <- .estimator_wtna_core(
      g, codes = codes, window_size = window_size, mode = mode,
      actor = NULL, wtna_method = wtna_method, type = "frequency"
    )$matrix
    mat <- mat[states, states, drop = FALSE]
    as.vector(t(mat))
  })

  if (length(rows) == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = nbins))
  }
  do.call(rbind, rows)
}


#' Pre-compute per-sequence counts from wide format (transition + co-occurrence)
#'
#' Single function handling both consecutive-pair transitions
#' (relative/frequency) and all-column-pair co-occurrences.
#' @noRd
.precompute_per_sequence_wide <- function(data, method, cols, id_col, states,
                                          begin_state = NULL,
                                          end_state = NULL,
                                          weighted = FALSE,
                                          concat = 1L) {
  # The estimator's own sequence preparation, so start/end markers and void
  # markers are placed exactly as in the network being resampled.
  mat <- .prepare_sequence_wide(data, id = id_col, cols = cols,
                                alphabet = states, concat = concat,
                                begin_state = begin_state,
                                end_state = end_state)$matrix
  nc <- ncol(mat)
  nr <- nrow(mat)
  n_states <- length(states)
  nbins <- n_states * n_states

  # Integer-encode the entire matrix once
  int_mat <- matrix(match(mat, states), nrow = nr, ncol = nc)

  result <- if (method %in% c("relative", "frequency")) {
    # Consecutive column pairs: (col[i], col[i+1])
    from_mat <- int_mat[, -nc, drop = FALSE]
    to_mat <- int_mat[, -1L, drop = FALSE]

    row_ids <- rep(seq_len(nr), times = nc - 1L)
    from_vec <- as.vector(from_mat)
    to_vec <- as.vector(to_mat)

    valid <- !is.na(from_vec) & !is.na(to_vec)
    row_ids <- row_ids[valid]
    from_vec <- from_vec[valid]
    to_vec <- to_vec[valid]

    pair_idx <- (from_vec - 1L) * n_states + to_vec
    combined_idx <- (row_ids - 1L) * nbins + pair_idx

    counts <- tabulate(combined_idx, nbins = nr * nbins)
    matrix(as.numeric(counts), nrow = nr, ncol = nbins, byrow = TRUE)

  } else {
    # Co-occurrence: all column pairs (i, j) where i < j
    result <- matrix(0, nrow = nr, ncol = nbins)

    # Loop kept on purpose: it accumulates into one nr x nbins matrix per
    # column pair; an apply over the nc * (nc - 1) / 2 pairs would hold one
    # such matrix per pair before summing them.
    for (i in seq_len(nc - 1L)) {
      col_i <- int_mat[, i]
      for (j in seq(i + 1L, nc)) {
        col_j <- int_mat[, j]
        valid <- !is.na(col_i) & !is.na(col_j)
        fi <- col_i[valid]
        tj <- col_j[valid]
        rows_valid <- which(valid)

        # Bidirectional: forward + reverse
        idx_fwd <- (fi - 1L) * n_states + tj
        idx_rev <- (tj - 1L) * n_states + fi
        combined <- c(
          (rows_valid - 1L) * nbins + idx_fwd,
          (rows_valid - 1L) * nbins + idx_rev
        )

        tab <- tabulate(combined, nbins = nr * nbins)
        result <- result + matrix(tab, nrow = nr, ncol = nbins, byrow = TRUE)
      }
    }

    # Self-pairs are double-counted by bidirectional approach
    diag_indices <- (seq_len(n_states) - 1L) * n_states + seq_len(n_states)
    result[, diag_indices] <- result[, diag_indices] / 2

    result
  }

  # Match the estimator's per-sequence weighting, including the length of
  # concatenated sequences and any start/end markers. Empty sequences
  # contribute zero rather than 0 / 0.
  if (isTRUE(weighted)) {
    seq_lengths <- rowSums(!is.na(mat))
    result <- result / pmax(seq_lengths, 1L)
  }
  result
}


# ---- Association path ----

#' Bootstrap association networks via full estimator calls
#' @noRd
.bootstrap_association <- function(data, estimator, params, states,
                                   scaling, threshold, iter,
                                   level, id_col, include_zeros = FALSE) {
  n <- nrow(data)
  n_states <- length(states)
  nbins <- n_states * n_states

  boot_flat <- vapply(seq_len(iter), function(i) {
    boot_data <- data[sample.int(n, n, replace = TRUE), , drop = FALSE]

    # For multilevel: apply decomposition to bootstrap sample
    if (!is.null(level) && !is.null(id_col) && !estimator$directed) {
      boot_data <- tryCatch( # nocov start
        .decompose_multilevel(boot_data, id_col = id_col, level = level),
        error = function(e) NULL
      )
      if (is.null(boot_data)) return(rep(NA_real_, nbins)) # nocov end
    }

    est <- tryCatch(
      do.call(estimator$fn, c(list(data = boot_data), params)),
      error = function(e) NULL
    )
    if (is.null(est)) return(rep(NA_real_, nbins))

    mat <- est$matrix
    # Align to expected states order
    if (!identical(rownames(mat), states)) {
      common <- intersect(states, rownames(mat))
      if (length(common) < n_states) return(rep(NA_real_, nbins))
      mat <- mat[states, states] # nocov
    }

    if (!is.null(scaling)) { # nocov start
      mat <- .apply_scaling(mat, scaling, include_zeros = include_zeros)
    }
    if (threshold > 0) mat[abs(mat) < threshold] <- 0
    as.vector(mat) # nocov end
  }, numeric(nbins))

  # A replicate that could not be estimated, or that lost a node, is a
  # column of NA left out of the statistics: counted and said, as the other
  # resampling verbs do.
  .warn_failed_resamples(sum(colSums(!is.na(boot_flat)) == 0L), iter,
                         "bootstrap_network()", below_threshold = FALSE)

  t(boot_flat)
}


# ---- Vectorized statistics ----

#' Compute bootstrap statistics from the iter x nbins matrix
#' @noRd
.compute_bootstrap_stats <- function(boot_matrices, original_matrix, states,
                                     directed, iter, ci_level, inference,
                                     consistency_range, edge_threshold,
                                     boundary = "inclusive",
                                     ci_method = "percentile") {
  n_states <- length(states)
  orig_flat <- as.vector(original_matrix)

  # Mean and SD (vectorized)
  weights_mean <- colMeans(boot_matrices, na.rm = TRUE)
  weights_sd <- apply(boot_matrices, 2, sd, na.rm = TRUE)

  # Percentile CIs (vectorized via apply)
  ci_lower <- apply(boot_matrices, 2, quantile,
                    probs = ci_level / 2, na.rm = TRUE)
  ci_upper <- apply(boot_matrices, 2, quantile,
                    probs = 1 - ci_level / 2, na.rm = TRUE)

  # Basic (Davison & Hinkley 1997, eq. 5.6) intervals: reflect the
  # percentile bounds around the observed estimate. Corrects first-order
  # bootstrap bias; can produce bounds outside the weight range when the
  # observed value is near a boundary (documented behavior).
  if (identical(ci_method, "basic")) {
    basic_lower <- 2 * orig_flat - ci_upper
    ci_upper <- 2 * orig_flat - ci_lower
    ci_lower <- basic_lower
  }

  # p-values via vectorized sweep
  valid_rows <- rowSums(is.na(boot_matrices)) == 0
  bm <- boot_matrices[valid_rows, , drop = FALSE]
  n_valid <- nrow(bm)

  if (n_valid < 1) {
    p_values <- rep(1, length(orig_flat))
  } else if (inference == "stability") {
    cr_low <- pmin(orig_flat * consistency_range[1],
                   orig_flat * consistency_range[2])
    cr_high <- pmax(orig_flat * consistency_range[1],
                    orig_flat * consistency_range[2])
    # boundary = "inclusive" (default, tna-compatible): <=, >= count an
    # iteration that lands EXACTLY on the band edge as outside the band.
    # boundary = "strict": <, > count only iterations strictly outside.
    cmp_lo <- if (identical(boundary, "strict")) "<"  else "<="
    cmp_hi <- if (identical(boundary, "strict")) ">"  else ">="
    below <- sweep(bm, 2, cr_low,  cmp_lo)
    above <- sweep(bm, 2, cr_high, cmp_hi)
    p_counts <- colSums(below | above)
    p_values <- (p_counts + 1) / (n_valid + 1)
  } else {
    below_thresh <- sweep(abs(bm), 2, edge_threshold, "<")
    p_counts <- colSums(below_thresh)
    p_values <- (p_counts + 1) / (n_valid + 1)
  }

  # Reshape all to n_states x n_states matrices
  to_mat <- function(x) {
    matrix(x, n_states, n_states, dimnames = list(states, states))
  }

  p_mat <- to_mat(p_values)
  sig_mask <- (p_mat < ci_level) * 1
  sig_mat <- original_matrix * sig_mask

  list(
    mean        = to_mat(weights_mean),
    sd          = to_mat(weights_sd),
    ci_lower    = to_mat(ci_lower),
    ci_upper    = to_mat(ci_upper),
    p_values    = p_mat,
    significant = sig_mat,
    cr_lower    = to_mat(pmin(orig_flat * consistency_range[1],
                              orig_flat * consistency_range[2])),
    cr_upper    = to_mat(pmax(orig_flat * consistency_range[1],
                              orig_flat * consistency_range[2]))
  )
}


#' Build long-format summary data frame from bootstrap stats
#' @noRd
.build_bootstrap_summary <- function(stats, original_matrix, states, directed,
                                     ci_level, inference) {
  n <- length(states)
  dt <- data.table::data.table(
    from     = rep(states, each = n),
    to       = rep(states, times = n),
    weight   = as.vector(t(original_matrix)),
    mean     = as.vector(t(stats$mean)),
    sd       = as.vector(t(stats$sd)),
    p_value  = as.vector(t(stats$p_values)),
    sig      = as.vector(t(stats$p_values)) < ci_level,
    ci_lower = as.vector(t(stats$ci_lower)),
    ci_upper = as.vector(t(stats$ci_upper))
  )
  if (inference == "stability") {
    dt[, cr_lower := as.vector(t(stats$cr_lower))]
    dt[, cr_upper := as.vector(t(stats$cr_upper))]
  }

  # Filter: keep only non-zero original edges (including self-loops)
  if (directed) {
    dt <- dt[weight != 0]
  } else {
    dt <- dt[weight != 0 & from <= to]
  }

  as.data.frame(dt)
}


# ---- Print helpers (shared with print.net_certainty) ----

# sprintf("%-20s") pads by bytes, so a label holding the multibyte arrow came
# out short; these pad by display width, and widen the column to the longest
# label so a long state name does not push the numbers out of line.
.edge_label_width <- function(labels, min_width = 20L) {
  max(min_width, nchar(labels, type = "width"))
}

.pad_label <- function(labels, width) {
  paste0(labels, strrep(" ", pmax(0L, width - nchar(labels, type = "width"))))
}

# Header of the top-edges table: the interval is labelled with the level the
# intervals were computed at, not a fixed 95%.
.cat_edge_table_header <- function(width, ci_pct) {
  cat(sprintf("  %s   Mean     %-6s          p\n",
              .pad_label("Edge", width), paste(ci_pct, "CI")))
  cat(sprintf("  %s\n", strrep("-", width + 27L)))
}


# ---- S3 Methods ----

#' @rdname bootstrap_network
#' @param ... For the `print()` and `summary()` methods: ignored.
#' @param object For the `summary()` methods: the result of
#'   `bootstrap_network()`.
#' @return The `print()` methods return their input invisibly. The
#'   `summary()` method of a `net_bootstrap` returns its per-edge table, the
#'   `summary` field described above. For a `net_bootstrap_group`
#'   it returns the tables of all groups stacked into one data frame, with a
#'   first column `group` naming the network each row comes from; the group
#'   object keeps that table from the start (its `"summary"` attribute); for
#'   a `wtna_boot_mixed` it
#'   returns a list of the `transition` and `cooccurrence` tables.
#' @export
print.net_bootstrap <- function(x, ...) {
  method_labels <- c(
    relative      = "Transition Network (relative)",
    frequency     = "Transition Network (frequency)",
    co_occurrence = "Co-occurrence Network",
    glasso        = "Partial Correlation (EBICglasso)",
    pcor          = "Partial Correlation",
    cor           = "Correlation Network",
    attention     = "Attention Network",
    wtna          = "Window TNA"
  )
  lbl_raw   <- method_labels[x$method]
  label     <- if (is.na(lbl_raw)) sprintf("Network (%s)", x$method) else unname(lbl_raw)
  dir_label <- if (x$original$directed) "directed" else "undirected"
  ci_pct    <- sprintf("%d%%", round((1 - x$ci_level) * 100))
  n_orig    <- x$original$n_edges
  n_sig     <- x$model$n_edges
  n_nonsig  <- n_orig - n_sig

  # Top significant edges first
  s <- x$summary
  if (!is.null(s) && nrow(s) > 0L) {
    sig_s <- s[s$sig, , drop = FALSE]
    if (nrow(sig_s) > 0L) {
      sig_s <- sig_s[order(abs(sig_s$mean), decreasing = TRUE), , drop = FALSE]
      top   <- head(sig_s, 5L)
      lbl <- if (!is.null(top$from) && !is.null(top$to)) {
        paste0(top$from, " \u2192 ", top$to)
      } else {
        as.character(top$edge)
      }
      width <- .edge_label_width(lbl)
      .cat_edge_table_header(width, ci_pct)
      stars <- ifelse(top$p_value < 0.001, "***",
                      ifelse(top$p_value < 0.01, "** ", "*  "))
      cat(sprintf("  %s  %6.3f  [%5.3f, %5.3f]  %s\n",
                  .pad_label(lbl, width), top$mean, top$ci_lower,
                  top$ci_upper, stars),
          sep = "")
      if (nrow(sig_s) > 5L)
        cat(sprintf("  ... and %d more significant edges\n", nrow(sig_s) - 5L))
      cat("\n")
    }
  }

  cat(sprintf("Bootstrap Network  [%s | %s]\n", label, dir_label))
  cat(sprintf("  Iterations : %d  |  Nodes : %d\n", x$iter, x$original$n_nodes))
  cat(sprintf("  Edges      : %d significant / %d total\n", n_sig, n_orig))
  cat(sprintf("  CI         : %s  |  Inference: %s", ci_pct, x$inference))
  if (x$inference == "stability") {
    cat(sprintf("  |  CR [%.2f, %.2f]", x$consistency_range[1L], x$consistency_range[2L]))
  } else {
    cat(sprintf("  |  Threshold: %g", x$edge_threshold))
  }
  cat("\n")
  if (!is.null(x$clustering)) {
    cl <- x$clustering
    if (length(x$actor) == 1L) {
      cat(sprintf("  Actor      : %s (%d actors)\n", x$actor, cl$n_actors))
    } else {
      cat(sprintf(paste0("  Actor      : %s (nested); resampled: %s ",
                         "(%d actors, outermost level)\n"),
                  paste(x$actor, collapse = " > "), x$actor[[1L]],
                  cl$n_actors))
    }
    cat(sprintf("  Nesting    : %s  |  Design effect (edges) = %.2f\n",
                .format_icc(cl$icc, cl$icc_ci_lower, cl$icc_ci_upper,
                            cl$n_sequences, cl$n_actors),
                cl$deff_edges))
  }

  invisible(x)
}


#' @rdname bootstrap_network
#' @export
summary.net_bootstrap <- function(object, ...) {
  object$summary
}

#' @rdname bootstrap_network
#' @export
print.net_bootstrap_group <- function(x, ...) {
  grp_names <- names(x)
  iter      <- x[[1L]]$iter
  ci_pct    <- sprintf("%d%%", round((1 - x[[1L]]$ci_level) * 100))

  # Per-group edge counts
  grp_stats <- vapply(grp_names, function(nm) {
    b <- x[[nm]]
    c(total = b$original$n_edges, sig = b$model$n_edges)
  }, numeric(2L))

  # Helper: make edge keys "from->to" from a summary df
  .edge_keys <- function(s) {
    if (is.null(s) || nrow(s) == 0L) return(character(0L))
    idx <- which(s$sig)
    if (length(idx) == 0L) return(character(0L))
    paste0(s$from[idx], "\u2192", s$to[idx])
  }

  # Shared significant edges (present in all groups)
  sig_keys  <- lapply(grp_names, function(nm) .edge_keys(x[[nm]]$summary))
  shared    <- Reduce(intersect, sig_keys)

  # Detect disjoint state spaces (e.g. mcml: macro carries cluster
  # names, each within-cluster carries different code names). When the
  # union of every group's FULL edge set has empty pairwise overlap,
  # the "shared edges" metric is meaningless -- suppress its line.
  full_keys <- lapply(grp_names, function(nm) {
    s <- x[[nm]]$summary
    if (is.null(s) || nrow(s) == 0L) character(0L)
    else paste0(s$from, "->", s$to)
  })
  has_shared_space <- length(Reduce(intersect, full_keys)) > 0L

  # Top shared edges table first
  if (length(shared) > 0L) {
    .mean_for_key <- function(nm, key) {
      s <- x[[nm]]$summary
      parts <- strsplit(key, "\u2192", fixed = TRUE)[[1L]]
      idx <- which(s$from == parts[1L] & s$to == parts[2L])
      if (length(idx) == 0L) NA_real_ else s$mean[idx[1L]]
    }
    shared_max <- vapply(shared, function(k) {
      max(abs(vapply(grp_names, .mean_for_key, numeric(1L), key = k)), na.rm = TRUE)
    }, numeric(1L))
    top_shared <- shared[order(shared_max, decreasing = TRUE)][seq_len(min(5L, length(shared)))]
    hdr <- paste(sprintf("%-8s", grp_names), collapse = "  ")
    width <- .edge_label_width(top_shared, min_width = 22L)
    cat(sprintf("  %s  %s\n", .pad_label("Edge", width), hdr))
    cat(sprintf("  %s\n", strrep("-", width + 2L + 10L * length(grp_names))))
    vals_str <- vapply(top_shared, function(k) {
      grp_vals <- vapply(grp_names, .mean_for_key, numeric(1L), key = k)
      paste(sprintf("%-8.3f", grp_vals), collapse = "  ")
    }, character(1L))
    cat(sprintf("  %s  %s\n", .pad_label(top_shared, width), vals_str),
        sep = "")
    if (length(shared) > 5L)
      cat(sprintf("  ... and %d more shared significant edges\n", length(shared) - 5L))
    cat("\n")
  }

  # certainty() on a group returns this class too; it draws nothing, so it
  # has a prior where a bootstrap has a number of iterations.
  is_certainty <- inherits(x, "net_certainty_group")
  if (is_certainty) {
    cat(sprintf("Grouped Certainty  [%d groups | Dirichlet(%.2f) prior | %s CI]\n",
                length(grp_names), x[[1L]]$prior, ci_pct))
  } else {
    cat(sprintf("Grouped Bootstrap  [%d groups | %d iterations | %s CI]\n",
                length(grp_names), iter, ci_pct))
  }
  sig_word <- if (is_certainty) "certain" else "sig"
  group_lines <- vapply(grp_names, function(nm) {
    cl <- x[[nm]]$clustering
    paste0(
      sprintf("  %-20s  %d %s / %d total",
              nm, grp_stats["sig", nm], sig_word, grp_stats["total", nm]),
      if (!is.null(cl)) {
        sprintf("  |  %s: %d actors, ICC %.3f, design effect %.2f",
                x[[nm]]$actor[[1L]], cl$n_actors, cl$icc, cl$deff_edges)
      },
      "\n")
  }, character(1L))
  cat(group_lines, sep = "")
  if (has_shared_space) {
    cat(sprintf("  Shared (all groups)   %d edges\n", length(shared)))
  }

  invisible(x)
}

#' @rdname bootstrap_network
#' @export
summary.net_bootstrap_group <- function(object, ...) {
  attr(object, "summary") %||% .bootstrap_group_table(object)
}

# The edge tables of a grouped bootstrap stacked into one, with a `group`
# column first.
#' @noRd
.bootstrap_group_table <- function(object) {
  do.call(rbind, lapply(names(object), function(nm) {
    df       <- object[[nm]]$summary
    df$group <- nm
    df[c("group", setdiff(names(df), "group"))]
  }))
}


#' @rdname bootstrap_network
#' @export
print.wtna_boot_mixed <- function(x, ...) {
  cat("Mixed Window TNA Bootstrap\n")
  cat("-- Transition --\n")
  print(x$transition)
  cat("-- Co-occurrence --\n")
  print(x$cooccurrence)
  invisible(x)
}


#' @rdname bootstrap_network
#' @export
summary.wtna_boot_mixed <- function(object, ...) {
  list(
    transition   = summary(object$transition),
    cooccurrence = summary(object$cooccurrence)
  )
}


# ---- Latency networks ----

# Sequences an edge's replies must come from before its stability is judged:
# with all its replies in one sequence, every resample that holds the edge
# gives it the same mean, so its resampled weights cannot vary.
.LATENCY_MIN_SEQUENCES <- 5L

#' Resample a latency network
#'
#' Sequences are resampled with replacement and every edge's weight is
#' recomputed from the replies of the drawn sequences, as the estimator
#' computes it: each reply weighs `1 - (1 - floor) * min(t / cap, 1)`, and a
#' pair's replies are averaged (or summed). A percentile cap that the network
#' resolved from its own replies is resolved again in every resample; a cap
#' in seconds, or one pooled over the groups of a grouped build, is held at
#' its value.
#' @param x A latency `netobject`.
#' @return List with `weights` (iter x n^2, column-major, after the network's
#'   scaling and threshold), `scores` (iter x n^2 mean speed scores
#'   `(w - floor) / (1 - floor)`, NA where a resample holds no reply on the
#'   edge; NULL for `aggregate = "sum"`), `n_sequences` (n^2: the sequences
#'   with a reply on each edge) and `cap` (`"resolved"` or `"held"`).
#' @noRd
.bootstrap_latency <- function(x, iter) {
  params <- x$params
  floor <- params$latency_floor
  aggregate <- params$latency_aggregate %||% "mean"
  states <- rownames(x$weights)
  k <- length(states)
  times <- x$latency_time
  times <- as.matrix(times[startsWith(names(times), .LATENCY_TIME_PREFIX)])
  pairs <- .latency_reply_pairs(x$data, times, alphabet = states)$pairs
  changed <- pairs$from != pairs$to
  cell <- (pairs$to - 1L) * k + pairs$from
  percentile <- is.numeric(params$latency_cap)
  if (percentile && !identical(params$latency_scope %||% "global", "global")) {
    stop(errorCondition(paste0(
      "bootstrap_network() resamples latency networks with one cap for all ",
      "replies; per-unit caps (latency_scope = \"", params$latency_scope,
      "\") are not resampled. Build the network with latency_scope = ",
      "\"global\" or a cap in seconds."),
      class = c("nestimate_latency_resampling", "nestimate_bad_input"),
      call = NULL))
  }
  # a cap the network resolved from its own replies is resolved again; one
  # pooled over the groups of a grouped build belongs to all the groups
  own_cap <- percentile && any(changed) && isTRUE(all.equal(
    stats::quantile(pairs$latency[changed], params$latency_cap,
                    names = FALSE, type = 7),
    x$latency_cap_seconds))
  include_zeros <- .scaling_includes_zeros("latency")
  rows_by_sequence <- split(seq_len(nrow(pairs)), pairs$row)
  n_seq <- nrow(x$data)
  replicate_once <- function() {
    drawn <- sample.int(n_seq, n_seq, replace = TRUE)
    rows <- unlist(rows_by_sequence[as.character(drawn)], use.names = FALSE)
    cap <- if (own_cap && any(changed[rows])) {
      stats::quantile(pairs$latency[rows][changed[rows]], params$latency_cap,
                      names = FALSE, type = 7)
    } else x$latency_cap_seconds
    # a resample whose replies are mostly instant resolves a cap of 0, which
    # build_network() refuses too: no reply speed can be scored against it
    if (!(cap > 0)) {
      stop(errorCondition(sprintf(paste0(
        "bootstrap_network(): a resample's %s percentile of the reply ",
        "latencies is 0 seconds, so its replies cannot be weighted. Build the ",
        "network with the cap in seconds (e.g. latency_cap = \"%gs\") to ",
        "resample it."), format(params$latency_cap),
        signif(x$latency_cap_seconds, 3)),
        class = c("nestimate_latency_zero_cap", "nestimate_bad_input"),
        call = NULL))
    }
    w <- 1 - (1 - floor) * pmin(pairs$latency[rows] / cap, 1)
    total <- numeric(k * k)
    sums <- rowsum(w, cell[rows], reorder = FALSE)
    total[as.integer(rownames(sums))] <- sums[, 1L]
    count <- tabulate(cell[rows], k * k)
    mean_w <- ifelse(count > 0, total / pmax(count, 1L), 0)
    mat <- matrix(if (identical(aggregate, "mean")) mean_w else total, k, k)
    c(.finish_replicate(mat, FALSE, x$scaling, include_zeros, x$threshold),
      ifelse(count > 0, (mean_w - floor) / (1 - floor), NA_real_))
  }
  draws <- t(vapply(seq_len(iter), function(i) replicate_once(),
                    numeric(2L * k * k)))
  first <- seq_len(k * k)
  list(weights = draws[, first, drop = FALSE],
       scores = if (identical(aggregate, "mean")) draws[, k * k + first, drop = FALSE],
       n_sequences = tabulate(cell[!duplicated(paste(cell, pairs$row))], k * k),
       cap = if (own_cap) "resolved" else "held")
}

#' The stability test of a mean latency network, on the speed score
#'
#' A mean reply weight lies in `[floor, 1]`, so a band around the weight
#' passes almost every edge. The test is run on the speed score
#' `(w - floor) / (1 - floor)`, which runs from 0 (every reply at or past the
#' cap) to 1 (every reply instant): an edge passes when its resampled score
#' stays inside `consistency_range` times the observed score. A resample that
#' holds no reply on the edge counts as outside. An edge whose replies come
#' from fewer than `.LATENCY_MIN_SEQUENCES` sequences gets a p-value of 1.
#' @return The `stats` list with `p_values`, `significant`, `cr_lower` and
#'   `cr_upper` replaced, `cv` (the coefficient of variation of the
#'   resampled scores) and `n_sequences` added.
#' @noRd
.latency_stability <- function(stats, boot, original, floor, consistency_range,
                               boundary, ci_level) {
  k <- nrow(original)
  observed <- (as.vector(original) - floor) / (1 - floor)
  lo <- pmin(observed * consistency_range[1], observed * consistency_range[2])
  hi <- pmax(observed * consistency_range[1], observed * consistency_range[2])
  scores <- boot$scores
  low <- if (identical(boundary, "strict")) sweep(scores, 2, lo, "<") else
    sweep(scores, 2, lo, "<=")
  high <- if (identical(boundary, "strict")) sweep(scores, 2, hi, ">") else
    sweep(scores, 2, hi, ">=")
  outside <- colSums(low | high | is.na(scores))
  p <- (outside + 1) / (nrow(scores) + 1)
  p[boot$n_sequences < .LATENCY_MIN_SEQUENCES] <- 1
  m <- colMeans(scores, na.rm = TRUE)
  s <- apply(scores, 2, stats::sd, na.rm = TRUE)
  cv <- ifelse(m > 0, s / m, Inf)
  cv[boot$n_sequences == 0L] <- NA_real_
  as_mat <- function(v) matrix(v, k, k, dimnames = dimnames(original))
  stats$p_values <- as_mat(p)
  stats$significant <- original * (stats$p_values < ci_level)
  stats$cr_lower <- as_mat(floor + (1 - floor) * lo)
  stats$cr_upper <- as_mat(floor + (1 - floor) * hi)
  stats$cv <- as_mat(cv)
  stats$n_sequences <- as_mat(boot$n_sequences)
  stats
}


#' The note that the latency bootstrap is experimental
#' @noRd
.latency_experimental_message <- function() {
  cond <- simpleMessage(paste0(
    "The latency bootstrap is experimental: tested in simulation, it ",
    "showed promising results, but it has not yet been tested widely on ",
    "real-life datasets.\n"))
  class(cond) <- c("nestimate_experimental", class(cond))
  message(cond)
}
