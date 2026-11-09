# ---- Multi-Order Generative Model (MOGen) ----
#
# Implements multi-order De Bruijn graph models for sequential data.
# Based on Scholtes (2017) "When is a Network a Network?" and
# Gote & Scholtes (2023) "Predicting variable-length paths using MOGen".
#
# Key idea: Build higher-order De Bruijn graphs at orders k=0,1,...,K,
# compute transition matrices, and use AIC/BIC/LRT to select the optimal
# Markov order for the data.

# ---------------------------------------------------------------------------
# Internal: Count k-grams and transitions
# ---------------------------------------------------------------------------

#' Count k-grams and transitions between consecutive k-grams
#'
#' At order k, nodes are k-tuples of states from trajectories.
#' Edges connect consecutive k-tuples (overlapping by k-1 states).
#'
#' @param trajectories List of character vectors.
#' @param k Integer order (k >= 1).
#' @return List with nodes, node_counts, edges (data.frame from/to/weight).
#' @noRd
.mogen_count_kgrams <- function(trajectories, k) {
  stopifnot(k >= 1L)

  results <- lapply(trajectories, function(traj) {
    n <- length(traj)
    if (n < k) return(list(nodes = character(0L), from = character(0L),
                            to = character(0L)))

    # k-grams via sliding window
    kgrams <- vapply(seq_len(n - k + 1L), function(i) {
      paste(traj[i:(i + k - 1L)], collapse = .HON_SEP)
    }, character(1L))

    # Transitions between consecutive k-grams
    if (length(kgrams) > 1L) {
      list(nodes = kgrams, from = kgrams[-length(kgrams)],
           to = kgrams[-1L])
    } else {
      list(nodes = kgrams, from = character(0L), to = character(0L))
    }
  })

  all_nodes <- unlist(lapply(results, `[[`, "nodes"))
  all_from <- unlist(lapply(results, `[[`, "from"))
  all_to <- unlist(lapply(results, `[[`, "to"))

  node_tab <- table(all_nodes)
  nodes <- names(node_tab)
  node_counts <- as.integer(node_tab)
  names(node_counts) <- nodes

  if (length(all_from) == 0L) {
    edges <- data.frame(from = character(0L), to = character(0L),
                        weight = integer(0L), stringsAsFactors = FALSE)
  } else {
    edge_keys <- paste(all_from, all_to, sep = "\x02")
    edge_tab <- table(edge_keys)
    edge_split <- strsplit(names(edge_tab), "\x02", fixed = TRUE)
    edges <- data.frame(
      from = vapply(edge_split, `[`, character(1L), 1L),
      to = vapply(edge_split, `[`, character(1L), 2L),
      weight = as.integer(edge_tab),
      stringsAsFactors = FALSE
    )
  }

  list(nodes = nodes, node_counts = node_counts, edges = edges)
}

# ---------------------------------------------------------------------------
# Internal: Build transition matrix
# ---------------------------------------------------------------------------

#' Build row-stochastic transition matrix from k-gram edge counts
#'
#' @param nodes Character vector of node names.
#' @param edges Data frame with from, to, weight columns.
#' @return Named square matrix (row-stochastic).
#' @noRd
.mogen_transition_matrix <- function(nodes, edges) {
  n <- length(nodes)
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(nodes, nodes))

  if (nrow(edges) > 0L) {
    idx <- cbind(match(edges$from, nodes), match(edges$to, nodes))
    mat[idx] <- edges$weight
  }

  row_sums <- rowSums(mat)
  nonzero <- row_sums > 0
  mat[nonzero, ] <- mat[nonzero, ] / row_sums[nonzero]
  mat
}

# ---------------------------------------------------------------------------
# Internal: Order-0 marginal distribution
# ---------------------------------------------------------------------------

#' Compute order-0 marginal distribution over states
#'
#' @param trajectories List of character vectors.
#' @return Named numeric vector of state probabilities.
#' @noRd
.mogen_marginal <- function(trajectories) {
  all_states <- unlist(trajectories)
  tab <- table(all_states)
  probs <- as.numeric(tab) / sum(tab)
  names(probs) <- names(tab)
  probs
}

# ---------------------------------------------------------------------------
# Internal: Log-likelihood computation
# ---------------------------------------------------------------------------

#' Compute log-likelihood of trajectories under k-th order model
#'
#' Uses hierarchical decomposition:
#' - Step 1: order-0 probability (state marginal)
#' - Steps 2..k: use increasing orders (order j for step j)
#' - Steps k+1..l: use order-k transitions
#'
#' @param trajectories List of character vectors.
#' @param k Integer model order.
#' @param trans_mats List of transition matrices (index j+1 = order j).
#'   Index 1 is order-0 marginal (named numeric vector).
#' @return Numeric scalar (log-likelihood).
#' @noRd
.mogen_log_likelihood <- function(trajectories, k, trans_mats) {
  log_eps <- log(.Machine$double.eps)
  p0 <- trans_mats[[1L]]

  sum(vapply(trajectories, function(traj) {
    n <- length(traj)
    if (n == 0L) return(0)

    # Order 0: initial state probability
    ll <- if (traj[1L] %in% names(p0) && p0[traj[1L]] > 0)
      log(p0[traj[1L]]) else log_eps

    if (n < 2L) return(ll)

    # Steps 2..n: use order min(step-1, k)
    step_lls <- vapply(2L:n, function(step) {
      order_used <- min(step - 1L, k)

      if (order_used == 0L) {
        # Order 0: each state drawn from marginal (no transition context)
        if (traj[step] %in% names(p0) && p0[traj[step]] > 0) # nocov start
          log(p0[traj[step]]) else log_eps # nocov end
      } else {
        src_start <- step - order_used
        src_key <- paste(traj[src_start:(step - 1L)], collapse = .HON_SEP)
        tgt_key <- paste(traj[(src_start + 1L):step], collapse = .HON_SEP)

        tm <- trans_mats[[order_used + 1L]]
        if (src_key %in% rownames(tm) && tgt_key %in% colnames(tm)) {
          p <- tm[src_key, tgt_key]
          if (p > 0) log(p) else log_eps
        } else {
          log_eps # nocov
        }
      }
    }, numeric(1L))

    ll + sum(step_lls)
  }, numeric(1L)))
}

# ---------------------------------------------------------------------------
# Internal: Degrees of freedom
# ---------------------------------------------------------------------------

#' Compute degrees of freedom for a row-stochastic transition matrix
#'
#' For each row with m non-zero entries: m - 1 free parameters.
#'
#' @param trans_mat Square transition matrix.
#' @return Integer (total DOF).
#' @noRd
.mogen_layer_dof <- function(trans_mat) {
  row_nonzero <- rowSums(trans_mat > 0)
  as.integer(sum(pmax(row_nonzero - 1L, 0L)))
}

#' Order selected by sequential likelihood-ratio tests
#'
#' The largest order whose test against the order below rejects at
#' \code{alpha} (Scholtes 2017), or 0 when none does.
#'
#' @param lrt_p p-values indexed by order 0..K (order 0 is \code{NA}).
#' @param alpha Significance threshold.
#' @return Integer order.
#' @noRd
.mogen_lrt_order <- function(lrt_p, alpha) {
  rejected <- which(!is.na(lrt_p) & lrt_p < alpha)
  if (length(rejected)) as.integer(max(rejected) - 1L) else 0L
}

#' Path-based degrees of freedom of each layer (Scholtes 2017, Eq. 7)
#'
#' The parameters of an order-k layer are the transitions of every path of
#' length k that the first-order topology allows, observed or not:
#' \eqn{\mathbf{1}^\top A^k \mathbf{1}}, less one per k-order node with a
#' possible successor (rows sum to one):
#' \eqn{\mathbf{1}^\top A^{k-1} d}, with \eqn{d_j = 1} when state j has an
#' outgoing transition. \eqn{A} is the binary matrix of observed first-order
#' transitions. Matches pathpy 2 \code{HigherOrderNetwork.dof_paths}.
#'
#' @param first_tm Order-1 transition matrix (states x states).
#' @param max_order Highest layer.
#' @return Numeric vector of length \code{max_order} (whole numbers, kept
#'   double because \eqn{s^k} overflows integers for large state spaces).
#' @noRd
.mogen_path_dofs <- function(first_tm, max_order) {
  A <- (first_tm > 0) * 1
  # sum(A^k) - colSums(A^(k-1)) . d  =  sum_i r_(k-1)[i] * (outdeg_i - d_i),
  # with r_k = 1' A^k (r_0 = 1'). Each state adds max(outdeg - 1, 0) per path
  # ending there: non-negative terms, so no large difference cancels and no
  # matrix power is stored (O(states) per layer).
  extra <- pmax(rowSums(A) - 1, 0)
  r0 <- rep(1, nrow(A))
  # r_0..r_(K-1) as a list. Reduce() returns the bare init on empty input and
  # unlists when every element has length one (a one-state model), so build
  # the K = 1 case directly and keep simplify = FALSE.
  r <- if (max_order > 1L) {
    Reduce(function(v, i) as.numeric(v %*% A), seq_len(max_order - 1L),
           accumulate = TRUE, init = r0, simplify = FALSE)
  } else {
    list(r0)
  }
  dofs <- vapply(r, function(v) sum(v * extra), numeric(1L))
  bad <- which(!is.finite(dofs))
  if (length(bad)) {
    stop(errorCondition(
      sprintf(paste0("path degrees of freedom exceed double precision at ",
                     "order %d; lower `max_order` or use `dof = \"observed\"`"),
              bad[1L]),
      class = "nestimate_dof_overflow", call = NULL
    ))
  }
  dofs
}

# ---------------------------------------------------------------------------
# Main function: build_mogen
# ---------------------------------------------------------------------------

#' Build Multi-Order Generative Model (MOGen)
#'
#' @description
#' `build_mogen()` fits Markov chains of every order from 0 to `max_order`
#' to a set of sequences and selects the order that best balances fit
#' against the number of parameters, by AIC, BIC or likelihood-ratio tests.
#' At order k the chain is a De Bruijn graph: its nodes are the k-tuples of
#' consecutive states seen in the data, and an edge joins two tuples that
#' overlap in k - 1 states, so each edge is an observed path of k + 1
#' states. Order 0 is the marginal distribution of the states. It is the
#' multi-order model of Scholtes (2017), and takes wide sequence data, a
#' list of state vectors, a `tna` model or a network built from sequences.
#'
#' The selected order is the number of earlier states the next state
#' depends on, as far as the sample can show. The number of parameters
#' grows quickly with the order (`summary()` lists them per order), so few
#' or short sequences favour low orders, and AIC and BIC can disagree:
#' `print()` shows both. Before reading a higher-order transition, check how
#' often its context was observed with [mogen_transitions()];
#' [markov_order_test()] tests the order against a permutation null, and
#' [path_dependence()] shows which contexts carry the extra information.
#'
#' @details
#' The log-likelihood of the order-k model scores the first state of each
#' trajectory by the marginal distribution, the j-th state by the order
#' j - 1 transitions while j <= k, and every later state by the order-k
#' transitions. The likelihood-ratio test of order k against order k - 1
#' refers twice the gain in log-likelihood to a chi-square distribution with
#' the order-k layer's degrees of freedom; its p-value is `NA` for order 0
#' and where the statistic or the degrees of freedom are not positive. With
#' `criterion = "lrt"` the selected order is the largest k whose test
#' rejects at `lrt_alpha` (Scholtes 2017), or 0 when none does. Every
#' p-value is kept in the result, so a non-monotone run of tests is visible.
#'
#' The model is the one pathpy 2 calls `MultiOrderModel`. It is not the
#' multi-order generative model of Gote et al. (2023), which adds start and
#' end states and selects by AIC only; the function name is kept for
#' compatibility.
#'
#' @param data The sequences: a data frame with one row per trajectory and
#'   one column per time point, a character matrix laid out the same way, a
#'   list of state vectors (one per trajectory), a `tna` model (its integer
#'   states are decoded to their labels) or a `netobject` built from
#'   sequences (its states are read as stored). A data frame needs position
#'   columns named as a numbered run (`T1`, `T2`, ...); any other frame
#'   raises an error of class `nestimate_bad_input` that says to declare the
#'   roles in [build_network()] and pass the network. A missing value ends a
#'   run of observed states, so each run is its own trajectory, and
#'   trajectories with fewer than two states are dropped; with none left,
#'   `build_mogen()` stops. A one-hot network, an association network
#'   (`cor`, `pcor`, `glasso`, `ising`, `mgm`), a `netobject_group` and a
#'   0/1 matrix raise an error of class `nestimate_bad_input`.
#' @param max_order The highest order fitted, a whole number of at least 1;
#'   a value such as `2.7` is an error, not truncated. Default `5`. A value
#'   at or above the length of the longest trajectory is lowered to that
#'   length minus 1, with a message.
#' @param criterion How the order is selected: `"aic"` (default) or
#'   `"bic"` takes the order with the smallest criterion, and `"lrt"` uses
#'   the sequential likelihood-ratio tests described in *Details*.
#' @param lrt_alpha The significance level of the likelihood-ratio tests,
#'   strictly between 0 and 1. Default `0.01`. It only affects the selection
#'   with `criterion = "lrt"`; the p-values are computed under every
#'   criterion.
#' @param dof How the parameters of each layer are counted. `"paths"`
#'   (default) counts every transition of the paths of length k that the
#'   observed first-order topology allows, as Scholtes (2017, Eq. 7) and
#'   pathpy do. `"observed"` counts only the transitions seen in the data
#'   (the rule before Nestimate 0.9.28); that count can fall as the order
#'   rises, and it under-penalises high orders on sparse data.
#' @return A `net_mogen` object (class `c("net_mogen", "cograph_network")`).
#'   `print()` shows the selected order and the AIC and BIC of every order,
#'   with each minimum starred. `summary()` prints the states and the
#'   selection, then returns, visibly, one row per order as a data frame
#'   (`order`, `layer_dof`, `cum_dof`, `loglik`, `aic`, `bic`, `best`
#'   (`"AIC"`, `"BIC"` or `"AIC+BIC"` on the order each criterion picks,
#'   empty elsewhere) and `selected` (`"<--"` on the chosen order)).
#'   `plot()` draws the criteria or the log-likelihood against the order and
#'   returns its input invisibly; `print()` also returns its input
#'   invisibly. [mogen_transitions()] gives the transitions of one
#'   order, one row per path, and [pathways()] the same as character
#'   strings. The object's fields:
#'   \describe{
#'     \item{optimal_order}{The selected order.}
#'     \item{criterion}{The criterion used for the selection.}
#'     \item{orders}{The orders fitted, 0 to `max_order` after any capping.}
#'     \item{aic, bic, log_likelihood}{Named numeric vectors, one value per
#'       order.}
#'     \item{dof}{Named numeric vector of the cumulative degrees of freedom
#'       of each model.}
#'     \item{layer_dof}{Named numeric vector of the degrees of freedom of
#'       each layer.}
#'     \item{dof_method}{`"paths"` or `"observed"`.}
#'     \item{lrt_p}{Named numeric vector: the likelihood-ratio p-value of
#'       order k against order k - 1.}
#'     \item{transition_matrices}{List of row-stochastic transition
#'       matrices (a context never followed by another state keeps a zero
#'       row), the first (order 0) held as the named vector of marginal
#'       probabilities.}
#'     \item{count_matrices}{List of the matching counts, in the same
#'       order, the first (order 0) as a named vector of state counts;
#'       [mogen_transitions()] reads them.}
#'     \item{states}{The distinct first-order states.}
#'     \item{n_paths}{Number of trajectories fitted.}
#'     \item{n_observations}{Total number of state observations.}
#'     \item{weights}{The transition matrix of the selected order, as the
#'       `cograph_network` weight matrix (a 1 x n matrix named `"marginal"`
#'       when the selected order is 0). Its dimnames are internal k-tuple
#'       keys (states joined by a non-printing separator), not arrow
#'       notation.}
#'     \item{nodes}{Data frame (`id`, `label`, `name`) of the nodes of the
#'       selected order.}
#'     \item{edges}{The `cograph_network` edge data frame, with integer
#'       `from` and `to` node indices and a numeric `weight`;
#'       [mogen_transitions()] gives the readable table.}
#'     \item{directed}{Always `TRUE`.}
#'     \item{n_nodes, n_edges}{Counts for the graph of the selected order.}
#'     \item{meta}{The `cograph_network` metadata list.}
#'     \item{node_groups}{Always `NULL`.}
#'   }
#'
#' @references
#' Scholtes, I. (2017). When is a Network a Network? Multi-Order Graphical
#' Model Selection in Pathways and Temporal Networks. \emph{KDD 2017}.
#'
#' Gote, C., Casiraghi, G., Schweitzer, F., & Scholtes, I. (2023). Predicting
#' variable-length paths in networked systems using multi-order generative
#' models. \emph{Applied Network Science}, 8, 68.
#'
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' mg <- build_mogen(seqs, max_order = 2)
#'
#' \donttest{
#' trajs <- list(c("A","B","C","D"), c("A","B","D","C"),
#'               c("B","C","D","A"), c("C","D","A","B"))
#' m <- build_mogen(trajs, max_order = 3)
#' print(m)
#' plot(m)
#' }
#'
#' @export
build_mogen <- function(data, max_order = 5L, criterion = c("aic", "bic", "lrt"),
                        lrt_alpha = 0.01, dof = c("paths", "observed")) {
  .refuse_onehot_network(data, "build_mogen()")
  .refuse_binary_matrix(data, "build_mogen()")
  if (is.data.frame(data)) {
    .refuse_undeclared_roles(data, "build_mogen()", via_network = TRUE,
                             whole_frame = TRUE)
  }
  data <- .coerce_sequence_input(data)
  criterion <- .ho_match_arg(criterion, "criterion")
  dof <- .ho_match_arg(dof, "dof")
  .assert(
    "'data' must be a data.frame or list" =
      is.data.frame(data) || is.list(data),
    "'max_order' must be a whole number" =
      is.numeric(max_order) && length(max_order) == 1L &&
      max_order == round(max_order),
    "'max_order' must be >= 1" = max_order >= 1L,
    "'lrt_alpha' must be in (0, 1)" = lrt_alpha > 0 && lrt_alpha < 1
  )
  max_order <- as.integer(max_order)

  trajectories <- .hon_parse_input(data, collapse_repeats = FALSE)
  if (length(trajectories) == 0L) {
    stop(errorCondition(
      "No valid trajectories (each must have at least 2 states)",
      class = "nestimate_bad_input", call = NULL))
  }

  # Cap max_order at max path length - 1
  max_path_len <- max(vapply(trajectories, length, integer(1L)))
  if (max_order >= max_path_len) {
    max_order <- max_path_len - 1L
    message(sprintf("max_order capped at %d (longest path = %d)", max_order,
                    max_path_len))
  }

  n_obs <- sum(vapply(trajectories, length, integer(1L)))

  # --- Order 0: marginal distribution ---
  marginal <- .mogen_marginal(trajectories)
  trans_mats <- vector("list", max_order + 1L)
  count_mats <- vector("list", max_order + 1L)
  trans_mats[[1L]] <- marginal
  count_mats[[1L]] <- marginal * sum(vapply(trajectories, length, integer(1L)))

  layer_dofs <- numeric(max_order + 1L)
  layer_dofs[1L] <- length(marginal) - 1L

  # --- Orders 1..max_order: De Bruijn graphs ---
  lapply(seq_len(max_order), function(k) {
    kg <- .mogen_count_kgrams(trajectories, k)
    tm <- .mogen_transition_matrix(kg$nodes, kg$edges)
    # Also store raw count matrix
    cm <- matrix(0L, nrow = length(kg$nodes), ncol = length(kg$nodes),
                 dimnames = list(kg$nodes, kg$nodes))
    if (nrow(kg$edges) > 0L) {
      idx <- cbind(match(kg$edges$from, kg$nodes),
                   match(kg$edges$to, kg$nodes))
      cm[idx] <- kg$edges$weight
    }
    trans_mats[[k + 1L]] <<- tm
    count_mats[[k + 1L]] <<- cm
    layer_dofs[k + 1L] <<- .mogen_layer_dof(tm)
    NULL
  })

  if (dof == "paths") {
    layer_dofs[-1L] <- .mogen_path_dofs(trans_mats[[2L]], max_order)
  }

  # --- Cumulative DOF and log-likelihoods ---
  cum_dofs <- cumsum(layer_dofs)

  logliks <- vapply(0L:max_order, function(k) {
    .mogen_log_likelihood(trajectories, k, trans_mats[seq_len(k + 1L)])
  }, numeric(1L))

  # Likelihood-ratio test of each order against the one below
  lrt_p <- c(NA_real_, vapply(seq_len(max_order), function(k) {
    x <- -2 * (logliks[k] - logliks[k + 1L])
    df_diff <- layer_dofs[k + 1L]
    if (df_diff > 0 && x > 0) {
      stats::pchisq(x, df = df_diff, lower.tail = FALSE)
    } else {
      NA_real_
    }
  }, numeric(1L)))

  # --- Information criteria ---
  aics <- 2 * cum_dofs - 2 * logliks
  bics <- log(n_obs) * cum_dofs - 2 * logliks

  order_names <- paste0("order_", 0L:max_order)
  names(aics) <- order_names
  names(bics) <- order_names
  names(logliks) <- order_names
  names(cum_dofs) <- order_names
  names(layer_dofs) <- order_names
  names(lrt_p) <- order_names

  # --- Select optimal order ---
  if (criterion == "aic") {
    optimal_order <- which.min(aics) - 1L
  } else if (criterion == "bic") {
    optimal_order <- which.min(bics) - 1L
  } else {
    optimal_order <- .mogen_lrt_order(lrt_p, lrt_alpha)
  }

  # Extract optimal-order transition matrix for cograph compatibility
  opt_mat <- trans_mats[[optimal_order + 1L]]
  if (is.null(dim(opt_mat))) {
    # Order 0: marginal vector -> 1xn matrix
    opt_mat <- matrix(opt_mat, nrow = 1,
                      dimnames = list("marginal", names(opt_mat)))
  }
  opt_nodes <- rownames(opt_mat) %||% names(marginal)
  cg <- .ho_cograph_fields(opt_mat, opt_nodes, method = "mogen")

  result <- list(
    optimal_order = as.integer(optimal_order),
    criterion = criterion,
    orders = 0L:max_order,
    aic = aics,
    bic = bics,
    log_likelihood = logliks,
    dof = cum_dofs,
    layer_dof = layer_dofs,
    dof_method = dof,
    lrt_p = lrt_p,
    transition_matrices = trans_mats,
    count_matrices = count_mats,
    states = names(marginal),
    n_paths = length(trajectories),
    n_observations = n_obs,
    weights = cg$weights,
    nodes = cg$nodes,
    edges = cg$edges,
    directed = TRUE,
    n_nodes = cg$n_nodes,
    n_edges = cg$n_edges,
    meta = cg$meta,
    node_groups = NULL
  )

  class(result) <- c("net_mogen", "cograph_network")
  result
}

# ---------------------------------------------------------------------------
# Utility: Extract transitions at any order
# ---------------------------------------------------------------------------

#' Extract Transition Table from a MOGen Model
#'
#' @description
#' `mogen_transitions()` returns the transitions of a fitted [build_mogen()]
#' model at one order, one row per path, most frequent first, with the
#' path written out in arrow notation, its observed count and its
#' transition probability. At order k each edge of the De Bruijn graph is a
#' path of k + 1 states: at order 2, the edge from node `"AI -> FAIL"` to
#' node `"FAIL -> SOLVE"` is the path `AI -> FAIL -> SOLVE`, read as the
#' probability of `SOLVE` after the context `AI -> FAIL`.
#'
#' A probability is the share of its context's observations that moved to
#' that state, so in a context seen a few times it rests on those few
#' observations: read it with its `count`, and use `min_count` to leave
#' rare transitions out.
#'
#' @param x A `net_mogen` object from [build_mogen()].
#' @param order The order whose transitions are returned, a whole number
#'   from 1 to the highest order fitted; a non-integer value is an error,
#'   not truncated. Default `NULL`, the order the model selected. Order 0
#'   has no transitions, so when the model selected order 0, give `order`
#'   explicitly; otherwise `mogen_transitions()` stops.
#' @param min_count The smallest observed count a transition needs to be
#'   listed. Default `1`, every observed transition. A value below 1 also
#'   lists the transitions that never occurred, with count and probability
#'   0: every pair of observed contexts where the second continues the first
#'   (at order 1, every pair of states).
#' @return A data frame with one row per listed transition, sorted by
#'   `count` (descending), with columns:
#'   \describe{
#'     \item{path}{The full sequence of states, such as
#'       `"AI -> FAIL -> SOLVE"`.}
#'     \item{count}{How often the transition was observed.}
#'     \item{probability}{The transition probability P(to | from), rounded
#'       to 4 decimal places.}
#'     \item{from}{The context: the k states the transition starts from.}
#'     \item{to}{The next state.}
#'   }
#'   When no transition reaches `min_count`, a data frame with the same
#'   columns and no rows.
#'
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' mg <- build_mogen(seqs, max_order = 2)
#' mogen_transitions(mg, order = 1)
#'
#' \donttest{
#' trajs <- list(c("A","B","C","D"), c("A","B","D","C"),
#'               c("B","C","D","A"), c("C","D","A","B"))
#' m <- build_mogen(trajs, max_order = 3)
#' mogen_transitions(m, order = 1)
#' }
#'
#' @export
mogen_transitions <- function(x, order = NULL, min_count = 1L) {
  .assert("'x' must be a net_mogen object from build_mogen()" =
            inherits(x, "net_mogen"))
  if (is.null(order)) order <- x$optimal_order
  .assert(
    "'order' must be a whole number" =
      is.numeric(order) && length(order) == 1L &&
      order == round(order)
  )
  order <- as.integer(order)
  min_count <- as.integer(min_count)
  .assert(
    "'order' must be >= 1" = order >= 1L,
    "'order' exceeds max tested order" = order <= max(x$orders)
  )

  tm <- x$transition_matrices[[order + 1L]]
  cm <- x$count_matrices[[order + 1L]]

  # Find edges with sufficient count
  idx <- which(cm >= min_count, arr.ind = TRUE)

  if (nrow(idx) == 0L) {
    return(data.frame(path = character(0L), count = integer(0L),
                      probability = numeric(0L), from = character(0L),
                      to = character(0L), stringsAsFactors = FALSE))
  }

  from_raw <- rownames(tm)[idx[, 1L]]
  to_raw <- colnames(tm)[idx[, 2L]]

  # Reconstruct full path and extract context/next_state
  parsed <- lapply(seq_len(nrow(idx)), function(i) {
    from_parts <- strsplit(from_raw[i], .HON_SEP, fixed = TRUE)[[1L]]
    to_parts <- strsplit(to_raw[i], .HON_SEP, fixed = TRUE)[[1L]]
    next_state <- to_parts[length(to_parts)]
    list(
      path = paste(c(from_parts, next_state), collapse = " -> "),
      from = paste(from_parts, collapse = " -> "),
      to = next_state,
      # the next context must continue this one (its first k - 1 states are
      # this context's last k - 1); always true at order 1
      follows = identical(from_parts[-1L], to_parts[-length(to_parts)])
    )
  })
  follows <- vapply(parsed, `[[`, logical(1L), "follows")
  parsed <- parsed[follows]
  idx <- idx[follows, , drop = FALSE]

  result <- data.frame(
    path = vapply(parsed, `[[`, character(1L), "path"),
    count = as.integer(cm[idx]),
    probability = round(tm[idx], 4),
    from = vapply(parsed, `[[`, character(1L), "from"),
    to = vapply(parsed, `[[`, character(1L), "to"),
    stringsAsFactors = FALSE
  )
  result <- result[order(-result$count), ]
  rownames(result) <- NULL
  result
}

# ---------------------------------------------------------------------------
# Utility: Count path frequencies
# ---------------------------------------------------------------------------

#' Count Path Frequencies in Trajectory Data
#'
#' @description
#' `path_counts()` counts how often each path of `k` consecutive states
#' (a k-gram) occurs across all trajectories, and returns one row per path,
#' most frequent first. Every position of a trajectory starts a path, so the
#' paths of one trajectory overlap. The counts are descriptive: they show
#' which runs of states dominate the data before a model such as
#' [build_mogen()] or [build_hon()] is fitted to them.
#'
#' @param data The trajectories: a list of state vectors, a data frame or
#'   character matrix with one row per trajectory and one column per time
#'   point, a `tna` model or a `netobject` built from sequences. A data
#'   frame needs position columns named as a numbered run (`T1`, `T2`, ...);
#'   any other frame, a one-hot or association network and a 0/1 matrix
#'   raise an error of class `nestimate_bad_input`. A missing value ends a
#'   run of observed states, and no path spans it.
#' @param k The number of states in each path, a whole number of at least
#'   2; a non-integer value is an error, not truncated. Default `2`, which
#'   counts single transitions; `3` counts two-step paths, and so on.
#' @param top The number of most frequent paths to keep. Default `NULL`,
#'   every path.
#' @return A data frame with one row per distinct path, sorted by `count`
#'   (descending), with columns `path` (the k states in arrow notation, such
#'   as `"A -> B"`), `count` and `proportion` (the path's share of all paths
#'   counted, rounded to 4 decimal places). When no trajectory holds `k`
#'   consecutive observed states, the data frame has no rows.
#'
#' @examples
#' trajs <- list(c("A","B","C","D"), c("A","B","D","C"))
#' path_counts(trajs, k = 2)
#'
#' \donttest{
#' path_counts(trajs, k = 3, top = 10)
#' }
#'
#' @export
path_counts <- function(data, k = 2L, top = NULL) {
  .refuse_onehot_network(data, "path_counts()")
  .refuse_binary_matrix(data, "path_counts()")
  if (is.data.frame(data)) {
    .refuse_undeclared_roles(data, "path_counts()", via_network = TRUE,
                             whole_frame = TRUE)
  }
  data <- .coerce_sequence_input(data)
  .assert(
    "'k' must be a whole number" =
      is.numeric(k) && length(k) == 1L && k == round(k)
  )
  k <- as.integer(k)
  .assert("'k' must be >= 2" = k >= 2L)

  # Normalize input
  if (is.data.frame(data)) {
    trajectories <- .sequence_rows(data)
  } else {
    trajectories <- lapply(data, as.character)
  }

  # a missing value ends a run: k-grams never span it (dropping the NA and
  # joining the two sides would count paths that were never observed)
  runs <- unlist(lapply(trajectories, function(traj) {
    observed <- !is.na(traj)
    run_id <- cumsum(!observed)
    unname(split(traj[observed], run_id[observed]))
  }), recursive = FALSE)
  all_grams <- unlist(lapply(runs, function(traj) {
    n <- length(traj)
    if (n < k) return(character(0L))
    vapply(seq_len(n - k + 1L), function(i) {
      paste(traj[i:(i + k - 1L)], collapse = " -> ")
    }, character(1L))
  }))

  tbl <- sort(table(all_grams), decreasing = TRUE)
  result <- data.frame(
    path = names(tbl),
    count = as.integer(tbl),
    proportion = round(as.numeric(tbl) / sum(tbl), 4),
    stringsAsFactors = FALSE
  )
  rownames(result) <- NULL

  if (!is.null(top)) result <- head(result, as.integer(top))
  result
}

# ---------------------------------------------------------------------------
# Utility: State frequencies
# ---------------------------------------------------------------------------

#' Compute State Frequencies from Trajectory Data
#'
#' @description
#' `state_frequencies()` counts how often each state occurs across all
#' trajectories, and returns one row per state, most frequent first.
#' Missing values are not counted.
#'
#' @param data The trajectories: a list of state vectors, a data frame or
#'   character matrix with one row per trajectory and one column per time
#'   point, a `tna` model or a `netobject` built from sequences. A data
#'   frame needs position columns named as a numbered run (`T1`, `T2`, ...);
#'   any other frame, a one-hot or association network and a 0/1 matrix
#'   raise an error of class `nestimate_bad_input`.
#' @return A data frame with one row per distinct state, sorted by `count`
#'   (descending), with columns `state`, `count` and `proportion` (the
#'   state's share of all observed states, rounded to 4 decimal places).
#'
#' @examples
#' trajs <- list(c("A","B","C"), c("A","B","A"))
#' state_frequencies(trajs)
#'
#' @export
state_frequencies <- function(data) {
  .refuse_onehot_network(data, "state_frequencies()")
  .refuse_binary_matrix(data, "state_frequencies()")
  if (is.data.frame(data)) {
    .refuse_undeclared_roles(data, "state_frequencies()", via_network = TRUE,
                             whole_frame = TRUE)
  }
  data <- .coerce_sequence_input(data)
  if (is.data.frame(data)) {
    all_states <- unlist(lapply(data, as.character), use.names = FALSE)
  } else {
    all_states <- unlist(lapply(data, as.character))
  }

  tbl <- sort(table(all_states), decreasing = TRUE)
  data.frame(
    state = names(tbl),
    count = as.integer(tbl),
    proportion = round(as.numeric(tbl) / sum(tbl), 4),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' @rdname build_mogen
#' @param x For the `print()` and `plot()` methods: a `net_mogen` object.
#' @param ... For the `plot()` method: further arguments passed to
#'   [graphics::plot()]. For the `print()` and `summary()` methods: ignored.
#' @param object For the `summary()` method: a `net_mogen` object.
#' @param type For the `plot()` method: `"ic"` (default) draws AIC and BIC
#'   against the order; `"likelihood"` draws the log-likelihood against the
#'   order. Both mark the selected order with a dashed line.
#' @export
print.net_mogen <- function(x, ...) {
  cat("Multi-Order Generative Model (MOGen)\n")
  cat(sprintf("  Optimal order:  %d (by %s)\n", x$optimal_order, x$criterion))
  cat(sprintf("  Orders tested:  0 to %d\n", max(x$orders)))
  cat(sprintf("  States:         %d\n", length(x$states)))
  cat(sprintf("  Paths:          %d (%d observations)\n",
              x$n_paths, x$n_observations))

  aic_min <- which.min(x$aic) - 1L
  bic_min <- which.min(x$bic) - 1L

  fmt_ic <- function(vals, opt_idx) {
    parts <- sprintf("%.1f", vals)
    parts[opt_idx] <- paste0(parts[opt_idx], "*")
    paste(parts, collapse = " | ")
  }

  cat(sprintf("  AIC:            %s\n", fmt_ic(x$aic, aic_min + 1L)))
  cat(sprintf("  BIC:            %s\n", fmt_ic(x$bic, bic_min + 1L)))
  if (aic_min == bic_min) {
    cat(sprintf("  (* = minimum;  AIC and BIC agree on order %d)\n", aic_min))
  } else {
    cat(sprintf("  (* = minimum;  AIC favors order %d, BIC favors order %d)\n",
                aic_min, bic_min))
  }
  invisible(x)
}

#' @rdname build_mogen
#' @export
summary.net_mogen <- function(object, ...) {
  cat("Multi-Order Generative Model (MOGen) Summary\n\n")
  cat(sprintf("  States: %s\n", paste(object$states, collapse = ", ")))
  cat(sprintf("  Paths: %d | Observations: %d\n",
              object$n_paths, object$n_observations))

  aic_min <- which.min(object$aic) - 1L
  bic_min <- which.min(object$bic) - 1L
  cat(sprintf("  Best by AIC: order %d  |  Best by BIC: order %d\n",
              aic_min, bic_min))
  cat(sprintf("  Selected:    order %d (by %s)\n\n",
              object$optimal_order, object$criterion))

  res <- data.frame(
    order     = object$orders,
    layer_dof = as.numeric(object$layer_dof),
    cum_dof   = as.numeric(object$dof),
    loglik    = round(object$log_likelihood, 2),
    aic       = round(object$aic, 2),
    bic       = round(object$bic, 2),
    stringsAsFactors = FALSE
  )

  res$best <- ""
  res$best[aic_min + 1L] <- "AIC"
  res$best[bic_min + 1L] <- if (res$best[bic_min + 1L] == "") "BIC" else
    paste0(res$best[bic_min + 1L], "+BIC")
  sel_idx <- object$optimal_order + 1L
  res$selected <- ""
  res$selected[sel_idx] <- "<--"

  res
}

#' @rdname build_mogen
#' @export
plot.net_mogen <- function(x, type = c("ic", "likelihood"), ...) {
  type <- .ho_match_arg(type, "type")

  if (type == "ic") {
    orders <- x$orders
    aic_vals <- x$aic
    bic_vals <- x$bic

    ylim <- range(c(aic_vals, bic_vals))
    plot(orders, aic_vals, type = "b", pch = 19, col = "steelblue",
         xlab = "Markov Order", ylab = "Information Criterion",
         main = "MOGen: Order Selection", ylim = ylim, xaxt = "n", ...)
    graphics::axis(1, at = orders)
    graphics::lines(orders, bic_vals, type = "b", pch = 17, col = "coral")
    graphics::legend("topright", legend = c("AIC", "BIC"),
                     col = c("steelblue", "coral"), pch = c(19, 17),
                     lty = 1, bty = "n")
    graphics::abline(v = x$optimal_order, lty = 2, col = "gray40")
    graphics::text(x$optimal_order, ylim[1], sprintf("k*=%d", x$optimal_order),
                   pos = 4, col = "gray40")
  } else {
    orders <- x$orders
    plot(orders, x$log_likelihood, type = "b", pch = 19, col = "steelblue",
         xlab = "Markov Order", ylab = "Log-Likelihood",
         main = "MOGen: Log-Likelihood by Order", xaxt = "n", ...)
    graphics::axis(1, at = orders)
    graphics::abline(v = x$optimal_order, lty = 2, col = "gray40")
  }

  invisible(x)
}
