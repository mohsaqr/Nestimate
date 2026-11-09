# Mixed Markov Model - latent-class mixture of Markov chains
#
# Discovers subgroups with different transition dynamics via EM.
# Uses precomputed per-sequence transition counts for speed.

# ---------------------------------------------------------------------------
# Core EM
# ---------------------------------------------------------------------------

#' Run one EM from a given initialization
#' @param counts N x K^2 per-sequence transition counts.
#'   Column j = (from-1)*K + to encodes the (from, to) transition pair.
#' @param init_state N-length integer vector of initial states (1..K)
#' @param n_comp Number of mixture components
#' @param max_iter Max iterations
#' @param tol Log-likelihood convergence tolerance
#' @param smooth Laplace smoothing
#' @param K Number of states
#' @param init_posterior Optional N x n_comp starting posterior. When NULL a
#'   uniform-random row-normalised posterior is drawn.
#' @param from_ind Accepted for call-site symmetry with the screen/refine
#'   driver and currently unused: the from-state grouping is recomputed
#'   internally as `rep(seq_len(K), each = K)`.
#' @param cov_df Optional data.frame of covariates. When supplied, the
#'   mixing proportions become covariate-dependent (softmax M-step) instead
#'   of a single vector.
#' @return List with P_all, init_all, pi_mix, posterior, ll, iterations,
#'   converged, and cov_beta (the softmax coefficients, NULL without cov_df)
#' @noRd
.mmm_em <- function(counts, init_state, n_comp, max_iter, tol, smooth, K,
                    init_posterior = NULL, from_ind = NULL, cov_df = NULL) {
  N <- nrow(counts)
  K2 <- K * K

  # Initialize posterior (N x M)
  if (is.null(init_posterior)) {
    post <- matrix(runif(N * n_comp), nrow = N, ncol = n_comp)
    post <- post / .rowSums(post, N, n_comp)
  } else {
    post <- init_posterior
  }

  # Index of "from" state for each K^2 pair (K^2-length integer vector).
  # pair j = (from-1)*K + to, so from = (j-1) %/% K + 1.
  from_idx <- rep(seq_len(K), each = K)

  # Initial state indicator (N x K, sparse: one 1 per row).
  # Used in the M-step only - the E-step indexes log(init_all) directly.
  init_ind <- matrix(0, N, K)
  valid_init <- !is.na(init_state)
  init_ind[cbind(which(valid_init), init_state[valid_init])] <- 1
  init_state_safe <- init_state
  init_state_safe[!valid_init] <- 1L

  ll_prev <- -Inf
  converged <- FALSE
  cov_beta <- NULL
  cov_wts <- NULL
  mstep_converged <- TRUE

  # Pre-compute design matrix for covariates (once)
  X_cov <- NULL
  if (!is.null(cov_df)) {
    X_cov <- stats::model.matrix(~ ., data = cov_df)
  }

  log_lik <- matrix(0, N, n_comp)

  for (iter in seq_len(max_iter)) {
    # ---- M-step ----
    # Compute transition table directly in K^2 x M shape so the E-step
    # can use it without a transpose. Pure-R BLAS one-shot:
    #   counts is N x K^2, post is N x M, so crossprod(counts, post)
    #   = t(counts) %*% post is K^2 x M.
    P_unnorm <- crossprod(counts, post) + smooth        # K^2 x M
    # Sum over from-state groups via 3D array view. P_unnorm has columns
    # indexed by j = (f-1)*K + t (column-major), so reshaping to
    # (K_to, K_from, M) and colSums-ing over the to-axis gives K_from x M.
    dim(P_unnorm) <- c(K, K, n_comp)
    from_sums_t <- colSums(P_unnorm)                    # K x M
    dim(P_unnorm) <- c(K2, n_comp)
    from_sums_t[from_sums_t == 0] <- 1
    # Expand K x M to K^2 x M by repeating each row K times in place.
    P_all <- P_unnorm / from_sums_t[from_idx, , drop = FALSE]   # K^2 x M

    # Initial state probabilities: K x M directly (skip the t() at the end).
    init_unnorm <- crossprod(init_ind, post) + smooth     # K x M
    init_sums <- .colSums(init_unnorm, K, n_comp)         # M
    init_sums[init_sums == 0] <- 1
    init_all <- init_unnorm / rep(init_sums, each = K)    # K x M

    # Mixing proportions
    if (is.null(cov_df)) {
      pi_mix <- .colMeans(post, N, n_comp)
      log_pi <- log(pi_mix + 1e-300)
    } else {
      sm <- .mmm_softmax_mstep(post, X_cov, wts = cov_wts)
      log_pi_mat <- sm$log_pi_mat
      cov_beta <- sm$beta
      cov_wts <- sm$wts
      mstep_converged <- sm$converged
      pi_mix <- .colMeans(exp(log_pi_mat), N, n_comp)
    }

    # ---- E-step ----
    # log_P: K^2 x M. log_init: K x M. log_pi: M (or N x M with covariates).
    log_P <- log(P_all + 1e-300)
    log_init <- log(init_all + 1e-300)

    # Transition log-likelihood: N x M = (N x K^2) %*% (K^2 x M)
    log_lik <- counts %*% log_P
    # Initial-state log-likelihood: each row i needs log_init[init_state[i], ].
    # init_ind %*% log_init reduces to indexed extraction.
    # A sequence with no observed first state carries no initial-state term
    # (uninformative), rather than being scored as starting in state 1.
    init_term <- log_init[init_state_safe, , drop = FALSE]
    init_term[!valid_init, ] <- 0
    log_lik <- log_lik + init_term
    # Mixing
    if (is.null(cov_df)) {
      log_lik <- log_lik + rep(log_pi, each = N)
    } else {
      log_lik <- log_lik + log_pi_mat
    }

    # Log-sum-exp row max: pmax.int is variadic, so a single C call covers
    # any n_comp. Special-cased only because column extraction is unrolled.
    if (n_comp == 2L) {
      log_max <- pmax.int(log_lik[, 1L], log_lik[, 2L])
    } else if (n_comp == 3L) {
      log_max <- pmax.int(log_lik[, 1L], log_lik[, 2L], log_lik[, 3L])
    } else if (n_comp == 4L) {
      log_max <- pmax.int(log_lik[, 1L], log_lik[, 2L], log_lik[, 3L], log_lik[, 4L])
    } else {
      log_max <- do.call(pmax.int, lapply(seq_len(n_comp), function(m) log_lik[, m]))
    }
    exp_lik <- exp(log_lik - log_max)
    row_sums <- .rowSums(exp_lik, N, n_comp)
    row_sums[row_sums == 0] <- 1e-300
    post <- exp_lik / row_sums

    ll <- sum(log_max + log(row_sums))

    if (abs(ll - ll_prev) < tol) {
      converged <- TRUE
      break
    }
    # For covariate models: also check relative LL change (handles
    # perfect separation where LL grows toward 0 but never converges)
    if (!is.null(cov_df) && abs(ll) > 0 &&
        abs((ll - ll_prev) / ll) < tol) {
      converged <- TRUE
      break
    }
    ll_prev <- ll
  }
  # EM convergence also needs the last covariate M-step to have converged.
  converged <- converged && mstep_converged

  list(
    P_all = P_all,
    init_all = init_all,
    pi_mix = pi_mix,
    posterior = post,
    ll = ll,
    iterations = iter,
    converged = converged,
    cov_beta = cov_beta
  )
}

# ---------------------------------------------------------------------------
# Initialization via k-means on count vectors
# ---------------------------------------------------------------------------

#' @noRd
.mmm_init_kmeans <- function(counts, n_comp) {
  N <- nrow(counts)
  if (N <= n_comp) {
    post <- matrix(0, nrow = N, ncol = n_comp)
    for (i in seq_len(N)) post[i, ((i - 1L) %% n_comp) + 1L] <- 1
    return(post)
  }

  km <- tryCatch(
    stats::kmeans(counts, centers = n_comp, nstart = 1L, iter.max = 20L),
    error = function(e) NULL
  )

  if (is.null(km)) {
    post <- matrix(runif(N * n_comp), nrow = N, ncol = n_comp)
    return(post / rowSums(post))
  }

  # Convert cluster assignments to soft posterior (vectorized)
  post <- matrix(0.01 / n_comp, nrow = N, ncol = n_comp)
  post[cbind(seq_len(N), km$cluster)] <- 0.99
  post / rowSums(post)
}

# ---------------------------------------------------------------------------
# Cluster quality metrics
# ---------------------------------------------------------------------------

#' Compute cluster quality metrics from posterior matrix
#' @noRd
.mmm_quality <- function(posterior, assignments, M) {
  N <- nrow(posterior)

  # Average Posterior Probability per class
  avepp <- vapply(seq_len(M), function(m) {
    idx <- which(assignments == m)
    if (length(idx) == 0L) return(NA_real_)
    mean(posterior[idx, m])
  }, numeric(1))

  # Max posterior per sequence - direct pmax.int beats `do.call(pmax,
  # as.data.frame(posterior))` which copies the matrix into a list.
  max_post <- if (M == 2L) {
    pmax.int(posterior[, 1L], posterior[, 2L])
  } else if (M == 3L) {
    pmax.int(posterior[, 1L], posterior[, 2L], posterior[, 3L])
  } else if (M == 4L) {
    pmax.int(posterior[, 1L], posterior[, 2L], posterior[, 3L], posterior[, 4L])
  } else {
    do.call(pmax.int, lapply(seq_len(M), function(m) posterior[, m]))
  }

  # Overall AvePP (mean of max posteriors)
  avepp_overall <- mean(max_post)

  # Entropy (normalized: 0 = perfect, 1 = random)
  ent_raw <- -sum(posterior * log(posterior + 1e-300))
  ent_max <- N * log(M)
  entropy <- if (ent_max > 0) ent_raw / ent_max else 0
  relative_entropy <- 1 - entropy

  # Classification error: the expected share of sequences misclassified by
  # modal assignment, sum_i (1 - max_k P(k | y_i)) / N (Vermunt & Magidson,
  # Latent GOLD 4.0 User's Guide, 2005, sec. on classification statistics).
  # The former share with max posterior < 0.5 was 0 by construction for
  # k = 2, where the largest posterior is always >= 0.5.
  classification_error <- mean(1 - max_post)

  # ICL classification entropy (vectorized via cbind indexing)
  class_ent <- -sum(log(posterior[cbind(seq_len(N), assignments)] + 1e-300))

  list(
    avepp = avepp,
    avepp_overall = avepp_overall,
    entropy = entropy,
    relative_entropy = relative_entropy,
    classification_error = classification_error,
    class_entropy = class_ent
  )
}


#' @noRd
.mmm_component_warnings <- function(models, assignments, k,
                                    tol = sqrt(.Machine$double.eps)) {
  sizes <- tabulate(assignments, nbins = k)
  empty <- which(sizes == 0L)
  if (length(empty) > 0L) {
    warning(
      "MMM produced empty hard-assignment cluster",
      if (length(empty) == 1L) " " else "s ",
      paste(empty, collapse = ", "),
      ". The fitted mixture components should be treated as unstable; ",
      "try fewer clusters, more starts, a fixed seed, or stronger sequence ",
      "signal.",
      call. = FALSE
    )
  }

  pairs <- utils::combn(seq_len(k), 2L, simplify = FALSE)
  duplicated <- vapply(pairs, function(pair) {
    max(abs(models[[pair[1L]]]$weights - models[[pair[2L]]]$weights),
        na.rm = TRUE) <= tol
  }, logical(1L))
  if (any(duplicated)) {
    dup_labels <- vapply(pairs[duplicated], function(pair) {
      paste(pair, collapse = "-")
    }, character(1L))
    warning(
      "MMM produced effectively identical transition components: ",
      paste(dup_labels, collapse = ", "),
      ". The requested k is not supported by the observed transition ",
      "structure in this fit.",
      call. = FALSE
    )
  }

  invisible(list(empty = empty, duplicated = pairs[duplicated]))
}


# ---------------------------------------------------------------------------
# Covariate M-step helpers
# ---------------------------------------------------------------------------

#' Covariate M-step: multinomial logit of the posterior on the covariates
#'
#' The exact M-step for covariate-dependent mixing: `nnet::multinom()` fitted
#' with the N x k posterior as soft class targets maximises
#' sum(post * log(pi)) over the softmax coefficients (first class is the
#' reference). Warm-started from the previous EM iteration's weights.
#'
#' @param post N x k posterior matrix
#' @param X N x p design matrix (with intercept column)
#' @param wts Network weights of the previous fit (`$wts`), or NULL
#' @return List with log_pi_mat (N x k), beta ((k-1) x p, columns named by
#'   `colnames(X)`), wts, and converged (FALSE when multinom hit `maxit`)
#' @noRd
.mmm_softmax_mstep <- function(post, X, wts = NULL) {
  # Standardised columns keep the optimiser on scale; the model is the same,
  # so a covariate shifted by a constant changes the intercept and nothing else.
  std <- .standardise_design(X)
  X_std <- std$X
  # nnet treats a supplied `Wts = NULL` as a weight vector of length 0, so
  # the warm start is passed only when there is one.
  fit <- if (is.null(wts)) {
    nnet::multinom(post ~ X_std - 1, maxit = 500L, trace = FALSE)
  } else {
    nnet::multinom(post ~ X_std - 1, Wts = wts, maxit = 500L, trace = FALSE)
  }
  beta <- stats::coef(fit)
  colnames(beta) <- colnames(X)
  beta <- beta %*% std$back
  list(log_pi_mat = log(stats::fitted(fit) + 1e-300), beta = beta,
       wts = fit$wts, converged = fit$convergence == 0L)
}


# ---------------------------------------------------------------------------
# Data decoding helper
# ---------------------------------------------------------------------------

#' Convert integer-encoded matrix (from tna/cograph) to character labels
#' @noRd
.mmm_decode_int_data <- function(raw_data, states) {
  raw_data <- as.data.frame(raw_data, stringsAsFactors = FALSE)
  if (is.integer(raw_data[[1]]) || is.numeric(raw_data[[1]])) {
    raw_data[] <- lapply(raw_data, function(col) {
      idx <- as.integer(col)
      ifelse(is.na(idx) | idx < 1L | idx > length(states), NA_character_,
             states[idx])
    })
  }
  raw_data
}

# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

#' Fit a Mixed Markov Model
#'
#' @description
#' `build_mmm()` fits a mixture of first-order Markov chains to a set of
#' sequences by expectation-maximisation: it looks for `k` subgroups whose
#' members move between states in different ways, each with its own
#' transition matrix and initial-state distribution. Every sequence gets a
#' posterior probability of belonging to each subgroup and is assigned to the
#' most probable one. Covariates can shape the mixing proportions during
#' estimation or be related to the subgroups afterwards (see
#' `covariate_effect`).
#'
#' EM finds a local maximum, so the fit is started `n_starts` times (one
#' start from k-means on the transition counts, the rest random), each is
#' screened with a short EM run, the most promising are continued for up to
#' `max_iter` iterations, and the one with the highest log-likelihood is
#' kept. `print()` reports whether it converged within
#' `max_iter`; a fit that did not is not a result. The information criteria
#' (BIC, AIC, ICL) compare fits with different `k` -- [compare_mmm()] fits a
#' range of them -- while average posterior probability (AvePP), entropy and
#' classification error describe how cleanly the sequences separate within
#' one fit. None of them picks `k` alone: also check that the subgroups'
#' transition matrices differ in ways that can be interpreted, and that
#' each subgroup holds enough sequences to estimate its matrix.
#'
#' @param data The sequences: a data frame with one row per sequence and one
#'   column per position, a `netobject` or bare `cograph_network` built from
#'   sequences, or a `tna` model (whose stored states are decoded to their
#'   labels). Columns take the same roles as in [build_network()]: an id or a
#'   covariate is not a position, and the columns set aside are named in a
#'   message and kept as metadata. Association networks (`cor`, `pcor`,
#'   `glasso`, `ising`, `mgm`) hold variables, not sequences, and raise an
#'   error of class `nestimate_bad_input`.
#' @param k The number of subgroups (mixture components), a whole number of
#'   at least 2. Default `2`.
#' @param n_starts How many starting configurations to try, a positive whole
#'   number. More starts make it likelier that the best local maximum is
#'   found. Default `50`.
#' @param max_iter The most EM iterations a start may run, a positive whole
#'   number. Default `200`.
#' @param tol The convergence tolerance: EM stops when the log-likelihood
#'   changes by less than this between iterations (with `covariate_effect =
#'   "em"`, also when its relative change does). A finite positive number;
#'   default `1e-6`.
#' @param smooth A pseudo-count added to every transition count before the
#'   matrices are estimated, so a transition unseen in a subgroup gets a small
#'   probability instead of zero. A finite non-negative number; default
#'   `0.01`.
#' @param seed An integer that makes the random starts reproducible.
#'   Default `NULL`. The session's own random stream is left as it was.
#' @param covariates Variables describing each sequence (such as a course or
#'   a grade), given as a column name, a character vector of names, a
#'   one-sided formula or a data frame, as in [build_clusters()]. For a
#'   `netobject`, names are looked up in its metadata first, so a typical call
#'   is `build_mmm(net, k = 3, covariates = "session_label")`. How they enter
#'   the model is set by `covariate_effect`. Default `NULL`.
#' @param covariate_effect How `covariates` enter the model. `"em"`
#'   (default) makes the mixing proportions depend on them during EM, so they
#'   shape the subgroups found and add `(k - 1) * p` parameters for `p`
#'   covariate columns; rows with a missing covariate are dropped first, with
#'   a warning. `"posthoc"` fits a plain mixture on every sequence and relates
#'   the covariates to the subgroups only afterwards, with a multinomial logit,
#'   so neither their values nor their missingness change the subgroups.
#'   Ignored without `covariates`.
#' @param estimator How the after-fit multinomial logit is estimated; it never
#'   affects the EM. `"auto"` (default) uses Firth's penalised likelihood
#'   (`brglm2::brmultinom()`, finite under separation) when a cell of the
#'   subgroup-by-covariate table has fewer than 5 sequences, and
#'   `nnet::multinom()` otherwise. `"firth"` and `"multinom"` force one of
#'   them (`"multinom"` warns about separation risk), and `"chisq"` runs
#'   descriptive tests instead of a logit. [build_clusters()] has the details.
#'   Ignored without `covariates`.
#'
#' @return A `net_mmm` object. `print()` shows the fit statistics,
#'   classification quality and subgroup sizes; [summary()] adds each
#'   subgroup's transition matrix and returns one row per subgroup (mixing
#'   proportion, sequences assigned, mean posterior, AvePP); `plot()` draws
#'   the posterior probabilities or the covariate effects;
#'   [build_network()] turns the fit into one network per subgroup (a
#'   `netobject_group`); [cluster_diagnostics()] and [session_ids()] give the
#'   classification and each sequence's subgroup. The object's fields:
#'   \describe{
#'     \item{data}{The full N-row sequence frame used for estimation.}
#'     \item{models}{List of \code{netobject}s, one per component. Each
#'       component carries the rows assigned to that component in its
#'       \code{$data} slot, while its transition matrix is the EM-estimated
#'       component transition matrix.}
#'     \item{k}{Number of components.}
#'     \item{mixing}{Numeric vector of mixing proportions.}
#'     \item{posterior}{N x k matrix of posterior probabilities.}
#'     \item{assignments}{Integer vector of hard assignments (1..k).}
#'     \item{quality}{List: \code{avepp} (per-class), \code{avepp_overall},
#'       \code{entropy}, \code{relative_entropy},
#'       \code{classification_error}, \code{class_entropy}.}
#'     \item{log_likelihood, BIC, AIC, ICL}{Model fit statistics.}
#'     \item{n_params}{Number of free parameters behind BIC/AIC/ICL. With
#'       \code{covariate_effect = "em"} it grows by \code{(k - 1) * p} for
#'       \code{p} covariate columns.}
#'     \item{iterations, converged}{EM iterations used by the retained fit
#'       and whether it met \code{tol}.}
#'     \item{states}{Character vector of state names.}
#'     \item{n_sequences}{Number of sequences actually fitted (rows with
#'       missing covariates are dropped under
#'       \code{covariate_effect = "em"}).}
#'     \item{covariates}{The post-hoc covariate analysis (list), or NULL
#'       when \code{covariates = NULL}.}
#'     \item{network_method, build_args, htna_partition}{Provenance kept
#'       from \code{netobject} / HTNA input so per-cluster networks can be
#'       rebuilt the same way; NULL otherwise.}
#'     \item{metadata}{The \code{netobject}'s per-sequence metadata, one row
#'       per fitted sequence in the row order of \code{posterior}, so
#'       \code{\link{session_ids}} can name each sequence; NULL for other
#'       input.}
#'   }
#'
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' mmm <- build_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
#' mmm
#' \donttest{
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
#'   V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
#' )
#' mmm <- build_mmm(seqs, k = 2, seed = 42)
#' print(mmm)
#' summary(mmm)
#' }
#'
#' @section Initial states:
#' Each sequence's initial state is the value in its first position column,
#' read as it is: the function does not look ahead to the first observed
#' state, and it does not convert missing-value symbols as [build_clusters()]
#' does. The states are the distinct non-missing values of all position
#' columns, so a placeholder such as `"*"` or `"%"` for an empty cell becomes
#' a state, and a sequence starting with it starts in that state. Recode such
#' placeholders to `NA` first: a missing first position gives the sequence
#' no initial state, which leaves its subgroup to the transitions alone. Or
#' trim leading empty cells so every first position holds an observed state.
#'
#' A one-hot network is different: its units are windows of codes, not
#' positions, so its initial state is the one code active in the first
#' window, and none when that window holds several codes or none.
#'
#' @references
#' Helske, S., & Helske, J. (2019). Mixture hidden Markov models for
#' sequence data: The seqHMM package in R. \emph{Journal of Statistical
#' Software}, 88(3), 1--32. \doi{10.18637/jss.v088.i03}
#'
#' Vermunt, J. K., & Magidson, J. (2005). \emph{Latent GOLD 4.0 User's
#' Guide}. Statistical Innovations.
#'
#' @seealso \code{\link{compare_mmm}}, \code{\link{build_network}}
#'
#' @export
build_mmm <- function(data,
                      k = 2L,
                      n_starts = 50L,
                      max_iter = 200L,
                      tol = 1e-6,
                      smooth = 0.01,
                      seed = NULL,
                      covariates = NULL,
                      covariate_effect = c("em", "posthoc"),
                      estimator = c("auto", "firth", "multinom", "chisq")) {

  estimator <- match.arg(estimator)
  covariate_effect <- match.arg(covariate_effect)
  .assert(
    "'k' must be a single numeric value" = is.numeric(k) && length(k) == 1L,
    "'n_starts' must be a single numeric value" =
      is.numeric(n_starts) && length(n_starts) == 1L,
    "'max_iter' must be a single numeric value" =
      is.numeric(max_iter) && length(max_iter) == 1L,
    "'tol' must be a single numeric value" =
      is.numeric(tol) && length(tol) == 1L,
    "'smooth' must be a single numeric value" =
      is.numeric(smooth) && length(smooth) == 1L
  )
  if (!is.finite(k) || k != floor(k) || k < 2L) {
    stop(errorCondition("'k' must be a whole finite number >= 2.",
                        class = "nestimate_bad_input", call = NULL))
  }
  if (!is.finite(n_starts) || n_starts != floor(n_starts) ||
      n_starts < 1L) {
    stop(errorCondition("'n_starts' must be a positive whole finite number.",
         class = "nestimate_bad_input", call = NULL))
  }
  if (!is.finite(max_iter) || max_iter != floor(max_iter) ||
      max_iter < 1L) {
    stop(errorCondition("'max_iter' must be a positive whole finite number.",
         class = "nestimate_bad_input", call = NULL))
  }
  if (!is.finite(tol) || tol <= 0) {
    stop(errorCondition("'tol' must be a finite positive number.",
                        class = "nestimate_bad_input", call = NULL))
  }
  if (!is.finite(smooth) || smooth < 0) {
    stop(errorCondition("'smooth' must be a finite non-negative number.",
                        class = "nestimate_bad_input", call = NULL))
  }
  k <- as.integer(k)
  n_starts <- as.integer(n_starts)
  max_iter <- as.integer(max_iter)

  # ---- Extract data and states ----
  htna_partition <- .capture_htna_partition(data)
  network_method <- NULL
  build_args     <- NULL
  metadata       <- NULL
  if (inherits(data, "tna") || inherits(data, "ftna")) {
    raw_data <- data$data
    states <- data$labels
    if (is.null(raw_data)) {
      stop(errorCondition(
        paste0(
          "tna model does not contain $data. ",
          "Rebuild with tna::tna(data) where data is kept."),
        class = "nestimate_bad_input", call = NULL)) # nocov end
    }
    # A windowed model (tna::import_onehot(), windowed sequences) counts within
    # windows whose layout its decoded data does not carry; fitted as plain
    # positions it would count across slots that are never adjacent.
    if (isTRUE(attr(data, "params")$windowed)) {
      stop(errorCondition(paste(
        "build_mmm(): this tna model is windowed, and its windows cannot be",
        "replayed. Build the network with build_network(format = \"onehot\",",
        "actor = ) or pass the sequences."),
        class = "nestimate_bad_input", call = NULL))
    }
    raw_data <- .mmm_decode_int_data(raw_data, states)
  } else if (inherits(data, "netobject_group")) {
    stop(errorCondition(paste(
      "build_mmm(): a netobject_group holds several networks, not one set of",
      "sequences. Fit each network, or pass the sequences of all groups."),
      class = "nestimate_bad_input", call = NULL))
  } else if (inherits(data, c("netobject", "cograph_network"))) {
    # A bare cograph_network is coerced first, so it meets the same refusals
    # and the same one-hot unit rule as a netobject. Its declared method (a
    # bare network keeps $method) rebuilds the clusters and decides the
    # association refusal; an inferred one is not a declaration.
    bare <- !inherits(data, "netobject")
    declared_method <- data$method %||% data$meta$tna$method
    if (bare && is.null(data$method) && identical(data$meta$source, "tna") &&
        is.null(.tna_conversion_counts(data))) {
      # cograph's conversion of a tna model keeps neither its type nor its
      # windows; its positions are read only when they reproduce its weights
      stop(errorCondition(paste(
        "build_mmm(): this network was converted from a tna model, and its",
        "sequences do not reproduce its weights (windows, scaling or other",
        "settings were lost). Pass the tna model itself or the sequences."),
        class = "nestimate_bad_input", call = NULL))
    }
    .mmm_refuse_association(list(method = declared_method))
    if (bare) data <- .as_netobject(data)
    .mmm_refuse_association(data)
    network_method <- if (bare) declared_method else data$method
    build_args     <- if (bare) NULL else data$build_args
    # the estimator's own settings (n_gram, start/end states, concatenation)
    # rebuild the clusters as the source network was built
    if (!bare && length(data$params) > 0L && is.null(build_args$params)) {
      build_args$params <- data$params
    }
    metadata       <- data$metadata
    raw_data <- data$data
    states <- data$nodes$label
    # The rows of a one-hot network are time points, not sequences.
    if (.is_onehot_params(data$params)) {
      raw_data <- .onehot_unit_sequences(data)
      build_args <- NULL
      metadata <- NULL
      # its units are windowed state sequences now, which the wtna methods
      # (they read 0/1 codes) cannot rebuild; their counting methods can
      if (!is.null(network_method)) {
        network_method <- .resampling_method(network_method, data$params)
      }
    }
  } else {
    raw_data <- as.data.frame(data)
    state_cols <- .select_state_cols(raw_data, id = NULL, cols = NULL)
    mat <- as.matrix(raw_data[, state_cols, drop = FALSE])
    states <- sort(unique(as.vector(mat)))
    states <- states[!is.na(states)]
  }

  n_states <- length(states)
  .assert("Need at least 2 states" = n_states >= 2L)

  # ---- Resolve covariates ----
  cov_df <- NULL
  if (!is.null(covariates)) {
    cov_resolved <- .resolve_covariates(covariates, data, nrow(raw_data))
    cov_df <- cov_resolved$cov_df
  }

  # Strip covariate columns from raw_data before state extraction
  if (!is.null(cov_df) && is.data.frame(raw_data)) {
    cov_col_names <- names(cov_df)
    keep <- setdiff(names(raw_data), cov_col_names)
    if (length(keep) > 0L) {
      raw_data <- .keep_window_attrs(raw_data[, keep, drop = FALSE], raw_data)
    }
    # Re-extract states after removing covariate columns
    state_cols_clean <- .select_state_cols(raw_data, id = NULL, cols = NULL)
    mat <- as.matrix(raw_data[, state_cols_clean, drop = FALSE])
    states <- sort(unique(as.vector(mat)))
    states <- states[!is.na(states)]
    n_states <- length(states)
    .assert("Need at least 2 states" = n_states >= 2L)
  }

  # Same roles as build_network(): an id or a covariate is not a position.
  if (is.data.frame(data)) {
    roles <- .split_sequence_columns(raw_data, "build_mmm")
    if (!is.null(roles$metadata)) {
      raw_data <- roles$sequences
      metadata <- roles$metadata
      states <- sort(unique(as.vector(as.matrix(raw_data))))
      n_states <- length(states)
      .assert("Need at least 2 states" = n_states >= 2L)
    }
  }

  # ---- Precompute per-sequence transition counts ----
  # Windowed rows (one-hot units, sequence windows) count within their windows,
  # as the network they came from does.
  counts <- if (isTRUE(attr(raw_data, "windowed"))) {
    .precompute_per_sequence_windowed(raw_data, method = "relative",
                                      params = list(), states = states)
  } else {
    .precompute_per_sequence_wide(
      raw_data, method = "relative", cols = NULL, id_col = NULL, states = states
    )
  }
  N <- nrow(counts)
  K2 <- n_states * n_states

  # Extract initial states from the FIRST sequence column. We do NOT scan
  # forward to the first non-NA position - see "Initial states" section in
  # the build_mmm() roxygen for the rationale and a workaround for
  # left-padded sequences. NA values here are treated downstream as an
  # uninformative initial distribution.
  state_cols <- .select_state_cols(raw_data, id = NULL, cols = NULL)
  first_col <- as.character(raw_data[[state_cols[1L]]])
  if (isTRUE(attr(raw_data, "windowed")) && !is.null(attr(raw_data, "codes"))) {
    # A one-hot window's slots are codes, not positions: the first state is
    # the one state observed in the first window, none when it holds several.
    width <- as.integer(attr(raw_data, "window_size") %||% 2L) *
      as.integer(attr(raw_data, "window_span") %||% 1L)
    first <- as.matrix(raw_data[, state_cols[seq_len(min(width,
                                                         length(state_cols)))],
                                drop = FALSE])
    first_col <- apply(first, 1L, \(row) {
      seen <- unique(row[!is.na(row)])
      if (length(seen) == 1L) as.character(seen) else NA_character_
    })
  }
  init_state <- match(first_col, states)

  # Handle covariate NAs: subset counts, init_state, cov_df. Done only in
  # "em" mode, so covariate missingness never changes which sequences get
  # clustered when covariates are meant to be post-hoc only. In "posthoc"
  # mode the after-fit logit (.run_covariate_analysis) does its own
  # complete-case subsetting on the full data instead.
  if (!is.null(cov_df) && covariate_effect == "em") {
    complete <- stats::complete.cases(cov_df)
    n_dropped <- sum(!complete)
    if (n_dropped > 0L) {
      warning(sprintf(
        "Dropped %d rows with NA covariates (%d remaining).",
        n_dropped, sum(complete)
      ), call. = FALSE)
      counts <- counts[complete, , drop = FALSE]
      init_state <- init_state[complete]
      cov_df <- cov_df[complete, , drop = FALSE]
      raw_data <- raw_data[complete, , drop = FALSE]
      if (!is.null(metadata)) metadata <- metadata[complete, , drop = FALSE]
      N <- nrow(counts)
    }
  }

  # em_cov_df is the covariate matrix the EM actually sees (covariate-
  # dependent mixing); NULL in "posthoc" mode so the mixture is plain.
  # cov_df remains the full set handed to the post-hoc covariate analysis.
  em_cov_df <- if (covariate_effect == "em") cov_df else NULL

  .assert("Need more sequences than components" = N > k)

  # ---- Run EM with screen-then-refine strategy ----
  if (!is.null(seed)) .local_seed(seed)

  # Pre-compute from-state indicator (shared across all EM runs)
  from_idx <- rep(seq_len(n_states), each = n_states)
  from_ind <- matrix(0, K2, n_states)
  from_ind[cbind(seq_len(K2), from_idx)] <- 1

  # Screen phase
  screen_iter <- min(25L, max_iter)
  n_refine <- max(3L, n_starts %/% 10L)

  # Collect all initializations (before forking for seed control)
  inits <- vector("list", n_starts)
  inits[[1L]] <- .mmm_init_kmeans(counts, k)
  for (s in seq_len(n_starts - 1L)) {
    inits[[s + 1L]] <- matrix(runif(N * k), nrow = N, ncol = k)
    inits[[s + 1L]] <- inits[[s + 1L]] / rowSums(inits[[s + 1L]])
  }

  # Screen: short EM runs (parallel on Unix)
  n_cores <- 1L
  if (.Platform$OS.type == "unix" && n_starts > 1L) {
    detected <- parallel::detectCores(logical = FALSE)
    n_cores <- min(if (is.finite(detected)) detected else 1L, n_starts)
    if (isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE")))) { # nocov start
      n_cores <- 1L
    } # nocov end
  }

  .run_screen <- function(init) {
    # Screen phase: no covariates (fast), just find good transition inits
    .mmm_em(counts, init_state, k, screen_iter, tol, smooth, n_states,
            init_posterior = init, from_ind = from_ind, cov_df = NULL)
  }

  if (n_cores > 1L) {
    screen_results <- parallel::mclapply(inits, .run_screen,
                                          mc.cores = n_cores)
  } else {
    screen_results <- lapply(inits, .run_screen)
  }

  # Handle failed parallel runs
  screen_lls <- vapply(screen_results, function(r) {
    if (is.null(r)) -Inf else r$ll # nocov
  }, numeric(1))
  screen_conv <- vapply(screen_results, function(r) {
    if (is.null(r)) FALSE else r$converged # nocov
  }, logical(1))

  # Refine: top candidates + any that converged early
  top_idx <- order(screen_lls, decreasing = TRUE)[seq_len(n_refine)]
  conv_idx <- which(screen_conv)
  refine_idx <- unique(c(top_idx, conv_idx))

  best <- list(ll = -Inf)
  for (idx in refine_idx) {
    if (is.null(screen_results[[idx]])) next
    # Always refine when covariates present (screen ran without them)
    if (screen_results[[idx]]$converged && is.null(em_cov_df)) {
      run <- screen_results[[idx]]
    } else {
      run <- .mmm_em(counts, init_state, k, max_iter, tol, smooth, n_states,
                     init_posterior = screen_results[[idx]]$posterior,
                     from_ind = from_ind, cov_df = em_cov_df)
    }
    if (run$ll > best$ll) best <- run
  }

  # Fallback: if all parallel screen runs returned NULL (e.g. macOS arm64 fork
  # failures), run one guaranteed sequential EM before touching best$P_all.
  if (is.null(best$P_all)) {
    run <- .mmm_em(counts, init_state, k, max_iter, tol, smooth, n_states,
                   init_posterior = inits[[1L]], from_ind = from_ind,
                   cov_df = em_cov_df)
    if (!is.null(run)) best <- run
  }
  if (is.null(best$P_all)) {
    stop("EM algorithm failed to produce a valid result in all starting configurations.",
         call. = FALSE)
  }

  # ---- Assignments ----
  # max.col is the vectorized row-wise argmax - avoids `apply` overhead.
  assignments <- max.col(best$posterior, ties.method = "first")

  # ---- Build netobjects for each component ----
  models <- lapply(seq_len(k), function(m) {
    P_vec <- best$P_all[, m]
    P_mat <- t(matrix(P_vec, nrow = n_states, ncol = n_states))
    dimnames(P_mat) <- list(states, states)

    edges <- .extract_edges_from_matrix(P_mat, directed = TRUE)

    nodes_df <- data.frame(
      id = seq_along(states), label = states, name = states,
      x = NA_real_, y = NA_real_, stringsAsFactors = FALSE
    )

    # Initial state probabilities from EM M-step (states x components matrix)
    init <- setNames(best$init_all[, m], states)
    model_data <- raw_data[assignments == m, , drop = FALSE]

    structure(list(
      data = model_data,
      weights = P_mat,
      nodes = nodes_df,
      edges = edges,
      directed = TRUE,
      method = "relative",
      params = list(),
      scaling = NULL,
      threshold = 0,
      n_nodes = n_states,
      n_edges = nrow(edges),
      level = NULL,
      initial = init,
      meta = list(source = "nestimate", layout = NULL,
                  tna = list(method = "relative")),
      node_groups = NULL
    ), class = c("netobject", "cograph_network"))
  })

  names(models) <- paste0("Cluster ", seq_len(k))

  # ---- Quality ----
  quality <- .mmm_quality(best$posterior, assignments, k)
  .mmm_component_warnings(models, assignments, k)

  # ---- Information criteria ----
  n_params <- k * n_states * (n_states - 1L) +
              k * (n_states - 1L) +
              # mixing: one free logit per non-reference class and design
              # column (a 4-level factor is 3 columns, not 1)
              if (is.null(em_cov_df)) (k - 1L) else
                (k - 1L) * ncol(stats::model.matrix(~ ., data = em_cov_df))
  BIC_val <- -2 * best$ll + n_params * log(N)
  AIC_val <- -2 * best$ll + 2 * n_params
  ICL_val <- BIC_val + 2 * quality$class_entropy

  # ---- Covariate output ----
  cov_result <- NULL
  if (!is.null(cov_df)) {
    # Post-hoc multinomial logit for SEs / ORs. Runs in both modes: in
    # "em" mode it describes the covariate-dependent mixing the EM fit; in
    # "posthoc" mode it is the only place covariates enter at all.
    cov_result <- .run_covariate_analysis(
      assignments, cov_df, paste(names(cov_df), collapse = " + "), k,
      estimator = estimator
    )
    # Store EM-estimated mixing beta when the EM produced one (em mode only).
    if (!is.null(best$cov_beta)) cov_result$beta <- best$cov_beta
  }

  structure(list(
    data = raw_data,
    models = models,
    k = k,
    mixing = best$pi_mix,
    posterior = best$posterior,
    assignments = assignments,
    quality = quality,
    log_likelihood = best$ll,
    BIC = BIC_val,
    AIC = AIC_val,
    ICL = ICL_val,
    n_params = n_params,
    iterations = best$iterations,
    converged = best$converged,
    states = states,
    n_sequences = N,
    covariates = cov_result,
    covariate_effect = if (is.null(cov_result)) NULL else covariate_effect,
    network_method = network_method,
    build_args     = build_args,
    metadata       = metadata,
    htna_partition = htna_partition
  ), class = "net_mmm")
}

# ---------------------------------------------------------------------------
# Compare multiple k values
# ---------------------------------------------------------------------------

#' Compare MMM fits across different k
#'
#' @description
#' `compare_mmm()` fits a mixed Markov model with [build_mmm()] for each
#' number of subgroups in `k` and returns one row per fit with its
#' log-likelihood, information criteria (AIC, BIC, ICL), overall average
#' posterior probability (AvePP), entropy and convergence. Every fit uses
#' the same data and the same [build_mmm()] settings.
#'
#' Lower AIC, BIC and ICL indicate more support relative to the other fits.
#' ICL adds a penalty for uncertain classification to BIC, so it favours
#' subgroups that separate cleanly (Biernacki, Celeux & Govaert, 2000);
#' AvePP and entropy describe that separation within each fit. The criteria
#' do not settle `k` on their own: check that every fit converged, that the
#' subgroups' transition matrices of the chosen fit differ in interpretable
#' ways (`summary()` of that fit prints them), and, when the criteria keep
#' falling at the largest `k`, extend the range.
#'
#' @param data The sequences, in any form [build_mmm()] accepts: a wide data
#'   frame, a `netobject` or `cograph_network` built from sequences, or a
#'   `tna` model.
#' @param k The numbers of subgroups to fit, a non-empty vector of whole
#'   numbers of at least 2; one model is fitted per value. Default `2:5`.
#' @param return_fits Whether to keep the fitted models. With `TRUE`, the
#'   `net_mmm` fits are attached to the result as its `"fits"` attribute, a
#'   list named by `k`, so the chosen model needs no second run of the EM.
#'   Default `FALSE`, which returns only the comparison table.
#' @inheritParams build_mmm
#' @param seed A seed for reproducible starts, applied afresh to each fit
#'   (every `k` starts from the same random state). Default `NULL`.
#' @param ... For the `print()`, `summary()` and `plot()` methods: not used;
#'   any extra argument raises an error.
#' @return An `mmm_compare` object, a data frame with one row per value of
#'   `k`. `print()` shows the table with a `best`
#'   column marking the lowest BIC (`"<-- BIC"`) and, when it falls on a
#'   different row, the lowest ICL (`"<-- ICL"`), and returns that marked
#'   table invisibly. `summary()` returns the same table as a plain data
#'   frame with a `best` column holding `"BIC"`, `"ICL"` or `""`; when one
#'   row has both minima it is marked `"BIC"` only. `plot()` draws BIC, AIC
#'   and ICL against `k`, and returns the `ggplot` invisibly after printing
#'   it. The table's columns are `k`, `log_likelihood`, `AIC`, `BIC`,
#'   `ICL`, `AvePP`, `Entropy` and `converged`; `return_fits = TRUE` adds the
#'   `"fits"` attribute.
#'
#' @references
#' Biernacki, C., Celeux, G., & Govaert, G. (2000). Assessing a mixture
#' model for clustering with the integrated completed likelihood.
#' \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence},
#' 22(7), 719--725. \doi{10.1109/34.865189}
#'
#' @seealso [build_mmm()], [cluster_mmm()]
#'
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' comp <- compare_mmm(seqs, k = 2:3, n_starts = 1, max_iter = 10, seed = 1)
#' comp
#' \donttest{
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
#'   V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
#' )
#' comp <- compare_mmm(seqs, k = 2:3, seed = 42)
#' print(comp)
#'
#' # Retain the fits so the chosen model needs no re-run; summary() marks
#' # the minimum-BIC and minimum-ICL rows in its `best` column.
#' comp_with_fits <- compare_mmm(seqs, k = 2:3, seed = 42, return_fits = TRUE)
#' summary(comp_with_fits)
#' }
#'
#' @export
compare_mmm <- function(data, k = 2:5, return_fits = FALSE, n_starts = 50L,
                        max_iter = 200L, tol = 1e-6, smooth = 0.01,
                        seed = NULL, covariates = NULL,
                        covariate_effect = c("em", "posthoc"),
                        estimator = c("auto", "firth", "multinom", "chisq")) {
  covariate_effect <- match.arg(covariate_effect)
  estimator <- match.arg(estimator)
  if (!is.logical(return_fits) || length(return_fits) != 1L ||
      is.na(return_fits)) {
    stop(errorCondition("'return_fits' must be TRUE or FALSE.",
                        class = "nestimate_bad_input", call = NULL))
  }
  if (!is.numeric(k) || !length(k) || any(!is.finite(k)) ||
      any(k != floor(k)) || any(k < 2L)) {
    stop(errorCondition(
      "'k' must be a non-empty vector of whole finite numbers >= 2.",
      class = "nestimate_bad_input", call = NULL))
  }
  k <- as.integer(k)
  fits <- lapply(k, function(m) {
    build_mmm(data, k = m, n_starts = n_starts, max_iter = max_iter,
              tol = tol, smooth = smooth, seed = seed, covariates = covariates,
              covariate_effect = covariate_effect, estimator = estimator)
  })
  rows <- lapply(seq_along(k), function(i) {
    fit <- fits[[i]]
    data.frame(
      k              = k[i],
      log_likelihood = fit$log_likelihood,
      AIC            = fit$AIC,
      BIC            = fit$BIC,
      ICL            = fit$ICL,
      AvePP          = fit$quality$avepp_overall,
      Entropy        = fit$quality$entropy,
      converged      = fit$converged,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  class(result) <- c("mmm_compare", "data.frame")
  if (return_fits) {
    names(fits) <- as.character(k)
    attr(result, "fits") <- fits
  }
  result
}


# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' @rdname build_mmm
#' @param x For the \code{print()} and \code{plot()} methods: a
#'   \code{net_mmm} fit, or the \code{net_mmm_clustering} attribute of the
#'   \code{netobject_group} built from one.
#' @param digits For the \code{print()} methods: decimal places of the
#'   statistics. Default \code{3}.
#' @param ... For the \code{print()}, \code{summary()} and \code{plot()}
#'   methods: not used; any extra argument raises an error.
#' @param object For the \code{summary()} method: a \code{net_mmm} fit.
#' @param type For the \code{plot()} methods: \code{"posterior"} (default), a
#'   histogram of each sequence's largest posterior probability, coloured by
#'   cluster; or \code{"covariates"}, a forest plot of the covariate effects,
#'   for a model fitted with \code{covariates}. On a
#'   \code{net_mmm_clustering}, \code{"predictors"} is an alias of
#'   \code{"covariates"}.
#' @param combined For the \code{plot()} methods with
#'   \code{type = "covariates"}: \code{TRUE} (default) draws the forest panels
#'   as one faceted plot; \code{FALSE} returns a list of ggplots, one per
#'   cluster.
#' @return The `print()` methods return their input invisibly. For a model
#'   fitted with `covariates`, `summary()` also prints the covariate profiles
#'   and effects and returns, invisibly, a `tidy_covariates` data frame of
#'   the effects (one row per subgroup other than the reference and per
#'   predictor term, so a factor gives one row per non-reference level; one
#'   row per covariate test with `estimator = "chisq"`), with the
#'   per-subgroup table as its `cluster_stats` attribute. The `plot()`
#'   methods print the plot and return the `ggplot` invisibly; with
#'   `type = "covariates"` and `combined = FALSE` they return, invisibly and
#'   without printing, a list of `ggplot`s named by subgroup.
#' @export
print.net_mmm <- function(x, digits = 3L, ...) {
  .net_mmm_check_unused_dots("print.net_mmm", ...)
  .net_mmm_check_digits(digits)
  digits <- as.integer(digits)
  k <- as.integer(x$k)
  n_total <- as.integer(x$n_sequences)

  cat("Mixed Markov Model\n")
  cat(sprintf("  Sequences: %d  |  Clusters: %d  |  States: %d\n",
              n_total, k, length(x$states)))
  cat(sprintf("  ICs: LL = %.*f  |  BIC = %.*f  |  AIC = %.*f  |  ICL = %.*f\n",
              digits, x$log_likelihood, digits, x$BIC,
              digits, x$AIC, digits, x$ICL))
  if (!is.null(x$quality)) {
    cat(sprintf(
      "  Quality: AvePP = %.*f  |  Entropy = %.*f  |  Class.Err = %.1f%%\n",
      digits, x$quality$avepp_overall, digits, x$quality$entropy,
      x$quality$classification_error * 100))
  }
  if (isFALSE(x$converged)) {
    cat(sprintf("  Status: did not converge in %d iterations\n",
                as.integer(x$iterations)))
  }

  cat("\n")
  sizes <- as.integer(tabulate(x$assignments, nbins = k))
  cols <- list(
    Cluster = sprintf("%d", seq_len(k)),
    N       = .fmt_size_pct(sizes, n_total),
    `Mix%`  = sprintf("%4.1f%%", as.numeric(x$mixing) * 100),
    AvePP   = sprintf(paste0("%.", digits, "f"),
                      as.numeric(x$quality$avepp))
  )
  cat(paste(.cluster_table_lines(cols), collapse = "\n"), "\n", sep = "")

  if (!is.null(x$covariates)) {
    cov_names <- setdiff(
      unique(x$covariates$coefficients$variable), "(Intercept)"
    )
    cat(sprintf("\n  Covariates: %s (%s, %d predictors)\n",
                paste(cov_names, collapse = ", "),
                .mmm_covariate_label(x$covariate_effect), length(cov_names)))
  }

  invisible(x)
}

# Label for how covariates entered the fit. Objects saved before
# covariate_effect was stored carry NULL; they were fitted with the EM default.
#' @noRd
.mmm_covariate_label <- function(effect) {
  if (identical(effect, "posthoc")) "post-hoc" else "integrated"
}

#' @rdname build_mmm
#' @export
summary.net_mmm <- function(object, ...) {
  .net_mmm_check_unused_dots("summary.net_mmm", ...)
  print(object)
  cat("\n")
  for (m in seq_len(object$k)) {
    n_in <- sum(object$assignments == m)
    cat(sprintf("--- Cluster %d (%.1f%%, n=%d) ---\n",
                m, object$mixing[m] * 100, n_in))
    print(round(object$models[[m]]$weights, 3))
    cat("\n")
  }

  # Covariate analysis
  if (!is.null(object$covariates)) {
    cov <- object$covariates
    .print_covariate_profiles(
      cov, if (identical(object$covariate_effect, "posthoc")) {
        "Covariate Analysis (post-hoc -- does not influence cluster membership)"
      } else {
        "Covariate Analysis (integrated into EM -- influences cluster membership)"
      }
    )
  }

  k <- object$k
  mean_posterior <- vapply(seq_len(k), function(m) {
    idx <- which(object$assignments == m)
    if (length(idx) == 0L) return(NA_real_)
    mean(object$posterior[idx, m])
  }, numeric(1L))

  component_stats <- data.frame(
    component      = seq_len(k),
    prior          = as.numeric(object$mixing),
    n_assigned     = as.integer(tabulate(object$assignments, nbins = k)),
    mean_posterior = mean_posterior,
    avepp          = as.numeric(object$quality$avepp),
    stringsAsFactors = FALSE,
    row.names      = NULL
  )

  if (!is.null(object$covariates)) {
    return(invisible(.tidy_covariates(object$covariates,
                                       cluster_stats = component_stats)))
  }
  component_stats
}

#' @rdname build_mmm
#' @export
plot.net_mmm <- function(x, type = c("posterior", "covariates"),
                          combined = TRUE, ...) {
  .net_mmm_check_unused_dots("plot.net_mmm", ...)
  type <- match.arg(type)
  .assert(is.logical(combined), length(combined) == 1L, !is.na(combined))

  if (type == "covariates") {
    if (is.null(x$covariates)) {
      stop(errorCondition(
        "No covariate analysis found. Run build_mmm() with covariates.",
        class = "nestimate_bad_input", call = NULL))
    }
    return(.plot_covariate_forest(
      x$covariates$coefficients,
      sprintf("Covariate Effects (ref: Cluster %s)",
              x$covariates$fit$reference_cluster),
      combined = combined
    ))
  }

  if (type == "posterior") {
    return(.plot_mmm_posterior(x))
  }
}

#' Plot posterior probability distribution
#' @noRd
.plot_mmm_posterior <- function(x) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required.", call. = FALSE) # nocov
  }

  max_post <- apply(x$posterior, 1L, max)
  df <- data.frame(
    cluster = factor(x$assignments),
    max_posterior = max_post
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = max_posterior, fill = cluster)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
    ggplot2::labs(x = "Max Posterior Probability", y = "Count",
                  title = "Classification Certainty", fill = "Cluster") +
    ggplot2::theme_minimal(base_size = 12)

  print(p)
  invisible(p)
}

#' @rdname compare_mmm
#' @param x For the `print()` and `plot()` methods: an `mmm_compare` table.
#' @param object For the `summary()` method: an `mmm_compare` table.
#' @export
print.mmm_compare <- function(x, ...) {
  .mmm_compare_check_unused_dots(...)
  cat("MMM Model Comparison\n\n")
  best_bic <- which.min(x$BIC)
  best_icl <- which.min(x$ICL)
  x$best <- ""
  x$best[best_bic] <- "<-- BIC"
  if (best_icl != best_bic) x$best[best_icl] <- paste(x$best[best_icl], "<-- ICL")
  print.data.frame(x, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' @rdname compare_mmm
#' @export
summary.mmm_compare <- function(object, ...) {
  .mmm_compare_check_unused_dots(...)
  best_bic <- which.min(object$BIC)
  best_icl <- which.min(object$ICL)
  best <- rep("", nrow(object))
  best[best_bic] <- "BIC"
  if (length(best_icl) && best_icl != best_bic) {
    best[best_icl] <- if (nzchar(best[best_icl]))
      paste(best[best_icl], "ICL", sep = "+") else "ICL"
  }
  out <- as.data.frame(object)
  out$best <- best
  row.names(out) <- NULL
  out
}

#' @rdname compare_mmm
#' @export
plot.mmm_compare <- function(x, ...) {
  .mmm_compare_check_unused_dots(...)
  if (!requireNamespace("ggplot2", quietly = TRUE)) { # nocov start
    stop("Package 'ggplot2' required.", call. = FALSE)
  } # nocov end

  df <- data.frame(
    k = rep(x$k, 3),
    value = c(x$BIC, x$AIC, x$ICL),
    criterion = rep(c("BIC", "AIC", "ICL"), each = nrow(x)),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = value, color = criterion)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_x_continuous(breaks = x$k) +
    ggplot2::labs(x = "k (components)", y = "Information Criterion",
                  title = "MMM Model Selection", color = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  print(p)
  invisible(p)
}

.net_mmm_check_unused_dots <- function(method, ...) {
  dots <- list(...)
  if (!length(dots)) {
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
  stop(errorCondition(
    paste0(
      method, "() got unsupported argument",
      if (length(dots) == 1L) ": " else "s: ",
      paste(dot_names, collapse = ", ")
    ),
    class = "nestimate_bad_input", call = NULL
  ))
}

.net_mmm_check_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits != floor(digits) || digits < 0L) {
    stop(errorCondition(
      "'digits' must be a single non-negative whole finite number.",
      class = "nestimate_bad_input", call = NULL))
  }
  invisible(TRUE)
}

.mmm_compare_check_unused_dots <- function(...) {
  dots <- list(...)
  if (!length(dots)) {
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
  stop(errorCondition(
    paste0(
      "mmm_compare method got unsupported argument",
      if (length(dots) == 1L) ": " else "s: ",
      paste(dot_names, collapse = ", ")
    ),
    class = "nestimate_bad_input", call = NULL
  ))
}

# ---------------------------------------------------------------------------
# MMM network materialization metadata
# ---------------------------------------------------------------------------

# Attach MMM clustering metadata (assignments, posterior, mixing, ICs, full
# sequence data) to a netobject_group built from a net_mmm. Used by
# build_network(net_mmm), including the network-building stage inside
# cluster_network(..., cluster_by = "mmm"). Every member has its own
# $weights/$nodes/$edges, and `attr(, "clustering")` carries N-row $data
# matching $assignments.
#' @noRd
.attach_mmm_clustering <- function(nets, mmm, full_data = NULL) {
  clustering_info <- mmm[setdiff(names(mmm), "models")]
  if (is.null(full_data) && !is.null(mmm$data)) {
    full_data <- mmm$data
  }
  if (is.null(full_data) && length(nets) > 0L) {
    full_data <- nets[[1L]]$data
  }
  if (!is.null(full_data)) clustering_info$data <- full_data
  class(clustering_info) <- "net_mmm_clustering"
  attr(nets, "clustering") <- clustering_info
  attr(nets, "group_col")  <- "cluster"
  .restore_htna_group(nets, mmm$htna_partition)
}

#' Cluster sequences using Mixed Markov Models
#'
#' @description
#' `cluster_mmm()` clusters sequences by fitting a mixture of first-order
#' Markov chains and returns the fitted `net_mmm` object: the same fit, with
#' the same arguments, as [build_mmm()]. The fit holds each sequence's
#' subgroup and posterior probabilities, the mixing proportions, the
#' information criteria and the subgroups' transition matrices.
#'
#' Pass the fit to [build_network()] to get one network per subgroup, or use
#' `cluster_network(..., cluster_by = "mmm")` to fit and build the networks
#' in one call. [build_mmm()] explains how to read the fit and what to check
#' before trusting it.
#'
#' @details
#' `cluster_mmm()` takes the arguments of [build_mmm()] plus `cluster_by`, so
#' its call has the same shape as [cluster_network()].
#'
#' @inheritParams build_mmm
#' @param cluster_by The clustering family. Only `"mmm"` (the default) is
#'   accepted; it is there so `cluster_mmm()` and [cluster_network()] share
#'   one call shape. Any other value raises an error that points to
#'   [cluster_network()].
#' @param ... Not used; any extra argument raises an error.
#' @return A `net_mmm` object, exactly as [build_mmm()] returns it: `print()`,
#'   [summary()] and `plot()` read it, [build_network()] turns it into a
#'   `netobject_group` with one network per subgroup, and
#'   [cluster_diagnostics()] gives the classification and, for sequences
#'   built from long data with `actor`/`session`, [session_ids()] gives each
#'   sequence's subgroup. [build_mmm()] lists its fields. For HTNA input,
#'   the actor partition is kept and restored when the fit is turned into
#'   networks with [build_network()] or [as_htna()].
#' @seealso [build_mmm()], [build_network()], and [cluster_network()] for
#'   fitting and building the per-cluster networks in one call.
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' fit <- cluster_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
#' fit
#' cluster_diagnostics(fit)
#' \donttest{
#' # Visualise with sequence_plot
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:3], 40, TRUE),
#'   V2 = sample(LETTERS[1:3], 40, TRUE),
#'   V3 = sample(LETTERS[1:3], 40, TRUE)
#' )
#' fit <- cluster_mmm(seqs, k = 2)
#' sequence_plot(fit, type = "index")
#' }
#' @export
cluster_mmm <- function(data, k = 2L, n_starts = 50L, max_iter = 200L,
                        tol = 1e-6, smooth = 0.01, seed = NULL,
                        covariates = NULL,
                        covariate_effect = c("em", "posthoc"),
                        estimator = c("auto", "firth", "multinom", "chisq"),
                        cluster_by = "mmm", ...) {
  estimator <- match.arg(estimator)
  covariate_effect <- match.arg(covariate_effect)
  dots <- list(...)
  if (length(dots) > 0L) {
    dot_names <- names(dots)
    dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
    stop(errorCondition(paste0("cluster_mmm() got unsupported argument",
         if (length(dots) == 1L) ": " else "s: ",
         paste(dot_names, collapse = ", ")),
         class = "nestimate_bad_input", call = NULL))
  }

  # cluster_by exists for API parity with cluster_network() so a single
  # surface argument toggles the clustering family. Only "mmm" is valid
  # here; anything else is a programming error worth catching loudly.
  if (!identical(as.character(cluster_by), "mmm")) {
    stop(errorCondition(
      paste0(
        "cluster_mmm() only supports cluster_by = \"mmm\". For other ",
        "clustering algorithms, use cluster_network(..., cluster_by = ...)."),
      class = "nestimate_bad_input", call = NULL))
  }
  build_mmm(data = data, k = k, n_starts = n_starts,
            max_iter = max_iter, tol = tol, smooth = smooth,
            seed = seed, covariates = covariates,
            covariate_effect = covariate_effect,
            estimator = estimator)
}

#' @rdname build_mmm
#' @export
print.net_mmm_clustering <- function(x, digits = 3L, ...) {
  .net_mmm_check_unused_dots("print.net_mmm_clustering", ...)
  .net_mmm_check_digits(digits)
  digits <- as.integer(digits)
  k <- as.integer(x$k)
  n_total <- as.integer(x$n_sequences)

  cat("MMM Clustering [k = ", k, "]\n", sep = "")
  cat(sprintf("  Sequences: %d  |  Clusters: %d\n", n_total, k))
  if (!is.null(x$quality)) {
    cat(sprintf(
      "  Quality: AvePP = %.*f  |  Entropy = %.*f  |  Class.Err = %.1f%%\n",
      digits, x$quality$avepp_overall, digits, x$quality$entropy,
      x$quality$classification_error * 100))
  }
  if (!is.null(x$BIC)) {
    cat(sprintf("  ICs: BIC = %.*f  |  AIC = %.*f  |  ICL = %.*f\n",
                digits, x$BIC, digits, x$AIC, digits, x$ICL))
  }

  cat("\n")
  sizes <- as.integer(tabulate(x$assignments, nbins = k))
  cols <- list(
    Cluster = sprintf("%d", seq_len(k)),
    N       = .fmt_size_pct(sizes, n_total),
    `Mix%`  = sprintf("%4.1f%%", as.numeric(x$mixing) * 100),
    AvePP   = sprintf(paste0("%.", digits, "f"),
                      as.numeric(x$quality$avepp))
  )
  cat(paste(.cluster_table_lines(cols), collapse = "\n"), "\n", sep = "")

  if (!is.null(x$covariates)) {
    cov_names <- setdiff(
      unique(x$covariates$coefficients$variable), "(Intercept)"
    )
    cat(sprintf("\n  Covariates: %s (%s, %d predictors)\n",
                paste(cov_names, collapse = ", "),
                .mmm_covariate_label(x$covariate_effect), length(cov_names)))
  }

  invisible(x)
}

#' @rdname build_mmm
#' @section Methods:
#' `print()` of a `net_mmm` shows a header with the numbers of sequences,
#' subgroups and states, the log-likelihood and the information criteria; a
#' quality line with the overall AvePP, the normalised entropy (0 when every
#' sequence is classified with certainty, 1 when all posteriors are equal)
#' and the classification error (the expected share of sequences
#' misclassified when each is assigned to its most likely subgroup, the mean
#' of one minus the largest posterior; Vermunt & Magidson, 2005); a line when the EM did not converge;
#' and one row per subgroup with its size, mixing proportion and AvePP. For a
#' model fitted with `covariates`, a last line names them, says whether they
#' entered the EM (integrated) or were related post hoc, and gives the number
#' of predictors. The layout follows `print()` of a `net_clustering` from
#' [build_clusters()], so dissimilarity-based and model-based clusterings read
#' alike.
#'
#' [build_network()] on a fit, and [cluster_network()] with
#' `cluster_by = "mmm"`, return a `netobject_group` whose `"clustering"`
#' attribute is a `net_mmm_clustering`: the fit without its
#' subgroup networks. Its `print()` shows the same quality line,
#' information criteria and subgroup table, without the state count, the
#' log-likelihood or the convergence line, and its `plot()` draws the same
#' two views as `plot()` of a `net_mmm`. The `"silhouette"`, `"mds"` and
#' `"heatmap"` views of a dissimilarity-based clustering need a distance
#' matrix, which an MMM clustering does not have, so asking for them raises
#' an error.
#' @export
plot.net_mmm_clustering <- function(x, type = c("posterior", "covariates",
                                                 "predictors"),
                                     combined = TRUE, ...) {
  .net_mmm_check_unused_dots("plot.net_mmm_clustering", ...)
  type <- match.arg(type)
  .assert(is.logical(combined), length(combined) == 1L, !is.na(combined))
  if (type == "predictors") type <- "covariates"

  if (type == "covariates") {
    if (is.null(x$covariates)) {
      stop(errorCondition(
        paste0(
          "No covariate analysis on this clustering. Re-run ",
          "cluster_mmm() (or build_mmm()) with covariates = ..."),
        class = "nestimate_bad_input", call = NULL))
    }
    return(.plot_covariate_forest(
      x$covariates$coefficients,
      sprintf("Covariate Effects (ref: Cluster %s)",
              x$covariates$fit$reference_cluster),
      combined = combined
    ))
  }

  # type == "posterior": .plot_mmm_posterior() reads only $posterior +
  # $assignments, both carried by net_mmm_clustering.
  .plot_mmm_posterior(x)
}

# The rows of an association network (cor, pcor, glasso, ising, mgm) are
# variables, so it has no transitions to fit a mixture of chains on.
.mmm_refuse_association <- function(net) {
  if (isTRUE(net$method %in% .association_methods)) {
    stop(errorCondition(sprintf(paste(
      "build_mmm(): this \"%s\" network holds variables, not sequences of",
      "states. Build it with a transition method or pass the sequences."),
      net$method), class = "nestimate_bad_input", call = NULL))
  }
  invisible(NULL)
}


#' Counts of a network converted from a tna model, when they can be trusted
#'
#' `cograph::as_cograph()` keeps a tna model's weights and encoded sequences,
#' but not its type, windows or settings. The sequences are read as plain
#' positions only when recounting them reproduces the weights, as counts or as
#' row probabilities; otherwise the conversion lost something and NULL is
#' returned for the caller to refuse.
#' @noRd
.tna_conversion_counts <- function(x) {
  w <- x$weights
  if (!is.matrix(w) || is.null(x$data) || is.null(rownames(w))) return(NULL)
  data <- .as_netobject(x)$data
  if (!is.data.frame(data) || nrow(data) == 0L) return(NULL)
  counts <- build_network(data, method = "frequency",
                          state_cols = names(data))$weights
  states <- rownames(w)
  if (!all(rownames(counts) %in% states)) return(NULL)
  full <- matrix(0, length(states), length(states),
                 dimnames = list(states, states))
  full[rownames(counts), colnames(counts)] <- counts
  sums <- rowSums(full)
  probs <- full / ifelse(sums > 0, sums, 1)
  same <- function(a) isTRUE(all.equal(unname(a), unname(w), tolerance = 1e-10))
  if (same(full) || same(probs)) full else NULL
}
