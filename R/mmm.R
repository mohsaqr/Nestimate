# Mixed Markov Model — latent-class mixture of Markov chains
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
#' @param from_ind Pre-computed K^2 x K indicator grouping by "from" state
#' @return List with P_all, init_all, pi_mix, posterior, ll, iterations, converged
#' @noRd
.mmm_em <- function(counts, init_state, n_comp, max_iter, tol, smooth, K,
                    init_posterior = NULL, from_ind = NULL, cov_df = NULL) {
  N <- nrow(counts)
  K2 <- K * K

  # Initialize posterior (N x M)
  if (is.null(init_posterior)) {
    post <- matrix(runif(N * n_comp), nrow = N, ncol = n_comp)
    post <- post / rowSums(post)
  } else {
    post <- init_posterior
  }

  # Pre-compute from-state indicator: K^2 x K binary matrix
  # from_ind[j, k] = 1 if pair j has from-state = k
  # j = (from-1)*K + to, so from = (j-1) %/% K + 1
  if (is.null(from_ind)) {
    from_idx <- rep(seq_len(K), each = K)
    from_ind <- matrix(0, K2, K)
    from_ind[cbind(seq_len(K2), from_idx)] <- 1
  }
  from_ind_t <- t(from_ind)  # K x K^2

  # Initial state indicator: N x K binary matrix
  init_ind <- matrix(0, N, K)
  valid_init <- !is.na(init_state)
  init_ind[cbind(which(valid_init), init_state[valid_init])] <- 1

  ll_prev <- -Inf
  converged <- FALSE
  cov_fit <- NULL
  cov_beta <- NULL

  # Pre-compute design matrix for covariates (once)
  X_cov <- NULL
  if (!is.null(cov_df)) {
    X_cov <- stats::model.matrix(~ ., data = cov_df)
  }

  for (iter in seq_len(max_iter)) {
    # ---- M-step ----
    # Transition probabilities: weighted counts grouped by from-state
    weighted_smooth <- crossprod(post, counts) + smooth  # M x K^2
    from_sums <- weighted_smooth %*% from_ind  # M x K
    from_sums[from_sums == 0] <- 1
    divisor <- from_sums %*% from_ind_t  # M x K^2
    P_all <- t(weighted_smooth / divisor)  # K^2 x M

    # Initial state probabilities: weighted first-state counts
    init_weighted <- crossprod(post, init_ind) + smooth  # M x K
    init_sums <- rowSums(init_weighted)
    init_sums[init_sums == 0] <- 1
    init_all <- t(init_weighted / init_sums)  # K x M

    # Mixing proportions
    if (is.null(cov_df)) {
      pi_mix <- .colMeans(post, N, n_comp)
      log_pi <- log(pi_mix + 1e-300)
    } else {
      sm <- .mmm_softmax_mstep_fast(post, X_cov, beta_prev = cov_beta,
                                      n_steps = 3L)
      log_pi_mat <- sm$log_pi_mat
      cov_beta <- sm$beta
      pi_mix <- colMeans(exp(log_pi_mat))
    }

    # ---- E-step ----
    # Transition log-likelihood: N x M
    log_lik <- counts %*% log(P_all + 1e-300)
    # Initial state log-likelihood: N x M
    log_lik <- log_lik + init_ind %*% log(init_all + 1e-300)
    # Mixing proportions
    if (is.null(cov_df)) {
      log_lik <- log_lik + rep(log_pi, each = N)
    } else {
      log_lik <- log_lik + log_pi_mat
    }

    # Log-sum-exp with vectorized row max
    if (n_comp == 2L) {
      log_max <- pmax(log_lik[, 1L], log_lik[, 2L])
    } else {
      log_max <- log_lik[, 1L]
      for (m in 2:n_comp) log_max <- pmax(log_max, log_lik[, m])
    }
    log_lik <- log_lik - log_max
    post <- exp(log_lik)
    row_sums <- rowSums(post)
    row_sums[row_sums == 0] <- 1e-300
    post <- post / row_sums

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

  # Max posterior per sequence (computed once, reused)
  max_post <- do.call(pmax, as.data.frame(posterior))

  # Overall AvePP (mean of max posteriors)
  avepp_overall <- mean(max_post)

  # Entropy (normalized: 0 = perfect, 1 = random)
  ent_raw <- -sum(posterior * log(posterior + 1e-300))
  ent_max <- N * log(M)
  entropy <- if (ent_max > 0) ent_raw / ent_max else 0
  relative_entropy <- 1 - entropy

  # Classification error (proportion with max posterior < 0.5)
  classification_error <- mean(max_post < 0.5)

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

# ---------------------------------------------------------------------------
# Covariate M-step helpers
# ---------------------------------------------------------------------------

#' Fast softmax regression M-step via Newton-Raphson
#'
#' Replaces nnet::multinom() inside the EM loop. Does a few Newton steps
#' on the pre-computed design matrix. Much faster than calling multinom()
#' at every EM iteration.
#'
#' @param post N x k posterior matrix
#' @param X N x p design matrix (with intercept column)
#' @param beta_prev (k-1) x p coefficient matrix from previous iteration (or NULL)
#' @param n_steps Number of Newton-Raphson steps (default: 5)
#' @return List with log_pi_mat (N x k) and beta ((k-1) x p)
#' @noRd
.mmm_softmax_mstep_fast <- function(post, X, beta_prev = NULL, n_steps = 5L) {
  N <- nrow(X)
  p <- ncol(X)
  k <- ncol(post)
  km1 <- k - 1L

  # Initialize beta
  if (is.null(beta_prev)) {
    beta <- matrix(0, km1, p)
  } else {
    beta <- beta_prev
  }

  for (step in seq_len(n_steps)) {
    # Compute softmax probabilities: N x k
    eta <- X %*% t(beta)  # N x (k-1)
    eta_full <- cbind(0, eta)  # N x k (reference class = 0)
    eta_max <- do.call(pmax, as.data.frame(eta_full))
    exp_eta <- exp(eta_full - eta_max)
    pi_mat <- exp_eta / rowSums(exp_eta)

    # Gradient for each non-reference class
    residual <- post[, -1L, drop = FALSE] - pi_mat[, -1L, drop = FALSE]  # N x (k-1)

    if (km1 == 1L) {
      # k=2: logistic regression Newton step
      grad <- crossprod(X, residual)  # p x 1
      W <- pi_mat[, 1L] * pi_mat[, 2L]  # N-length
      H <- crossprod(X * W, X) + diag(1e-6, p)
      H_inv <- tryCatch(solve(H), error = function(e) NULL)
      if (is.null(H_inv)) break # nocov
      beta <- beta + t(H_inv %*% grad)
    } else {
      # General k: block Newton step
      beta_vec <- as.numeric(t(beta))
      grad_vec <- as.numeric(t(crossprod(X, residual)))  # km1*p

      delta <- numeric(km1 * p)
      for (m in seq_len(km1)) {
        W <- pi_mat[, m + 1L] * (1 - pi_mat[, m + 1L])
        H <- crossprod(X * W, X) + diag(1e-6, p)
        H_inv <- tryCatch(solve(H), error = function(e) NULL)
        if (is.null(H_inv)) next # nocov
        idx <- ((m - 1L) * p + 1L):(m * p)
        delta[idx] <- H_inv %*% grad_vec[idx]
      }
      beta <- matrix(beta_vec + delta, km1, p, byrow = TRUE)
    }
  }

  log_pi_mat <- log(pi_mat + 1e-300)
  list(log_pi_mat = log_pi_mat, beta = beta)
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
#' Discovers latent subgroups with different transition dynamics using
#' Expectation-Maximization. Each mixture component has its own transition
#' matrix. Sequences are probabilistically assigned to components.
#'
#' @param data A data.frame (wide format), \code{netobject}, or
#'   \code{tna} model. For tna objects, extracts the stored data.
#' @param k Integer. Number of mixture components. Default: 2.
#' @param n_starts Integer. Number of random restarts. Default: 50.
#' @param max_iter Integer. Maximum EM iterations per start. Default: 200.
#' @param tol Numeric. Convergence tolerance. Default: 1e-6.
#' @param smooth Numeric. Laplace smoothing constant. Default: 0.01.
#' @param seed Integer or NULL. Random seed.
#' @param covariates Optional. Covariates integrated into the EM algorithm
#'   to model covariate-dependent mixing proportions. Accepts formula,
#'   character vector, string, or data.frame (same forms as
#'   \code{\link{cluster_data}}). Unlike the post-hoc analysis in
#'   \code{cluster_data()}, these covariates directly influence cluster
#'   membership during estimation. Requires the \pkg{nnet} package.
#'
#' @return An object of class \code{net_mmm} with components:
#'   \describe{
#'     \item{models}{List of \code{netobject}s, one per component.}
#'     \item{k}{Number of components.}
#'     \item{mixing}{Numeric vector of mixing proportions.}
#'     \item{posterior}{N x k matrix of posterior probabilities.}
#'     \item{assignments}{Integer vector of hard assignments (1..k).}
#'     \item{quality}{List: \code{avepp} (per-class), \code{avepp_overall},
#'       \code{entropy}, \code{relative_entropy},
#'       \code{classification_error}.}
#'     \item{log_likelihood, BIC, AIC, ICL}{Model fit statistics.}
#'     \item{states}{Character vector of state names.}
#'   }
#'
#' @examples
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
                      covariates = NULL) {

  k <- as.integer(k)
  n_starts <- as.integer(n_starts)
  max_iter <- as.integer(max_iter)
  stopifnot(
    "'k' must be >= 2" = k >= 2L,
    "'n_starts' must be >= 1" = n_starts >= 1L,
    "'max_iter' must be >= 1" = max_iter >= 1L,
    "'tol' must be > 0" = tol > 0,
    "'smooth' must be >= 0" = smooth >= 0
  )

  # ---- Extract data and states ----
  if (inherits(data, "tna") || inherits(data, "ftna")) {
    raw_data <- data$data
    states <- data$labels
    if (is.null(raw_data)) {
      stop("tna model does not contain $data. ", # nocov start
           "Rebuild with tna::tna(data) where data is kept.", call. = FALSE) # nocov end
    }
    raw_data <- .mmm_decode_int_data(raw_data, states)
  } else if (inherits(data, "cograph_network") && !inherits(data, "netobject")) {
    # Pure cograph_network (from cograph package): coerce to decode integer data
    data <- .as_netobject(data)
    raw_data <- data$data
    states <- data$nodes$label
  } else if (inherits(data, "netobject")) {
    raw_data <- data$data
    states <- data$nodes$label
  } else {
    raw_data <- as.data.frame(data)
    state_cols <- .select_state_cols(raw_data, id = NULL, cols = NULL)
    mat <- as.matrix(raw_data[, state_cols, drop = FALSE])
    states <- sort(unique(as.vector(mat)))
    states <- states[!is.na(states)]
  }

  n_states <- length(states)
  stopifnot("Need at least 2 states" = n_states >= 2L)

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
    if (length(keep) > 0L) raw_data <- raw_data[, keep, drop = FALSE]
    # Re-extract states after removing covariate columns
    state_cols_clean <- .select_state_cols(raw_data, id = NULL, cols = NULL)
    mat <- as.matrix(raw_data[, state_cols_clean, drop = FALSE])
    states <- sort(unique(as.vector(mat)))
    states <- states[!is.na(states)]
    n_states <- length(states)
    stopifnot("Need at least 2 states" = n_states >= 2L)
  }

  # ---- Precompute per-sequence transition counts ----
  counts <- .precompute_per_sequence_wide(
    raw_data, method = "relative", cols = NULL, id_col = NULL, states = states
  )
  N <- nrow(counts)
  K2 <- n_states * n_states

  # Extract initial states (first non-NA state per sequence)
  state_cols <- .select_state_cols(raw_data, id = NULL, cols = NULL)
  first_col <- as.character(raw_data[[state_cols[1L]]])
  init_state <- match(first_col, states)

  # Handle covariate NAs: subset counts, init_state, cov_df
  if (!is.null(cov_df)) {
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
      N <- nrow(counts)
    }
  }

  stopifnot("Need more sequences than components" = N > k)

  # ---- Run EM with screen-then-refine strategy ----
  if (!is.null(seed)) set.seed(seed)

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
    n_cores <- min(parallel::detectCores(logical = FALSE) %||% 1L, n_starts)
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
    if (screen_results[[idx]]$converged && is.null(cov_df)) {
      run <- screen_results[[idx]]
    } else {
      run <- .mmm_em(counts, init_state, k, max_iter, tol, smooth, n_states,
                     init_posterior = screen_results[[idx]]$posterior,
                     from_ind = from_ind, cov_df = cov_df)
    }
    if (run$ll > best$ll) best <- run
  }

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

    structure(list(
      data = raw_data,
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

  names(models) <- paste0("Component_", seq_len(k))

  # ---- Assignments & quality ----
  assignments <- apply(best$posterior, 1L, which.max)
  quality <- .mmm_quality(best$posterior, assignments, k)

  # ---- Information criteria ----
  n_params <- k * n_states * (n_states - 1L) +
              k * (n_states - 1L) +
              if (is.null(cov_df)) (k - 1L) else (k - 1L) * (ncol(cov_df) + 1L)
  BIC_val <- -2 * best$ll + n_params * log(N)
  AIC_val <- -2 * best$ll + 2 * n_params
  ICL_val <- BIC_val + 2 * quality$class_entropy

  # ---- Covariate output ----
  cov_result <- NULL
  if (!is.null(cov_df) && !is.null(best$cov_beta)) {
    # Run final nnet::multinom once for SEs and proper inference
    cov_result <- .run_covariate_analysis(
      assignments, cov_df, paste(names(cov_df), collapse = " + "), k
    )
    # Store EM-estimated beta
    cov_result$beta <- best$cov_beta
  }

  structure(list(
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
    covariates = cov_result
  ), class = "net_mmm")
}

# ---------------------------------------------------------------------------
# Compare multiple k values
# ---------------------------------------------------------------------------

#' Compare MMM fits across different k
#'
#' @param data Data frame, netobject, or tna model.
#' @param k Integer vector of component counts. Default: 2:5.
#' @param ... Arguments passed to \code{\link{build_mmm}}.
#'
#' @return A \code{mmm_compare} data frame with BIC, AIC, ICL, AvePP,
#'   entropy per k.
#'
#' @examples
#' \donttest{
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
#'   V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
#' )
#' comp <- compare_mmm(seqs, k = 2:3, seed = 42)
#' print(comp)
#' }
#'
#' @export
compare_mmm <- function(data, k = 2:5, ...) {
  results <- lapply(k, function(m) {
    fit <- build_mmm(data, k = m, ...)
    data.frame(
      k = m,
      log_likelihood = fit$log_likelihood,
      AIC = fit$AIC,
      BIC = fit$BIC,
      ICL = fit$ICL,
      AvePP = fit$quality$avepp_overall,
      Entropy = fit$quality$entropy,
      converged = fit$converged,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, results)
  class(result) <- c("mmm_compare", "data.frame")
  result
}


# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print Method for net_mmm
#'
#' @param x A \code{net_mmm} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
print.net_mmm <- function(x, ...) {
  cat("Mixed Markov Model\n")
  cat(sprintf("  k = %d | %d sequences | %d states\n",
              x$k, x$n_sequences, length(x$states)))
  cat(sprintf("  LL = %.1f | BIC = %.1f | ICL = %.1f\n",
              x$log_likelihood, x$BIC, x$ICL))

  # Cluster table
  cat("\n  Cluster  Size  Mix%%   AvePP\n")
  cat("  " , strrep("-", 30), "\n", sep = "")
  for (m in seq_len(x$k)) {
    n_in <- sum(x$assignments == m)
    cat(sprintf("  %7d  %4d  %4.1f%%  %.3f\n",
                m, n_in, x$mixing[m] * 100, x$quality$avepp[m]))
  }
  cat(sprintf("\n  Overall AvePP = %.3f | Entropy = %.3f | Class.Err = %.1f%%\n",
              x$quality$avepp_overall, x$quality$entropy,
              x$quality$classification_error * 100))
  if (!is.null(x$covariates)) {
    cov_names <- setdiff(
      unique(x$covariates$coefficients$variable), "(Intercept)"
    )
    cat(sprintf("  Covariates:    %s (integrated, %d predictors)\n",
                paste(cov_names, collapse = ", "), length(cov_names)))
  }
  invisible(x)
}

#' Summary Method for net_mmm
#'
#' @param object A \code{net_mmm} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
summary.net_mmm <- function(object, ...) {
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
      cov, "Covariate Analysis (integrated into EM -- influences cluster membership)"
    )
  }

  invisible(object)
}

#' Plot Method for net_mmm
#'
#' @param x A \code{net_mmm} object.
#' @param type Character. Plot type: \code{"posterior"} (default) or \code{"covariates"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{ggplot} object, invisibly.
#'
#' @export
plot.net_mmm <- function(x, type = c("posterior", "covariates"), ...) {
  type <- match.arg(type)

  if (type == "covariates") {
    if (is.null(x$covariates)) {
      stop("No covariate analysis found. Run build_mmm() with covariates.",
           call. = FALSE)
    }
    return(.plot_covariate_forest(
      x$covariates$coefficients,
      sprintf("Covariate Effects (ref: Cluster %s)",
              x$covariates$fit$reference_cluster)
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
    ggplot2::theme_minimal()

  print(p)
  invisible(p)
}

#' Print Method for mmm_compare
#'
#' @param x An \code{mmm_compare} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
print.mmm_compare <- function(x, ...) {
  cat("MMM Model Comparison\n\n")
  best_bic <- which.min(x$BIC)
  best_icl <- which.min(x$ICL)
  x$best <- ""
  x$best[best_bic] <- "<-- BIC"
  if (best_icl != best_bic) x$best[best_icl] <- paste(x$best[best_icl], "<-- ICL")
  print.data.frame(x, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' Plot Method for mmm_compare
#'
#' @param x An \code{mmm_compare} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{ggplot} object, invisibly.
#'
#' @export
plot.mmm_compare <- function(x, ...) {
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
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  print(p)
  invisible(p)
}

# ---------------------------------------------------------------------------
# cluster_mmm — convenience alias for build_mmm
# ---------------------------------------------------------------------------

#' Cluster sequences using Mixed Markov Models
#'
#' Convenience alias for \code{\link{build_mmm}}. Fits a mixture of Markov
#' chains to sequence data and returns per-component transition networks with
#' EM-fitted initial state probabilities.
#'
#' Use \code{\link{build_network}} on the result to extract per-cluster
#' networks with any estimation method, or use \code{\link{cluster_network}}
#' for a one-shot clustering + network call.
#'
#' @inheritParams build_mmm
#' @return A \code{net_mmm} object. See \code{\link{build_mmm}} for details.
#' @seealso \code{\link{build_mmm}}, \code{\link{cluster_network}}
#' @export
cluster_mmm <- build_mmm
