# ---- Sequence Clustering ----

# Available metrics and methods (private constants)
.clustering_metrics <- c(
  "hamming", "osa", "lv", "dl", "lcs", "qgram", "cosine", "jaccard", "jw"
)
.clustering_methods <- c(
  "pam", "ward.D2", "ward.D", "complete", "average",
  "single", "mcquitty", "median", "centroid"
)

# ==============================================================================
# 1. Sequence Encoding
# ==============================================================================

#' Encode sequences as integer matrix
#'
#' @param data Data frame or matrix of sequences (rows = sequences, cols = time)
#' @param na_syms Character vector of symbols to treat as NA
#' @return List with int_mat, len, states, n_states
#' @noRd
.encode_sequences <- function(data, na_syms = c("*", "%")) {
  mat <- as.matrix(data)
  storage.mode(mat) <- "character"
  mat[mat %in% na_syms] <- NA
  obs <- !is.na(mat)
  last_obs <- max.col(obs, ties.method = "last")
  all_na <- rowSums(obs) == 0L
  last_obs[all_na] <- 0L
  states <- sort(unique(mat[obs]))
  n_states <- length(states)
  int_mat <- matrix(match(mat, states), nrow = nrow(mat), ncol = ncol(mat))
  list(int_mat = int_mat, len = last_obs, states = states, n_states = n_states)
}

# ==============================================================================
# 2. Distance Functions (R fallback)
# ==============================================================================
# The loops in this section are deliberate. The edit distances (Levenshtein,
# OSA, Damerau-Levenshtein, LCS) are dynamic programmes: every cell of the
# table is read from cells filled just before it. Jaro matching is greedy and
# sequential (a matched position cannot be matched again), and the
# Jaro-Winkler prefix stops at the first mismatch.

#' @noRd
.hamming_dist_r <- function(x, y, n, m, weights, ...) {
  sum((x != y) * weights)
}

#' @noRd
.levenshtein_dist_r <- function(x, y, n, m, ...) {
  d <- matrix(0L, n + 1L, m + 1L)
  d[, 1L] <- 0:n
  d[1L, ] <- 0:m
  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      cost <- as.integer(x[i] != y[j])
      d[i + 1L, j + 1L] <- min(
        d[i, j + 1L] + 1L,
        d[i + 1L, j] + 1L,
        d[i, j] + cost
      )
    }
  }
  d[n + 1L, m + 1L]
}

#' @noRd
.osa_dist_r <- function(x, y, n, m, ...) {

  d <- matrix(0L, n + 1L, m + 1L)
  d[, 1L] <- 0:n
  d[1L, ] <- 0:m
  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      cost <- as.integer(x[i] != y[j])
      d[i + 1L, j + 1L] <- min(
        d[i, j + 1L] + 1L,
        d[i + 1L, j] + 1L,
        d[i, j] + cost
      )
      if (i > 1L && j > 1L && x[i] == y[j - 1L] && x[i - 1L] == y[j]) {
        d[i + 1L, j + 1L] <- min(d[i + 1L, j + 1L], d[i - 1L, j - 1L] + 1L)
      }
    }
  }
  d[n + 1L, m + 1L]
}

#' @noRd
.dl_dist_r <- function(x, y, n, m, ...) {
  if (n == 0L) return(m)
  if (m == 0L) return(n)
  x_sub <- x[seq_len(n)]
  y_sub <- y[seq_len(m)]
  alphabet <- unique(c(x_sub, y_sub))
  last <- stats::setNames(rep(0L, length(alphabet)), as.character(alphabet))
  d <- matrix(0L, n + 2L, m + 2L)
  maxdist <- n + m
  d[1L, ] <- maxdist
  d[, 1L] <- maxdist
  d[seq_len(n + 1L) + 1L, 2L] <- seq_len(n + 1L) - 1L
  d[2L, seq_len(m + 1L) + 1L] <- seq_len(m + 1L) - 1L
  for (i in seq_len(n)) {
    db <- 0L
    for (j in seq_len(m)) {
      k <- last[[as.character(y[j])]]
      if (is.null(k)) k <- 0L # nocov
      l <- db
      if (x[i] == y[j]) {
        cost <- 0L
        db <- j
      } else {
        cost <- 1L
      }
      d[i + 2L, j + 2L] <- min(
        d[i + 1L, j + 2L] + 1L,
        d[i + 2L, j + 1L] + 1L,
        d[i + 1L, j + 1L] + cost,
        d[k + 1L, l + 1L] + (i - k - 1L) + (j - l - 1L) + 1L
      )
    }
    last[[as.character(x[i])]] <- i
  }
  d[n + 2L, m + 2L]
}

#' @noRd
.lcs_dist_r <- function(x, y, n, m, ...) {
  lcs_len <- matrix(0L, n + 1L, m + 1L)
  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      if (x[i] == y[j]) {
        lcs_len[i + 1L, j + 1L] <- lcs_len[i, j] + 1L
      } else {
        lcs_len[i + 1L, j + 1L] <- max(lcs_len[i, j + 1L], lcs_len[i + 1L, j])
      }
    }
  }
  n + m - 2L * lcs_len[n + 1L, m + 1L]
}

#' Extract q-gram frequency profile from integer sequence
#' @noRd
.get_qgram_r <- function(x, n, q) {
  if (n < q) return(integer()) # nocov start
  ng <- n - q + 1L
  grams <- vapply(seq_len(ng), function(i) {
    paste0(x[i:(i + q - 1L)], collapse = "\x01")
  }, character(1L))
  tab <- table(grams)
  counts <- as.integer(tab)
  names(counts) <- names(tab)
  counts # nocov end
}

#' @noRd
.qgram_dist_r <- function(x, y, qx, qy, ...) {
  grams <- union(names(qx), names(qy))
  x_all <- qx[grams]
  y_all <- qy[grams]
  x_all[is.na(x_all)] <- 0L
  y_all[is.na(y_all)] <- 0L
  sum(abs(x_all - y_all))
}

#' @noRd
.cosine_dist_r <- function(x, y, qx, qy, ...) {
  grams <- union(names(qx), names(qy))
  x_all <- qx[grams]
  y_all <- qy[grams]
  x_all[is.na(x_all)] <- 0L
  y_all[is.na(y_all)] <- 0L
  num <- sum(x_all * y_all)
  den <- sqrt(sum(x_all^2)) * sqrt(sum(y_all^2))
  if (den == 0) return(1)
  1 - num / den
}

#' @noRd
.jaccard_dist_r <- function(qx, qy, ...) {
  qx_names <- names(qx)
  qy_names <- names(qy)
  if (length(qx_names) == 0L && length(qy_names) == 0L) return(0)
  1 - length(intersect(qx_names, qy_names)) / length(union(qx_names, qy_names))
}

#' @noRd
.jaro_dist_r <- function(x, y, n, m, ...) {
  if (n == 0L && m == 0L) return(0) # nocov start
  if (n == 0L || m == 0L) return(1) # nocov end
  match_dist <- max(floor(max(n, m) / 2) - 1L, 0L)
  x_match <- rep(FALSE, n)
  y_match <- rep(FALSE, m)
  matches <- 0L
  for (i in seq_len(n)) {
    start <- max(1L, i - match_dist)
    end <- min(m, i + match_dist)
    if (start > end) next # nocov
    for (j in start:end) {
      if (!y_match[j] && x[i] == y[j]) {
        x_match[i] <- TRUE
        y_match[j] <- TRUE
        matches <- matches + 1L
        break
      }
    }
  }
  if (matches == 0L) return(1) # nocov
  trans <- 0L
  k <- 1L
  for (i in seq_len(n)) {
    if (x_match[i]) {
      while (k <= m && !y_match[k]) k <- k + 1L
      if (k <= m && x[i] != y[k]) trans <- trans + 1L
      k <- k + 1L
    }
  }
  trans <- trans / 2L
  jaro_sim <- (matches / n + matches / m + (matches - trans) / matches) / 3
  1 - jaro_sim
}

#' @noRd
.jw_dist_r <- function(x, y, n, m, p = 0.1, max_l = 4L, ...) {
  jaro_sim <- 1 - .jaro_dist_r(x, y, n, m)
  if (jaro_sim == 0) return(1) # nocov
  pref <- min(n, m, max_l)
  l <- 0L
  for (k in seq_len(pref)) {
    if (x[k] == y[k]) {
      l <- l + 1L
    } else {
      break
    }
  }
  p <- min(p, 0.25)
  jw <- jaro_sim + l * p * (1 - jaro_sim)
  1 - jw
}

# ==============================================================================
# 3. Vectorized Distance Matrix (hamming, qgram-family)
# ==============================================================================

#' Hamming distance matrix via column-wise outer() - no pair loop
#' @noRd
.hamming_matrix_outer <- function(mat_s, weights) {
  n <- nrow(mat_s)
  k <- ncol(mat_s)
  d <- matrix(0, n, n)
  # Loops on purpose: one n x n mismatch matrix is added at a time, instead of
  # holding an n x n x k array.
  if (length(weights) == 1L) {
    # Unweighted: just count mismatches per column
    for (c in seq_len(k)) {
      d <- d + outer(mat_s[, c], mat_s[, c], "!=")
    }
  } else {
    for (c in seq_len(k)) {
      d <- d + outer(mat_s[, c], mat_s[, c], "!=") * weights[c]
    }
  }
  stats::as.dist(d)
}

#' Build term-document matrix of q-gram counts (rows = seqs, cols = qgrams)
#' @noRd
.build_qgram_tdm <- function(mat_s, len, q) {
  n <- nrow(mat_s)
  # Collect all q-gram profiles
  profiles <- lapply(seq_len(n), function(i) {
    if (len[i] < q) return(NULL)
    grams <- vapply(seq_len(len[i] - q + 1L), function(j) {
      paste0(mat_s[i, j:(j + q - 1L)], collapse = "\x01")
    }, character(1L))
    table(grams)
  })
  gram_names <- as.character(unlist(lapply(profiles, names)))
  # Columns in order of first appearance.
  all_grams <- unique(gram_names)
  # Build matrix: rows = sequences, cols = unique qgrams
  tdm <- matrix(0L, nrow = n, ncol = length(all_grams))
  colnames(tdm) <- all_grams
  tdm[cbind(rep(seq_len(n), lengths(profiles)),
            match(gram_names, all_grams))] <-
    as.integer(unlist(lapply(profiles, as.integer)))
  tdm
}

#' Q-gram distance matrix via term-document matrix subtraction
#' @noRd
.qgram_matrix_tdm <- function(tdm) {
  n <- nrow(tdm)
  # |a_i - b_i| summed over all qgrams = Manhattan distance between rows
  # Expand via outer difference per column, then sum absolute values
  d <- matrix(0, n, n)
  # A loop on purpose: one n x n matrix per q-gram, not an n x n x grams array.
  for (c in seq_len(ncol(tdm))) {
    d <- d + abs(outer(tdm[, c], tdm[, c], "-"))
  }
  stats::as.dist(d)
}

#' Cosine distance matrix via matrix multiply
#' @noRd
.cosine_matrix_tdm <- function(tdm) {
  # cosine_sim = (A %*% t(A)) / (norm_a * norm_b)
  norms <- sqrt(rowSums(tdm^2))
  norms[norms == 0] <- 1  # avoid division by zero (returns dist=1 for zero vecs)
  dot <- tcrossprod(tdm)
  sim <- dot / outer(norms, norms)
  # Clamp to [0, 1] for numerical safety
  sim[sim > 1] <- 1
  sim[sim < 0] <- 0
  # Handle zero-norm rows: distance = 1
  zero_rows <- which(rowSums(tdm^2) == 0)
  d <- 1 - sim
  diag(d) <- 0
  if (length(zero_rows) > 0L) {
    d[zero_rows, ] <- 1
    d[, zero_rows] <- 1
    d[cbind(zero_rows, zero_rows)] <- 0
  }
  stats::as.dist(d)
}

#' Jaccard distance matrix via set intersection/union on qgram presence
#' @noRd
.jaccard_matrix_tdm <- function(tdm) {
  # Binary presence matrix
  bin <- (tdm > 0L) * 1L
  # intersection = A %*% t(B), union = |A| + |B| - intersection
  isect <- tcrossprod(bin)
  sizes <- rowSums(bin)
  union_mat <- outer(sizes, sizes, "+") - isect
  # Avoid 0/0 (both empty -> distance 0)
  d <- ifelse(union_mat == 0, 0, 1 - isect / union_mat)
  diag(d) <- 0
  stats::as.dist(d)
}

# ==============================================================================
# 4. Pairwise Distance Matrix (DP distances - pair loop unavoidable)
# ==============================================================================

#' Compute pairwise dissimilarity matrix using pure R
#' @noRd
.dissimilarity_matrix_r <- function(enc, dissimilarity, lambda,
                                    q = 2L, p = 0.1) {
  int_mat <- enc$int_mat
  len <- enc$len
  n <- nrow(int_mat)
  k <- ncol(int_mat)
  sentinel <- enc$n_states + 1L

  # Full-width matrix with sentinel for NA
  mat_s <- int_mat
  mat_s[is.na(mat_s)] <- sentinel

  # --- Vectorized fast paths ---
  if (dissimilarity == "hamming") {
    weights <- 1
    if (lambda > 0) {
      weights <- exp(-lambda * seq(0, k - 1))
      weights <- weights / max(weights)
    }
    return(.hamming_matrix_outer(mat_s, weights))
  }

  if (dissimilarity %in% c("qgram", "cosine", "jaccard")) {
    tdm <- .build_qgram_tdm(mat_s, len, q)
    if (dissimilarity == "qgram") return(.qgram_matrix_tdm(tdm))
    if (dissimilarity == "cosine") return(.cosine_matrix_tdm(tdm))
    return(.jaccard_matrix_tdm(tdm))
  }

  # --- Pair loop for DP distances (osa, lv, dl, lcs, jw) ---
  d_fun <- switch(dissimilarity,
    osa = .osa_dist_r,
    lv = .levenshtein_dist_r,
    dl = .dl_dist_r,
    lcs = .lcs_dist_r,
    jw = .jw_dist_r
  )

  d <- matrix(0, n, n)
  pairs <- combn(n, 2L)

  dists <- vapply(seq_len(ncol(pairs)), function(idx) {
    i <- pairs[1L, idx]
    j <- pairs[2L, idx]
    if (dissimilarity == "jw") {
      d_fun(x = mat_s[i, ], y = mat_s[j, ],
            n = len[i], m = len[j], p = p)
    } else {
      d_fun(x = mat_s[i, ], y = mat_s[j, ],
            n = len[i], m = len[j])
    }
  }, numeric(1L))

  d[lower.tri(d)] <- dists
  d <- d + t(d)
  stats::as.dist(d)
}

# ==============================================================================
# 4. stringdist Fast Path
# ==============================================================================

#' Compute pairwise dissimilarity matrix using stringdist (C-level)
#' @noRd
.dissimilarity_matrix_stringdist <- function(enc, dissimilarity, lambda,
                                             q = 2L, p = 0.1) {
  int_mat <- enc$int_mat
  len <- enc$len
  n_states <- enc$n_states
  k <- ncol(int_mat)

  base_chars <- c(letters, LETTERS, as.character(0:9))
  sentinel_char <- base_chars[n_states + 1L]

  # Build character matrix (sentinel for NA)
  chr_mat <- matrix(sentinel_char, nrow = nrow(int_mat), ncol = k)
  obs <- !is.na(int_mat)
  chr_mat[obs] <- base_chars[int_mat[obs]]

  if (dissimilarity == "hamming") {
    # Weighted hamming not supported by stringdist - fall back to R
    if (lambda > 0) {
      return(.dissimilarity_matrix_r(enc, dissimilarity, lambda, q, p)) # nocov
    }
    # Full width strings (same length, NAs are sentinel chars)
    strings <- vapply(seq_len(nrow(chr_mat)), function(i) {
      paste0(chr_mat[i, ], collapse = "")
    }, character(1L))
    stringdist::stringdistmatrix(strings, method = "hamming")
  } else {
    # Truncate to effective length per sequence
    strings <- vapply(seq_len(nrow(chr_mat)), function(i) {
      if (len[i] == 0L) return("") # nocov
      paste0(chr_mat[i, seq_len(len[i])], collapse = "")
    }, character(1L))
    if (dissimilarity == "jw") {
      stringdist::stringdistmatrix(strings, method = "jw", p = p)
    } else if (dissimilarity %in% c("qgram", "cosine", "jaccard")) {
      # stringdist returns NaN/NA when a sequence is shorter than q
      # (no q-grams -> zero-norm vector -> cosine undefined). Replace
      # with distance 1 (maximally dissimilar from everything else,
      # 0 to self), matching the R fast path's zero-norm convention.
      d <- stringdist::stringdistmatrix(strings, method = dissimilarity, q = q)
      m <- as.matrix(d)
      if (anyNA(m)) {
        m[is.na(m)] <- 1
        diag(m) <- 0
      }
      stats::as.dist(m)
    } else {
      stringdist::stringdistmatrix(strings, method = dissimilarity)
    }
  }
}

# ==============================================================================
# 5. Input Extraction
# ==============================================================================

#' One windowed sequence per unit of a one-hot network
#'
#' The rows of a one-hot network are time points. Treating each as a sequence
#' would cluster strings of 0 and 1 and break the windows apart.
#' @noRd
.onehot_unit_sequences <- function(net) {
  data <- .resampling_transition_data(net$data, net$metadata, net$params)
  actor <- .param_get(net$params, "actor")
  if (is.null(actor)) {
    stop(errorCondition(
      paste0("This one-hot network has no `actor` or `session`, so it is one ",
             "sequence and there is nothing to cluster."),
      class = "nestimate_bad_input", call = NULL))
  }
  data$.unit <- .observed_group_id(data[actor], context = "actor/session")
  .prepare_onehot(data, cols = .param_get(net$params, "codes"), actor = ".unit",
                 window_size = .param_get(net$params, "window_size", 1L),
                 window_type = .param_get(net$params, "mode", "non-overlapping"))
}


#' Extract wide-format sequence data from various object types
#'
#' Supports: data.frame, matrix, netobject, tna, cograph_network.
#' For tna objects, decodes integer-encoded data using stored labels.
#' Rejects association-method netobjects (cor/pcor/glasso/ising) since
#' those contain numeric observations, not categorical sequences.
#' @param data Input object
#' @return A data frame of wide-format sequences
#' @noRd
.extract_sequence_data <- function(data) {
  # A bare cograph_network is coerced first, so it meets the same method and
  # one-hot checks as a netobject.
  if (inherits(data, "cograph_network") && !inherits(data, "netobject")) {
    data <- .as_netobject(data)
  }
  # --- netobject (from Nestimate) ---
  if (inherits(data, "netobject")) {
    seq_methods <- c("relative", "frequency", "co_occurrence", "attention")
    if (!is.null(data$method) && !(data$method %in% seq_methods)) {
      stop(errorCondition(paste0(sprintf(
        "build_clusters() requires sequence data, but this netobject uses method '%s'. ",
        data$method
      ), "Pass the original sequence data frame instead."),
      class = "nestimate_bad_input", call = NULL))
    }
    if (is.null(data$data)) { # nocov start
      stop(errorCondition(
        "netobject has no $data. Pass the original sequence data.",
        class = "nestimate_bad_input", call = NULL))
    } # nocov end
    if (.is_onehot_params(data$params)) return(.onehot_unit_sequences(data))
    return(as.data.frame(data$data, stringsAsFactors = FALSE))
  }

  # --- tna model (integer-encoded with $labels) ---
  if (inherits(data, "tna")) {
    if (is.null(data$data) || is.null(data$labels)) { # nocov start
      stop(errorCondition(
        "tna object missing $data or $labels.",
        class = "nestimate_bad_input", call = NULL))
    } # nocov end
    # Decode integer matrix back to state names
    decoded <- matrix(data$labels[data$data], nrow = nrow(data$data),
                      ncol = ncol(data$data))
    colnames(decoded) <- colnames(data$data)
    return(as.data.frame(decoded, stringsAsFactors = FALSE))
  }

  # --- cograph_network ---
  if (inherits(data, "cograph_network")) {
    if (is.null(data$data)) { # nocov start
      stop(errorCondition(
        "cograph_network has no $data. Pass the original sequence data.",
        class = "nestimate_bad_input", call = NULL))
    } # nocov end
    return(as.data.frame(data$data, stringsAsFactors = FALSE))
  }

  # --- raw data.frame or matrix ---
  if (is.data.frame(data) || is.matrix(data)) {
    return(as.data.frame(data, stringsAsFactors = FALSE))
  }

  stop(errorCondition(
    paste0(
      "Unsupported input type: ",
      paste(class(data), collapse = ", "),
      ". Expected data.frame, matrix, netobject, tna, or cograph_network."),
    class = "nestimate_bad_input", call = NULL))
}

#' Split a wide frame into sequence positions and set-aside covariates
#'
#' Applies the role rules of \code{build_network()} so a verb that reads wide
#' frames never treats an id or a covariate as a sequence position.
#' @param data A wide data frame.
#' @param verb Name of the calling verb, for the message.
#' @return List with `sequences` (data frame of positions) and `metadata`
#'   (data frame of the other columns, or NULL when there are none).
#' @noRd
.split_sequence_columns <- function(data, verb) {
  .refuse_undeclared_roles(data, paste0(verb, "()"), via_network = TRUE,
                           reads_onehot = TRUE)
  selected <- .build_sequence_cols(data)
  if (is.null(selected) || setequal(selected, names(data))) {
    return(list(sequences = data, metadata = NULL))
  }
  aside <- setdiff(names(data), selected)
  if (length(selected) < 2L) {
    stop(errorCondition(
      sprintf(paste0("%s() found fewer than two sequence columns: `%s` were ",
                     "read as covariates. Pass a frame holding only the ",
                     "sequence columns."),
              verb, paste(aside, collapse = "`, `")),
      class = "nestimate_bad_input", call = NULL
    ))
  }
  message(sprintf("%s(): treating `%s` as covariates, not sequence positions.",
                  verb, paste(aside, collapse = "`, `")))
  list(sequences = .keep_window_attrs(data[, selected, drop = FALSE], data),
       metadata = data[, aside, drop = FALSE])
}

# ==============================================================================
# 6. Main Function
# ==============================================================================

#' Cluster Sequences by Dissimilarity
#'
#' @description
#' `build_clusters()` splits a set of sequences into `k` groups of similar
#' sequences. It computes the dissimilarity between every pair of sequences
#' (the Hamming distance, an edit distance, a q-gram distance or the
#' Jaro-Winkler distance) and partitions that matrix with PAM (partitioning
#' around medoids) or a hierarchical method. It takes wide sequence data, a
#' `netobject` or bare `cograph_network` built from sequences, or a `tna`
#' model. A one-hot network built with `"relative"`, `"frequency"` or
#' `"co_occurrence"` is clustered by unit, each actor's or session's windows
#' forming one sequence; a `wtna` network is refused (see `data`).
#' Covariates can be related to the clusters after the fit (see
#' `covariates`); they never change the assignments.
#'
#' Every dissimilarity and method returns `k` groups, whether or not the
#' sequences fall into groups. The average silhouette width (Rousseeuw,
#' 1987), which `print()` reports, compares each sequence's mean
#' dissimilarity to its own cluster with that to the nearest other cluster:
#' it runs from -1 to 1, values near 1 meaning well-separated clusters and
#' values near 0 overlapping ones. It describes one solution and is not a
#' test. Compare values of `k`, dissimilarities and methods with
#' [cluster_choice()], check the classification with [cluster_diagnostics()],
#' and build one network per cluster with [build_network()] (or
#' [cluster_network()] in one call) to see how the groups differ.
#'
#' @details
#' With the \pkg{stringdist} package installed, the distances are computed
#' by `stringdist::stringdistmatrix()` in C, many times faster than the base
#' R implementation for the edit distances. The base R implementation is
#' used when \pkg{stringdist} is absent, for weighted Hamming, and when the
#' data hold 62 or more states.
#'
#' @param data The sequences: a data frame or matrix with one row per
#'   sequence and one column per position, holding state names; a
#'   `netobject` or bare `cograph_network` built from sequences with method
#'   `"relative"`, `"frequency"`, `"co_occurrence"` or `"attention"`, whose
#'   stored sequences are used; or a `tna` model, whose integer-coded
#'   sequences are decoded to their labels. Columns of a data frame take the
#'   same roles as in [build_network()]: an id or a covariate is not a
#'   position, and the columns set aside are named in a message and kept as
#'   metadata. A one-hot network needs an `actor` or `session` to define its
#'   units. A one-hot network without one, and a network of any other
#'   method, raise an error of class `nestimate_bad_input`.
#' @param k The number of clusters, a whole number from 2 to one less than
#'   the number of sequences. Required.
#' @param dissimilarity How two sequences are compared. `"hamming"`
#'   (default) counts the positions where they differ. `"lv"` (Levenshtein)
#'   counts the insertions, deletions and substitutions that turn one into
#'   the other; `"osa"` (optimal string alignment) and `"dl"`
#'   (Damerau-Levenshtein) also count a swap of two adjacent states as one
#'   edit; `"lcs"` counts the states left unpaired by their longest common
#'   subsequence. `"qgram"`, `"cosine"` and `"jaccard"` compare the runs of
#'   `q` consecutive states (q-grams) the sequences contain: the summed
#'   absolute difference of their counts, one minus the cosine similarity of
#'   the counts, and one minus the share of distinct q-grams they have in
#'   common. `"jw"` is the Jaro-Winkler distance, which weighs agreement in
#'   the first states more (see `p`).
#' @param method The clustering algorithm. `"pam"` (default) partitions
#'   around medoids: it chooses `k` sequences as centres and assigns every
#'   sequence to the nearest one (Kaufman & Rousseeuw, 1990). `"ward.D2"`,
#'   `"ward.D"`, `"complete"`, `"average"`, `"single"`, `"mcquitty"`,
#'   `"median"` and `"centroid"` build a tree with [stats::hclust()] using
#'   that linkage and cut it into `k` groups.
#' @param na_syms The symbols read as missing values, besides `NA`. Default
#'   `c("*", "%")`. A missing cell is not dropped from a comparison: it is
#'   treated as one more state, so two missing cells at the same position
#'   match and a missing cell against an observed state is a mismatch (a
#'   distance of 1 for Hamming). This is the usual rule for aligned
#'   sequences, because dropping missing cells pair by pair would change
#'   each pair's effective length and break the metric. For every
#'   dissimilarity except Hamming, a sequence ends at its last observed
#'   state, so trailing missing cells are not compared. To treat missing
#'   cells differently, drop or recode them before the call.
#' @param weighted Whether to weight the Hamming distance by position, so
#'   that a difference early in the sequences counts more than a late one
#'   (see `lambda`). Default `FALSE`. With any dissimilarity other than
#'   `"hamming"`, `TRUE` raises an error.
#' @param lambda The decay rate of the position weights when
#'   `weighted = TRUE`: position \eqn{t} gets weight
#'   \eqn{e^{-\lambda (t - 1)}}{exp(-lambda * (t - 1))}, so a larger value
#'   puts more of the distance on the earlier positions. A finite
#'   non-negative number; default `1`. Ignored, and stored as `0`, when
#'   `weighted = FALSE`.
#' @param seed An integer used as the random seed for clustering; the
#'   session's own random stream is left as it was. Default `NULL`. PAM and the
#'   hierarchical methods are deterministic, so it does not change the
#'   clusters; it is stored with the result.
#' @param q The length of the q-grams compared by `"qgram"`, `"cosine"` and
#'   `"jaccard"`, a positive whole number. Default `2`. Ignored by the other
#'   dissimilarities.
#' @param p The Winkler prefix weight of `"jw"`: how much agreement in the
#'   first states (up to four) lowers the Jaro distance, from 0 (the plain
#'   Jaro distance) to 0.25. Default `0.1`. Ignored by the other
#'   dissimilarities.
#' @param covariates Variables describing each sequence (such as an
#'   achievement level or a course), related to the clusters after they are
#'   found; they never change the assignments. Give a column name, a
#'   character vector of names, a one-sided formula such as
#'   `~ Age + Gender` (or the string `"Age + Gender"`), or a data frame with
#'   one row per sequence, all of whose columns are used. For a `netobject`
#'   or `cograph_network`, names are looked up in its metadata first and
#'   then among the columns of its data, so a typical call is
#'   `build_clusters(net, k = 3, covariates = "Achiever")`; for a data frame
#'   or matrix, among its columns. A `tna` model takes only the data-frame
#'   form. Named covariate columns are taken out of the sequence data before
#'   the distances are computed, and so are the columns of a supplied data
#'   frame that match a sequence column in name and values; a sequence
#'   column that shares only a name is kept. For a network, a name found in
#'   its metadata leaves a same-named data column in place as a sequence
#'   position. Rows with a missing covariate
#'   are left out of the analysis with a warning, and a constant covariate
#'   raises an error. How the analysis is run is set by `estimator`.
#'   Default `NULL` (none).
#' @param estimator How the relation between the covariates and the
#'   clusters is estimated; it never affects the clustering. `"auto"`
#'   (default) fits a multinomial logit of cluster membership on the
#'   covariates: with Firth's bias-reduced penalised likelihood
#'   (`brglm2::brmultinom()`, finite under quasi-complete separation) when a
#'   cell of the cluster-by-covariate table of a categorical covariate holds
#'   fewer than 5 sequences, and with `nnet::multinom()` otherwise.
#'   `"firth"` forces the Firth fit, which is about 200 times slower on
#'   well-conditioned data. `"multinom"` forces classical maximum likelihood
#'   and warns that rare cells can give extreme odds ratios with degenerate
#'   confidence intervals. `"chisq"` runs descriptive tests instead of a
#'   logit, as WeightedCluster does: for each categorical covariate a
#'   chi-square test with Cramer's V and standardised adjusted residuals,
#'   for each numeric one a Kruskal-Wallis test with eta-squared. Ignored
#'   without `covariates`.
#' @param ... Not used: an extra argument raises an error, here and in the
#'   `print()`, `summary()` and `plot()` methods for a `net_clustering`.
#'   `print()` of a `tidy_covariates` table ignores it.
#' @return A `net_clustering` object. `print()` shows the method,
#'   dissimilarity and average silhouette width, then one row per cluster
#'   with its size (count and percentage), mean within-cluster distance and,
#'   for PAM, its medoid; after a covariate analysis it adds a line
#'   naming the covariates. [summary()] prints the method, dissimilarity and
#'   average silhouette width and a per-cluster table of sizes and mean
#'   within-cluster distances, followed by the covariate analysis when there
#'   is one, and returns a per-cluster table, or the covariate table after a
#'   covariate analysis (see the methods' return value below). `plot()` draws
#'   the silhouette widths, a map of the distances, the distance matrix or
#'   the covariate effects. [build_network()] turns the clustering into one
#'   network per cluster (a `netobject_group`), and
#'   [cluster_diagnostics()] and [session_ids()] give the classification
#'   quality and each sequence's cluster. The object's fields:
#' \describe{
#'   \item{data}{The sequence positions that were clustered, one row per
#'     sequence: covariates and set-aside columns are removed, and a
#'     one-hot network gives one row per unit.}
#'   \item{k}{Number of clusters.}
#'   \item{assignments}{Integer vector of cluster memberships (1 to
#'     \code{k}), one per sequence in row order.}
#'   \item{silhouette}{Average silhouette width over all sequences.}
#'   \item{sizes}{Named integer vector of cluster sizes.}
#'   \item{method}{Clustering method used.}
#'   \item{dissimilarity}{Dissimilarity used.}
#'   \item{distance}{The pairwise dissimilarities, a \code{dist} object.}
#'   \item{medoids}{For PAM, the row numbers of the \code{k} medoid
#'     sequences; NULL for the hierarchical methods.}
#'   \item{seed}{Seed given (or NULL).}
#'   \item{weighted}{Whether weighted Hamming was used.}
#'   \item{lambda}{Decay rate used (0 when not weighted).}
#'   \item{covariates}{The covariate analysis (a list; see
#'     \code{estimator}), or NULL when \code{covariates = NULL}.}
#'   \item{network_method, build_args}{For \code{netobject} input, the
#'     source network's method and stored build arguments, so per-cluster
#'     networks can be rebuilt the same way. NULL otherwise, and
#'     \code{build_args} is NULL for a one-hot network.}
#'   \item{metadata}{For \code{netobject} input, its per-sequence metadata,
#'     one row per clustered sequence, so \code{\link{session_ids}} can name
#'     each sequence. For a data frame, the columns set aside as covariates
#'     (see \code{data}). NULL otherwise, and for a one-hot network.}
#'   \item{htna_partition}{For HTNA input, the preserved node-to-actor
#'     partition used to restore HTNA children when networks are built.}
#' }
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
#'                    V3 = c("C","A","B","C","B"))
#' cl <- build_clusters(seqs, k = 2)
#' cl
#' \donttest{
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:3], 20, TRUE), V2 = sample(LETTERS[1:3], 20, TRUE),
#'   V3 = sample(LETTERS[1:3], 20, TRUE), V4 = sample(LETTERS[1:3], 20, TRUE)
#' )
#' cl <- build_clusters(seqs, k = 2)
#' print(cl)
#' summary(cl)
#' }
#'
#' @references
#' Kaufman, L., & Rousseeuw, P. J. (1990). \emph{Finding Groups in Data: An
#' Introduction to Cluster Analysis}. Wiley.
#'
#' Rousseeuw, P. J. (1987). Silhouettes: A graphical aid to the
#' interpretation and validation of cluster analysis. \emph{Journal of
#' Computational and Applied Mathematics}, 20, 53--65.
#' \doi{10.1016/0377-0427(87)90125-7}
#'
#' @importFrom utils combn
#' @export
build_clusters <- function(data, k, dissimilarity = "hamming", method = "pam",
                         na_syms = c("*", "%"), weighted = FALSE, lambda = 1,
                         seed = NULL, q = 2L, p = 0.1,
                         covariates = NULL,
                         estimator = c("auto", "firth", "multinom", "chisq"),
                         ...) {
  estimator <- match.arg(estimator)
  dots <- list(...)
  if (length(dots) > 0L) {
    dot_names <- names(dots)
    dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
    stop(errorCondition(paste0("build_clusters() got unsupported argument",
         if (length(dots) == 1L) ": " else "s: ",
         paste(dot_names, collapse = ", ")),
         class = "nestimate_bad_input", call = NULL))
  }

  # --- Extract sequence data from supported objects ---
  raw_data <- data
  onehot_net <- inherits(raw_data, "netobject") &&
    .is_onehot_params(raw_data$params)
  data <- .extract_sequence_data(data)
  # Resolve post-hoc covariates before constructing distances: named columns
  # in the sequence frame must not become additional sequence positions.
  cov_resolved <- .resolve_covariates(covariates, raw_data, nrow(data))
  if (!is.null(cov_resolved)) {
    cov_cols <- names(cov_resolved$cov_df)
    if (is.data.frame(covariates)) {
      # A supplied frame may be a subset of the raw input. Exclude those
      # columns too, but preserve unrelated positions sharing only a name.
      overlap <- intersect(cov_cols, names(data))
      cov_cols <- overlap[vapply(overlap, function(nm) {
        identical(as.character(data[[nm]]), as.character(covariates[[nm]]))
      }, logical(1L))]
    }
    # Metadata takes precedence in .resolve_covariates(); a same-named state
    # column remains a position when the requested covariate came from there.
    if (!is.data.frame(covariates) &&
        inherits(raw_data, c("netobject", "cograph_network"))) {
      cov_cols <- setdiff(cov_cols, names(raw_data$metadata))
    }
    data <- .keep_window_attrs(data[, setdiff(names(data), cov_cols), drop = FALSE],
                               data)
    if (ncol(data) == 0L) {
      stop(errorCondition(
        "No sequence columns remain after excluding post-hoc covariates.",
        class = "nestimate_bad_input", call = NULL))
    }
  }
  # Same roles as build_network(): an id or a covariate is not a position.
  set_aside <- NULL
  if (is.data.frame(raw_data)) {
    roles <- .split_sequence_columns(data, "build_clusters")
    data <- roles$sequences
    set_aside <- roles$metadata
  }

  # --- Input validation ---
  # Use named-condition stopifnot for type checks so the failure message
  # names the offending argument instead of dumping the predicate.
  .assert(
    "'k' must be a single numeric value"            = is.numeric(k) && length(k) == 1L,
    "'dissimilarity' must be a single character"    = is.character(dissimilarity) && length(dissimilarity) == 1L,
    "'method' must be a single character"           = is.character(method) && length(method) == 1L,
    "'na_syms' must be a character vector"          = is.character(na_syms),
    "'weighted' must be a single logical"           = is.logical(weighted) && length(weighted) == 1L,
    "'lambda' must be a single numeric"             = is.numeric(lambda) && length(lambda) == 1L,
    "'q' must be a single numeric"                  = is.numeric(q) && length(q) == 1L,
    "'p' must be a single numeric"                  = is.numeric(p) && length(p) == 1L
  )
  if (!is.finite(k) || k != floor(k)) {
    stop(errorCondition(
      "'k' must be a whole finite number.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (is.na(weighted)) {
    stop(errorCondition(
      "'weighted' must be TRUE or FALSE, not NA.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.finite(lambda) || lambda < 0) {
    stop(errorCondition(
      "'lambda' must be a finite non-negative number.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.finite(q) || q != floor(q) || q < 1) {
    stop(errorCondition(
      "'q' must be a positive whole number.",
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.finite(p) || p < 0 || p > 0.25) {
    stop(errorCondition(
      "'p' must be a finite number between 0 and 0.25.",
      class = "nestimate_bad_input", call = NULL))
  }
  k <- as.integer(k)
  q <- as.integer(q)
  n <- nrow(data)
  if (k < 2L) {
    stop(errorCondition(sprintf(
      "'k' must be at least 2 (got k = %d). Clustering with one group is not meaningful.",
      k
    ), class = "nestimate_bad_input", call = NULL))
  }
  if (k > n - 1L) {
    stop(errorCondition(sprintf(
      paste0("'k' (= %d) must be <= n - 1 where n = %d sequences. ",
             "Either reduce k or supply more sequences (need at least %d)."),
      k, n, k + 1L
    ), class = "nestimate_bad_input", call = NULL))
  }

  dissimilarity <- match.arg(dissimilarity, .clustering_metrics)
  method <- match.arg(method, .clustering_methods)

  if (weighted && dissimilarity != "hamming") {
    stop(errorCondition(
      "Weighting is only supported for Hamming distance.",
      class = "nestimate_bad_input", call = NULL))
  }
  lambda <- if (weighted) lambda else 0

  if (!is.null(seed)) .local_seed(seed)

  # --- Encode sequences ---
  enc <- .encode_sequences(data, na_syms)

  if (enc$n_states == 0L) {
    stop(errorCondition(
      paste0(
        "No observed sequence states were found after applying na_syms = ",
        sprintf("c(%s)", paste(shQuote(na_syms), collapse = ", ")),
        ". Check that the input data contains at least one non-missing ",
        "state cell, or adjust na_syms to exclude symbols that should ",
        "be treated as observed states."),
      class = "nestimate_bad_input", call = NULL))
  }

  # --- Compute distance matrix ---
  use_stringdist <- requireNamespace("stringdist", quietly = TRUE) &&
    enc$n_states < 62L &&  # need room for sentinel char
    !(dissimilarity == "hamming" && lambda > 0)  # weighted hamming not in stringdist

  if (use_stringdist) {
    dist_mat <- .dissimilarity_matrix_stringdist(
      enc, dissimilarity, lambda, q, p
    )
  } else {
    dist_mat <- .dissimilarity_matrix_r(enc, dissimilarity, lambda, q, p)
  }

  # --- Cluster ---
  medoids <- NULL
  if (method == "pam") {
    clust_result <- cluster::pam(dist_mat, diss = TRUE, k = k)
    assignments <- clust_result$clustering
    silhouette_score <- clust_result$silinfo$avg.width
    medoids <- clust_result$id.med
  } else {
    hc <- stats::hclust(dist_mat, method = method)
    assignments <- stats::cutree(hc, k = k)
    sil <- cluster::silhouette(assignments, dist = dist_mat)
    silhouette_score <- mean(sil[, 3L])
  }

  sizes <- c(table(assignments))

  # --- Post-hoc covariate analysis ---
  cov_result <- NULL
  if (!is.null(cov_resolved)) {
    cov_result <- .run_covariate_analysis(
      assignments, cov_resolved$cov_df, cov_resolved$rhs, k,
      estimator = estimator
    )
  }

  result <- structure(
    list(
      data = data,
      k = k,
      assignments = assignments,
      silhouette = silhouette_score,
      sizes = sizes,
      method = method,
      dissimilarity = dissimilarity,
      distance = dist_mat,
      medoids = medoids,
      seed = seed,
      weighted = weighted,
      lambda = lambda,
      covariates = cov_result,
      network_method = if (inherits(raw_data, "netobject")) raw_data$method else NULL,
      build_args     = if (inherits(raw_data, "netobject") && !onehot_net) raw_data$build_args else NULL,
      metadata       = if (onehot_net) NULL else if (inherits(raw_data, "netobject")) raw_data$metadata else set_aside,
      htna_partition = .capture_htna_partition(raw_data)
    ),
    class = "net_clustering"
  )
  result
}

#' Cluster sequence data (deprecated alias)
#'
#' `cluster_data()` is the old name of [build_clusters()], renamed in
#' Nestimate 0.4.3. It is kept so older tutorials and the historical pkgdown
#' reference still run: every call raises a deprecation warning and passes
#' its arguments to [build_clusters()] unchanged.
#'
#' @param ... Arguments passed unchanged to [build_clusters()].
#' @return The `net_clustering` object that [build_clusters()] returns.
#' @seealso \code{\link{build_clusters}}, \code{\link{cluster_network}},
#'   \code{\link{cluster_mmm}}.
#' @keywords internal
#' @export
cluster_data <- function(...) {
  .Deprecated("build_clusters",
              package = "Nestimate",
              msg = paste0(
                "cluster_data() was renamed to build_clusters() in ",
                "Nestimate 0.4.3. Please use build_clusters() instead."))
  build_clusters(...)
}


# ==============================================================================
# 6. S3 Methods
# ==============================================================================

# ------------------------------------------------------------------------------
# Shared formatter for cluster summary tables.
#
# `cols` is a named list of vectors (length k each). The first column drives
# row count; all other columns must match that length. Returns a character
# vector of formatted lines, indented and headed by a separator rule. Used by
# print.net_clustering / print.net_mmm / print.net_mmm_clustering /
# print.netobject_group so a netobject_group built by cluster_network() and
# one built from net_mmm print with the same shape.
# ------------------------------------------------------------------------------
.cluster_table_lines <- function(cols, indent = "  ", sep = "  ") {
  stopifnot(is.list(cols), length(cols) >= 1L, !is.null(names(cols)))
  k <- length(cols[[1L]])
  if (k == 0L) return(character(0L))

  names_v <- names(cols)
  # Column widths from header + cell content (cells are already character).
  widths <- vapply(seq_along(cols), function(j) {
    cells <- as.character(cols[[j]])
    max(nchar(names_v[[j]]), if (length(cells)) max(nchar(cells)) else 0L)
  }, integer(1L))

  ncol_v <- length(cols)
  # `print(data.frame)` left-aligns; we keep the same convention but pad
  # only up to (not past) the *last* column so trailing whitespace on the
  # right edge is gone.
  pad <- function(s, w, last) if (last) s else formatC(s, width = w, flag = "-")

  build_row <- function(cells) {
    parts <- vapply(seq_len(ncol_v), function(j) {
      pad(cells[[j]], widths[[j]], last = j == ncol_v)
    }, character(1L))
    paste(parts, collapse = sep)
  }

  header <- build_row(as.list(names_v))
  rows   <- vapply(seq_len(k), function(i) {
    build_row(lapply(cols, function(col) as.character(col[[i]])))
  }, character(1L))

  paste0(indent, c(header, rows))
}

# Helper: format size + percentage as a "30 (50.0%)" string.
.fmt_size_pct <- function(sizes, total = sum(sizes)) {
  pct <- if (isTRUE(total > 0)) sizes / total * 100 else rep(0, length(sizes))
  sprintf("%d (%4.1f%%)", as.integer(sizes), pct)
}

# Per-cluster mean within-cluster distance. Singletons return 0 (no
# within-distance defined). Used by print.net_clustering,
# summary.net_clustering, cluster_diagnostics.net_clustering, and
# cluster_choice's per-row builder.
.per_cluster_within_dist <- function(distance, assignments, k = max(assignments)) {
  if (is.null(distance)) return(rep(NA_real_, k))
  dmat <- as.matrix(distance)
  vapply(seq_len(k), function(cl) {
    members <- which(assignments == cl)
    if (length(members) > 1L) {
      sub <- dmat[members, members]
      mean(sub[lower.tri(sub)])
    } else {
      0
    }
  }, numeric(1L))
}

#' @rdname build_clusters
#' @param x For the `print()` and `plot()` methods: a `net_clustering`; for
#'   `print()`, also the `tidy_covariates` table that `summary()` returns.
#' @param digits For `print()` of a `net_clustering`: decimal places of the
#'   silhouette width and the within-cluster distances, a non-negative whole
#'   number. Default `3`.
#' @param object For `summary()`: a `net_clustering`.
#' @param type For `plot()`: the view. `"silhouette"` (default) draws each
#'   sequence's silhouette width as a bar, grouped by cluster, with the
#'   average as a dashed line; `"mds"` places the sequences on the first two
#'   dimensions of a classical multidimensional scaling of the distances,
#'   coloured by cluster; `"heatmap"` draws the distance matrix with the
#'   sequences ordered by cluster; `"predictors"` draws a forest plot of the
#'   covariates' odds ratios (log scale) for each cluster against the
#'   reference cluster. `"predictors"` needs a fit with `covariates` and a
#'   logit estimator; with `estimator = "chisq"`, which has no odds ratios,
#'   it raises an error that points to [summary()] for the tests.
#' @param combined For `plot()` with `type = "predictors"`: `TRUE` (default)
#'   draws one plot with a panel per non-reference cluster; `FALSE` returns
#'   a list of plots, one per non-reference cluster. Ignored by the other
#'   views.
#' @return For the `print()` methods: the input, invisibly. `print()` of a
#'   `tidy_covariates` table prints a header naming the estimator, then the
#'   table. For `summary()`: a data frame with one row per cluster and
#'   columns `cluster`, `size` and `mean_within_dist`, returned visibly
#'   after the printed summary. For a clustering fitted with `covariates`,
#'   the printed summary adds the covariate profiles and the logit or test
#'   tables, and a `tidy_covariates` data frame is returned invisibly
#'   instead: one row per cluster and predictor (`cluster`, `predictor`,
#'   `OR`, `CI_low`, `CI_high`, `p`, `sig`) for a logit, or one row per
#'   covariate (`variable`, `type`, `statistic`, `df`, `p`, `effect`,
#'   `effect_kind`, `sig`) for `estimator = "chisq"`, with the per-cluster
#'   table, fit statistics and profiles attached as attributes. For
#'   `plot()`: a `ggplot` object, drawn and returned invisibly; with
#'   `type = "predictors"` and `combined = FALSE`, a list of `ggplot`
#'   objects, one per non-reference cluster and named by it, returned
#'   invisibly without drawing.
#' @export
print.net_clustering <- function(x, digits = 3L, ...) {
  .net_clustering_check_unused_dots("print.net_clustering", ...)
  .net_clustering_check_digits(digits)
  digits <- as.integer(digits)
  k <- as.integer(x$k)
  sizes <- as.integer(x$sizes)
  n_total <- sum(sizes)

  diss_lab <- if (isTRUE(x$weighted)) {
    sprintf("%s (weighted, lambda = %g)", x$dissimilarity, x$lambda)
  } else {
    x$dissimilarity
  }

  cat(sprintf("Sequence Clustering [%s]\n", x$method))
  cat(sprintf("  Sequences: %d  |  Clusters: %d\n", n_total, k))
  cat(sprintf("  Dissimilarity: %s\n", diss_lab))
  if (!is.null(x$silhouette)) {
    cat(sprintf("  Quality: silhouette = %.*f\n",
                digits, x$silhouette))
  }

  # Per-cluster mean within-cluster distance, when distance matrix is present.
  mean_within <- if (!is.null(x$distance) && !is.null(x$assignments) && k > 0L) {
    .per_cluster_within_dist(x$distance, x$assignments, k)
  } else {
    rep(NA_real_, k)
  }

  cat("\n")
  cols <- list(
    Cluster = sprintf("%d", seq_len(k)),
    N       = .fmt_size_pct(sizes, n_total)
  )
  if (!all(is.na(mean_within))) {
    cols[["Mean within-dist"]] <- sprintf(paste0("%.", digits, "f"),
                                          mean_within)
  }
  if (!is.null(x$medoids) && length(x$medoids) == k) {
    cols[["Medoid"]] <- as.character(x$medoids)
  }
  cat(paste(.cluster_table_lines(cols), collapse = "\n"), "\n", sep = "")

  if (!is.null(x$covariates)) {
    cat(.covariate_line(x$covariates, "post-hoc"))
  }

  invisible(x)
}

# The covariates of a post-hoc or EM covariate analysis, as the user named
# them. They are read from the profiles, which every estimator fills; the
# coefficient table holds design columns (one per dummy level) and a
# chi-square analysis has none.
#' @noRd
.covariate_names <- function(cov) {
  unique(c(cov$profiles$numeric$variable, cov$profiles$categorical$variable))
}

# The "Covariates:" line shared by the print methods of net_clustering,
# net_mmm, net_mmm_clustering and net_cluster_diagnostics.
#' @noRd
.covariate_line <- function(cov, label) {
  cov_names <- .covariate_names(cov)
  sprintf("\n  Covariates: %s (%s, %d %s)\n",
          paste(cov_names, collapse = ", "), label, length(cov_names),
          if (length(cov_names) == 1L) "covariate" else "covariates")
}

# The forest plot draws odds ratios, which a chi-square covariate analysis
# does not estimate.
#' @noRd
.require_covariate_coefficients <- function(cov, type) {
  if (!is.null(cov$coefficients)) return(invisible(TRUE))
  stop(errorCondition(
    paste0(
      "plot(..., type = \"", type, "\") draws the odds ratios of a logit ",
      "covariate analysis (estimator = \"auto\", \"firth\" or ",
      "\"multinom\"). This analysis used estimator = \"",
      cov$fit$method %||% "chisq", "\", which estimates none; summary() ",
      "shows its tests."),
    class = "nestimate_bad_input", call = NULL))
}

#' @rdname build_clusters
#' @export
summary.net_clustering <- function(object, ...) {
  .net_clustering_check_unused_dots("summary.net_clustering", ...)
  k <- object$k
  assignments <- object$assignments
  sizes <- as.integer(tabulate(assignments, nbins = k))
  cluster_stats <- data.frame(
    cluster          = seq_len(k),
    size             = sizes,
    mean_within_dist = .per_cluster_within_dist(object$distance,
                                                 assignments, k),
    stringsAsFactors = FALSE
  )
  rownames(cluster_stats) <- NULL

  cat("Sequence Clustering Summary\n")
  cat("  Method:       ", object$method, "\n")
  cat("  Dissimilarity:", object$dissimilarity, "\n")
  cat("  Silhouette:   ", round(object$silhouette, 4), "\n\n")
  cat("Per-cluster statistics:\n")
  print(cluster_stats, row.names = FALSE)

  if (!is.null(object$covariates)) {
    cat("\n")
    .print_covariate_profiles(object$covariates)
    cat("\nNote: Covariates are post-hoc and do not influence cluster assignments.\n")
    # Return a tidy data.frame of the inferential output. Cluster sizes,
    # fit stats, residuals, and the full untidied list are attached as
    # attributes so power users still have access without spelunking
    # through `$covariates`.
    return(invisible(.tidy_covariates(object$covariates,
                                       cluster_stats = cluster_stats)))
  }

  cluster_stats
}

#' @rdname build_clusters
#' @import ggplot2
#' @export
plot.net_clustering <- function(x, type = c("silhouette", "mds", "heatmap",
                                             "predictors"),
                                 combined = TRUE, ...) {
  .net_clustering_check_unused_dots("plot.net_clustering", ...)
  type <- match.arg(type)
  .assert(is.logical(combined), length(combined) == 1L, !is.na(combined))

  if (type == "predictors") {
    if (is.null(x$covariates)) {
      stop(errorCondition(
        "No covariate analysis found. Run build_clusters() with covariates.",
        class = "nestimate_bad_input", call = NULL))
    }
    .require_covariate_coefficients(x$covariates, "predictors")
    return(.plot_covariate_forest(
      x$covariates$coefficients,
      sprintf("Covariate Effects on Cluster Membership (ref: Cluster %s)",
              x$covariates$fit$reference_cluster),
      combined = combined
    ))
  }

  if (type == "silhouette") {
    sil <- cluster::silhouette(x$assignments, dist = x$distance)
    sil_df <- data.frame(
      obs = seq_len(nrow(sil)),
      cluster = factor(sil[, 1L]),
      width = sil[, 3L]
    )
    # Order by cluster then width within cluster
    sil_df <- sil_df[order(sil_df$cluster, -sil_df$width), ]
    sil_df$rank <- seq_len(nrow(sil_df))

    p <- ggplot(sil_df, aes(x = .data$rank, y = .data$width,
                             fill = .data$cluster)) +
      geom_col(width = 1) +
      coord_flip() +
      labs(x = "Observation", y = "Silhouette Width",
           title = sprintf("Silhouette Plot (avg = %.3f)", x$silhouette),
           fill = "Cluster") +
      theme_minimal(base_size = 12) +
      geom_hline(yintercept = x$silhouette, linetype = "dashed", alpha = 0.5)

  } else if (type == "mds") {
    mds <- stats::cmdscale(x$distance, k = 2L)
    mds_df <- data.frame(
      dim1 = mds[, 1L],
      dim2 = mds[, 2L],
      cluster = factor(x$assignments)
    )

    p <- ggplot(mds_df, aes(x = .data$dim1, y = .data$dim2,
                             colour = .data$cluster)) +
      geom_point(size = 2.5) +
      labs(x = "MDS Dimension 1", y = "MDS Dimension 2",
           title = sprintf("MDS Projection (k = %d)", x$k),
           colour = "Cluster") +
      theme_minimal(base_size = 12)

  } else {
    # Heatmap - order by cluster assignment
    ord <- order(x$assignments)
    dist_mat <- as.matrix(x$distance)[ord, ord]
    n <- nrow(dist_mat)

    # Build long-format data frame
    idx <- expand.grid(row = seq_len(n), col = seq_len(n))
    idx$distance <- as.vector(dist_mat)
    idx$row <- factor(idx$row, levels = seq_len(n))
    idx$col <- factor(idx$col, levels = rev(seq_len(n)))

    p <- ggplot(idx, aes(x = .data$row, y = .data$col,
                          fill = .data$distance)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "steelblue") +
      labs(x = "Sequence", y = "Sequence",
           title = "Distance Matrix (ordered by cluster)",
           fill = "Distance") +
      theme_minimal(base_size = 12) +
      theme(axis.text = element_blank(), axis.ticks = element_blank())
  }

  print(p)
  invisible(p)
}

.net_clustering_check_unused_dots <- function(method, ...) {
  dots <- list(...)
  if (!length(dots)) {
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
  stop(errorCondition(paste0(
    method, "() got unsupported argument",
    if (length(dots) == 1L) ": " else "s: ",
    paste(dot_names, collapse = ", ")),
    class = "nestimate_bad_input", call = NULL)
  )
}

.net_clustering_check_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits != floor(digits) || digits < 0L) {
    stop(errorCondition(
      "'digits' must be a single non-negative whole finite number.",
      class = "nestimate_bad_input", call = NULL))
  }
  invisible(TRUE)
}

# ==============================================================================
# 7. Covariate Analysis (post-hoc)
# ==============================================================================

#' Resolve covariates into a data.frame
#'
#' Handles five input forms: formula, character vector, string, data.frame,
#' NULL. Looks up column names in metadata or raw data as needed.
#' @param covariates User-supplied covariates argument
#' @param raw_data Original input to build_clusters (before extraction)
#' @param n Number of sequences
#' @return A list with \code{cov_df} (data.frame) and \code{formula} (formula),
#'   or NULL if covariates is NULL.
#' @noRd
.resolve_covariates <- function(covariates, raw_data, n) {
  if (is.null(covariates)) return(NULL)

  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop(errorCondition(paste0( # nocov start
      "The 'nnet' package is required for covariate analysis. ",
      "Install it with install.packages('nnet')."),
      class = "nestimate_missing_package", call = NULL)) # nocov end
  }

  # --- Covariates supplied as a data.frame: use all columns ---
  if (is.data.frame(covariates)) {
    if (nrow(covariates) != n) {
      stop(errorCondition(sprintf(
        "Covariate data.frame has %d rows but sequence data has %d rows.",
        nrow(covariates), n
      ), class = "nestimate_bad_input", call = NULL))
    }
    rhs <- paste(names(covariates), collapse = " + ")
    return(list(cov_df = covariates, rhs = rhs))
  }

  # --- Parse into column names ---
  if (inherits(covariates, "formula")) {
    rhs_text <- deparse(covariates[[length(covariates)]])
    cov_names <- all.vars(covariates)
  } else if (is.character(covariates) && length(covariates) == 1L &&
             grepl("+", covariates, fixed = TRUE)) {
    # String like "Age + Gender"
    rhs_text <- covariates
    cov_names <- all.vars(stats::as.formula(paste("~", covariates)))
  } else if (is.character(covariates)) {
    # Character vector like c("Age", "Gender")
    cov_names <- covariates
    rhs_text <- paste(cov_names, collapse = " + ")
  } else {
    stop(errorCondition(
      "covariates must be a formula, character vector, string, or data.frame.",
      class = "nestimate_bad_input", call = NULL))
  }

  if (length(cov_names) == 0L) {
    stop(errorCondition(
      "No covariate variables specified.",
      class = "nestimate_bad_input", call = NULL))
  }

  # --- Look up columns in metadata or raw data ---
  cov_df <- NULL

  # Try netobject / cograph_network: walk $metadata first, then $data columns.
  # Both classes carry the same slot shape, so we resolve from either.
  if (inherits(raw_data, "netobject") ||
      inherits(raw_data, "cograph_network")) {
    avail_df <- NULL
    avail_names <- character(0L)
    if (!is.null(raw_data$metadata) && is.data.frame(raw_data$metadata)) {
      avail_df <- raw_data$metadata
      avail_names <- names(raw_data$metadata)
    }
    if (!is.null(raw_data$data) && is.data.frame(raw_data$data)) {
      extra <- setdiff(names(raw_data$data), avail_names)
      if (length(extra) > 0L) {
        avail_names <- c(avail_names, extra)
        avail_df <- if (is.null(avail_df)) {
          raw_data$data[, extra, drop = FALSE]
        } else {
          cbind(avail_df, raw_data$data[, extra, drop = FALSE])
        }
      }
    }
    if (length(avail_names) == 0L) {
      stop(errorCondition(
        paste0(
          "Network object has no $metadata or $data slot to resolve ",
          "column-name covariates against. Pass covariates as a ",
          "data.frame instead."),
        class = "nestimate_bad_input", call = NULL))
    }
    if (all(cov_names %in% avail_names)) {
      cov_df <- avail_df[, cov_names, drop = FALSE]
    } else {
      missing <- setdiff(cov_names, avail_names)
      stop(errorCondition(sprintf(
        "Covariate columns not found in network object: %s. Available: %s",
        paste(missing, collapse = ", "),
        paste(avail_names, collapse = ", ")
      ), class = "nestimate_bad_input", call = NULL))
    }
  }

  # tna objects: require data.frame covariates (no canonical metadata slot).
  if (is.null(cov_df) && inherits(raw_data, "tna")) {
    stop(errorCondition(
      paste0(
        "Column-name covariates are not supported for tna inputs. ",
        "Pass covariates as a data.frame instead."),
      class = "nestimate_bad_input", call = NULL))
  }

  # Try raw data.frame columns
  if (is.null(cov_df) && (is.data.frame(raw_data) || is.matrix(raw_data))) {
    avail <- colnames(raw_data)
    if (all(cov_names %in% avail)) {
      cov_df <- as.data.frame(raw_data[, cov_names, drop = FALSE],
                               stringsAsFactors = FALSE)
    } else {
      missing <- setdiff(cov_names, avail)
      stop(errorCondition(sprintf(
        "Covariate columns not found: %s. Available: %s",
        paste(missing, collapse = ", "),
        paste(avail, collapse = ", ")
      ), class = "nestimate_bad_input", call = NULL))
    }
  }

  if (is.null(cov_df)) {
    stop(errorCondition( # nocov start
      "Could not resolve covariate columns from the input data.",
      class = "nestimate_bad_input", call = NULL)) # nocov end
  }

  if (nrow(cov_df) != n) { # nocov start
    stop(errorCondition(sprintf(
      "Covariate data has %d rows but sequence data has %d rows.",
      nrow(cov_df), n
    ), class = "nestimate_bad_input", call = NULL)) # nocov end
  }

  list(cov_df = cov_df, rhs = rhs_text)
}

#' Compute per-cluster descriptive profiles
#' @noRd
.compute_cluster_profiles <- function(cov_df, assignments, k) {
  n_total <- length(assignments)
  num_cols <- vapply(cov_df, is.numeric, logical(1L))
  cat_cols <- !num_cols

  # --- Numeric profiles ---
  numeric_profile <- NULL
  if (any(num_cols)) {
    num_names <- names(cov_df)[num_cols]
    rows <- lapply(seq_len(k), function(cl) {
      mask <- assignments == cl
      cl_n <- sum(mask)
      cl_pct <- round(100 * cl_n / n_total, 1)
      lapply(num_names, function(vname) {
        vals <- cov_df[[vname]][mask]
        vals <- vals[!is.na(vals)]
        data.frame(
          cluster = cl, n = cl_n, pct = cl_pct, variable = vname,
          mean = mean(vals), sd = stats::sd(vals),
          median = stats::median(vals),
          stringsAsFactors = FALSE, row.names = NULL
        )
      })
    })
    numeric_profile <- do.call(rbind, unlist(rows, recursive = FALSE))
  }

  # --- Categorical profiles ---
  categorical_profile <- NULL
  if (any(cat_cols)) {
    cat_names <- names(cov_df)[cat_cols]
    # One row per cluster x variable x level, in that order.
    rows <- lapply(seq_len(k), function(cl) {
      mask <- assignments == cl
      cl_n <- sum(mask)
      lapply(cat_names, function(vname) {
        tab <- table(cov_df[[vname]][mask], useNA = "no")
        lapply(names(tab), function(lev) {
          data.frame(
            cluster = cl, n = cl_n, variable = vname,
            level = lev, count = as.integer(tab[[lev]]),
            pct = round(100 * as.integer(tab[[lev]]) / cl_n, 1),
            stringsAsFactors = FALSE, row.names = NULL
          )
        })
      })
    })
    categorical_profile <- do.call(
      rbind, unlist(unlist(rows, recursive = FALSE), recursive = FALSE))
  }

  list(numeric = numeric_profile, categorical = categorical_profile)
}

#' Run covariate analysis on cluster assignments
#'
#' Dispatches on `estimator`:
#' \itemize{
#'   \item `"auto"` (default) - resolves to `"firth"` when
#'     `.detect_separation_risk()` finds a cluster x covariate cell with
#'     fewer than 5 observations, and to `"multinom"` otherwise. The
#'     resolved choice and its reason are stored on the `profiles` object
#'     as the `auto_choice` / `auto_reason` attributes.
#'   \item `"firth"` - Firth's penalised multinomial logit via
#'     `brglm2::brmultinom`. Bias-reduced; estimates stay finite under
#'     quasi-complete separation. AIC/BIC/McFadden derived from the same
#'     log-likelihood scale as the ML fit, so they remain comparable.
#'   \item `"multinom"` - classical ML multinomial via `nnet::multinom`.
#'     Identical to the pre-2026-05-10 behaviour. A `warning()` is emitted
#'     because rare-cell pathologies (separation, infinite ORs) are silent
#'     in this estimator.
#'   \item `"chisq"` - descriptive only, WeightedCluster-style: per-covariate
#'     test of cluster association (chi-square + Cramer's V for factors,
#'     Kruskal-Wallis + eta^2 for numerics). No ORs; a separate residuals
#'     table reports per-cell standardized adjusted residuals so the user
#'     can read which (cluster, level) cells drive the association.
#' }
#'
#' Output shape:
#' \itemize{
#'   \item firth/multinom: `list(profiles, coefficients, fit, model)` -
#'     `coefficients` carries the OR table, `fit` carries AIC/BIC/McFadden.
#'   \item chisq: `list(profiles, tests, residuals, fit)` - `tests` is a
#'     per-covariate row, `residuals` is a tidy frame of standardized
#'     adjusted residuals (factor covariates only).
#' }
#'
#' @noRd
## Detect quasi-complete-separation risk in cluster x covariate cross-tabs.
## Returns TRUE when any categorical-covariate cell has fewer than
## `min_cell` observations -- the case where classical multinom diverges
## and Firth's penalty is needed.
##
## Continuous covariates are not checked here (separation in continuous
## covariates is rare and the user can force estimator = "firth"
## explicitly if known).
.detect_separation_risk <- function(fit_df, min_cell = 5L) {
  cat_cols <- vapply(setdiff(names(fit_df), "cluster"), function(nm) {
    is.factor(fit_df[[nm]]) || is.character(fit_df[[nm]])
  }, logical(1L))
  if (!any(cat_cols)) return(FALSE)
  cat_names <- names(cat_cols)[cat_cols]
  any(vapply(cat_names, function(nm) {
    tab <- table(fit_df$cluster, fit_df[[nm]])
    min(tab) < min_cell
  }, logical(1L)))
}

#' @noRd
.run_covariate_analysis <- function(assignments, cov_df, rhs, k,
                                     estimator = c("auto", "firth",
                                                   "multinom", "chisq")) {
  estimator <- match.arg(estimator)

  # --- Prepare data (shared across all three estimators) ---
  fit_df <- cov_df
  to_factor <- names(fit_df)[vapply(fit_df, function(col) {
    is.character(col) || is.ordered(col)
  }, logical(1L))]
  fit_df[to_factor] <- lapply(fit_df[to_factor], factor, ordered = FALSE)
  complete <- stats::complete.cases(fit_df)
  n_dropped <- sum(!complete)
  if (n_dropped > 0L) {
    warning(sprintf(
      "Dropped %d rows with NA covariates (%d remaining).",
      n_dropped, sum(complete)
    ), call. = FALSE)
    fit_df <- fit_df[complete, , drop = FALSE]
    assignments <- assignments[complete]
  }
  fit_df$cluster <- factor(assignments)
  cov_cols <- setdiff(names(fit_df), "cluster")
  constant <- cov_cols[vapply(cov_cols, function(nm) {
    length(unique(fit_df[[nm]][!is.na(fit_df[[nm]])])) < 2L
  }, logical(1L))]
  if (length(constant) > 0L) {
    stop(errorCondition(
      sprintf("Covariate '%s' is constant. Remove it.", constant[1L]),
      class = "nestimate_bad_input", call = NULL))
  }
  profiles_df <- if (n_dropped > 0L) cov_df[complete, , drop = FALSE] else cov_df
  profiles <- .compute_cluster_profiles(profiles_df, assignments, k)

  if (estimator == "chisq") {
    return(.covariate_chisq(fit_df, profiles))
  }

  # "auto": pick firth only when cluster x covariate cross-tab has cells
  # < 5 (separation risk where multinom would diverge); use multinom
  # otherwise (235-330x faster, identical answer when MLE is finite).
  if (estimator == "auto") {
    risk <- .detect_separation_risk(fit_df, min_cell = 5L)
    estimator <- if (risk) "firth" else "multinom"
    attr(profiles, "auto_choice") <- estimator
    attr(profiles, "auto_reason") <- if (risk) {
      "rare/zero cells detected in cluster x covariate cross-tab"
    } else {
      "no rare cells (all cluster x covariate cells >= 5)"
    }
  }

  .covariate_logit(fit_df, rhs, profiles, estimator,
                    raw_assignments = assignments)
}


# Factor + numeric per-covariate descriptive tests (WeightedCluster style).
# No multinomial fit; no ORs. Produces a per-covariate table + a tidy
# standardized-residuals frame for the factor covariates so users can
# read which (cluster, level) cells drive each association.
#' @noRd
.covariate_chisq <- function(fit_df, profiles) {
  cov_names <- setdiff(names(fit_df), "cluster")
  cluster <- fit_df$cluster

  # One covariate at a time: its test row and, for a factor, its residuals.
  per_covariate <- lapply(cov_names, function(nm) {
    x <- fit_df[[nm]]
    if (is.factor(x)) {
      tab <- table(cluster = cluster, level = x)
      # Suppress chi-sq's small-expected-count warning; we report the
      # statistic regardless and let users read the residuals.
      ct <- suppressWarnings(stats::chisq.test(tab))
      n_total <- sum(tab)
      n_cells <- prod(dim(tab))
      # Cramer's V on a (r x c) table: sqrt(chi^2 / (n * min(r-1, c-1))).
      v <- sqrt(as.numeric(ct$statistic) /
                  (n_total * min(nrow(tab) - 1L, ncol(tab) - 1L)))
      test_row <- data.frame(
        variable = nm, type = "factor",
        statistic = unname(ct$statistic),
        df = unname(ct$parameter),
        p_value = ct$p.value,
        effect_size = v,
        effect_label = "Cramer's V",
        stringsAsFactors = FALSE
      )
      # stdres is the residual divided by its standard error -> approx N(0,1)
      # under H_0; |stdres| > 2 flags a cell with significant departure.
      sr <- suppressWarnings(ct$stdres)
      sr_df <- as.data.frame(as.table(sr), stringsAsFactors = FALSE)
      names(sr_df) <- c("cluster", "level", "std_residual")
      sr_df$variable <- nm
      sr_df$observed <- as.numeric(as.vector(tab))
      sr_df$expected <- as.numeric(as.vector(ct$expected))
      list(test = test_row,
           resid = sr_df[, c("variable", "cluster", "level", "observed",
                             "expected", "std_residual"),
                         drop = FALSE])
    } else if (is.numeric(x)) {
      kw <- suppressWarnings(stats::kruskal.test(x ~ cluster))
      # eta-squared for KW: H / (n - 1) bounded to [0, 1]; sometimes
      # reported as (H - g + 1) / (n - g). Use the simpler, more common
      # form so it matches `effectsize::rank_epsilon_squared`'s scale.
      n_total <- length(x)
      eta_sq <- as.numeric(kw$statistic) / (n_total - 1L)
      eta_sq <- max(0, min(1, eta_sq))
      list(test = data.frame(
        variable = nm, type = "numeric",
        statistic = unname(kw$statistic),
        df = unname(kw$parameter),
        p_value = kw$p.value,
        effect_size = eta_sq,
        effect_label = "eta^2 (KW)",
        stringsAsFactors = FALSE
      ), resid = NULL)
    } else {
      stop(errorCondition(paste0(sprintf(
        "Covariate '%s' has unsupported type '%s' for estimator = 'chisq'. ",
        nm, class(x)[1L]),
        "Convert to factor or numeric."),
        class = "nestimate_bad_input", call = NULL))
    }
  })
  test_rows <- lapply(per_covariate, `[[`, "test")
  resid_rows <- Filter(Negate(is.null), lapply(per_covariate, `[[`, "resid"))

  tests <- do.call(rbind, test_rows)
  if (!is.null(tests)) {
    tests$sig <- ifelse(is.na(tests$p_value), "",
      ifelse(tests$p_value < 0.001, "***",
      ifelse(tests$p_value < 0.01,  "**",
      ifelse(tests$p_value < 0.05,  "*",
      ifelse(tests$p_value < 0.1,   ".", "")))))
  }
  residuals_df <- if (length(resid_rows)) do.call(rbind, resid_rows) else NULL

  list(
    profiles = profiles,
    tests = tests,
    residuals = residuals_df,
    fit = list(method = "chisq",
               n_factor_tests = sum(vapply(cov_names,
                 function(n) is.factor(fit_df[[n]]), logical(1L))),
               n_numeric_tests = sum(vapply(cov_names,
                 function(n) is.numeric(fit_df[[n]]), logical(1L)))),
    model = NULL
  )
}


# Multinomial logit, dispatched on estimator. Firth via brglm2::brmultinom
# (default); ML via nnet::multinom. Output shape is shared so downstream
# print/summary code does not need to branch.
#' @noRd
.covariate_logit <- function(fit_df, rhs, profiles, estimator,
                              raw_assignments) {
  # If estimator was selected via auto (attribute set in caller), suppress
  # the multinom safety warning -- auto already verified no separation risk.
  via_auto <- isTRUE(identical(attr(profiles, "auto_choice"), "multinom"))
  if (estimator == "multinom") {
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop(errorCondition(paste0( # nocov start
        "Package 'nnet' is required for estimator = 'multinom'. ",
        "Install with: install.packages('nnet')"),
        class = "nestimate_missing_package", call = NULL)) # nocov end
    }
    if (!via_auto) warning(
      "estimator = 'multinom' uses nnet::multinom (classical ML). ",
      "Coefficient estimates can diverge under quasi-complete separation ",
      "(rare-cell or zero-cell covariate levels), producing astronomical ",
      "ORs with degenerate confidence intervals. The default ",
      "estimator = 'auto' falls back to firth (brglm2) when separation ",
      "risk is detected. See Heinze & Schemper (2002).",
      call. = FALSE
    )
    design <- stats::model.matrix(stats::as.formula(paste("~", rhs)), fit_df)
    fit <- .multinom_standardised(fit_df$cluster, design)
    null_fit <- nnet::multinom(cluster ~ 1, data = fit_df, trace = FALSE)
  } else {
    # estimator == "firth"
    if (!requireNamespace("brglm2", quietly = TRUE)) {
      stop(errorCondition(paste0( # nocov start
        "Package 'brglm2' is required for estimator = 'firth'. ",
        "Install with: install.packages('brglm2'). ",
        "Pass estimator = 'multinom' to use nnet::multinom instead, ",
        "or estimator = 'chisq' for descriptive-only output."),
        class = "nestimate_missing_package", call = NULL)) # nocov end
    }
    fml <- stats::as.formula(paste("cluster ~", rhs))
    # Not wrapped in suppressWarnings(): a convergence warning from the Firth
    # fit is information about the reported odds ratios.
    fit <- brglm2::brmultinom(fml, data = fit_df, ref = 1L)
    null_fit <- brglm2::brmultinom(cluster ~ 1, data = fit_df, ref = 1L)
  }

  estimates <- if (estimator == "multinom") fit$estimates else summary(fit)
  coefs <- estimates$coefficients
  ses <- estimates$standard.errors
  if (estimator == "multinom") fit <- fit$model
  if (is.null(ses)) {
    # brmultinom stores SE on the vcov diagonal
    ses <- matrix(sqrt(diag(stats::vcov(fit))),
                  nrow = nrow(coefs), byrow = FALSE,
                  dimnames = dimnames(coefs))
  }
  if (is.null(dim(coefs))) {
    coefs <- matrix(coefs, nrow = 1L,
                    dimnames = list(levels(fit_df$cluster)[2L],
                                     names(coefs)))
    ses <- matrix(ses, nrow = 1L,
                  dimnames = list(levels(fit_df$cluster)[2L],
                                   names(ses)))
  }
  n_params <- length(coefs)
  min_cl <- min(table(raw_assignments))
  if (min_cl < n_params) {
    warning(sprintf(
      "Smallest cluster has %d observations but model has %d parameters. ",
      min_cl, n_params),
      "Estimates may be unreliable.", call. = FALSE)
  }

  # One row per coefficient, clusters outermost.
  rows <- lapply(seq_len(nrow(coefs)), function(i) {
    lapply(seq_len(ncol(coefs)), function(j) {
      est <- coefs[i, j]; se <- ses[i, j]
      z_val <- est / se
      p_val <- 2 * (1 - stats::pnorm(abs(z_val)))
      or <- exp(est); ci_lo <- exp(est - 1.96 * se); ci_hi <- exp(est + 1.96 * se)
      sig <- if (is.na(p_val)) "" else if (p_val < 0.001) "***" else
        if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else
        if (p_val < 0.1) "." else ""
      data.frame(
        cluster = rownames(coefs)[i],
        variable = colnames(coefs)[j],
        estimate = est, std_error = se,
        odds_ratio = or, ci_lower = ci_lo, ci_upper = ci_hi,
        z = z_val, p = p_val, sig = sig,
        stringsAsFactors = FALSE, row.names = NULL
      )
    })
  })
  coef_df <- do.call(rbind, unlist(rows, recursive = FALSE))

  null_dev <- tryCatch(-2 * as.numeric(stats::logLik(null_fit)),
                       error = function(e) NA_real_)
  model_dev <- tryCatch(-2 * as.numeric(stats::logLik(fit)),
                        error = function(e) NA_real_)
  mcfadden_r2 <- if (is.na(null_dev) || is.na(model_dev)) NA_real_
                 else 1 - model_dev / null_dev

  list(
    profiles = profiles,
    coefficients = coef_df,
    fit = list(
      method = estimator,
      aic = tryCatch(stats::AIC(fit), error = function(e) NA_real_),
      bic = tryCatch(stats::BIC(fit), error = function(e) NA_real_),
      deviance = model_dev,
      mcfadden_r2 = mcfadden_r2,
      reference_cluster = levels(fit_df$cluster)[1L],
      converged = if (estimator == "multinom") fit$convergence == 0L else NA
    ),
    model = fit
  )
}

#' Centre and scale a design matrix, and the map back to its original scale
#'
#' nnet's optimiser and its numeric Hessian fail on columns far from unit
#' scale (a year, an income). Fitting on standardised columns keeps the same
#' model: the likelihood and the fitted values do not change.
#' @param X Design matrix with named columns.
#' @return List with `X` (standardised) and `back`, the matrix with
#'   `original = standardised %*% back` for coefficients in rows.
#' @noRd
.standardise_design <- function(X) {
  is_intercept <- colnames(X) == "(Intercept)"
  center <- colMeans(X)
  center[is_intercept] <- 0
  # Without an intercept, centring would change the model.
  if (!any(is_intercept)) center[] <- 0
  spread <- apply(X, 2L, stats::sd)
  spread[is_intercept | !is.finite(spread) | spread == 0] <- 1
  back <- diag(1 / spread, nrow = ncol(X))
  dimnames(back) <- list(colnames(X), colnames(X))
  if (any(is_intercept)) back[, is_intercept] <- back[, is_intercept] - center / spread
  list(X = sweep(sweep(X, 2L, center, "-"), 2L, spread, "/"), back = back)
}

#' Multinomial logit on standardised columns, reported on the original scale
#'
#' @param y Factor outcome.
#' @param design Design matrix of the original covariates.
#' @return List with `model` (the nnet fit, on the standardised columns) and
#'   `estimates` (list of `coefficients` and `standard.errors`, one row per
#'   non-reference level, on the original scale). Warns with class
#'   `nestimate_no_converge` when the optimiser stopped early.
#' @noRd
.multinom_standardised <- function(y, design) {
  std <- .standardise_design(design)
  X <- std$X
  fit <- nnet::multinom(y ~ X - 1, Hess = TRUE, maxit = 1000L, trace = FALSE)
  .warn_multinom_convergence(fit)
  levels_kept <- levels(y)[-1L]
  p <- ncol(design)
  beta <- matrix(stats::coef(fit), nrow = length(levels_kept), ncol = p,
                 dimnames = list(levels_kept, colnames(design)))
  covariance <- stats::vcov(fit)
  # vcov() holds one p x p block per non-reference level, in level order.
  ses <- vapply(seq_along(levels_kept), \(level) {
    block <- (level - 1L) * p + seq_len(p)
    sqrt(diag(t(std$back) %*% covariance[block, block, drop = FALSE] %*%
                std$back))
  }, numeric(p))
  ses <- matrix(ses, nrow = length(levels_kept), ncol = p, byrow = TRUE,
                dimnames = dimnames(beta))
  list(model = fit,
       estimates = list(coefficients = beta %*% std$back,
                        standard.errors = ses))
}

#' Warn when nnet::multinom stopped before converging
#' @param fit A fit with a `convergence` field (0 = converged).
#' @return `fit`, invisibly.
#' @noRd
.warn_multinom_convergence <- function(fit) {
  if (!identical(as.integer(fit$convergence), 0L)) {
    warning(warningCondition(
      paste0("nnet::multinom did not converge within its iteration limit, ",
             "so the covariate estimates are not reliable."),
      class = "nestimate_no_converge", call = NULL
    ))
  }
  invisible(fit)
}

#' Print covariate profiles and coefficient table (shared by build_clusters/mmm)
#' @noRd
.print_covariate_profiles <- function(cov, header = "Post-hoc Covariate Analysis (does not influence cluster membership)") {
  cat(header, "\n\n")

  # Numeric profiles
  if (!is.null(cov$profiles$numeric)) {
    cat("Cluster Profiles (numeric):\n")
    np <- cov$profiles$numeric
    vars <- unique(np$variable)
    rows <- lapply(unique(np$cluster), function(cl) {
      sub <- np[np$cluster == cl, , drop = FALSE]
      row <- data.frame(
        Cluster = cl,
        `N (%)` = sprintf("%d (%.0f%%)", sub$n[1L], sub$pct[1L]),
        stringsAsFactors = FALSE, check.names = FALSE
      )
      cells <- unlist(lapply(vars, function(v) {
        sv <- sub[sub$variable == v, ]
        stats::setNames(
          list(sprintf("%.2f (%.2f)", sv$mean, sv$sd),
               sprintf("%.2f", sv$median)),
          c(sprintf("%s Mean (SD)", v), sprintf("%s Median", v)))
      }), recursive = FALSE)
      row[names(cells)] <- cells
      row
    })
    print(do.call(rbind, rows), row.names = FALSE, right = FALSE)
    cat("\n")
  }

  # Categorical profiles
  if (!is.null(cov$profiles$categorical)) {
    cat("Cluster Profiles (categorical):\n")
    cp <- cov$profiles$categorical
    # Pre-compute all variable=level combos for consistent columns
    all_cols <- unique(paste0(cp$variable, "=", cp$level))
    rows <- lapply(unique(cp$cluster), function(cl) {
      sub <- cp[cp$cluster == cl, , drop = FALSE]
      row <- data.frame(
        Cluster = cl, N = sub$n[1L],
        stringsAsFactors = FALSE, check.names = FALSE
      )
      cells <- lapply(all_cols, function(col) {
        parts <- strsplit(col, "=", fixed = TRUE)[[1L]]
        v <- parts[1L]
        lev <- paste(parts[-1L], collapse = "=")  # handle = in level names
        slev <- sub[sub$variable == v & sub$level == lev, ]
        if (nrow(slev) > 0L) {
          sprintf("%d (%.0f%%)", slev$count, slev$pct)
        } else {
          "0 (0%)"
        }
      })
      row[sprintf("%s N(%%)", all_cols)] <- cells
      row
    })
    print(do.call(rbind, rows), row.names = FALSE, right = FALSE)
    cat("\n")
  }

  # Estimator-specific section: ORs (firth/multinom) or test table (chisq).
  if (identical(cov$fit$method, "chisq")) {
    cat("Per-Covariate Cluster Association Tests:\n")
    tt <- cov$tests
    disp <- data.frame(
      Variable = tt$variable,
      Type = tt$type,
      Statistic = sprintf("%.3f", tt$statistic),
      df = ifelse(is.na(tt$df), "", as.character(tt$df)),
      p = ifelse(tt$p_value < 0.001, "<0.001",
                 sprintf("%.3f", tt$p_value)),
      Effect = paste0(sprintf("%.3f", tt$effect_size), " (", tt$effect_label, ")"),
      Sig = tt$sig,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    print(disp, row.names = FALSE, right = FALSE)

    if (!is.null(cov$residuals) && nrow(cov$residuals) > 0L) {
      cat("\nStandardized Adjusted Residuals (|r| > 2 flags significant cells):\n")
      rs <- cov$residuals
      rs$std_residual <- sprintf("%.2f", rs$std_residual)
      rs$observed <- as.integer(rs$observed)
      rs$expected <- sprintf("%.1f", rs$expected)
      print(rs, row.names = FALSE, right = FALSE)
    }
    cat("\nMethod: chi-square (categorical) / Kruskal-Wallis (numeric); ",
        "descriptive only.\n", sep = "")
  } else {
    cat(sprintf("Predictors of Membership (estimator = %s, reference: Cluster %s):\n",
                cov$fit$method %||% "multinom",
                cov$fit$reference_cluster))
    coef_display <- cov$coefficients[
      cov$coefficients$variable != "(Intercept)", , drop = FALSE]
    disp <- data.frame(
      Cluster = coef_display$cluster,
      Variable = coef_display$variable,
      OR = sprintf("%.2f", coef_display$odds_ratio),
      `95% CI` = sprintf("[%.2f, %.2f]", coef_display$ci_lower,
                          coef_display$ci_upper),
      p = ifelse(coef_display$p < 0.001, "<0.001",
                 sprintf("%.3f", coef_display$p)),
      Sig = coef_display$sig,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    print(disp, row.names = FALSE, right = FALSE)
    cat(sprintf("\nModel: AIC = %.1f | BIC = %.1f | McFadden R-squared = %.2f\n",
                cov$fit$aic, cov$fit$bic, cov$fit$mcfadden_r2))
  }
}

# Tidy a covariate-analysis result list into a single data.frame the user
# can hand to write.csv / kable / dplyr without further reshaping. Returns
# different shapes per estimator (logit vs chisq) - the column set names
# the relevant statistics for each path so the meaning is unambiguous.
#
# Attributes carry the metadata that didn't fit a tidy column (cluster
# sizes, fit stats, residuals frame for chisq). `print.tidy_covariates`
# renders the same human-readable view that `.print_covariate_profiles()`
# does, so the tidy frame remains REPL-friendly.
#' @noRd
.tidy_covariates <- function(cov, cluster_stats = NULL) {
  if (is.null(cov)) return(NULL)
  method <- cov$fit$method %||% "multinom"

  if (method == "chisq") {
    out <- cov$tests
    if (!is.null(out)) {
      out <- out[, c("variable", "type", "statistic", "df", "p_value",
                     "effect_size", "effect_label", "sig"),
                  drop = FALSE]
      names(out)[names(out) == "p_value"] <- "p"
      names(out)[names(out) == "effect_label"] <- "effect_kind"
      names(out)[names(out) == "effect_size"] <- "effect"
      rownames(out) <- NULL
    }
  } else {
    co <- cov$coefficients
    co <- co[co$variable != "(Intercept)", , drop = FALSE]
    out <- data.frame(
      cluster   = co$cluster,
      predictor = co$variable,
      OR        = co$odds_ratio,
      CI_low    = co$ci_lower,
      CI_high   = co$ci_upper,
      p         = co$p,
      sig       = co$sig,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  attr(out, "estimator")       <- method
  attr(out, "fit")             <- cov$fit
  attr(out, "profiles")        <- cov$profiles
  attr(out, "residuals")       <- cov$residuals
  attr(out, "cluster_stats")   <- cluster_stats
  attr(out, "covariates_full") <- cov
  class(out) <- c("tidy_covariates", "data.frame")
  out
}

#' @rdname build_clusters
#' @export
print.tidy_covariates <- function(x, ...) {
  cat(sprintf("# Tidy covariate analysis (estimator = %s)\n",
              attr(x, "estimator") %||% "unknown"))
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}


#' Odds ratio forest plot (shared by build_clusters/mmm)
#' @noRd
.plot_covariate_forest <- function(coef_df, title, combined = TRUE) {
  coef_df <- coef_df[coef_df$variable != "(Intercept)", , drop = FALSE]
  coef_df$cluster <- factor(coef_df$cluster)
  coef_df$significant <- coef_df$p < 0.05

  base_p <- function(d, ttl) {
    ggplot(d, aes(
      x = .data$odds_ratio, y = .data$variable,
      xmin = .data$ci_lower, xmax = .data$ci_upper,
      colour = .data$significant
    )) +
      geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
      geom_pointrange(size = 0.6) +
      scale_colour_manual(
        values = c("TRUE" = "#D62728", "FALSE" = "#7F7F7F"),
        labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
        name = NULL,
        drop = FALSE
      ) +
      scale_x_log10() +
      labs(x = "Odds Ratio (log scale)", y = NULL, title = ttl) +
      theme_minimal(base_size = 12)
  }

  if (!combined) {
    clusters <- levels(coef_df$cluster)
    plots <- lapply(clusters, function(cl) {
      sub <- coef_df[coef_df$cluster == cl, , drop = FALSE]
      base_p(sub, sprintf("%s -- Cluster %s", title, cl))
    })
    names(plots) <- clusters
    return(invisible(plots))
  }

  p <- base_p(coef_df, title) +
    facet_wrap(~ .data$cluster, labeller = label_both)
  print(p)
  invisible(p)
}

# ==============================================================================
# 8. build_network dispatch for net_clustering
# ==============================================================================

#' Build per-cluster networks from clustering results
#'
#' Called by \code{build_network()} when data is a \code{net_clustering}
#' object. Splits sequence data by cluster assignment and builds a network
#' for each cluster.
#'
#' @param x A \code{net_clustering} object.
#' @param method Network estimation method (default: "relative").
#' @param ... Passed to \code{build_network()}.
#' @return A \code{netobject_group}. When \code{data} is an HTNA network,
#'   the result is an \code{htna_group}; every child retains the original
#'   node-to-actor partition and HTNA class, while the clustering diagnostics
#'   remain attached to the outer object.
#' @noRd
.build_network_clustering <- function(x, method = "relative", ...) {
  assignments <- x$assignments
  seq_data <- x$data
  k <- x$k

  # Merge stored build_args with caller's ...; caller takes precedence.
  # Drop long-format-only column-name args (actor/time/session/order/action)
  # because x$data is already in wide format (T1..Tn columns) at this point,
  # so forwarding e.g. actor = "session_id" would fail validation.
  long_only <- c("actor", "time", "session", "order", "action",
                 "state_cols", "metadata_cols")
  dots <- list(...)
  build_args <- if (!is.null(x$build_args))
    modifyList(x$build_args, dots) else dots
  build_args <- build_args[setdiff(names(build_args), long_only)]

  nets <- lapply(seq_len(k), function(cl) {
    rows <- assignments == cl
    sub <- seq_data[rows, , drop = FALSE]
    net <- do.call(build_network, c(list(data = sub, method = method),
                                    build_args))
    .attach_set_aside_metadata(net, x$metadata, rows)
  })

  names(nets) <- paste("Cluster", seq_len(k))
  attr(nets, "group_col") <- "cluster"
  attr(nets, "clustering") <- x
  class(nets) <- "netobject_group"
  .restore_htna_group(nets, x$htna_partition)
}

# A network's estimator settings live in `$params` (an attention network's
# `lambda`, an n-gram's `n_gram`), not in `$build_args`; without them the
# per-cluster rebuild silently fell back to the defaults. The caller's own
# `params`, and a method other than the network's, take precedence.
#' @noRd
.inherit_estimator_params <- function(build_dots, net) {
  same_method <- identical(.resolve_method_alias(build_dots$method),
                           .resolve_method_alias(net$method))
  if (is.null(build_dots$params) && same_method && length(net$params) > 0L) {
    build_dots$params <- net$params
  }
  build_dots
}

# Add the columns set aside at fitting time to one component's network. The
# rebuild can find metadata of its own (a position that is empty in this
# component), so it joins what is there instead of giving way to it.
#' @noRd
.attach_set_aside_metadata <- function(net, metadata, rows) {
  if (is.null(metadata)) return(net)
  new_cols <- setdiff(names(metadata), names(net$metadata))
  if (length(new_cols) == 0L) return(net)
  extra <- metadata[rows, new_cols, drop = FALSE]
  net$metadata <- if (is.null(net$metadata)) extra else
    cbind(net$metadata, extra)
  net
}

# ==============================================================================
# 9. cluster_network - one-shot clustering + network estimation
# ==============================================================================

#' Cluster data and build per-cluster networks in one step
#'
#' @description
#' `cluster_network()` clusters a set of sequences and estimates one network
#' per cluster in a single call. It runs [build_clusters()] (dissimilarity
#' clustering with PAM or a hierarchical method) or, with
#' `cluster_by = "mmm"`, [build_mmm()] (a mixture of Markov chains), then
#' [build_network()] on each cluster's sequences. The network method is
#' the `method` given in `...`; without one, a `netobject` input passes on
#' its own method and stored build arguments, so the per-cluster networks
#' are of the same type as the input, and other input gives `"relative"`
#' transition networks.
#'
#' The clusters were formed to make the sequences within each one alike, so
#' their networks differ partly by construction. Check the classification
#' with [cluster_diagnostics()] and the choice of `k` with
#' [cluster_choice()] or [compare_mmm()], and see how far each cluster's
#' edges move under resampling with [bootstrap_network()] before reading
#' differences between the networks.
#'
#' @param data The sequences: a data frame, matrix, `netobject`, bare
#'   `cograph_network` or `tna` model, as accepted by [build_clusters()] (or
#'   by [build_mmm()] when `cluster_by = "mmm"`). A `wtna` network is
#'   accepted only with `cluster_by = "mmm"`; dissimilarity clustering
#'   refuses it with an error of class `nestimate_bad_input`.
#' @param k The number of clusters, a whole number from 2 to one less than
#'   the number of sequences. Required.
#' @param cluster_by The clustering algorithm. `"pam"` (default) and the
#'   hierarchical methods `"ward.D2"`, `"ward.D"`, `"complete"`,
#'   `"average"`, `"single"`, `"mcquitty"`, `"median"` and `"centroid"` are
#'   passed to [build_clusters()] as its `method`; `"mmm"` fits a mixed
#'   Markov model with [build_mmm()] instead.
#' @param dissimilarity The distance between sequences for dissimilarity
#'   clustering, one of those [build_clusters()] accepts. Default
#'   `"hamming"`. Giving it with `cluster_by = "mmm"` raises an error.
#' @param ... Arguments for the two stages, routed by name. For
#'   dissimilarity clustering, [build_clusters()]'s `na_syms`, `weighted`,
#'   `lambda`, `seed`, `q`, `p`, `covariates` and `estimator` go to the
#'   clustering; with `cluster_by = "mmm"`, [build_mmm()]'s `n_starts`,
#'   `max_iter`, `tol`, `smooth`, `seed`, `covariates`, `covariate_effect`
#'   and `estimator` go to the mixture fit. Everything else is meant for the
#'   network stage and must be an argument of [build_network()] or of a
#'   registered estimator; any other name raises an error before either
#'   stage runs. A name both stages accept goes to the clustering: under
#'   dissimilarity clustering a `lambda` given here is the Hamming decay of
#'   the clustering and never reaches an attention network, and the same
#'   holds for `weighted`. When `data` is a `netobject`, its stored build
#'   arguments are reused for the per-cluster networks, and a value given in
#'   the call replaces the stored one. The *Details* section lists which
#'   network arguments take effect.
#' @details
#' `method`, the estimators' own arguments (such as `n_gram` or
#' `direction`) and [build_network()]'s `threshold`, `scaling`, `start`,
#' `end` and `params` reach the per-cluster networks. A `netobject` input
#' of the same method also passes on its estimator settings (`$params`, such
#' as an attention network's `lambda`) unless `params` is given. With
#' `cluster_by = "mmm"` and the default relative networks, which are the
#' fitted chains, `scaling` and `threshold` are applied to the fitted
#' weights, and an estimator argument, `start`, `end` or a new `params`
#' raises an error of class `nestimate_bad_input`. For a one-hot `netobject`
#' the stored build arguments are not used. Under `"mmm"`, rebuilding a
#' network built from
#' long-format data with any method other than `"relative"` currently raises
#' an error, because its column names are passed on to the wide per-cluster
#' data.
#' @return A `netobject_group`: one `netobject` per cluster, named
#'   `"Cluster 1"`, `"Cluster 2"`, and so on, with the clustering (a
#'   `net_clustering`, or for `"mmm"` a `net_mmm_clustering`) kept as
#'   its `"clustering"` attribute. `print()` names the clustering method and
#'   shows each cluster's network size, weight range and share of the
#'   sequences; [cluster_diagnostics()] reads the clustering, and
#'   [summary()] (the networks' metrics, a column per cluster),
#'   [bootstrap_network()] and [net_centrality()] take the group as it is.
#' @seealso \code{\link{build_clusters}}, \code{\link{cluster_mmm}},
#'   \code{\link{build_network}}
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
#'                    V3 = c("C","A","B","C","B"))
#' grp <- cluster_network(seqs, k = 2)
#' grp
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 50, TRUE), V2 = sample(LETTERS[1:4], 50, TRUE),
#'   V3 = sample(LETTERS[1:4], 50, TRUE), V4 = sample(LETTERS[1:4], 50, TRUE)
#' )
#' # Default: PAM clustering, relative (transition) networks
#' grp <- cluster_network(seqs, k = 3)
#'
#'
#' # MMM-based clustering
#' grp <- cluster_network(seqs, k = 2, cluster_by = "mmm")
#' }
#' @export
cluster_network <- function(data, k, cluster_by = "pam",
                             dissimilarity = "hamming", ...) {
  if (missing(k)) {
    stop(errorCondition(
      "argument \"k\" is missing, with no default",
      class = "nestimate_bad_input", call = NULL))
  }
  caller_dots <- list(...)

  # Coerce cograph_network to netobject so the build_args lift below works.
  if (inherits(data, "cograph_network") && !inherits(data, "netobject")) {
    data <- .as_netobject(data)
  }

  if (identical(cluster_by, "mmm")) {
    if (!missing(dissimilarity)) {
      stop(errorCondition(paste0(
        "`dissimilarity` is only valid for distance clustering. ",
        "Omit it when `cluster_by = \"mmm\"`."),
        class = "nestimate_bad_input", call = NULL)
      )
    }

    # Split caller's `...` into MMM EM hyperparameters vs build_network args
    # so a single cluster_network() call controls both stages without a
    # silently-dropped `n_starts`/`max_iter`/`seed`. Anything else falls
    # through to build_network() unchanged.
    # Taken from build_mmm()'s own formals so a new build_mmm() argument
    # (covariate_effect was missed) is routed without editing this list.
    mmm_arg_names  <- setdiff(names(formals(build_mmm)), c("data", "k"))
    mmm_args   <- caller_dots[intersect(names(caller_dots), mmm_arg_names)]
    build_dots <- caller_dots[setdiff(names(caller_dots), mmm_arg_names)]
    # A one-hot network's units are fitted as windowed state sequences, so the
    # fit carries the method that rebuilds them (build_mmm()); its own wtna
    # method and one-hot build arguments would read those as 0/1 codes.
    onehot <- inherits(data, "netobject") && .is_onehot_params(data$params)
    if (inherits(data, "netobject") && !onehot) {
      if (!is.null(data$build_args))
        build_dots <- modifyList(data$build_args, build_dots)
      if (is.null(build_dots$method))
        build_dots$method <- data$method
    }
    .validate_cluster_network_build_args(build_dots)
    mmm_fit <- do.call(build_mmm, c(list(data = data, k = k), mmm_args))
    if (onehot && is.null(build_dots$method)) {
      build_dots$method <- mmm_fit$network_method %||% "relative"
    }
    return(do.call(build_network, c(list(data = mmm_fit), build_dots)))
  }

  # Split caller's `...` into build_clusters() args vs build_network() args.
  # Without this split, weighted/lambda/seed/q/p/na_syms/covariates are
  # silently dropped because they reach build_network() instead of the
  # clusterer (audit_clustering finding #1). Splitting on caller_dots only -
  # never on data$build_args - keeps network-history args (e.g. an `atna`
  # net's `lambda` decay) routed to build_network where they belong.
  cluster_arg_names <- c("na_syms", "weighted", "lambda", "seed",
                         "q", "p", "covariates", "estimator")
  cluster_args <- caller_dots[intersect(names(caller_dots), cluster_arg_names)]
  build_dots   <- caller_dots[setdiff(names(caller_dots), cluster_arg_names)]

  if (inherits(data, "netobject")) {
    if (!is.null(data$build_args))
      build_dots <- modifyList(data$build_args, build_dots)
    if (is.null(build_dots$method))
      build_dots$method <- data$method
    build_dots <- .inherit_estimator_params(build_dots, data)
  }
  .validate_cluster_network_build_args(build_dots)

  cls <- do.call(build_clusters,
                 c(list(data = data, k = k, method = cluster_by,
                        dissimilarity = dissimilarity),
                   cluster_args))
  do.call(build_network, c(list(data = cls), build_dots))
}

.validate_cluster_network_build_args <- function(args) {
  arg_names <- names(args)
  if (is.null(arg_names) || !length(arg_names)) return(invisible(TRUE))

  build_formals <- setdiff(names(formals(build_network)), c("data", "..."))
  estimator_names <- list_estimators()$name
  estimator_formals <- unique(unlist(lapply(estimator_names, function(nm) {
    entry <- get_estimator(nm)
    if (is.null(entry$fn)) return(character())
    setdiff(names(formals(entry$fn)), c("data", "..."))
  }), use.names = FALSE))
  accepted <- unique(c(build_formals, estimator_formals))

  unknown <- setdiff(arg_names, accepted)
  if (length(unknown)) {
    stop(errorCondition(paste0(
      "Unknown argument(s) in cluster_network(): ",
      paste(unknown, collapse = ", "), ". ",
      "Pass estimator-specific options through `params = list(...)` ",
      "when using custom estimators."),
      class = "nestimate_bad_input", call = NULL)
    )
  }

  invisible(TRUE)
}
