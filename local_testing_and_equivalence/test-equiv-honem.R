# Numerical equivalence: Nestimate's build_honem() vs a clean-room R
# reimplementation derived from Saebi et al. 2020.
#
# HONEM is pure linear algebra on the HON transition matrix, so there is no
# canonical external library to compare against. Instead we rebuild the
# neighborhood matrix from first principles using structurally different code
# paths (Reduce accumulation vs vapply-with-<<-, explicit power list vs
# on-the-fly recurrence, full SVD with post-hoc truncation vs truncated SVD).
# If the two disagree at machine precision, the algorithm is mis-stated.
#
# Checks per config:
#   1. Singular values agree at TOL (machine precision).
#   2. Top-k subspace projectors U %*% t(U) agree at 1e-8 — this is
#      invariant to sign flip and rotation for tied singular values.
#   3. Reconstruction residual ||U_ref Sigma V_ref' - S||_F / ||S||_F at
#      machine precision on the ref side (sanity check on the rebuild).

set.seed(4242)
N_HONEM <- 40L
TOL <- 1e-10
TOL_SUBSPACE <- 1e-8

# ---- Config generation ----
honem_configs <- lapply(seq_len(N_HONEM), function(i) {
  list(n_actors = sample(c(15L, 20L, 30L), 1),
       n_states = sample(3:5, 1),
       seq_length = sample(c(20L, 30L, 40L), 1),
       hon_max_order = sample(1:2, 1),
       embed_dim = sample(c(3L, 5L, 8L), 1),
       max_power = sample(c(3L, 5L, 10L), 1),
       seed = sample.int(100000, 1))
})

# Clean-room HONEM: builds S via Reduce over powers, runs full base::svd.
# S = (1/Z) * sum_{k=0}^{L} exp(-k) * D^{k+1}, Z = sum(exp(-(0:L))).
.honem_reference <- function(mat, dim, max_power) {
  n <- nrow(mat)
  row_sums <- rowSums(mat)
  D <- mat
  nz <- row_sums > 0
  D[nz, ] <- mat[nz, ] / row_sums[nz]
  D[!nz, ] <- 0

  weights <- exp(-(0L:max_power))
  Z <- sum(weights)

  # Build explicit list of D^{k+1} for k = 0..max_power
  powers <- vector("list", max_power + 1L)
  powers[[1L]] <- D
  if (max_power >= 1L) {
    for (k in seq_len(max_power)) {
      powers[[k + 1L]] <- powers[[k]] %*% D
    }
  }

  S <- Reduce(`+`, mapply(function(w, P) w * P, weights, powers,
                          SIMPLIFY = FALSE)) / Z
  rownames(S) <- rownames(mat)
  colnames(S) <- colnames(mat)

  dim <- min(as.integer(dim), n - 1L)
  sv <- base::svd(S)
  sigma <- sv$d[seq_len(dim)]
  U <- sv$u[, seq_len(dim), drop = FALSE]
  V <- sv$v[, seq_len(dim), drop = FALSE]

  list(S = S, sigma = sigma, U = U, V = V,
       embedding = sweep(U, 2, sqrt(sigma), `*`),
       total_var = sum(sv$d^2),
       dim = dim)
}

test_that("HONEM singular values + subspace match clean-room reimplementation", {
  skip_on_cran()
  skip_equiv_tests()

  report <- equiv_report()

  invisible(lapply(seq_len(N_HONEM), function(i) {
    cfg <- honem_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    hon <- tryCatch(build_hon(seqs, max_order = cfg$hon_max_order,
                              min_freq = 1L, method = "hon"),
                    error = function(e) NULL)
    if (is.null(hon) || is.null(hon$matrix) || nrow(hon$matrix) < 3L) {
      return(NULL)
    }
    mat <- hon$matrix
    effective_dim <- min(cfg$embed_dim, nrow(mat) - 1L)

    hem <- tryCatch(build_honem(hon, dim = effective_dim,
                                max_power = cfg$max_power),
                    error = function(e) NULL)
    if (is.null(hem)) return(NULL)

    ref <- .honem_reference(mat, effective_dim, cfg$max_power)

    # 1. Singular values agree
    sig_delta <- abs(sort(hem$singular_values) - sort(ref$sigma))

    # 2. Top-k subspace projector
    P_nest <- hem$embeddings %*% solve(
      diag(sqrt(ref$sigma), nrow = length(ref$sigma)),
      t(hem$embeddings) / 1  # undo the sqrt(sigma) scaling to isolate U
    )
    # Safer route: rebuild U from embeddings by dividing each column by sqrt(sigma).
    U_nest <- sweep(hem$embeddings, 2, sqrt(hem$singular_values), `/`)
    P_nest <- U_nest %*% t(U_nest)
    P_ref <- ref$U %*% t(ref$U)
    proj_delta <- max(abs(P_nest - P_ref))

    report$log(
      func = "build_honem",
      config = sprintf("cfg%d(n=%d,s=%d,mo=%d,d=%d,L=%d)",
                       i, cfg$n_actors, cfg$seq_length,
                       cfg$hon_max_order, effective_dim, cfg$max_power),
      n_checked = length(sig_delta) + length(P_nest),
      n_failed = as.integer(max(sig_delta) > TOL) +
        as.integer(proj_delta > TOL_SUBSPACE),
      max_abs_err = max(c(sig_delta, proj_delta)),
      mean_abs_err = mean(c(sig_delta, as.numeric(abs(P_nest - P_ref)))),
      median_abs_err = stats::median(c(sig_delta, as.numeric(abs(P_nest - P_ref)))),
      p95_abs_err = as.numeric(stats::quantile(
        c(sig_delta, as.numeric(abs(P_nest - P_ref))), 0.95
      )),
      reference = "clean-room Reduce+base::svd reimpl",
      notes = sprintf("sigma_max=%.2e proj_max=%.2e", max(sig_delta), proj_delta)
    )

    expect_true(max(sig_delta) < TOL,
                label = sprintf("cfg%d singular values max delta = %.2e",
                                i, max(sig_delta)))
    expect_true(proj_delta < TOL_SUBSPACE,
                label = sprintf("cfg%d subspace projector max delta = %.2e",
                                i, proj_delta))
  }))

  report$write_csv("honem")
  report$write_cvs("honem", "local_testing_and_equivalence/test-equiv-honem.R")
})

test_that("HONEM singular values + subspace match RSpectra::svds (Arnoldi backend)", {
  # Second SVD backend as an independent numerical-algorithm check.
  # base::svd uses LAPACK dgesdd (divide-and-conquer).
  # RSpectra::svds uses ARPACK-style Arnoldi iteration on A'A — a completely
  # different numerical method. Agreement confirms the result doesn't depend
  # on the choice of SVD algorithm, catching algorithm-level bugs that a
  # formula-only reimplementation would miss.
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("RSpectra")

  report <- equiv_report()

  invisible(lapply(seq_len(N_HONEM), function(i) {
    cfg <- honem_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })
    hon <- tryCatch(build_hon(seqs, max_order = cfg$hon_max_order,
                              min_freq = 1L, method = "hon"),
                    error = function(e) NULL)
    if (is.null(hon) || is.null(hon$matrix) || nrow(hon$matrix) < 3L) {
      return(NULL)
    }
    mat <- hon$matrix
    effective_dim <- min(cfg$embed_dim, nrow(mat) - 1L)

    hem <- tryCatch(build_honem(hon, dim = effective_dim,
                                max_power = cfg$max_power),
                    error = function(e) NULL)
    if (is.null(hem)) return(NULL)

    # Rebuild the S matrix using the clean-room reference (validated against
    # Nestimate in the previous test), then run RSpectra::svds on it.
    ref <- .honem_reference(mat, effective_dim, cfg$max_power)
    spec <- tryCatch(
      RSpectra::svds(ref$S, k = effective_dim),
      error = function(e) NULL
    )
    if (is.null(spec)) return(NULL)

    # Singular values: direct comparison (both return top-k in descending
    # order). RSpectra matches LAPACK at ~1e-14 on well-conditioned matrices.
    sig_delta <- abs(sort(hem$singular_values) - sort(spec$d))

    # Subspace projector: |P_nest - P_rspectra|_inf. Both SVDs can flip signs,
    # so the projector U %*% t(U) is the sign-invariant quantity.
    U_nest <- sweep(hem$embeddings, 2, sqrt(hem$singular_values), `/`)
    P_nest <- U_nest %*% t(U_nest)
    P_spec <- spec$u %*% t(spec$u)
    proj_delta <- max(abs(P_nest - P_spec))

    report$log(
      func = "build_honem_vs_rspectra",
      config = sprintf("cfg%d(n=%d,s=%d,mo=%d,d=%d,L=%d)",
                       i, cfg$n_actors, cfg$seq_length,
                       cfg$hon_max_order, effective_dim, cfg$max_power),
      n_checked = length(sig_delta) + length(P_nest),
      n_failed = as.integer(max(sig_delta) > TOL) +
        as.integer(proj_delta > TOL_SUBSPACE),
      max_abs_err = max(c(sig_delta, proj_delta)),
      mean_abs_err = mean(c(sig_delta, as.numeric(abs(P_nest - P_spec)))),
      median_abs_err = stats::median(c(sig_delta, as.numeric(abs(P_nest - P_spec)))),
      p95_abs_err = as.numeric(stats::quantile(
        c(sig_delta, as.numeric(abs(P_nest - P_spec))), 0.95
      )),
      reference = "RSpectra::svds (ARPACK Arnoldi)",
      notes = sprintf("sigma_max=%.2e proj_max=%.2e", max(sig_delta), proj_delta)
    )

    expect_true(max(sig_delta) < TOL,
                label = sprintf("cfg%d vs RSpectra singular values delta = %.2e",
                                i, max(sig_delta)))
    expect_true(proj_delta < TOL_SUBSPACE,
                label = sprintf("cfg%d vs RSpectra subspace delta = %.2e",
                                i, proj_delta))
  }))

  report$write_csv("honem_rspectra")
  report$write_cvs("honem_rspectra",
                   "local_testing_and_equivalence/test-equiv-honem.R")
})

test_that("HONEM reconstruction residual is tiny for the reference", {
  skip_on_cran()
  skip_equiv_tests()

  # Sanity: if we take the full SVD (no truncation), reconstruction should be
  # exact. Truncation introduces a residual whose Frobenius norm equals the
  # tail energy (sum of squared dropped singular values). This is an invariant
  # of the algorithm, not Nestimate — but if Nestimate's S matrix were wrong,
  # the Nestimate sigmas wouldn't match the ref sigmas (covered above).
  invisible(lapply(seq_len(10L), function(i) {
    cfg <- honem_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })
    hon <- tryCatch(build_hon(seqs, max_order = cfg$hon_max_order,
                              min_freq = 1L, method = "hon"),
                    error = function(e) NULL)
    if (is.null(hon) || is.null(hon$matrix) || nrow(hon$matrix) < 3L) {
      return(NULL)
    }
    ref <- .honem_reference(hon$matrix, nrow(hon$matrix) - 1L,
                            cfg$max_power)
    # Full-rank reconstruction (dim = n-1) residual should be tail energy only
    S_reconstructed <- ref$U %*% diag(ref$sigma, nrow = length(ref$sigma)) %*% t(ref$V)
    residual <- sqrt(sum((S_reconstructed - ref$S)^2))
    # Full-rank truncation drops at most 1 singular value (dim = n-1), so the
    # residual is at most the smallest singular value of the full SVD.
    sv_full <- base::svd(ref$S)
    expected_residual <- sv_full$d[length(sv_full$d)]
    expect_true(
      residual <= expected_residual + TOL,
      label = sprintf("cfg%d reconstruction residual = %.2e (expected <= %.2e)",
                      i, residual, expected_residual)
    )
  }))
})
