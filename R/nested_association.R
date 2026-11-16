# ---- Association tests for nested (clustered) contingency tables ----
#
# Engine behind the clustered residuals of mosaic_plot() and the clustered
# tests of mosaic_analysis(). Every method reads one R x C x K count array:
# the contingency table of each of K clusters (sequences, actors, teams).
# Pooling the K tables and running Pearson's chi-square treats every count as
# an independent draw; the methods below correct for the dependence of counts
# that share a cluster. Each returns per-cell residuals (drawn as the mosaic
# fill) and one global test.
#
#   rao_scott   design-based second-order Rao-Scott F (survey::svychisq's
#               default), residuals deflated by each cell's design effect
#   donner_klar Pearson X2 divided by the cluster-size inflation factor
#               1 + (m - 1) rho; one margin must be constant within clusters
#   cmh         generalised Cochran-Mantel-Haenszel, clusters as strata
#   gee         marginal model of the response given the predictor
#               (identity link, independence working correlation), Wald X2
#               with the Liang-Zeger cluster-robust covariance
#   glmm        log-linear Poisson mixed model with cluster random effects
#               (lme4), Wald X2 on the row-by-column interaction


.nested_methods <- c("rao_scott", "donner_klar", "cmh", "gee", "glmm")

.nested_labels <- c(
  rao_scott   = "Rao-Scott",
  donner_klar = "Donner-Klar",
  cmh         = "Cochran-Mantel-Haenszel",
  gee         = "GEE Wald",
  glmm        = "GLMM Wald"
)


#' Clustered association test and residuals for an R x C x K count array
#'
#' @param counts Numeric array, rows x columns x clusters, whole numbers.
#' @param method One of `.nested_methods`.
#' @param levels For `glmm`: a data frame with one row per cluster (third
#'   dimension) and one column per nesting level, outermost first. `NULL`
#'   means one level.
#' @param predictor `"rows"` or `"cols"`: which margin `gee` and `glmm`
#'   condition on (the other is the response).
#' @return A list: `residuals` (R x C matrix, `NA` where undefined; empty
#'   rows and columns are `NA`) and `test` (one-row data frame).
#' @noRd
.nested_association <- function(counts, method = .nested_methods,
                                levels = NULL,
                                predictor = c("rows", "cols")) {
  method    <- match.arg(method)
  predictor <- match.arg(predictor)
  .assert(
    "`counts` must be a numeric array of three dimensions" =
      is.array(counts) && is.numeric(counts) && length(dim(counts)) == 3L,
    "`counts` must hold finite, non-negative whole numbers" =
      all(is.finite(counts)) && all(counts >= 0) &&
      all(abs(counts - round(counts)) <= 1e-8)
  )
  d <- dim(counts)
  dimnames(counts) <- list(
    dimnames(counts)[[1L]] %||% as.character(seq_len(d[1L])),
    dimnames(counts)[[2L]] %||% as.character(seq_len(d[2L])),
    dimnames(counts)[[3L]] %||% as.character(seq_len(d[3L])))
  if (is.null(levels)) {
    levels <- data.frame(cluster = dimnames(counts)[[3L]],
                         stringsAsFactors = FALSE)
  }
  .assert("`levels` must have one row per cluster" =
            is.data.frame(levels) && nrow(levels) == d[3L] &&
            ncol(levels) >= 1L)

  tab <- rowSums(counts, dims = 2L)
  keep_r <- rowSums(tab) > 0
  keep_c <- colSums(tab) > 0
  if (sum(keep_r) < 2L || sum(keep_c) < 2L) {
    stop(errorCondition(
      "A clustered test needs at least two non-empty rows and two non-empty columns.",
      class = "nestimate_bad_input", call = NULL))
  }
  sub <- counts[keep_r, keep_c, , drop = FALSE]
  keep_k <- apply(sub, 3L, sum) > 0
  sub <- sub[, , keep_k, drop = FALSE]
  levels <- levels[keep_k, , drop = FALSE]
  if (dim(sub)[3L] < 2L) {
    stop(errorCondition(
      sprintf("A clustered test needs at least two clusters with counts; found %d.",
              dim(sub)[3L]),
      class = "nestimate_bad_input", call = NULL))
  }

  # Only the GLMM models nesting levels; the other methods use the outermost
  # level, the ultimate cluster, which carries the dependence of every level
  # below it.
  if (!identical(method, "glmm") && ncol(levels) > 1L) {
    outer <- as.character(levels[[1L]])
    flat <- rowsum(t(matrix(sub, prod(dim(sub)[1:2]), dim(sub)[3L])), outer,
                   reorder = FALSE)
    sub <- array(t(flat), c(dim(sub)[1:2], nrow(flat)),
                 dimnames = c(dimnames(sub)[1:2], list(rownames(flat))))
    levels <- data.frame(cluster = rownames(flat), stringsAsFactors = FALSE)
  }

  flip <- identical(predictor, "cols") && method %in% c("gee", "glmm")
  work <- if (flip) aperm(sub, c(2L, 1L, 3L)) else sub
  fit <- switch(method,
    rao_scott   = .nest_rao_scott(work),
    donner_klar = .nest_donner_klar(work),
    cmh         = .nest_cmh(work),
    gee         = .nest_gee(work),
    glmm        = .nest_glmm(work, levels)
  )
  res_sub <- if (flip) t(fit$residuals) else fit$residuals

  residuals <- matrix(NA_real_, nrow(tab), ncol(tab), dimnames = dimnames(tab))
  residuals[keep_r, keep_c] <- res_sub
  test <- data.frame(
    test          = .nested_labels[[method]],
    statistic     = fit$statistic,
    df            = fit$df,
    df2           = fit$df2 %||% NA_real_,
    p_value       = fit$p_value,
    clusters      = fit$clusters %||% dim(sub)[3L],
    design_effect = fit$design_effect %||% NA_real_,
    icc           = fit$icc %||% NA_real_,
    stringsAsFactors = FALSE
  )
  list(residuals = residuals, test = test)
}


# Element [i, j, k] = m[i, k] (an R x K matrix spread over C columns).
#' @noRd
.nest_spread_rows <- function(m, n_cols) {
  aperm(array(m, c(nrow(m), ncol(m), n_cols)), c(1L, 3L, 2L))
}

# Pearson X2 and adjusted residuals of a pooled table without empty margins.
#' @noRd
.nest_pearson <- function(tab) {
  N <- sum(tab)
  r <- rowSums(tab) / N
  cc <- colSums(tab) / N
  E <- outer(rowSums(tab), colSums(tab)) / N
  list(x2 = sum((tab - E)^2 / E),
       stdres = (tab - E) / sqrt(E * outer(1 - r, 1 - cc)),
       expected = E, r = r, c = cc, N = N)
}

# A ratio that is undefined (zero or non-finite denominator) is NA.
#' @noRd
.nest_safe_div <- function(num, den) {
  out <- num / den
  out[!is.finite(out) | !(den > 0)] <- NA_real_
  out
}


# ---- Rao-Scott ----
# Rao & Scott (1981, 1984); the computation follows survey::svychisq(
# statistic = "F") for an unweighted one-stage cluster sample with
# with-replacement variance: the generalised design effects are the
# eigenvalues of Delta, X2 / tr(Delta) is referred to F(d0, d0 * (K - 1)).
# Per cell: the adjusted residual divided by the square root of the design
# effect of p_ij - p_i. p_.j (cluster-robust variance over simple-random-
# sampling variance, both by linearisation). With one observation per
# cluster the two variances coincide and the residual is chisq.test()'s.
#' @noRd
.nest_rao_scott <- function(n) {
  R <- dim(n)[1L]; C <- dim(n)[2L]; K <- dim(n)[3L]
  tab <- rowSums(n, dims = 2L)
  pr <- .nest_pearson(tab)
  N <- pr$N
  p <- as.vector(tab) / N
  m <- apply(n, 3L, sum)

  # Cluster totals of the linearised cell proportions, one column per cluster.
  Z <- (matrix(n, R * C, K) - outer(p, m)) / N
  Z <- Z - rowMeans(Z)
  V <- tcrossprod(Z) * K / (K - 1)

  grid <- expand.grid(rows = seq_len(R), cols = seq_len(C))
  X1  <- stats::model.matrix(~ factor(rows) + factor(cols), grid)
  X12 <- stats::model.matrix(~ factor(rows) * factor(cols), grid)
  Cmat <- qr.resid(qr(X1), X12[, -seq_len(R + C - 1L), drop = FALSE])
  ip <- ifelse(p == 0, 0, 1 / p)
  ICm <- Cmat * ip
  denom <- crossprod(Cmat, ICm) / N
  numr  <- crossprod(ICm, V %*% ICm)
  # The generalised design effects are the eigenvalues of denom^-1 numr.
  # survey::svychisq() takes them from solve(denom, numr), which stops when
  # empty cells make denom singular. Here they come from the symmetric form
  # T' numr T with T = V+ L+^(-1/2) on the eigenvectors of denom with
  # positive eigenvalues: identical traces when denom is invertible, the
  # design effects of the estimable contrasts when it is not.
  ed <- eigen((denom + t(denom)) / 2, symmetric = TRUE)
  keep <- ed$values > max(ed$values) * 1e-10
  tmat <- ed$vectors[, keep, drop = FALSE] %*%
    diag(1 / sqrt(ed$values[keep]), sum(keep))
  smat <- crossprod(tmat, numr %*% tmat)
  tr_delta <- sum(diag(smat))
  if (!(is.finite(tr_delta) && tr_delta > 0)) {
    stop(errorCondition(paste0(
      "Rao-Scott: the design effects are zero or undefined (no variation ",
      "between actors in the tested contrasts); the test cannot be formed."),
      class = "nestimate_bad_input", call = NULL))
  }
  d0 <- tr_delta^2 / sum(smat * t(smat))
  nu <- K - 1
  f_stat <- pr$x2 / tr_delta

  # Per-cell design effects of d_ij = p_ij - p_i. p_.j.
  r <- pr$r; cc <- pr$c
  h <- tab / N - 2 * outer(r, cc)
  nik <- apply(n, c(1L, 3L), sum)
  njk <- apply(n, c(2L, 3L), sum)
  cc_nik <- .nest_spread_rows(nik, C) *
    aperm(array(cc, c(C, R, K)), c(2L, 1L, 3L))
  z <- (n - cc_nik - outer(r, njk) - outer(h, m)) / N
  v_clus <- rowSums(z^2, dims = 2L) * K / (K - 1)
  ni <- rowSums(tab); nj <- colSums(tab)
  ni_m <- matrix(ni, R, C); nj_m <- matrix(nj, R, C, byrow = TRUE)
  r_m <- matrix(r, R, C); c_m <- matrix(cc, R, C, byrow = TRUE)
  s_srs <- tab * (1 - c_m - r_m - h)^2 +
    (ni_m - tab) * (c_m + h)^2 +
    (nj_m - tab) * (r_m + h)^2 +
    (N - ni_m - nj_m + tab) * h^2
  v_srs <- s_srs / (N * (N - 1))
  deff <- .nest_safe_div(v_clus, v_srs)
  # A design effect of (numerically) zero leaves the residual undefined, not
  # infinite: no actor-to-actor variation in that contrast.
  deff[!is.na(deff) & deff < 1e-8] <- NA_real_

  list(residuals = pr$stdres / sqrt(deff),
       statistic = f_stat, df = d0, df2 = d0 * nu,
       p_value = stats::pf(f_stat, d0, d0 * nu, lower.tail = FALSE),
       design_effect = tr_delta / sum(keep))
}


# ---- Donner-Klar ----
# Donner, Birkett & Buck (1981); Donner & Klar (1994, 2000). The clusters are
# nested in the levels of one margin (the groups). Pearson X2 with each
# group's terms divided by its inflation factor
# C_g = sum_k m_k (1 + (m_k - 1) rho) / N_g, which is 1 + (m - 1) rho for
# equal cluster sizes. rho is the ANOVA intracluster correlation of the
# category indicators, pooled within groups; a negative estimate is reported
# and used as 0. Residuals: adjusted residuals divided by sqrt(C_g).
#' @noRd
.nest_donner_klar <- function(n) {
  rows_per_k <- colSums(apply(n, c(1L, 3L), sum) > 0)
  cols_per_k <- colSums(apply(n, c(2L, 3L), sum) > 0)
  if (all(rows_per_k == 1L)) {
    fit <- .nest_donner_klar_groups(n)
    return(fit)
  }
  if (all(cols_per_k == 1L)) {
    fit <- .nest_donner_klar_groups(aperm(n, c(2L, 1L, 3L)))
    fit$residuals <- t(fit$residuals)
    return(fit)
  }
  stop(errorCondition(paste0(
    "Donner-Klar needs one variable that is constant within every cluster ",
    "(clusters nested in its levels, as in a cluster-randomised design); ",
    "here both variables vary within clusters. Use \"rao_scott\", \"cmh\", ",
    "\"gee\" or \"glmm\"."),
    class = "nestimate_bad_input", call = NULL))
}

#' @noRd
.nest_donner_klar_groups <- function(n) {
  G <- dim(n)[1L]; C <- dim(n)[2L]; K <- dim(n)[3L]
  tab <- rowSums(n, dims = 2L)
  pr <- .nest_pearson(tab)
  N <- pr$N
  group_k <- apply(apply(n, c(1L, 3L), sum) > 0, 2L, which)
  y <- apply(n, c(2L, 3L), sum)
  m <- colSums(y)
  if (K - G < 1L || N - K < 1L) {
    stop(errorCondition(paste0(
      "Donner-Klar needs more clusters than groups and at least one cluster ",
      "with two or more observations, to estimate the intracluster correlation."),
      class = "nestimate_bad_input", call = NULL))
  }
  n_g <- rowSums(tab)
  sq_k <- colSums(y^2) / m
  ssw <- sum(m - sq_k)
  ssc <- sum(sq_k) - sum(rowSums(tab^2) / n_g)
  msc <- ssc / (K - G)
  msw <- ssw / (N - K)
  m0 <- (N - sum(m^2 / n_g[group_k])) / (K - G)
  rho_den <- msc + (m0 - 1) * msw
  if (!(rho_den > 0)) {
    stop(errorCondition(
      "Donner-Klar: the intracluster correlation is undefined (no variation within or between clusters).",
      class = "nestimate_bad_input", call = NULL))
  }
  rho <- (msc - msw) / rho_den
  rho_used <- max(rho, 0)
  infl <- m * (1 + (m - 1) * rho_used)
  c_g <- as.vector(rowsum(infl, group_k, reorder = TRUE)) / n_g
  x2a <- sum((tab - pr$expected)^2 / (c_g * pr$expected))
  df <- (G - 1) * (C - 1)
  list(residuals = pr$stdres / sqrt(c_g),
       statistic = x2a, df = df,
       p_value = stats::pchisq(x2a, df, lower.tail = FALSE),
       design_effect = sum(infl) / N, icc = rho)
}


# ---- Cochran-Mantel-Haenszel ----
# Mantel & Haenszel (1959); generalised to R x C by Landis, Heyman & Koch
# (1978). Each cluster is a stratum, so the test asks whether the variables
# are associated within clusters. stats::mantelhaen.test() (no continuity
# correction); clusters with fewer than two observations carry no
# information and are left out. Residuals: sum over strata of O - E divided
# by the square root of the summed hypergeometric variances.
#' @noRd
.nest_cmh <- function(n) {
  m_all <- apply(n, 3L, sum)
  n <- n[, , m_all >= 2, drop = FALSE]
  K <- dim(n)[3L]
  R <- dim(n)[1L]; C <- dim(n)[2L]
  nik <- apply(n, c(1L, 3L), sum)
  njk <- apply(n, c(2L, 3L), sum)
  informative <- colSums(nik > 0) >= 2L & colSums(njk > 0) >= 2L
  if (K < 2L) {
    stop(errorCondition(sprintf(paste0(
      "Cochran-Mantel-Haenszel needs at least two actors with two or more ",
      "observations; found %d."), K),
      class = "nestimate_bad_input", call = NULL))
  }
  if (!any(informative)) {
    stop(errorCondition(paste0(
      "Cochran-Mantel-Haenszel stratifies by cluster, so both variables must ",
      "vary within clusters; here no cluster holds two levels of each. Use ",
      "\"donner_klar\", \"rao_scott\", \"gee\" or \"glmm\"."),
      class = "nestimate_bad_input", call = NULL))
  }
  ht <- withCallingHandlers(
    stats::mantelhaen.test(n, correct = FALSE),
    error = function(e) {
      stop(errorCondition(paste0(
        "Cochran-Mantel-Haenszel failed: ", conditionMessage(e),
        ". With many sparse clusters the generalised statistic can be ",
        "singular; use \"rao_scott\"."),
        class = "nestimate_bad_input", call = NULL))
    })
  m <- colSums(nik)
  tab <- rowSums(n, dims = 2L)
  expected <- (nik * rep(1 / m, each = R)) %*% t(njk)
  a <- nik * (rep(m, each = R) - nik)
  a <- a * rep(1 / (m^2 * (m - 1)), each = R)
  b <- njk * (rep(m, each = C) - njk)
  v <- a %*% t(b)
  list(residuals = .nest_safe_div(tab - expected, sqrt(v)),
       statistic = as.numeric(ht$statistic),
       df = as.numeric(ht$parameter),
       p_value = ht$p.value,
       clusters = sum(informative))
}


# ---- GEE ----
# Liang & Zeger (1986). Marginal model for the response (columns) given the
# predictor (rows) with identity link and independence working
# correlation: the estimates are the row-conditional proportions and the
# covariance is the cluster-robust sandwich (no small-sample factor, as
# geepack). Wald X2 for H0: every row has the same conditional distribution,
# with a pseudo-inverse and df = its rank. Residuals: row-conditional minus
# pooled proportion, over its cluster-robust standard error.
#' @noRd
.nest_gee <- function(n) {
  R <- dim(n)[1L]; C <- dim(n)[2L]; K <- dim(n)[3L]
  tab <- rowSums(n, dims = 2L)
  N <- sum(tab)
  ni <- rowSums(tab)
  pc <- tab / ni
  nik <- apply(n, c(1L, 3L), sum)
  njk <- apply(n, c(2L, 3L), sum)
  m <- colSums(nik)
  U <- (n - .nest_spread_rows(nik, C) * as.vector(pc)) / ni

  # Parameters: P(column j | row i) for j < C, row index fastest.
  u_theta <- t(matrix(U[, -C, , drop = FALSE], R * (C - 1L), K))
  v_theta <- crossprod(u_theta)
  theta <- as.vector(pc[, -C, drop = FALSE])
  L <- kronecker(diag(C - 1L), cbind(diag(R - 1L), -1))
  lt <- as.vector(L %*% theta)
  lvl <- L %*% v_theta %*% t(L)
  eig <- eigen((lvl + t(lvl)) / 2, symmetric = TRUE)
  pos <- eig$values > max(eig$values) * 1e-10
  # A rank-deficient covariance (fewer actors than the test's degrees of
  # freedom, or cells that never occur) makes the pseudo-inverse Wald depend
  # on which row is the reference and which column is dropped (a factor of
  # 20 under relabelling was found), so the test is refused.
  if (sum(pos) < nrow(lvl)) {
    stop(errorCondition(sprintf(paste0(
      "GEE: the cluster-robust covariance has rank %d, below the %d degrees ",
      "of freedom of the test (fewer actors than degrees of freedom, or cells ",
      "that never occur), so the Wald test would depend on the order of the ",
      "rows and columns. Use \"rao_scott\"."), sum(pos), nrow(lvl)),
      class = "nestimate_bad_input", call = NULL))
  }
  proj <- crossprod(eig$vectors[, pos, drop = FALSE], lt)
  wald <- sum(proj^2 / eig$values[pos])
  df <- sum(pos)

  pj <- colSums(tab) / N
  u_pool <- aperm(array(njk - outer(pj, m), c(C, K, R)), c(3L, 1L, 2L)) / N
  ue <- U - u_pool
  v_e <- rowSums(ue^2, dims = 2L)
  e <- pc - matrix(pj, R, C, byrow = TRUE)
  list(residuals = .nest_safe_div(e, sqrt(v_e)),
       statistic = wald, df = df,
       p_value = stats::pchisq(wald, df, lower.tail = FALSE))
}


# ---- GLMM ----
# Log-linear Poisson mixed model on the cluster-by-cell counts (Agresti 2013,
# ch. 13): count ~ row * col + random intercepts for cluster-by-row and
# cluster-by-column at every nesting level, which absorb each cluster's own
# row and column composition. Cells of a cluster whose predictor row is empty
# are left out (that cluster has no response distribution for the row). The
# Wald X2 tests the (R - 1)(C - 1) interaction terms. Fitted with nAGQ = 0
# (fixed effects from the penalised IRLS step): the default two-stage fit
# failed to converge in 268 s on 200 teams x 81 cells, nAGQ = 0 took 3 s.
# Residuals use pooled independence contrasts and a model-free sandwich
# variance, clustering influences at the outermost actor level.
#' @noRd
.nest_glmm <- function(n, levels) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop(errorCondition(
      "The \"glmm\" method needs the lme4 package; install it or use \"gee\" or \"rao_scott\".",
      class = "nestimate_missing_package", call = NULL))
  }
  R <- dim(n)[1L]; C <- dim(n)[2L]; K <- dim(n)[3L]
  nik <- apply(n, c(1L, 3L), sum)
  ids <- Reduce(function(acc, j) {
    c(acc, list(paste(acc[[length(acc)]], as.character(levels[[j]]),
                      sep = "\r")))
  }, seq_along(levels)[-1L],
  list(as.character(levels[[1L]])))
  cell <- data.frame(
    count = as.vector(n),
    row   = factor(rep(seq_len(R), times = C * K)),
    col   = factor(rep(rep(seq_len(C), each = R), times = K)),
    k     = rep(seq_len(K), each = R * C)
  )
  lvl_names <- paste0("lv", seq_along(ids))
  cell[lvl_names] <- lapply(ids, function(v) factor(v[cell$k]))
  cell <- cell[as.vector(.nest_spread_rows(nik, C)) > 0, , drop = FALSE]
  random <- paste0("(1 | ", lvl_names, ":row) + (1 | ", lvl_names, ":col)",
                   collapse = " + ")
  fit_glmm <- function(fixed) {
    .nest_classed(
      lme4::glmer(stats::as.formula(paste("count ~", fixed, "+", random)),
                  data = cell, family = stats::poisson,
                  contrasts = list(row = "contr.sum", col = "contr.sum"),
                  nAGQ = 0L),
      "GLMM fit")
  }
  fit <- fit_glmm("row * col")
  beta_all <- lme4::fixef(fit)
  idx <- grep("^row[0-9]+:col[0-9]+$", names(beta_all))
  beta <- beta_all[idx]
  vb <- as.matrix(stats::vcov(fit))[idx, idx, drop = FALSE]
  wald <- as.numeric(crossprod(beta, .nest_classed(solve(vb, beta),
                                                   "GLMM Wald test")))
  df <- length(idx)

  # Residuals: observed minus expected under independence, O_ij - n_i. p_.j,
  # over its cluster-robust standard error. Each actor's influence is
  # psi_ijk = (n_ijk - n_ik p_j) - (n_i. / N)(n_.jk - m_k p_j) (the second
  # term carries the estimation of p_j), and the variance is
  # K / (K - 1) * sum_k (psi_k - mean psi)^2. Model-based expected counts
  # were tried and rejected: conditional on the actor-by-column effects they
  # absorbed the association (no tile beyond +-2 at 91% power), and the
  # population-level ones left the residual SD at 0.75 or 1.26 under the
  # null depending on the variance used.
  N <- sum(n)
  nik_all <- apply(n, c(1L, 3L), sum)
  njk <- apply(n, c(2L, 3L), sum)
  m_k <- colSums(njk)
  p_j <- rowSums(njk) / N
  ni <- rowSums(nik_all)
  e <- n - .nest_spread_rows(nik_all, C) *
    aperm(array(p_j, c(C, R, K)), c(2L, 1L, 3L))
  g <- njk - outer(p_j, m_k)
  psi <- matrix(e, R * C, K) -
    (rep(ni, times = C) / N) * g[rep(seq_len(C), each = R), , drop = FALSE]
  num <- rowSums(matrix(e, R * C, K))
  # Inner actors in the same outer unit are dependent. Sum their influences
  # before squaring, preserving the cross-products within each outer unit.
  outer_psi <- t(rowsum(t(psi), as.character(levels[[1L]]), reorder = FALSE))
  n_outer <- ncol(outer_psi)
  v <- if (n_outer > 1L) {
    rowSums((outer_psi - rowMeans(outer_psi))^2) * n_outer / (n_outer - 1)
  } else rep(NA_real_, R * C)
  list(residuals = matrix(.nest_safe_div(num, sqrt(v)), R, C),
       statistic = wald, df = df,
       p_value = stats::pchisq(wald, df, lower.tail = FALSE))
}

# Evaluate `expr`; an error from the fitting code becomes a classed
# nestimate_bad_input naming the step, so mosaic_test() can report it.
#' @noRd
.nest_classed <- function(expr, what) {
  withCallingHandlers(expr, error = function(e) {
    stop(errorCondition(sprintf("%s failed: %s", what, conditionMessage(e)),
                        class = "nestimate_bad_input", call = NULL))
  })
}


# ---- Per-cluster tables for the mosaics ----

# One cluster id per sequence row of `net`, one column per level of `actor`
# (outermost first). NULL: every sequence is its own cluster.
#' @noRd
.mosaic_cluster_levels <- function(net, actor, n_rows, label) {
  if (is.null(actor)) {
    return(data.frame(sequence = paste(label, seq_len(n_rows), sep = "\r"),
                      stringsAsFactors = FALSE))
  }
  .assert("`actor` must be a character vector of column names" =
            is.character(actor) && length(actor) >= 1L &&
            !anyNA(actor) && !anyDuplicated(actor))
  ids <- lapply(actor, function(col) {
    .permutation_block_ids(net, col, n_rows, label)
  })
  stats::setNames(as.data.frame(ids, stringsAsFactors = FALSE,
                                col.names = paste0("lv", seq_along(ids))),
                  actor)
}

# Sum rows of `counts` (one per sequence) within clusters; the key joins all
# levels, so a cluster is one combination of them.
#' @noRd
.mosaic_sum_clusters <- function(counts, levels) {
  key <- do.call(paste, c(unname(as.list(levels)), sep = "\r"))
  first <- !duplicated(key)
  list(counts = rowsum(counts, key, reorder = FALSE),
       levels = levels[first, , drop = FALSE],
       key = key[first])
}

# The R x C x K array behind a network mosaic, oriented like its table `tab`.
#' @noRd
.mosaic_cluster_array <- function(x, source_class, actor, tab) {
  if (source_class %in% c("netobject", "htna")) {
    out <- .mosaic_network_clusters(x, actor)
  } else if (identical(source_class, "netobject_group")) {
    out <- .mosaic_group_clusters(x, actor, colnames(tab))
  } else {
    stop(errorCondition(sprintf(
      "Clustered residuals are not available for a %s.", source_class),
      class = "nestimate_bad_input", call = NULL))
  }
  counts <- out$counts[rownames(tab), colnames(tab), , drop = FALSE]
  pooled <- rowSums(counts, dims = 2L)
  if (!isTRUE(all.equal(as.vector(pooled), as.vector(unclass(tab)),
                        tolerance = 1e-10))) {
    stop(errorCondition(sprintf(paste0(
      "The per-sequence counts (total %g) do not add up to the mosaic's ",
      "table (total %g), so the clustered residuals would describe a ",
      "different table. Use residuals = \"permutation\" or \"asymptotic\"."),
      sum(pooled), sum(tab)),
      class = "nestimate_bad_input", call = NULL))
  }
  list(counts = counts, levels = out$levels, predictor = out$predictor)
}

# Transitions (or co-occurrences) per sequence, counted as the estimator
# counts them, then summed within clusters. Mosaic rows are the targets and
# columns the sources, so the source ("cols") is the predictor.
#' @noRd
.mosaic_network_clusters <- function(x, actor) {
  method <- .resampling_method(.resolve_method_alias(x$method %||% ""),
                               x$params)
  if (!method %in% c("relative", "frequency", "co_occurrence") ||
      !is.data.frame(x$data) || nrow(x$data) == 0L) {
    stop(errorCondition(sprintf(paste0(
      "Clustered residuals count the transitions of each sequence, which ",
      "needs a relative, frequency or co_occurrence network that keeps its ",
      "sequences; this one is method \"%s\"%s."),
      x$method %||% "unknown",
      if (is.data.frame(x$data)) "" else " without sequence data"),
      class = "nestimate_bad_input", call = NULL))
  }
  states <- rownames(x$weights)
  per <- .precompute_per_sequence(
    .resampling_transition_data(x$data, x$metadata, x$params),
    method, x$params, states)
  lv <- .mosaic_cluster_levels(x, actor, nrow(per), "the network")
  agg <- .mosaic_sum_clusters(per, lv)
  s <- length(states)
  # A row of `per` is as.vector(t(from x to)): column-major (to x from).
  counts <- array(t(agg$counts), c(s, s, nrow(agg$counts)),
                  dimnames = list(states, states, agg$key))
  list(counts = counts, levels = agg$levels, predictor = "cols")
}

# State counts per sequence of every group, summed within clusters. A cluster
# id shared by two groups is one cluster holding both. Rows are the groups,
# the predictor.
#' @noRd
.mosaic_group_clusters <- function(x, actor, states) {
  group_names <- names(x) %||% paste0("Group ", seq_along(x))
  parts <- lapply(seq_along(x), function(i) {
    g <- x[[i]]
    label <- sprintf("group '%s'", group_names[i])
    if (!is.data.frame(g$data) && !is.matrix(g$data)) {
      stop(errorCondition(sprintf(
        "Clustered residuals count the states of each sequence; %s keeps no sequence data.",
        label), class = "nestimate_bad_input", call = NULL))
    }
    m <- as.matrix(g$data)
    vals <- as.character(m)
    row_id <- rep(seq_len(nrow(m)), times = ncol(m))
    ok <- !is.na(vals) & nzchar(vals)
    per <- unclass(table(factor(row_id[ok], levels = seq_len(nrow(m))),
                         factor(vals[ok], levels = states)))
    agg <- .mosaic_sum_clusters(per, .mosaic_cluster_levels(g, actor,
                                                            nrow(m), label))
    data.frame(group = group_names[i],
               key = rep(agg$key, times = length(states)),
               state = rep(states, each = length(agg$key)),
               count = as.vector(agg$counts),
               agg$levels[rep(seq_along(agg$key), times = length(states)), ,
                          drop = FALSE],
               stringsAsFactors = FALSE, check.names = FALSE,
               row.names = NULL)
  })
  long <- do.call(rbind, parts)
  keys <- unique(long$key)
  long$group <- factor(long$group, levels = group_names)
  long$state <- factor(long$state, levels = states)
  long$key <- factor(long$key, levels = keys)
  counts <- unclass(stats::xtabs(count ~ group + state + key, data = long))
  attr(counts, "call") <- NULL
  level_cols <- setdiff(names(long), c("group", "key", "state", "count"))
  levels <- long[match(keys, as.character(long$key)), level_cols,
                 drop = FALSE]
  rownames(levels) <- NULL
  list(counts = counts, levels = levels, predictor = "rows")
}


# One line naming the clustered test, for a plot subtitle or a print.
#' @noRd
.nested_test_label <- function(test) {
  stat <- if (is.na(test$df2)) {
    sprintf("X2(%s) = %.2f", format(round(test$df, 2)), test$statistic)
  } else {
    sprintf("F(%.2f, %.1f) = %.2f", test$df, test$df2, test$statistic)
  }
  extra <- c(if (!is.na(test$clusters))
               sprintf("%d actors", as.integer(test$clusters)),
             if (!is.na(test$design_effect))
               sprintf("design effect %.2f", test$design_effect),
             if (!is.na(test$icc)) sprintf("ICC %.3f", test$icc))
  line <- if (is.na(test$statistic)) {   # Fisher's exact: a p-value only
    sprintf("%s: p %s", test$test, .nested_format_p(test$p_value))
  } else {
    sprintf("%s: %s, p %s", test$test, stat, .nested_format_p(test$p_value))
  }
  if (length(extra)) paste0(line, "; ", paste(extra, collapse = ", ")) else
    line
}

# report = TRUE: the test as the one-row table mosaic_test() returns, printed.
#' @noRd
.mosaic_report <- function(test, method, n) {
  row <- data.frame(method = method, test = test$test,
                    statistic = test$statistic, df = test$df,
                    df2 = test$df2, p_value = test$p_value, n = n,
                    clusters = as.integer(test$clusters),
                    design_effect = test$design_effect, icc = test$icc,
                    stringsAsFactors = FALSE)
  print(row, row.names = FALSE)
  invisible(row)
}

#' @noRd
.nested_format_p <- function(p) {
  txt <- format.pval(p, digits = 3, eps = 1e-16)
  if (startsWith(txt, "<")) paste("<", sub("^<\\s*", "", txt)) else
    paste("=", txt)
}


# Resolve `nested` / `nested_method` for a mosaic call: the method to run, or
# NULL. Arguments the call does not use raise nestimate_ignored_argument:
# `residuals` (or `test`) when nested, `actor` and `nested_method` when not.
#' @noRd
.resolve_nesting <- function(nested, nested_method, other_given, method_given,
                             actor, other = "residuals") {
  .assert("`nested` must be TRUE or FALSE" =
            is.logical(nested) && length(nested) == 1L && !is.na(nested))
  ignored <- function(arg, why) {
    warning(warningCondition(sprintf("`%s` is ignored %s.", arg, why),
                             class = "nestimate_ignored_argument",
                             call = NULL))
  }
  if (nested) {
    if (other_given) ignored(other, "when nested = TRUE")
    return(nested_method)
  }
  if (!is.null(actor)) ignored("actor", "unless nested = TRUE")
  if (method_given) ignored("nested_method", "unless nested = TRUE")
  NULL
}


# Residuals and test line of one mosaic. `method` (from .resolve_nesting())
# fits the nested test on `counts`; with `method` NULL and report = TRUE the
# test is Pearson's on `tab` and the residuals stay the caller's. With
# report = TRUE the test is also printed as a one-row table.
#' @noRd
.mosaic_fit_report <- function(tab, method, report, counts = NULL,
                               levels = NULL, predictor = "rows") {
  .assert("`report` must be TRUE or FALSE" =
            is.logical(report) && length(report) == 1L && !is.na(report))
  if (is.null(method) && !report) return(list(stdres = NULL, subtitle = NULL))
  if (is.null(method)) {
    fit <- .mosaic_test_pearson(unclass(tab))
    fit$residuals <- NULL
    method <- "chisq"
  } else {
    fit <- .nested_association(counts, method, levels = levels,
                               predictor = predictor)
  }
  if (report) .mosaic_report(fit$test, method, sum(tab))
  list(stdres = fit$residuals,
       subtitle = if (report) .nested_test_label(fit$test))
}
