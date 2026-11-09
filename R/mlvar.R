# ---- Multilevel Vector Autoregression (mlVAR) ----

#' Build a Multilevel Vector Autoregression (mlVAR) network
#'
#' @description
#' `build_mlvar()` fits a multilevel vector autoregression to intensive
#' longitudinal panel data (ESM or EMA: many people, each measured on many
#' occasions) and returns its three networks. The *temporal* network is
#' directed and holds the fixed-effect lagged regression coefficients: the
#' weight in row `i`, column `j` is the coefficient of variable `i` at an
#' earlier occasion on variable `j` now, so each edge runs from the earlier
#' variable to the one it predicts.
#' The *contemporaneous* network is undirected: the partial correlations of
#' the within-person residuals at the same occasion. The *between-subjects*
#' network is undirected too: the partial correlations of the person means.
#' Each variable is regressed, in a mixed model, on every variable
#' `lag` occasions earlier and on the other variables' person means; the
#' temporal and contemporaneous effects are fixed and only the intercepts are
#' random.
#'
#' The temporal and contemporaneous edges describe within-person relations
#' averaged over people; the between-subjects edges describe how people who
#' differ in one variable on average differ in another (Epskamp et al.,
#' 2018). Only random intercepts are estimated, so the networks do not show
#' how people differ in their dynamics. A temporal edge is a lagged
#' predictive association, and a contemporaneous edge has no temporal order.
#' [coefs()] gives each temporal coefficient with its standard error,
#' p-value and confidence interval; check them before interpreting any
#' single temporal edge.
#'
#' @details Estimation is done by [idiographic::fit_mlvar()], called with
#' `estimator = "lmer"`, `temporal = "fixed"` and `contemporaneous =
#' "fixed"`. It follows mlVAR's lmer path:
#' \enumerate{
#'   \item Drop rows with NA in id/day/beep and, with `standardize = TRUE`,
#'         grand-mean standardise each variable.
#'   \item Expand the per-(id, day) beep grid and right-join the original
#'         values, giving the augmented panel (`augData`).
#'   \item Add within-person lagged predictors (`L1_*`) and person-mean
#'         predictors (`PM_*`).
#'   \item For each outcome variable fit
#'         `lmer(y ~ within + between-except-own-PM + (1 | id))` with
#'         `REML = FALSE`. Collect the fixed-effect temporal matrix `B`, the
#'         between-effect matrix `Gamma`, the random-intercept SDs (`mu_SD`)
#'         and the lmer residual SDs.
#'   \item Contemporaneous network:
#'         `cor2pcor(D %*% cov2cor(cor(resid)) %*% D)`.
#'   \item Between-subjects network:
#'         `cor2pcor(pseudoinverse(forcePositive(D (I - Gamma))))`.
#' }
#'
#' The result agrees with `mlVAR::mlVAR()` to machine precision (largest
#' difference below 1e-10) on 25 real ESM datasets from `openesm` and 20
#' simulated configurations, and exactly with the implementation Nestimate
#' carried before 0.9.0 on all layers, coefficients and observation counts,
#' across settings of `lag`, `standardize`, `day` and `beep`.
#'
#' The temporal weight matrix has predictors in rows and outcomes in columns,
#' the transpose of mlVAR's `Beta` and the same matrix as
#' `mlVAR::getNet(, "temporal")`. A Nestimate network reads its weight matrix
#' from row to column, so the edges of the temporal network (its edge list
#' and [extract_edges()]) run from the predictor at `t - lag` to the outcome
#' at `t`. Nestimate 1.0.9 and earlier did not transpose it, so every
#' temporal edge pointed from the outcome to its predictor.
#'
#' When the data carry no between-person variance (a random-intercept SD of
#' zero), the between-subjects network cannot be estimated. Since 0.9.0
#' \pkg{idiographic} then raises a warning and the zero matrix is returned by
#' convention; before 0.9.0 the zero matrix was returned silently.
#'
#' @section Dispatch limitation:
#' `plot()` draws the three networks side by side, or one with `layer`,
#' with `cograph::splot()` when \pkg{cograph} is installed and in base
#' graphics otherwise; each constituent can also be drawn on its own
#' (`cograph::splot()` on the temporal network). The constituents are
#' wrappers around their weight matrices and carry no
#' data, so the verbs that resample or reread the data,
#' [bootstrap_network()], [network_reliability()] and
#' [centrality_stability()], raise an error on the fit and on each
#' constituent, and [certainty()] refuses them as not transition networks.
#' [net_centrality()] and `summary()` read the weights and work on each
#' constituent. [coefs()] gives the model output.
#'
#' @param data The panel data: a data frame with one row per measurement
#'   occasion. A column named in `vars`, `id`, `day` or `beep` that is not in
#'   it raises an error.
#' @param vars The columns to model, at least two.
#' @param id The column identifying the person each row belongs to.
#' @param day The column identifying the day (or session) of each row. When
#'   given, lagged pairs are formed only within the same day. Default `NULL`.
#' @param beep The column numbering the measurement occasions. When `NULL`
#'   (default), the row position within each person (and day) is used.
#' @param lag How many occasions back the predictors of the temporal network
#'   are, a whole number of at least 1: `lag = 2` regresses each variable on
#'   every variable two occasions earlier (that lag only, not lags 1 and 2).
#'   Default `1`.
#' @param standardize Whether to grand-mean centre each variable and divide
#'   it by its pooled SD before the panel is augmented, as
#'   `mlVAR::mlVAR(scale = TRUE)` does. Default `FALSE`, matching
#'   `mlVAR::mlVAR(scale = FALSE)` (mlVAR's own default is `TRUE`); `FALSE`
#'   is the only setting for which agreement with mlVAR has been validated.
#'
#' @return A `net_mlvar` object (class `c("net_mlvar", "netobject_group")`):
#'   a named list of three networks, with the model-level results stored as
#'   attributes. `print()` gives the number of subjects, observations and
#'   variables, the lag, and how many temporal edges are significant;
#'   [coefs()] returns the temporal coefficients, one row per outcome and
#'   predictor; `summary()` prints the three weight matrices and the
#'   significant temporal edges and returns the same table as [coefs()],
#'   visibly, so at the console the matrices are followed by the table.
#'   \pkg{idiographic} registers its own `print()` and `summary()` methods
#'   for `net_mlvar`; Nestimate loads it first and registers its methods
#'   after, so these two are the ones that run. Each
#'   network is a `netobject` (class `c("netobject", "cograph_network")`)
#'   built from its weight matrix, without data; the *Dispatch limitation*
#'   section lists the verbs that do not work on it. The parts:
#'   \describe{
#'     \item{`temporal`}{Directed network of the `d x d` matrix of
#'       fixed-effect lagged coefficients: the weight in row `i`, column `j`
#'       is the effect of variable `i` at `t - lag` on variable `j` at `t`,
#'       so edges run from predictor to outcome. Method
#'       `"mlvar_temporal"`.}
#'     \item{`contemporaneous`}{Undirected `d x d` network of the partial
#'       correlations of the within-person lmer residuals. Method
#'       `"mlvar_contemporaneous"`.}
#'     \item{`between`}{Undirected `d x d` network of the partial
#'       correlations of the person means, derived from `D (I - Gamma)`.
#'       Method `"mlvar_between"`.}
#'     \item{coefficient table (attribute `coefs`)}{The data frame [coefs()]
#'       returns: one row per `(outcome, predictor)` pair, with columns
#'       `outcome`, `predictor`, `beta`, `se`, `t`, `p`, `ci_lower`,
#'       `ci_upper` and `significant`.}
#'     \item{attribute `n_obs`}{The number of rows of the augmented panel
#'       after rows with missing values are removed.}
#'     \item{attribute `n_subjects`}{The number of subjects remaining.}
#'     \item{attribute `lag`}{The lag used.}
#'     \item{attribute `standardize`}{Whether the variables were
#'       standardised before augmentation.}
#'   }
#'   The `print()` method returns its input invisibly.
#'
#' @examples
#' # A three-variable ESM panel: 20 people x 20 beeps. `tired` is driven by
#' # `happy` one beep earlier, so the temporal network should recover it.
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   set.seed(1)
#'   n_beep <- 20
#'   ar1 <- function(n, phi) as.numeric(stats::filter(stats::rnorm(n), phi,
#'                                                    method = "recursive"))
#'   panel <- do.call(rbind, lapply(seq_len(20), function(i) {
#'     happy <- ar1(n_beep, 0.4)
#'     data.frame(
#'       id    = i,
#'       beep  = seq_len(n_beep),
#'       happy = happy + stats::rnorm(1),
#'       calm  = ar1(n_beep, 0.3) + stats::rnorm(1),
#'       tired = 0.5 * c(0, happy[-n_beep]) + stats::rnorm(n_beep) +
#'               stats::rnorm(1)
#'     )
#'   }))
#'   fit <- build_mlvar(panel, vars = c("happy", "calm", "tired"),
#'                      id = "id", beep = "beep")
#'   fit
#'   coefs(fit)
#'   summary(fit)
#' }
#'
#' @references
#' Epskamp, S., Waldorp, L. J., M\enc{õ}{o}ttus, R., & Borsboom, D. (2018).
#' The Gaussian graphical model in cross-sectional and time-series data.
#' \emph{Multivariate Behavioral Research}, 53(4), 453--480.
#' \doi{10.1080/00273171.2018.1454823}
#'
#' @seealso [build_network()]
#' @importFrom idiographic fit_mlvar
#' @export
build_mlvar <- function(data, vars, id,
                        day = NULL, beep = NULL,
                        lag = 1L,
                        standardize = FALSE) {
  # ---- Input validation ----
  .assert(
    is.data.frame(data),
    is.character(vars), length(vars) >= 2L,
    is.character(id), length(id) == 1L,
    is.numeric(lag), length(lag) == 1L, lag >= 1L,
    is.logical(standardize), length(standardize) == 1L
  )

  required <- c(vars, id)
  if (!is.null(day))  required <- c(required, day)
  if (!is.null(beep)) required <- c(required, beep)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(errorCondition(paste0("Columns not found in data: ",
         paste(missing_cols, collapse = ", ")),
                        class = "nestimate_bad_input", call = NULL))
  }

  lag <- as.integer(lag)

  # The importFrom(idiographic, fit_mlvar) above is load-bearing: it loads
  # idiographic's namespace before Nestimate registers its S3 methods, so
  # Nestimate's print/summary.net_mlvar are registered last and win. Without
  # it idiographic loaded lazily here and replaced them on the first fit.
  # Estimation is owned by idiographic; verified exactly equal (max abs
  # delta == 0 on all three layers, coefs, and n_obs) to the former
  # in-package lmer pipeline across lag/standardize/day/beep configs.
  fit <- idiographic::fit_mlvar(
    data, vars = vars, id = id, day = day, beep = beep,
    lags = lag,
    estimator = "lmer", temporal = "fixed", contemporaneous = "fixed",
    scale = standardize
  )

  stopifnot(
    "idiographic::fit_mlvar() did not return the three network layers" =
      all(c("temporal", "contemporaneous", "between") %in% names(fit)),
    "idiographic::fit_mlvar() result carries no coefficient table" =
      is.data.frame(attr(fit, "coefs"))
  )

  # Re-wrap each layer through Nestimate's own constructor so the return
  # contract (fields, meta, class) is byte-identical to the pre-delegation
  # object regardless of idiographic's internal netobject flavour.
  # idiographic (like mlVAR's Beta) stores the temporal matrix as outcome x
  # predictor. A netobject reads its weights from row to column, so the
  # matrix is transposed to predictor x outcome: the edge then runs from the
  # lagged predictor to the outcome, as mlVAR::getNet(, "temporal") and
  # plot.mlVAR() draw it (both use t(Beta)).
  nets <- list(
    temporal        = .wrap_netobject(t(fit$temporal$weights),
                                      method   = "mlvar_temporal",
                                      directed = TRUE),
    contemporaneous = .wrap_netobject(fit$contemporaneous$weights,
                                      method   = "mlvar_contemporaneous",
                                      directed = FALSE),
    between         = .wrap_netobject(fit$between$weights,
                                      method   = "mlvar_between",
                                      directed = FALSE)
  )

  # Model-level metadata lives in attributes so the list stays a pure
  # netobject_group (each element is a netobject). Use coefs(fit) to
  # retrieve the tidy coefs data.frame.
  # Two-sided normal p-values from the upper tail, 2 * pnorm(-|t|):
  # idiographic's 2 * (1 - pnorm(|t|)) underflows to exactly 0 once |t|
  # exceeds about 8.3.
  coef_table <- attr(fit, "coefs")
  coef_table$p <- 2 * stats::pnorm(-abs(coef_table$t))
  coef_table$significant <- coef_table$p < 0.05
  attr(nets, "coefs")       <- coef_table
  attr(nets, "n_obs")       <- attr(fit, "n_obs")
  attr(nets, "n_subjects")  <- attr(fit, "n_subjects")
  attr(nets, "lag")         <- lag
  attr(nets, "standardize") <- standardize
  attr(nets, "group_col")   <- "network_type"

  class(nets) <- c("net_mlvar", "netobject_group")
  nets
}

#' Tidy coefficients from a fitted mlvar model
#'
#' @description
#' `coefs()` returns the temporal coefficients of a [build_mlvar()] fit as a
#' data frame, one row per `(outcome, predictor)` pair: the fixed-effect
#' lagged regression coefficients that make up the temporal network, with
#' their standard errors and tests.
#'
#' Only the within-person (temporal) coefficients are tabulated. The
#' between-subjects effects that enter the between-subjects network through
#' the `D (I - Gamma)` transformation are not given as a table of their own.
#' The p-values are not adjusted for the number of coefficients, and
#' `significant` is `p < 0.05` for each coefficient on its own.
#'
#' @param x A fitted model; only a `net_mlvar` from [build_mlvar()] has a
#'   method, and any other object raises an error.
#' @param ... Not used.
#' @return A data frame with one row per `(outcome, predictor)` pair and the
#'   columns `outcome` (the variable predicted), `predictor` (the variable at
#'   the earlier occasion), `beta` (the fixed-effect coefficient), `se` (its
#'   standard error), `t` (`beta / se`), `p` (the two-sided p-value from the
#'   normal distribution, `2 * pnorm(-abs(t))`, which stays above zero for
#'   large `t`), `ci_lower` and `ci_upper` (the 95% Wald interval,
#'   `beta` plus or minus `qnorm(0.975)` standard errors) and `significant`
#'   (`p < 0.05`).
#' @inherit build_mlvar examples
#' @export
coefs <- function(x, ...) {
  UseMethod("coefs")
}

#' @rdname coefs
#' @export
coefs.net_mlvar <- function(x, ...) {
  attr(x, "coefs")
}

#' @rdname coefs
#' @export
coefs.default <- function(x, ...) {
  stop(errorCondition(paste0("No coefs() method for object of class '",
       class(x)[1], "'"),
                      class = "nestimate_bad_input", call = NULL))
}

# ---- S3 methods --------------------------------------------------------

#' @rdname build_mlvar
#' @param x For the `print()` method: a `net_mlvar`.
#' @param ... For the `print()` and `summary()` methods: not used.
#' @param object For the `summary()` method: a `net_mlvar`.
#' @export
print.net_mlvar <- function(x, ...) {
  coef_df <- attr(x, "coefs")
  d <- nrow(x$temporal$weights)
  n_sig <- sum(coef_df$significant, na.rm = TRUE)
  n_tot <- nrow(coef_df)

  cat(sprintf("mlVAR result: %d subjects, %d observations, %d variables (lag %d)\n",
              attr(x, "n_subjects"),
              attr(x, "n_obs"),
              d,
              attr(x, "lag")))
  cat(sprintf("  $temporal        %d x %d directed    (%d/%d edges significant at p<0.05)\n",
              d, d, n_sig, n_tot))
  cat(sprintf("  $contemporaneous %d x %d undirected\n", d, d))
  cat(sprintf("  $between         %d x %d undirected\n", d, d))
  cat("Access: fit$temporal, fit$contemporaneous, fit$between (each a netobject);",
      "coefs(fit) for tidy coefs.\n")
  invisible(x)
}

#' @rdname build_mlvar
#' @section Methods:
#' * `summary.net_mlvar()` prints the three weight matrices and the
#'   significant temporal edges, then returns the coefficient table of
#'   [coefs()].
#' @export
summary.net_mlvar <- function(object, ...) {
  coef_df <- attr(object, "coefs")
  vars    <- rownames(object$temporal$weights)
  d       <- length(vars)

  cat("=== mlVAR Summary ===\n")
  cat("Subjects:",      attr(object, "n_subjects"),
      " | Observations:", attr(object, "n_obs"),
      " | Variables:",  d,
      " | Lag:",        attr(object, "lag"), "\n")
  cat("Standardized:",  attr(object, "standardize"), "\n\n")

  cat("--- Temporal Network (rows: predictor at t - lag; columns: outcome) ---\n")
  print(round(object$temporal$weights, 4))
  cat("\n")

  sig_rows <- coef_df[coef_df$significant, , drop = FALSE]
  if (nrow(sig_rows) > 0L) {
    cat("Significant temporal edges (p < 0.05):\n")
    sig_print <- sig_rows[, c("outcome", "predictor",
                              "beta", "se", "t", "p"), drop = FALSE]
    sig_print$beta <- round(sig_print$beta, 4)
    sig_print$se   <- round(sig_print$se, 4)
    sig_print$t    <- round(sig_print$t, 3)
    sig_print$p    <- signif(sig_print$p, 3)
    rownames(sig_print) <- NULL
    print(sig_print)
  } else {
    cat("No significant temporal edges at p < 0.05.\n")
  }
  cat("\n")

  cat("--- Contemporaneous Network ---\n")
  print(round(object$contemporaneous$weights, 4))
  cat("\n")

  cat("--- Between-Subjects Network ---\n")
  print(round(object$between$weights, 4))
  cat("\n")

  coef_df
}

#' @rdname build_mlvar
#' @param layer For the `plot()` method: `"temporal"`, `"contemporaneous"`
#'   or `"between"` to draw that network alone. Default `NULL` draws all
#'   three side by side.
#' @export
plot.net_mlvar <- function(x, layer = NULL, ...) {
  # Nestimate draws its own layers: idiographic's plot.net_mlvar transposes
  # the temporal matrix (it stores outcome x predictor), which would reverse
  # the predictor -> outcome edges Nestimate stores.
  layers <- c("temporal", "contemporaneous", "between")
  if (!is.null(layer)) layers <- match.arg(layer, layers)
  if (length(layers) > 1L) {
    op <- graphics::par(mfrow = c(1L, length(layers)), mar = c(2, 2, 3, 2))
    on.exit(graphics::par(op), add = TRUE, after = FALSE)
  }
  lapply(layers, \(nm) .mlvar_draw_layer(x[[nm]], nm, ...))
  invisible(x)
}

# One mlVAR layer: cograph::splot() when installed, otherwise the base
# graphics circle shared with compare_networks().
.mlvar_draw_layer <- function(net, nm, ...) {
  ttl <- paste0(toupper(substring(nm, 1L, 1L)), substring(nm, 2L))
  if (requireNamespace("cograph", quietly = TRUE)) {
    cograph::splot(net, title = ttl, ...)
  } else {
    .cn_draw_network_base(net$weights, ttl, signed = TRUE,
                          directed = isTRUE(net$directed))
  }
  invisible(NULL)
}
