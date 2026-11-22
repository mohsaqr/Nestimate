# Layer 4, recovery of known parameters: latency networks. Message logs are drawn with a known structure: who
# replies to whom follows a Markov chain over three senders, and every pair's reply time is lognormal with its own
# median (fast, medium and slow ties). The mean reply weight of each edge then has a closed form, and
# bootstrap_network()'s percentile intervals must cover it at the nominal rate; an edge its stability test calls
# stable must be within 25% of its true speed score. The study behind the design is in the package's
# development notes (tmp/latency-resampling, not shipped).

# ---- Criteria, fixed before the first replicate is drawn ---------------------------------------------------------
BAND_LEVEL <- 0.999
CI_LEVEL <- 0.05
R_COVERAGE <- 200L
BOOT_ITER <- 200L
N_THREADS <- 300L
FLOOR <- 0.5
CAP <- 0.95
SDLOG <- 0.8
# An edge called stable is within 25% of its true speed score in at least this share of the calls.
STABLE_VALID <- 0.95

senders <- c("A", "B", "C")
# who replies: row = previous sender, column = responder
reply_chain <- matrix(c(0.1, 0.6, 0.3,
                        0.5, 0.1, 0.4,
                        0.4, 0.5, 0.1), 3, 3, byrow = TRUE, dimnames = list(senders, senders))
# median reply time of responder (row) to previous sender (column), seconds
reply_median <- matrix(c(30, 10, 120,
                         10, 30, 30,
                        120, 30, 30), 3, 3, byrow = TRUE, dimnames = list(senders, senders))

lt_log <- function(n_threads) {
  do.call(rbind, lapply(seq_len(n_threads), function(i) {
    len <- sample(10:30, 1L)
    s <- integer(len)
    s[1L] <- sample.int(3L, 1L)
    gap <- numeric(len)
    for (t in seq_len(len)[-1L]) {  # each sender depends on the one before: sequential
      s[t] <- sample.int(3L, 1L, prob = reply_chain[s[t - 1L], ])
      gap[t] <- exp(log(reply_median[s[t], s[t - 1L]]) + SDLOG * stats::rnorm(1L))
    }
    data.frame(thread = sprintf("t%04d", i), sender = senders[s], at = 86400 * i + cumsum(gap))
  }))
}

# expected number of a-replies-to-b per thread (responder a, previous b), averaged over the lengths
lt_reply_counts <- function() {
  per_length <- lapply(10:30, function(len) {
    dist <- rep(1 / 3, 3); total <- matrix(0, 3, 3)
    for (t in seq_len(len - 1L)) {  # propagation of the chain: sequential
      total <- total + dist * reply_chain
      dist <- as.vector(dist %*% reply_chain)
    }
    t(total)
  })
  Reduce(`+`, per_length) / length(per_length)
}

# the population cap: the CAP-quantile of the sender-change reply times, a mixture over the off-diagonal pairs
lt_population_cap <- function() {
  counts <- lt_reply_counts(); off <- row(counts) != col(counts)
  weights <- counts[off] / sum(counts[off]); mu <- log(reply_median[off])
  stats::uniroot(function(x) sum(weights * stats::pnorm((log(x) - mu) / SDLOG)) - CAP, c(1e-3, 1e6),
                 tol = 1e-10)$root
}

# the true mean reply weight of every edge at cap `cap`
lt_true_weights <- function(cap) {
  mu <- log(reply_median)
  capped <- exp(mu + SDLOG^2 / 2) * stats::pnorm((log(cap) - mu - SDLOG^2) / SDLOG) +
    cap * (1 - stats::pnorm((log(cap) - mu) / SDLOG))
  out <- 1 - (1 - FLOOR) * capped / cap
  dimnames(out) <- dimnames(reply_median)
  out
}

test_that("bootstrap_network() intervals of a latency network cover the true mean reply weights", {
  band <- rc_binomial_band(R_COVERAGE, 1 - CI_LEVEL, BAND_LEVEL)
  truth <- lt_true_weights(lt_population_cap())
  true_score <- (truth - FLOOR) / (1 - FLOOR)
  results <- rc_replicate(seq_len(R_COVERAGE), \(r) {
    log <- rc_with_seed(9100L + r, lt_log(N_THREADS))
    net <- build_network(log, method = "latency", action = "sender", time = "at", session = "thread",
                         time_threshold = FALSE, latency_floor = FLOOR, latency_cap = CAP)
    s <- summary(bootstrap_network(net, iter = BOOT_ITER, ci_level = CI_LEVEL, seed = r))
    edge <- cbind(s$from, s$to)
    score <- (s$weight - FLOOR) / (1 - FLOOR)
    list(covered = stats::setNames(s$ci_lower <= truth[edge] & truth[edge] <= s$ci_upper, paste(s$from, s$to)),
         stable_within = abs(score - true_score[edge]) <= 0.25 * true_score[edge],
         stable = s$sig)
  })
  coverage <- rowMeans(simplify2array(lapply(results, \(x) x$covered)))
  expect_reference(coverage, stats::setNames(rep(1 - CI_LEVEL, length(coverage)), names(coverage)),
                   reference = "independent: closed-form mean reply weight of lognormal reply times",
                   covers = "bootstrap_network()",
                   quantity = sprintf("latency network: coverage of the 95%% percentile interval, %d draws, R = %d",
                                      BOOT_ITER, R_COVERAGE),
                   tolerance = band[["tolerance"]],
                   note = sprintf("99.9%% binomial range [%.3f, %.3f]; %d threads, cap at the %.2f quantile",
                                  band[["lower"]], band[["upper"]], N_THREADS, CAP))
  calls <- unlist(lapply(results, \(x) x$stable_within[x$stable]))
  expect_reference(c(stable_within_25pct = mean(calls)), c(stable_within_25pct = 1),
                   reference = "independent: closed-form mean reply weight of lognormal reply times",
                   covers = "bootstrap_network()",
                   quantity = sprintf("latency network: share of stable calls within 25%% of the true speed score, %d calls",
                                      length(calls)),
                   tolerance = 1 - STABLE_VALID)
})
