# build_latency(): reply weights computed directly from their definition on a small message log. Within each
# actor's messages in time order (a gap above time_threshold starts a new sequence), each message is a reply from
# its sender to the sender of the message before it, weighted 1 - (1 - floor) * min(t / cap, 1) for a reply that
# took t seconds; the weights of a pair are summed or averaged.

latency_log <- function() {
  start <- as.POSIXct("2024-03-01 10:00:00", tz = "UTC")
  data.frame(
    team = rep(c("t1", "t2"), c(7L, 6L)),
    who = c("A", "B", "A", "C", "C", "B", "A", "B", "C", "B", "A", "A", "C"),
    at = format(start + c(0, 20, 50, 130, 140, 400, 2000, 0, 5, 65, 95, 300, 330), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

# The reference: reply pairs and their latencies, from the definition above.
latency_reference <- function(data, floor, cap, aggregate, time_threshold = 900) {
  data$secs <- as.numeric(as.POSIXct(data$at, tz = "UTC"))
  replies <- do.call(rbind, lapply(split(data, data$team), \(team) {
    team <- team[order(team$secs), ]
    gap <- c(Inf, diff(team$secs))
    sequence_id <- cumsum(gap > time_threshold)
    keep <- gap <= time_threshold
    data.frame(from = team$who[keep], to = team$who[which(keep) - 1L], t = gap[keep],
               seq = sequence_id[keep], stringsAsFactors = FALSE)
  }))
  replies$w <- 1 - (1 - floor) * pmin(replies$t / cap, 1)
  senders <- sort(unique(data$who))
  out <- matrix(0, length(senders), length(senders), dimnames = list(senders, senders))
  agg <- stats::aggregate(w ~ from + to, data = replies, FUN = if (aggregate == "sum") sum else mean)
  out[cbind(agg$from, agg$to)] <- agg$w
  out
}

test_that("build_latency() weights equal the reply weights computed from their definition", {
  data <- latency_log()
  cases <- list(list(floor = 0.5, cap = 60, aggregate = "sum"), list(floor = 0.2, cap = 120, aggregate = "sum"),
                list(floor = 0, cap = 30, aggregate = "mean"), list(floor = 1, cap = 60, aggregate = "sum"))
  invisible(lapply(cases, \(case) {
    net <- build_latency(data, action = "who", time = "at", actor = "team", floor = case$floor,
                         cap = sprintf("%ds", as.integer(case$cap)), aggregate = case$aggregate)
    expect_reference(extract_transition_matrix(net), latency_reference(data, case$floor, case$cap, case$aggregate),
                     reference = "independent: reply weights from the definition, base R",
                     covers = c("build_latency()", "build_network()"),
                     quantity = sprintf("floor %g, cap %gs, %s", case$floor, case$cap, case$aggregate),
                     tolerance = 1e-12)
  }))
})
