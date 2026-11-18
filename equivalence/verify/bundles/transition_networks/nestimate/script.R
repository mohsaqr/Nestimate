library(Nestimate)
sequences <- read.csv("../data/sequences.csv")
relative <- extract_edges(build_network(sequences, method = "relative"), include_self = TRUE)
frequency <- extract_edges(build_network(sequences, method = "frequency"), include_self = TRUE)
attention <- extract_edges(build_network(sequences, method = "attention"), include_self = TRUE)
result <- rbind(transform(relative, kind = "relative"),
                transform(frequency, kind = "frequency"),
                transform(attention, kind = "attention"))
result <- transform(result, quantity = paste0(kind, ": ", from, " -> ", to), value = weight)
result
