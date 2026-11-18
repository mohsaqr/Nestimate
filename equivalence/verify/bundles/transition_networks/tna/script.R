library(tna)
sequences <- read.csv("../data/sequences.csv")
edges_of <- function(model, kind) {
  edges <- as.data.frame(as.table(model$weights), stringsAsFactors = FALSE)
  names(edges) <- c("from", "to", "value")
  edges <- edges[edges$value != 0, ]
  data.frame(quantity = paste0(kind, ": ", edges$from, " -> ", edges$to), kind = kind,
             value = edges$value)
}
result <- rbind(edges_of(tna(sequences), "relative"),
                edges_of(ftna(sequences), "frequency"),
                edges_of(atna(sequences), "attention"))
result
