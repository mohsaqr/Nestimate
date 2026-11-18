library(tna)
sequences <- read.csv("../data/sequences.csv")
model <- tna(sequences)
result <- as.data.frame(centralities(model, measures = c("OutStrength", "InStrength", "ClosenessIn",
                                                         "ClosenessOut", "Closeness", "Betweenness")))
result
