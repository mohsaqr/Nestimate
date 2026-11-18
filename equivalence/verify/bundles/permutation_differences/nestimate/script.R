library(Nestimate)
group_a <- read.csv("../data/sequences.csv")
group_b <- read.csv("../data/sequences_b.csv")
test <- permutation(build_network(group_a, method = "relative"),
                    build_network(group_b, method = "relative"), iter = 200, seed = 1)
test
summary(test)
result <- subset(transform(summary(test), edge = paste(from, "->", to)),
                 select = c(edge, weight_x, weight_y, diff))
result
