library(tna)
group_a <- read.csv("../data/sequences.csv")
group_b <- read.csv("../data/sequences_b.csv")
model_a <- tna(group_a)
model_b <- tna(group_b)
set.seed(1)
test <- permutation_test(model_a, model_b, iter = 200)
as.data.frame(test$edges$stats)
tested <- which(model_a$weights != 0 | model_b$weights != 0, arr.ind = TRUE)
result <- data.frame(edge = paste(model_a$labels[tested[, "row"]], "->", model_a$labels[tested[, "col"]]),
                     weight_x = model_a$weights[tested], weight_y = model_b$weights[tested],
                     diff = test$edges$diffs_true[tested])
result
