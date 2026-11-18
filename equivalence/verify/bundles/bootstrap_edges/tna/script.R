library(tna)
sequences <- read.csv("../data/sequences.csv")
model <- tna(sequences)
set.seed(1)
boot <- bootstrap(model, iter = 200)
edges <- boot$summary
at <- cbind(match(edges$from, model$labels), match(edges$to, model$labels))
result <- data.frame(edge = paste(edges$from, "->", edges$to), weight = edges$weight,
                     mean = boot$weights_mean[at], sd = boot$weights_sd[at],
                     p_value = edges$p_value, sig = edges$sig,
                     ci_lower = edges$ci_lower, ci_upper = edges$ci_upper,
                     cr_lower = edges$cr_lower, cr_upper = edges$cr_upper)
result
