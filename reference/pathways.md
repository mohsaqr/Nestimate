# Extract Pathways from Higher-Order Network Objects

Extracts higher-order pathway strings suitable for
[`cograph::plot_simplicial()`](https://sonsoles.me/cograph/reference/plot_simplicial.html).
Each pathway represents a multi-step dependency: source states lead to a
target state.

For `net_hon`: extracts edges where the source node is higher-order
(order \> 1), i.e., the transitions that differ from first-order Markov.

For `net_hypa`: extracts anomalous paths (over- or under-represented
relative to the hypergeometric null model).

For `net_mogen`: extracts all transitions at the optimal order (or a
specified order).

## Usage

``` r
pathways(x, ...)

# S3 method for class 'net_hon'
pathways(x, min_count = 1L, min_prob = 0, top = NULL, order = NULL, ...)

# S3 method for class 'net_hypa'
pathways(x, type = "all", ...)

# S3 method for class 'netobject'
pathways(x, ho_method = c("hon", "hypa"), ...)

# S3 method for class 'net_association_rules'
pathways(x, top = NULL, min_lift = NULL, min_confidence = NULL, ...)

# S3 method for class 'net_link_prediction'
pathways(x, method = NULL, top = 10L, evidence = TRUE, max_evidence = 3L, ...)

# S3 method for class 'net_mogen'
pathways(x, order = NULL, min_count = 1L, min_prob = 0, top = NULL, ...)
```

## Arguments

- x:

  A higher-order network object (`net_hon`, `net_hypa`, or `net_mogen`).

- ...:

  Additional arguments.

- min_count:

  Integer. Minimum transition count to include (default: 1).

- min_prob:

  Numeric. Minimum transition probability to include (default: 0).

- top:

  Integer or NULL. Return only the top N pathways ranked by count
  (default: NULL = all).

- order:

  Integer or NULL. Markov order to extract. Default: optimal order from
  model selection.

- type:

  Character. Which anomalies to include: `"all"` (default), `"over"`, or
  `"under"`.

- ho_method:

  Character. Higher-order method: `"hon"` (default) or `"hypa"`.

- min_lift:

  Numeric or NULL. Additional lift filter applied on top of the object's
  original threshold (default: NULL).

- min_confidence:

  Numeric or NULL. Additional confidence filter (default: NULL).

- method:

  Character or NULL. Which prediction method to use. Default: first
  method in the object.

- evidence:

  Logical. If TRUE, include common neighbor evidence nodes in each
  pathway. Default: TRUE.

- max_evidence:

  Integer. Maximum number of evidence nodes per pathway (default: 3).

## Value

A character vector of pathway strings in arrow notation (e.g.
`"A B -> C"`), suitable for
[`cograph::plot_simplicial()`](https://sonsoles.me/cograph/reference/plot_simplicial.html).

A character vector of pathway strings.

A character vector of pathway strings.

A character vector of pathway strings.

A character vector of pathway strings.

A character vector of pathway strings.

A character vector of pathway strings.

## Methods (by class)

- `pathways(net_hon)`: Extract higher-order pathways from HON

- `pathways(net_hypa)`: Extract anomalous pathways from HYPA

- `pathways(netobject)`: Extract pathways from a netobject

  Builds a Higher-Order Network (HON) from the netobject's sequence data
  and returns the higher-order pathways. Requires that the netobject was
  built from sequence data (has `$data`).

- `pathways(net_association_rules)`: Extract pathways from association
  rules

  Converts association rules `{A, B} => {C}` into pathway strings
  `"A B -> C"` suitable for
  [`cograph::plot_simplicial()`](https://sonsoles.me/cograph/reference/plot_simplicial.html).
  Antecedent items become source nodes; consequent items become the
  target.

- `pathways(net_link_prediction)`: Extract pathways from link
  predictions

  Converts predicted links into pathway strings for
  [`cograph::plot_simplicial()`](https://sonsoles.me/cograph/reference/plot_simplicial.html).
  When `evidence = TRUE` (default), each predicted edge `A -> B` is
  enriched with common neighbors that structurally support the
  prediction, producing `"A cn1 cn2 -> B"`.

- `pathways(net_mogen)`: Extract transition pathways from MOGen

## Examples

``` r
# \donttest{
seqs <- list(c("A","B","C","D"), c("A","B","C","A"))
hon <- build_hon(seqs, max_order = 3)
pw <- pathways(hon)
# }

trans <- list(c("A","B","C"), c("A","B"), c("B","C","D"), c("A","C","D"))
rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.3,
                           min_lift = 0)
pathways(rules)
#> [1] "D -> C" "A -> B" "A -> C" "B -> C"

seqs <- data.frame(
  V1 = sample(LETTERS[1:5], 50, TRUE),
  V2 = sample(LETTERS[1:5], 50, TRUE),
  V3 = sample(LETTERS[1:5], 50, TRUE)
)
net <- build_network(seqs, method = "relative")
pred <- predict_links(net, methods = "common_neighbors")
pathways(pred)
#> character(0)
```
