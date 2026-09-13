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

  A higher-order network object (`net_hon`, `net_hypa`, `net_mogen`), a
  `netobject`, a `net_association_rules` or a `net_link_prediction`.

- ...:

  Additional arguments passed on to the method.

- min_count:

  Integer. Minimum transition count to include (default: 1). Filters
  noise from rare observations. Used by the `net_hon` and `net_mogen`
  methods.

- min_prob:

  Numeric. Minimum transition probability to include (default: 0).
  Useful for filtering weak transitions. Used by the `net_hon` and
  `net_mogen` methods.

- top:

  Integer or NULL. Keep only the top N pathways: ranked by count for
  `net_hon` and `net_mogen`, by lift then confidence for
  `net_association_rules`, and by score for `net_link_prediction`.
  Default `NULL` (all) everywhere except `net_link_prediction`, whose
  default is `10`.

- order:

  Integer or NULL. For `net_hon`, keep only pathways whose source is of
  this order; default `NULL` = every order above 1. For `net_mogen`, the
  Markov order to extract; default `NULL` = the optimal order from model
  selection.

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
Every method returns this shape, and `character(0)` when nothing
survives its filters.

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
  target. Rules whose antecedent and consequent share the same item set
  describe the same simplex, so only the highest-lift rule per item set
  is returned; `top` is applied after that de-duplication.

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
#> [1] "C A B D -> E"
```
