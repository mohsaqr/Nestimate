# Extract Pathways from Higher-Order Network Objects

`pathways()` turns a higher-order result into pathway strings in arrow
notation, `"A B -> C"`: the states before the arrow lead, together, to
the state after it. The strings are what
[`cograph::plot_simplicial()`](https://sonsoles.me/cograph/reference/plot_simplicial.html)
draws, one simplex per pathway. Which pathways are returned depends on
the input:

- a `net_hon` gives its higher-order edges, those whose source node has
  order above 1, i.e. the transitions that differ from a first-order
  Markov chain, ranked by count;

- a `net_hypa` gives its anomalous paths, over- or under-represented
  relative to the hypergeometric null model, the over-represented ones
  first;

- a `net_mogen` gives every transition at the optimal order (or the
  order given), ranked by count;

- a `netobject` built from sequences is first turned into a HON or a
  HYPA, and that result's pathways are returned;

- a `net_association_rules` gives its rules `{A, B} => {C}` as
  `"A B -> C"`: antecedent items become the sources, consequent items
  the target. Rules whose antecedent and consequent together hold the
  same items describe the same simplex, so only the highest-lift rule
  per item set is kept; the rules are ranked by lift, then confidence;

- a `net_link_prediction` gives its predicted links, ranked by score;
  with `evidence = TRUE` each predicted edge `A -> B` is extended with
  the common neighbours that support it, `"A cn1 cn2 -> B"`.

The strings keep only which pathways passed the filters, not the counts,
probabilities, scores or tests behind them; those stay in the source
object (for a MOGen,
[`mogen_transitions()`](https://pak.dynasite.org/Nestimate/reference/mogen_transitions.md)
gives them per order).

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

  The object to read: a `net_hon`, `net_hypa` or `net_mogen` (from
  [`build_hon()`](https://pak.dynasite.org/Nestimate/reference/build_hon.md),
  [`build_hypa()`](https://pak.dynasite.org/Nestimate/reference/build_hypa.md)
  or
  [`build_mogen()`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)),
  a `netobject` built from sequences, a `net_association_rules` (from
  [`association_rules()`](https://pak.dynasite.org/Nestimate/reference/association_rules.md))
  or a `net_link_prediction` (from
  [`predict_links()`](https://pak.dynasite.org/Nestimate/reference/predict_links.md)).
  A `netobject` of an association method (whose rows are variables) or a
  one-hot `netobject` raises an error of class `nestimate_bad_input`; a
  network without its sequences stops `pathways()`.

- ...:

  For the `netobject` method: passed to the `pathways()` method of the
  HON or HYPA built from it, so `min_count`, `min_prob`, `top` and
  `order` apply to a HON and `type` to a HYPA. Ignored by the other
  methods.

- min_count:

  For a HON or a MOGen: the smallest transition count kept, which
  filters out rare observations. Default `1` (all).

- min_prob:

  For a HON or a MOGen: the smallest transition probability kept.
  Default `0` (all).

- top:

  How many pathways to keep after ranking: by count for a HON or a
  MOGen, by lift then confidence for association rules (after the item
  sets are de-duplicated), and by score for link predictions. Default
  `NULL` (all), except for link predictions, whose default is `10`. Not
  used for a HYPA.

- order:

  For a HON: keep only pathways whose source has this order. Default
  `NULL`, every order above 1; order 1 returns nothing. For a MOGen: the
  Markov order to extract. Default `NULL`, the optimal order from model
  selection; an order below 1 returns `character(0)`.

- type:

  For a HYPA: which anomalies to include, `"all"` (default), `"over"`
  (over-represented paths) or `"under"` (under-represented ones).

- ho_method:

  For a `netobject`: the higher-order model built from its sequences,
  `"hon"` (default,
  [`build_hon()`](https://pak.dynasite.org/Nestimate/reference/build_hon.md))
  or `"hypa"`
  ([`build_hypa()`](https://pak.dynasite.org/Nestimate/reference/build_hypa.md)),
  each with its default settings. The network must carry the sequences
  it was built from.

- min_lift:

  For association rules: the smallest lift kept, on top of the threshold
  the rules were mined with. Default `NULL` (no further filter).

- min_confidence:

  For association rules: the smallest confidence kept, on top of the
  mining threshold. Default `NULL` (no further filter).

- method:

  For link predictions: the prediction method whose links are returned.
  Default `NULL`, the first method in the object. A method the object
  does not hold returns `character(0)`.

- evidence:

  For link predictions: whether to extend each predicted link `A -> B`
  with its common neighbours, the nodes with an edge from `A` and an
  edge to `B`, giving `"A cn1 cn2 -> B"`. Default `TRUE`. Without a
  stored adjacency matrix the links are returned as `"A -> B"`.

- max_evidence:

  For link predictions with `evidence = TRUE`: the most common
  neighbours added to one link; when there are more, those with the
  largest summed weight of the edges from `A` and to `B` are kept.
  Default `3`.

## Value

A character vector of pathway strings in arrow notation (e.g.
`"A B -> C"`), one per pathway and in the ranking order described above,
ready for
[`cograph::plot_simplicial()`](https://sonsoles.me/cograph/reference/plot_simplicial.html).
Every method returns this shape, and `character(0)` when nothing passes
its filters.

## Methods (by class)

- `pathways(net_hon)`: Higher-order edges of a HON.

- `pathways(net_hypa)`: Anomalous paths of a HYPA.

- `pathways(netobject)`: Higher-order pathways of a network's sequences.

- `pathways(net_association_rules)`: Association rules as pathways.

- `pathways(net_link_prediction)`: Predicted links as pathways.

- `pathways(net_mogen)`: Transitions of a MOGen at one order.

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
