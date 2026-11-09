# Discover Association Rules from Sequential or Transaction Data

`association_rules()` finds association rules of the form {A} -\> {B}
with the Apriori algorithm (Agrawal & Srikant, 1994): sets of items that
occur together in at least `min_support` of the transactions, and the
rules among them that pass `min_confidence` and `min_lift`. It reads a
list of transactions, a wide data frame or matrix of states, a 0/1
incidence table, or a network: the sequences behind a `netobject` or
`tna` model become transactions of their distinct states, and a one-hot
network's windows become transactions of their active codes.

Support is how common the whole rule is, confidence how often B occurs
when A does, and lift how much more often A and B occur together than
they would if they were independent (see *Metrics*). A sequence read as
a transaction loses its order, so a rule says that states occur in the
same sequence, not that one follows the other;
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
models order. The rules are not tested for significance, and lowering
the thresholds quickly multiplies them, so check whether a rule of
interest stays when `min_support` changes.

## Usage

``` r
association_rules(
  x,
  min_support = 0.1,
  min_confidence = 0.5,
  min_lift = 1,
  max_length = 5L,
  format = c("auto", "binary")
)

# S3 method for class 'net_association_rules'
print(x, ...)

# S3 method for class 'net_association_rules'
summary(object, ...)

# S3 method for class 'net_association_rules'
plot(x, ...)
```

## Arguments

- x:

  The transactions, in one of these forms:

  `netobject`

  :   Each sequence in its data becomes a transaction of its distinct
      states. In a one-hot network the codes are the items and each
      window is one transaction. A network without data, and an
      association network (`cor`, `pcor`, `glasso`, `ising`, `mgm`),
      whose rows are variables, raise an error; a `netobject_group`
      raises one of class `nestimate_bad_input` (pass one group's
      network). A bare `cograph_network` is read as a `netobject`, a
      `tna` model by its sequence data, and an `mcml` by its
      between-cluster network, so its items are the clusters.

  list

  :   A plain list holding one character vector of items per
      transaction.

  data frame

  :   Each row is a transaction of its distinct non-missing states. The
      frame must be a single column or a numbered run of position
      columns (`T1`, `T2`, ...) with nothing beside it; otherwise it
      raises an error of class `nestimate_bad_input`; declare the
      columns with `build_network(state_cols = )` and pass the network
      instead. A column called `action` and repeated column names are
      refused the same way. A frame of 0/1 or logical columns is a
      binary table, one column per item, read with `format = "binary"`.

  matrix

  :   A character matrix is read like a character data frame. A numeric
      or logical matrix is a binary transaction table (rows are
      transactions, columns are items, named `I1`, `I2`, ... when the
      matrix has no column names).

  Binary input must hold only 0/1 or `TRUE`/`FALSE`: a missing cell or
  any other value raises an error of class `nestimate_bad_input`. For
  the [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods: a
  `net_association_rules` object.

- min_support:

  The smallest share of transactions an itemset must occur in to be
  frequent, above 0 and at most 1. Default `0.1`. It is applied as a
  count, `ceiling(min_support * n_transactions)`.

- min_confidence:

  The smallest confidence a rule must have, between 0 and 1. Default
  `0.5`.

- min_lift:

  The smallest lift a rule must have, at least 0. Default `1`, which
  keeps rules whose items occur together at least as often as
  independence predicts.

- max_length:

  The largest itemset size searched, at least 2. Default `5`.

- format:

  How 0/1 columns are read. Columns of only 0 and 1 could be items or
  recoded positions, so with `"auto"` (default) they raise an error of
  class `nestimate_bad_input`. `"binary"` declares them items, one row
  per transaction; it raises the same error on a data frame with any
  other value.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods:
  ignored.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_association_rules` object.

## Value

A `net_association_rules` object.
[`print()`](https://rdrr.io/r/base/print.html) shows the counts, the
thresholds and the ten rules with the highest lift;
[`summary()`](https://rdrr.io/r/base/summary.html) returns the rules
table, one row per rule;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the rules
by support and confidence. The object's fields:

- rules:

  Data frame, one row per rule, ordered by decreasing lift and then
  confidence, with columns `antecedent` and `consequent` (the itemsets
  as comma-separated strings), `support`, `confidence`, `lift`,
  `conviction`, `count` (transactions holding the whole rule) and
  `n_transactions`.

- frequent:

  Data frame, one row per frequent itemset, with columns `itemset`,
  `size`, `support` and `count`.

- frequent_itemsets:

  List of the frequent itemsets per size k, each entry a list of
  `items`, `count` and `support`.

- items:

  Character vector of the frequent single items the mining ran on (all
  items when none reaches `min_support`).

- n_transactions:

  Number of transactions.

- n_rules:

  Number of rules.

- params:

  List of `min_support`, `min_confidence`, `min_lift` and `max_length`.

The [`print()`](https://rdrr.io/r/base/print.html) method returns its
input invisibly. [`summary()`](https://rdrr.io/r/base/summary.html)
returns the rules data frame described above.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints a
scatter plot of support against confidence, with point size and colour
showing lift, and returns the `ggplot` object invisibly; with no rules
it says so in a message and returns `NULL` invisibly.

## Details

### Algorithm

Level-wise Apriori (Agrawal & Srikant, 1994) with the full pruning step:
after the join step generates the candidates of size k, every subset of
size k - 1 is checked to be frequent before support is counted, which
keeps the search fast at k \>= 4. Support of pairs is counted with
[`crossprod()`](https://rdrr.io/r/base/crossprod.html), of larger
itemsets by logical indexing. Every split of a frequent itemset into a
non-empty antecedent and consequent is a candidate rule.

### Metrics

- support:

  P(A and B): the share of transactions holding both the antecedent and
  the consequent.

- confidence:

  P(B \| A): the share of transactions holding the antecedent that also
  hold the consequent.

- lift:

  P(A and B) / (P(A) P(B)). Values above 1 indicate positive
  association, below 1 negative association.

- conviction:

  (1 - P(B)) / (1 - confidence), the departure from independence in the
  direction of the rule (Brin et al., 1997); higher means a stronger
  implication, and a rule with confidence 1 has conviction `Inf`.

## References

Agrawal, R. & Srikant, R. (1994). Fast algorithms for mining association
rules. In *Proc. 20th VLDB Conference*, 487–499.

Brin, S., Motwani, R., Ullman, J. D. & Tsur, S. (1997). Dynamic itemset
counting and implication rules for market basket data. In *Proc. ACM
SIGMOD*, 255–264. (lift and conviction)

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
[`predict_links`](https://pak.dynasite.org/Nestimate/reference/predict_links.md)

## Examples

``` r
# From a list of transactions
trans <- list(
  c("plan", "discuss", "execute"),
  c("plan", "research", "analyze"),
  c("discuss", "execute", "reflect"),
  c("plan", "discuss", "execute", "reflect"),
  c("research", "analyze", "reflect")
)
rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.5)
print(rules)
#> Association Rules  [24 rules | 6 items | 5 transactions]
#>   Support >= 0.30  |  Confidence >= 0.50  |  Lift >= 1.00
#> 
#>   Top rules (by lift):
#>     1. analyze -> research  (sup=0.400 conf=1.000 lift=2.50)
#>     2. research -> analyze  (sup=0.400 conf=1.000 lift=2.50)
#>     3. discuss -> execute  (sup=0.600 conf=1.000 lift=1.67)
#>     4. execute -> discuss  (sup=0.600 conf=1.000 lift=1.67)
#>     5. discuss, plan -> execute  (sup=0.400 conf=1.000 lift=1.67)
#>     6. execute, plan -> discuss  (sup=0.400 conf=1.000 lift=1.67)
#>     7. discuss, reflect -> execute  (sup=0.400 conf=1.000 lift=1.67)
#>     8. execute, reflect -> discuss  (sup=0.400 conf=1.000 lift=1.67)
#>     9. discuss -> execute, plan  (sup=0.400 conf=0.667 lift=1.67)
#>     10. execute -> discuss, plan  (sup=0.400 conf=0.667 lift=1.67)
#>     ... and 14 more rules

# From a netobject (sequences as transactions)
seqs <- data.frame(
  V1 = sample(LETTERS[1:5], 50, TRUE),
  V2 = sample(LETTERS[1:5], 50, TRUE),
  V3 = sample(LETTERS[1:5], 50, TRUE)
)
net <- build_network(seqs, method = "relative")
rules <- association_rules(net, min_support = 0.1)
```
