# Discover Association Rules from Sequential or Transaction Data

Discovers association rules using the Apriori algorithm with proper
candidate pruning. Accepts `netobject` (extracts sequences as
transactions), data frames, lists, or binary matrices.

Support counting is vectorized via
[`crossprod()`](https://rdrr.io/r/base/crossprod.html) for 2-itemsets
and logical matrix indexing for k-itemsets.

## Usage

``` r
association_rules(
  x,
  min_support = 0.1,
  min_confidence = 0.5,
  min_lift = 1,
  max_length = 5L
)
```

## Arguments

- x:

  Input data. Accepts:

  netobject

  :   Uses `$data` sequences — each sequence becomes a transaction of
      its unique states.

  list

  :   Each element is a character vector of items (one transaction).

  data.frame

  :   Wide format: each row is a transaction, character columns are item
      occurrences. Or a binary matrix (0/1).

  matrix

  :   Binary transaction matrix (rows = transactions, columns = items).

- min_support:

  Numeric. Minimum support threshold. Default: 0.1.

- min_confidence:

  Numeric. Minimum confidence threshold. Default: 0.5.

- min_lift:

  Numeric. Minimum lift threshold. Default: 1.0.

- max_length:

  Integer. Maximum itemset size. Default: 5.

## Value

An object of class `"net_association_rules"` containing:

- rules:

  Data frame with columns: antecedent (list), consequent (list),
  support, confidence, lift, conviction, count, n_transactions.

- frequent_itemsets:

  List of frequent itemsets per level k.

- items:

  Character vector of all items.

- n_transactions:

  Integer.

- n_rules:

  Integer.

- params:

  List of min_support, min_confidence, min_lift, max_length.

## Details

### Algorithm

Uses level-wise Apriori (Agrawal & Srikant, 1994) with the full pruning
step: after the join step generates k-candidates, all (k-1)-subsets are
verified as frequent before support counting. This is critical for
efficiency at k \>= 4.

### Metrics

- support:

  P(A and B). Fraction of transactions containing both antecedent and
  consequent.

- confidence:

  P(B \| A). Fraction of antecedent transactions that also contain the
  consequent.

- lift:

  P(A and B) / (P(A) \* P(B)). Values \> 1 indicate positive
  association; \< 1 indicate negative association.

- conviction:

  (1 - P(B)) / (1 - confidence). Measures departure from independence.
  Higher = stronger implication.

## References

Agrawal, R. & Srikant, R. (1994). Fast algorithms for mining association
rules. In *Proc. 20th VLDB Conference*, 487–499.

## See also

[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md),
[`predict_links`](https://mohsaqr.github.io/Nestimate/reference/predict_links.md)

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
