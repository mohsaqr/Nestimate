# Extract Initial Probabilities from Model

`extract_initial_probs()` returns the initial-state probabilities of a
network: the share of sequences that start in each state. It reads a
`netobject` (or any list) that stores them as `initial`,
`initial_probs`, `inits` or `initial_probabilities`, a `tna`-family
object that stores them as `initial`, `initial_probs` or `inits` (models
from the tna package store `inits`), or an `mcml` (one vector per
layer).

Of the
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
methods, the transition methods (`"relative"`, `"frequency"`,
`"attention"`, `"ngram"`, `"gap"` and `"reverse"`) store initial
probabilities, as do
[`wtna()`](https://pak.dynasite.org/Nestimate/reference/wtna.md)
transition networks. For co-occurrence and association networks, and for
a model that keeps them under another name, the function returns a
uniform distribution over the states of the weight matrix, with a
warning: those values are a placeholder, not an estimate.

## Usage

``` r
extract_initial_probs(model)

# S3 method for class 'nest_initial_probs'
summary(object, ...)
```

## Arguments

- model:

  The network: a `netobject`, a `tna`-family model, an `mcml`, or a list
  holding an `initial` element. Input from which neither initial
  probabilities nor a weight matrix can be read raises an error.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `nest_initial_probs`.

- ...:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method:
  ignored.

## Value

For a single network: a named numeric vector of initial-state
probabilities that sums to 1 (stored values are rescaled when their sum
differs from 1 by more than 1e-6), of class
`c("nest_initial_probs", "numeric")`.
[`summary()`](https://rdrr.io/r/base/summary.html) returns a data frame
with one row per state and the columns `state` and `prob`, sorted by
decreasing probability.

For an `mcml`: a named list of plain numeric vectors, `macro` first and
then one per cluster, taken as stored in each layer (`inits`).

## See also

[`extract_transition_matrix`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md)
for extracting the transition matrix,
[`extract_edges`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
for extracting an edge list.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
net <- build_network(seqs, method = "relative")
init_probs <- extract_initial_probs(net)
print(init_probs)
#>         A         B         C 
#> 0.6666667 0.3333333 0.0000000 
#> attr(,"class")
#> [1] "nest_initial_probs" "numeric"           
```
