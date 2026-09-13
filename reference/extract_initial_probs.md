# Extract Initial Probabilities from Model

Extract the initial state probability vector from a TNA model object.

## Usage

``` r
extract_initial_probs(model)
```

## Arguments

- model:

  A `netobject`, TNA model object, `mcml` object, or a list containing
  an `initial` element.

## Value

For a single network: a named numeric vector of initial state
probabilities summing to 1, of class
`c("nest_initial_probs", "numeric")`. The class stamp only adds a
[`summary()`](https://rdrr.io/r/base/summary.html) method returning a
tidy `state`/`prob` data frame.

For an `mcml` object: a named list of such vectors, with `macro` first
and then one element per cluster (taken from each layer's `$inits`,
unstamped).

## Details

Initial probabilities represent the probability of starting a sequence
in each state. If the model doesn't have explicit initial probabilities,
this function falls back to a uniform distribution over the states of
the transition matrix and warns.

## See also

[`extract_transition_matrix`](https://saqr.me/Nestimate/reference/extract_transition_matrix.md)
for extracting the transition matrix,
[`extract_edges`](https://saqr.me/Nestimate/reference/extract_edges.md)
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
