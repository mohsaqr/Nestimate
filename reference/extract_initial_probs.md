# Extract Initial Probabilities from Model

Extract the initial state probability vector from a TNA model object.

## Usage

``` r
extract_initial_probs(model)
```

## Arguments

- model:

  A TNA model object or a list containing an 'initial' element.

## Value

A named numeric vector of initial state probabilities.

## Details

Initial probabilities represent the probability of starting a sequence
in each state. If the model doesn't have explicit initial probabilities,
this function attempts to estimate them from the data or use uniform
probabilities.

## See also

[`extract_transition_matrix`](https://mohsaqr.github.io/Nestimate/reference/extract_transition_matrix.md)
for extracting the transition matrix,
[`extract_edges`](https://mohsaqr.github.io/Nestimate/reference/extract_edges.md)
for extracting an edge list.

## Examples

``` r
# \donttest{
seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
net <- build_network(seqs, method = "relative")
init_probs <- extract_initial_probs(net)
print(init_probs)
#>         A         B         C 
#> 0.6666667 0.3333333 0.0000000 
# }
```
