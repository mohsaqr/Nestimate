# Extract Transition Matrix from Model

Extract the transition probability matrix from a TNA model object.

## Usage

``` r
extract_transition_matrix(model, type = c("raw", "scaled"))
```

## Arguments

- model:

  A TNA model object or a list containing a 'weights' element.

- type:

  Character. Type of matrix to return:

  "raw"

  :   The raw weight matrix as stored in the model.

  "scaled"

  :   Row-normalized to ensure rows sum to 1.

  Default: "raw".

## Value

A square numeric matrix with row and column names as state names.

## Details

TNA models store transition weights in different locations depending on
the model type. This function handles the extraction automatically.

For "scaled" type, each row is divided by its sum to create valid
transition probabilities. This is useful when the original weights don't
sum to 1.

## See also

[`extract_initial_probs`](https://mohsaqr.github.io/Nestimate/reference/extract_initial_probs.md)
for extracting initial probabilities,
[`extract_edges`](https://mohsaqr.github.io/Nestimate/reference/extract_edges.md)
for extracting an edge list.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
net <- build_network(seqs, method = "relative")
trans_mat <- extract_transition_matrix(net)
print(trans_mat)
#>   A         B         C
#> A 0 0.3333333 0.6666667
#> B 1 0.0000000 0.0000000
#> C 0 1.0000000 0.0000000
#> attr(,"class")
#> [1] "nest_transition_matrix" "matrix"                 "array"                 
```
