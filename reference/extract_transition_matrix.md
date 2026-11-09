# Extract Transition Matrix from Model

`extract_transition_matrix()` returns the weight matrix of a network:
for a transition network, the transition probabilities or counts between
states, with rows as the states moved from and columns as the states
moved to. It reads a `netobject` (or any list with a `weights` element),
a `tna`, `ftna`, `ctna` or `atna` model, an `mcml` (one matrix per
layer) or a bare matrix, and can rescale the rows to sum to 1.

## Usage

``` r
extract_transition_matrix(model, type = c("raw", "scaled"))

# S3 method for class 'nest_transition_matrix'
summary(object, ...)
```

## Arguments

- model:

  The network: a `netobject`, a `tna`-family model, an `mcml`, a list
  holding the matrix as `weights` (or, failing that, `transition_matrix`
  or `transition`), or a weight matrix. Anything else raises an error.

- type:

  Which matrix to return. `"raw"` (default) returns the weights as
  stored. `"scaled"` divides each row by its sum, so the rows of a count
  matrix become transition probabilities; a row that sums to zero is
  left as it is.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `nest_transition_matrix`.

- ...:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method:
  ignored.

## Value

For a single network: a square numeric matrix with the model's state
names (if any) as row and column names, of class
`c("nest_transition_matrix", "matrix", "array")`. It behaves as an
ordinary matrix in arithmetic and indexing;
[`summary()`](https://rdrr.io/r/base/summary.html) returns a data frame
with one row per non-zero entry (self-transitions included) and the
columns `from`, `to` and `weight`, sorted by decreasing absolute weight.

For an `mcml`: a named list of such matrices, `macro` first and then one
per cluster.

## Details

Models store their weights in different places; the function looks in
each in turn. `"scaled"` is useful when the stored weights are counts or
otherwise do not sum to 1 per row.

## See also

[`extract_initial_probs`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md)
for extracting initial probabilities,
[`extract_edges`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
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
