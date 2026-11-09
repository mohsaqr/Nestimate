# Window-based Transition Network Analysis

`wtna()` builds a network from one-hot data, where each row is a time
point and each code column holds 1 when that code is active, by counting
codes across windows of consecutive rows. It returns a directed
transition network (codes active in one window, followed by codes active
in the next), an undirected co-occurrence network (codes active in the
same window), or both. With `actor`, each actor's rows are counted on
their own and the counts summed, so no window spans two actors. The
networks are the ones `build_network(method = "wtna")` and
`build_network(method = "wtna_cooccurrence")` estimate.

The weights count pairs of active codes, not events: a row with several
active codes adds to several edges at once, and a larger window pairs
more rows. Compare networks built with the same `window_size` and
`mode`. How far the weights move under resampling of actors is checked
with
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md),
and two networks are compared with
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md);
both need `actor`, since rows without it form a single unit.

## Usage

``` r
wtna(
  data,
  method = c("transition", "cooccurrence", "both"),
  type = c("frequency", "relative"),
  codes = NULL,
  window_size = 3L,
  mode = c("non-overlapping", "overlapping"),
  actor = NULL
)

# S3 method for class 'wtna_mixed'
print(x, ...)
```

## Arguments

- data:

  A data frame or matrix whose rows are time points, in order, holding
  the one-hot code columns and, optionally, the `actor` column. A column
  named in `codes` may hold only 0, 1 or missing; a missing cell counts
  as inactive, and any other value raises an error of class
  `nestimate_bad_input`. Under the default `codes = NULL` such a column
  is not read as a code.

- method:

  The network to build: `"transition"` (default), directed, from the
  codes of one window to those of the next; `"cooccurrence"`,
  undirected, between the codes of the same window; or `"both"`, which
  returns the two networks together as a `wtna_mixed`.

- type:

  The weights: `"frequency"` (default) keeps the counts, and
  `"relative"` divides each row by its sum (a row without counts stays
  zero). On a co-occurrence network `"relative"` gives an asymmetric
  matrix, the co-occurrence conditional on the row's code, so keep
  `"frequency"` when symmetric counts are needed.

- codes:

  The one-hot columns, as names, as column positions, or as a single
  `"first:last"` string naming a range of adjacent columns. Default
  `NULL`, every numeric or logical column holding only 0, 1 and missing
  values, except `actor`. At least two codes are needed, and a name that
  is not a column stops `wtna()`.

- window_size:

  The number of consecutive rows in a window, a whole number of at
  least 1. Default `3`. With `1` there is no windowing: transitions are
  counted from each row to the next, and co-occurrence within each row,
  and `mode` is ignored.

- mode:

  How windows are formed when `window_size` is above 1.
  `"non-overlapping"` (default) cuts the rows into consecutive, separate
  blocks of `window_size` rows. `"overlapping"` first replaces the rows
  by one row per sliding window of `window_size` rows (step 1), marking
  a code active when it is active in any of them, so n rows give n -
  `window_size` + 1 windows, and then counts those in blocks of
  `window_size`, as `tna::import_onehot(window_type = "sliding")` does.

- actor:

  The column, or columns, naming who each row belongs to. The counts are
  computed per actor (per combination of values, for several columns)
  and summed. Default `NULL`, all rows as one sequence. A missing
  identifier raises an error of class `nestimate_bad_input`, and a name
  that is not a column stops `wtna()`. With several columns, the
  initial-state distribution groups the rows by the first column only.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: a
  `wtna_mixed` object.

- ...:

  For the [`print()`](https://rdrr.io/r/base/print.html) method:
  ignored.

## Value

For `method = "transition"` or `"cooccurrence"`: a `netobject` (class
`c("netobject", "cograph_network")`, see
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md))
with `method` set to `"wtna_transition"` or `"wtna_cooccurrence"`.
[`print()`](https://rdrr.io/r/base/print.html) shows the weights and,
for transitions, the initial probabilities;
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
gives one row per edge,
[`extract_initial_probs()`](https://pak.dynasite.org/Nestimate/reference/extract_initial_probs.md)
the initial-state distribution of a transition network,
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws it,
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
resamples its actors, and
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
compares it with a second network by permuting actors. It is directed
only for transitions, holds the input as given in `data`, and records
the windowing settings (`type`, `window_size`, `mode`, `codes`, `actor`)
in `params`. A transition network also stores the initial-state
distribution in `initial`: for each actor (the first `actor` column when
there are several), the codes active in its first row with an active
code share its weight equally, and the shares are averaged over actors
(without `actor`, the first row with an active code in the whole data).

For `method = "both"`: a `wtna_mixed` object, a list of the `transition`
and `cooccurrence` networks (each a `netobject` as above) and
`method = "wtna_both"`. Its
[`print()`](https://rdrr.io/r/base/print.html) gives the node and edge
counts of each network and returns its input invisibly;
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
and
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
test both networks.

## Details

**Transitions**: with `window_size = 1`, `crossprod(X[-n, ], X[-1, ])`
counts how often code i is active in row t and code j in row t + 1. With
larger windows, every row of a window is paired with every row of the
next window, which is the outer product of the two windows' column sums;
a final, shorter window is counted as the successor of the last full
one.

**Co-occurrence**: with `window_size = 1`, `crossprod(X)` counts the
codes active in the same row, so the diagonal is how often each code is
active. With larger windows, every row of a window is paired with every
row of the same window, itself included, which is the outer product of
the window's column sums.

**Windows**: in non-overlapping windows a code counts once for each row
it is active in. These counts follow tna's windowed algorithm.

**Per actor**: with `actor`, the counts are computed per group and
summed.

## See also

[`build_network`](https://pak.dynasite.org/Nestimate/reference/build_network.md)

## Examples

``` r
oh <- matrix(c(1,0,0, 0,1,0, 0,0,1, 1,0,0), nrow = 4, byrow = TRUE,
             dimnames = list(NULL, c("A","B","C")))
w <- wtna(oh)

# \donttest{
# Simple one-hot data
df <- data.frame(
  A = c(1, 0, 1, 0, 1),
  B = c(0, 1, 0, 1, 0),
  C = c(0, 0, 1, 0, 0)
)

# Transition network
net <- wtna(df)
print(net)
#> Network (method: wtna_transition) [directed]
#>   Weights: [1.000, 2.000]  |  mean: 1.333
#> 
#>   Weight matrix:
#>     A B C
#>   A 2 2 0
#>   B 1 1 0
#>   C 1 1 0 
#> 
#>   Initial probabilities:
#>   A             1.000  ████████████████████████████████████████
#>   B             0.000  
#>   C             0.000  

# Both networks
nets <- wtna(df, method = "both")
print(nets)   # the transition and the co-occurrence network
#> Mixed Window TNA (transition + co-occurrence)
#> -- Transition (directed) --
#>   Nodes: 3  |  Edges: 6
#> -- Co-occurrence (undirected) --
#>   Nodes: 3  |  Edges: 6

# With windowing
net <- wtna(df, window_size = 2, mode = "non-overlapping")
# }
```
