# Estimate a window transition network analysis (WTNA) model

`build_network(method = "wtna")`, or for short `build_wtna()`, estimates
a window transition network analysis (WTNA) model from one-hot codes.
The rows are cut into windows of `window_size` rows, and the edge from
code A to code B counts the times B is active in the window after one
where A is active.

Windows do not cross actors (`actor`) or sessions (`session`).
`build_network(method = "wtna_cooccurrence")` counts the codes active in
the same window instead;
[`wtna()`](https://pak.dynasite.org/Nestimate/reference/wtna.md) returns
both networks at once.

## Usage

``` r
build_wtna(
  data,
  codes = NULL,
  window_size = 3L,
  mode = c("non-overlapping", "overlapping"),
  type = c("frequency", "relative"),
  actor = NULL,
  session = NULL,
  metadata_cols = NULL,
  group = NULL,
  scaling = NULL,
  threshold = 0,
  labels = NULL,
  params = list(),
  ...
)
```

## Arguments

- data:

  One-hot data: one row per time point and one column of 0 and 1 per
  code.

- codes:

  The code columns (one per state, 1 when the state is active). Default
  `NULL`: every column outside the other roles.

- window_size:

  How many consecutive rows form one window. Default `3`.

- mode:

  `"non-overlapping"` (default) or `"overlapping"` windows: overlapping
  ones slide by one row, as
  `tna::import_onehot(window_type = "sliding")`, and are counted in
  blocks of `window_size`.

- type:

  `"frequency"` (default) keeps the counts; `"relative"` divides each
  row by its sum.

- actor:

  The column naming who each row belongs to; windows do not cross
  actors. Default `NULL`.

- session:

  The column identifying sessions; windows do not cross sessions.
  Default `NULL`.

- metadata_cols:

  Columns to keep as the network's metadata, never read as codes.
  Default `NULL`.

- group:

  A column to split the data by: one network per value, returned as a
  `netobject_group`. Default `NULL`.

- scaling:

  Rescaling applied to the weights after estimation, in the order given:
  `"minmax"` rescales the weights to run from 0 to 1, over every cell
  for the transition methods (tna's rule, so with no zero cell the
  smallest weight becomes 0) and over the non-zero weights, into (0,
  1\], for the other methods, `"max"` divides by the largest absolute
  weight, `"rank"` replaces the non-zero weights by their ranks, and
  `"normalize"` makes each row's absolute weights sum to 1; an
  undirected network whose matrix it makes asymmetric becomes directed,
  so both of a pair's values are kept as edges. Zero weights stay zero
  under every option. Combine them as in `c("rank", "minmax")`. Default
  `NULL` (none).

- threshold:

  Weights whose absolute value is below this are set to zero. Default
  `0` (none).

- labels:

  New display labels for the nodes, applied after estimation: a
  two-column data frame (`name`, `label`), a named character vector
  `c(name = "label")` or a named list. It rewrites the node labels and
  the weight matrix's dimnames; nodes not listed keep their names.
  Default `NULL`.

- params:

  Settings for the estimator, as a named list. The network keeps them.
  Default [`list()`](https://rdrr.io/r/base/list.html).

- ...:

  Further arguments for
  [`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).
  `method` cannot be given (error of class `nestimate_bad_input`).

## Value

A `netobject` (a `netobject_group` with `group`), as
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
returns it.

## See also

[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
for every method;
[`wtna()`](https://pak.dynasite.org/Nestimate/reference/wtna.md) for the
transition and co-occurrence networks together;
[`action_to_onehot()`](https://pak.dynasite.org/Nestimate/reference/action_to_onehot.md)
to turn an event log into one-hot codes.

## Examples

``` r
codes <- c("Reading", "Video", "Forum", "Quiz", "Coding", "Review")
build_network(learning_activities, method = "wtna", codes = codes,
              actor = "student")
#> Window TNA (transitions) [directed]
#>   Weights: [2493.000, 4314.000]  |  mean: 3025.056
#> 
#>   Strongest edges (10 of 36):
#>    from    to      weight
#>    Video   Video   4314  
#>    Reading Reading 4297  
#>    Review  Review  4206  
#>    Coding  Coding  4005  
#>    Forum   Forum   3881  
#>    Quiz    Quiz    3695  
#>    Review  Video   3113  
#>    Coding  Video   3094  
#>    Video   Review  3018  
#>    Reading Review  3008   
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix.
build_wtna(learning_activities, codes = codes, actor = "student")
#> Window TNA (transitions) [directed]
#>   Weights: [2493.000, 4314.000]  |  mean: 3025.056
#> 
#>   Strongest edges (10 of 36):
#>    from    to      weight
#>    Video   Video   4314  
#>    Reading Reading 4297  
#>    Review  Review  4206  
#>    Coding  Coding  4005  
#>    Forum   Forum   3881  
#>    Quiz    Quiz    3695  
#>    Review  Video   3113  
#>    Coding  Video   3094  
#>    Video   Review  3018  
#>    Reading Review  3008   
#> 
#>   extract_edges(include_self = TRUE) lists every edge; extract_transition_matrix() gives the weight matrix.
```
