# Extract Transition Table from a MOGen Model

`mogen_transitions()` returns the transitions of a fitted
[`build_mogen()`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
model at one order, one row per path, most frequent first, with the path
written out in arrow notation, its observed count and its transition
probability. At order k each edge of the De Bruijn graph is a path of
k + 1 states: at order 2, the edge from node `"AI -> FAIL"` to node
`"FAIL -> SOLVE"` is the path `AI -> FAIL -> SOLVE`, read as the
probability of `SOLVE` after the context `AI -> FAIL`.

A probability is the share of its context's observations that moved to
that state, so in a context seen a few times it rests on those few
observations: read it with its `count`, and use `min_count` to leave
rare transitions out.

## Usage

``` r
mogen_transitions(x, order = NULL, min_count = 1L)
```

## Arguments

- x:

  A `net_mogen` object from
  [`build_mogen()`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md).

- order:

  The order whose transitions are returned, a whole number from 1 to the
  highest order fitted; a non-integer value is an error, not truncated.
  Default `NULL`, the order the model selected. Order 0 has no
  transitions, so when the model selected order 0, give `order`
  explicitly; otherwise `mogen_transitions()` stops.

- min_count:

  The smallest observed count a transition needs to be listed. Default
  `1`, every observed transition. A value below 1 also lists the
  transitions that never occurred, with count and probability 0: every
  pair of observed contexts where the second continues the first (at
  order 1, every pair of states).

## Value

A data frame with one row per listed transition, sorted by `count`
(descending), with columns:

- path:

  The full sequence of states, such as `"AI -> FAIL -> SOLVE"`.

- count:

  How often the transition was observed.

- probability:

  The transition probability P(to \| from), rounded to 4 decimal places.

- from:

  The context: the k states the transition starts from.

- to:

  The next state.

When no transition reaches `min_count`, a data frame with the same
columns and no rows.

## Examples

``` r
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
mg <- build_mogen(seqs, max_order = 2)
mogen_transitions(mg, order = 1)
#>     path count probability from to
#> 1 B -> C     3      1.0000    B  C
#> 2 A -> B     2      1.0000    A  B
#> 3 C -> D     2      0.6667    C  D
#> 4 C -> A     1      0.3333    C  A
#> 5 D -> A     1      1.0000    D  A

# \donttest{
trajs <- list(c("A","B","C","D"), c("A","B","D","C"),
              c("B","C","D","A"), c("C","D","A","B"))
m <- build_mogen(trajs, max_order = 3)
mogen_transitions(m, order = 1)
#>     path count probability from to
#> 1 A -> B     3      1.0000    A  B
#> 2 C -> D     3      1.0000    C  D
#> 3 D -> A     2      0.6667    D  A
#> 4 B -> C     2      0.6667    B  C
#> 5 D -> C     1      0.3333    D  C
#> 6 B -> D     1      0.3333    B  D
# }
```
