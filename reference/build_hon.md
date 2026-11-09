# Build a Higher-Order Network (HON)

`build_hon()` builds a higher-order network from sequences with the
BuildHON algorithm (Xu, Wickramarathne & Chawla, 2016) or its faster
variant BuildHON+ (Saebi et al., 2020), and returns a `net_hon`. A node
is a state together with as much of its preceding history as changes
where the sequence goes next: `"B"` is a first-order node, and
`"A -> B"` is state B reached from A, kept as a node of its own when the
distribution of the next state after A then B differs enough from the
one after B alone (see *Details*), or when it begins a longer history
that is kept, such as `"A -> B -> C"`. It reads sequences in a wide data
frame or character matrix, a list of state vectors, or a `netobject`,
`cograph_network` or `tna` built from sequences.

An edge weight is the probability of moving from the source node to the
target node. A higher-order node marks a dependence on earlier states
that a first-order network, with one node per state, averages away.
Whether a history is extended is decided by a rule, not a test: the
Kullback-Leibler divergence between the longer and the shorter history's
next-state distributions must exceed a threshold that grows with the
order and shrinks with how often the longer history was observed, so the
same divergence passes more easily for a history observed more often,
and nodes resting on a few observations are noisy (raise `min_freq` to
set them aside). To test whether the sequences need memory beyond the
first order, use
[`markov_order_test()`](https://pak.dynasite.org/Nestimate/reference/markov_order_test.md);
[`build_mogen()`](https://pak.dynasite.org/Nestimate/reference/build_mogen.md)
selects an order by likelihood, and
[`pathways()`](https://pak.dynasite.org/Nestimate/reference/pathways.md)
lists the higher-order paths found here.

## Usage

``` r
build_hon(
  data,
  max_order = 5L,
  min_freq = 1L,
  collapse_repeats = FALSE,
  method = "hon+"
)

# S3 method for class 'net_hon'
print(x, ...)

# S3 method for class 'net_hon'
summary(object, ...)
```

## Arguments

- data:

  The sequences. A data frame or character matrix has one row per
  sequence and one column per position; a data frame of several columns
  must hold only position columns named as a numbered run (`T1`, `T2`,
  ...), or be a sequence table returned by a Nestimate verb such as
  [`extract_sequences()`](https://pak.dynasite.org/Nestimate/reference/extract_sequences.md),
  and any other frame raises an error of class `nestimate_bad_input`
  (declare its positions with `build_network(state_cols = )` and pass
  the network). A plain list holds one vector of states per sequence. A
  `netobject`, a bare `cograph_network` or a `tna` built from sequences
  contributes its stored sequences, with a `tna`'s numeric state codes
  converted to their labels. Values are read as character. A missing
  value ends a run of observed states, so a sequence with an internal
  gap becomes two sequences, and sequences with fewer than two states
  are dropped; when none is left, an error is raised. A one-hot or
  association network, a `netobject_group`, a classed list, and a matrix
  or data frame of only 0 and 1 raise an error of class
  `nestimate_bad_input`.

- max_order:

  The highest order a node may reach, that is, how many states of
  history it carries. Default `5`. Lower orders result where the data do
  not support a longer history;
  [`print()`](https://rdrr.io/r/base/print.html) reports the highest
  order reached beside the one requested. Must be at least 1; a
  fractional value is truncated to a whole number.

- min_freq:

  The minimum count for a transition to be used. Transitions observed
  fewer than `min_freq` times are treated as zero, also in the totals
  that set the divergence threshold. Default `1`. Must be at least 1; a
  fractional value is truncated to a whole number.

- collapse_repeats:

  Whether to collapse each run of a repeated state into one state before
  counting (after splitting at missing values). Default `FALSE`.

- method:

  The algorithm. `"hon+"` (default) is BuildHON+, which counts a longer
  history only when it is needed and skips an extension when even the
  largest divergence the history could reach would not pass the
  threshold. `"hon"` is the original BuildHON, which first counts the
  transitions from every history of up to `max_order` states. Both stop
  at `max_order`.

- x:

  For the [`print()`](https://rdrr.io/r/base/print.html) method: a
  `net_hon`.

- ...:

  Ignored by the [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods.

- object:

  For the [`summary()`](https://rdrr.io/r/base/summary.html) method: a
  `net_hon`.

## Value

A `net_hon` (class `c("net_hon", "cograph_network")`). Print it for the
number of nodes, edges and sequences and the orders reached;
[`summary()`](https://rdrr.io/r/base/summary.html) prints the same with
the first-order states and the number of nodes of each order, and
returns the network's edge table (integer `from`/`to` node indices and
`weight`).
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
gives the weight of each ordered pair of nodes by node name (pairs with
zero weight are included at its default `threshold = 0`, and self-loops
only with `include_self = TRUE`),
[`pathways()`](https://pak.dynasite.org/Nestimate/reference/pathways.md)
the higher-order paths,
[`build_honem()`](https://pak.dynasite.org/Nestimate/reference/build_honem.md)
an embedding of the nodes and `build_simplicial(type = "pathway")` a
simplicial complex; it draws with
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).
The [`print()`](https://rdrr.io/r/base/print.html) method returns its
input invisibly. The fields:

- weights, matrix:

  The weighted adjacency matrix (rows are sources, columns targets)
  under two names: `weights` is the `cograph_network` field, `matrix`
  the higher-order name kept for compatibility. Rows and columns carry
  the node names in arrow notation (e.g. `"A -> B"`).

- ho_edges:

  The higher-order edge table, one row per edge: `path` (the full state
  sequence, e.g. `"A -> B -> C"`), `from` (the source node, its
  history), `to` (the next state), `count` (how often the transition was
  observed), `probability` (the transition probability), and
  `from_order` and `to_order` (the orders of the source and target
  nodes).

- edges:

  The `cograph_network` edge table: one row per non-zero cell of
  `weights`, with integer `from`/`to` indices into `nodes` and a numeric
  `weight`. It is not the arrow-notation table, which is `ho_edges`.

- nodes:

  A data frame with columns `id`, `label` and `name`, one row per node;
  `label` and `name` hold the node names in arrow notation.

- n_nodes:

  The number of nodes.

- n_edges:

  The number of rows of `ho_edges`.

- first_order_states:

  The distinct states of the sequences analysed (after splitting at
  missing values, collapsing repeats and dropping short sequences),
  sorted.

- max_order_requested:

  The `max_order` used.

- max_order_observed:

  The highest order reached.

- min_freq:

  The `min_freq` used.

- n_trajectories:

  The number of sequences analysed, after splitting at missing values
  and dropping short ones.

- directed:

  Always `TRUE`.

- meta:

  The `cograph_network` metadata list (`source`, `layout`, and
  `tna$method = "hon"`).

- node_groups:

  Always `NULL`.

## Details

**Node names**: a first-order node is the state itself, `"A"`. A
second-order node for "came from A, now at B" is `"A -> B"`; a
third-order one is `"A -> B -> C"`, and so on.

**Algorithm**:

1.  Count the transitions from every history of up to `max_order` states
    (subsequences of up to `max_order + 1` states); BuildHON+ counts the
    longer histories only when they are needed.

2.  Turn the counts into next-state distributions, setting counts below
    `min_freq` to zero.

3.  For each first-order source, extend the history one state at a time
    and keep the longer history when the Kullback-Leibler divergence
    (base 2) of its distribution from the currently kept one exceeds
    `order / log2(1 + n)`, where `order` is the new order and `n` the
    number of transitions counted from the longer history (after setting
    counts below `min_freq` to zero).

4.  Build the network from the kept histories, rewiring edges so the
    higher-order nodes are connected to each other.

The implementation follows the reference Python code of the two papers
(pyHON and pyHON+).

## References

Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
higher-order dependencies in networks. *Science Advances*, 2(5),
e1600028.

Saebi, M., Xu, J., Kaplan, L. M., Ribeiro, B., & Chawla, N. V. (2020).
Efficient modeling of higher-order dependencies in networks: from
algorithm to application for anomaly detection. *EPJ Data Science*,
9(1), 15.

## Examples

``` r
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
hon <- build_hon(seqs, max_order = 2)

# \donttest{
# From list of trajectories
trajs <- list(
  c("A", "B", "C", "D", "A"),
  c("A", "B", "D", "C", "A"),
  c("A", "B", "C", "D", "A")
)
hon <- build_hon(trajs, max_order = 3, min_freq = 1)
print(hon)
#> Higher-Order Network (HON)
#>   Nodes:        4 (4 first-order states)
#>   Edges:        7
#>   Max order:    1 (requested 3)
#>   Min freq:     1
#>   Trajectories: 3
summary(hon)
#> Higher-Order Network (HON) Summary
#>   Nodes: 4 | Edges: 7 | Trajectories: 3
#>   First-order states: A, B, C, D
#>   Max order observed: 1 (requested: 3)
#>   Min frequency: 1
#>   Node order distribution:
#>     Order 1: 4 nodes
#>   from to    weight
#> 1    3  1 0.3333333
#> 2    4  1 0.6666667
#> 3    1  2 1.0000000
#> 4    2  3 0.6666667
#> 5    4  3 0.3333333
#> 6    2  4 0.3333333
#> 7    3  4 0.6666667

# From data.frame (rows = trajectories)
df <- data.frame(T1 = c("A", "A"), T2 = c("B", "B"),
                 T3 = c("C", "D"), T4 = c("D", "C"))
hon <- build_hon(df, max_order = 2)
# }
```
