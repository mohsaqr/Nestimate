# Build Multi-Order Generative Model (MOGen)

Constructs higher-order De Bruijn graphs from sequential trajectory data
and selects the optimal Markov order using AIC, BIC, or likelihood ratio
tests.

## Usage

``` r
build_mogen(
  data,
  max_order = 5L,
  criterion = c("aic", "bic", "lrt"),
  lrt_alpha = 0.01
)
```

## Arguments

- data:

  A data.frame (rows = trajectories, columns = time points), a list of
  character/numeric vectors (one per trajectory), a `tna` object, or a
  `netobject` with sequence data. For `tna`/`netobject`, numeric state
  IDs are automatically converted to label names.

- max_order:

  Integer. Maximum Markov order to test (default 5). Must be a whole
  number; a non-integer value (e.g. `2.7`) is an error rather than being
  silently truncated. A `max_order` at or above the longest trajectory
  is capped at (longest path - 1) with a message.

- criterion:

  Character. Model selection criterion: `"aic"` (default), `"bic"`, or
  `"lrt"` (likelihood ratio test).

- lrt_alpha:

  Numeric. Significance threshold for LRT (default 0.01).

## Value

An object of class `c("net_mogen", "cograph_network")` with components:

- optimal_order:

  Selected optimal Markov order.

- criterion:

  Which criterion was used for selection.

- orders:

  Integer vector of tested orders (0 to max_order, after any capping).

- aic:

  Named numeric vector of AIC values per order.

- bic:

  Named numeric vector of BIC values per order.

- log_likelihood:

  Named numeric vector of log-likelihoods.

- dof:

  Named integer vector of cumulative DOF per model.

- layer_dof:

  Named integer vector of per-layer DOF.

- transition_matrices:

  List of row-stochastic transition matrices (index 1 = order 0, held as
  the marginal named numeric vector).

- count_matrices:

  List of the matching raw count matrices, same indexing; read by
  [`mogen_transitions()`](https://saqr.me/Nestimate/reference/mogen_transitions.md).

- states:

  Unique first-order states.

- n_paths:

  Number of trajectories.

- n_observations:

  Total number of state observations.

- weights:

  `cograph_network` weight matrix: the transition matrix of the selected
  optimal order (a 1 x n matrix named `"marginal"` when the optimal
  order is 0). Its dimnames are the internal k-gram keys (states joined
  by a non-printing separator), not arrow notation.

- nodes:

  data.frame (`id`, `label`, `name`) of the optimal-order De Bruijn
  nodes.

- edges:

  `cograph_network` edge data.frame with integer `from`/`to` node
  indices and a numeric `weight`. The readable arrow-notation table is
  [`mogen_transitions()`](https://saqr.me/Nestimate/reference/mogen_transitions.md).

- directed:

  Logical. Always `TRUE`.

- n_nodes, n_edges:

  Counts for the optimal-order graph.

- meta:

  `cograph_network` metadata list.

- node_groups:

  Always `NULL`.

## Details

At order k, nodes are k-tuples of states and edges represent transitions
between overlapping k-tuples. The model tests increasingly complex
Markov orders and selects the one that best balances fit and parsimony.

## References

Scholtes, I. (2017). When is a Network a Network? Multi-Order Graphical
Model Selection in Pathways and Temporal Networks. *KDD 2017*.

Gote, C., Casiraghi, G., Schweitzer, F., & Scholtes, I. (2023).
Predicting variable-length paths in networked systems using multi-order
generative models. *Applied Network Science*, 8, 68.

## Examples

``` r
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
mg <- build_mogen(seqs, max_order = 2)

# \donttest{
trajs <- list(c("A","B","C","D"), c("A","B","D","C"),
              c("B","C","D","A"), c("C","D","A","B"))
m <- build_mogen(trajs, max_order = 3)
print(m)
#> Multi-Order Generative Model (MOGen)
#>   Optimal order:  2 (by aic)
#>   Orders tested:  0 to 3
#>   States:         4
#>   Paths:          4 (16 observations)
#>   AIC:            50.4 | 28.7 | 26.7* | 26.7
#>   BIC:            52.7 | 32.6 | 31.3* | 31.3
#>   (* = minimum;  AIC and BIC agree on order 2)
plot(m)

# }
```
