# Nestimate Objects for cograph

Nestimate outputs S3 objects that cograph can plot directly. This document describes the structure of each object type.

## netobject / cograph_network

Class: `c("netobject", "cograph_network")`

Returned by `build_network()`. The dual class means cograph functions (`splot()`, `communities()`, etc.) work directly — no conversion needed.

```
$weights        numeric matrix (p x p), named rows/cols = node labels
$nodes          data.frame: id (int), label (chr), name (chr), x (num), y (num)
$edges          data.frame: from (int node ID), to (int node ID), weight (num)
$directed       logical
$meta           list: source, layout, tna (list with $method)
$node_groups    NULL or factor vector
$method         character: "relative", "frequency", "co_occurrence", "cor", "pcor", "glasso", "ising", "attention"
$n_nodes        integer
$n_edges        integer
```

Node IDs in `$edges$from` / `$edges$to` are integers matching `$nodes$id`. Labels are in `$nodes$label`.

### Variants

- **`netobject_ml`** — from `build_network(..., level = "both")`. Contains `$between` and `$within`, each a `netobject`.
- **`netobject_group`** — from grouped `build_network()` or `cluster_data()`. Named list of `netobject`s, one per group.

## net_bootstrap

Class: `"net_bootstrap"`

Returned by `bootstrap_network()`. All matrices are p x p with named rows/cols matching node labels.

```
$original       netobject/cograph_network — the input network
$model          netobject/cograph_network — pruned network (non-significant edges zeroed)
$mean           numeric matrix — mean bootstrap weights
$sd             numeric matrix — standard deviation of bootstrap weights
$p_values       numeric matrix — edge-level p-values (proportion of bootstrap samples where edge is zero)
$significant    numeric matrix — 1 where significant, 0 otherwise
$ci_lower       numeric matrix — lower confidence interval bound
$ci_upper       numeric matrix — upper confidence interval bound
$cr_lower       numeric matrix — lower consistency range bound
$cr_upper       numeric matrix — upper consistency range bound
$summary        data.frame: from, to, observed, mean, sd, ci_lower, ci_upper, p_value, significant
$method         character — estimation method used
$iter           integer — number of bootstrap iterations
$ci_level       numeric — confidence level (default 0.95)
$inference      character — "stability" or "threshold"
$edge_threshold numeric — threshold used (if inference = "threshold")
```

**For cograph plotting**: use `$original` for the observed network, `$model` for the pruned (significant-only) network. Both are `cograph_network` objects ready for `splot()`.

## net_permutation

Class: `"net_permutation"`

Returned by `permutation_test()`. Compares two networks edge by edge. All matrices are p x p with named rows/cols.

```
$x              netobject/cograph_network — first network
$y              netobject/cograph_network — second network
$diff           numeric matrix — raw difference (x - y) per edge
$diff_sig       numeric matrix — difference with non-significant edges zeroed
$p_values       numeric matrix — permutation p-values per edge
$effect_size    numeric matrix — Cohen's d effect size per edge
$summary        data.frame: from, to, weight_x, weight_y, diff, effect_size, p_value, sig
$method         character — estimation method
$iter           integer — number of permutations
$alpha          numeric — significance level
$paired         logical — paired or unpaired test
$adjust         character — p-value adjustment method
```

**For cograph plotting**: use `$x` and `$y` as the two networks. Use `$diff_sig` as a weight matrix for a difference network (significant differences only).

## boot_glasso

Class: `"boot_glasso"`

Returned by `boot_glasso()`. Specialized bootstrap for EBICglasso networks.

```
$pcor           numeric matrix — estimated partial correlation network
$ci             list of matrices — edge CI bounds ($lower, $upper per edge)
$cs_results     list — centrality stability results (CS-coefficient per measure)
$edge_diff_p    numeric matrix or NULL — edge difference test p-values
$cent_diff_p    list or NULL — centrality difference test p-values
$iter           integer
$alpha          numeric
```

**For cograph plotting**: use `$pcor` directly as a weight matrix for `splot()`.

## Common patterns

All weight matrices share the same structure:
- Square numeric matrix, p x p
- Named rows and columns (node labels)
- Directed networks: `mat[i,j]` = edge from node i to node j
- Undirected networks: symmetric matrix

All `netobject` / `cograph_network` objects can be passed directly to cograph's `splot()` — the `$weights`, `$nodes`, `$edges`, and `$directed` fields are in cograph's expected format.
