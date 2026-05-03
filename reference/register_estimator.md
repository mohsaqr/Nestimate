# Register a Network Estimator

Register a custom or built-in network estimator function by name.
Estimators registered here can be used by
[`estimate_network`](https://mohsaqr.github.io/Nestimate/reference/estimate_network.md)
via the `method` parameter.

## Usage

``` r
register_estimator(name, fn, description, directed)
```

## Arguments

- name:

  Character. Unique name for the estimator (e.g. `"relative"`,
  `"glasso"`).

- fn:

  Function. The estimator function. Must accept `data` as its first
  argument and `...` for additional parameters. Must return a list with
  at least: `matrix` (square numeric matrix), `nodes` (character
  vector), `directed` (logical).

- description:

  Character. Short description of the estimator.

- directed:

  Logical. Whether the estimator produces directed networks.

## Value

Invisible `NULL`.

## See also

[`get_estimator`](https://mohsaqr.github.io/Nestimate/reference/get_estimator.md),
[`list_estimators`](https://mohsaqr.github.io/Nestimate/reference/list_estimators.md),
[`remove_estimator`](https://mohsaqr.github.io/Nestimate/reference/remove_estimator.md),
[`estimate_network`](https://mohsaqr.github.io/Nestimate/reference/estimate_network.md)

## Examples

``` r
my_fn <- function(data, ...) {
  m <- cor(data)
  diag(m) <- 0
  list(matrix = m, nodes = colnames(m), directed = FALSE)
}
register_estimator("my_cor", my_fn, "Custom correlation", directed = FALSE)
df <- data.frame(A = rnorm(20), B = rnorm(20), C = rnorm(20))
net <- build_network(df, method = "my_cor")
remove_estimator("my_cor")
```
