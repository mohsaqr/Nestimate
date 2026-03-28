# Compute Node Predictability

Computes the proportion of variance explained (R\\^2\\) for each node in
the network, following Haslbeck & Waldorp (2018).

For `method = "glasso"` or `"pcor"`, predictability is computed
analytically from the precision matrix: \$\$R^2_j = 1 - 1 /
\Omega\_{jj}\$\$ where \\\Omega\\ is the precision (inverse correlation)
matrix.

For `method = "cor"`, predictability is the multiple R\\^2\\ from
regressing each node on its network neighbors (nodes with non-zero
edges).

## Usage

``` r
predictability(object, ...)

# S3 method for class 'netobject'
predictability(object, ...)

# S3 method for class 'netobject_ml'
predictability(object, ...)
```

## Arguments

- object:

  A `netobject` or `netobject_ml` object.

- ...:

  Additional arguments (ignored).

## Value

For `netobject`: a named numeric vector of R\\^2\\ values (one per node,
between 0 and 1).

For `netobject_ml`: a list with elements `$between` and `$within`, each
a named numeric vector.

A named numeric vector of predictability values per node.

A list with `within` and `between` predictability vectors.

## References

Haslbeck, J. M. B., & Waldorp, L. J. (2018). How well do network models
predict observations? On the importance of predictability in network
models. *Behavior Research Methods*, 50(2), 853–861.
[doi:10.3758/s13428-017-0910-x](https://doi.org/10.3758/s13428-017-0910-x)

## Examples

``` r
# \donttest{
set.seed(42)
mat <- matrix(rnorm(60), ncol = 4)
colnames(mat) <- LETTERS[1:4]
net <- build_network(as.data.frame(mat), method = "glasso")
predictability(net)
#> A B C D 
#> 0 0 0 0 
# }
```
