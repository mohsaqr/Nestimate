# Summary Method for boot_glasso

Summary Method for boot_glasso

## Usage

``` r
# S3 method for class 'boot_glasso'
summary(object, type = "edges", ...)
```

## Arguments

- object:

  A `boot_glasso` object.

- type:

  Character. Summary type: `"edges"` (default), `"centrality"`, `"cs"`,
  `"predictability"`, or `"all"`.

- ...:

  Additional arguments (ignored).

## Value

A data frame or list of data frames depending on `type`.

## Examples

``` r
set.seed(1)
dat <- as.data.frame(matrix(rnorm(60), ncol = 3))
bg <- boot_glasso(dat, iter = 10, cs_iter = 5, centrality = "strength")
summary(bg, type = "edges")
#>       edge weight   ci_lower  ci_upper inclusion
#> 1 V1 -- V2      0 -0.4571610 0.3104397       0.5
#> 2 V1 -- V3      0  0.0000000 0.5580800       0.8
#> 3 V2 -- V3      0 -0.4873358 0.0000000       0.6
# \donttest{
set.seed(42)
mat <- matrix(rnorm(60), ncol = 4)
colnames(mat) <- LETTERS[1:4]
boot <- boot_glasso(as.data.frame(mat), iter = 20, cs_iter = 10,
  centrality = "strength", seed = 42)
summary(boot, type = "edges")
#>     edge weight    ci_lower  ci_upper inclusion
#> 1 A -- B      0 -0.11158199 0.2237173      0.35
#> 2 A -- C      0 -0.35185504 0.0000000      0.45
#> 3 B -- C      0 -0.30946880 0.2065199      0.25
#> 4 A -- D      0  0.00000000 0.3570903      0.50
#> 5 B -- D      0 -0.15048304 0.3105727      0.35
#> 6 C -- D      0 -0.03724162 0.1787669      0.15
# }
```
