# Validate a netobject / cograph_network against the shared schema

Enforces the structural contract that both Nestimate netobjects and
psychnet objects must satisfy to be interchangeable across the package
boundary. This is the single place that says what "a network object"
means, so a drift on either side (a renamed field, a mistyped edge
column) fails loudly here rather than mis-rendering three layers
downstream.

## Usage

``` r
validate_netobject(x)
```

## Arguments

- x:

  An object expected to satisfy the `cograph_network` contract.

## Value

Invisibly `TRUE` if `x` conforms; otherwise stops with the full list of
violations.

## Details

The contract is deliberately the *shared* subset: the `$nodes` `x`/`y`
layout columns and the Nestimate pipeline fields (`$data`, `$level`,
...) are not required, and `$edges` endpoints may be either integer node
indices (Nestimate) or character labels (psychnet).

## See also

[`as_netobject`](https://saqr.me/Nestimate/reference/as_netobject.md)

## Examples

``` r
net <- build_cor(data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50)))
validate_netobject(net)
```
