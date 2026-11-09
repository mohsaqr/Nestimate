# Cluster sequence data (deprecated alias)

`cluster_data()` is the old name of
[`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
renamed in Nestimate 0.4.3. It is kept so older tutorials and the
historical pkgdown reference still run: every call raises a deprecation
warning and passes its arguments to
[`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
unchanged.

## Usage

``` r
cluster_data(...)
```

## Arguments

- ...:

  Arguments passed unchanged to
  [`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md).

## Value

The `net_clustering` object that
[`build_clusters()`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
returns.

## See also

[`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
[`cluster_network`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md),
[`cluster_mmm`](https://pak.dynasite.org/Nestimate/reference/cluster_mmm.md).
