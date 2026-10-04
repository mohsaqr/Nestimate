# Cluster sequence data (deprecated alias)

Renamed to
[`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md)
in Nestimate 0.4.3. This thin wrapper is preserved so the function name
in older tutorials and the historical pkgdown reference continues to
work; it issues a one-shot deprecation warning and forwards every
argument unchanged.

## Usage

``` r
cluster_data(...)
```

## Arguments

- ...:

  Passed verbatim to
  [`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md).

## Value

The `net_clustering` object returned by
[`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md).

## See also

[`build_clusters`](https://pak.dynasite.org/Nestimate/reference/build_clusters.md),
[`cluster_network`](https://pak.dynasite.org/Nestimate/reference/cluster_network.md),
[`cluster_mmm`](https://pak.dynasite.org/Nestimate/reference/cluster_mmm.md).
