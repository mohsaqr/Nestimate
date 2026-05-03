# Build a Co-occurrence Network

Constructs an undirected co-occurrence network from various input
formats. Entities that appear together in the same transaction,
document, or record are connected, with edge weights reflecting raw
counts or a similarity measure. Argument names follow the citenets
convention.

## Usage

``` r
cooccurrence(
  data,
  field = NULL,
  by = NULL,
  sep = NULL,
  similarity = c("none", "jaccard", "cosine", "inclusion", "association", "dice",
    "equivalence", "relative"),
  threshold = 0,
  min_occur = 1L,
  diagonal = TRUE,
  top_n = NULL,
  ...
)
```

## Arguments

- data:

  Input data. Accepts:

  - A `data.frame` with a delimited column (`field` + `sep`).

  - A `data.frame` in long/bipartite format (`field` + `by`).

  - A binary (0/1) `data.frame` or `matrix` (auto-detected).

  - A wide sequence `data.frame` or `matrix` (non-binary).

  - A `list` of character vectors (each element is a transaction).

- field:

  Character. The entity column — determines what the nodes are. For
  delimited format, a single column whose values are split by `sep`. For
  long/bipartite format, the item column. For multi-column delimited, a
  vector of column names whose split values are pooled per row.

- by:

  Character or `NULL`. What links the nodes. For long/bipartite format,
  the grouping column (e.g., `"paper_id"`, `"session_id"`). Each unique
  value of `by` defines one transaction. If `NULL` (default), entities
  co-occur within the same row/document.

- sep:

  Character or `NULL`. Separator for splitting delimited fields (e.g.,
  `";"`, `","`). Default `NULL`.

- similarity:

  Character. Similarity measure applied to the raw co-occurrence counts.
  One of:

  `"none"`

  :   Raw co-occurrence counts.

  `"jaccard"`

  :   \\C\_{ij} / (f_i + f_j - C\_{ij})\\.

  `"cosine"`

  :   \\C\_{ij} / \sqrt{f_i \cdot f_j}\\ (Salton's cosine).

  `"inclusion"`

  :   \\C\_{ij} / \min(f_i, f_j)\\ (Simpson coefficient).

  `"association"`

  :   \\C\_{ij} / (f_i \cdot f_j)\\ (association strength /
      probabilistic affinity index; van Eck & Waltman, 2009).

  `"dice"`

  :   \\2 C\_{ij} / (f_i + f_j)\\.

  `"equivalence"`

  :   \\C\_{ij}^2 / (f_i \cdot f_j)\\ (Salton's cosine squared).

  `"relative"`

  :   Row-normalized: each row sums to 1.

- threshold:

  Numeric. Minimum edge weight to retain. Edges below this value are set
  to zero. Applied *after* similarity normalization. Default 0.

- min_occur:

  Integer. Minimum entity frequency (number of transactions an entity
  must appear in). Entities below this threshold are dropped before
  computing co-occurrence. Default 1 (keep all).

- diagonal:

  Logical. If `TRUE` (default), the diagonal of the co-occurrence matrix
  is kept (item self-co-occurrence = item frequency). If `FALSE`, the
  diagonal is zeroed.

- top_n:

  Integer or `NULL`. If specified, return only the top `top_n` edges by
  weight. Default `NULL` (all edges).

- ...:

  Currently unused.

## Value

A `netobject` (undirected) with `method = "co_occurrence_fn"`. The
`$weights` matrix contains similarity (or raw) co-occurrence values. The
`$params` list stores the similarity method, threshold, and the number
of transactions.

## Details

Six input formats are supported, auto-detected from the combination of
`field`, `by`, and `sep`:

1.  **Delimited**: `field` + `sep` (single column). Each cell is split
    by `sep`, trimmed, and de-duplicated per row.

2.  **Multi-column delimited**: `field` (vector) + `sep`. Values from
    multiple columns are split, pooled, and de-duplicated per row.

3.  **Long bipartite**: `field` + `by`. Groups by `by`; unique values of
    `field` within each group form a transaction.

4.  **Binary matrix**: No `field`/`by`/`sep`, all values 0/1. Columns
    are items, rows are transactions.

5.  **Wide sequence**: No `field`/`by`/`sep`, non-binary. Unique values
    across each row form a transaction.

6.  **List**: A plain list of character vectors.

The pipeline converts all formats into a list of character vectors
(transactions), optionally filters by `min_occur`, builds a binary
transaction matrix, computes `crossprod(B)` for the raw co-occurrence
counts, normalizes via the chosen `similarity`, then applies `threshold`
and `top_n` filtering.

## References

van Eck, N. J., & Waltman, L. (2009). How to normalize co-occurrence
data? An analysis of some well-known similarity measures. *Journal of
the American Society for Information Science and Technology*, 60(8),
1635–1651.

## See also

[`build_cna`](https://mohsaqr.github.io/Nestimate/reference/build_cna.md)
for sequence-positional co-occurrence via
[`build_network()`](https://mohsaqr.github.io/Nestimate/reference/build_network.md).

## Examples

``` r
# Delimited field (e.g., keyword co-occurrence)
df <- data.frame(
  id = 1:4,
  keywords = c("network; graph", "graph; matrix; network",
               "matrix; algebra", "network; algebra; graph")
)
net <- cooccurrence(df, field = "keywords", sep = ";")

# Long/bipartite
long_df <- data.frame(
  paper = c(1, 1, 1, 2, 2, 3, 3),
  keyword = c("network", "graph", "matrix", "graph", "algebra",
              "network", "algebra")
)
net <- cooccurrence(long_df, field = "keyword", by = "paper")

# List of transactions
transactions <- list(c("A", "B"), c("B", "C"), c("A", "B", "C"))
net <- cooccurrence(transactions, similarity = "jaccard")

# Binary matrix
bin <- matrix(c(1,0,1, 1,1,0, 0,1,1), nrow = 3, byrow = TRUE,
              dimnames = list(NULL, c("X", "Y", "Z")))
net <- cooccurrence(bin)
```
