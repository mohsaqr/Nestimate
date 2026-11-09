# Build a Co-occurrence Network

`cooccurrence()` builds an undirected network of items that appear
together in the same transaction (a document, a session, a record): an
edge joins two items that occur together, weighted by how often they do,
or by a similarity measure that adjusts that count for how often each
item occurs overall. With `similarity = "relative"` the network is
directed. It reads delimited columns (keywords separated by `";"`), a
long table of items per transaction, a 0/1 incidence table, wide
sequences or a list of transactions; the *Details* say how each is
recognised.

Raw counts grow with how frequent the two items are, so the most
frequent items dominate a count network. The similarity measures correct
for that in different ways (van Eck & Waltman, 2009); compare a few
before reading which pairs are strongest. The network is descriptive: it
keeps no source data, so
[`bootstrap_network()`](https://pak.dynasite.org/Nestimate/reference/bootstrap_network.md)
and
[`permutation()`](https://pak.dynasite.org/Nestimate/reference/permutation.md)
cannot resample it.

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
  ...,
  format = c("auto", "binary")
)
```

## Arguments

- data:

  The transactions, in one of these forms:

  - a data frame with a delimited column (`field` and `sep`);

  - a data frame in long (bipartite) form, one row per item in a
    transaction (`field` and `by`);

  - a 0/1 or `TRUE`/`FALSE` data frame or matrix, one row per
    transaction and one column per item, declared with
    `format = "binary"`;

  - a wide sequence data frame or matrix with other values, whose rows
    are the transactions;

  - a plain list of character vectors, one per transaction.

  A network object (`netobject`, `netobject_group`, `cograph_network`,
  `tna`) raises an error of class `nestimate_bad_input`; pass the data
  it was built from.

- field:

  The column holding the items, which become the nodes. With `sep`, its
  values are split by `sep`; several column names pool the split values
  of each row. With `by`, it is the item column of a long table. Default
  `NULL`. Setting `by` or `sep` without `field`, or `field` without
  either, raises an error of class `nestimate_bad_input`, as does naming
  a column that `data` does not have, or setting any of the three for a
  list.

- by:

  The column that defines the transactions of a long table: each of its
  values is one transaction, holding the distinct `field` values of its
  rows (for example `"paper_id"` or `"session_id"`). Default `NULL`,
  where items co-occur within the same row. Ignored when `sep` is set,
  which reads `field` as delimited.

- sep:

  The separator that splits delimited `field` values, such as `";"` or
  `","`, matched literally; the split values are trimmed. Default
  `NULL`. See `field` for the combinations that raise an error.

- similarity:

  How the co-occurrence count \\C\_{ij}\\ of items \\i\\ and \\j\\, with
  item frequencies \\f_i\\ and \\f_j\\, is turned into an edge weight.
  Default `"none"`.

  `"none"`

  :   Raw co-occurrence counts.

  `"jaccard"`

  :   \\C\_{ij} / (f_i + f_j - C\_{ij})\\.

  `"cosine"`

  :   \\C\_{ij} / \sqrt{f_i \cdot f_j}\\ (Salton's cosine).

  `"inclusion"`

  :   \\C\_{ij} / \min(f_i, f_j)\\ (Simpson coefficient).

  `"association"`

  :   \\C\_{ij} / (f_i \cdot f_j)\\ (association strength, or
      probabilistic affinity index; van Eck & Waltman, 2009).

  `"dice"`

  :   \\2 C\_{ij} / (f_i + f_j)\\.

  `"equivalence"`

  :   \\C\_{ij}^2 / (f_i \cdot f_j)\\ (Salton's cosine squared).

  `"relative"`

  :   Each row divided by its sum (the diagonal included when
      `diagonal = TRUE`), so each row sums to 1: the weight from `i` to
      `j` is the share of `i`'s co-occurrences that involve `j`. The
      matrix is no longer symmetric, so the network is directed and both
      directions of a pair are edges.

  A zero denominator is replaced by 1, so an item without co-occurrences
  gets weight 0 rather than `NaN`.

- threshold:

  The smallest weight kept: weights below it, on the scale of
  `similarity`, are set to zero. Default `0` (none).

- min_occur:

  The fewest transactions an item must appear in to be kept. Rarer items
  are dropped before the co-occurrences are counted, and transactions
  left empty are dropped with them. Default `1` (keep all).

- diagonal:

  Whether the diagonal is kept. With `TRUE` (default) it holds each
  item's frequency (the number of transactions it appears in),
  transformed by `similarity`; with `FALSE` it is set to zero. Item
  frequencies for the similarity measures are the same either way.

- top_n:

  How many of the strongest off-diagonal edges to keep: pairs in an
  undirected network, directed cells with `similarity = "relative"`.
  Edges tied with the `top_n`-th weight are all kept, so a few more may
  remain, and the diagonal is not affected. Default `NULL` (all edges).

- ...:

  Not used.

- format:

  How 0/1 columns are read. Columns of only 0 and 1 (or `TRUE` and
  `FALSE`) could be items or recoded positions, so with `"auto"`
  (default) they raise an error of class `nestimate_bad_input`.
  `"binary"` declares them items, one row per transaction, where a
  missing cell counts as no occurrence; it raises the same error on data
  with any other value.

## Value

A `netobject` (class `c("netobject", "cograph_network")`), undirected
(directed with `similarity = "relative"`), with
`method = "co_occurrence_fn"`. Print it for an overview;
[`extract_edges()`](https://pak.dynasite.org/Nestimate/reference/extract_edges.md)
gives one row per edge,
[`extract_transition_matrix()`](https://pak.dynasite.org/Nestimate/reference/extract_transition_matrix.md)
the weight matrix (one row and column per item kept),
[`net_centrality()`](https://pak.dynasite.org/Nestimate/reference/net_centrality.md)
the node centralities, and
[`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html)
draws it. It is built without
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md)
and holds no source data (`data` is `NULL`). Its `params` record
`similarity`, `threshold`, `min_occur`, `diagonal`, `top_n`,
`n_transactions` and `n_items`.

## Details

Six input forms are recognised from `data` and from which of `field`,
`by` and `sep` are set:

1.  **Delimited**: `field` (one column) and `sep`. Each cell is split by
    `sep`, trimmed, and its distinct values form a transaction.

2.  **Multi-column delimited**: `field` (several columns) and `sep`. The
    split values of a row's columns are pooled into one transaction.

3.  **Long bipartite**: `field` and `by`. The distinct `field` values
    within each `by` group form a transaction.

4.  **Binary**: none of the three, all values 0/1 or logical, and
    `format = "binary"`. Columns are items, rows are transactions; a
    missing cell is no occurrence. Without `format = "binary"` such data
    is refused, since 0/1 columns could also be recoded positions.

5.  **Wide sequence**: none of the three, other values. The distinct
    states of each row form a transaction, leaving out missing values
    and the void markers (`"%"`, `"*"`, `""`, `"NA"`, `"NaN"`). A data
    frame must be a single column or a numbered run of position columns
    (`T1`, `T2`, ...) with nothing beside it; any other frame raises an
    error of class `nestimate_bad_input`, and a column called `action`
    and repeated column names are refused the same way.

6.  **List**: a plain list of character vectors.

To count co-occurrence of states within sequences through
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md),
with its resampling verbs, use `build_network(method = "co_occurrence")`
or
[`build_cna()`](https://pak.dynasite.org/Nestimate/reference/build_cna.md)
instead.

All forms are turned into a list of transactions, which `min_occur`
filters. The raw counts are `crossprod(B)` of the binary transaction
matrix `B`; `similarity` is applied to them, then `threshold`, then
`top_n`.

## References

van Eck, N. J., & Waltman, L. (2009). How to normalize co-occurrence
data? An analysis of some well-known similarity measures. *Journal of
the American Society for Information Science and Technology*, 60(8),
1635–1651.

## See also

[`build_cna`](https://pak.dynasite.org/Nestimate/reference/build_cna.md)
for sequence-positional co-occurrence via
[`build_network()`](https://pak.dynasite.org/Nestimate/reference/build_network.md).

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
net <- cooccurrence(bin, format = "binary")
```
