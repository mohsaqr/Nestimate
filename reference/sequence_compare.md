# Compare Subsequence Patterns Between Groups

Extracts all k-gram patterns (subsequences of length k) from sequences
in each group, computes standardized residuals against the independence
model, and optionally runs a permutation or chi-square test of group
differences.

## Usage

``` r
sequence_compare(
  x,
  group = NULL,
  sub = 3:5,
  min_freq = 5L,
  test = c("permutation", "chisq", "none"),
  iter = 1000L,
  adjust = "fdr"
)
```

## Arguments

- x:

  A `netobject_group` (from grouped `build_network`), a `netobject`
  (requires `group`), or a wide-format `data.frame` (requires `group`).

- group:

  Character or vector. Column name or vector of group labels. Not needed
  for `netobject_group`.

- sub:

  Integer vector. Pattern lengths to analyze. Default: `3:5`.

- min_freq:

  Integer. Minimum frequency in each group for a pattern to be included.
  Default: 5.

- test:

  Character. Inference method: one of `"permutation"` (default),
  `"chisq"`, or `"none"`. See Details.

- iter:

  Integer. Permutation iterations. Only used when
  `test = "permutation"`. Default: 1000.

- adjust:

  Character. P-value correction method (see
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)). Default:
  `"fdr"`.

## Value

An object of class `"net_sequence_comparison"` containing:

- patterns:

  Tidy data.frame. Always present: `pattern`, `length`, `freq_<group>`,
  `prop_<group>`, `resid_<group>`. If `test = "permutation"`:
  `effect_size`, `p_value`. If `test = "chisq"`: `statistic`, `p_value`.

- groups:

  Character vector of group names.

- n_patterns:

  Integer. Number of patterns passing min_freq.

- params:

  List of sub, min_freq, test, iter, adjust.

## Details

Standardized residuals are always computed from a 2xG contingency table
of (this pattern vs. everything else) using the textbook formula
`(o - e) / sqrt(e * (1 - r/N) * (1 - c/N))`. They describe how much each
group's count for a given pattern deviates from expectation under
independence, scaled to be approximately N(0,1) under the null.

The optional `test` argument chooses an inference method:

- `"permutation"`:

  Shuffles group labels across sequences and recomputes a per-pattern
  statistic (row-wise Euclidean residual norm). Answers: "is this
  pattern's distribution associated with group membership at the *actor*
  level?" Respects the sequence as the unit of analysis; can be
  underpowered when the number of sequences is small.

- `"chisq"`:

  Runs `chisq.test` on the 2xG table per pattern. Answers: "do the group
  *streams* generate this pattern at different rates?" Treats each
  k-gram occurrence as an event; fast and powerful even with few
  sequences, but the iid assumption it makes is optimistic when
  sequences are strongly autocorrelated.

- `"none"`:

  Skip inference. Only residuals, frequencies, and proportions are
  returned.

P-values are adjusted once across all patterns (not per-pattern) using
any method supported by
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). The default is
`"fdr"` (Benjamini-Hochberg).

## Examples

``` r
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 60, TRUE),
  V2 = sample(LETTERS[1:4], 60, TRUE),
  V3 = sample(LETTERS[1:4], 60, TRUE),
  V4 = sample(LETTERS[1:4], 60, TRUE)
)
grp <- rep(c("X", "Y"), 30)
net <- build_network(seqs, method = "relative")
res <- sequence_compare(net, group = grp, sub = 2:3, test = "chisq")
```
