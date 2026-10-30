# Extract Sequences in Wide Format

Returns the sequences behind any Nestimate input as a wide data frame,
one row per sequence and one column per position. Accepts the output of
[`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md),
a long event log, a wide data frame or character matrix, a `netobject`,
a grouped `netobject_group`, or a `tna` model. Model objects return
their stored, decoded states.

## Usage

``` r
extract_sequences(
  x,
  cols = NULL,
  actor = NULL,
  action = NULL,
  time = NULL,
  meta = FALSE,
  ...
)
```

## Arguments

- x:

  Sequence data: a `nestimate_data` object from
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md);
  a long event log (supply `action`); a wide `data.frame` or character
  matrix; a `netobject`, `netobject_group`, `cograph_network` or `tna`
  built from sequence data.

- cols:

  Character vector of state columns for wide data frames. When `NULL`,
  numeric columns of a mixed data frame are treated as covariates and
  the remaining columns as states. Ignored for model objects.

- actor, action, time:

  Column names of a long event log, passed to
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md).
  Supplying `action` marks `x` as long.

- meta:

  Logical. When `TRUE`, the covariate columns that belong to each
  sequence (actor, session, group, metadata) are placed before the state
  columns. Default `FALSE`.

- ...:

  Further arguments passed to
  [`prepare()`](https://pak.dynasite.org/Nestimate/reference/prepare.md)
  for long input (e.g. `time_threshold`, `order`).

## Value

A `data.frame` with one row per sequence. State columns are character,
with `NA` after a sequence ends; their names are the position names of
the input (`T1`, `T2`, ... for prepared or long data). With
`meta = TRUE` the covariate columns come first. A grouped network adds
its grouping column. Input that has no sequence data raises
`nestimate_bad_input`.

## See also

[`prepare`](https://pak.dynasite.org/Nestimate/reference/prepare.md),
[`sequence_patterns`](https://pak.dynasite.org/Nestimate/reference/sequence_patterns.md)

## Examples

``` r
extract_sequences(trajectories)
#>              1          2          3          4          5          6
#> 1       Active Disengaged Disengaged Disengaged     Active     Active
#> 2      Average    Average    Average    Average    Average    Average
#> 3      Average     Active     Active     Active     Active     Active
#> 4       Active     Active     Active     Active     Active    Average
#> 5       Active     Active     Active    Average     Active    Average
#> 6      Average    Average Disengaged    Average Disengaged Disengaged
#> 7       Active    Average    Average    Average    Average    Average
#> 8      Average Disengaged    Average Disengaged Disengaged Disengaged
#> 9       Active     Active     Active     Active     Active     Active
#> 10     Average    Average    Average    Average    Average    Average
#> 11     Average    Average    Average     Active    Average     Active
#> 12     Average    Average    Average    Average    Average Disengaged
#> 13      Active     Active     Active     Active     Active     Active
#> 14      Active     Active     Active     Active     Active     Active
#> 15     Average    Average Disengaged Disengaged    Average Disengaged
#> 16     Average     Active     Active    Average    Average    Average
#> 17      Active Disengaged     Active Disengaged Disengaged Disengaged
#> 18     Average    Average    Average    Average    Average     Active
#> 19      Active    Average Disengaged    Average Disengaged    Average
#> 20     Average    Average    Average    Average    Average    Average
#> 21      Active Disengaged    Average Disengaged    Average Disengaged
#> 22     Average    Average     Active     Active       <NA>       <NA>
#> 23     Average Disengaged Disengaged     Active    Average    Average
#> 24      Active     Active     Active    Average    Average    Average
#> 25      Active     Active     Active     Active     Active     Active
#> 26      Active     Active     Active     Active     Active     Active
#> 27  Disengaged    Average     Active    Average     Active    Average
#> 28  Disengaged    Average    Average    Average Disengaged    Average
#> 29  Disengaged    Average    Average    Average    Average     Active
#> 30  Disengaged Disengaged       <NA>       <NA>       <NA>       <NA>
#> 31  Disengaged    Average     Active    Average    Average    Average
#> 32      Active     Active     Active     Active     Active     Active
#> 33  Disengaged       <NA>       <NA>       <NA>       <NA>       <NA>
#> 34      Active     Active     Active     Active     Active    Average
#> 35  Disengaged Disengaged Disengaged    Average    Average    Average
#> 36     Average    Average    Average    Average Disengaged    Average
#> 37      Active     Active    Average     Active    Average    Average
#> 38     Average    Average Disengaged    Average Disengaged     Active
#> 39      Active    Average     Active     Active     Active     Active
#> 40     Average     Active     Active     Active    Average     Active
#> 41     Average Disengaged Disengaged Disengaged Disengaged     Active
#> 42      Active Disengaged     Active     Active     Active     Active
#> 43      Active    Average    Average    Average    Average    Average
#> 44     Average    Average    Average    Average    Average    Average
#> 45     Average Disengaged Disengaged Disengaged Disengaged Disengaged
#> 46     Average Disengaged Disengaged Disengaged Disengaged Disengaged
#> 47      Active     Active     Active     Active     Active    Average
#> 48  Disengaged Disengaged Disengaged    Average Disengaged Disengaged
#> 49     Average     Active    Average    Average    Average     Active
#> 50      Active     Active     Active     Active     Active     Active
#> 51     Average Disengaged Disengaged Disengaged Disengaged Disengaged
#> 52  Disengaged    Average    Average    Average    Average    Average
#> 53     Average    Average Disengaged Disengaged Disengaged     Active
#> 54  Disengaged Disengaged Disengaged Disengaged Disengaged    Average
#> 55      Active    Average    Average     Active     Active    Average
#> 56     Average    Average    Average    Average    Average    Average
#> 57     Average Disengaged Disengaged    Average    Average Disengaged
#> 58     Average    Average    Average    Average    Average    Average
#> 59  Disengaged    Average    Average    Average Disengaged    Average
#> 60     Average Disengaged Disengaged Disengaged    Average Disengaged
#> 61     Average    Average Disengaged    Average    Average    Average
#> 62  Disengaged     Active    Average    Average    Average    Average
#> 63  Disengaged    Average    Average    Average Disengaged    Average
#> 64      Active Disengaged Disengaged Disengaged       <NA>       <NA>
#> 65     Average    Average    Average    Average    Average Disengaged
#> 66      Active     Active     Active     Active     Active     Active
#> 67     Average    Average    Average    Average    Average    Average
#> 68      Active     Active     Active     Active     Active     Active
#> 69     Average    Average Disengaged    Average    Average    Average
#> 70  Disengaged    Average Disengaged Disengaged Disengaged Disengaged
#> 71  Disengaged    Average    Average     Active     Active    Average
#> 72  Disengaged    Average Disengaged Disengaged Disengaged Disengaged
#> 73      Active    Average    Average     Active     Active    Average
#> 74  Disengaged    Average    Average    Average    Average    Average
#> 75      Active     Active     Active     Active     Active     Active
#> 76      Active     Active    Average    Average     Active    Average
#> 77        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 78      Active     Active     Active    Average    Average     Active
#> 79      Active     Active     Active    Average     Active     Active
#> 80     Average     Active     Active     Active     Active     Active
#> 81     Average     Active    Average     Active    Average     Active
#> 82      Active     Active     Active     Active     Active     Active
#> 83     Average    Average    Average    Average    Average Disengaged
#> 84      Active    Average    Average    Average     Active    Average
#> 85     Average     Active     Active     Active    Average     Active
#> 86      Active     Active     Active     Active     Active     Active
#> 87      Active     Active     Active    Average     Active     Active
#> 88     Average     Active    Average     Active     Active Disengaged
#> 89      Active     Active     Active     Active     Active     Active
#> 90     Average     Active     Active     Active    Average     Active
#> 91     Average    Average    Average     Active    Average    Average
#> 92      Active     Active     Active     Active     Active     Active
#> 93     Average    Average    Average Disengaged    Average    Average
#> 94      Active     Active     Active     Active     Active     Active
#> 95  Disengaged    Average Disengaged     Active    Average    Average
#> 96  Disengaged    Average    Average    Average Disengaged    Average
#> 97  Disengaged       <NA>       <NA>       <NA>       <NA>       <NA>
#> 98  Disengaged       <NA>       <NA>       <NA>       <NA>       <NA>
#> 99     Average Disengaged Disengaged Disengaged       <NA>       <NA>
#> 100 Disengaged    Average     Active    Average    Average    Average
#> 101 Disengaged Disengaged Disengaged Disengaged Disengaged Disengaged
#> 102 Disengaged       <NA>       <NA>       <NA>       <NA>       <NA>
#> 103 Disengaged Disengaged       <NA>       <NA>       <NA>       <NA>
#> 104       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 105 Disengaged Disengaged Disengaged    Average       <NA>       <NA>
#> 106 Disengaged       <NA>       <NA>       <NA>       <NA>       <NA>
#> 107     Active    Average    Average     Active     Active    Average
#> 108     Active    Average    Average     Active     Active    Average
#> 109 Disengaged    Average    Average    Average Disengaged    Average
#> 110    Average    Average    Average    Average    Average Disengaged
#> 111     Active    Average Disengaged    Average Disengaged    Average
#> 112     Active    Average    Average     Active     Active    Average
#> 113    Average    Average Disengaged Disengaged Disengaged     Active
#> 114 Disengaged    Average     Active    Average     Active    Average
#> 115    Average     Active     Active     Active     Active     Active
#> 116     Active    Average    Average    Average    Average    Average
#> 117     Active Disengaged    Average Disengaged    Average Disengaged
#> 118 Disengaged    Average    Average    Average    Average     Active
#> 119     Active     Active     Active     Active     Active    Average
#> 120     Active     Active     Active     Active     Active    Average
#> 121 Disengaged    Average     Active    Average    Average    Average
#> 122     Active Disengaged     Active     Active     Active     Active
#> 123    Average    Average Disengaged    Average Disengaged Disengaged
#> 124     Active     Active     Active     Active     Active     Active
#> 125    Average     Active     Active     Active    Average     Active
#> 126     Active     Active     Active     Active     Active    Average
#> 127    Average    Average    Average    Average    Average Disengaged
#> 128 Disengaged Disengaged Disengaged    Average Disengaged Disengaged
#> 129 Disengaged    Average     Active    Average     Active    Average
#> 130    Average    Average    Average    Average Disengaged    Average
#> 131    Average     Active    Average    Average    Average     Active
#> 132    Average     Active     Active    Average    Average    Average
#> 133     Active     Active     Active     Active     Active     Active
#> 134    Average    Average    Average    Average    Average Disengaged
#> 135     Active     Active     Active     Active     Active     Active
#> 136     Active Disengaged     Active     Active     Active     Active
#> 137     Active     Active    Average    Average     Active    Average
#> 138 Disengaged    Average     Active    Average    Average    Average
#>              7          8          9         10         11         12
#> 1       Active     Active    Average     Active    Average     Active
#> 2      Average    Average    Average     Active    Average    Average
#> 3       Active     Active     Active    Average     Active     Active
#> 4      Average     Active     Active     Active    Average     Active
#> 5       Active     Active     Active Disengaged     Active     Active
#> 6      Average    Average    Average Disengaged    Average Disengaged
#> 7      Average    Average    Average Disengaged    Average    Average
#> 8       Active     Active     Active     Active     Active     Active
#> 9      Average     Active    Average    Average     Active     Active
#> 10      Active    Average    Average Disengaged    Average    Average
#> 11      Active     Active     Active    Average    Average     Active
#> 12     Average Disengaged    Average    Average Disengaged Disengaged
#> 13     Average    Average     Active     Active     Active     Active
#> 14      Active     Active     Active     Active     Active     Active
#> 15     Average Disengaged    Average    Average Disengaged Disengaged
#> 16     Average    Average    Average     Active    Average    Average
#> 17  Disengaged     Active Disengaged    Average    Average     Active
#> 18      Active    Average    Average    Average Disengaged    Average
#> 19      Active     Active    Average    Average    Average    Average
#> 20      Active    Average    Average    Average Disengaged    Average
#> 21     Average    Average    Average    Average    Average    Average
#> 22        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 23     Average     Active    Average     Active    Average    Average
#> 24     Average    Average    Average       <NA>       <NA>       <NA>
#> 25     Average     Active     Active    Average    Average    Average
#> 26      Active     Active Disengaged     Active     Active     Active
#> 27      Active    Average Disengaged     Active     Active    Average
#> 28  Disengaged    Average Disengaged    Average    Average    Average
#> 29  Disengaged    Average Disengaged    Average    Average     Active
#> 30        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 31  Disengaged    Average Disengaged    Average    Average    Average
#> 32      Active     Active     Active     Active     Active     Active
#> 33        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 34  Disengaged    Average Disengaged    Average     Active    Average
#> 35     Average       <NA>       <NA>       <NA>       <NA>       <NA>
#> 36  Disengaged    Average Disengaged     Active    Average    Average
#> 37  Disengaged     Active     Active     Active     Active     Active
#> 38     Average     Active    Average     Active    Average    Average
#> 39      Active     Active     Active    Average    Average     Active
#> 40     Average     Active     Active    Average    Average     Active
#> 41      Active Disengaged    Average    Average     Active    Average
#> 42      Active     Active     Active     Active     Active     Active
#> 43     Average    Average     Active     Active     Active     Active
#> 44     Average    Average    Average     Active    Average    Average
#> 45     Average Disengaged    Average Disengaged    Average Disengaged
#> 46     Average Disengaged    Average Disengaged Disengaged Disengaged
#> 47     Average    Average    Average     Active     Active    Average
#> 48  Disengaged Disengaged Disengaged Disengaged Disengaged Disengaged
#> 49      Active Disengaged Disengaged     Active    Average Disengaged
#> 50      Active Disengaged Disengaged    Average    Average     Active
#> 51  Disengaged Disengaged Disengaged Disengaged Disengaged Disengaged
#> 52     Average Disengaged Disengaged    Average    Average     Active
#> 53     Average Disengaged Disengaged    Average    Average Disengaged
#> 54  Disengaged     Active Disengaged Disengaged Disengaged    Average
#> 55     Average Disengaged Disengaged    Average    Average Disengaged
#> 56     Average Disengaged Disengaged Disengaged    Average Disengaged
#> 57  Disengaged Disengaged Disengaged Disengaged Disengaged    Average
#> 58     Average    Average Disengaged    Average Disengaged     Active
#> 59     Average     Active    Average    Average Disengaged    Average
#> 60  Disengaged    Average    Average Disengaged Disengaged    Average
#> 61  Disengaged     Active     Active    Average Disengaged     Active
#> 62      Active     Active    Average    Average    Average    Average
#> 63     Average    Average    Average    Average    Average    Average
#> 64        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 65     Average    Average    Average Disengaged Disengaged    Average
#> 66     Average     Active    Average     Active    Average     Active
#> 67  Disengaged    Average    Average Disengaged    Average    Average
#> 68      Active     Active    Average    Average    Average     Active
#> 69     Average    Average    Average    Average Disengaged    Average
#> 70  Disengaged    Average    Average Disengaged Disengaged Disengaged
#> 71     Average    Average    Average    Average    Average    Average
#> 72  Disengaged    Average    Average Disengaged Disengaged    Average
#> 73     Average    Average    Average     Active     Active     Active
#> 74     Average     Active     Active     Active    Average     Active
#> 75      Active     Active     Active     Active     Active     Active
#> 76      Active     Active    Average     Active     Active     Active
#> 77        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 78      Active    Average     Active     Active     Active     Active
#> 79      Active     Active     Active     Active     Active     Active
#> 80      Active     Active     Active     Active     Active     Active
#> 81      Active    Average     Active     Active     Active     Active
#> 82      Active     Active     Active     Active     Active     Active
#> 83     Average     Active    Average    Average    Average     Active
#> 84     Average    Average     Active    Average    Average     Active
#> 85      Active    Average    Average     Active     Active     Active
#> 86      Active     Active     Active     Active     Active     Active
#> 87      Active     Active     Active     Active     Active     Active
#> 88     Average    Average    Average    Average    Average       <NA>
#> 89      Active     Active     Active     Active     Active     Active
#> 90      Active    Average     Active    Average     Active     Active
#> 91     Average    Average    Average    Average     Active     Active
#> 92      Active     Active     Active     Active     Active     Active
#> 93      Active     Active    Average    Average    Average     Active
#> 94      Active     Active     Active     Active     Active     Active
#> 95      Active    Average    Average     Active    Average    Average
#> 96  Disengaged    Average Disengaged    Average Disengaged    Average
#> 97        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 98        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 99        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 100    Average    Average    Average    Average    Average    Average
#> 101 Disengaged    Average       <NA>       <NA>       <NA>       <NA>
#> 102       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 103       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 104       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 105       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 106       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 107    Average Disengaged Disengaged    Average    Average Disengaged
#> 108    Average    Average    Average     Active     Active     Active
#> 109 Disengaged    Average Disengaged    Average    Average    Average
#> 110    Average     Active    Average    Average    Average     Active
#> 111     Active     Active    Average    Average    Average    Average
#> 112    Average Disengaged Disengaged    Average    Average Disengaged
#> 113    Average Disengaged Disengaged    Average    Average Disengaged
#> 114     Active    Average Disengaged     Active     Active    Average
#> 115     Active     Active     Active     Active     Active     Active
#> 116    Average    Average     Active     Active     Active     Active
#> 117    Average    Average    Average    Average    Average    Average
#> 118 Disengaged    Average Disengaged    Average    Average     Active
#> 119    Average     Active     Active     Active    Average     Active
#> 120    Average    Average    Average     Active     Active    Average
#> 121 Disengaged    Average Disengaged    Average    Average    Average
#> 122     Active     Active     Active     Active     Active     Active
#> 123    Average    Average    Average Disengaged    Average Disengaged
#> 124     Active     Active     Active     Active     Active     Active
#> 125    Average     Active     Active    Average    Average     Active
#> 126    Average     Active     Active     Active    Average     Active
#> 127    Average    Average    Average Disengaged Disengaged    Average
#> 128 Disengaged Disengaged Disengaged Disengaged Disengaged Disengaged
#> 129     Active    Average Disengaged     Active     Active    Average
#> 130 Disengaged    Average Disengaged     Active    Average    Average
#> 131     Active Disengaged Disengaged     Active    Average Disengaged
#> 132    Average    Average    Average     Active    Average    Average
#> 133    Average     Active     Active    Average    Average    Average
#> 134    Average    Average    Average Disengaged Disengaged    Average
#> 135    Average     Active    Average    Average     Active     Active
#> 136     Active     Active     Active     Active     Active     Active
#> 137     Active     Active    Average     Active     Active     Active
#> 138    Average    Average    Average    Average    Average    Average
#>             13         14         15
#> 1         <NA>       <NA>       <NA>
#> 2      Average    Average    Average
#> 3       Active     Active     Active
#> 4      Average     Active    Average
#> 5       Active    Average    Average
#> 6      Average    Average    Average
#> 7      Average    Average    Average
#> 8       Active     Active     Active
#> 9       Active     Active     Active
#> 10     Average    Average    Average
#> 11      Active     Active    Average
#> 12     Average    Average    Average
#> 13      Active     Active     Active
#> 14      Active     Active    Average
#> 15  Disengaged Disengaged    Average
#> 16     Average    Average    Average
#> 17      Active     Active     Active
#> 18     Average Disengaged     Active
#> 19      Active    Average     Active
#> 20  Disengaged     Active    Average
#> 21     Average    Average    Average
#> 22        <NA>       <NA>       <NA>
#> 23     Average    Average    Average
#> 24        <NA>       <NA>       <NA>
#> 25      Active     Active    Average
#> 26      Active    Average     Active
#> 27     Average Disengaged Disengaged
#> 28     Average    Average Disengaged
#> 29      Active    Average     Active
#> 30        <NA>       <NA>       <NA>
#> 31     Average    Average Disengaged
#> 32      Active     Active     Active
#> 33        <NA>       <NA>       <NA>
#> 34     Average    Average Disengaged
#> 35        <NA>       <NA>       <NA>
#> 36     Average Disengaged     Active
#> 37      Active     Active    Average
#> 38     Average    Average    Average
#> 39     Average    Average     Active
#> 40     Average     Active    Average
#> 41      Active Disengaged     Active
#> 42      Active     Active     Active
#> 43      Active     Active    Average
#> 44     Average    Average    Average
#> 45  Disengaged Disengaged     Active
#> 46  Disengaged Disengaged    Average
#> 47     Average    Average    Average
#> 48  Disengaged Disengaged Disengaged
#> 49  Disengaged Disengaged Disengaged
#> 50      Active     Active     Active
#> 51  Disengaged Disengaged Disengaged
#> 52      Active     Active     Active
#> 53  Disengaged Disengaged Disengaged
#> 54     Average    Average    Average
#> 55  Disengaged Disengaged Disengaged
#> 56  Disengaged Disengaged Disengaged
#> 57     Average    Average    Average
#> 58      Active     Active     Active
#> 59     Average    Average    Average
#> 60     Average    Average    Average
#> 61      Active     Active     Active
#> 62     Average    Average    Average
#> 63     Average    Average    Average
#> 64        <NA>       <NA>       <NA>
#> 65     Average    Average    Average
#> 66      Active     Active     Active
#> 67     Average    Average    Average
#> 68      Active     Active     Active
#> 69     Average    Average    Average
#> 70     Average Disengaged       <NA>
#> 71     Average    Average    Average
#> 72     Average    Average    Average
#> 73     Average Disengaged     Active
#> 74     Average Disengaged    Average
#> 75     Average    Average    Average
#> 76      Active Disengaged     Active
#> 77        <NA>       <NA>       <NA>
#> 78     Average    Average     Active
#> 79      Active     Active     Active
#> 80      Active     Active     Active
#> 81      Active     Active     Active
#> 82      Active     Active     Active
#> 83     Average Disengaged    Average
#> 84     Average    Average    Average
#> 85     Average     Active     Active
#> 86      Active     Active     Active
#> 87      Active     Active    Average
#> 88        <NA>       <NA>       <NA>
#> 89      Active     Active    Average
#> 90      Active    Average     Active
#> 91      Active    Average    Average
#> 92      Active    Average    Average
#> 93     Average Disengaged     Active
#> 94      Active     Active     Active
#> 95      Active     Active    Average
#> 96     Average    Average       <NA>
#> 97        <NA>       <NA>       <NA>
#> 98        <NA>       <NA>       <NA>
#> 99        <NA>       <NA>       <NA>
#> 100    Average    Average    Average
#> 101       <NA>       <NA>       <NA>
#> 102       <NA>       <NA>       <NA>
#> 103       <NA>       <NA>       <NA>
#> 104       <NA>       <NA>       <NA>
#> 105       <NA>       <NA>       <NA>
#> 106       <NA>       <NA>       <NA>
#> 107 Disengaged Disengaged Disengaged
#> 108    Average Disengaged     Active
#> 109    Average    Average Disengaged
#> 110    Average Disengaged    Average
#> 111     Active    Average     Active
#> 112 Disengaged Disengaged Disengaged
#> 113 Disengaged Disengaged Disengaged
#> 114    Average Disengaged Disengaged
#> 115     Active     Active     Active
#> 116     Active     Active    Average
#> 117    Average    Average    Average
#> 118     Active    Average     Active
#> 119    Average     Active    Average
#> 120    Average    Average    Average
#> 121    Average    Average Disengaged
#> 122     Active     Active     Active
#> 123    Average    Average    Average
#> 124     Active     Active     Active
#> 125    Average     Active    Average
#> 126    Average     Active    Average
#> 127    Average    Average    Average
#> 128 Disengaged Disengaged Disengaged
#> 129    Average Disengaged Disengaged
#> 130    Average Disengaged     Active
#> 131 Disengaged Disengaged Disengaged
#> 132    Average    Average    Average
#> 133     Active     Active    Average
#> 134    Average    Average    Average
#> 135     Active     Active     Active
#> 136     Active     Active     Active
#> 137     Active Disengaged     Active
#> 138    Average    Average    Average

extract_sequences(group_regulation_long, actor = "Actor",
                  action = "Action", time = "Time", meta = TRUE)
#>      .session_id .session_label Actor Achiever Group Course         T1
#> 1           1 s1           1 s1     1     High     1      A   cohesion
#> 2           2 s1           2 s1     2     High     1      A       plan
#> 3           3 s1           3 s1     3     High     1      A  consensus
#> 4           4 s1           4 s1     4     High     1      A    monitor
#> 5           5 s1           5 s1     5     High     1      A    discuss
#> 6           6 s1           6 s1     6     High     1      A       plan
#> 7           7 s1           7 s1     7     High     1      A       plan
#> 8           8 s1           8 s1     8     High     1      A       plan
#> 9           9 s1           9 s1     9     High     1      A    discuss
#> 10         10 s1          10 s1    10     High     1      A    emotion
#> 11         11 s1          11 s1    11     High     2      A    monitor
#> 12         12 s1          12 s1    12     High     2      A  consensus
#> 13         13 s1          13 s1    13     High     2      A    discuss
#> 14         14 s1          14 s1    14     High     2      A    monitor
#> 15         15 s1          15 s1    15     High     2      A      adapt
#> 16         16 s1          16 s1    16     High     2      A       plan
#> 17         17 s1          17 s1    17     High     2      A    discuss
#> 18         18 s1          18 s1    18     High     2      A   cohesion
#> 19         19 s1          19 s1    19     High     2      A    monitor
#> 20         20 s1          20 s1    20     High     2      A   cohesion
#> 21         21 s1          21 s1    21     High     3      A       plan
#> 22         22 s1          22 s1    22     High     3      A    monitor
#> 23         23 s1          23 s1    23     High     3      A   cohesion
#> 24         24 s1          24 s1    24     High     3      A       plan
#> 25         25 s1          25 s1    25     High     3      A       plan
#> 26         26 s1          26 s1    26     High     3      A       plan
#> 27         27 s1          27 s1    27     High     3      A  synthesis
#> 28         28 s1          28 s1    28     High     3      A    monitor
#> 29         29 s1          29 s1    29     High     3      A  consensus
#> 30         30 s1          30 s1    30     High     3      A    emotion
#> 31         31 s1          31 s1    31     High     4      A    emotion
#> 32         32 s1          32 s1    32     High     4      A    monitor
#> 33         33 s1          33 s1    33     High     4      A  consensus
#> 34         34 s1          34 s1    34     High     4      A       plan
#> 35         35 s1          35 s1    35     High     4      A    monitor
#> 36         36 s1          36 s1    36     High     4      A    discuss
#> 37         37 s1          37 s1    37     High     4      A    discuss
#> 38         38 s1          38 s1    38     High     4      A    discuss
#> 39         39 s1          39 s1    39     High     4      A   cohesion
#> 40         40 s1          40 s1    40     High     4      A  consensus
#> 41         41 s1          41 s1    41     High     5      A  consensus
#> 42         42 s1          42 s1    42     High     5      A    emotion
#> 43         43 s1          43 s1    43     High     5      A  consensus
#> 44         44 s1          44 s1    44     High     5      A    discuss
#> 45         45 s1          45 s1    45     High     5      A    emotion
#> 46         46 s1          46 s1    46     High     5      A    monitor
#> 47         47 s1          47 s1    47     High     5      A    discuss
#> 48         48 s1          48 s1    48     High     5      A       plan
#> 49         49 s1          49 s1    49     High     5      A    discuss
#> 50         50 s1          50 s1    50     High     5      A  consensus
#> 51         51 s1          51 s1    51     High     6      A  consensus
#> 52         52 s1          52 s1    52     High     6      A   cohesion
#> 53         53 s1          53 s1    53     High     6      A    discuss
#> 54         54 s1          54 s1    54     High     6      A       plan
#> 55         55 s1          55 s1    55     High     6      A       plan
#> 56         56 s1          56 s1    56     High     6      A    monitor
#> 57         57 s1          57 s1    57     High     6      A    monitor
#> 58         58 s1          58 s1    58     High     6      A    emotion
#> 59         59 s1          59 s1    59     High     6      A  consensus
#> 60         60 s1          60 s1    60     High     6      A    emotion
#> 61         61 s1          61 s1    61     High     7      A       plan
#> 62         62 s1          62 s1    62     High     7      A  synthesis
#> 63         63 s1          63 s1    63     High     7      A    emotion
#> 64         64 s1          64 s1    64     High     7      A    discuss
#> 65         65 s1          65 s1    65     High     7      A  consensus
#> 66         66 s1          66 s1    66     High     7      A    monitor
#> 67         67 s1          67 s1    67     High     7      A   cohesion
#> 68         68 s1          68 s1    68     High     7      A  consensus
#> 69         69 s1          69 s1    69     High     7      A       plan
#> 70         70 s1          70 s1    70     High     7      A    monitor
#> 71         71 s1          71 s1    71     High     8      A    discuss
#> 72         72 s1          72 s1    72     High     8      A       plan
#> 73         73 s1          73 s1    73     High     8      A   cohesion
#> 74         74 s1          74 s1    74     High     8      A    discuss
#> 75         75 s1          75 s1    75     High     8      A       plan
#> 76         76 s1          76 s1    76     High     8      A   cohesion
#> 77         77 s1          77 s1    77     High     8      A   cohesion
#> 78         78 s1          78 s1    78     High     8      A    emotion
#> 79         79 s1          79 s1    79     High     8      A    monitor
#> 80         80 s1          80 s1    80     High     8      A       plan
#> 81         81 s1          81 s1    81     High     9      A       plan
#> 82         82 s1          82 s1    82     High     9      A    discuss
#> 83         83 s1          83 s1    83     High     9      A    emotion
#> 84         84 s1          84 s1    84     High     9      A   cohesion
#> 85         85 s1          85 s1    85     High     9      A  consensus
#> 86         86 s1          86 s1    86     High     9      A   cohesion
#> 87         87 s1          87 s1    87     High     9      A       plan
#> 88         88 s1          88 s1    88     High     9      A    monitor
#> 89         89 s1          89 s1    89     High     9      A    monitor
#> 90         90 s1          90 s1    90     High     9      A       plan
#> 91         91 s1          91 s1    91     High    10      A  consensus
#> 92         92 s1          92 s1    92     High    10      A    discuss
#> 93         93 s1          93 s1    93     High    10      A  consensus
#> 94         94 s1          94 s1    94     High    10      A       plan
#> 95         95 s1          95 s1    95     High    10      A    discuss
#> 96         96 s1          96 s1    96     High    10      A   cohesion
#> 97         97 s1          97 s1    97     High    10      A    monitor
#> 98         98 s1          98 s1    98     High    10      A   cohesion
#> 99         99 s1          99 s1    99     High    10      A      adapt
#> 100       100 s1         100 s1   100     High    10      A       plan
#> 101       101 s1         101 s1   101     High    11      A    monitor
#> 102       102 s1         102 s1   102     High    11      A    emotion
#> 103       103 s1         103 s1   103     High    11      A    discuss
#> 104       104 s1         104 s1   104     High    11      A    emotion
#> 105       105 s1         105 s1   105     High    11      A  consensus
#> 106       106 s1         106 s1   106     High    11      A    monitor
#> 107       107 s1         107 s1   107     High    11      A    monitor
#> 108       108 s1         108 s1   108     High    11      A    emotion
#> 109       109 s1         109 s1   109     High    11      A    discuss
#> 110       110 s1         110 s1   110     High    11      A    monitor
#> 111       111 s1         111 s1   111     High    12      A    emotion
#> 112       112 s1         112 s1   112     High    12      A   cohesion
#> 113       113 s1         113 s1   113     High    12      A       plan
#> 114       114 s1         114 s1   114     High    12      A  consensus
#> 115       115 s1         115 s1   115     High    12      A       plan
#> 116       116 s1         116 s1   116     High    12      A  consensus
#> 117       117 s1         117 s1   117     High    12      A    monitor
#> 118       118 s1         118 s1   118     High    12      A  consensus
#> 119       119 s1         119 s1   119     High    12      A    discuss
#> 120       120 s1         120 s1   120     High    12      A    monitor
#> 121       121 s1         121 s1   121     High    13      A  synthesis
#> 122       122 s1         122 s1   122     High    13      A    discuss
#> 123       123 s1         123 s1   123     High    13      A   cohesion
#> 124       124 s1         124 s1   124     High    13      A  consensus
#> 125       125 s1         125 s1   125     High    13      A    emotion
#> 126       126 s1         126 s1   126     High    13      A    emotion
#> 127       127 s1         127 s1   127     High    13      A  consensus
#> 128       128 s1         128 s1   128     High    13      A       plan
#> 129       129 s1         129 s1   129     High    13      A       plan
#> 130       130 s1         130 s1   130     High    13      A    emotion
#> 131       131 s1         131 s1   131     High    14      A    emotion
#> 132       132 s1         132 s1   132     High    14      A    emotion
#> 133       133 s1         133 s1   133     High    14      A    monitor
#> 134       134 s1         134 s1   134     High    14      A    discuss
#> 135       135 s1         135 s1   135     High    14      A       plan
#> 136       136 s1         136 s1   136     High    14      A       plan
#> 137       137 s1         137 s1   137     High    14      A  consensus
#> 138       138 s1         138 s1   138     High    14      A   cohesion
#> 139       139 s1         139 s1   139     High    14      A  consensus
#> 140       140 s1         140 s1   140     High    14      A       plan
#> 141       141 s1         141 s1   141     High    15      A       plan
#> 142       142 s1         142 s1   142     High    15      A       plan
#> 143       143 s1         143 s1   143     High    15      A    monitor
#> 144       144 s1         144 s1   144     High    15      A    emotion
#> 145       145 s1         145 s1   145     High    15      A       plan
#> 146       146 s1         146 s1   146     High    15      A  consensus
#> 147       147 s1         147 s1   147     High    15      A  consensus
#> 148       148 s1         148 s1   148     High    15      A       plan
#> 149       149 s1         149 s1   149     High    15      A       plan
#> 150       150 s1         150 s1   150     High    15      A    emotion
#> 151       151 s1         151 s1   151     High    16      A    discuss
#> 152       152 s1         152 s1   152     High    16      A  consensus
#> 153       153 s1         153 s1   153     High    16      A  consensus
#> 154       154 s1         154 s1   154     High    16      A    discuss
#> 155       155 s1         155 s1   155     High    16      A       plan
#> 156       156 s1         156 s1   156     High    16      A       plan
#> 157       157 s1         157 s1   157     High    16      A    emotion
#> 158       158 s1         158 s1   158     High    16      A  consensus
#> 159       159 s1         159 s1   159     High    16      A  consensus
#> 160       160 s1         160 s1   160     High    16      A    emotion
#> 161       161 s1         161 s1   161     High    17      A  synthesis
#> 162       162 s1         162 s1   162     High    17      A    discuss
#> 163       163 s1         163 s1   163     High    17      A    monitor
#> 164       164 s1         164 s1   164     High    17      A   cohesion
#> 165       165 s1         165 s1   165     High    17      A   cohesion
#> 166       166 s1         166 s1   166     High    17      A       plan
#> 167       167 s1         167 s1   167     High    17      A    discuss
#> 168       168 s1         168 s1   168     High    17      A    discuss
#> 169       169 s1         169 s1   169     High    17      A    emotion
#> 170       170 s1         170 s1   170     High    17      A       plan
#> 171       171 s1         171 s1   171     High    18      A    emotion
#> 172       172 s1         172 s1   172     High    18      A    emotion
#> 173       173 s1         173 s1   173     High    18      A    monitor
#> 174       174 s1         174 s1   174     High    18      A       plan
#> 175       175 s1         175 s1   175     High    18      A       plan
#> 176       176 s1         176 s1   176     High    18      A  consensus
#> 177       177 s1         177 s1   177     High    18      A       plan
#> 178       178 s1         178 s1   178     High    18      A    emotion
#> 179       179 s1         179 s1   179     High    18      A       plan
#> 180       180 s1         180 s1   180     High    18      A  consensus
#> 181       181 s1         181 s1   181     High    19      A   cohesion
#> 182       182 s1         182 s1   182     High    19      A   cohesion
#> 183       183 s1         183 s1   183     High    19      A    monitor
#> 184       184 s1         184 s1   184     High    19      A   cohesion
#> 185       185 s1         185 s1   185     High    19      A  consensus
#> 186       186 s1         186 s1   186     High    19      A    emotion
#> 187       187 s1         187 s1   187     High    19      A  consensus
#> 188       188 s1         188 s1   188     High    19      A  consensus
#> 189       189 s1         189 s1   189     High    19      A       plan
#> 190       190 s1         190 s1   190     High    19      A    discuss
#> 191       191 s1         191 s1   191     High    20      A    discuss
#> 192       192 s1         192 s1   192     High    20      A    discuss
#> 193       193 s1         193 s1   193     High    20      A    monitor
#> 194       194 s1         194 s1   194     High    20      A  consensus
#> 195       195 s1         195 s1   195     High    20      A  consensus
#> 196       196 s1         196 s1   196     High    20      A  consensus
#> 197       197 s1         197 s1   197     High    20      A    discuss
#> 198       198 s1         198 s1   198     High    20      A    emotion
#> 199       199 s1         199 s1   199     High    20      A    monitor
#> 200       200 s1         200 s1   200     High    20      A  consensus
#> 201       201 s1         201 s1   201     High    21      A    monitor
#> 202       202 s1         202 s1   202     High    21      A    discuss
#> 203       203 s1         203 s1   203     High    21      A    discuss
#> 204       204 s1         204 s1   204     High    21      A   cohesion
#> 205       205 s1         205 s1   205     High    21      A    monitor
#> 206       206 s1         206 s1   206     High    21      A  consensus
#> 207       207 s1         207 s1   207     High    21      A       plan
#> 208       208 s1         208 s1   208     High    21      A    emotion
#> 209       209 s1         209 s1   209     High    21      A  consensus
#> 210       210 s1         210 s1   210     High    21      A    monitor
#> 211       211 s1         211 s1   211     High    22      A       plan
#> 212       212 s1         212 s1   212     High    22      A  consensus
#> 213       213 s1         213 s1   213     High    22      A    monitor
#> 214       214 s1         214 s1   214     High    22      A  synthesis
#> 215       215 s1         215 s1   215     High    22      A    emotion
#> 216       216 s1         216 s1   216     High    22      A    monitor
#> 217       217 s1         217 s1   217     High    22      A    discuss
#> 218       218 s1         218 s1   218     High    22      A       plan
#> 219       219 s1         219 s1   219     High    22      A       plan
#> 220       220 s1         220 s1   220     High    22      A       plan
#> 221       221 s1         221 s1   221     High    23      A    emotion
#> 222       222 s1         222 s1   222     High    23      A    discuss
#> 223       223 s1         223 s1   223     High    23      A       plan
#> 224       224 s1         224 s1   224     High    23      A    discuss
#> 225       225 s1         225 s1   225     High    23      A    emotion
#> 226       226 s1         226 s1   226     High    23      A  consensus
#> 227       227 s1         227 s1   227     High    23      A       plan
#> 228       228 s1         228 s1   228     High    23      A  consensus
#> 229       229 s1         229 s1   229     High    23      A   cohesion
#> 230       230 s1         230 s1   230     High    23      A  consensus
#> 231       231 s1         231 s1   231     High    24      A  consensus
#> 232       232 s1         232 s1   232     High    24      A       plan
#> 233       233 s1         233 s1   233     High    24      A    discuss
#> 234       234 s1         234 s1   234     High    24      A  consensus
#> 235       235 s1         235 s1   235     High    24      A    discuss
#> 236       236 s1         236 s1   236     High    24      A    discuss
#> 237       237 s1         237 s1   237     High    24      A    discuss
#> 238       238 s1         238 s1   238     High    24      A  consensus
#> 239       239 s1         239 s1   239     High    24      A       plan
#> 240       240 s1         240 s1   240     High    24      A  consensus
#> 241       241 s1         241 s1   241     High    25      A  consensus
#> 242       242 s1         242 s1   242     High    25      A  consensus
#> 243       243 s1         243 s1   243     High    25      A    monitor
#> 244       244 s1         244 s1   244     High    25      A  synthesis
#> 245       245 s1         245 s1   245     High    25      A  consensus
#> 246       246 s1         246 s1   246     High    25      A       plan
#> 247       247 s1         247 s1   247     High    25      A       plan
#> 248       248 s1         248 s1   248     High    25      A    emotion
#> 249       249 s1         249 s1   249     High    25      A    emotion
#> 250       250 s1         250 s1   250     High    25      A       plan
#> 251       251 s1         251 s1   251     High    26      A    monitor
#> 252       252 s1         252 s1   252     High    26      A   cohesion
#> 253       253 s1         253 s1   253     High    26      A    discuss
#> 254       254 s1         254 s1   254     High    26      A    emotion
#> 255       255 s1         255 s1   255     High    26      A  synthesis
#> 256       256 s1         256 s1   256     High    26      A    emotion
#> 257       257 s1         257 s1   257     High    26      A  consensus
#> 258       258 s1         258 s1   258     High    26      A   cohesion
#> 259       259 s1         259 s1   259     High    26      A    monitor
#> 260       260 s1         260 s1   260     High    26      A    emotion
#> 261       261 s1         261 s1   261     High    27      A  consensus
#> 262       262 s1         262 s1   262     High    27      A    emotion
#> 263       263 s1         263 s1   263     High    27      A       plan
#> 264       264 s1         264 s1   264     High    27      A  consensus
#> 265       265 s1         265 s1   265     High    27      A  consensus
#> 266       266 s1         266 s1   266     High    27      A    discuss
#> 267       267 s1         267 s1   267     High    27      A   cohesion
#> 268       268 s1         268 s1   268     High    27      A    discuss
#> 269       269 s1         269 s1   269     High    27      A    discuss
#> 270       270 s1         270 s1   270     High    27      A    emotion
#> 271       271 s1         271 s1   271     High    28      A  consensus
#> 272       272 s1         272 s1   272     High    28      A    discuss
#> 273       273 s1         273 s1   273     High    28      A   cohesion
#> 274       274 s1         274 s1   274     High    28      A       plan
#> 275       275 s1         275 s1   275     High    28      A coregulate
#> 276       276 s1         276 s1   276     High    28      A    discuss
#> 277       277 s1         277 s1   277     High    28      A    discuss
#> 278       278 s1         278 s1   278     High    28      A  consensus
#> 279       279 s1         279 s1   279     High    28      A  consensus
#> 280       280 s1         280 s1   280     High    28      A    emotion
#> 281       281 s1         281 s1   281     High    29      A       plan
#> 282       282 s1         282 s1   282     High    29      A    monitor
#> 283       283 s1         283 s1   283     High    29      A    emotion
#> 284       284 s1         284 s1   284     High    29      A       plan
#> 285       285 s1         285 s1   285     High    29      A       plan
#> 286       286 s1         286 s1   286     High    29      A    discuss
#> 287       287 s1         287 s1   287     High    29      A    discuss
#> 288       288 s1         288 s1   288     High    29      A  consensus
#> 289       289 s1         289 s1   289     High    29      A    monitor
#> 290       290 s1         290 s1   290     High    29      A    monitor
#> 291       291 s1         291 s1   291     High    30      A    monitor
#> 292       292 s1         292 s1   292     High    30      A  consensus
#> 293       293 s1         293 s1   293     High    30      A  consensus
#> 294       294 s1         294 s1   294     High    30      A    discuss
#> 295       295 s1         295 s1   295     High    30      A  consensus
#> 296       296 s1         296 s1   296     High    30      A    emotion
#> 297       297 s1         297 s1   297     High    30      A    monitor
#> 298       298 s1         298 s1   298     High    30      A    monitor
#> 299       299 s1         299 s1   299     High    30      A    emotion
#> 300       300 s1         300 s1   300     High    30      A    emotion
#> 301       301 s1         301 s1   301     High    31      A    discuss
#> 302       302 s1         302 s1   302     High    31      A    emotion
#> 303       303 s1         303 s1   303     High    31      A    emotion
#> 304       304 s1         304 s1   304     High    31      A    discuss
#> 305       305 s1         305 s1   305     High    31      A  consensus
#> 306       306 s1         306 s1   306     High    31      A    discuss
#> 307       307 s1         307 s1   307     High    31      A    discuss
#> 308       308 s1         308 s1   308     High    31      A   cohesion
#> 309       309 s1         309 s1   309     High    31      A    emotion
#> 310       310 s1         310 s1   310     High    31      A       plan
#> 311       311 s1         311 s1   311     High    32      A  consensus
#> 312       312 s1         312 s1   312     High    32      A    emotion
#> 313       313 s1         313 s1   313     High    32      A  consensus
#> 314       314 s1         314 s1   314     High    32      A    emotion
#> 315       315 s1         315 s1   315     High    32      A  consensus
#> 316       316 s1         316 s1   316     High    32      A    emotion
#> 317       317 s1         317 s1   317     High    32      A  consensus
#> 318       318 s1         318 s1   318     High    32      A  synthesis
#> 319       319 s1         319 s1   319     High    32      A       plan
#> 320       320 s1         320 s1   320     High    32      A   cohesion
#> 321       321 s1         321 s1   321     High    33      A    emotion
#> 322       322 s1         322 s1   322     High    33      A       plan
#> 323       323 s1         323 s1   323     High    33      A  consensus
#> 324       324 s1         324 s1   324     High    33      A    emotion
#> 325       325 s1         325 s1   325     High    33      A       plan
#> 326       326 s1         326 s1   326     High    33      A    monitor
#> 327       327 s1         327 s1   327     High    33      A      adapt
#> 328       328 s1         328 s1   328     High    33      A    monitor
#> 329       329 s1         329 s1   329     High    33      A  consensus
#> 330       330 s1         330 s1   330     High    33      A    emotion
#> 331       331 s1         331 s1   331     High    34      A    discuss
#> 332       332 s1         332 s1   332     High    34      A       plan
#> 333       333 s1         333 s1   333     High    34      A    monitor
#> 334       334 s1         334 s1   334     High    34      A  consensus
#> 335       335 s1         335 s1   335     High    34      A    emotion
#> 336       336 s1         336 s1   336     High    34      A    emotion
#> 337       337 s1         337 s1   337     High    34      A  consensus
#> 338       338 s1         338 s1   338     High    34      A    discuss
#> 339       339 s1         339 s1   339     High    34      A       plan
#> 340       340 s1         340 s1   340     High    34      A  consensus
#> 341       341 s1         341 s1   341     High    35      A    discuss
#> 342       342 s1         342 s1   342     High    35      A    emotion
#> 343       343 s1         343 s1   343     High    35      A  synthesis
#> 344       344 s1         344 s1   344     High    35      A    discuss
#> 345       345 s1         345 s1   345     High    35      A  consensus
#> 346       346 s1         346 s1   346     High    35      A  synthesis
#> 347       347 s1         347 s1   347     High    35      A    monitor
#> 348       348 s1         348 s1   348     High    35      A    emotion
#> 349       349 s1         349 s1   349     High    35      A       plan
#> 350       350 s1         350 s1   350     High    35      A       plan
#> 351       351 s1         351 s1   351     High    36      A coregulate
#> 352       352 s1         352 s1   352     High    36      A  consensus
#> 353       353 s1         353 s1   353     High    36      A    discuss
#> 354       354 s1         354 s1   354     High    36      A    discuss
#> 355       355 s1         355 s1   355     High    36      A    discuss
#> 356       356 s1         356 s1   356     High    36      A    emotion
#> 357       357 s1         357 s1   357     High    36      A  consensus
#> 358       358 s1         358 s1   358     High    36      A  consensus
#> 359       359 s1         359 s1   359     High    36      A    discuss
#> 360       360 s1         360 s1   360     High    36      A    monitor
#> 361       361 s1         361 s1   361     High    37      A  consensus
#> 362       362 s1         362 s1   362     High    37      A  synthesis
#> 363       363 s1         363 s1   363     High    37      A  consensus
#> 364       364 s1         364 s1   364     High    37      A       plan
#> 365       365 s1         365 s1   365     High    37      A    monitor
#> 366       366 s1         366 s1   366     High    37      A  consensus
#> 367       367 s1         367 s1   367     High    37      A       plan
#> 368       368 s1         368 s1   368     High    37      A   cohesion
#> 369       369 s1         369 s1   369     High    37      A    emotion
#> 370       370 s1         370 s1   370     High    37      A    emotion
#> 371       371 s1         371 s1   371     High    38      A       plan
#> 372       372 s1         372 s1   372     High    38      A  synthesis
#> 373       373 s1         373 s1   373     High    38      A       plan
#> 374       374 s1         374 s1   374     High    38      A    emotion
#> 375       375 s1         375 s1   375     High    38      A    emotion
#> 376       376 s1         376 s1   376     High    38      A  synthesis
#> 377       377 s1         377 s1   377     High    38      A      adapt
#> 378       378 s1         378 s1   378     High    38      A  consensus
#> 379       379 s1         379 s1   379     High    38      A  consensus
#> 380       380 s1         380 s1   380     High    38      A       plan
#> 381       381 s1         381 s1   381     High    39      A    discuss
#> 382       382 s1         382 s1   382     High    39      A    discuss
#> 383       383 s1         383 s1   383     High    39      A  consensus
#> 384       384 s1         384 s1   384     High    39      A    emotion
#> 385       385 s1         385 s1   385     High    39      A    discuss
#> 386       386 s1         386 s1   386     High    39      A    monitor
#> 387       387 s1         387 s1   387     High    39      A  consensus
#> 388       388 s1         388 s1   388     High    39      A    monitor
#> 389       389 s1         389 s1   389     High    39      A    emotion
#> 390       390 s1         390 s1   390     High    39      A    discuss
#> 391       391 s1         391 s1   391     High    40      A   cohesion
#> 392       392 s1         392 s1   392     High    40      A    monitor
#> 393       393 s1         393 s1   393     High    40      A  consensus
#> 394       394 s1         394 s1   394     High    40      A   cohesion
#> 395       395 s1         395 s1   395     High    40      A  consensus
#> 396       396 s1         396 s1   396     High    40      A    monitor
#> 397       397 s1         397 s1   397     High    40      A    discuss
#> 398       398 s1         398 s1   398     High    40      A    emotion
#> 399       399 s1         399 s1   399     High    40      A  consensus
#> 400       400 s1         400 s1   400     High    40      A    monitor
#> 401       401 s1         401 s1   401     High    41      A   cohesion
#> 402       402 s1         402 s1   402     High    41      A    discuss
#> 403       403 s1         403 s1   403     High    41      A       plan
#> 404       404 s1         404 s1   404     High    41      A    monitor
#> 405       405 s1         405 s1   405     High    41      A  consensus
#> 406       406 s1         406 s1   406     High    41      A    monitor
#> 407       407 s1         407 s1   407     High    41      A    discuss
#> 408       408 s1         408 s1   408     High    41      A    discuss
#> 409       409 s1         409 s1   409     High    41      A    discuss
#> 410       410 s1         410 s1   410     High    41      A       plan
#> 411       411 s1         411 s1   411     High    42      A       plan
#> 412       412 s1         412 s1   412     High    42      A    discuss
#> 413       413 s1         413 s1   413     High    42      A       plan
#> 414       414 s1         414 s1   414     High    42      A    discuss
#> 415       415 s1         415 s1   415     High    42      A    monitor
#> 416       416 s1         416 s1   416     High    42      A       plan
#> 417       417 s1         417 s1   417     High    42      A  consensus
#> 418       418 s1         418 s1   418     High    42      A    monitor
#> 419       419 s1         419 s1   419     High    42      A    monitor
#> 420       420 s1         420 s1   420     High    42      A   cohesion
#> 421       421 s1         421 s1   421     High    43      A       plan
#> 422       422 s1         422 s1   422     High    43      A       plan
#> 423       423 s1         423 s1   423     High    43      A    monitor
#> 424       424 s1         424 s1   424     High    43      A    discuss
#> 425       425 s1         425 s1   425     High    43      A    monitor
#> 426       426 s1         426 s1   426     High    43      A    emotion
#> 427       427 s1         427 s1   427     High    43      A    monitor
#> 428       428 s1         428 s1   428     High    43      A       plan
#> 429       429 s1         429 s1   429     High    43      A       plan
#> 430       430 s1         430 s1   430     High    43      A    monitor
#> 431       431 s1         431 s1   431     High    44      A  consensus
#> 432       432 s1         432 s1   432     High    44      A       plan
#> 433       433 s1         433 s1   433     High    44      A       plan
#> 434       434 s1         434 s1   434     High    44      A    emotion
#> 435       435 s1         435 s1   435     High    44      A       plan
#> 436       436 s1         436 s1   436     High    44      A    emotion
#> 437       437 s1         437 s1   437     High    44      A  synthesis
#> 438       438 s1         438 s1   438     High    44      A    discuss
#> 439       439 s1         439 s1   439     High    44      A    discuss
#> 440       440 s1         440 s1   440     High    44      A  consensus
#> 441       441 s1         441 s1   441     High    45      A    monitor
#> 442       442 s1         442 s1   442     High    45      A       plan
#> 443       443 s1         443 s1   443     High    45      A    discuss
#> 444       444 s1         444 s1   444     High    45      A    discuss
#> 445       445 s1         445 s1   445     High    45      A    emotion
#> 446       446 s1         446 s1   446     High    45      A   cohesion
#> 447       447 s1         447 s1   447     High    45      A  consensus
#> 448       448 s1         448 s1   448     High    45      A    emotion
#> 449       449 s1         449 s1   449     High    45      A  consensus
#> 450       450 s1         450 s1   450     High    45      A       plan
#> 451       451 s1         451 s1   451     High    46      A    monitor
#> 452       452 s1         452 s1   452     High    46      A    discuss
#> 453       453 s1         453 s1   453     High    46      A  consensus
#> 454       454 s1         454 s1   454     High    46      A       plan
#> 455       455 s1         455 s1   455     High    46      A  consensus
#> 456       456 s1         456 s1   456     High    46      A   cohesion
#> 457       457 s1         457 s1   457     High    46      A    emotion
#> 458       458 s1         458 s1   458     High    46      A  consensus
#> 459       459 s1         459 s1   459     High    46      A    discuss
#> 460       460 s1         460 s1   460     High    46      A    emotion
#> 461       461 s1         461 s1   461     High    47      A    emotion
#> 462       462 s1         462 s1   462     High    47      A       plan
#> 463       463 s1         463 s1   463     High    47      A  consensus
#> 464       464 s1         464 s1   464     High    47      A    emotion
#> 465       465 s1         465 s1   465     High    47      A       plan
#> 466       466 s1         466 s1   466     High    47      A    emotion
#> 467       467 s1         467 s1   467     High    47      A  consensus
#> 468       468 s1         468 s1   468     High    47      A  consensus
#> 469       469 s1         469 s1   469     High    47      A  consensus
#> 470       470 s1         470 s1   470     High    47      A       plan
#> 471       471 s1         471 s1   471     High    48      A    emotion
#> 472       472 s1         472 s1   472     High    48      A    discuss
#> 473       473 s1         473 s1   473     High    48      A    monitor
#> 474       474 s1         474 s1   474     High    48      A    emotion
#> 475       475 s1         475 s1   475     High    48      A coregulate
#> 476       476 s1         476 s1   476     High    48      A  synthesis
#> 477       477 s1         477 s1   477     High    48      A       plan
#> 478       478 s1         478 s1   478     High    48      A  consensus
#> 479       479 s1         479 s1   479     High    48      A    discuss
#> 480       480 s1         480 s1   480     High    48      A  synthesis
#> 481       481 s1         481 s1   481     High    49      A    discuss
#> 482       482 s1         482 s1   482     High    49      A       plan
#> 483       483 s1         483 s1   483     High    49      A    discuss
#> 484       484 s1         484 s1   484     High    49      A    discuss
#> 485       485 s1         485 s1   485     High    49      A  consensus
#> 486       486 s1         486 s1   486     High    49      A  consensus
#> 487       487 s1         487 s1   487     High    49      A    emotion
#> 488       488 s1         488 s1   488     High    49      A    monitor
#> 489       489 s1         489 s1   489     High    49      A    discuss
#> 490       490 s1         490 s1   490     High    49      A    monitor
#> 491       491 s1         491 s1   491     High    50      A       plan
#> 492       492 s1         492 s1   492     High    50      A  consensus
#> 493       493 s1         493 s1   493     High    50      A    emotion
#> 494       494 s1         494 s1   494     High    50      A    monitor
#> 495       495 s1         495 s1   495     High    50      A    monitor
#> 496       496 s1         496 s1   496     High    50      A   cohesion
#> 497       497 s1         497 s1   497     High    50      A  consensus
#> 498       498 s1         498 s1   498     High    50      A  consensus
#> 499       499 s1         499 s1   499     High    50      A    discuss
#> 500       500 s1         500 s1   500     High    50      A   cohesion
#> 501       501 s1         501 s1   501     High    51      A    emotion
#> 502       502 s1         502 s1   502     High    51      A       plan
#> 503       503 s1         503 s1   503     High    51      A  consensus
#> 504       504 s1         504 s1   504     High    51      A  consensus
#> 505       505 s1         505 s1   505     High    51      A  consensus
#> 506       506 s1         506 s1   506     High    51      A  consensus
#> 507       507 s1         507 s1   507     High    51      A    discuss
#> 508       508 s1         508 s1   508     High    51      A    discuss
#> 509       509 s1         509 s1   509     High    51      A    discuss
#> 510       510 s1         510 s1   510     High    51      A       plan
#> 511       511 s1         511 s1   511     High    52      A    emotion
#> 512       512 s1         512 s1   512     High    52      A  consensus
#> 513       513 s1         513 s1   513     High    52      A    emotion
#> 514       514 s1         514 s1   514     High    52      A    monitor
#> 515       515 s1         515 s1   515     High    52      A    discuss
#> 516       516 s1         516 s1   516     High    52      A    discuss
#> 517       517 s1         517 s1   517     High    52      A       plan
#> 518       518 s1         518 s1   518     High    52      A    emotion
#> 519       519 s1         519 s1   519     High    52      A    monitor
#> 520       520 s1         520 s1   520     High    52      A    discuss
#> 521       521 s1         521 s1   521     High    53      A    emotion
#> 522       522 s1         522 s1   522     High    53      A  consensus
#> 523       523 s1         523 s1   523     High    53      A  consensus
#> 524       524 s1         524 s1   524     High    53      A    emotion
#> 525       525 s1         525 s1   525     High    53      A       plan
#> 526       526 s1         526 s1   526     High    53      A       plan
#> 527       527 s1         527 s1   527     High    53      A    discuss
#> 528       528 s1         528 s1   528     High    53      A       plan
#> 529       529 s1         529 s1   529     High    53      A    discuss
#> 530       530 s1         530 s1   530     High    53      A    monitor
#> 531       531 s1         531 s1   531     High    54      A       plan
#> 532       532 s1         532 s1   532     High    54      A    emotion
#> 533       533 s1         533 s1   533     High    54      A  consensus
#> 534       534 s1         534 s1   534     High    54      A       plan
#> 535       535 s1         535 s1   535     High    54      A    emotion
#> 536       536 s1         536 s1   536     High    54      A    emotion
#> 537       537 s1         537 s1   537     High    54      A    emotion
#> 538       538 s1         538 s1   538     High    54      A    discuss
#> 539       539 s1         539 s1   539     High    54      A       plan
#> 540       540 s1         540 s1   540     High    54      A    discuss
#> 541       541 s1         541 s1   541     High    55      A       plan
#> 542       542 s1         542 s1   542     High    55      A    discuss
#> 543       543 s1         543 s1   543     High    55      A    monitor
#> 544       544 s1         544 s1   544     High    55      A    monitor
#> 545       545 s1         545 s1   545     High    55      A      adapt
#> 546       546 s1         546 s1   546     High    55      A    emotion
#> 547       547 s1         547 s1   547     High    55      A  consensus
#> 548       548 s1         548 s1   548     High    55      A  consensus
#> 549       549 s1         549 s1   549     High    55      A   cohesion
#> 550       550 s1         550 s1   550     High    55      A       plan
#> 551       551 s1         551 s1   551     High    56      A       plan
#> 552       552 s1         552 s1   552     High    56      A    discuss
#> 553       553 s1         553 s1   553     High    56      A       plan
#> 554       554 s1         554 s1   554     High    56      A       plan
#> 555       555 s1         555 s1   555     High    56      A  consensus
#> 556       556 s1         556 s1   556     High    56      A    discuss
#> 557       557 s1         557 s1   557     High    56      A    emotion
#> 558       558 s1         558 s1   558     High    56      A   cohesion
#> 559       559 s1         559 s1   559     High    56      A  consensus
#> 560       560 s1         560 s1   560     High    56      A    emotion
#> 561       561 s1         561 s1   561     High    57      A  consensus
#> 562       562 s1         562 s1   562     High    57      A    emotion
#> 563       563 s1         563 s1   563     High    57      A    emotion
#> 564       564 s1         564 s1   564     High    57      A    discuss
#> 565       565 s1         565 s1   565     High    57      A    discuss
#> 566       566 s1         566 s1   566     High    57      A    monitor
#> 567       567 s1         567 s1   567     High    57      A    emotion
#> 568       568 s1         568 s1   568     High    57      A      adapt
#> 569       569 s1         569 s1   569     High    57      A  consensus
#> 570       570 s1         570 s1   570     High    57      A       plan
#> 571       571 s1         571 s1   571     High    58      A  synthesis
#> 572       572 s1         572 s1   572     High    58      A   cohesion
#> 573       573 s1         573 s1   573     High    58      A    monitor
#> 574       574 s1         574 s1   574     High    58      A    emotion
#> 575       575 s1         575 s1   575     High    58      A    monitor
#> 576       576 s1         576 s1   576     High    58      A  consensus
#> 577       577 s1         577 s1   577     High    58      A       plan
#> 578       578 s1         578 s1   578     High    58      A    emotion
#> 579       579 s1         579 s1   579     High    58      A  consensus
#> 580       580 s1         580 s1   580     High    58      A    emotion
#> 581       581 s1         581 s1   581     High    59      A  consensus
#> 582       582 s1         582 s1   582     High    59      A    emotion
#> 583       583 s1         583 s1   583     High    59      A    discuss
#> 584       584 s1         584 s1   584     High    59      A   cohesion
#> 585       585 s1         585 s1   585     High    59      A  consensus
#> 586       586 s1         586 s1   586     High    59      A    discuss
#> 587       587 s1         587 s1   587     High    59      A    discuss
#> 588       588 s1         588 s1   588     High    59      A       plan
#> 589       589 s1         589 s1   589     High    59      A  consensus
#> 590       590 s1         590 s1   590     High    59      A  synthesis
#> 591       591 s1         591 s1   591     High    60      A    monitor
#> 592       592 s1         592 s1   592     High    60      A    monitor
#> 593       593 s1         593 s1   593     High    60      A  consensus
#> 594       594 s1         594 s1   594     High    60      A    emotion
#> 595       595 s1         595 s1   595     High    60      A  consensus
#> 596       596 s1         596 s1   596     High    60      A    discuss
#> 597       597 s1         597 s1   597     High    60      A    monitor
#> 598       598 s1         598 s1   598     High    60      A    monitor
#> 599       599 s1         599 s1   599     High    60      A  consensus
#> 600       600 s1         600 s1   600     High    60      A   cohesion
#> 601       601 s1         601 s1   601     High    61      A  consensus
#> 602       602 s1         602 s1   602     High    61      A   cohesion
#> 603       603 s1         603 s1   603     High    61      A       plan
#> 604       604 s1         604 s1   604     High    61      A    monitor
#> 605       605 s1         605 s1   605     High    61      A    discuss
#> 606       606 s1         606 s1   606     High    61      A    discuss
#> 607       607 s1         607 s1   607     High    61      A  consensus
#> 608       608 s1         608 s1   608     High    61      A    monitor
#> 609       609 s1         609 s1   609     High    61      A    emotion
#> 610       610 s1         610 s1   610     High    61      A    monitor
#> 611       611 s1         611 s1   611     High    62      A    emotion
#> 612       612 s1         612 s1   612     High    62      A    emotion
#> 613       613 s1         613 s1   613     High    62      A  consensus
#> 614       614 s1         614 s1   614     High    62      A    emotion
#> 615       615 s1         615 s1   615     High    62      A  consensus
#> 616       616 s1         616 s1   616     High    62      A    discuss
#> 617       617 s1         617 s1   617     High    62      A  consensus
#> 618       618 s1         618 s1   618     High    62      A    discuss
#> 619       619 s1         619 s1   619     High    62      A    discuss
#> 620       620 s1         620 s1   620     High    62      A    monitor
#> 621       621 s1         621 s1   621     High    63      A  consensus
#> 622       622 s1         622 s1   622     High    63      A    emotion
#> 623       623 s1         623 s1   623     High    63      A    discuss
#> 624       624 s1         624 s1   624     High    63      A  consensus
#> 625       625 s1         625 s1   625     High    63      A    emotion
#> 626       626 s1         626 s1   626     High    63      A   cohesion
#> 627       627 s1         627 s1   627     High    63      A  consensus
#> 628       628 s1         628 s1   628     High    63      A  consensus
#> 629       629 s1         629 s1   629     High    63      A    monitor
#> 630       630 s1         630 s1   630     High    63      A       plan
#> 631       631 s1         631 s1   631     High    64      A    discuss
#> 632       632 s1         632 s1   632     High    64      A    discuss
#> 633       633 s1         633 s1   633     High    64      A    emotion
#> 634       634 s1         634 s1   634     High    64      A       plan
#> 635       635 s1         635 s1   635     High    64      A    emotion
#> 636       636 s1         636 s1   636     High    64      A    monitor
#> 637       637 s1         637 s1   637     High    64      A       plan
#> 638       638 s1         638 s1   638     High    64      A    emotion
#> 639       639 s1         639 s1   639     High    64      A  consensus
#> 640       640 s1         640 s1   640     High    64      A  synthesis
#> 641       641 s1         641 s1   641     High    65      A       plan
#> 642       642 s1         642 s1   642     High    65      A   cohesion
#> 643       643 s1         643 s1   643     High    65      A    monitor
#> 644       644 s1         644 s1   644     High    65      A       plan
#> 645       645 s1         645 s1   645     High    65      A      adapt
#> 646       646 s1         646 s1   646     High    65      A   cohesion
#> 647       647 s1         647 s1   647     High    65      A   cohesion
#> 648       648 s1         648 s1   648     High    65      A  synthesis
#> 649       649 s1         649 s1   649     High    65      A       plan
#> 650       650 s1         650 s1   650     High    65      A    monitor
#> 651       651 s1         651 s1   651     High    66      A    emotion
#> 652       652 s1         652 s1   652     High    66      A  consensus
#> 653       653 s1         653 s1   653     High    66      A       plan
#> 654       654 s1         654 s1   654     High    66      A  consensus
#> 655       655 s1         655 s1   655     High    66      A       plan
#> 656       656 s1         656 s1   656     High    66      A       plan
#> 657       657 s1         657 s1   657     High    66      A       plan
#> 658       658 s1         658 s1   658     High    66      A    discuss
#> 659       659 s1         659 s1   659     High    66      A    monitor
#> 660       660 s1         660 s1   660     High    66      A    monitor
#> 661       661 s1         661 s1   661     High    67      A      adapt
#> 662       662 s1         662 s1   662     High    67      A       plan
#> 663       663 s1         663 s1   663     High    67      A  consensus
#> 664       664 s1         664 s1   664     High    67      A      adapt
#> 665       665 s1         665 s1   665     High    67      A   cohesion
#> 666       666 s1         666 s1   666     High    67      A    discuss
#> 667       667 s1         667 s1   667     High    67      A       plan
#> 668       668 s1         668 s1   668     High    67      A  consensus
#> 669       669 s1         669 s1   669     High    67      A    emotion
#> 670       670 s1         670 s1   670     High    67      A  consensus
#> 671       671 s1         671 s1   671     High    68      A   cohesion
#> 672       672 s1         672 s1   672     High    68      A    emotion
#> 673       673 s1         673 s1   673     High    68      A       plan
#> 674       674 s1         674 s1   674     High    68      A    discuss
#> 675       675 s1         675 s1   675     High    68      A   cohesion
#> 676       676 s1         676 s1   676     High    68      A    emotion
#> 677       677 s1         677 s1   677     High    68      A   cohesion
#> 678       678 s1         678 s1   678     High    68      A    discuss
#> 679       679 s1         679 s1   679     High    68      A   cohesion
#> 680       680 s1         680 s1   680     High    68      A    discuss
#> 681       681 s1         681 s1   681     High    69      A    emotion
#> 682       682 s1         682 s1   682     High    69      A    monitor
#> 683       683 s1         683 s1   683     High    69      A    discuss
#> 684       684 s1         684 s1   684     High    69      A       plan
#> 685       685 s1         685 s1   685     High    69      A       plan
#> 686       686 s1         686 s1   686     High    69      A  consensus
#> 687       687 s1         687 s1   687     High    69      A    monitor
#> 688       688 s1         688 s1   688     High    69      A       plan
#> 689       689 s1         689 s1   689     High    69      A  consensus
#> 690       690 s1         690 s1   690     High    69      A    emotion
#> 691       691 s1         691 s1   691     High    70      A    emotion
#> 692       692 s1         692 s1   692     High    70      A    emotion
#> 693       693 s1         693 s1   693     High    70      A  consensus
#> 694       694 s1         694 s1   694     High    70      A  consensus
#> 695       695 s1         695 s1   695     High    70      A    discuss
#> 696       696 s1         696 s1   696     High    70      A    discuss
#> 697       697 s1         697 s1   697     High    70      A    monitor
#> 698       698 s1         698 s1   698     High    70      A    monitor
#> 699       699 s1         699 s1   699     High    70      A    discuss
#> 700       700 s1         700 s1   700     High    70      A  consensus
#> 701       701 s1         701 s1   701     High    71      A    emotion
#> 702       702 s1         702 s1   702     High    71      A    emotion
#> 703       703 s1         703 s1   703     High    71      A  consensus
#> 704       704 s1         704 s1   704     High    71      A    discuss
#> 705       705 s1         705 s1   705     High    71      A    discuss
#> 706       706 s1         706 s1   706     High    71      A  consensus
#> 707       707 s1         707 s1   707     High    71      A    discuss
#> 708       708 s1         708 s1   708     High    71      A  consensus
#> 709       709 s1         709 s1   709     High    71      A    monitor
#> 710       710 s1         710 s1   710     High    71      A    discuss
#> 711       711 s1         711 s1   711     High    72      A       plan
#> 712       712 s1         712 s1   712     High    72      A   cohesion
#> 713       713 s1         713 s1   713     High    72      A    emotion
#> 714       714 s1         714 s1   714     High    72      A    monitor
#> 715       715 s1         715 s1   715     High    72      A       plan
#> 716       716 s1         716 s1   716     High    72      A    monitor
#> 717       717 s1         717 s1   717     High    72      A       plan
#> 718       718 s1         718 s1   718     High    72      A       plan
#> 719       719 s1         719 s1   719     High    72      A       plan
#> 720       720 s1         720 s1   720     High    72      A  consensus
#> 721       721 s1         721 s1   721     High    73      A coregulate
#> 722       722 s1         722 s1   722     High    73      A    emotion
#> 723       723 s1         723 s1   723     High    73      A    discuss
#> 724       724 s1         724 s1   724     High    73      A    discuss
#> 725       725 s1         725 s1   725     High    73      A    discuss
#> 726       726 s1         726 s1   726     High    73      A       plan
#> 727       727 s1         727 s1   727     High    73      A    emotion
#> 728       728 s1         728 s1   728     High    73      A   cohesion
#> 729       729 s1         729 s1   729     High    73      A   cohesion
#> 730       730 s1         730 s1   730     High    73      A    discuss
#> 731       731 s1         731 s1   731     High    74      A  consensus
#> 732       732 s1         732 s1   732     High    74      A    discuss
#> 733       733 s1         733 s1   733     High    74      A    emotion
#> 734       734 s1         734 s1   734     High    74      A  consensus
#> 735       735 s1         735 s1   735     High    74      A       plan
#> 736       736 s1         736 s1   736     High    74      A    discuss
#> 737       737 s1         737 s1   737     High    74      A  consensus
#> 738       738 s1         738 s1   738     High    74      A    monitor
#> 739       739 s1         739 s1   739     High    74      A    monitor
#> 740       740 s1         740 s1   740     High    74      A  consensus
#> 741       741 s1         741 s1   741     High    75      A coregulate
#> 742       742 s1         742 s1   742     High    75      A      adapt
#> 743       743 s1         743 s1   743     High    75      A   cohesion
#> 744       744 s1         744 s1   744     High    75      A    discuss
#> 745       745 s1         745 s1   745     High    75      A  consensus
#> 746       746 s1         746 s1   746     High    75      A    monitor
#> 747       747 s1         747 s1   747     High    75      A    emotion
#> 748       748 s1         748 s1   748     High    75      A    emotion
#> 749       749 s1         749 s1   749     High    75      A       plan
#> 750       750 s1         750 s1   750     High    75      A   cohesion
#> 751       751 s1         751 s1   751     High    76      A  consensus
#> 752       752 s1         752 s1   752     High    76      A    emotion
#> 753       753 s1         753 s1   753     High    76      A  consensus
#> 754       754 s1         754 s1   754     High    76      A    discuss
#> 755       755 s1         755 s1   755     High    76      A    discuss
#> 756       756 s1         756 s1   756     High    76      A    emotion
#> 757       757 s1         757 s1   757     High    76      A    emotion
#> 758       758 s1         758 s1   758     High    76      A  consensus
#> 759       759 s1         759 s1   759     High    76      A       plan
#> 760       760 s1         760 s1   760     High    76      A       plan
#> 761       761 s1         761 s1   761     High    77      A       plan
#> 762       762 s1         762 s1   762     High    77      A       plan
#> 763       763 s1         763 s1   763     High    77      A       plan
#> 764       764 s1         764 s1   764     High    77      A  consensus
#> 765       765 s1         765 s1   765     High    77      A    emotion
#> 766       766 s1         766 s1   766     High    77      A       plan
#> 767       767 s1         767 s1   767     High    77      A   cohesion
#> 768       768 s1         768 s1   768     High    77      A   cohesion
#> 769       769 s1         769 s1   769     High    77      A    emotion
#> 770       770 s1         770 s1   770     High    77      A    monitor
#> 771       771 s1         771 s1   771     High    78      A    discuss
#> 772       772 s1         772 s1   772     High    78      A  consensus
#> 773       773 s1         773 s1   773     High    78      A  synthesis
#> 774       774 s1         774 s1   774     High    78      A    emotion
#> 775       775 s1         775 s1   775     High    78      A   cohesion
#> 776       776 s1         776 s1   776     High    78      A   cohesion
#> 777       777 s1         777 s1   777     High    78      A      adapt
#> 778       778 s1         778 s1   778     High    78      A    discuss
#> 779       779 s1         779 s1   779     High    78      A    discuss
#> 780       780 s1         780 s1   780     High    78      A    emotion
#> 781       781 s1         781 s1   781     High    79      A   cohesion
#> 782       782 s1         782 s1   782     High    79      A   cohesion
#> 783       783 s1         783 s1   783     High    79      A    monitor
#> 784       784 s1         784 s1   784     High    79      A  consensus
#> 785       785 s1         785 s1   785     High    79      A    discuss
#> 786       786 s1         786 s1   786     High    79      A    monitor
#> 787       787 s1         787 s1   787     High    79      A  consensus
#> 788       788 s1         788 s1   788     High    79      A    discuss
#> 789       789 s1         789 s1   789     High    79      A    discuss
#> 790       790 s1         790 s1   790     High    79      A  consensus
#> 791       791 s1         791 s1   791     High    80      A    discuss
#> 792       792 s1         792 s1   792     High    80      A    emotion
#> 793       793 s1         793 s1   793     High    80      A   cohesion
#> 794       794 s1         794 s1   794     High    80      A   cohesion
#> 795       795 s1         795 s1   795     High    80      A    monitor
#> 796       796 s1         796 s1   796     High    80      A       plan
#> 797       797 s1         797 s1   797     High    80      A    discuss
#> 798       798 s1         798 s1   798     High    80      A    discuss
#> 799       799 s1         799 s1   799     High    80      A  consensus
#> 800       800 s1         800 s1   800     High    80      A    monitor
#> 801       801 s1         801 s1   801     High    81      A       plan
#> 802       802 s1         802 s1   802     High    81      A    emotion
#> 803       803 s1         803 s1   803     High    81      A    emotion
#> 804       804 s1         804 s1   804     High    81      A  consensus
#> 805       805 s1         805 s1   805     High    81      A  consensus
#> 806       806 s1         806 s1   806     High    81      A    emotion
#> 807       807 s1         807 s1   807     High    81      A       plan
#> 808       808 s1         808 s1   808     High    81      A    monitor
#> 809       809 s1         809 s1   809     High    81      A    discuss
#> 810       810 s1         810 s1   810     High    81      A    emotion
#> 811       811 s1         811 s1   811     High    82      A    monitor
#> 812       812 s1         812 s1   812     High    82      A    discuss
#> 813       813 s1         813 s1   813     High    82      A       plan
#> 814       814 s1         814 s1   814     High    82      A       plan
#> 815       815 s1         815 s1   815     High    82      A  consensus
#> 816       816 s1         816 s1   816     High    82      A    discuss
#> 817       817 s1         817 s1   817     High    82      A       plan
#> 818       818 s1         818 s1   818     High    82      A    discuss
#> 819       819 s1         819 s1   819     High    82      A       plan
#> 820       820 s1         820 s1   820     High    82      A  consensus
#> 821       821 s1         821 s1   821     High    83      A  consensus
#> 822       822 s1         822 s1   822     High    83      A       plan
#> 823       823 s1         823 s1   823     High    83      A   cohesion
#> 824       824 s1         824 s1   824     High    83      A    discuss
#> 825       825 s1         825 s1   825     High    83      A   cohesion
#> 826       826 s1         826 s1   826     High    83      A  consensus
#> 827       827 s1         827 s1   827     High    83      A    emotion
#> 828       828 s1         828 s1   828     High    83      A    monitor
#> 829       829 s1         829 s1   829     High    83      A   cohesion
#> 830       830 s1         830 s1   830     High    83      A    discuss
#> 831       831 s1         831 s1   831     High    84      A  consensus
#> 832       832 s1         832 s1   832     High    84      A  consensus
#> 833       833 s1         833 s1   833     High    84      A  consensus
#> 834       834 s1         834 s1   834     High    84      A    monitor
#> 835       835 s1         835 s1   835     High    84      A   cohesion
#> 836       836 s1         836 s1   836     High    84      A       plan
#> 837       837 s1         837 s1   837     High    84      A       plan
#> 838       838 s1         838 s1   838     High    84      A    emotion
#> 839       839 s1         839 s1   839     High    84      A    monitor
#> 840       840 s1         840 s1   840     High    84      A       plan
#> 841       841 s1         841 s1   841     High    85      A    monitor
#> 842       842 s1         842 s1   842     High    85      A       plan
#> 843       843 s1         843 s1   843     High    85      A    discuss
#> 844       844 s1         844 s1   844     High    85      A    emotion
#> 845       845 s1         845 s1   845     High    85      A  consensus
#> 846       846 s1         846 s1   846     High    85      A    discuss
#> 847       847 s1         847 s1   847     High    85      A  consensus
#> 848       848 s1         848 s1   848     High    85      A    discuss
#> 849       849 s1         849 s1   849     High    85      A   cohesion
#> 850       850 s1         850 s1   850     High    85      A   cohesion
#> 851       851 s1         851 s1   851     High    86      A    discuss
#> 852       852 s1         852 s1   852     High    86      A   cohesion
#> 853       853 s1         853 s1   853     High    86      A  consensus
#> 854       854 s1         854 s1   854     High    86      A    monitor
#> 855       855 s1         855 s1   855     High    86      A  consensus
#> 856       856 s1         856 s1   856     High    86      A    discuss
#> 857       857 s1         857 s1   857     High    86      A    monitor
#> 858       858 s1         858 s1   858     High    86      A       plan
#> 859       859 s1         859 s1   859     High    86      A  consensus
#> 860       860 s1         860 s1   860     High    86      A    discuss
#> 861       861 s1         861 s1   861     High    87      A    discuss
#> 862       862 s1         862 s1   862     High    87      A  consensus
#> 863       863 s1         863 s1   863     High    87      A       plan
#> 864       864 s1         864 s1   864     High    87      A    discuss
#> 865       865 s1         865 s1   865     High    87      A  consensus
#> 866       866 s1         866 s1   866     High    87      A  consensus
#> 867       867 s1         867 s1   867     High    87      A    discuss
#> 868       868 s1         868 s1   868     High    87      A  consensus
#> 869       869 s1         869 s1   869     High    87      A    monitor
#> 870       870 s1         870 s1   870     High    87      A  consensus
#> 871       871 s1         871 s1   871     High    88      A    emotion
#> 872       872 s1         872 s1   872     High    88      A  synthesis
#> 873       873 s1         873 s1   873     High    88      A  consensus
#> 874       874 s1         874 s1   874     High    88      A       plan
#> 875       875 s1         875 s1   875     High    88      A    discuss
#> 876       876 s1         876 s1   876     High    88      A  synthesis
#> 877       877 s1         877 s1   877     High    88      A       plan
#> 878       878 s1         878 s1   878     High    88      A       plan
#> 879       879 s1         879 s1   879     High    88      A       plan
#> 880       880 s1         880 s1   880     High    88      A       plan
#> 881       881 s1         881 s1   881     High    89      A       plan
#> 882       882 s1         882 s1   882     High    89      A    emotion
#> 883       883 s1         883 s1   883     High    89      A    discuss
#> 884       884 s1         884 s1   884     High    89      A  consensus
#> 885       885 s1         885 s1   885     High    89      A    monitor
#> 886       886 s1         886 s1   886     High    89      A       plan
#> 887       887 s1         887 s1   887     High    89      A    discuss
#> 888       888 s1         888 s1   888     High    89      A  consensus
#> 889       889 s1         889 s1   889     High    89      A    monitor
#> 890       890 s1         890 s1   890     High    89      A    monitor
#> 891       891 s1         891 s1   891     High    90      A  consensus
#> 892       892 s1         892 s1   892     High    90      A       plan
#> 893       893 s1         893 s1   893     High    90      A  consensus
#> 894       894 s1         894 s1   894     High    90      A  consensus
#> 895       895 s1         895 s1   895     High    90      A    discuss
#> 896       896 s1         896 s1   896     High    90      A  consensus
#> 897       897 s1         897 s1   897     High    90      A    emotion
#> 898       898 s1         898 s1   898     High    90      A    monitor
#> 899       899 s1         899 s1   899     High    90      A    discuss
#> 900       900 s1         900 s1   900     High    90      A  consensus
#> 901       901 s1         901 s1   901     High    91      B  consensus
#> 902       902 s1         902 s1   902     High    91      B       plan
#> 903       903 s1         903 s1   903     High    91      B       plan
#> 904       904 s1         904 s1   904     High    91      B    discuss
#> 905       905 s1         905 s1   905     High    91      B    emotion
#> 906       906 s1         906 s1   906     High    91      B       plan
#> 907       907 s1         907 s1   907     High    91      B    monitor
#> 908       908 s1         908 s1   908     High    91      B    emotion
#> 909       909 s1         909 s1   909     High    91      B       plan
#> 910       910 s1         910 s1   910     High    91      B    discuss
#> 911       911 s1         911 s1   911     High    92      B    discuss
#> 912       912 s1         912 s1   912     High    92      B    emotion
#> 913       913 s1         913 s1   913     High    92      B       plan
#> 914       914 s1         914 s1   914     High    92      B    emotion
#> 915       915 s1         915 s1   915     High    92      B    discuss
#> 916       916 s1         916 s1   916     High    92      B   cohesion
#> 917       917 s1         917 s1   917     High    92      B    emotion
#> 918       918 s1         918 s1   918     High    92      B  consensus
#> 919       919 s1         919 s1   919     High    92      B    discuss
#> 920       920 s1         920 s1   920     High    92      B  consensus
#> 921       921 s1         921 s1   921     High    93      B       plan
#> 922       922 s1         922 s1   922     High    93      B    discuss
#> 923       923 s1         923 s1   923     High    93      B    emotion
#> 924       924 s1         924 s1   924     High    93      B  consensus
#> 925       925 s1         925 s1   925     High    93      B       plan
#> 926       926 s1         926 s1   926     High    93      B    emotion
#> 927       927 s1         927 s1   927     High    93      B       plan
#> 928       928 s1         928 s1   928     High    93      B   cohesion
#> 929       929 s1         929 s1   929     High    93      B    discuss
#> 930       930 s1         930 s1   930     High    93      B    emotion
#> 931       931 s1         931 s1   931     High    94      B    emotion
#> 932       932 s1         932 s1   932     High    94      B  consensus
#> 933       933 s1         933 s1   933     High    94      B  consensus
#> 934       934 s1         934 s1   934     High    94      B    emotion
#> 935       935 s1         935 s1   935     High    94      B    emotion
#> 936       936 s1         936 s1   936     High    94      B  consensus
#> 937       937 s1         937 s1   937     High    94      B    discuss
#> 938       938 s1         938 s1   938     High    94      B    emotion
#> 939       939 s1         939 s1   939     High    94      B       plan
#> 940       940 s1         940 s1   940     High    94      B   cohesion
#> 941       941 s1         941 s1   941     High    95      B    emotion
#> 942       942 s1         942 s1   942     High    95      B   cohesion
#> 943       943 s1         943 s1   943     High    95      B  consensus
#> 944       944 s1         944 s1   944     High    95      B       plan
#> 945       945 s1         945 s1   945     High    95      B    monitor
#> 946       946 s1         946 s1   946     High    95      B    emotion
#> 947       947 s1         947 s1   947     High    95      B  consensus
#> 948       948 s1         948 s1   948     High    95      B    emotion
#> 949       949 s1         949 s1   949     High    95      B    monitor
#> 950       950 s1         950 s1   950     High    95      B       plan
#> 951       951 s1         951 s1   951     High    96      B    monitor
#> 952       952 s1         952 s1   952     High    96      B   cohesion
#> 953       953 s1         953 s1   953     High    96      B    emotion
#> 954       954 s1         954 s1   954     High    96      B  consensus
#> 955       955 s1         955 s1   955     High    96      B  consensus
#> 956       956 s1         956 s1   956     High    96      B    emotion
#> 957       957 s1         957 s1   957     High    96      B    discuss
#> 958       958 s1         958 s1   958     High    96      B       plan
#> 959       959 s1         959 s1   959     High    96      B    emotion
#> 960       960 s1         960 s1   960     High    96      B       plan
#> 961       961 s1         961 s1   961     High    97      B       plan
#> 962       962 s1         962 s1   962     High    97      B       plan
#> 963       963 s1         963 s1   963     High    97      B    emotion
#> 964       964 s1         964 s1   964     High    97      B   cohesion
#> 965       965 s1         965 s1   965     High    97      B    emotion
#> 966       966 s1         966 s1   966     High    97      B    emotion
#> 967       967 s1         967 s1   967     High    97      B       plan
#> 968       968 s1         968 s1   968     High    97      B   cohesion
#> 969       969 s1         969 s1   969     High    97      B  consensus
#> 970       970 s1         970 s1   970     High    97      B  consensus
#> 971       971 s1         971 s1   971     High    98      B  consensus
#> 972       972 s1         972 s1   972     High    98      B    emotion
#> 973       973 s1         973 s1   973     High    98      B       plan
#> 974       974 s1         974 s1   974     High    98      B    discuss
#> 975       975 s1         975 s1   975     High    98      B  consensus
#> 976       976 s1         976 s1   976     High    98      B    emotion
#> 977       977 s1         977 s1   977     High    98      B       plan
#> 978       978 s1         978 s1   978     High    98      B  consensus
#> 979       979 s1         979 s1   979     High    98      B       plan
#> 980       980 s1         980 s1   980     High    98      B    emotion
#> 981       981 s1         981 s1   981     High    99      B    monitor
#> 982       982 s1         982 s1   982     High    99      B    discuss
#> 983       983 s1         983 s1   983     High    99      B    discuss
#> 984       984 s1         984 s1   984     High    99      B       plan
#> 985       985 s1         985 s1   985     High    99      B  consensus
#> 986       986 s1         986 s1   986     High    99      B   cohesion
#> 987       987 s1         987 s1   987     High    99      B  consensus
#> 988       988 s1         988 s1   988     High    99      B    discuss
#> 989       989 s1         989 s1   989     High    99      B    monitor
#> 990       990 s1         990 s1   990     High    99      B       plan
#> 991       991 s1         991 s1   991     High   100      B  consensus
#> 992       992 s1         992 s1   992     High   100      B       plan
#> 993       993 s1         993 s1   993     High   100      B    emotion
#> 994       994 s1         994 s1   994     High   100      B    monitor
#> 995       995 s1         995 s1   995     High   100      B  consensus
#> 996       996 s1         996 s1   996     High   100      B      adapt
#> 997       997 s1         997 s1   997     High   100      B       plan
#> 998       998 s1         998 s1   998     High   100      B    monitor
#> 999       999 s1         999 s1   999     High   100      B    discuss
#> 1000     1000 s1        1000 s1  1000     High   100      B    discuss
#> 1001     1001 s1        1001 s1  1001      Low   101      B   cohesion
#> 1002     1002 s1        1002 s1  1002      Low   101      B    discuss
#> 1003     1003 s1        1003 s1  1003      Low   101      B    discuss
#> 1004     1004 s1        1004 s1  1004      Low   101      B   cohesion
#> 1005     1005 s1        1005 s1  1005      Low   101      B    emotion
#> 1006     1006 s1        1006 s1  1006      Low   101      B    emotion
#> 1007     1007 s1        1007 s1  1007      Low   101      B  consensus
#> 1008     1008 s1        1008 s1  1008      Low   101      B    discuss
#> 1009     1009 s1        1009 s1  1009      Low   101      B       plan
#> 1010     1010 s1        1010 s1  1010      Low   101      B  consensus
#> 1011     1011 s1        1011 s1  1011      Low   102      B    discuss
#> 1012     1012 s1        1012 s1  1012      Low   102      B  synthesis
#> 1013     1013 s1        1013 s1  1013      Low   102      B  consensus
#> 1014     1014 s1        1014 s1  1014      Low   102      B  consensus
#> 1015     1015 s1        1015 s1  1015      Low   102      B  consensus
#> 1016     1016 s1        1016 s1  1016      Low   102      B    emotion
#> 1017     1017 s1        1017 s1  1017      Low   102      B    monitor
#> 1018     1018 s1        1018 s1  1018      Low   102      B  consensus
#> 1019     1019 s1        1019 s1  1019      Low   102      B    discuss
#> 1020     1020 s1        1020 s1  1020      Low   102      B  consensus
#> 1021     1021 s1        1021 s1  1021      Low   103      B    emotion
#> 1022     1022 s1        1022 s1  1022      Low   103      B    discuss
#> 1023     1023 s1        1023 s1  1023      Low   103      B    monitor
#> 1024     1024 s1        1024 s1  1024      Low   103      B    discuss
#> 1025     1025 s1        1025 s1  1025      Low   103      B    discuss
#> 1026     1026 s1        1026 s1  1026      Low   103      B    discuss
#> 1027     1027 s1        1027 s1  1027      Low   103      B    discuss
#> 1028     1028 s1        1028 s1  1028      Low   103      B       plan
#> 1029     1029 s1        1029 s1  1029      Low   103      B       plan
#> 1030     1030 s1        1030 s1  1030      Low   103      B  consensus
#> 1031     1031 s1        1031 s1  1031      Low   104      B   cohesion
#> 1032     1032 s1        1032 s1  1032      Low   104      B    discuss
#> 1033     1033 s1        1033 s1  1033      Low   104      B    emotion
#> 1034     1034 s1        1034 s1  1034      Low   104      B       plan
#> 1035     1035 s1        1035 s1  1035      Low   104      B    discuss
#> 1036     1036 s1        1036 s1  1036      Low   104      B  consensus
#> 1037     1037 s1        1037 s1  1037      Low   104      B    discuss
#> 1038     1038 s1        1038 s1  1038      Low   104      B  consensus
#> 1039     1039 s1        1039 s1  1039      Low   104      B    emotion
#> 1040     1040 s1        1040 s1  1040      Low   104      B  consensus
#> 1041     1041 s1        1041 s1  1041      Low   105      B       plan
#> 1042     1042 s1        1042 s1  1042      Low   105      B       plan
#> 1043     1043 s1        1043 s1  1043      Low   105      B    discuss
#> 1044     1044 s1        1044 s1  1044      Low   105      B coregulate
#> 1045     1045 s1        1045 s1  1045      Low   105      B  consensus
#> 1046     1046 s1        1046 s1  1046      Low   105      B       plan
#> 1047     1047 s1        1047 s1  1047      Low   105      B       plan
#> 1048     1048 s1        1048 s1  1048      Low   105      B    discuss
#> 1049     1049 s1        1049 s1  1049      Low   105      B       plan
#> 1050     1050 s1        1050 s1  1050      Low   105      B    discuss
#> 1051     1051 s1        1051 s1  1051      Low   106      B  consensus
#> 1052     1052 s1        1052 s1  1052      Low   106      B  consensus
#> 1053     1053 s1        1053 s1  1053      Low   106      B    emotion
#> 1054     1054 s1        1054 s1  1054      Low   106      B    discuss
#> 1055     1055 s1        1055 s1  1055      Low   106      B    emotion
#> 1056     1056 s1        1056 s1  1056      Low   106      B  consensus
#> 1057     1057 s1        1057 s1  1057      Low   106      B       plan
#> 1058     1058 s1        1058 s1  1058      Low   106      B    emotion
#> 1059     1059 s1        1059 s1  1059      Low   106      B  consensus
#> 1060     1060 s1        1060 s1  1060      Low   106      B  consensus
#> 1061     1061 s1        1061 s1  1061      Low   107      B    monitor
#> 1062     1062 s1        1062 s1  1062      Low   107      B coregulate
#> 1063     1063 s1        1063 s1  1063      Low   107      B  consensus
#> 1064     1064 s1        1064 s1  1064      Low   107      B coregulate
#> 1065     1065 s1        1065 s1  1065      Low   107      B       plan
#> 1066     1066 s1        1066 s1  1066      Low   107      B  consensus
#> 1067     1067 s1        1067 s1  1067      Low   107      B  synthesis
#> 1068     1068 s1        1068 s1  1068      Low   107      B    monitor
#> 1069     1069 s1        1069 s1  1069      Low   107      B       plan
#> 1070     1070 s1        1070 s1  1070      Low   107      B       plan
#> 1071     1071 s1        1071 s1  1071      Low   108      B  consensus
#> 1072     1072 s1        1072 s1  1072      Low   108      B  consensus
#> 1073     1073 s1        1073 s1  1073      Low   108      B       plan
#> 1074     1074 s1        1074 s1  1074      Low   108      B       plan
#> 1075     1075 s1        1075 s1  1075      Low   108      B    monitor
#> 1076     1076 s1        1076 s1  1076      Low   108      B    monitor
#> 1077     1077 s1        1077 s1  1077      Low   108      B    monitor
#> 1078     1078 s1        1078 s1  1078      Low   108      B       plan
#> 1079     1079 s1        1079 s1  1079      Low   108      B    discuss
#> 1080     1080 s1        1080 s1  1080      Low   108      B  consensus
#> 1081     1081 s1        1081 s1  1081      Low   109      B    discuss
#> 1082     1082 s1        1082 s1  1082      Low   109      B    discuss
#> 1083     1083 s1        1083 s1  1083      Low   109      B coregulate
#> 1084     1084 s1        1084 s1  1084      Low   109      B  consensus
#> 1085     1085 s1        1085 s1  1085      Low   109      B  consensus
#> 1086     1086 s1        1086 s1  1086      Low   109      B    emotion
#> 1087     1087 s1        1087 s1  1087      Low   109      B      adapt
#> 1088     1088 s1        1088 s1  1088      Low   109      B    monitor
#> 1089     1089 s1        1089 s1  1089      Low   109      B  consensus
#> 1090     1090 s1        1090 s1  1090      Low   109      B   cohesion
#> 1091     1091 s1        1091 s1  1091      Low   110      B       plan
#> 1092     1092 s1        1092 s1  1092      Low   110      B  consensus
#> 1093     1093 s1        1093 s1  1093      Low   110      B       plan
#> 1094     1094 s1        1094 s1  1094      Low   110      B    emotion
#> 1095     1095 s1        1095 s1  1095      Low   110      B  consensus
#> 1096     1096 s1        1096 s1  1096      Low   110      B    emotion
#> 1097     1097 s1        1097 s1  1097      Low   110      B    emotion
#> 1098     1098 s1        1098 s1  1098      Low   110      B  consensus
#> 1099     1099 s1        1099 s1  1099      Low   110      B  consensus
#> 1100     1100 s1        1100 s1  1100      Low   110      B    emotion
#> 1101     1101 s1        1101 s1  1101      Low   111      B       plan
#> 1102     1102 s1        1102 s1  1102      Low   111      B       plan
#> 1103     1103 s1        1103 s1  1103      Low   111      B       plan
#> 1104     1104 s1        1104 s1  1104      Low   111      B       plan
#> 1105     1105 s1        1105 s1  1105      Low   111      B    monitor
#> 1106     1106 s1        1106 s1  1106      Low   111      B  consensus
#> 1107     1107 s1        1107 s1  1107      Low   111      B coregulate
#> 1108     1108 s1        1108 s1  1108      Low   111      B    discuss
#> 1109     1109 s1        1109 s1  1109      Low   111      B    discuss
#> 1110     1110 s1        1110 s1  1110      Low   111      B       plan
#> 1111     1111 s1        1111 s1  1111      Low   112      B    discuss
#> 1112     1112 s1        1112 s1  1112      Low   112      B    discuss
#> 1113     1113 s1        1113 s1  1113      Low   112      B   cohesion
#> 1114     1114 s1        1114 s1  1114      Low   112      B       plan
#> 1115     1115 s1        1115 s1  1115      Low   112      B    emotion
#> 1116     1116 s1        1116 s1  1116      Low   112      B  consensus
#> 1117     1117 s1        1117 s1  1117      Low   112      B  consensus
#> 1118     1118 s1        1118 s1  1118      Low   112      B    emotion
#> 1119     1119 s1        1119 s1  1119      Low   112      B    monitor
#> 1120     1120 s1        1120 s1  1120      Low   112      B       plan
#> 1121     1121 s1        1121 s1  1121      Low   113      B coregulate
#> 1122     1122 s1        1122 s1  1122      Low   113      B coregulate
#> 1123     1123 s1        1123 s1  1123      Low   113      B  consensus
#> 1124     1124 s1        1124 s1  1124      Low   113      B coregulate
#> 1125     1125 s1        1125 s1  1125      Low   113      B    discuss
#> 1126     1126 s1        1126 s1  1126      Low   113      B       plan
#> 1127     1127 s1        1127 s1  1127      Low   113      B    monitor
#> 1128     1128 s1        1128 s1  1128      Low   113      B    monitor
#> 1129     1129 s1        1129 s1  1129      Low   113      B    monitor
#> 1130     1130 s1        1130 s1  1130      Low   113      B       plan
#> 1131     1131 s1        1131 s1  1131      Low   114      B  consensus
#> 1132     1132 s1        1132 s1  1132      Low   114      B  consensus
#> 1133     1133 s1        1133 s1  1133      Low   114      B  consensus
#> 1134     1134 s1        1134 s1  1134      Low   114      B       plan
#> 1135     1135 s1        1135 s1  1135      Low   114      B    monitor
#> 1136     1136 s1        1136 s1  1136      Low   114      B  consensus
#> 1137     1137 s1        1137 s1  1137      Low   114      B    discuss
#> 1138     1138 s1        1138 s1  1138      Low   114      B  consensus
#> 1139     1139 s1        1139 s1  1139      Low   114      B    emotion
#> 1140     1140 s1        1140 s1  1140      Low   114      B  consensus
#> 1141     1141 s1        1141 s1  1141      Low   115      B       plan
#> 1142     1142 s1        1142 s1  1142      Low   115      B    emotion
#> 1143     1143 s1        1143 s1  1143      Low   115      B    discuss
#> 1144     1144 s1        1144 s1  1144      Low   115      B    discuss
#> 1145     1145 s1        1145 s1  1145      Low   115      B       plan
#> 1146     1146 s1        1146 s1  1146      Low   115      B  consensus
#> 1147     1147 s1        1147 s1  1147      Low   115      B  consensus
#> 1148     1148 s1        1148 s1  1148      Low   115      B       plan
#> 1149     1149 s1        1149 s1  1149      Low   115      B       plan
#> 1150     1150 s1        1150 s1  1150      Low   115      B    emotion
#> 1151     1151 s1        1151 s1  1151      Low   116      B   cohesion
#> 1152     1152 s1        1152 s1  1152      Low   116      B    emotion
#> 1153     1153 s1        1153 s1  1153      Low   116      B       plan
#> 1154     1154 s1        1154 s1  1154      Low   116      B  consensus
#> 1155     1155 s1        1155 s1  1155      Low   116      B    discuss
#> 1156     1156 s1        1156 s1  1156      Low   116      B  consensus
#> 1157     1157 s1        1157 s1  1157      Low   116      B    emotion
#> 1158     1158 s1        1158 s1  1158      Low   116      B   cohesion
#> 1159     1159 s1        1159 s1  1159      Low   116      B  consensus
#> 1160     1160 s1        1160 s1  1160      Low   116      B  consensus
#> 1161     1161 s1        1161 s1  1161      Low   117      B  consensus
#> 1162     1162 s1        1162 s1  1162      Low   117      B    discuss
#> 1163     1163 s1        1163 s1  1163      Low   117      B coregulate
#> 1164     1164 s1        1164 s1  1164      Low   117      B    emotion
#> 1165     1165 s1        1165 s1  1165      Low   117      B  consensus
#> 1166     1166 s1        1166 s1  1166      Low   117      B    discuss
#> 1167     1167 s1        1167 s1  1167      Low   117      B    monitor
#> 1168     1168 s1        1168 s1  1168      Low   117      B  synthesis
#> 1169     1169 s1        1169 s1  1169      Low   117      B    monitor
#> 1170     1170 s1        1170 s1  1170      Low   117      B    discuss
#> 1171     1171 s1        1171 s1  1171      Low   118      B    discuss
#> 1172     1172 s1        1172 s1  1172      Low   118      B       plan
#> 1173     1173 s1        1173 s1  1173      Low   118      B       plan
#> 1174     1174 s1        1174 s1  1174      Low   118      B    monitor
#> 1175     1175 s1        1175 s1  1175      Low   118      B    discuss
#> 1176     1176 s1        1176 s1  1176      Low   118      B  synthesis
#> 1177     1177 s1        1177 s1  1177      Low   118      B       plan
#> 1178     1178 s1        1178 s1  1178      Low   118      B    emotion
#> 1179     1179 s1        1179 s1  1179      Low   118      B       plan
#> 1180     1180 s1        1180 s1  1180      Low   118      B    discuss
#> 1181     1181 s1        1181 s1  1181      Low   119      B   cohesion
#> 1182     1182 s1        1182 s1  1182      Low   119      B       plan
#> 1183     1183 s1        1183 s1  1183      Low   119      B    monitor
#> 1184     1184 s1        1184 s1  1184      Low   119      B  consensus
#> 1185     1185 s1        1185 s1  1185      Low   119      B    monitor
#> 1186     1186 s1        1186 s1  1186      Low   119      B    monitor
#> 1187     1187 s1        1187 s1  1187      Low   119      B  consensus
#> 1188     1188 s1        1188 s1  1188      Low   119      B    monitor
#> 1189     1189 s1        1189 s1  1189      Low   119      B    monitor
#> 1190     1190 s1        1190 s1  1190      Low   119      B       plan
#> 1191     1191 s1        1191 s1  1191      Low   120      B       plan
#> 1192     1192 s1        1192 s1  1192      Low   120      B  consensus
#> 1193     1193 s1        1193 s1  1193      Low   120      B    discuss
#> 1194     1194 s1        1194 s1  1194      Low   120      B    discuss
#> 1195     1195 s1        1195 s1  1195      Low   120      B    monitor
#> 1196     1196 s1        1196 s1  1196      Low   120      B       plan
#> 1197     1197 s1        1197 s1  1197      Low   120      B  consensus
#> 1198     1198 s1        1198 s1  1198      Low   120      B    monitor
#> 1199     1199 s1        1199 s1  1199      Low   120      B    monitor
#> 1200     1200 s1        1200 s1  1200      Low   120      B    discuss
#> 1201     1201 s1        1201 s1  1201      Low   121      B  consensus
#> 1202     1202 s1        1202 s1  1202      Low   121      B       plan
#> 1203     1203 s1        1203 s1  1203      Low   121      B    monitor
#> 1204     1204 s1        1204 s1  1204      Low   121      B       plan
#> 1205     1205 s1        1205 s1  1205      Low   121      B    emotion
#> 1206     1206 s1        1206 s1  1206      Low   121      B    monitor
#> 1207     1207 s1        1207 s1  1207      Low   121      B       plan
#> 1208     1208 s1        1208 s1  1208      Low   121      B  synthesis
#> 1209     1209 s1        1209 s1  1209      Low   121      B       plan
#> 1210     1210 s1        1210 s1  1210      Low   121      B       plan
#> 1211     1211 s1        1211 s1  1211      Low   122      B       plan
#> 1212     1212 s1        1212 s1  1212      Low   122      B       plan
#> 1213     1213 s1        1213 s1  1213      Low   122      B    discuss
#> 1214     1214 s1        1214 s1  1214      Low   122      B       plan
#> 1215     1215 s1        1215 s1  1215      Low   122      B  consensus
#> 1216     1216 s1        1216 s1  1216      Low   122      B    emotion
#> 1217     1217 s1        1217 s1  1217      Low   122      B    discuss
#> 1218     1218 s1        1218 s1  1218      Low   122      B  consensus
#> 1219     1219 s1        1219 s1  1219      Low   122      B    emotion
#> 1220     1220 s1        1220 s1  1220      Low   122      B    emotion
#> 1221     1221 s1        1221 s1  1221      Low   123      B    emotion
#> 1222     1222 s1        1222 s1  1222      Low   123      B    monitor
#> 1223     1223 s1        1223 s1  1223      Low   123      B    monitor
#> 1224     1224 s1        1224 s1  1224      Low   123      B    discuss
#> 1225     1225 s1        1225 s1  1225      Low   123      B    emotion
#> 1226     1226 s1        1226 s1  1226      Low   123      B    emotion
#> 1227     1227 s1        1227 s1  1227      Low   123      B    emotion
#> 1228     1228 s1        1228 s1  1228      Low   123      B       plan
#> 1229     1229 s1        1229 s1  1229      Low   123      B  synthesis
#> 1230     1230 s1        1230 s1  1230      Low   123      B  consensus
#> 1231     1231 s1        1231 s1  1231      Low   124      B    emotion
#> 1232     1232 s1        1232 s1  1232      Low   124      B    monitor
#> 1233     1233 s1        1233 s1  1233      Low   124      B  consensus
#> 1234     1234 s1        1234 s1  1234      Low   124      B  consensus
#> 1235     1235 s1        1235 s1  1235      Low   124      B    discuss
#> 1236     1236 s1        1236 s1  1236      Low   124      B       plan
#> 1237     1237 s1        1237 s1  1237      Low   124      B coregulate
#> 1238     1238 s1        1238 s1  1238      Low   124      B    monitor
#> 1239     1239 s1        1239 s1  1239      Low   124      B    emotion
#> 1240     1240 s1        1240 s1  1240      Low   124      B  consensus
#> 1241     1241 s1        1241 s1  1241      Low   125      B    emotion
#> 1242     1242 s1        1242 s1  1242      Low   125      B    monitor
#> 1243     1243 s1        1243 s1  1243      Low   125      B   cohesion
#> 1244     1244 s1        1244 s1  1244      Low   125      B coregulate
#> 1245     1245 s1        1245 s1  1245      Low   125      B      adapt
#> 1246     1246 s1        1246 s1  1246      Low   125      B  consensus
#> 1247     1247 s1        1247 s1  1247      Low   125      B    monitor
#> 1248     1248 s1        1248 s1  1248      Low   125      B    emotion
#> 1249     1249 s1        1249 s1  1249      Low   125      B    discuss
#> 1250     1250 s1        1250 s1  1250      Low   125      B    monitor
#> 1251     1251 s1        1251 s1  1251      Low   126      B  consensus
#> 1252     1252 s1        1252 s1  1252      Low   126      B       plan
#> 1253     1253 s1        1253 s1  1253      Low   126      B       plan
#> 1254     1254 s1        1254 s1  1254      Low   126      B       plan
#> 1255     1255 s1        1255 s1  1255      Low   126      B  consensus
#> 1256     1256 s1        1256 s1  1256      Low   126      B  consensus
#> 1257     1257 s1        1257 s1  1257      Low   126      B    discuss
#> 1258     1258 s1        1258 s1  1258      Low   126      B       plan
#> 1259     1259 s1        1259 s1  1259      Low   126      B    monitor
#> 1260     1260 s1        1260 s1  1260      Low   126      B    monitor
#> 1261     1261 s1        1261 s1  1261      Low   127      B    emotion
#> 1262     1262 s1        1262 s1  1262      Low   127      B  consensus
#> 1263     1263 s1        1263 s1  1263      Low   127      B    emotion
#> 1264     1264 s1        1264 s1  1264      Low   127      B       plan
#> 1265     1265 s1        1265 s1  1265      Low   127      B    emotion
#> 1266     1266 s1        1266 s1  1266      Low   127      B    monitor
#> 1267     1267 s1        1267 s1  1267      Low   127      B    discuss
#> 1268     1268 s1        1268 s1  1268      Low   127      B    discuss
#> 1269     1269 s1        1269 s1  1269      Low   127      B    monitor
#> 1270     1270 s1        1270 s1  1270      Low   127      B  consensus
#> 1271     1271 s1        1271 s1  1271      Low   128      B       plan
#> 1272     1272 s1        1272 s1  1272      Low   128      B  consensus
#> 1273     1273 s1        1273 s1  1273      Low   128      B       plan
#> 1274     1274 s1        1274 s1  1274      Low   128      B    discuss
#> 1275     1275 s1        1275 s1  1275      Low   128      B   cohesion
#> 1276     1276 s1        1276 s1  1276      Low   128      B    discuss
#> 1277     1277 s1        1277 s1  1277      Low   128      B       plan
#> 1278     1278 s1        1278 s1  1278      Low   128      B  synthesis
#> 1279     1279 s1        1279 s1  1279      Low   128      B       plan
#> 1280     1280 s1        1280 s1  1280      Low   128      B  consensus
#> 1281     1281 s1        1281 s1  1281      Low   129      B       plan
#> 1282     1282 s1        1282 s1  1282      Low   129      B    monitor
#> 1283     1283 s1        1283 s1  1283      Low   129      B    discuss
#> 1284     1284 s1        1284 s1  1284      Low   129      B  synthesis
#> 1285     1285 s1        1285 s1  1285      Low   129      B    monitor
#> 1286     1286 s1        1286 s1  1286      Low   129      B    monitor
#> 1287     1287 s1        1287 s1  1287      Low   129      B    discuss
#> 1288     1288 s1        1288 s1  1288      Low   129      B    discuss
#> 1289     1289 s1        1289 s1  1289      Low   129      B    discuss
#> 1290     1290 s1        1290 s1  1290      Low   129      B       plan
#> 1291     1291 s1        1291 s1  1291      Low   130      B       plan
#> 1292     1292 s1        1292 s1  1292      Low   130      B  consensus
#> 1293     1293 s1        1293 s1  1293      Low   130      B       plan
#> 1294     1294 s1        1294 s1  1294      Low   130      B       plan
#> 1295     1295 s1        1295 s1  1295      Low   130      B    discuss
#> 1296     1296 s1        1296 s1  1296      Low   130      B       plan
#> 1297     1297 s1        1297 s1  1297      Low   130      B       plan
#> 1298     1298 s1        1298 s1  1298      Low   130      B  consensus
#> 1299     1299 s1        1299 s1  1299      Low   130      B    monitor
#> 1300     1300 s1        1300 s1  1300      Low   130      B    monitor
#> 1301     1301 s1        1301 s1  1301      Low   131      B    discuss
#> 1302     1302 s1        1302 s1  1302      Low   131      B    monitor
#> 1303     1303 s1        1303 s1  1303      Low   131      B    monitor
#> 1304     1304 s1        1304 s1  1304      Low   131      B  consensus
#> 1305     1305 s1        1305 s1  1305      Low   131      B    discuss
#> 1306     1306 s1        1306 s1  1306      Low   131      B    monitor
#> 1307     1307 s1        1307 s1  1307      Low   131      B  consensus
#> 1308     1308 s1        1308 s1  1308      Low   131      B   cohesion
#> 1309     1309 s1        1309 s1  1309      Low   131      B    monitor
#> 1310     1310 s1        1310 s1  1310      Low   131      B  consensus
#> 1311     1311 s1        1311 s1  1311      Low   132      B       plan
#> 1312     1312 s1        1312 s1  1312      Low   132      B       plan
#> 1313     1313 s1        1313 s1  1313      Low   132      B  consensus
#> 1314     1314 s1        1314 s1  1314      Low   132      B  consensus
#> 1315     1315 s1        1315 s1  1315      Low   132      B   cohesion
#> 1316     1316 s1        1316 s1  1316      Low   132      B       plan
#> 1317     1317 s1        1317 s1  1317      Low   132      B   cohesion
#> 1318     1318 s1        1318 s1  1318      Low   132      B  consensus
#> 1319     1319 s1        1319 s1  1319      Low   132      B  consensus
#> 1320     1320 s1        1320 s1  1320      Low   132      B    monitor
#> 1321     1321 s1        1321 s1  1321      Low   133      B    monitor
#> 1322     1322 s1        1322 s1  1322      Low   133      B    monitor
#> 1323     1323 s1        1323 s1  1323      Low   133      B       plan
#> 1324     1324 s1        1324 s1  1324      Low   133      B  consensus
#> 1325     1325 s1        1325 s1  1325      Low   133      B    emotion
#> 1326     1326 s1        1326 s1  1326      Low   133      B coregulate
#> 1327     1327 s1        1327 s1  1327      Low   133      B  consensus
#> 1328     1328 s1        1328 s1  1328      Low   133      B       plan
#> 1329     1329 s1        1329 s1  1329      Low   133      B coregulate
#> 1330     1330 s1        1330 s1  1330      Low   133      B    monitor
#> 1331     1331 s1        1331 s1  1331      Low   134      B       plan
#> 1332     1332 s1        1332 s1  1332      Low   134      B    discuss
#> 1333     1333 s1        1333 s1  1333      Low   134      B  consensus
#> 1334     1334 s1        1334 s1  1334      Low   134      B    emotion
#> 1335     1335 s1        1335 s1  1335      Low   134      B    monitor
#> 1336     1336 s1        1336 s1  1336      Low   134      B    emotion
#> 1337     1337 s1        1337 s1  1337      Low   134      B    emotion
#> 1338     1338 s1        1338 s1  1338      Low   134      B    monitor
#> 1339     1339 s1        1339 s1  1339      Low   134      B    emotion
#> 1340     1340 s1        1340 s1  1340      Low   134      B  consensus
#> 1341     1341 s1        1341 s1  1341      Low   135      B    discuss
#> 1342     1342 s1        1342 s1  1342      Low   135      B    emotion
#> 1343     1343 s1        1343 s1  1343      Low   135      B       plan
#> 1344     1344 s1        1344 s1  1344      Low   135      B    discuss
#> 1345     1345 s1        1345 s1  1345      Low   135      B   cohesion
#> 1346     1346 s1        1346 s1  1346      Low   135      B       plan
#> 1347     1347 s1        1347 s1  1347      Low   135      B    discuss
#> 1348     1348 s1        1348 s1  1348      Low   135      B    discuss
#> 1349     1349 s1        1349 s1  1349      Low   135      B    discuss
#> 1350     1350 s1        1350 s1  1350      Low   135      B  consensus
#> 1351     1351 s1        1351 s1  1351      Low   136      B    monitor
#> 1352     1352 s1        1352 s1  1352      Low   136      B  consensus
#> 1353     1353 s1        1353 s1  1353      Low   136      B    monitor
#> 1354     1354 s1        1354 s1  1354      Low   136      B    discuss
#> 1355     1355 s1        1355 s1  1355      Low   136      B    discuss
#> 1356     1356 s1        1356 s1  1356      Low   136      B       plan
#> 1357     1357 s1        1357 s1  1357      Low   136      B       plan
#> 1358     1358 s1        1358 s1  1358      Low   136      B    discuss
#> 1359     1359 s1        1359 s1  1359      Low   136      B    emotion
#> 1360     1360 s1        1360 s1  1360      Low   136      B    discuss
#> 1361     1361 s1        1361 s1  1361      Low   137      B       plan
#> 1362     1362 s1        1362 s1  1362      Low   137      B       plan
#> 1363     1363 s1        1363 s1  1363      Low   137      B       plan
#> 1364     1364 s1        1364 s1  1364      Low   137      B       plan
#> 1365     1365 s1        1365 s1  1365      Low   137      B    monitor
#> 1366     1366 s1        1366 s1  1366      Low   137      B      adapt
#> 1367     1367 s1        1367 s1  1367      Low   137      B coregulate
#> 1368     1368 s1        1368 s1  1368      Low   137      B    monitor
#> 1369     1369 s1        1369 s1  1369      Low   137      B  consensus
#> 1370     1370 s1        1370 s1  1370      Low   137      B  consensus
#> 1371     1371 s1        1371 s1  1371      Low   138      B    discuss
#> 1372     1372 s1        1372 s1  1372      Low   138      B  consensus
#> 1373     1373 s1        1373 s1  1373      Low   138      B   cohesion
#> 1374     1374 s1        1374 s1  1374      Low   138      B  consensus
#> 1375     1375 s1        1375 s1  1375      Low   138      B  consensus
#> 1376     1376 s1        1376 s1  1376      Low   138      B coregulate
#> 1377     1377 s1        1377 s1  1377      Low   138      B    monitor
#> 1378     1378 s1        1378 s1  1378      Low   138      B    discuss
#> 1379     1379 s1        1379 s1  1379      Low   138      B    discuss
#> 1380     1380 s1        1380 s1  1380      Low   138      B  consensus
#> 1381     1381 s1        1381 s1  1381      Low   139      B  consensus
#> 1382     1382 s1        1382 s1  1382      Low   139      B    discuss
#> 1383     1383 s1        1383 s1  1383      Low   139      B  consensus
#> 1384     1384 s1        1384 s1  1384      Low   139      B    discuss
#> 1385     1385 s1        1385 s1  1385      Low   139      B  consensus
#> 1386     1386 s1        1386 s1  1386      Low   139      B       plan
#> 1387     1387 s1        1387 s1  1387      Low   139      B    discuss
#> 1388     1388 s1        1388 s1  1388      Low   139      B  consensus
#> 1389     1389 s1        1389 s1  1389      Low   139      B    monitor
#> 1390     1390 s1        1390 s1  1390      Low   139      B       plan
#> 1391     1391 s1        1391 s1  1391      Low   140      B    emotion
#> 1392     1392 s1        1392 s1  1392      Low   140      B       plan
#> 1393     1393 s1        1393 s1  1393      Low   140      B       plan
#> 1394     1394 s1        1394 s1  1394      Low   140      B   cohesion
#> 1395     1395 s1        1395 s1  1395      Low   140      B    emotion
#> 1396     1396 s1        1396 s1  1396      Low   140      B       plan
#> 1397     1397 s1        1397 s1  1397      Low   140      B    emotion
#> 1398     1398 s1        1398 s1  1398      Low   140      B    emotion
#> 1399     1399 s1        1399 s1  1399      Low   140      B    monitor
#> 1400     1400 s1        1400 s1  1400      Low   140      B    discuss
#> 1401     1401 s1        1401 s1  1401      Low   141      B    emotion
#> 1402     1402 s1        1402 s1  1402      Low   141      B coregulate
#> 1403     1403 s1        1403 s1  1403      Low   141      B    discuss
#> 1404     1404 s1        1404 s1  1404      Low   141      B    discuss
#> 1405     1405 s1        1405 s1  1405      Low   141      B  consensus
#> 1406     1406 s1        1406 s1  1406      Low   141      B    monitor
#> 1407     1407 s1        1407 s1  1407      Low   141      B    emotion
#> 1408     1408 s1        1408 s1  1408      Low   141      B  consensus
#> 1409     1409 s1        1409 s1  1409      Low   141      B   cohesion
#> 1410     1410 s1        1410 s1  1410      Low   141      B  consensus
#> 1411     1411 s1        1411 s1  1411      Low   142      B    monitor
#> 1412     1412 s1        1412 s1  1412      Low   142      B       plan
#> 1413     1413 s1        1413 s1  1413      Low   142      B  consensus
#> 1414     1414 s1        1414 s1  1414      Low   142      B    emotion
#> 1415     1415 s1        1415 s1  1415      Low   142      B    monitor
#> 1416     1416 s1        1416 s1  1416      Low   142      B    discuss
#> 1417     1417 s1        1417 s1  1417      Low   142      B   cohesion
#> 1418     1418 s1        1418 s1  1418      Low   142      B    monitor
#> 1419     1419 s1        1419 s1  1419      Low   142      B    discuss
#> 1420     1420 s1        1420 s1  1420      Low   142      B    monitor
#> 1421     1421 s1        1421 s1  1421      Low   143      B   cohesion
#> 1422     1422 s1        1422 s1  1422      Low   143      B       plan
#> 1423     1423 s1        1423 s1  1423      Low   143      B       plan
#> 1424     1424 s1        1424 s1  1424      Low   143      B    discuss
#> 1425     1425 s1        1425 s1  1425      Low   143      B  consensus
#> 1426     1426 s1        1426 s1  1426      Low   143      B  consensus
#> 1427     1427 s1        1427 s1  1427      Low   143      B    emotion
#> 1428     1428 s1        1428 s1  1428      Low   143      B    monitor
#> 1429     1429 s1        1429 s1  1429      Low   143      B    monitor
#> 1430     1430 s1        1430 s1  1430      Low   143      B coregulate
#> 1431     1431 s1        1431 s1  1431      Low   144      B    monitor
#> 1432     1432 s1        1432 s1  1432      Low   144      B  consensus
#> 1433     1433 s1        1433 s1  1433      Low   144      B coregulate
#> 1434     1434 s1        1434 s1  1434      Low   144      B  consensus
#> 1435     1435 s1        1435 s1  1435      Low   144      B  consensus
#> 1436     1436 s1        1436 s1  1436      Low   144      B    monitor
#> 1437     1437 s1        1437 s1  1437      Low   144      B   cohesion
#> 1438     1438 s1        1438 s1  1438      Low   144      B  consensus
#> 1439     1439 s1        1439 s1  1439      Low   144      B  consensus
#> 1440     1440 s1        1440 s1  1440      Low   144      B  consensus
#> 1441     1441 s1        1441 s1  1441      Low   145      B       plan
#> 1442     1442 s1        1442 s1  1442      Low   145      B       plan
#> 1443     1443 s1        1443 s1  1443      Low   145      B coregulate
#> 1444     1444 s1        1444 s1  1444      Low   145      B    emotion
#> 1445     1445 s1        1445 s1  1445      Low   145      B  synthesis
#> 1446     1446 s1        1446 s1  1446      Low   145      B  consensus
#> 1447     1447 s1        1447 s1  1447      Low   145      B    monitor
#> 1448     1448 s1        1448 s1  1448      Low   145      B    discuss
#> 1449     1449 s1        1449 s1  1449      Low   145      B    emotion
#> 1450     1450 s1        1450 s1  1450      Low   145      B       plan
#> 1451     1451 s1        1451 s1  1451      Low   146      B       plan
#> 1452     1452 s1        1452 s1  1452      Low   146      B    discuss
#> 1453     1453 s1        1453 s1  1453      Low   146      B    emotion
#> 1454     1454 s1        1454 s1  1454      Low   146      B    monitor
#> 1455     1455 s1        1455 s1  1455      Low   146      B    emotion
#> 1456     1456 s1        1456 s1  1456      Low   146      B    discuss
#> 1457     1457 s1        1457 s1  1457      Low   146      B  consensus
#> 1458     1458 s1        1458 s1  1458      Low   146      B       plan
#> 1459     1459 s1        1459 s1  1459      Low   146      B    discuss
#> 1460     1460 s1        1460 s1  1460      Low   146      B  consensus
#> 1461     1461 s1        1461 s1  1461      Low   147      B       plan
#> 1462     1462 s1        1462 s1  1462      Low   147      B       plan
#> 1463     1463 s1        1463 s1  1463      Low   147      B    discuss
#> 1464     1464 s1        1464 s1  1464      Low   147      B  consensus
#> 1465     1465 s1        1465 s1  1465      Low   147      B    emotion
#> 1466     1466 s1        1466 s1  1466      Low   147      B    monitor
#> 1467     1467 s1        1467 s1  1467      Low   147      B       plan
#> 1468     1468 s1        1468 s1  1468      Low   147      B    discuss
#> 1469     1469 s1        1469 s1  1469      Low   147      B    emotion
#> 1470     1470 s1        1470 s1  1470      Low   147      B  consensus
#> 1471     1471 s1        1471 s1  1471      Low   148      B       plan
#> 1472     1472 s1        1472 s1  1472      Low   148      B    discuss
#> 1473     1473 s1        1473 s1  1473      Low   148      B    discuss
#> 1474     1474 s1        1474 s1  1474      Low   148      B    emotion
#> 1475     1475 s1        1475 s1  1475      Low   148      B    emotion
#> 1476     1476 s1        1476 s1  1476      Low   148      B    emotion
#> 1477     1477 s1        1477 s1  1477      Low   148      B       plan
#> 1478     1478 s1        1478 s1  1478      Low   148      B    emotion
#> 1479     1479 s1        1479 s1  1479      Low   148      B  consensus
#> 1480     1480 s1        1480 s1  1480      Low   148      B       plan
#> 1481     1481 s1        1481 s1  1481      Low   149      B  consensus
#> 1482     1482 s1        1482 s1  1482      Low   149      B    monitor
#> 1483     1483 s1        1483 s1  1483      Low   149      B  consensus
#> 1484     1484 s1        1484 s1  1484      Low   149      B       plan
#> 1485     1485 s1        1485 s1  1485      Low   149      B    monitor
#> 1486     1486 s1        1486 s1  1486      Low   149      B       plan
#> 1487     1487 s1        1487 s1  1487      Low   149      B       plan
#> 1488     1488 s1        1488 s1  1488      Low   149      B    monitor
#> 1489     1489 s1        1489 s1  1489      Low   149      B   cohesion
#> 1490     1490 s1        1490 s1  1490      Low   149      B    emotion
#> 1491     1491 s1        1491 s1  1491      Low   150      B    discuss
#> 1492     1492 s1        1492 s1  1492      Low   150      B      adapt
#> 1493     1493 s1        1493 s1  1493      Low   150      B    discuss
#> 1494     1494 s1        1494 s1  1494      Low   150      B    emotion
#> 1495     1495 s1        1495 s1  1495      Low   150      B  consensus
#> 1496     1496 s1        1496 s1  1496      Low   150      B    emotion
#> 1497     1497 s1        1497 s1  1497      Low   150      B  consensus
#> 1498     1498 s1        1498 s1  1498      Low   150      B       plan
#> 1499     1499 s1        1499 s1  1499      Low   150      B    emotion
#> 1500     1500 s1        1500 s1  1500      Low   150      B  consensus
#> 1501     1501 s1        1501 s1  1501      Low   151      B  consensus
#> 1502     1502 s1        1502 s1  1502      Low   151      B       plan
#> 1503     1503 s1        1503 s1  1503      Low   151      B    emotion
#> 1504     1504 s1        1504 s1  1504      Low   151      B  consensus
#> 1505     1505 s1        1505 s1  1505      Low   151      B    monitor
#> 1506     1506 s1        1506 s1  1506      Low   151      B       plan
#> 1507     1507 s1        1507 s1  1507      Low   151      B    emotion
#> 1508     1508 s1        1508 s1  1508      Low   151      B   cohesion
#> 1509     1509 s1        1509 s1  1509      Low   151      B       plan
#> 1510     1510 s1        1510 s1  1510      Low   151      B    monitor
#> 1511     1511 s1        1511 s1  1511      Low   152      B    monitor
#> 1512     1512 s1        1512 s1  1512      Low   152      B  consensus
#> 1513     1513 s1        1513 s1  1513      Low   152      B    emotion
#> 1514     1514 s1        1514 s1  1514      Low   152      B    monitor
#> 1515     1515 s1        1515 s1  1515      Low   152      B    discuss
#> 1516     1516 s1        1516 s1  1516      Low   152      B       plan
#> 1517     1517 s1        1517 s1  1517      Low   152      B  consensus
#> 1518     1518 s1        1518 s1  1518      Low   152      B    discuss
#> 1519     1519 s1        1519 s1  1519      Low   152      B       plan
#> 1520     1520 s1        1520 s1  1520      Low   152      B  consensus
#> 1521     1521 s1        1521 s1  1521      Low   153      B       plan
#> 1522     1522 s1        1522 s1  1522      Low   153      B    emotion
#> 1523     1523 s1        1523 s1  1523      Low   153      B    discuss
#> 1524     1524 s1        1524 s1  1524      Low   153      B    monitor
#> 1525     1525 s1        1525 s1  1525      Low   153      B    discuss
#> 1526     1526 s1        1526 s1  1526      Low   153      B  consensus
#> 1527     1527 s1        1527 s1  1527      Low   153      B       plan
#> 1528     1528 s1        1528 s1  1528      Low   153      B    discuss
#> 1529     1529 s1        1529 s1  1529      Low   153      B    discuss
#> 1530     1530 s1        1530 s1  1530      Low   153      B    monitor
#> 1531     1531 s1        1531 s1  1531      Low   154      B    discuss
#> 1532     1532 s1        1532 s1  1532      Low   154      B coregulate
#> 1533     1533 s1        1533 s1  1533      Low   154      B    discuss
#> 1534     1534 s1        1534 s1  1534      Low   154      B    discuss
#> 1535     1535 s1        1535 s1  1535      Low   154      B    monitor
#> 1536     1536 s1        1536 s1  1536      Low   154      B    discuss
#> 1537     1537 s1        1537 s1  1537      Low   154      B    emotion
#> 1538     1538 s1        1538 s1  1538      Low   154      B  synthesis
#> 1539     1539 s1        1539 s1  1539      Low   154      B  consensus
#> 1540     1540 s1        1540 s1  1540      Low   154      B  consensus
#> 1541     1541 s1        1541 s1  1541      Low   155      B    monitor
#> 1542     1542 s1        1542 s1  1542      Low   155      B    emotion
#> 1543     1543 s1        1543 s1  1543      Low   155      B  consensus
#> 1544     1544 s1        1544 s1  1544      Low   155      B  consensus
#> 1545     1545 s1        1545 s1  1545      Low   155      B       plan
#> 1546     1546 s1        1546 s1  1546      Low   155      B  consensus
#> 1547     1547 s1        1547 s1  1547      Low   155      B    emotion
#> 1548     1548 s1        1548 s1  1548      Low   155      B       plan
#> 1549     1549 s1        1549 s1  1549      Low   155      B  consensus
#> 1550     1550 s1        1550 s1  1550      Low   155      B       plan
#> 1551     1551 s1        1551 s1  1551      Low   156      B   cohesion
#> 1552     1552 s1        1552 s1  1552      Low   156      B       plan
#> 1553     1553 s1        1553 s1  1553      Low   156      B   cohesion
#> 1554     1554 s1        1554 s1  1554      Low   156      B      adapt
#> 1555     1555 s1        1555 s1  1555      Low   156      B    discuss
#> 1556     1556 s1        1556 s1  1556      Low   156      B       plan
#> 1557     1557 s1        1557 s1  1557      Low   156      B    monitor
#> 1558     1558 s1        1558 s1  1558      Low   156      B    discuss
#> 1559     1559 s1        1559 s1  1559      Low   156      B    emotion
#> 1560     1560 s1        1560 s1  1560      Low   156      B    monitor
#> 1561     1561 s1        1561 s1  1561      Low   157      B    emotion
#> 1562     1562 s1        1562 s1  1562      Low   157      B    emotion
#> 1563     1563 s1        1563 s1  1563      Low   157      B    discuss
#> 1564     1564 s1        1564 s1  1564      Low   157      B       plan
#> 1565     1565 s1        1565 s1  1565      Low   157      B coregulate
#> 1566     1566 s1        1566 s1  1566      Low   157      B    monitor
#> 1567     1567 s1        1567 s1  1567      Low   157      B    discuss
#> 1568     1568 s1        1568 s1  1568      Low   157      B  consensus
#> 1569     1569 s1        1569 s1  1569      Low   157      B       plan
#> 1570     1570 s1        1570 s1  1570      Low   157      B    emotion
#> 1571     1571 s1        1571 s1  1571      Low   158      B    discuss
#> 1572     1572 s1        1572 s1  1572      Low   158      B  consensus
#> 1573     1573 s1        1573 s1  1573      Low   158      B    discuss
#> 1574     1574 s1        1574 s1  1574      Low   158      B    emotion
#> 1575     1575 s1        1575 s1  1575      Low   158      B       plan
#> 1576     1576 s1        1576 s1  1576      Low   158      B    monitor
#> 1577     1577 s1        1577 s1  1577      Low   158      B       plan
#> 1578     1578 s1        1578 s1  1578      Low   158      B       plan
#> 1579     1579 s1        1579 s1  1579      Low   158      B    monitor
#> 1580     1580 s1        1580 s1  1580      Low   158      B    discuss
#> 1581     1581 s1        1581 s1  1581      Low   159      B    discuss
#> 1582     1582 s1        1582 s1  1582      Low   159      B       plan
#> 1583     1583 s1        1583 s1  1583      Low   159      B    discuss
#> 1584     1584 s1        1584 s1  1584      Low   159      B    discuss
#> 1585     1585 s1        1585 s1  1585      Low   159      B       plan
#> 1586     1586 s1        1586 s1  1586      Low   159      B  consensus
#> 1587     1587 s1        1587 s1  1587      Low   159      B       plan
#> 1588     1588 s1        1588 s1  1588      Low   159      B       plan
#> 1589     1589 s1        1589 s1  1589      Low   159      B   cohesion
#> 1590     1590 s1        1590 s1  1590      Low   159      B    discuss
#> 1591     1591 s1        1591 s1  1591      Low   160      B    discuss
#> 1592     1592 s1        1592 s1  1592      Low   160      B       plan
#> 1593     1593 s1        1593 s1  1593      Low   160      B       plan
#> 1594     1594 s1        1594 s1  1594      Low   160      B    discuss
#> 1595     1595 s1        1595 s1  1595      Low   160      B  consensus
#> 1596     1596 s1        1596 s1  1596      Low   160      B    emotion
#> 1597     1597 s1        1597 s1  1597      Low   160      B  consensus
#> 1598     1598 s1        1598 s1  1598      Low   160      B    emotion
#> 1599     1599 s1        1599 s1  1599      Low   160      B       plan
#> 1600     1600 s1        1600 s1  1600      Low   160      B    monitor
#> 1601     1601 s1        1601 s1  1601      Low   161      C coregulate
#> 1602     1602 s1        1602 s1  1602      Low   161      C       plan
#> 1603     1603 s1        1603 s1  1603      Low   161      C    emotion
#> 1604     1604 s1        1604 s1  1604      Low   161      C  consensus
#> 1605     1605 s1        1605 s1  1605      Low   161      C  consensus
#> 1606     1606 s1        1606 s1  1606      Low   161      C      adapt
#> 1607     1607 s1        1607 s1  1607      Low   161      C  consensus
#> 1608     1608 s1        1608 s1  1608      Low   161      C    discuss
#> 1609     1609 s1        1609 s1  1609      Low   161      C    discuss
#> 1610     1610 s1        1610 s1  1610      Low   161      C  consensus
#> 1611     1611 s1        1611 s1  1611      Low   162      C       plan
#> 1612     1612 s1        1612 s1  1612      Low   162      C  consensus
#> 1613     1613 s1        1613 s1  1613      Low   162      C    emotion
#> 1614     1614 s1        1614 s1  1614      Low   162      C    discuss
#> 1615     1615 s1        1615 s1  1615      Low   162      C       plan
#> 1616     1616 s1        1616 s1  1616      Low   162      C       plan
#> 1617     1617 s1        1617 s1  1617      Low   162      C       plan
#> 1618     1618 s1        1618 s1  1618      Low   162      C  consensus
#> 1619     1619 s1        1619 s1  1619      Low   162      C    monitor
#> 1620     1620 s1        1620 s1  1620      Low   162      C   cohesion
#> 1621     1621 s1        1621 s1  1621      Low   163      C    emotion
#> 1622     1622 s1        1622 s1  1622      Low   163      C    monitor
#> 1623     1623 s1        1623 s1  1623      Low   163      C       plan
#> 1624     1624 s1        1624 s1  1624      Low   163      C       plan
#> 1625     1625 s1        1625 s1  1625      Low   163      C       plan
#> 1626     1626 s1        1626 s1  1626      Low   163      C       plan
#> 1627     1627 s1        1627 s1  1627      Low   163      C       plan
#> 1628     1628 s1        1628 s1  1628      Low   163      C  consensus
#> 1629     1629 s1        1629 s1  1629      Low   163      C   cohesion
#> 1630     1630 s1        1630 s1  1630      Low   163      C    monitor
#> 1631     1631 s1        1631 s1  1631      Low   164      C    discuss
#> 1632     1632 s1        1632 s1  1632      Low   164      C    monitor
#> 1633     1633 s1        1633 s1  1633      Low   164      C    discuss
#> 1634     1634 s1        1634 s1  1634      Low   164      C   cohesion
#> 1635     1635 s1        1635 s1  1635      Low   164      C       plan
#> 1636     1636 s1        1636 s1  1636      Low   164      C   cohesion
#> 1637     1637 s1        1637 s1  1637      Low   164      C  consensus
#> 1638     1638 s1        1638 s1  1638      Low   164      C       plan
#> 1639     1639 s1        1639 s1  1639      Low   164      C       plan
#> 1640     1640 s1        1640 s1  1640      Low   164      C    monitor
#> 1641     1641 s1        1641 s1  1641      Low   165      C    discuss
#> 1642     1642 s1        1642 s1  1642      Low   165      C       plan
#> 1643     1643 s1        1643 s1  1643      Low   165      C       plan
#> 1644     1644 s1        1644 s1  1644      Low   165      C    emotion
#> 1645     1645 s1        1645 s1  1645      Low   165      C  consensus
#> 1646     1646 s1        1646 s1  1646      Low   165      C    monitor
#> 1647     1647 s1        1647 s1  1647      Low   165      C       plan
#> 1648     1648 s1        1648 s1  1648      Low   165      C    emotion
#> 1649     1649 s1        1649 s1  1649      Low   165      C    discuss
#> 1650     1650 s1        1650 s1  1650      Low   165      C  consensus
#> 1651     1651 s1        1651 s1  1651      Low   166      C    emotion
#> 1652     1652 s1        1652 s1  1652      Low   166      C    discuss
#> 1653     1653 s1        1653 s1  1653      Low   166      C       plan
#> 1654     1654 s1        1654 s1  1654      Low   166      C    monitor
#> 1655     1655 s1        1655 s1  1655      Low   166      C       plan
#> 1656     1656 s1        1656 s1  1656      Low   166      C       plan
#> 1657     1657 s1        1657 s1  1657      Low   166      C       plan
#> 1658     1658 s1        1658 s1  1658      Low   166      C  consensus
#> 1659     1659 s1        1659 s1  1659      Low   166      C    discuss
#> 1660     1660 s1        1660 s1  1660      Low   166      C coregulate
#> 1661     1661 s1        1661 s1  1661      Low   167      C    monitor
#> 1662     1662 s1        1662 s1  1662      Low   167      C    discuss
#> 1663     1663 s1        1663 s1  1663      Low   167      C    discuss
#> 1664     1664 s1        1664 s1  1664      Low   167      C  consensus
#> 1665     1665 s1        1665 s1  1665      Low   167      C      adapt
#> 1666     1666 s1        1666 s1  1666      Low   167      C    emotion
#> 1667     1667 s1        1667 s1  1667      Low   167      C    monitor
#> 1668     1668 s1        1668 s1  1668      Low   167      C    emotion
#> 1669     1669 s1        1669 s1  1669      Low   167      C  synthesis
#> 1670     1670 s1        1670 s1  1670      Low   167      C    monitor
#> 1671     1671 s1        1671 s1  1671      Low   168      C       plan
#> 1672     1672 s1        1672 s1  1672      Low   168      C       plan
#> 1673     1673 s1        1673 s1  1673      Low   168      C    emotion
#> 1674     1674 s1        1674 s1  1674      Low   168      C  consensus
#> 1675     1675 s1        1675 s1  1675      Low   168      C  consensus
#> 1676     1676 s1        1676 s1  1676      Low   168      C    monitor
#> 1677     1677 s1        1677 s1  1677      Low   168      C    emotion
#> 1678     1678 s1        1678 s1  1678      Low   168      C    emotion
#> 1679     1679 s1        1679 s1  1679      Low   168      C    discuss
#> 1680     1680 s1        1680 s1  1680      Low   168      C       plan
#> 1681     1681 s1        1681 s1  1681      Low   169      C    discuss
#> 1682     1682 s1        1682 s1  1682      Low   169      C    monitor
#> 1683     1683 s1        1683 s1  1683      Low   169      C    monitor
#> 1684     1684 s1        1684 s1  1684      Low   169      C       plan
#> 1685     1685 s1        1685 s1  1685      Low   169      C       plan
#> 1686     1686 s1        1686 s1  1686      Low   169      C    monitor
#> 1687     1687 s1        1687 s1  1687      Low   169      C       plan
#> 1688     1688 s1        1688 s1  1688      Low   169      C    discuss
#> 1689     1689 s1        1689 s1  1689      Low   169      C       plan
#> 1690     1690 s1        1690 s1  1690      Low   169      C   cohesion
#> 1691     1691 s1        1691 s1  1691      Low   170      C    emotion
#> 1692     1692 s1        1692 s1  1692      Low   170      C  consensus
#> 1693     1693 s1        1693 s1  1693      Low   170      C    discuss
#> 1694     1694 s1        1694 s1  1694      Low   170      C    discuss
#> 1695     1695 s1        1695 s1  1695      Low   170      C    monitor
#> 1696     1696 s1        1696 s1  1696      Low   170      C  consensus
#> 1697     1697 s1        1697 s1  1697      Low   170      C       plan
#> 1698     1698 s1        1698 s1  1698      Low   170      C  synthesis
#> 1699     1699 s1        1699 s1  1699      Low   170      C    emotion
#> 1700     1700 s1        1700 s1  1700      Low   170      C  consensus
#> 1701     1701 s1        1701 s1  1701      Low   171      C    monitor
#> 1702     1702 s1        1702 s1  1702      Low   171      C  consensus
#> 1703     1703 s1        1703 s1  1703      Low   171      C    monitor
#> 1704     1704 s1        1704 s1  1704      Low   171      C  consensus
#> 1705     1705 s1        1705 s1  1705      Low   171      C       plan
#> 1706     1706 s1        1706 s1  1706      Low   171      C    monitor
#> 1707     1707 s1        1707 s1  1707      Low   171      C  consensus
#> 1708     1708 s1        1708 s1  1708      Low   171      C       plan
#> 1709     1709 s1        1709 s1  1709      Low   171      C       plan
#> 1710     1710 s1        1710 s1  1710      Low   171      C       plan
#> 1711     1711 s1        1711 s1  1711      Low   172      C  consensus
#> 1712     1712 s1        1712 s1  1712      Low   172      C       plan
#> 1713     1713 s1        1713 s1  1713      Low   172      C       plan
#> 1714     1714 s1        1714 s1  1714      Low   172      C  consensus
#> 1715     1715 s1        1715 s1  1715      Low   172      C      adapt
#> 1716     1716 s1        1716 s1  1716      Low   172      C    discuss
#> 1717     1717 s1        1717 s1  1717      Low   172      C      adapt
#> 1718     1718 s1        1718 s1  1718      Low   172      C    discuss
#> 1719     1719 s1        1719 s1  1719      Low   172      C       plan
#> 1720     1720 s1        1720 s1  1720      Low   172      C coregulate
#> 1721     1721 s1        1721 s1  1721      Low   173      C      adapt
#> 1722     1722 s1        1722 s1  1722      Low   173      C  consensus
#> 1723     1723 s1        1723 s1  1723      Low   173      C  consensus
#> 1724     1724 s1        1724 s1  1724      Low   173      C       plan
#> 1725     1725 s1        1725 s1  1725      Low   173      C       plan
#> 1726     1726 s1        1726 s1  1726      Low   173      C    monitor
#> 1727     1727 s1        1727 s1  1727      Low   173      C    discuss
#> 1728     1728 s1        1728 s1  1728      Low   173      C       plan
#> 1729     1729 s1        1729 s1  1729      Low   173      C  consensus
#> 1730     1730 s1        1730 s1  1730      Low   173      C    emotion
#> 1731     1731 s1        1731 s1  1731      Low   174      C    emotion
#> 1732     1732 s1        1732 s1  1732      Low   174      C   cohesion
#> 1733     1733 s1        1733 s1  1733      Low   174      C    discuss
#> 1734     1734 s1        1734 s1  1734      Low   174      C  consensus
#> 1735     1735 s1        1735 s1  1735      Low   174      C    discuss
#> 1736     1736 s1        1736 s1  1736      Low   174      C    discuss
#> 1737     1737 s1        1737 s1  1737      Low   174      C  consensus
#> 1738     1738 s1        1738 s1  1738      Low   174      C  consensus
#> 1739     1739 s1        1739 s1  1739      Low   174      C       plan
#> 1740     1740 s1        1740 s1  1740      Low   174      C  consensus
#> 1741     1741 s1        1741 s1  1741      Low   175      C    discuss
#> 1742     1742 s1        1742 s1  1742      Low   175      C    monitor
#> 1743     1743 s1        1743 s1  1743      Low   175      C    discuss
#> 1744     1744 s1        1744 s1  1744      Low   175      C    monitor
#> 1745     1745 s1        1745 s1  1745      Low   175      C    emotion
#> 1746     1746 s1        1746 s1  1746      Low   175      C       plan
#> 1747     1747 s1        1747 s1  1747      Low   175      C  consensus
#> 1748     1748 s1        1748 s1  1748      Low   175      C    monitor
#> 1749     1749 s1        1749 s1  1749      Low   175      C    discuss
#> 1750     1750 s1        1750 s1  1750      Low   175      C  consensus
#> 1751     1751 s1        1751 s1  1751      Low   176      C    discuss
#> 1752     1752 s1        1752 s1  1752      Low   176      C  consensus
#> 1753     1753 s1        1753 s1  1753      Low   176      C       plan
#> 1754     1754 s1        1754 s1  1754      Low   176      C    monitor
#> 1755     1755 s1        1755 s1  1755      Low   176      C    discuss
#> 1756     1756 s1        1756 s1  1756      Low   176      C    monitor
#> 1757     1757 s1        1757 s1  1757      Low   176      C       plan
#> 1758     1758 s1        1758 s1  1758      Low   176      C coregulate
#> 1759     1759 s1        1759 s1  1759      Low   176      C    discuss
#> 1760     1760 s1        1760 s1  1760      Low   176      C coregulate
#> 1761     1761 s1        1761 s1  1761      Low   177      C  consensus
#> 1762     1762 s1        1762 s1  1762      Low   177      C    monitor
#> 1763     1763 s1        1763 s1  1763      Low   177      C  consensus
#> 1764     1764 s1        1764 s1  1764      Low   177      C    discuss
#> 1765     1765 s1        1765 s1  1765      Low   177      C  consensus
#> 1766     1766 s1        1766 s1  1766      Low   177      C       plan
#> 1767     1767 s1        1767 s1  1767      Low   177      C    discuss
#> 1768     1768 s1        1768 s1  1768      Low   177      C       plan
#> 1769     1769 s1        1769 s1  1769      Low   177      C    monitor
#> 1770     1770 s1        1770 s1  1770      Low   177      C   cohesion
#> 1771     1771 s1        1771 s1  1771      Low   178      C       plan
#> 1772     1772 s1        1772 s1  1772      Low   178      C  consensus
#> 1773     1773 s1        1773 s1  1773      Low   178      C    emotion
#> 1774     1774 s1        1774 s1  1774      Low   178      C  consensus
#> 1775     1775 s1        1775 s1  1775      Low   178      C    emotion
#> 1776     1776 s1        1776 s1  1776      Low   178      C  consensus
#> 1777     1777 s1        1777 s1  1777      Low   178      C    monitor
#> 1778     1778 s1        1778 s1  1778      Low   178      C    monitor
#> 1779     1779 s1        1779 s1  1779      Low   178      C    monitor
#> 1780     1780 s1        1780 s1  1780      Low   178      C  consensus
#> 1781     1781 s1        1781 s1  1781      Low   179      C    monitor
#> 1782     1782 s1        1782 s1  1782      Low   179      C    discuss
#> 1783     1783 s1        1783 s1  1783      Low   179      C coregulate
#> 1784     1784 s1        1784 s1  1784      Low   179      C  synthesis
#> 1785     1785 s1        1785 s1  1785      Low   179      C    emotion
#> 1786     1786 s1        1786 s1  1786      Low   179      C  consensus
#> 1787     1787 s1        1787 s1  1787      Low   179      C  consensus
#> 1788     1788 s1        1788 s1  1788      Low   179      C  consensus
#> 1789     1789 s1        1789 s1  1789      Low   179      C       plan
#> 1790     1790 s1        1790 s1  1790      Low   179      C    discuss
#> 1791     1791 s1        1791 s1  1791      Low   180      C    discuss
#> 1792     1792 s1        1792 s1  1792      Low   180      C    monitor
#> 1793     1793 s1        1793 s1  1793      Low   180      C    monitor
#> 1794     1794 s1        1794 s1  1794      Low   180      C    monitor
#> 1795     1795 s1        1795 s1  1795      Low   180      C      adapt
#> 1796     1796 s1        1796 s1  1796      Low   180      C       plan
#> 1797     1797 s1        1797 s1  1797      Low   180      C    discuss
#> 1798     1798 s1        1798 s1  1798      Low   180      C    discuss
#> 1799     1799 s1        1799 s1  1799      Low   180      C    discuss
#> 1800     1800 s1        1800 s1  1800      Low   180      C       plan
#> 1801     1801 s1        1801 s1  1801      Low   181      C       plan
#> 1802     1802 s1        1802 s1  1802      Low   181      C    emotion
#> 1803     1803 s1        1803 s1  1803      Low   181      C   cohesion
#> 1804     1804 s1        1804 s1  1804      Low   181      C    emotion
#> 1805     1805 s1        1805 s1  1805      Low   181      C    discuss
#> 1806     1806 s1        1806 s1  1806      Low   181      C  consensus
#> 1807     1807 s1        1807 s1  1807      Low   181      C    monitor
#> 1808     1808 s1        1808 s1  1808      Low   181      C  consensus
#> 1809     1809 s1        1809 s1  1809      Low   181      C   cohesion
#> 1810     1810 s1        1810 s1  1810      Low   181      C    emotion
#> 1811     1811 s1        1811 s1  1811      Low   182      C    discuss
#> 1812     1812 s1        1812 s1  1812      Low   182      C       plan
#> 1813     1813 s1        1813 s1  1813      Low   182      C    discuss
#> 1814     1814 s1        1814 s1  1814      Low   182      C    monitor
#> 1815     1815 s1        1815 s1  1815      Low   182      C    monitor
#> 1816     1816 s1        1816 s1  1816      Low   182      C       plan
#> 1817     1817 s1        1817 s1  1817      Low   182      C       plan
#> 1818     1818 s1        1818 s1  1818      Low   182      C coregulate
#> 1819     1819 s1        1819 s1  1819      Low   182      C       plan
#> 1820     1820 s1        1820 s1  1820      Low   182      C  consensus
#> 1821     1821 s1        1821 s1  1821      Low   183      C  consensus
#> 1822     1822 s1        1822 s1  1822      Low   183      C  consensus
#> 1823     1823 s1        1823 s1  1823      Low   183      C       plan
#> 1824     1824 s1        1824 s1  1824      Low   183      C       plan
#> 1825     1825 s1        1825 s1  1825      Low   183      C  consensus
#> 1826     1826 s1        1826 s1  1826      Low   183      C       plan
#> 1827     1827 s1        1827 s1  1827      Low   183      C       plan
#> 1828     1828 s1        1828 s1  1828      Low   183      C    monitor
#> 1829     1829 s1        1829 s1  1829      Low   183      C  consensus
#> 1830     1830 s1        1830 s1  1830      Low   183      C    emotion
#> 1831     1831 s1        1831 s1  1831      Low   184      C  synthesis
#> 1832     1832 s1        1832 s1  1832      Low   184      C    emotion
#> 1833     1833 s1        1833 s1  1833      Low   184      C    monitor
#> 1834     1834 s1        1834 s1  1834      Low   184      C    monitor
#> 1835     1835 s1        1835 s1  1835      Low   184      C  consensus
#> 1836     1836 s1        1836 s1  1836      Low   184      C    emotion
#> 1837     1837 s1        1837 s1  1837      Low   184      C    discuss
#> 1838     1838 s1        1838 s1  1838      Low   184      C       plan
#> 1839     1839 s1        1839 s1  1839      Low   184      C    discuss
#> 1840     1840 s1        1840 s1  1840      Low   184      C  consensus
#> 1841     1841 s1        1841 s1  1841      Low   185      C    monitor
#> 1842     1842 s1        1842 s1  1842      Low   185      C  consensus
#> 1843     1843 s1        1843 s1  1843      Low   185      C    emotion
#> 1844     1844 s1        1844 s1  1844      Low   185      C  consensus
#> 1845     1845 s1        1845 s1  1845      Low   185      C       plan
#> 1846     1846 s1        1846 s1  1846      Low   185      C  consensus
#> 1847     1847 s1        1847 s1  1847      Low   185      C       plan
#> 1848     1848 s1        1848 s1  1848      Low   185      C       plan
#> 1849     1849 s1        1849 s1  1849      Low   185      C    monitor
#> 1850     1850 s1        1850 s1  1850      Low   185      C       plan
#> 1851     1851 s1        1851 s1  1851      Low   186      C       plan
#> 1852     1852 s1        1852 s1  1852      Low   186      C    emotion
#> 1853     1853 s1        1853 s1  1853      Low   186      C  consensus
#> 1854     1854 s1        1854 s1  1854      Low   186      C    discuss
#> 1855     1855 s1        1855 s1  1855      Low   186      C    monitor
#> 1856     1856 s1        1856 s1  1856      Low   186      C    emotion
#> 1857     1857 s1        1857 s1  1857      Low   186      C    discuss
#> 1858     1858 s1        1858 s1  1858      Low   186      C       plan
#> 1859     1859 s1        1859 s1  1859      Low   186      C    discuss
#> 1860     1860 s1        1860 s1  1860      Low   186      C       plan
#> 1861     1861 s1        1861 s1  1861      Low   187      C  consensus
#> 1862     1862 s1        1862 s1  1862      Low   187      C  consensus
#> 1863     1863 s1        1863 s1  1863      Low   187      C       plan
#> 1864     1864 s1        1864 s1  1864      Low   187      C       plan
#> 1865     1865 s1        1865 s1  1865      Low   187      C    emotion
#> 1866     1866 s1        1866 s1  1866      Low   187      C    emotion
#> 1867     1867 s1        1867 s1  1867      Low   187      C    discuss
#> 1868     1868 s1        1868 s1  1868      Low   187      C    emotion
#> 1869     1869 s1        1869 s1  1869      Low   187      C    emotion
#> 1870     1870 s1        1870 s1  1870      Low   187      C    monitor
#> 1871     1871 s1        1871 s1  1871      Low   188      C  consensus
#> 1872     1872 s1        1872 s1  1872      Low   188      C    monitor
#> 1873     1873 s1        1873 s1  1873      Low   188      C    emotion
#> 1874     1874 s1        1874 s1  1874      Low   188      C       plan
#> 1875     1875 s1        1875 s1  1875      Low   188      C    emotion
#> 1876     1876 s1        1876 s1  1876      Low   188      C    discuss
#> 1877     1877 s1        1877 s1  1877      Low   188      C  consensus
#> 1878     1878 s1        1878 s1  1878      Low   188      C    monitor
#> 1879     1879 s1        1879 s1  1879      Low   188      C    monitor
#> 1880     1880 s1        1880 s1  1880      Low   188      C       plan
#> 1881     1881 s1        1881 s1  1881      Low   189      C  consensus
#> 1882     1882 s1        1882 s1  1882      Low   189      C    monitor
#> 1883     1883 s1        1883 s1  1883      Low   189      C       plan
#> 1884     1884 s1        1884 s1  1884      Low   189      C       plan
#> 1885     1885 s1        1885 s1  1885      Low   189      C coregulate
#> 1886     1886 s1        1886 s1  1886      Low   189      C       plan
#> 1887     1887 s1        1887 s1  1887      Low   189      C       plan
#> 1888     1888 s1        1888 s1  1888      Low   189      C   cohesion
#> 1889     1889 s1        1889 s1  1889      Low   189      C    monitor
#> 1890     1890 s1        1890 s1  1890      Low   189      C  consensus
#> 1891     1891 s1        1891 s1  1891      Low   190      C       plan
#> 1892     1892 s1        1892 s1  1892      Low   190      C       plan
#> 1893     1893 s1        1893 s1  1893      Low   190      C    emotion
#> 1894     1894 s1        1894 s1  1894      Low   190      C  consensus
#> 1895     1895 s1        1895 s1  1895      Low   190      C   cohesion
#> 1896     1896 s1        1896 s1  1896      Low   190      C       plan
#> 1897     1897 s1        1897 s1  1897      Low   190      C  consensus
#> 1898     1898 s1        1898 s1  1898      Low   190      C  consensus
#> 1899     1899 s1        1899 s1  1899      Low   190      C    monitor
#> 1900     1900 s1        1900 s1  1900      Low   190      C    discuss
#> 1901     1901 s1        1901 s1  1901      Low   191      C    discuss
#> 1902     1902 s1        1902 s1  1902      Low   191      C    emotion
#> 1903     1903 s1        1903 s1  1903      Low   191      C       plan
#> 1904     1904 s1        1904 s1  1904      Low   191      C  consensus
#> 1905     1905 s1        1905 s1  1905      Low   191      C    discuss
#> 1906     1906 s1        1906 s1  1906      Low   191      C    monitor
#> 1907     1907 s1        1907 s1  1907      Low   191      C    discuss
#> 1908     1908 s1        1908 s1  1908      Low   191      C    monitor
#> 1909     1909 s1        1909 s1  1909      Low   191      C  consensus
#> 1910     1910 s1        1910 s1  1910      Low   191      C       plan
#> 1911     1911 s1        1911 s1  1911      Low   192      C    monitor
#> 1912     1912 s1        1912 s1  1912      Low   192      C    monitor
#> 1913     1913 s1        1913 s1  1913      Low   192      C    emotion
#> 1914     1914 s1        1914 s1  1914      Low   192      C  consensus
#> 1915     1915 s1        1915 s1  1915      Low   192      C       plan
#> 1916     1916 s1        1916 s1  1916      Low   192      C  consensus
#> 1917     1917 s1        1917 s1  1917      Low   192      C  consensus
#> 1918     1918 s1        1918 s1  1918      Low   192      C  consensus
#> 1919     1919 s1        1919 s1  1919      Low   192      C    monitor
#> 1920     1920 s1        1920 s1  1920      Low   192      C    emotion
#> 1921     1921 s1        1921 s1  1921      Low   193      C    monitor
#> 1922     1922 s1        1922 s1  1922      Low   193      C    discuss
#> 1923     1923 s1        1923 s1  1923      Low   193      C    discuss
#> 1924     1924 s1        1924 s1  1924      Low   193      C       plan
#> 1925     1925 s1        1925 s1  1925      Low   193      C    emotion
#> 1926     1926 s1        1926 s1  1926      Low   193      C       plan
#> 1927     1927 s1        1927 s1  1927      Low   193      C       plan
#> 1928     1928 s1        1928 s1  1928      Low   193      C  consensus
#> 1929     1929 s1        1929 s1  1929      Low   193      C    monitor
#> 1930     1930 s1        1930 s1  1930      Low   193      C  consensus
#> 1931     1931 s1        1931 s1  1931      Low   194      C    discuss
#> 1932     1932 s1        1932 s1  1932      Low   194      C       plan
#> 1933     1933 s1        1933 s1  1933      Low   194      C  consensus
#> 1934     1934 s1        1934 s1  1934      Low   194      C  synthesis
#> 1935     1935 s1        1935 s1  1935      Low   194      C  consensus
#> 1936     1936 s1        1936 s1  1936      Low   194      C    emotion
#> 1937     1937 s1        1937 s1  1937      Low   194      C    emotion
#> 1938     1938 s1        1938 s1  1938      Low   194      C       plan
#> 1939     1939 s1        1939 s1  1939      Low   194      C    monitor
#> 1940     1940 s1        1940 s1  1940      Low   194      C    discuss
#> 1941     1941 s1        1941 s1  1941      Low   195      C   cohesion
#> 1942     1942 s1        1942 s1  1942      Low   195      C    discuss
#> 1943     1943 s1        1943 s1  1943      Low   195      C  consensus
#> 1944     1944 s1        1944 s1  1944      Low   195      C    emotion
#> 1945     1945 s1        1945 s1  1945      Low   195      C  consensus
#> 1946     1946 s1        1946 s1  1946      Low   195      C    emotion
#> 1947     1947 s1        1947 s1  1947      Low   195      C    monitor
#> 1948     1948 s1        1948 s1  1948      Low   195      C coregulate
#> 1949     1949 s1        1949 s1  1949      Low   195      C    emotion
#> 1950     1950 s1        1950 s1  1950      Low   195      C    emotion
#> 1951     1951 s1        1951 s1  1951      Low   196      C  consensus
#> 1952     1952 s1        1952 s1  1952      Low   196      C       plan
#> 1953     1953 s1        1953 s1  1953      Low   196      C       plan
#> 1954     1954 s1        1954 s1  1954      Low   196      C    monitor
#> 1955     1955 s1        1955 s1  1955      Low   196      C    monitor
#> 1956     1956 s1        1956 s1  1956      Low   196      C coregulate
#> 1957     1957 s1        1957 s1  1957      Low   196      C    monitor
#> 1958     1958 s1        1958 s1  1958      Low   196      C   cohesion
#> 1959     1959 s1        1959 s1  1959      Low   196      C       plan
#> 1960     1960 s1        1960 s1  1960      Low   196      C    monitor
#> 1961     1961 s1        1961 s1  1961      Low   197      C coregulate
#> 1962     1962 s1        1962 s1  1962      Low   197      C    monitor
#> 1963     1963 s1        1963 s1  1963      Low   197      C  consensus
#> 1964     1964 s1        1964 s1  1964      Low   197      C    discuss
#> 1965     1965 s1        1965 s1  1965      Low   197      C    discuss
#> 1966     1966 s1        1966 s1  1966      Low   197      C    discuss
#> 1967     1967 s1        1967 s1  1967      Low   197      C  consensus
#> 1968     1968 s1        1968 s1  1968      Low   197      C  consensus
#> 1969     1969 s1        1969 s1  1969      Low   197      C coregulate
#> 1970     1970 s1        1970 s1  1970      Low   197      C  consensus
#> 1971     1971 s1        1971 s1  1971      Low   198      C    emotion
#> 1972     1972 s1        1972 s1  1972      Low   198      C    emotion
#> 1973     1973 s1        1973 s1  1973      Low   198      C   cohesion
#> 1974     1974 s1        1974 s1  1974      Low   198      C    emotion
#> 1975     1975 s1        1975 s1  1975      Low   198      C    emotion
#> 1976     1976 s1        1976 s1  1976      Low   198      C    emotion
#> 1977     1977 s1        1977 s1  1977      Low   198      C    monitor
#> 1978     1978 s1        1978 s1  1978      Low   198      C    emotion
#> 1979     1979 s1        1979 s1  1979      Low   198      C       plan
#> 1980     1980 s1        1980 s1  1980      Low   198      C  consensus
#> 1981     1981 s1        1981 s1  1981      Low   199      C  synthesis
#> 1982     1982 s1        1982 s1  1982      Low   199      C    monitor
#> 1983     1983 s1        1983 s1  1983      Low   199      C    discuss
#> 1984     1984 s1        1984 s1  1984      Low   199      C    emotion
#> 1985     1985 s1        1985 s1  1985      Low   199      C       plan
#> 1986     1986 s1        1986 s1  1986      Low   199      C    monitor
#> 1987     1987 s1        1987 s1  1987      Low   199      C  consensus
#> 1988     1988 s1        1988 s1  1988      Low   199      C    monitor
#> 1989     1989 s1        1989 s1  1989      Low   199      C       plan
#> 1990     1990 s1        1990 s1  1990      Low   199      C    monitor
#> 1991     1991 s1        1991 s1  1991      Low   200      C    emotion
#> 1992     1992 s1        1992 s1  1992      Low   200      C    discuss
#> 1993     1993 s1        1993 s1  1993      Low   200      C    discuss
#> 1994     1994 s1        1994 s1  1994      Low   200      C    monitor
#> 1995     1995 s1        1995 s1  1995      Low   200      C  consensus
#> 1996     1996 s1        1996 s1  1996      Low   200      C    monitor
#> 1997     1997 s1        1997 s1  1997      Low   200      C       plan
#> 1998     1998 s1        1998 s1  1998      Low   200      C       plan
#> 1999     1999 s1        1999 s1  1999      Low   200      C  consensus
#> 2000     2000 s1        2000 s1  2000      Low   200      C       plan
#>              T2         T3         T4         T5         T6         T7
#> 1     consensus    discuss  synthesis      adapt  consensus       plan
#> 2       emotion  consensus    discuss  synthesis      adapt    emotion
#> 3    coregulate    monitor  consensus       plan    emotion  consensus
#> 4       emotion       plan    discuss  synthesis  consensus    discuss
#> 5       emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 6          plan  consensus       plan       plan       plan       plan
#> 7       discuss coregulate       <NA>       <NA>       <NA>       <NA>
#> 8       emotion  consensus    discuss  consensus       plan  consensus
#> 9     consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 10     cohesion    discuss  synthesis       <NA>       <NA>       <NA>
#> 11         plan       plan  consensus    discuss  consensus coregulate
#> 12    consensus       plan   cohesion  consensus    discuss  synthesis
#> 13    consensus       plan       plan    emotion    monitor       plan
#> 14      monitor       <NA>       <NA>       <NA>       <NA>       <NA>
#> 15      emotion       plan       plan    monitor       plan  consensus
#> 16    consensus    discuss    emotion  consensus   cohesion coregulate
#> 17    consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 18   coregulate       plan    emotion    discuss  consensus coregulate
#> 19    consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 20    consensus       plan       plan    emotion  consensus       plan
#> 21      emotion    discuss  consensus coregulate    emotion  consensus
#> 22      emotion  consensus       plan  consensus  consensus       plan
#> 23    consensus    monitor       plan       plan  consensus       plan
#> 24         plan       plan       <NA>       <NA>       <NA>       <NA>
#> 25         plan       plan   cohesion    emotion  consensus       plan
#> 26   coregulate    emotion       <NA>       <NA>       <NA>       <NA>
#> 27    consensus       plan    emotion  consensus    discuss coregulate
#> 28      emotion  consensus       plan  consensus       plan  consensus
#> 29      discuss  synthesis       <NA>       <NA>       <NA>       <NA>
#> 30         plan       plan       plan  consensus    discuss    discuss
#> 31    consensus       plan    monitor    monitor  consensus coregulate
#> 32      emotion   cohesion  consensus       plan  consensus  consensus
#> 33      emotion    discuss    discuss  synthesis  consensus coregulate
#> 34    consensus       plan  consensus coregulate    monitor       plan
#> 35    consensus       plan    discuss  consensus  consensus    discuss
#> 36      discuss  consensus    discuss  synthesis  consensus      adapt
#> 37    consensus  consensus    discuss   cohesion  consensus coregulate
#> 38    synthesis  consensus       plan    discuss    discuss  consensus
#> 39    consensus       plan    emotion   cohesion  consensus   cohesion
#> 40      emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 41      emotion  consensus  consensus  consensus    discuss  consensus
#> 42     cohesion       plan  consensus    emotion    discuss coregulate
#> 43    synthesis      adapt       <NA>       <NA>       <NA>       <NA>
#> 44      discuss  consensus       plan       plan       plan  consensus
#> 45    consensus coregulate    emotion       plan  consensus    discuss
#> 46         plan coregulate  consensus       plan coregulate    emotion
#> 47   coregulate    discuss      adapt  consensus    discuss   cohesion
#> 48         plan      adapt  consensus    discuss  consensus       plan
#> 49    consensus    emotion  consensus       plan    discuss    emotion
#> 50         plan    discuss  consensus       plan       plan       plan
#> 51         plan  consensus coregulate       plan    monitor  consensus
#> 52    consensus       plan    emotion   cohesion    emotion  consensus
#> 53      emotion   cohesion    emotion  consensus coregulate    discuss
#> 54    consensus    discuss    discuss  consensus      adapt       <NA>
#> 55    consensus    discuss    emotion   cohesion       plan       plan
#> 56     cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 57      monitor  consensus       plan       plan    monitor    discuss
#> 58     cohesion  consensus    discuss    discuss  consensus       plan
#> 59   coregulate       plan    emotion    discuss  consensus       plan
#> 60     cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 61      monitor  consensus       <NA>       <NA>       <NA>       <NA>
#> 62    consensus    discuss  consensus       <NA>       <NA>       <NA>
#> 63      emotion  consensus       plan       plan    monitor    discuss
#> 64      discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 65      emotion  consensus       plan       plan       plan       <NA>
#> 66         plan       plan    discuss  consensus    discuss    discuss
#> 67    consensus       plan    discuss  consensus  consensus       plan
#> 68         plan       plan    emotion    discuss    emotion   cohesion
#> 69      discuss   cohesion  consensus    discuss  synthesis  consensus
#> 70      emotion    emotion    emotion  consensus       <NA>       <NA>
#> 71      discuss  consensus    discuss    discuss    discuss  consensus
#> 72      emotion   cohesion  consensus       plan       plan    emotion
#> 73         plan       plan  consensus    discuss    discuss       <NA>
#> 74    consensus    monitor    discuss  synthesis       plan       plan
#> 75    consensus    discuss  synthesis  consensus    discuss       plan
#> 76   coregulate    emotion   cohesion       <NA>       <NA>       <NA>
#> 77    consensus       plan  consensus       plan    emotion    discuss
#> 78    consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 79    consensus       plan    monitor       plan       plan    emotion
#> 80         plan       plan       <NA>       <NA>       <NA>       <NA>
#> 81    consensus    emotion    emotion   cohesion    emotion  consensus
#> 82        adapt  consensus coregulate    discuss  consensus       <NA>
#> 83    consensus    discuss coregulate       plan  consensus    discuss
#> 84    consensus   cohesion  consensus       plan       plan       plan
#> 85   coregulate    discuss  consensus       plan    emotion   cohesion
#> 86         plan       plan       <NA>       <NA>       <NA>       <NA>
#> 87      discuss  synthesis       plan  consensus coregulate       plan
#> 88      discuss  consensus    discuss  consensus    emotion  consensus
#> 89      monitor       plan  consensus       <NA>       <NA>       <NA>
#> 90         plan       plan       plan    monitor       plan    emotion
#> 91   coregulate    emotion   cohesion  consensus       plan       plan
#> 92    consensus    emotion    discuss  consensus       plan  consensus
#> 93      discuss  consensus       plan  synthesis  consensus    discuss
#> 94    consensus       plan coregulate    monitor    discuss       <NA>
#> 95      emotion  consensus       plan    discuss       plan       plan
#> 96    consensus    discuss  consensus    discuss    discuss    discuss
#> 97         plan  consensus       plan  consensus  consensus coregulate
#> 98         plan  consensus       plan  consensus    monitor    discuss
#> 99    consensus       plan       plan    emotion   cohesion  consensus
#> 100   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 101     discuss  consensus    discuss  synthesis  consensus       plan
#> 102   consensus coregulate    discuss  consensus    discuss    discuss
#> 103   synthesis      adapt       <NA>       <NA>       <NA>       <NA>
#> 104    cohesion  consensus       plan       <NA>       <NA>       <NA>
#> 105        plan  consensus    emotion coregulate    discuss  synthesis
#> 106  coregulate  consensus    emotion  consensus coregulate coregulate
#> 107        plan    emotion   cohesion  consensus  consensus  consensus
#> 108   consensus       plan       plan    monitor    discuss    emotion
#> 109   consensus    discuss  consensus  consensus       plan   cohesion
#> 110        plan       plan    emotion  consensus       plan       plan
#> 111    cohesion  consensus    discuss  consensus       plan    discuss
#> 112   consensus   cohesion       plan       plan       plan       plan
#> 113     emotion  consensus       plan    discuss  consensus    emotion
#> 114  coregulate       plan    emotion   cohesion  consensus       plan
#> 115   consensus       plan    emotion   cohesion  consensus coregulate
#> 116   consensus coregulate    emotion   cohesion    emotion   cohesion
#> 117     discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 118        plan    emotion  consensus    discuss  consensus  consensus
#> 119   consensus    discuss    discuss  consensus  consensus       plan
#> 120     discuss  consensus coregulate       plan       plan    emotion
#> 121   consensus    emotion       <NA>       <NA>       <NA>       <NA>
#> 122   consensus    discuss   cohesion    monitor    discuss    discuss
#> 123     emotion  consensus coregulate    discuss       plan    monitor
#> 124     emotion   cohesion  consensus coregulate    discuss       <NA>
#> 125    cohesion  consensus coregulate    emotion   cohesion       plan
#> 126    cohesion       plan    monitor coregulate    discuss    emotion
#> 127        plan       plan  consensus    discuss  consensus       plan
#> 128     emotion    discuss    monitor    discuss  consensus coregulate
#> 129        plan  consensus    discuss    emotion  consensus       plan
#> 130   consensus       plan    emotion  consensus  consensus    monitor
#> 131    cohesion       plan       plan    emotion    emotion  consensus
#> 132   consensus       plan       plan  synthesis  consensus  consensus
#> 133   consensus       plan  consensus  consensus       <NA>       <NA>
#> 134     discuss  consensus       plan    monitor  consensus       plan
#> 135   consensus    discuss  consensus       plan coregulate       plan
#> 136     discuss    monitor    discuss    discuss    emotion  consensus
#> 137     discuss  synthesis  consensus       plan       plan       plan
#> 138   consensus  consensus    discuss coregulate    emotion   cohesion
#> 139  coregulate      adapt    monitor  consensus    discuss    discuss
#> 140        plan    discuss    emotion  consensus   cohesion  consensus
#> 141     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 142     monitor  consensus   cohesion coregulate    monitor    discuss
#> 143  coregulate    discuss  consensus coregulate  consensus    discuss
#> 144    cohesion coregulate  consensus    discuss    emotion   cohesion
#> 145   consensus coregulate       plan  consensus coregulate    emotion
#> 146  coregulate    discuss coregulate    discuss    emotion  consensus
#> 147        plan  consensus       plan       plan       plan    emotion
#> 148     discuss  consensus       plan  consensus    emotion  consensus
#> 149  coregulate   cohesion coregulate       <NA>       <NA>       <NA>
#> 150   consensus    discuss    discuss    discuss    discuss       plan
#> 151    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 152  coregulate    monitor       plan       plan       plan       plan
#> 153     discuss  synthesis       <NA>       <NA>       <NA>       <NA>
#> 154   synthesis       plan    emotion    emotion    emotion  consensus
#> 155    cohesion  consensus  consensus       plan       plan  consensus
#> 156   consensus       plan    monitor    discuss       plan       plan
#> 157     monitor  consensus       plan  consensus    discuss   cohesion
#> 158     discuss    emotion    emotion   cohesion       plan    emotion
#> 159  coregulate    emotion       plan  consensus coregulate    emotion
#> 160    cohesion  consensus coregulate       plan       plan    emotion
#> 161   consensus  consensus  consensus    monitor    discuss coregulate
#> 162     emotion    discuss    discuss  consensus       plan  consensus
#> 163     discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 164   consensus  synthesis  consensus       plan       plan  consensus
#> 165   consensus    emotion  consensus    discuss  consensus    discuss
#> 166     emotion   cohesion    emotion      adapt  consensus    discuss
#> 167   consensus       plan  consensus       plan       plan    emotion
#> 168   consensus    discuss  synthesis  consensus coregulate    monitor
#> 169    cohesion  consensus    discuss  consensus coregulate    discuss
#> 170        plan       plan  consensus       <NA>       <NA>       <NA>
#> 171     monitor  consensus       plan  consensus    discuss       plan
#> 172        plan   cohesion    emotion    discuss  consensus coregulate
#> 173     emotion    discuss    discuss  consensus    discuss  consensus
#> 174   consensus       plan  consensus    monitor  consensus  consensus
#> 175   consensus       plan  consensus coregulate    emotion  consensus
#> 176        plan   cohesion    emotion   cohesion    emotion  consensus
#> 177     emotion   cohesion  consensus coregulate    discuss  consensus
#> 178   consensus    emotion  consensus  consensus coregulate    discuss
#> 179   consensus  consensus       <NA>       <NA>       <NA>       <NA>
#> 180   synthesis  consensus    emotion  consensus       plan    emotion
#> 181     monitor    discuss  consensus       plan  consensus coregulate
#> 182    cohesion   cohesion       plan  consensus   cohesion  consensus
#> 183     emotion       plan    emotion   cohesion  consensus  consensus
#> 184   consensus       plan  consensus coregulate   cohesion  consensus
#> 185     emotion      adapt  consensus coregulate       plan       plan
#> 186    cohesion  consensus       plan       plan       <NA>       <NA>
#> 187     emotion   cohesion    emotion   cohesion  consensus    discuss
#> 188        plan    discuss  consensus  consensus coregulate    discuss
#> 189   consensus       plan       plan       <NA>       <NA>       <NA>
#> 190   consensus       plan       plan  consensus       plan    discuss
#> 191   consensus    discuss  consensus       plan  consensus  consensus
#> 192   consensus    monitor       <NA>       <NA>       <NA>       <NA>
#> 193     discuss    discuss  synthesis       <NA>       <NA>       <NA>
#> 194     discuss  consensus       plan       plan    emotion       <NA>
#> 195        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 196        plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 197   synthesis   cohesion coregulate    monitor       plan    monitor
#> 198   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 199  coregulate    discuss  consensus    monitor   cohesion  consensus
#> 200     discuss  consensus       plan       plan       plan       plan
#> 201     monitor       plan       plan    monitor    discuss  consensus
#> 202     emotion    discuss  consensus       plan    emotion   cohesion
#> 203    cohesion    emotion       <NA>       <NA>       <NA>       <NA>
#> 204   consensus  consensus    discuss  consensus coregulate    discuss
#> 205     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 206     discuss    monitor       plan       plan  consensus    discuss
#> 207   consensus       plan       plan  consensus       plan    emotion
#> 208    cohesion  consensus       plan       plan    discuss       <NA>
#> 209        plan  consensus       plan coregulate    monitor       plan
#> 210     discuss   cohesion   cohesion coregulate       plan   cohesion
#> 211    cohesion   cohesion    emotion    monitor       plan       plan
#> 212        plan  consensus coregulate       plan       plan  consensus
#> 213  coregulate    emotion       <NA>       <NA>       <NA>       <NA>
#> 214   consensus       plan   cohesion  consensus       plan       plan
#> 215   consensus    emotion    monitor    discuss       <NA>       <NA>
#> 216     discuss  consensus coregulate    emotion   cohesion       <NA>
#> 217       adapt  consensus coregulate    discuss    discuss  consensus
#> 218     discuss    discuss  consensus coregulate coregulate       plan
#> 219        plan       plan    discuss  consensus       plan       plan
#> 220     emotion   cohesion  consensus  consensus       <NA>       <NA>
#> 221     emotion   cohesion  consensus coregulate    discuss  consensus
#> 222   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 223        plan  consensus  consensus    emotion  consensus       plan
#> 224  coregulate    monitor    discuss  consensus    discuss       plan
#> 225     discuss    discuss      adapt    emotion   cohesion    emotion
#> 226    cohesion  consensus  consensus  consensus    discuss  consensus
#> 227        plan  consensus       plan    discuss  synthesis      adapt
#> 228  coregulate    discuss  consensus       plan    discuss    monitor
#> 229   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 230        plan   cohesion  consensus    discuss  consensus       plan
#> 231   consensus    discuss  consensus    discuss   cohesion    emotion
#> 232   consensus coregulate       plan  consensus  consensus  consensus
#> 233  coregulate    monitor       plan    emotion    discuss  synthesis
#> 234        plan  consensus       plan       plan  consensus    emotion
#> 235  coregulate       plan       plan  consensus  consensus  consensus
#> 236     discuss  consensus       plan coregulate       plan       plan
#> 237   consensus       plan    discuss    emotion    monitor       <NA>
#> 238     emotion  consensus coregulate    emotion   cohesion coregulate
#> 239        plan    emotion  consensus       plan    monitor       plan
#> 240   consensus    monitor      adapt   cohesion  consensus       plan
#> 241        plan    emotion   cohesion  consensus    emotion    discuss
#> 242     emotion       plan    monitor    emotion  consensus       plan
#> 243     emotion  consensus       plan       plan    monitor  consensus
#> 244   consensus       plan    emotion   cohesion  consensus    discuss
#> 245        plan    emotion  consensus coregulate    emotion   cohesion
#> 246   consensus       plan       plan  consensus    discuss  consensus
#> 247     emotion   cohesion    emotion  consensus  consensus       plan
#> 248    cohesion  consensus    discuss    discuss  consensus    emotion
#> 249    cohesion    emotion    discuss  synthesis  consensus       plan
#> 250        plan  consensus coregulate  consensus       plan    emotion
#> 251   consensus    discuss    discuss    discuss   cohesion  consensus
#> 252   consensus       plan    emotion   cohesion    emotion    discuss
#> 253     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 254    cohesion   cohesion  synthesis    discuss  consensus coregulate
#> 255   consensus    discuss   cohesion  consensus    discuss    discuss
#> 256    cohesion   cohesion       plan  consensus       plan  consensus
#> 257        plan       plan    emotion   cohesion  consensus coregulate
#> 258   consensus       plan       plan    discuss   cohesion  consensus
#> 259     monitor    discuss  consensus  consensus    discuss  consensus
#> 260     discuss    emotion       plan    emotion   cohesion    discuss
#> 261     discuss  synthesis       <NA>       <NA>       <NA>       <NA>
#> 262   consensus  consensus coregulate    discuss  synthesis  consensus
#> 263        plan    discuss    discuss  synthesis  consensus       plan
#> 264     discuss  consensus    discuss  consensus  consensus    monitor
#> 265     monitor  consensus coregulate       <NA>       <NA>       <NA>
#> 266   consensus       plan  consensus    discuss   cohesion  consensus
#> 267   consensus    emotion  consensus       plan  consensus    emotion
#> 268    cohesion  consensus       plan  consensus       plan    monitor
#> 269   consensus  consensus       plan    emotion   cohesion  consensus
#> 270   consensus    discuss  synthesis       plan       plan       plan
#> 271  coregulate       plan       plan    discuss       <NA>       <NA>
#> 272     emotion  consensus       plan   cohesion  consensus    discuss
#> 273   consensus  consensus    discuss  consensus    emotion    discuss
#> 274     monitor    discuss       <NA>       <NA>       <NA>       <NA>
#> 275   synthesis  consensus coregulate    discuss  consensus  consensus
#> 276  coregulate   cohesion  consensus    discuss  synthesis  consensus
#> 277    cohesion  consensus       plan    monitor    discuss    discuss
#> 278        plan    emotion  consensus  consensus  consensus       plan
#> 279   consensus    discuss coregulate    emotion  consensus    emotion
#> 280    cohesion       plan    emotion    discuss coregulate       <NA>
#> 281     emotion       plan  consensus coregulate       plan       plan
#> 282     discuss    emotion  consensus       plan       plan  consensus
#> 283    cohesion  consensus    discuss  consensus    emotion   cohesion
#> 284        plan  consensus       plan    monitor    discuss coregulate
#> 285     emotion   cohesion  consensus       plan       plan       plan
#> 286   consensus  consensus coregulate       plan  consensus  consensus
#> 287    cohesion  consensus    discuss  consensus       plan    emotion
#> 288  coregulate    emotion    discuss    discuss  synthesis  consensus
#> 289     discuss  consensus   cohesion coregulate    discuss  consensus
#> 290   consensus       plan  consensus coregulate       plan       <NA>
#> 291  coregulate  consensus coregulate    discuss  synthesis  consensus
#> 292        plan  consensus coregulate  consensus    discuss  consensus
#> 293        plan       plan    discuss  consensus coregulate    emotion
#> 294   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 295        plan    emotion  consensus       plan    emotion  consensus
#> 296        plan       plan    emotion   cohesion       plan    emotion
#> 297        plan  consensus  consensus       plan  consensus    discuss
#> 298     emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 299    cohesion  consensus coregulate       plan    discuss  consensus
#> 300    cohesion coregulate    emotion   cohesion  consensus    discuss
#> 301    cohesion  consensus    discuss  consensus       plan      adapt
#> 302     emotion  consensus    discuss    emotion  consensus       plan
#> 303    cohesion    emotion   cohesion       plan  consensus    discuss
#> 304     emotion   cohesion  consensus    emotion  consensus  consensus
#> 305        plan coregulate    emotion  consensus       plan  consensus
#> 306    cohesion    discuss  consensus    emotion  consensus    discuss
#> 307     emotion   cohesion  consensus  consensus    discuss   cohesion
#> 308     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 309  coregulate       plan  consensus coregulate    discuss  synthesis
#> 310     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 311     discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 312    cohesion    emotion    emotion  consensus    discuss      adapt
#> 313        plan  consensus  consensus    discuss  consensus       plan
#> 314    cohesion  consensus    discuss coregulate       <NA>       <NA>
#> 315     emotion   cohesion    emotion   cohesion       plan  consensus
#> 316    cohesion  consensus    discuss  synthesis       plan  consensus
#> 317        plan       plan  consensus    monitor       plan       plan
#> 318   consensus    discuss   cohesion    emotion   cohesion       plan
#> 319        plan       plan    emotion  consensus coregulate       <NA>
#> 320        plan  consensus coregulate    emotion coregulate       plan
#> 321   consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 322        plan    monitor       <NA>       <NA>       <NA>       <NA>
#> 323  coregulate  consensus       plan       plan    monitor  consensus
#> 324    cohesion    emotion  consensus coregulate    emotion  consensus
#> 325        plan       plan    emotion  consensus    monitor       plan
#> 326        plan       plan    discuss       <NA>       <NA>       <NA>
#> 327   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 328     discuss  consensus    monitor   cohesion  consensus    discuss
#> 329        plan coregulate    monitor   cohesion       plan  consensus
#> 330    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 331     emotion    discuss coregulate    discuss    discuss  consensus
#> 332     emotion       plan    emotion       plan    monitor       plan
#> 333        plan    emotion  consensus  consensus    monitor       plan
#> 334        plan  consensus  consensus coregulate       plan  consensus
#> 335   consensus       plan       plan  consensus    emotion   cohesion
#> 336    cohesion    emotion  consensus    monitor       plan  consensus
#> 337  coregulate  consensus coregulate    monitor    discuss  consensus
#> 338   consensus    discuss    emotion       plan  consensus coregulate
#> 339     emotion  consensus    discuss    emotion  consensus    discuss
#> 340        plan    emotion    discuss  consensus  consensus    discuss
#> 341   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 342   consensus       plan  consensus coregulate       plan    emotion
#> 343     emotion  consensus  consensus    emotion  consensus    emotion
#> 344     discuss   cohesion    discuss    discuss  consensus coregulate
#> 345        plan       plan  consensus       plan    discuss    discuss
#> 346        plan       plan  consensus coregulate    discuss  consensus
#> 347        plan       plan       plan    monitor       plan    emotion
#> 348   consensus    emotion   cohesion  consensus       plan       plan
#> 349        plan       plan       <NA>       <NA>       <NA>       <NA>
#> 350     emotion  consensus    discuss  consensus       <NA>       <NA>
#> 351        plan       plan   cohesion  consensus coregulate    monitor
#> 352        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 353     emotion   cohesion  consensus  consensus       plan       plan
#> 354   consensus    discuss coregulate    discuss coregulate       plan
#> 355   synthesis  consensus   cohesion  consensus coregulate coregulate
#> 356   consensus       plan  consensus       plan       <NA>       <NA>
#> 357        plan  consensus  consensus coregulate    discuss    monitor
#> 358        plan       plan    emotion  consensus  consensus    monitor
#> 359     emotion   cohesion  consensus    discuss      adapt       plan
#> 360   consensus       plan  consensus    discuss  consensus    discuss
#> 361  coregulate    emotion    emotion   cohesion    monitor  consensus
#> 362   consensus    discuss  consensus    monitor    discuss  consensus
#> 363     discuss  consensus       plan   cohesion  consensus coregulate
#> 364        plan   cohesion  consensus    discuss    discuss  synthesis
#> 365   consensus coregulate    emotion  consensus       plan   cohesion
#> 366        plan       plan    emotion  consensus       plan       plan
#> 367   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 368   consensus       plan  consensus    emotion  consensus  consensus
#> 369     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 370     monitor    discuss  consensus       plan  consensus    discuss
#> 371        plan    discuss  consensus       plan    discuss coregulate
#> 372       adapt  consensus       <NA>       <NA>       <NA>       <NA>
#> 373   consensus       plan       plan       plan    monitor    discuss
#> 374    cohesion  consensus    discuss  consensus  consensus       plan
#> 375   consensus    emotion  consensus       plan       plan    discuss
#> 376   consensus       plan       plan    emotion    emotion       <NA>
#> 377    cohesion  consensus    discuss  consensus       plan  consensus
#> 378     discuss coregulate   cohesion  consensus       plan       plan
#> 379     discuss    discuss  synthesis      adapt  consensus    monitor
#> 380   consensus coregulate    discuss       <NA>       <NA>       <NA>
#> 381       adapt  consensus  consensus    monitor coregulate    emotion
#> 382   consensus    discuss  synthesis  consensus    discuss    emotion
#> 383     emotion   cohesion  consensus  consensus       <NA>       <NA>
#> 384    cohesion  consensus  consensus  consensus       plan       plan
#> 385     emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 386     discuss  consensus       plan    emotion   cohesion  consensus
#> 387        plan  consensus       plan coregulate    discuss  consensus
#> 388     discuss    discuss  consensus    discuss    discuss  consensus
#> 389  coregulate  synthesis       <NA>       <NA>       <NA>       <NA>
#> 390     emotion    discuss  synthesis  consensus       plan       plan
#> 391   consensus    discuss  consensus    discuss  consensus  consensus
#> 392        plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 393   consensus       plan  consensus    monitor  consensus    discuss
#> 394     monitor    discuss       <NA>       <NA>       <NA>       <NA>
#> 395        plan    emotion  consensus    discuss  consensus    monitor
#> 396  coregulate    monitor    discuss   cohesion  consensus    emotion
#> 397     monitor      adapt  consensus       plan    monitor  consensus
#> 398    cohesion  consensus       plan       plan  consensus       plan
#> 399        plan  consensus       plan       plan       plan    monitor
#> 400     emotion   cohesion  consensus    discuss  consensus coregulate
#> 401   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 402   consensus coregulate       plan  consensus       plan   cohesion
#> 403   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 404     discuss  consensus coregulate      adapt  consensus coregulate
#> 405        plan       plan  consensus    monitor coregulate    emotion
#> 406     discuss  synthesis       <NA>       <NA>       <NA>       <NA>
#> 407     discuss  consensus       plan       plan       plan  consensus
#> 408    cohesion       plan    emotion   cohesion    emotion  consensus
#> 409   consensus       plan       plan       plan  consensus       <NA>
#> 410        plan       plan       plan       plan       plan       plan
#> 411     discuss  synthesis coregulate    emotion       plan       plan
#> 412     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 413   consensus coregulate  consensus       plan  consensus       plan
#> 414     discuss  consensus       plan       plan    emotion    emotion
#> 415     discuss      adapt  consensus    monitor  consensus    discuss
#> 416     emotion       plan       plan    discuss    emotion  consensus
#> 417        plan  consensus coregulate  consensus  consensus    discuss
#> 418   consensus       plan   cohesion       plan    emotion    discuss
#> 419        plan  consensus coregulate    discuss  consensus  consensus
#> 420   consensus       plan       plan    emotion  consensus  consensus
#> 421     emotion  consensus coregulate    discuss  consensus       plan
#> 422   consensus    discuss  consensus       plan    monitor  consensus
#> 423   consensus    discuss    discuss  consensus       plan       plan
#> 424  coregulate    emotion   cohesion   cohesion  consensus    emotion
#> 425     discuss  consensus    discuss    emotion       <NA>       <NA>
#> 426     discuss  consensus       plan       plan       plan  consensus
#> 427  coregulate    discuss    discuss  consensus  synthesis  consensus
#> 428     monitor       plan    emotion   cohesion coregulate       plan
#> 429        plan   cohesion  consensus coregulate       plan    emotion
#> 430        plan       plan    emotion   cohesion  consensus    discuss
#> 431     discuss    emotion coregulate    discuss    discuss coregulate
#> 432     discuss  consensus    discuss  consensus       plan  consensus
#> 433  coregulate  consensus       plan    emotion   cohesion       plan
#> 434   consensus  consensus    emotion  consensus    discuss  consensus
#> 435    cohesion  consensus       plan       plan       plan    discuss
#> 436  coregulate       plan       plan    emotion  consensus coregulate
#> 437   consensus       plan       plan    monitor    emotion  consensus
#> 438   consensus       plan    monitor    discuss  synthesis   cohesion
#> 439    cohesion  consensus    discuss    discuss    discuss coregulate
#> 440        plan       plan       plan       plan       plan       plan
#> 441        plan    discuss       <NA>       <NA>       <NA>       <NA>
#> 442        plan    emotion  consensus    discuss       <NA>       <NA>
#> 443     discuss    emotion       plan  consensus  consensus coregulate
#> 444    cohesion  consensus    discuss  consensus    discuss       plan
#> 445    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 446   consensus    emotion   cohesion coregulate    emotion  consensus
#> 447        plan  consensus    emotion  consensus coregulate    emotion
#> 448     discuss   cohesion   cohesion  consensus       plan  consensus
#> 449        plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 450     emotion  consensus       plan  consensus coregulate    discuss
#> 451     emotion       plan  consensus    monitor    discuss   cohesion
#> 452   consensus coregulate       plan       plan    emotion   cohesion
#> 453     discuss    discuss    emotion       <NA>       <NA>       <NA>
#> 454   consensus       plan    monitor    discuss    emotion   cohesion
#> 455        plan  consensus       plan    emotion   cohesion  consensus
#> 456   consensus    emotion  consensus    discuss    discuss coregulate
#> 457     discuss       plan    monitor    discuss   cohesion  consensus
#> 458  coregulate    emotion    discuss  synthesis  consensus  consensus
#> 459     discuss  consensus       plan    emotion   cohesion       plan
#> 460     monitor    discuss  consensus       plan       plan    monitor
#> 461   consensus coregulate    emotion  consensus    discuss      adapt
#> 462        plan  consensus       plan       plan       plan       plan
#> 463  coregulate       plan    discuss    discuss  consensus       plan
#> 464     monitor  consensus coregulate    discuss  consensus  consensus
#> 465     emotion  consensus    monitor coregulate       <NA>       <NA>
#> 466        plan       plan       plan  consensus      adapt       <NA>
#> 467     discuss    discuss   cohesion  consensus  consensus    discuss
#> 468        plan  consensus    discuss    emotion    discuss  synthesis
#> 469        plan       plan  consensus       plan       plan    discuss
#> 470     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 471   consensus       plan  consensus  consensus   cohesion  consensus
#> 472  coregulate    emotion   cohesion    emotion   cohesion  consensus
#> 473     emotion  consensus coregulate       plan       <NA>       <NA>
#> 474     monitor coregulate       plan    emotion    discuss   cohesion
#> 475     emotion    discuss    discuss  consensus       plan  consensus
#> 476   consensus       plan       plan       plan    emotion       plan
#> 477     monitor    discuss    discuss  consensus       plan  consensus
#> 478   consensus    discuss    discuss  synthesis    discuss    discuss
#> 479   consensus    discuss    emotion       plan   cohesion  consensus
#> 480   consensus  synthesis    emotion       plan   cohesion  consensus
#> 481  coregulate   cohesion       <NA>       <NA>       <NA>       <NA>
#> 482     monitor    discuss coregulate    emotion  consensus       plan
#> 483   consensus       plan  consensus coregulate    emotion    discuss
#> 484   consensus       plan  consensus    discuss    discuss  consensus
#> 485    cohesion       plan       plan  consensus    discuss    emotion
#> 486        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 487     emotion   cohesion  consensus    discuss  consensus    discuss
#> 488     discuss  consensus    discuss   cohesion       <NA>       <NA>
#> 489   consensus coregulate    discuss       <NA>       <NA>       <NA>
#> 490     discuss   cohesion  consensus coregulate    discuss  synthesis
#> 491        plan    discuss       <NA>       <NA>       <NA>       <NA>
#> 492        plan       plan    emotion    emotion  consensus       plan
#> 493   consensus    emotion  consensus  consensus       plan  consensus
#> 494        plan  consensus    discuss  consensus       plan  consensus
#> 495   consensus    discuss      adapt    emotion  consensus    emotion
#> 496        plan    emotion       plan       plan  consensus       <NA>
#> 497     discuss    emotion   cohesion  consensus    emotion    monitor
#> 498    cohesion  consensus       plan    emotion  consensus    discuss
#> 499   consensus    discuss  synthesis  consensus       plan       plan
#> 500     emotion coregulate    emotion   cohesion  consensus       plan
#> 501    cohesion  consensus       plan       plan       plan       plan
#> 502    cohesion  consensus       plan  consensus    discuss  consensus
#> 503   consensus coregulate       plan    emotion   cohesion    emotion
#> 504        plan       plan       plan  consensus       plan    emotion
#> 505  coregulate    emotion       <NA>       <NA>       <NA>       <NA>
#> 506     monitor       plan  consensus       plan  consensus       plan
#> 507   consensus coregulate       plan  consensus       plan       plan
#> 508   consensus       plan  consensus    discuss    emotion  consensus
#> 509   consensus coregulate   cohesion       plan  consensus coregulate
#> 510     emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 511        plan       plan    emotion   cohesion       plan       plan
#> 512     monitor    discuss  consensus       plan    emotion  consensus
#> 513     emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 514     discuss  consensus       plan       plan  consensus       plan
#> 515   consensus coregulate    discuss       <NA>       <NA>       <NA>
#> 516   synthesis    discuss  consensus       plan       plan  consensus
#> 517        plan       plan  consensus       plan    emotion  consensus
#> 518     discuss      adapt  consensus coregulate    discuss    discuss
#> 519     discuss    emotion    discuss  consensus       plan  consensus
#> 520     discuss  consensus    discuss  consensus coregulate    emotion
#> 521    cohesion  consensus    emotion   cohesion    emotion  consensus
#> 522        plan    emotion   cohesion   cohesion  consensus  synthesis
#> 523        plan  consensus       plan    emotion       plan    emotion
#> 524     monitor    emotion   cohesion coregulate       <NA>       <NA>
#> 525   synthesis  consensus       <NA>       <NA>       <NA>       <NA>
#> 526        plan  consensus       plan       <NA>       <NA>       <NA>
#> 527    cohesion  consensus    discuss   cohesion       plan  consensus
#> 528     monitor       plan       plan  consensus       plan  consensus
#> 529   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 530     discuss  consensus       plan    emotion  consensus       <NA>
#> 531        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 532     emotion   cohesion    emotion   cohesion  consensus       plan
#> 533        plan       plan    emotion  consensus    monitor       plan
#> 534        plan       plan       plan       plan    emotion  consensus
#> 535    cohesion   cohesion    monitor    discuss    emotion coregulate
#> 536   consensus  consensus   cohesion  consensus coregulate       plan
#> 537    cohesion   cohesion  consensus       plan       plan       <NA>
#> 538   synthesis    discuss    discuss   cohesion    discuss      adapt
#> 539     emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 540   consensus  consensus    discuss coregulate    discuss  synthesis
#> 541   consensus       plan  consensus coregulate    emotion    discuss
#> 542  coregulate    discuss  consensus    discuss    discuss  consensus
#> 543     emotion   cohesion coregulate  consensus       plan  consensus
#> 544        plan   cohesion    emotion       plan    discuss    monitor
#> 545     emotion   cohesion    discuss    discuss  consensus  consensus
#> 546    cohesion  consensus       plan       plan  consensus       plan
#> 547        plan       plan    emotion coregulate    emotion  consensus
#> 548  coregulate       plan       plan    discuss    discuss  consensus
#> 549   consensus       plan   cohesion  consensus    monitor  consensus
#> 550        plan coregulate       plan       plan   cohesion coregulate
#> 551        plan       plan  consensus    discuss  consensus       plan
#> 552   consensus       plan  consensus       plan       plan       plan
#> 553  coregulate    emotion   cohesion  consensus    discuss  consensus
#> 554     discuss  consensus       plan       plan  consensus    emotion
#> 555  coregulate    monitor  consensus    discuss  consensus    discuss
#> 556   synthesis  consensus coregulate    discuss  synthesis   cohesion
#> 557  coregulate    emotion    monitor  consensus coregulate    monitor
#> 558  coregulate    emotion    discuss  consensus    emotion   cohesion
#> 559        plan coregulate    emotion  consensus  consensus    discuss
#> 560     emotion  consensus       plan  consensus    emotion  consensus
#> 561  coregulate    discuss  consensus coregulate       plan    discuss
#> 562   consensus    discuss  consensus       plan  consensus       plan
#> 563  coregulate  consensus    emotion    discuss  consensus    discuss
#> 564   synthesis       <NA>       <NA>       <NA>       <NA>       <NA>
#> 565     discuss  synthesis  consensus       plan    monitor  synthesis
#> 566        plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 567   consensus  consensus       plan    monitor    discuss    emotion
#> 568    cohesion  consensus    discuss coregulate       plan  consensus
#> 569        plan  consensus       plan  consensus       plan  consensus
#> 570   consensus  consensus coregulate  consensus    monitor       plan
#> 571   consensus    discuss    discuss  consensus    monitor    discuss
#> 572     monitor    emotion  consensus    discuss  consensus    discuss
#> 573   consensus    emotion  consensus  consensus       plan    monitor
#> 574   consensus    discuss  consensus       plan    discuss    discuss
#> 575   consensus    discuss  consensus    discuss  consensus       plan
#> 576  coregulate    discuss    monitor    discuss      adapt  consensus
#> 577     discuss    emotion  consensus    discuss  consensus coregulate
#> 578    cohesion   cohesion       plan    discuss  consensus       plan
#> 579     discuss coregulate    emotion  consensus coregulate       plan
#> 580   consensus    discuss  synthesis      adapt  consensus       plan
#> 581        plan  consensus coregulate       plan  consensus       <NA>
#> 582   consensus    discuss    discuss  synthesis  consensus       plan
#> 583   consensus    discuss    discuss    emotion   cohesion       plan
#> 584    cohesion coregulate       plan    emotion   cohesion       <NA>
#> 585  coregulate       plan       <NA>       <NA>       <NA>       <NA>
#> 586    cohesion  consensus coregulate    monitor      adapt    emotion
#> 587   synthesis      adapt   cohesion  consensus    discuss    monitor
#> 588     emotion  consensus    emotion    discuss  consensus       plan
#> 589     discuss    discuss    discuss   cohesion       plan  consensus
#> 590   consensus  consensus   cohesion  consensus       plan    emotion
#> 591     discuss    discuss  consensus    discuss    discuss  consensus
#> 592        plan    emotion  consensus       plan       plan  consensus
#> 593        plan    emotion  consensus    discuss  synthesis  consensus
#> 594    cohesion coregulate   cohesion  consensus       plan  consensus
#> 595  coregulate    discuss    monitor       <NA>       <NA>       <NA>
#> 596   consensus    discuss  consensus  synthesis  consensus       plan
#> 597        plan       plan       plan    discuss   cohesion  consensus
#> 598        plan  consensus       plan  consensus       plan    emotion
#> 599     discuss  consensus   cohesion  consensus       plan       plan
#> 600   synthesis       <NA>       <NA>       <NA>       <NA>       <NA>
#> 601  coregulate    monitor       plan       plan    monitor       plan
#> 602   consensus  consensus    emotion   cohesion    emotion       plan
#> 603    cohesion    discuss  consensus       plan    monitor    discuss
#> 604     discuss  synthesis  consensus    discuss    discuss  consensus
#> 605     emotion    discuss  synthesis      adapt    emotion  consensus
#> 606     emotion  consensus       plan       plan  consensus coregulate
#> 607        plan  consensus    emotion    discuss    discuss    emotion
#> 608  coregulate    discuss    emotion    discuss  consensus    discuss
#> 609   consensus  consensus    discuss   cohesion coregulate    emotion
#> 610   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 611     discuss  consensus       plan       plan       plan    emotion
#> 612     discuss    emotion    emotion       <NA>       <NA>       <NA>
#> 613     emotion    discuss  consensus       plan    monitor       plan
#> 614  coregulate    monitor   cohesion      adapt       plan  consensus
#> 615        plan    emotion    discuss  consensus coregulate       plan
#> 616   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 617  coregulate    emotion    emotion    emotion    emotion  consensus
#> 618   consensus    discuss coregulate    emotion  consensus       plan
#> 619     emotion coregulate    emotion coregulate    emotion   cohesion
#> 620        plan    monitor    emotion       plan  consensus    discuss
#> 621        plan  consensus       plan    emotion  consensus  consensus
#> 622    cohesion       plan       plan       plan  consensus       plan
#> 623     emotion   cohesion       plan       plan    monitor    monitor
#> 624  coregulate       plan    emotion   cohesion  consensus       plan
#> 625    cohesion  consensus  consensus    discuss    emotion  consensus
#> 626   consensus       plan  consensus      adapt  consensus       <NA>
#> 627        plan    emotion   cohesion       <NA>       <NA>       <NA>
#> 628  coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 629     discuss  synthesis  consensus    emotion    emotion       <NA>
#> 630        plan  consensus       plan    emotion  consensus       <NA>
#> 631   consensus       plan       plan    emotion coregulate       plan
#> 632  coregulate    emotion  consensus  consensus    monitor    discuss
#> 633   consensus    discuss    discuss    discuss    discuss  consensus
#> 634   consensus       plan       plan  consensus    discuss  consensus
#> 635     discuss  consensus       plan    emotion   cohesion  consensus
#> 636     monitor    discuss  synthesis  consensus    emotion  consensus
#> 637     monitor  consensus       <NA>       <NA>       <NA>       <NA>
#> 638     discuss  consensus       plan       plan       plan       plan
#> 639        plan       plan    discuss    emotion       <NA>       <NA>
#> 640   consensus  consensus    discuss  consensus       plan  consensus
#> 641        plan       plan       plan       plan       plan  consensus
#> 642     emotion    emotion   cohesion  consensus       plan   cohesion
#> 643     discuss coregulate  consensus coregulate  consensus       plan
#> 644     emotion  consensus coregulate  consensus       plan    emotion
#> 645     emotion  consensus       plan       plan       plan  consensus
#> 646   consensus coregulate    discuss  consensus       plan   cohesion
#> 647     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 648        plan  consensus    emotion  consensus    monitor  consensus
#> 649    cohesion  consensus    discuss   cohesion   cohesion coregulate
#> 650     discuss    emotion  consensus       plan    emotion  consensus
#> 651    cohesion coregulate      adapt  consensus   cohesion  consensus
#> 652        plan       plan  synthesis      adapt       <NA>       <NA>
#> 653   consensus    discuss coregulate       plan    emotion   cohesion
#> 654       adapt   cohesion       plan  consensus coregulate    emotion
#> 655   consensus       plan  consensus       plan  consensus  consensus
#> 656     emotion  consensus       plan       plan    discuss   cohesion
#> 657     emotion   cohesion  consensus    emotion   cohesion coregulate
#> 658    cohesion    emotion   cohesion    discuss  consensus coregulate
#> 659   consensus       plan       plan    emotion   cohesion  consensus
#> 660        plan       plan       plan  consensus       plan    discuss
#> 661   consensus       plan  consensus    discuss       <NA>       <NA>
#> 662   consensus       plan       plan    emotion  consensus coregulate
#> 663     discuss    discuss       <NA>       <NA>       <NA>       <NA>
#> 664   consensus    emotion    discuss  synthesis      adapt  consensus
#> 665   consensus coregulate       plan    monitor    discuss coregulate
#> 666   synthesis  consensus       plan  consensus  consensus       plan
#> 667        plan       plan       <NA>       <NA>       <NA>       <NA>
#> 668     discuss  consensus    discuss  consensus       plan    emotion
#> 669    cohesion       plan    emotion  consensus       plan  consensus
#> 670        plan       plan    discuss  consensus       plan       plan
#> 671   consensus    discuss    monitor    discuss    emotion    emotion
#> 672    cohesion  consensus coregulate    monitor    discuss       <NA>
#> 673     emotion  consensus    discuss    discuss  consensus    discuss
#> 674    cohesion coregulate      adapt   cohesion coregulate    discuss
#> 675       adapt    emotion  consensus       plan       plan  consensus
#> 676   consensus       plan  consensus coregulate    discuss  consensus
#> 677        plan    emotion       plan       plan    emotion       plan
#> 678     discuss    discuss  consensus    discuss  consensus       plan
#> 679   consensus coregulate      adapt  consensus       plan    monitor
#> 680   consensus coregulate   cohesion    emotion   cohesion       plan
#> 681     discuss coregulate       plan    emotion  consensus  consensus
#> 682     discuss  consensus    discuss  synthesis      adapt  consensus
#> 683     emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 684        plan       plan       plan  consensus       plan       plan
#> 685        plan       plan    monitor       plan    emotion   cohesion
#> 686        plan  consensus coregulate    discuss  consensus coregulate
#> 687   consensus coregulate   cohesion  consensus    discuss  consensus
#> 688        plan  consensus  consensus    emotion  consensus       plan
#> 689     discuss  consensus coregulate  synthesis    emotion   cohesion
#> 690    cohesion    emotion   cohesion       plan   cohesion  consensus
#> 691    cohesion    emotion    discuss    discuss  consensus coregulate
#> 692     discuss    discuss  consensus    emotion  consensus       plan
#> 693   consensus    discuss coregulate    emotion   cohesion  consensus
#> 694     discuss   cohesion       plan    emotion   cohesion  consensus
#> 695     discuss  consensus    emotion  consensus       plan  consensus
#> 696   synthesis  consensus coregulate       plan  consensus       plan
#> 697     discuss    emotion  consensus  consensus    discuss  consensus
#> 698   synthesis       plan    emotion   cohesion    emotion  consensus
#> 699     emotion  consensus       plan  consensus       <NA>       <NA>
#> 700  coregulate       plan       plan  consensus    discuss    discuss
#> 701    cohesion      adapt       <NA>       <NA>       <NA>       <NA>
#> 702    cohesion  consensus       plan  consensus       plan       plan
#> 703     monitor    discuss    emotion    emotion   cohesion  consensus
#> 704    cohesion  consensus       <NA>       <NA>       <NA>       <NA>
#> 705   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 706        plan       plan       plan    monitor       plan       plan
#> 707   consensus    emotion       plan       plan  consensus       plan
#> 708        plan    emotion   cohesion  consensus       plan  consensus
#> 709   consensus    emotion  consensus coregulate    discuss  consensus
#> 710   synthesis       plan    emotion    discuss  consensus coregulate
#> 711        plan  consensus    discuss  synthesis      adapt  consensus
#> 712   consensus  consensus    discuss  consensus    discuss  consensus
#> 713   consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 714   consensus       plan       plan    emotion  consensus coregulate
#> 715    cohesion       plan   cohesion    emotion  consensus       plan
#> 716  coregulate       plan       plan  consensus coregulate    discuss
#> 717     discuss    discuss  consensus coregulate    monitor    discuss
#> 718  coregulate       plan  consensus    discuss       <NA>       <NA>
#> 719  coregulate    discuss       <NA>       <NA>       <NA>       <NA>
#> 720     discuss      adapt  consensus    discuss   cohesion       plan
#> 721     discuss   cohesion  consensus       plan  consensus       plan
#> 722    cohesion    discuss  synthesis  consensus  synthesis    discuss
#> 723   consensus       plan       plan  consensus    discuss  consensus
#> 724   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 725   consensus    discuss   cohesion  consensus       plan       plan
#> 726   consensus    discuss    emotion  consensus       plan  consensus
#> 727   consensus    discuss  synthesis      adapt       <NA>       <NA>
#> 728     discuss  consensus    discuss    discuss    emotion coregulate
#> 729        plan  consensus coregulate       plan       plan   cohesion
#> 730   consensus    discuss coregulate       <NA>       <NA>       <NA>
#> 731     discuss  synthesis  consensus   cohesion  consensus       plan
#> 732     emotion    emotion   cohesion  consensus       plan  consensus
#> 733     discuss  consensus       plan       plan  consensus       plan
#> 734        plan coregulate       plan       plan  consensus coregulate
#> 735        plan       plan  consensus       plan  consensus       plan
#> 736   consensus       plan   cohesion  consensus       plan    emotion
#> 737  coregulate    emotion  consensus       plan  consensus coregulate
#> 738     discuss coregulate    discuss    discuss    discuss  synthesis
#> 739     discuss   cohesion  consensus       plan       <NA>       <NA>
#> 740     discuss  consensus    discuss    discuss    discuss  consensus
#> 741     discuss  consensus      adapt  consensus    discuss  consensus
#> 742     emotion   cohesion  consensus coregulate  consensus    discuss
#> 743   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 744   consensus coregulate       plan    emotion  consensus    monitor
#> 745  coregulate       plan  consensus   cohesion       plan  consensus
#> 746        plan  consensus    discuss    emotion  consensus    discuss
#> 747    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 748        plan       plan       plan  consensus    monitor   cohesion
#> 749   consensus coregulate  consensus       plan    emotion  consensus
#> 750   consensus       plan    discuss  consensus coregulate  synthesis
#> 751        plan  consensus    emotion   cohesion  consensus    discuss
#> 752     emotion       plan       <NA>       <NA>       <NA>       <NA>
#> 753  coregulate   cohesion  consensus       plan    emotion  consensus
#> 754   consensus coregulate  consensus    discuss  consensus       <NA>
#> 755  coregulate  consensus    emotion   cohesion    emotion    discuss
#> 756     emotion coregulate    emotion   cohesion  consensus    monitor
#> 757    cohesion       plan    discuss  consensus       plan       plan
#> 758        plan       plan  consensus coregulate    monitor       plan
#> 759   consensus coregulate    discuss    discuss    discuss  consensus
#> 760        plan       plan       plan    monitor    discuss    emotion
#> 761     monitor       <NA>       <NA>       <NA>       <NA>       <NA>
#> 762        plan  consensus       plan       plan       plan       <NA>
#> 763   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 764        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 765        plan       plan   cohesion       <NA>       <NA>       <NA>
#> 766   consensus       plan  consensus  consensus       plan  consensus
#> 767     discuss  consensus    discuss    emotion    discuss    emotion
#> 768        plan    emotion       plan       plan    emotion coregulate
#> 769    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 770        plan    monitor    discuss coregulate   cohesion  consensus
#> 771   consensus    discuss  synthesis       plan  consensus coregulate
#> 772     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 773     emotion    discuss  consensus       plan    monitor    discuss
#> 774   consensus  consensus    discuss  consensus       plan  consensus
#> 775        plan  consensus coregulate       plan    monitor   cohesion
#> 776   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 777    cohesion  consensus    emotion   cohesion    monitor    emotion
#> 778       adapt   cohesion   cohesion  consensus    monitor    monitor
#> 779     discuss  consensus       plan  consensus       plan       plan
#> 780   consensus  consensus       <NA>       <NA>       <NA>       <NA>
#> 781        plan       plan  consensus       plan    monitor    discuss
#> 782     discuss  synthesis  consensus    discuss   cohesion coregulate
#> 783     emotion       plan    emotion  consensus coregulate       plan
#> 784     discuss  consensus  consensus    discuss    discuss  consensus
#> 785   consensus coregulate       plan  consensus    discuss   cohesion
#> 786     discuss    discuss       <NA>       <NA>       <NA>       <NA>
#> 787   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 788   consensus coregulate    monitor    emotion  consensus       plan
#> 789     discuss    discuss coregulate  synthesis  consensus       plan
#> 790        plan    monitor    discuss  consensus       <NA>       <NA>
#> 791   synthesis      adapt    monitor coregulate    emotion  synthesis
#> 792   consensus coregulate  consensus coregulate  consensus coregulate
#> 793   consensus       plan       plan   cohesion  consensus   cohesion
#> 794   consensus       plan       plan       plan    discuss  consensus
#> 795   consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 796        plan    monitor    emotion       <NA>       <NA>       <NA>
#> 797    cohesion  consensus       plan       plan       plan    emotion
#> 798     discuss  consensus    emotion       <NA>       <NA>       <NA>
#> 799        plan    monitor    discuss  consensus    discuss    discuss
#> 800     discuss    emotion   cohesion  consensus    discuss   cohesion
#> 801        plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 802   consensus    discuss      adapt   cohesion  consensus    discuss
#> 803    cohesion  consensus       plan    monitor    discuss  consensus
#> 804        plan       plan       plan  consensus   cohesion coregulate
#> 805     monitor       plan    emotion  consensus       plan       plan
#> 806    cohesion  consensus       plan coregulate       plan  consensus
#> 807        plan  consensus       plan    discuss  consensus coregulate
#> 808        plan       plan  synthesis  consensus       plan    emotion
#> 809   consensus    discuss coregulate    discuss    discuss  consensus
#> 810    cohesion    monitor    discuss   cohesion    emotion       <NA>
#> 811  coregulate    emotion    discuss    discuss  consensus coregulate
#> 812   consensus       plan  consensus    discuss    discuss    emotion
#> 813     monitor    discuss  consensus       plan  consensus       plan
#> 814   consensus       plan       plan    emotion   cohesion   cohesion
#> 815   consensus    discuss    emotion  consensus       plan  consensus
#> 816   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 817        plan       plan       plan  consensus       plan    discuss
#> 818   consensus       plan       plan    emotion    discuss  consensus
#> 819   consensus  consensus       plan  consensus       plan       plan
#> 820        plan    emotion   cohesion    emotion  consensus    discuss
#> 821        plan    emotion    discuss  consensus    emotion   cohesion
#> 822        plan       plan  synthesis  consensus       plan       plan
#> 823  coregulate    discuss coregulate    emotion   cohesion       plan
#> 824   consensus coregulate       plan    discuss  synthesis      adapt
#> 825   consensus    discuss   cohesion  consensus       plan    emotion
#> 826       adapt  consensus       plan       plan    monitor       plan
#> 827   consensus       plan       plan       plan   cohesion  consensus
#> 828   consensus       plan    monitor    discuss  consensus    discuss
#> 829   consensus  consensus    discuss    emotion  consensus    discuss
#> 830    cohesion   cohesion       plan    emotion    monitor    discuss
#> 831  coregulate    discuss      adapt   cohesion  consensus coregulate
#> 832        plan       plan    monitor  consensus  consensus coregulate
#> 833     emotion   cohesion       plan  consensus    emotion   cohesion
#> 834        plan    emotion    emotion    discuss    discuss   cohesion
#> 835     emotion    discuss   cohesion  consensus    discuss    discuss
#> 836   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 837   consensus       plan       plan      adapt   cohesion       plan
#> 838    cohesion coregulate    emotion       <NA>       <NA>       <NA>
#> 839   consensus  synthesis       <NA>       <NA>       <NA>       <NA>
#> 840     emotion  consensus    discuss  consensus    monitor    discuss
#> 841        plan  consensus    emotion    emotion  consensus coregulate
#> 842        plan    discuss coregulate       <NA>       <NA>       <NA>
#> 843   consensus    discuss    discuss  consensus    discuss  consensus
#> 844   consensus    monitor       <NA>       <NA>       <NA>       <NA>
#> 845     discuss    emotion  consensus       <NA>       <NA>       <NA>
#> 846     discuss    discuss  consensus coregulate       plan    monitor
#> 847        plan       plan       plan  consensus    discuss  consensus
#> 848     emotion   cohesion  consensus    discuss    emotion  consensus
#> 849  coregulate       plan  consensus    discuss    emotion    discuss
#> 850        plan  consensus    discuss    discuss  synthesis  consensus
#> 851   synthesis       plan    emotion   cohesion   cohesion coregulate
#> 852   consensus  consensus coregulate    emotion    discuss  consensus
#> 853     monitor    discuss  consensus    monitor coregulate  synthesis
#> 854  coregulate       plan       <NA>       <NA>       <NA>       <NA>
#> 855    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 856   synthesis coregulate    emotion   cohesion    emotion   cohesion
#> 857  coregulate  consensus  consensus    discuss    emotion   cohesion
#> 858     emotion   cohesion   cohesion       plan       plan  consensus
#> 859        plan    discuss  consensus       <NA>       <NA>       <NA>
#> 860   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 861   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 862  coregulate coregulate       plan    discuss coregulate    discuss
#> 863        plan       plan  consensus    discuss  consensus  consensus
#> 864   consensus       plan  consensus       plan       plan  consensus
#> 865     discuss  consensus       plan    emotion  consensus    emotion
#> 866     discuss  consensus    discuss  consensus    discuss       <NA>
#> 867   consensus   cohesion    emotion  consensus    discuss      adapt
#> 868  coregulate       plan       plan    emotion   cohesion  consensus
#> 869     discuss  synthesis  consensus    emotion   cohesion  consensus
#> 870     discuss  consensus    discuss  consensus    discuss    emotion
#> 871     discuss  consensus    discuss  consensus       plan  consensus
#> 872        plan  consensus coregulate       <NA>       <NA>       <NA>
#> 873  coregulate    discuss  synthesis  consensus coregulate       plan
#> 874        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 875   consensus    discuss  consensus       <NA>       <NA>       <NA>
#> 876       adapt  consensus coregulate    emotion  consensus       plan
#> 877     emotion   cohesion    emotion       plan       plan    emotion
#> 878     discuss  consensus    discuss    emotion   cohesion    discuss
#> 879   consensus    monitor    discuss  consensus    discuss  consensus
#> 880   consensus  consensus coregulate       plan       <NA>       <NA>
#> 881   consensus    discuss  consensus       plan       plan       plan
#> 882    cohesion  consensus       plan       plan    emotion  consensus
#> 883    cohesion   cohesion    emotion  consensus       plan  consensus
#> 884        plan  consensus       plan       plan    monitor coregulate
#> 885     discuss  consensus       plan    emotion  consensus  consensus
#> 886   consensus       plan       plan    discuss  consensus    discuss
#> 887  coregulate       plan    discuss  consensus coregulate       plan
#> 888     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 889     discuss  consensus       plan       plan    emotion   cohesion
#> 890        plan   cohesion       <NA>       <NA>       <NA>       <NA>
#> 891     discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 892     discuss    discuss  consensus    discuss  consensus       plan
#> 893     discuss  consensus    discuss  consensus    discuss   cohesion
#> 894        plan    discuss    discuss    discuss coregulate       plan
#> 895     emotion   cohesion coregulate    emotion    discuss  synthesis
#> 896     discuss  synthesis  consensus   cohesion  consensus       <NA>
#> 897     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 898     discuss  consensus    discuss      adapt  consensus       plan
#> 899     discuss    discuss  synthesis    emotion       <NA>       <NA>
#> 900  coregulate       plan  synthesis      adapt  consensus       plan
#> 901     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 902     emotion   cohesion       plan  consensus       plan       plan
#> 903   consensus coregulate    emotion  consensus coregulate    emotion
#> 904   consensus coregulate  consensus    emotion       <NA>       <NA>
#> 905   consensus   cohesion  consensus  consensus    emotion       plan
#> 906        plan       plan    discuss  consensus       plan       plan
#> 907        plan    monitor    emotion    monitor       <NA>       <NA>
#> 908   consensus coregulate    discuss   cohesion       <NA>       <NA>
#> 909   consensus    discuss  consensus    discuss    discuss       plan
#> 910   consensus       plan    emotion  consensus    discuss  consensus
#> 911   synthesis  consensus    emotion   cohesion       plan       plan
#> 912   consensus coregulate      adapt  consensus    discuss  consensus
#> 913     monitor    discuss  synthesis  consensus       <NA>       <NA>
#> 914   consensus coregulate       plan  consensus coregulate       plan
#> 915     discuss coregulate       <NA>       <NA>       <NA>       <NA>
#> 916     discuss  synthesis  consensus coregulate coregulate    monitor
#> 917   consensus  consensus coregulate   cohesion  consensus coregulate
#> 918   consensus    emotion    monitor    discuss       <NA>       <NA>
#> 919   consensus       plan       plan       plan       plan    emotion
#> 920        plan    emotion  consensus coregulate       plan       plan
#> 921     emotion   cohesion coregulate    discuss   cohesion  consensus
#> 922     monitor    emotion   cohesion  consensus    discuss    discuss
#> 923    cohesion  consensus    discuss    discuss  consensus    discuss
#> 924  coregulate coregulate    monitor       <NA>       <NA>       <NA>
#> 925   consensus       plan       plan       plan    emotion       plan
#> 926     discuss  consensus    discuss coregulate       plan       plan
#> 927     emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 928        plan    monitor    emotion  consensus  consensus       plan
#> 929     emotion    discuss  consensus       plan    monitor       plan
#> 930    cohesion  consensus    discuss   cohesion       plan       plan
#> 931    cohesion  consensus       plan    emotion  consensus  consensus
#> 932        plan  consensus    emotion    emotion  consensus  consensus
#> 933    cohesion  consensus    discuss    discuss  consensus       plan
#> 934        plan  consensus    discuss  consensus    discuss       plan
#> 935     emotion    discuss  consensus       plan  consensus   cohesion
#> 936     discuss  consensus  consensus  consensus       plan coregulate
#> 937  coregulate  synthesis    discuss  consensus       <NA>       <NA>
#> 938    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 939     emotion  consensus coregulate    discuss  consensus coregulate
#> 940   consensus       plan   cohesion  consensus       plan  consensus
#> 941  coregulate       plan  consensus  synthesis    emotion  consensus
#> 942   consensus  consensus    emotion    discuss       <NA>       <NA>
#> 943  coregulate  consensus    discuss  consensus       plan    emotion
#> 944   consensus    discuss  consensus       plan  consensus  consensus
#> 945     discuss    emotion  consensus       plan   cohesion  consensus
#> 946   consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 947        plan  consensus  consensus  consensus       plan    emotion
#> 948     monitor    discuss  synthesis  consensus    discuss coregulate
#> 949        plan  consensus       plan    emotion   cohesion  consensus
#> 950   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 951     emotion  consensus       plan   cohesion  consensus  consensus
#> 952        plan       plan       plan    emotion  consensus       plan
#> 953    cohesion    emotion   cohesion  consensus coregulate    discuss
#> 954     monitor    discuss   cohesion  consensus       plan    emotion
#> 955   consensus       plan  consensus    discuss    emotion   cohesion
#> 956   consensus    discuss  synthesis       plan    emotion   cohesion
#> 957   consensus       plan    emotion   cohesion       <NA>       <NA>
#> 958        plan       plan  consensus coregulate  synthesis       <NA>
#> 959   consensus    discuss    emotion   cohesion      adapt   cohesion
#> 960     discuss  consensus       plan    discuss  synthesis   cohesion
#> 961   consensus    emotion    emotion   cohesion  consensus       plan
#> 962   consensus    discuss    discuss  consensus coregulate    monitor
#> 963    cohesion   cohesion       plan    discuss    emotion   cohesion
#> 964   consensus  synthesis       <NA>       <NA>       <NA>       <NA>
#> 965    cohesion       plan       plan  consensus       plan    monitor
#> 966    cohesion  consensus  consensus       plan  consensus       plan
#> 967        plan    emotion  consensus       plan    monitor    discuss
#> 968        plan       plan       plan       <NA>       <NA>       <NA>
#> 969        plan    emotion  consensus  consensus    discuss  consensus
#> 970     monitor       plan    emotion  consensus    discuss    discuss
#> 971   synthesis       plan       plan coregulate       plan       plan
#> 972        plan    discuss  consensus coregulate       plan  consensus
#> 973     emotion    discuss  consensus       plan    monitor       plan
#> 974   synthesis       plan    emotion coregulate    monitor  consensus
#> 975        plan  consensus coregulate       plan    discuss coregulate
#> 976    cohesion       plan       plan       plan    emotion coregulate
#> 977    cohesion   cohesion  consensus    discuss      adapt  consensus
#> 978     discuss  consensus    discuss  consensus  consensus  consensus
#> 979    cohesion  consensus    discuss   cohesion       <NA>       <NA>
#> 980     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 981    cohesion  consensus    discuss coregulate coregulate       <NA>
#> 982   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 983     discuss  consensus       plan  consensus    discuss      adapt
#> 984     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 985  coregulate    emotion  consensus       plan    discuss coregulate
#> 986        plan  consensus       plan       plan    emotion   cohesion
#> 987        plan       plan coregulate    emotion       plan       plan
#> 988     discuss  consensus       plan   cohesion  consensus coregulate
#> 989   synthesis    discuss  consensus coregulate  consensus       plan
#> 990   consensus coregulate    monitor    discuss  consensus coregulate
#> 991     discuss  consensus       plan    emotion       <NA>       <NA>
#> 992   consensus    discuss  consensus    discuss    emotion    discuss
#> 993    cohesion  consensus  consensus coregulate       plan coregulate
#> 994     discuss  consensus  consensus    discuss      adapt  consensus
#> 995     discuss  consensus       plan  consensus coregulate       plan
#> 996     discuss   cohesion       plan       plan    discuss  synthesis
#> 997     emotion  consensus  consensus  consensus coregulate    discuss
#> 998     discuss  synthesis  consensus       plan       plan    monitor
#> 999   consensus    monitor    monitor       plan  consensus       plan
#> 1000    discuss  consensus       plan   cohesion  consensus    discuss
#> 1001  consensus       plan       plan    monitor       plan  consensus
#> 1002      adapt   cohesion  consensus    discuss    emotion   cohesion
#> 1003    emotion   cohesion  consensus coregulate coregulate       plan
#> 1004       plan  consensus       plan  consensus    discuss    discuss
#> 1005   cohesion    emotion       plan    monitor    discuss    emotion
#> 1006   cohesion  consensus       plan       plan       plan       plan
#> 1007       plan       plan       plan       <NA>       <NA>       <NA>
#> 1008    discuss  consensus    discuss      adapt  consensus  consensus
#> 1009       plan       plan  consensus    discuss  consensus    emotion
#> 1010    emotion  consensus    discuss    monitor    discuss  synthesis
#> 1011 coregulate    discuss    discuss    emotion       <NA>       <NA>
#> 1012 coregulate  consensus       plan       plan       plan    emotion
#> 1013       plan       plan       plan       plan       plan       <NA>
#> 1014    monitor    discuss  consensus    discuss  consensus       plan
#> 1015    monitor    discuss  consensus  consensus coregulate       plan
#> 1016   cohesion    emotion   cohesion  consensus    discuss    discuss
#> 1017    discuss    emotion   cohesion coregulate    monitor       plan
#> 1018    discuss      adapt  consensus       plan    emotion  consensus
#> 1019  consensus       plan    monitor       plan       plan  consensus
#> 1020 coregulate   cohesion  consensus       plan  consensus coregulate
#> 1021       plan       plan   cohesion coregulate       plan       plan
#> 1022  consensus       plan    emotion       <NA>       <NA>       <NA>
#> 1023 coregulate  consensus       plan       plan  consensus    emotion
#> 1024  synthesis    discuss      adapt  consensus  consensus       plan
#> 1025    discuss coregulate  consensus       plan       plan  consensus
#> 1026  synthesis      adapt       <NA>       <NA>       <NA>       <NA>
#> 1027  synthesis      adapt  consensus coregulate coregulate   cohesion
#> 1028    monitor    emotion  consensus       plan  consensus    monitor
#> 1029  consensus  consensus       plan       plan       plan       plan
#> 1030    discuss  synthesis coregulate       plan  consensus    emotion
#> 1031  consensus    emotion  consensus coregulate    emotion    discuss
#> 1032  consensus       plan    discuss    discuss    emotion  consensus
#> 1033  consensus coregulate    discuss    emotion   cohesion coregulate
#> 1034       plan       plan  consensus       plan  consensus       plan
#> 1035      adapt    emotion  consensus       plan  consensus    monitor
#> 1036 coregulate    monitor  consensus       plan       plan       plan
#> 1037    discuss    discuss coregulate       plan    monitor       <NA>
#> 1038    discuss    emotion       plan  consensus coregulate    discuss
#> 1039       plan  consensus coregulate    discuss   cohesion  consensus
#> 1040       plan  consensus    discuss  consensus       plan  consensus
#> 1041    monitor       plan    monitor    discuss    discuss      adapt
#> 1042   cohesion       plan  consensus coregulate       plan    emotion
#> 1043      adapt  consensus coregulate    discuss  synthesis      adapt
#> 1044       plan  consensus    monitor       plan  consensus       plan
#> 1045       plan    monitor       plan  consensus       <NA>       <NA>
#> 1046       plan       plan    emotion       plan  consensus       plan
#> 1047  consensus       plan       plan       plan       <NA>       <NA>
#> 1048    discuss      adapt   cohesion coregulate    discuss  consensus
#> 1049 coregulate  synthesis      adapt   cohesion    discuss    discuss
#> 1050  synthesis      adapt    emotion   cohesion coregulate    discuss
#> 1051 coregulate  consensus coregulate       <NA>       <NA>       <NA>
#> 1052       plan    discuss  consensus    monitor       plan   cohesion
#> 1053  consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 1054   cohesion  consensus       plan    monitor      adapt   cohesion
#> 1055 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1056  consensus    discuss  consensus  consensus coregulate    emotion
#> 1057  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1058       plan       plan    monitor       plan       plan    discuss
#> 1059       plan       plan       plan  consensus  consensus      adapt
#> 1060       plan    discuss      adapt coregulate  consensus       plan
#> 1061    discuss coregulate   cohesion  consensus    discuss  consensus
#> 1062  consensus       plan       plan       <NA>       <NA>       <NA>
#> 1063 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1064    emotion       plan  consensus    discuss  synthesis  consensus
#> 1065  consensus coregulate  consensus       plan    emotion  consensus
#> 1066       plan    emotion       plan       plan       <NA>       <NA>
#> 1067  consensus    discuss  consensus       plan       <NA>       <NA>
#> 1068       plan       plan  consensus coregulate      adapt   cohesion
#> 1069    emotion   cohesion coregulate  consensus       plan       plan
#> 1070    discuss      adapt  consensus coregulate    emotion  consensus
#> 1071       plan       plan    monitor       plan       plan  consensus
#> 1072       plan  consensus       plan  consensus coregulate    discuss
#> 1073       plan       plan    emotion  consensus coregulate coregulate
#> 1074       plan   cohesion  consensus       <NA>       <NA>       <NA>
#> 1075    emotion   cohesion       plan       plan    emotion    monitor
#> 1076    discuss    emotion       plan       plan       plan    emotion
#> 1077 coregulate    emotion   cohesion  consensus       plan    emotion
#> 1078  consensus coregulate    discuss    discuss      adapt  consensus
#> 1079 coregulate  consensus       plan       plan    emotion   cohesion
#> 1080       plan coregulate    emotion       plan       plan       plan
#> 1081      adapt  consensus coregulate  consensus coregulate    emotion
#> 1082    discuss    monitor    discuss  synthesis    emotion    emotion
#> 1083    discuss   cohesion  consensus coregulate       plan  consensus
#> 1084       plan  consensus       plan  consensus       plan       plan
#> 1085    discuss    discuss  consensus       plan    emotion       <NA>
#> 1086    emotion   cohesion       plan       plan       plan       plan
#> 1087  consensus       plan    discuss    discuss       <NA>       <NA>
#> 1088    discuss  consensus       plan       plan  consensus coregulate
#> 1089       plan    discuss    discuss    emotion   cohesion coregulate
#> 1090  consensus       plan       plan  consensus    monitor    monitor
#> 1091    emotion   cohesion  consensus       plan  consensus    discuss
#> 1092       plan    monitor    discuss    monitor  consensus    discuss
#> 1093       plan    emotion  consensus coregulate       plan       plan
#> 1094   cohesion coregulate    discuss  synthesis  consensus       plan
#> 1095 coregulate    discuss  consensus    emotion   cohesion  consensus
#> 1096  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1097   cohesion  consensus       plan    monitor       plan       <NA>
#> 1098 coregulate    emotion       <NA>       <NA>       <NA>       <NA>
#> 1099       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1100  consensus    discuss  synthesis      adapt    emotion    emotion
#> 1101       plan  consensus  consensus       plan       <NA>       <NA>
#> 1102    emotion       plan       <NA>       <NA>       <NA>       <NA>
#> 1103    discuss    monitor    discuss  synthesis   cohesion  consensus
#> 1104  consensus    monitor  consensus coregulate    discuss    discuss
#> 1105       plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 1106    discuss  consensus  consensus    discuss  synthesis      adapt
#> 1107  consensus       plan  consensus       plan  consensus    monitor
#> 1108    discuss    emotion       plan  consensus   cohesion coregulate
#> 1109    discuss    emotion    emotion coregulate  consensus       plan
#> 1110   cohesion  consensus       plan    emotion  consensus    emotion
#> 1111      adapt   cohesion  consensus coregulate    emotion    discuss
#> 1112    discuss    discuss       <NA>       <NA>       <NA>       <NA>
#> 1113    monitor  consensus    discuss    emotion  consensus  consensus
#> 1114  consensus       plan       plan       plan       plan       plan
#> 1115   cohesion coregulate       plan    discuss    emotion   cohesion
#> 1116  consensus coregulate   cohesion  consensus    emotion       plan
#> 1117 coregulate    discuss  synthesis  consensus       plan    discuss
#> 1118   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1119    discuss      adapt coregulate  synthesis    monitor    monitor
#> 1120    discuss      adapt   cohesion  consensus       plan       plan
#> 1121       plan       plan    emotion coregulate    discuss  synthesis
#> 1122    discuss      adapt   cohesion  consensus       plan  consensus
#> 1123       plan  consensus       plan       plan       plan    monitor
#> 1124    monitor    discuss  consensus       plan    monitor       plan
#> 1125  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1126       plan       plan  consensus       plan  consensus    discuss
#> 1127    discuss    emotion  consensus       plan       plan       plan
#> 1128    discuss  synthesis  consensus coregulate    emotion  consensus
#> 1129   cohesion  consensus    discuss    emotion    monitor    discuss
#> 1130    monitor       plan       plan       plan       plan    monitor
#> 1131       plan       plan  consensus    discuss    discuss coregulate
#> 1132       plan  consensus       plan  consensus       plan    discuss
#> 1133    discuss      adapt  consensus    discuss   cohesion  consensus
#> 1134       plan       plan  consensus       plan    monitor  consensus
#> 1135       plan  consensus coregulate  consensus       plan       plan
#> 1136    discuss      adapt  consensus coregulate coregulate      adapt
#> 1137       plan  consensus       plan  consensus       plan    monitor
#> 1138    discuss    discuss    discuss    discuss  consensus    discuss
#> 1139   cohesion  consensus    emotion    emotion  consensus    discuss
#> 1140       plan       plan       plan       plan       plan  consensus
#> 1141    emotion  consensus    discuss coregulate   cohesion       plan
#> 1142   cohesion coregulate  consensus    discuss  consensus    discuss
#> 1143    emotion   cohesion  consensus       plan       plan       plan
#> 1144 coregulate  consensus       plan    discuss  synthesis  consensus
#> 1145       plan    emotion    emotion   cohesion       plan       plan
#> 1146       plan       plan       plan       plan       plan    monitor
#> 1147    monitor    discuss    discuss   cohesion coregulate    emotion
#> 1148    discuss coregulate  synthesis  consensus       plan       plan
#> 1149       plan       plan       plan       plan    emotion   cohesion
#> 1150    emotion    discuss    emotion    discuss      adapt coregulate
#> 1151    discuss  consensus    discuss  consensus       plan       plan
#> 1152 coregulate coregulate       plan       plan       plan       plan
#> 1153       plan  consensus    monitor       plan  consensus       plan
#> 1154 coregulate    monitor  consensus    emotion  consensus coregulate
#> 1155    emotion   cohesion  consensus coregulate    monitor coregulate
#> 1156       plan  consensus       plan  consensus    emotion       <NA>
#> 1157    emotion   cohesion coregulate       plan  consensus       <NA>
#> 1158    discuss    emotion    emotion  consensus    discuss coregulate
#> 1159    discuss    discuss       <NA>       <NA>       <NA>       <NA>
#> 1160       plan  consensus       plan       plan  consensus    monitor
#> 1161       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1162 coregulate  consensus  consensus    discuss    discuss  synthesis
#> 1163  consensus    discuss  consensus coregulate    emotion   cohesion
#> 1164   cohesion  consensus       <NA>       <NA>       <NA>       <NA>
#> 1165    emotion   cohesion  consensus       plan    emotion   cohesion
#> 1166  consensus       plan       plan  consensus coregulate    discuss
#> 1167       plan  consensus       plan       plan  consensus coregulate
#> 1168  consensus       plan    monitor    discuss  consensus    discuss
#> 1169       plan       plan       plan   cohesion  consensus    discuss
#> 1170  consensus coregulate    emotion   cohesion coregulate    emotion
#> 1171      adapt  consensus       plan       plan  consensus    discuss
#> 1172       plan    emotion coregulate       <NA>       <NA>       <NA>
#> 1173    monitor    discuss    discuss    discuss  consensus  consensus
#> 1174   cohesion  consensus       <NA>       <NA>       <NA>       <NA>
#> 1175  consensus coregulate   cohesion  consensus    monitor  consensus
#> 1176      adapt  consensus       plan  consensus       plan  consensus
#> 1177    discuss coregulate    emotion   cohesion  consensus  consensus
#> 1178    discuss    discuss  consensus coregulate   cohesion       plan
#> 1179       plan  consensus       plan  consensus coregulate  consensus
#> 1180    emotion    discuss  synthesis    emotion  consensus       plan
#> 1181  consensus    monitor    emotion    discuss  synthesis    emotion
#> 1182    discuss    discuss      adapt  consensus       plan       <NA>
#> 1183    discuss coregulate       plan       plan       plan       <NA>
#> 1184       plan    emotion       plan coregulate       plan  consensus
#> 1185 coregulate   cohesion    emotion   cohesion    discuss    discuss
#> 1186    discuss   cohesion    emotion   cohesion coregulate    emotion
#> 1187  consensus    discuss    discuss    discuss  consensus coregulate
#> 1188    emotion   cohesion  consensus       plan    discuss       <NA>
#> 1189      adapt  consensus coregulate       plan       plan    emotion
#> 1190    discuss      adapt   cohesion  consensus  consensus       plan
#> 1191       plan    emotion    discuss    discuss   cohesion   cohesion
#> 1192       plan    discuss    discuss coregulate    discuss      adapt
#> 1193      adapt    emotion       plan   cohesion    emotion   cohesion
#> 1194 coregulate    emotion   cohesion       <NA>       <NA>       <NA>
#> 1195  consensus  consensus  consensus coregulate    discuss  synthesis
#> 1196       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1197       plan  consensus coregulate    discuss    discuss  synthesis
#> 1198       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1199  consensus       plan       plan       plan       plan  consensus
#> 1200  consensus  consensus       plan       plan coregulate       plan
#> 1201       plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 1202  consensus       plan       plan  consensus coregulate       plan
#> 1203   cohesion  consensus coregulate    discuss    discuss    monitor
#> 1204       plan    monitor    emotion  consensus    discuss  synthesis
#> 1205   cohesion  consensus       plan    emotion    monitor  consensus
#> 1206    emotion  consensus    monitor       plan  consensus  consensus
#> 1207    discuss      adapt   cohesion  consensus  consensus       plan
#> 1208  consensus coregulate    discuss      adapt       plan    emotion
#> 1209  consensus  consensus       plan coregulate   cohesion  consensus
#> 1210    emotion  consensus coregulate       plan       plan       plan
#> 1211       plan    discuss      adapt  consensus    discuss  consensus
#> 1212  consensus  consensus  consensus       plan  consensus       plan
#> 1213  consensus       plan    emotion   cohesion  consensus       plan
#> 1214    monitor    discuss       <NA>       <NA>       <NA>       <NA>
#> 1215    discuss    discuss    discuss  consensus       plan       plan
#> 1216  consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 1217  synthesis    discuss       <NA>       <NA>       <NA>       <NA>
#> 1218       plan       plan  consensus coregulate  consensus    discuss
#> 1219  consensus       plan       plan       plan       plan  consensus
#> 1220    monitor       plan    monitor   cohesion       plan  consensus
#> 1221       plan       plan       plan       plan       plan       plan
#> 1222    discuss    emotion   cohesion       plan       plan  consensus
#> 1223    discuss coregulate    discuss    emotion       plan       plan
#> 1224  consensus    discuss  synthesis  consensus       plan    monitor
#> 1225  consensus coregulate    discuss    discuss  consensus       <NA>
#> 1226   cohesion       plan       plan  consensus    emotion   cohesion
#> 1227    emotion  consensus       plan  consensus       plan       plan
#> 1228  consensus       plan       plan  consensus       plan  consensus
#> 1229  consensus       plan  consensus  synthesis      adapt coregulate
#> 1230 coregulate    discuss       <NA>       <NA>       <NA>       <NA>
#> 1231   cohesion coregulate    discuss   cohesion  consensus coregulate
#> 1232    discuss      adapt       <NA>       <NA>       <NA>       <NA>
#> 1233  synthesis    emotion    discuss    discuss  consensus    emotion
#> 1234       plan    monitor       plan  consensus    discuss       plan
#> 1235    discuss coregulate    discuss  consensus    discuss    discuss
#> 1236    emotion   cohesion    emotion    monitor   cohesion  consensus
#> 1237       plan       plan    emotion       plan       <NA>       <NA>
#> 1238       plan    emotion   cohesion  consensus    discuss    emotion
#> 1239   cohesion       plan    monitor  consensus       plan       plan
#> 1240       plan       plan       plan       plan       plan       plan
#> 1241  consensus    discuss coregulate    emotion   cohesion  consensus
#> 1242    discuss  consensus coregulate    emotion  consensus       plan
#> 1243    discuss    emotion coregulate    emotion    discuss coregulate
#> 1244    discuss  synthesis      adapt  consensus       plan    emotion
#> 1245   cohesion coregulate       plan    monitor    discuss  synthesis
#> 1246       plan       plan       plan       plan       plan    discuss
#> 1247    discuss  consensus       plan    monitor       plan    emotion
#> 1248  consensus       plan       plan    discuss    emotion  consensus
#> 1249    monitor    discuss    discuss  synthesis  consensus       plan
#> 1250       plan       plan       plan   cohesion    emotion   cohesion
#> 1251  consensus  consensus coregulate    discuss    discuss    discuss
#> 1252       plan       plan       plan  consensus    emotion   cohesion
#> 1253       plan  consensus       plan       plan       plan       plan
#> 1254  consensus    discuss    emotion    discuss  synthesis  consensus
#> 1255  consensus    monitor    discuss coregulate       plan       plan
#> 1256    discuss    discuss      adapt  consensus       plan       plan
#> 1257  synthesis      adapt    monitor    discuss    discuss    discuss
#> 1258    monitor       plan       <NA>       <NA>       <NA>       <NA>
#> 1259    discuss  synthesis       <NA>       <NA>       <NA>       <NA>
#> 1260  consensus    discuss    emotion   cohesion    discuss  consensus
#> 1261   cohesion coregulate    discuss  synthesis  consensus  consensus
#> 1262    discuss coregulate    monitor       plan    emotion    discuss
#> 1263    emotion       plan       plan       <NA>       <NA>       <NA>
#> 1264    emotion   cohesion       plan       plan    monitor       plan
#> 1265    monitor  consensus  consensus       plan       plan       plan
#> 1266    discuss  synthesis  consensus       <NA>       <NA>       <NA>
#> 1267    emotion  consensus    monitor  consensus  consensus       plan
#> 1268 coregulate    emotion   cohesion    emotion  consensus       plan
#> 1269       plan       plan       plan    monitor  consensus       plan
#> 1270    monitor    discuss  synthesis      adapt    emotion    monitor
#> 1271       plan       plan  consensus       plan       plan       plan
#> 1272       plan    emotion  consensus  consensus       plan  consensus
#> 1273       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1274  synthesis  consensus    discuss  consensus    discuss    emotion
#> 1275 coregulate       plan       plan       plan       plan  consensus
#> 1276      adapt  consensus    discuss  synthesis  consensus       plan
#> 1277       plan       plan       plan  consensus coregulate       plan
#> 1278  consensus       plan  consensus       plan       plan       plan
#> 1279       plan       plan       plan  consensus       plan       plan
#> 1280       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1281       plan    discuss    emotion  consensus    discuss    emotion
#> 1282    emotion       plan       <NA>       <NA>       <NA>       <NA>
#> 1283      adapt   cohesion    monitor    discuss    emotion    monitor
#> 1284      adapt coregulate       plan  consensus       plan    emotion
#> 1285       plan    discuss    discuss    emotion  consensus       plan
#> 1286    emotion coregulate    monitor coregulate    discuss    discuss
#> 1287      adapt    discuss  synthesis      adapt  consensus    discuss
#> 1288  synthesis      adapt   cohesion       plan       plan       plan
#> 1289    discuss  synthesis      adapt  consensus    discuss  synthesis
#> 1290  consensus       plan       plan       plan    emotion       <NA>
#> 1291       plan  consensus       plan    monitor       plan       <NA>
#> 1292    discuss coregulate       <NA>       <NA>       <NA>       <NA>
#> 1293    discuss coregulate       plan coregulate    emotion    discuss
#> 1294    monitor  consensus       plan  consensus  consensus coregulate
#> 1295  synthesis      adapt  consensus       plan    monitor  consensus
#> 1296       plan  consensus    emotion       plan  consensus       plan
#> 1297    discuss  synthesis      adapt   cohesion  consensus       plan
#> 1298       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1299    discuss coregulate    discuss      adapt   cohesion  consensus
#> 1300       plan       plan       plan    emotion       plan    emotion
#> 1301    discuss      adapt  consensus  consensus       plan       plan
#> 1302    emotion  consensus  consensus       plan       plan       plan
#> 1303    discuss    discuss    discuss      adapt  consensus       plan
#> 1304       plan  consensus       plan       plan  consensus    monitor
#> 1305    discuss  synthesis coregulate       plan       plan    emotion
#> 1306 coregulate  consensus    discuss  consensus    monitor  consensus
#> 1307  consensus       plan  consensus       plan       plan       plan
#> 1308    discuss      adapt  consensus       plan    emotion    discuss
#> 1309    discuss      adapt   cohesion  consensus       plan  consensus
#> 1310 coregulate       plan  consensus       plan  consensus    discuss
#> 1311       plan  consensus coregulate       plan    discuss    emotion
#> 1312       plan       plan       plan  consensus       <NA>       <NA>
#> 1313       plan       plan       plan       plan       <NA>       <NA>
#> 1314    discuss  synthesis       <NA>       <NA>       <NA>       <NA>
#> 1315  consensus       plan    discuss  consensus       plan       plan
#> 1316       plan  consensus       plan  consensus  consensus       <NA>
#> 1317       plan  consensus       plan       plan  consensus  consensus
#> 1318  consensus coregulate       plan    discuss      adapt   cohesion
#> 1319    discuss  synthesis      adapt    emotion   cohesion  consensus
#> 1320    discuss  consensus coregulate    emotion   cohesion  consensus
#> 1321       plan  consensus    emotion   cohesion       plan    emotion
#> 1322  consensus       plan       plan    emotion   cohesion  consensus
#> 1323  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1324    discuss      adapt   cohesion coregulate  consensus       plan
#> 1325       plan       plan   cohesion       plan       plan       plan
#> 1326    discuss    discuss    discuss  consensus       plan       <NA>
#> 1327       plan    discuss  synthesis    emotion       plan       plan
#> 1328  consensus    discuss    emotion  consensus  consensus    emotion
#> 1329    monitor    emotion    emotion coregulate    discuss    emotion
#> 1330  consensus       plan  consensus      adapt  consensus       plan
#> 1331    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1332    discuss    emotion  consensus coregulate       plan       plan
#> 1333  consensus  consensus coregulate  consensus coregulate    emotion
#> 1334   cohesion  consensus       plan  consensus    emotion    discuss
#> 1335 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1336  consensus coregulate       plan    monitor coregulate    discuss
#> 1337       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1338       plan    emotion  consensus  consensus       plan    monitor
#> 1339  consensus coregulate coregulate    emotion  consensus coregulate
#> 1340       plan       plan    emotion    emotion  consensus       plan
#> 1341  consensus  consensus       plan  consensus    discuss  consensus
#> 1342   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1343    monitor    discuss  consensus  consensus    emotion   cohesion
#> 1344  synthesis      adapt  consensus coregulate       plan    monitor
#> 1345  consensus       plan    monitor  consensus    emotion       plan
#> 1346       plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 1347  consensus coregulate    emotion       plan    emotion coregulate
#> 1348 coregulate       plan       plan   cohesion  consensus coregulate
#> 1349    emotion   cohesion    emotion       <NA>       <NA>       <NA>
#> 1350 coregulate    discuss    emotion       <NA>       <NA>       <NA>
#> 1351       plan  consensus       plan  consensus      adapt  consensus
#> 1352       plan  consensus       plan    monitor       plan       plan
#> 1353   cohesion  consensus    discuss  consensus       plan       plan
#> 1354 coregulate    discuss    emotion   cohesion  consensus coregulate
#> 1355  consensus       plan       plan       plan    emotion    emotion
#> 1356       plan  consensus    emotion   cohesion  consensus  consensus
#> 1357       plan       plan  consensus       plan coregulate   cohesion
#> 1358  consensus coregulate  consensus    discuss      adapt  consensus
#> 1359       plan       plan  consensus    monitor  synthesis  consensus
#> 1360  synthesis coregulate  consensus       plan       plan  consensus
#> 1361       plan  consensus coregulate    emotion   cohesion    monitor
#> 1362       plan  consensus    monitor  consensus    discuss    discuss
#> 1363    emotion   cohesion       plan       plan  consensus    emotion
#> 1364       plan       plan       plan       plan       plan  consensus
#> 1365       plan    emotion   cohesion    discuss       <NA>       <NA>
#> 1366    monitor   cohesion  consensus       <NA>       <NA>       <NA>
#> 1367    emotion    monitor coregulate       plan  consensus coregulate
#> 1368    discuss coregulate       <NA>       <NA>       <NA>       <NA>
#> 1369 coregulate    discuss    emotion coregulate   cohesion  consensus
#> 1370 coregulate  consensus coregulate    discuss  consensus       plan
#> 1371  synthesis      adapt  consensus  consensus       plan  consensus
#> 1372    monitor    emotion       <NA>       <NA>       <NA>       <NA>
#> 1373    emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1374       plan       plan       plan  consensus       plan       plan
#> 1375    discuss  consensus       plan    monitor coregulate       <NA>
#> 1376    emotion       plan    monitor    discuss    discuss  consensus
#> 1377       plan       plan       plan  consensus coregulate       plan
#> 1378      adapt  consensus       <NA>       <NA>       <NA>       <NA>
#> 1379  synthesis       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1380       plan       plan  consensus  consensus       plan  consensus
#> 1381    discuss  consensus    emotion   cohesion  consensus       plan
#> 1382    monitor    discuss   cohesion coregulate    emotion  consensus
#> 1383       plan coregulate       plan       plan  consensus    monitor
#> 1384  consensus       plan  consensus coregulate    discuss    discuss
#> 1385 coregulate    discuss coregulate    discuss    emotion   cohesion
#> 1386       plan    emotion  consensus    discuss      adapt       <NA>
#> 1387      adapt  consensus coregulate    discuss    monitor coregulate
#> 1388    discuss      adapt  consensus    emotion       <NA>       <NA>
#> 1389   cohesion       plan  consensus       plan  consensus       plan
#> 1390  consensus coregulate  consensus coregulate       <NA>       <NA>
#> 1391  consensus    emotion       plan    emotion       plan  consensus
#> 1392       plan    discuss  consensus       plan       <NA>       <NA>
#> 1393       plan       plan    discuss      adapt   cohesion    discuss
#> 1394 coregulate    monitor       plan       plan  consensus coregulate
#> 1395  consensus coregulate  synthesis  consensus    discuss coregulate
#> 1396    emotion       plan       plan       <NA>       <NA>       <NA>
#> 1397  consensus       plan  consensus       plan    emotion  consensus
#> 1398  consensus  consensus coregulate    emotion       plan    monitor
#> 1399    discuss    emotion   cohesion  consensus       plan       plan
#> 1400  synthesis  consensus    discuss coregulate       plan    emotion
#> 1401   cohesion    monitor  consensus coregulate    discuss  consensus
#> 1402       plan    discuss       <NA>       <NA>       <NA>       <NA>
#> 1403  consensus       plan       plan       plan    discuss  consensus
#> 1404      adapt   cohesion coregulate    discuss  synthesis  consensus
#> 1405 coregulate       plan       plan       plan       <NA>       <NA>
#> 1406    discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 1407    monitor      adapt  consensus       plan       plan       plan
#> 1408  synthesis    discuss  synthesis    emotion   cohesion  consensus
#> 1409  consensus    emotion   cohesion  consensus       plan  consensus
#> 1410 coregulate    discuss      adapt  consensus       plan       plan
#> 1411    discuss  consensus       plan       <NA>       <NA>       <NA>
#> 1412       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1413  consensus       plan  consensus coregulate    emotion   cohesion
#> 1414  consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 1415  consensus    discuss  consensus    discuss   cohesion coregulate
#> 1416    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1417       plan  consensus    monitor    discuss  synthesis    emotion
#> 1418  consensus  consensus    emotion       plan  consensus coregulate
#> 1419  consensus       plan       plan    emotion       plan       plan
#> 1420    emotion    discuss  consensus       plan       plan  consensus
#> 1421    emotion    monitor       <NA>       <NA>       <NA>       <NA>
#> 1422    emotion    emotion   cohesion  consensus       plan  consensus
#> 1423       plan    emotion coregulate    discuss  consensus       plan
#> 1424  consensus  consensus coregulate    emotion  consensus    discuss
#> 1425 coregulate    emotion       plan  consensus       plan  consensus
#> 1426 coregulate       plan       plan  consensus       plan  consensus
#> 1427   cohesion  consensus    monitor    discuss  synthesis  consensus
#> 1428 coregulate    emotion       plan       <NA>       <NA>       <NA>
#> 1429  consensus coregulate       plan       plan       plan  consensus
#> 1430    discuss    discuss   cohesion  consensus coregulate    discuss
#> 1431    emotion    discuss    discuss  synthesis      adapt    discuss
#> 1432       plan       plan       plan       plan       plan       plan
#> 1433    emotion   cohesion  consensus       plan       plan       plan
#> 1434 coregulate    monitor    discuss   cohesion    discuss    discuss
#> 1435       plan  consensus coregulate       plan  consensus    discuss
#> 1436    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1437  consensus       plan    discuss  consensus    emotion    emotion
#> 1438 coregulate  consensus  consensus       plan       plan  consensus
#> 1439    discuss      adapt  consensus  consensus    discuss  consensus
#> 1440       plan    emotion    discuss   cohesion  consensus       plan
#> 1441       plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 1442    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1443       plan  consensus    emotion   cohesion coregulate    discuss
#> 1444   cohesion  consensus   cohesion       plan       plan       plan
#> 1445    monitor   cohesion  consensus       plan    monitor  consensus
#> 1446       plan       plan       plan  consensus coregulate       plan
#> 1447    discuss    emotion  consensus  consensus       plan  consensus
#> 1448      adapt    discuss    monitor  consensus    discuss  synthesis
#> 1449       plan    emotion  consensus coregulate    discuss    discuss
#> 1450       plan  consensus coregulate    discuss  consensus coregulate
#> 1451  consensus       plan  consensus       plan       plan   cohesion
#> 1452      adapt    emotion   cohesion       plan       plan       plan
#> 1453  consensus       plan       plan    emotion  consensus    discuss
#> 1454       plan       plan    emotion   cohesion coregulate    emotion
#> 1455   cohesion    discuss    emotion coregulate    discuss    discuss
#> 1456    discuss  synthesis    discuss  consensus       plan    emotion
#> 1457    discuss  consensus coregulate    discuss    emotion  consensus
#> 1458       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1459    emotion  consensus       plan    discuss    emotion  consensus
#> 1460    monitor       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1461    discuss    monitor  synthesis      adapt   cohesion  consensus
#> 1462       plan  consensus    discuss  consensus       plan  consensus
#> 1463    discuss  synthesis  consensus    monitor    discuss  consensus
#> 1464       plan       plan       plan       <NA>       <NA>       <NA>
#> 1465    discuss  consensus coregulate   cohesion    emotion    discuss
#> 1466       plan       plan    monitor    discuss  synthesis  consensus
#> 1467       plan  consensus    discuss    discuss    discuss    emotion
#> 1468    discuss  consensus       plan    discuss    discuss    emotion
#> 1469  consensus    emotion   cohesion    discuss    discuss       <NA>
#> 1470       plan       plan  consensus       plan    monitor       plan
#> 1471    monitor    discuss  consensus    discuss  consensus       plan
#> 1472  synthesis  consensus  consensus       plan       plan  consensus
#> 1473 coregulate    monitor    discuss    discuss    discuss  synthesis
#> 1474    discuss    discuss      adapt    monitor       plan    discuss
#> 1475   cohesion  consensus       plan    discuss      adapt   cohesion
#> 1476   cohesion  consensus coregulate    discuss      adapt   cohesion
#> 1477       plan       plan       plan       plan  consensus coregulate
#> 1478  consensus  consensus       <NA>       <NA>       <NA>       <NA>
#> 1479      adapt   cohesion    discuss    discuss    discuss    discuss
#> 1480       plan    emotion coregulate    discuss coregulate    discuss
#> 1481 coregulate       plan    emotion       <NA>       <NA>       <NA>
#> 1482  consensus       plan       plan       plan       plan       plan
#> 1483    emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 1484    emotion coregulate    monitor    emotion   cohesion  consensus
#> 1485    discuss coregulate    discuss  consensus  consensus coregulate
#> 1486    emotion  consensus coregulate  consensus    discuss    discuss
#> 1487    discuss    discuss  consensus coregulate       plan       plan
#> 1488    discuss    emotion  consensus       plan    emotion       plan
#> 1489    monitor    discuss  consensus  consensus       plan  consensus
#> 1490   cohesion       plan   cohesion       plan    emotion   cohesion
#> 1491 coregulate  consensus       <NA>       <NA>       <NA>       <NA>
#> 1492  consensus coregulate coregulate    emotion    discuss  synthesis
#> 1493    discuss    discuss       <NA>       <NA>       <NA>       <NA>
#> 1494   cohesion       plan    emotion  consensus    discuss  synthesis
#> 1495       plan       plan  consensus  synthesis    emotion       plan
#> 1496   cohesion    discuss coregulate       plan  consensus    monitor
#> 1497       plan       plan       plan   cohesion    discuss      adapt
#> 1498  consensus    monitor    discuss    discuss   cohesion coregulate
#> 1499    discuss  synthesis  consensus coregulate    emotion   cohesion
#> 1500 coregulate    emotion    monitor  consensus       plan  consensus
#> 1501    emotion   cohesion  consensus       plan  consensus       plan
#> 1502  consensus       plan  consensus       plan       plan  consensus
#> 1503  consensus       plan       plan       plan       plan    discuss
#> 1504 coregulate   cohesion    monitor    emotion   cohesion coregulate
#> 1505  consensus    emotion  consensus       plan       plan       plan
#> 1506       plan    monitor    monitor       plan       plan       plan
#> 1507   cohesion    discuss  synthesis   cohesion       plan       plan
#> 1508 coregulate       plan       plan  consensus       plan       plan
#> 1509  consensus       plan    discuss    discuss       plan  consensus
#> 1510       plan  consensus coregulate    discuss  synthesis      adapt
#> 1511  consensus    discuss  synthesis    emotion       <NA>       <NA>
#> 1512   cohesion  consensus       plan       plan  consensus  consensus
#> 1513   cohesion       plan       <NA>       <NA>       <NA>       <NA>
#> 1514       plan       plan    emotion   cohesion       <NA>       <NA>
#> 1515      adapt  consensus    monitor    discuss      adapt   cohesion
#> 1516    monitor       plan    emotion  consensus    monitor    discuss
#> 1517 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1518  consensus       plan    discuss    discuss  synthesis  consensus
#> 1519    monitor  consensus    discuss coregulate  consensus coregulate
#> 1520    discuss  consensus    monitor    discuss  synthesis  consensus
#> 1521  consensus       plan  consensus    discuss    discuss       <NA>
#> 1522  consensus coregulate    emotion    emotion       <NA>       <NA>
#> 1523  synthesis  consensus    discuss  synthesis  consensus       plan
#> 1524    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1525  synthesis    monitor    discuss    emotion  consensus       plan
#> 1526       plan    emotion   cohesion       plan    emotion  consensus
#> 1527       plan       plan  consensus    monitor       plan       plan
#> 1528  consensus       plan    emotion   cohesion  consensus coregulate
#> 1529    discuss  synthesis       <NA>       <NA>       <NA>       <NA>
#> 1530       plan       plan   cohesion  consensus    monitor       plan
#> 1531  consensus    discuss      adapt    emotion       plan    monitor
#> 1532  consensus       plan    emotion  consensus       plan  consensus
#> 1533    emotion   cohesion  consensus       plan       plan    discuss
#> 1534      adapt  consensus    monitor  consensus       plan  consensus
#> 1535    discuss  synthesis  consensus coregulate  consensus    discuss
#> 1536  consensus    discuss  consensus    discuss  consensus  consensus
#> 1537    emotion   cohesion    monitor       plan  consensus       plan
#> 1538      adapt    emotion    monitor   cohesion coregulate       <NA>
#> 1539 coregulate    discuss coregulate       plan  consensus       plan
#> 1540       plan  consensus       plan  consensus    discuss  consensus
#> 1541   cohesion coregulate    discuss       <NA>       <NA>       <NA>
#> 1542       plan       plan  consensus    discuss    discuss coregulate
#> 1543       plan       plan  consensus       <NA>       <NA>       <NA>
#> 1544    discuss    discuss      adapt   cohesion  consensus coregulate
#> 1545  consensus coregulate    discuss    emotion   cohesion    emotion
#> 1546       plan       plan    monitor       plan    emotion    emotion
#> 1547  consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 1548    monitor  consensus coregulate  consensus       plan       <NA>
#> 1549       plan       plan       plan       plan    emotion    emotion
#> 1550       plan  consensus       plan   cohesion coregulate       plan
#> 1551    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1552       plan       plan  consensus       plan       plan    monitor
#> 1553 coregulate    emotion  consensus       plan       plan       <NA>
#> 1554   cohesion  consensus coregulate    discuss    discuss  synthesis
#> 1555    emotion coregulate  consensus  consensus coregulate    emotion
#> 1556    emotion  consensus       plan    monitor    discuss       <NA>
#> 1557   cohesion  consensus    monitor    discuss       plan  consensus
#> 1558  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1559       plan  consensus coregulate       plan    emotion    emotion
#> 1560    discuss    emotion   cohesion    discuss  consensus coregulate
#> 1561  consensus    emotion    monitor    discuss  consensus       plan
#> 1562  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1563  consensus    discuss      adapt  consensus    discuss  consensus
#> 1564  consensus       plan       plan       <NA>       <NA>       <NA>
#> 1565    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1566    emotion       plan  consensus       plan    monitor    discuss
#> 1567 coregulate    monitor    discuss   cohesion  consensus       plan
#> 1568       plan    emotion   cohesion    monitor    discuss    discuss
#> 1569       plan    monitor  consensus coregulate    emotion    discuss
#> 1570  consensus   cohesion coregulate    discuss    discuss coregulate
#> 1571      adapt   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1572       plan       plan  consensus coregulate       plan       plan
#> 1573    discuss  synthesis      adapt  consensus       plan  consensus
#> 1574   cohesion coregulate       <NA>       <NA>       <NA>       <NA>
#> 1575    emotion coregulate    discuss coregulate    monitor    discuss
#> 1576    discuss   cohesion  consensus coregulate  consensus       plan
#> 1577  consensus       plan    emotion   cohesion  consensus    discuss
#> 1578  synthesis      adapt  consensus       plan       plan       plan
#> 1579       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1580 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1581      adapt   cohesion       plan    monitor  consensus       plan
#> 1582 coregulate    emotion    emotion coregulate    discuss  consensus
#> 1583      adapt   cohesion       plan       plan       plan       plan
#> 1584    monitor    discuss    emotion  consensus  consensus coregulate
#> 1585    monitor    discuss       <NA>       <NA>       <NA>       <NA>
#> 1586    discuss    discuss    discuss coregulate       plan  consensus
#> 1587    emotion  consensus       plan       plan    emotion    emotion
#> 1588  consensus       plan       plan  consensus       plan  consensus
#> 1589  consensus coregulate  consensus       plan       plan  consensus
#> 1590  consensus       plan  consensus coregulate       plan       plan
#> 1591  consensus       plan    discuss  consensus coregulate       plan
#> 1592    monitor       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1593       plan  consensus       plan    monitor    discuss      adapt
#> 1594 coregulate    emotion    monitor  synthesis      adapt  consensus
#> 1595 coregulate    discuss    monitor coregulate coregulate    discuss
#> 1596  consensus       plan    emotion   cohesion  consensus       plan
#> 1597       plan       plan  consensus  consensus coregulate    discuss
#> 1598       plan       plan       plan       plan  consensus coregulate
#> 1599    discuss  synthesis   cohesion  consensus coregulate    discuss
#> 1600  consensus       plan    discuss    discuss      adapt  consensus
#> 1601    discuss coregulate    discuss  synthesis  consensus    discuss
#> 1602       plan  consensus coregulate  consensus       plan    emotion
#> 1603  consensus    discuss  synthesis coregulate       plan       plan
#> 1604       plan       plan  consensus coregulate    discuss    emotion
#> 1605 coregulate    discuss  synthesis    discuss  consensus    emotion
#> 1606   cohesion  consensus       plan    monitor    discuss  consensus
#> 1607    discuss    discuss  consensus       plan       plan       plan
#> 1608    discuss      adapt   cohesion  consensus       plan  consensus
#> 1609  consensus       plan  consensus       plan    emotion  consensus
#> 1610 coregulate       plan    emotion    discuss    monitor    discuss
#> 1611  consensus       plan       plan  consensus       plan       plan
#> 1612    monitor       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1613   cohesion coregulate  consensus coregulate    discuss  consensus
#> 1614 coregulate    discuss      adapt    emotion  consensus       plan
#> 1615       plan    monitor       plan  consensus    monitor  consensus
#> 1616       plan  consensus       plan    emotion   cohesion  consensus
#> 1617    monitor    emotion  consensus  consensus  consensus    monitor
#> 1618       plan    emotion coregulate    discuss    emotion    emotion
#> 1619   cohesion    emotion       <NA>       <NA>       <NA>       <NA>
#> 1620  consensus    discuss coregulate    discuss    discuss  consensus
#> 1621    emotion  consensus       plan       <NA>       <NA>       <NA>
#> 1622    discuss    discuss    discuss    discuss  consensus       plan
#> 1623  consensus       plan       plan  consensus    emotion  consensus
#> 1624       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1625       plan  consensus coregulate       plan  consensus      adapt
#> 1626    monitor    discuss      adapt    monitor    discuss    monitor
#> 1627  consensus coregulate coregulate    monitor  consensus       plan
#> 1628   cohesion  consensus    discuss  consensus coregulate    discuss
#> 1629  consensus       plan  consensus   cohesion  consensus    discuss
#> 1630    discuss    monitor       plan  consensus       plan       plan
#> 1631      adapt  consensus       plan       plan       plan       plan
#> 1632  consensus coregulate       plan  consensus    discuss  synthesis
#> 1633  consensus    discuss    monitor    emotion   cohesion coregulate
#> 1634    monitor    discuss  consensus       plan       plan  consensus
#> 1635       plan  consensus       plan    monitor    discuss    discuss
#> 1636  consensus       plan   cohesion  consensus coregulate    discuss
#> 1637 coregulate    monitor       plan    discuss  synthesis    emotion
#> 1638    emotion   cohesion    emotion  consensus       plan       plan
#> 1639  consensus       plan       plan    emotion  consensus       plan
#> 1640       plan       plan    monitor       plan  consensus  consensus
#> 1641  synthesis  consensus       plan       plan    emotion  consensus
#> 1642  consensus coregulate   cohesion  consensus       plan    monitor
#> 1643       plan       plan       plan  consensus coregulate  consensus
#> 1644  consensus       plan    discuss  synthesis    discuss    emotion
#> 1645       plan   cohesion  consensus       plan  consensus  consensus
#> 1646       plan       plan  consensus    discuss coregulate    discuss
#> 1647  consensus       plan  consensus    emotion   cohesion       plan
#> 1648  consensus       plan       plan       plan       plan    emotion
#> 1649    discuss    discuss  consensus       plan    discuss       <NA>
#> 1650       plan    emotion  consensus       plan    monitor coregulate
#> 1651   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1652    emotion       plan       plan  consensus coregulate      adapt
#> 1653  consensus       plan  consensus    discuss    emotion    emotion
#> 1654    discuss    emotion  consensus    emotion  consensus    discuss
#> 1655       plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 1656       plan    discuss  consensus       plan    emotion       plan
#> 1657    emotion  consensus coregulate       plan       plan    monitor
#> 1658 coregulate    discuss  consensus       plan  consensus    discuss
#> 1659  synthesis      adapt    emotion    emotion       plan  consensus
#> 1660    discuss    discuss   cohesion  consensus    discuss    emotion
#> 1661    discuss      adapt    monitor       plan    monitor    discuss
#> 1662  consensus       plan       plan       plan       plan    monitor
#> 1663    discuss    discuss coregulate    discuss  synthesis    emotion
#> 1664       plan   cohesion    emotion       plan       plan    discuss
#> 1665  consensus    discuss    emotion  consensus       plan       plan
#> 1666   cohesion  consensus coregulate    discuss coregulate       plan
#> 1667       plan  consensus coregulate       plan    emotion  consensus
#> 1668 coregulate    discuss coregulate       <NA>       <NA>       <NA>
#> 1669      adapt  consensus    discuss    discuss    discuss       <NA>
#> 1670    discuss  synthesis  consensus coregulate       <NA>       <NA>
#> 1671    monitor    discuss   cohesion    emotion   cohesion  consensus
#> 1672  consensus coregulate    emotion  consensus       plan    discuss
#> 1673   cohesion  consensus coregulate    discuss  synthesis      adapt
#> 1674  consensus    discuss  synthesis      adapt   cohesion       <NA>
#> 1675       plan coregulate  consensus coregulate    discuss    discuss
#> 1676    discuss  synthesis    emotion   cohesion       <NA>       <NA>
#> 1677   cohesion coregulate    discuss    emotion       plan  consensus
#> 1678    emotion   cohesion       plan       plan    emotion coregulate
#> 1679      adapt    monitor       <NA>       <NA>       <NA>       <NA>
#> 1680    emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 1681  consensus       plan       plan    emotion   cohesion       plan
#> 1682    discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 1683      adapt  consensus  consensus       plan       plan    discuss
#> 1684       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1685       plan       plan       plan  consensus       plan       <NA>
#> 1686    discuss    emotion       plan  consensus  consensus    discuss
#> 1687  consensus       plan       plan  consensus    discuss    discuss
#> 1688    emotion    emotion    emotion   cohesion    discuss  consensus
#> 1689  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1690 coregulate       plan       plan  consensus       plan       plan
#> 1691    discuss  consensus       plan       plan    monitor   cohesion
#> 1692       plan  consensus    emotion       <NA>       <NA>       <NA>
#> 1693    discuss    discuss coregulate       plan    emotion  consensus
#> 1694  synthesis coregulate    discuss  consensus       plan       plan
#> 1695       plan       plan       plan    emotion   cohesion    emotion
#> 1696       plan       plan       plan    monitor       plan       plan
#> 1697    discuss  synthesis      adapt   cohesion       plan       plan
#> 1698      adapt  consensus coregulate    discuss  synthesis  consensus
#> 1699  consensus    discuss       plan  consensus    emotion    emotion
#> 1700 coregulate   cohesion coregulate  synthesis    discuss    discuss
#> 1701    discuss  consensus  consensus       plan       plan  consensus
#> 1702       plan       plan  consensus       plan       <NA>       <NA>
#> 1703  consensus       plan       plan  consensus coregulate       plan
#> 1704       plan    discuss    discuss  consensus       plan  consensus
#> 1705  consensus  consensus coregulate       plan       plan       plan
#> 1706    discuss  synthesis   cohesion    emotion   cohesion  consensus
#> 1707       plan    monitor coregulate       plan       plan       plan
#> 1708    discuss    discuss    discuss  consensus coregulate    monitor
#> 1709       plan  consensus    emotion   cohesion       <NA>       <NA>
#> 1710    monitor coregulate       <NA>       <NA>       <NA>       <NA>
#> 1711       plan  consensus    discuss    discuss    discuss    discuss
#> 1712       plan       plan    monitor   cohesion coregulate    discuss
#> 1713       plan       plan    emotion  consensus  consensus       plan
#> 1714 coregulate    discuss    discuss coregulate    discuss    discuss
#> 1715    emotion    discuss  synthesis      adapt    monitor    discuss
#> 1716  consensus       plan       plan       plan    monitor    discuss
#> 1717  consensus  synthesis    discuss coregulate    emotion    monitor
#> 1718      adapt  consensus  consensus       plan    monitor    discuss
#> 1719    emotion coregulate    discuss  synthesis      adapt   cohesion
#> 1720    discuss  consensus       plan  consensus    discuss  consensus
#> 1721  consensus    discuss    emotion   cohesion  consensus       <NA>
#> 1722    emotion   cohesion       plan  consensus       plan    emotion
#> 1723       plan    emotion       plan       plan    emotion       <NA>
#> 1724  consensus coregulate    discuss    discuss  consensus coregulate
#> 1725    discuss    discuss       <NA>       <NA>       <NA>       <NA>
#> 1726       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1727    emotion  consensus  consensus  consensus    emotion   cohesion
#> 1728       plan       plan    discuss    discuss    discuss      adapt
#> 1729 coregulate    discuss  consensus       plan  consensus coregulate
#> 1730    emotion       plan  consensus       plan    monitor       plan
#> 1731    monitor  consensus    monitor    emotion   cohesion  consensus
#> 1732    monitor    emotion  consensus       plan  consensus coregulate
#> 1733   cohesion  consensus coregulate       plan       plan       plan
#> 1734    discuss    discuss  synthesis  consensus coregulate    discuss
#> 1735    discuss      adapt  consensus coregulate    discuss    discuss
#> 1736    emotion   cohesion  consensus    emotion    monitor   cohesion
#> 1737       plan    emotion    emotion   cohesion    discuss    discuss
#> 1738      adapt  consensus       plan       plan       plan       plan
#> 1739       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1740  synthesis      adapt    monitor    discuss   cohesion  consensus
#> 1741      adapt       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1742    discuss    discuss      adapt  consensus coregulate    discuss
#> 1743 coregulate   cohesion  consensus    emotion    monitor    discuss
#> 1744   cohesion coregulate    monitor    discuss  synthesis  consensus
#> 1745   cohesion    discuss       <NA>       <NA>       <NA>       <NA>
#> 1746       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1747 coregulate    monitor    discuss  consensus coregulate  consensus
#> 1748       plan    discuss      adapt  consensus    emotion  consensus
#> 1749    discuss       plan       plan  consensus coregulate    discuss
#> 1750       plan  consensus       plan  consensus    monitor       plan
#> 1751 coregulate    discuss    emotion  consensus  consensus coregulate
#> 1752 coregulate       plan   cohesion  consensus       plan       plan
#> 1753    emotion    monitor    discuss coregulate      adapt    monitor
#> 1754  consensus       plan    emotion       plan    monitor  consensus
#> 1755  consensus       plan  consensus coregulate       plan  consensus
#> 1756       plan       plan    discuss  consensus    emotion  consensus
#> 1757    discuss    discuss  consensus    emotion       plan    discuss
#> 1758    discuss    discuss coregulate       plan       plan       plan
#> 1759  consensus    emotion  consensus coregulate   cohesion  consensus
#> 1760    discuss  consensus       plan  consensus    discuss coregulate
#> 1761    emotion  consensus       plan  consensus       plan       plan
#> 1762      adapt    emotion   cohesion    monitor  consensus coregulate
#> 1763       plan  consensus       plan       plan  consensus       plan
#> 1764    emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1765 coregulate    emotion    emotion  consensus    discuss  synthesis
#> 1766    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1767   cohesion  consensus coregulate    discuss    emotion   cohesion
#> 1768  consensus       plan    discuss  synthesis   cohesion  consensus
#> 1769    discuss    emotion  consensus coregulate  consensus    discuss
#> 1770 coregulate coregulate       plan    monitor coregulate  consensus
#> 1771  consensus    emotion    monitor       plan       plan    emotion
#> 1772       plan    discuss    emotion       <NA>       <NA>       <NA>
#> 1773   cohesion coregulate    discuss    emotion  consensus coregulate
#> 1774       plan  consensus  consensus       plan       plan    discuss
#> 1775  consensus    discuss coregulate    discuss      adapt    discuss
#> 1776       plan       plan       plan    monitor    discuss  consensus
#> 1777  consensus    monitor    discuss    emotion       plan    discuss
#> 1778 coregulate    monitor       <NA>       <NA>       <NA>       <NA>
#> 1779    discuss    discuss  synthesis  consensus    discuss    discuss
#> 1780    discuss  consensus    emotion  consensus    discuss    discuss
#> 1781    discuss    discuss    discuss      adapt  consensus coregulate
#> 1782  synthesis  consensus  consensus       plan       plan  consensus
#> 1783    discuss    discuss    emotion    monitor    discuss  synthesis
#> 1784      adapt  consensus       plan       plan  consensus coregulate
#> 1785       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1786       plan  consensus       plan       plan       plan  consensus
#> 1787       plan       plan       plan       <NA>       <NA>       <NA>
#> 1788 coregulate       plan    monitor  consensus       plan       plan
#> 1789       plan       plan  consensus    discuss       <NA>       <NA>
#> 1790  synthesis coregulate   cohesion    emotion  consensus coregulate
#> 1791      adapt    emotion    emotion  consensus       plan  consensus
#> 1792    discuss      adapt   cohesion  consensus    emotion  consensus
#> 1793  consensus       plan    emotion    discuss       <NA>       <NA>
#> 1794    discuss    discuss  consensus       plan       plan  consensus
#> 1795    emotion  consensus   cohesion    emotion       plan    discuss
#> 1796  consensus       plan       plan coregulate    discuss  synthesis
#> 1797      adapt   cohesion  consensus       <NA>       <NA>       <NA>
#> 1798  synthesis      adapt    discuss  synthesis  consensus       plan
#> 1799 coregulate    discuss    emotion  consensus       plan    discuss
#> 1800       plan    discuss    emotion  consensus       plan       plan
#> 1801  consensus       plan    emotion  consensus  consensus       plan
#> 1802   cohesion    discuss      adapt    emotion       plan  consensus
#> 1803    monitor       plan    emotion   cohesion coregulate coregulate
#> 1804 coregulate    monitor    discuss  synthesis    discuss    discuss
#> 1805    monitor    emotion   cohesion       plan   cohesion       plan
#> 1806       plan  consensus       plan    discuss      adapt  consensus
#> 1807    discuss  synthesis    monitor    discuss    emotion   cohesion
#> 1808    monitor       plan    emotion    discuss    emotion  consensus
#> 1809  consensus       plan  consensus    discuss    emotion  consensus
#> 1810  consensus coregulate    discuss  consensus       plan  consensus
#> 1811    emotion       plan  consensus    monitor       plan    emotion
#> 1812       plan       plan  consensus coregulate  consensus coregulate
#> 1813    monitor coregulate       plan       plan    emotion   cohesion
#> 1814    discuss  synthesis      adapt  consensus       plan       plan
#> 1815  consensus       plan       plan    monitor coregulate   cohesion
#> 1816  consensus       plan    monitor      adapt   cohesion    monitor
#> 1817       plan    emotion    emotion    emotion    monitor  consensus
#> 1818    discuss  consensus    emotion       plan    emotion   cohesion
#> 1819  consensus       plan  consensus       plan       plan  consensus
#> 1820 coregulate  consensus coregulate   cohesion    emotion  consensus
#> 1821       plan       plan  consensus       <NA>       <NA>       <NA>
#> 1822       plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 1823    discuss      adapt   cohesion  consensus       plan  consensus
#> 1824       plan    emotion   cohesion    emotion       plan    monitor
#> 1825    emotion  consensus       plan  consensus coregulate       plan
#> 1826    monitor   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1827       plan  consensus       plan    discuss      adapt  consensus
#> 1828    discuss coregulate  consensus coregulate    emotion coregulate
#> 1829    emotion       plan   cohesion  consensus    emotion       plan
#> 1830  consensus    discuss    discuss      adapt  consensus coregulate
#> 1831  consensus  consensus    emotion   cohesion    emotion    monitor
#> 1832       plan  consensus       plan  consensus       plan    emotion
#> 1833   cohesion    emotion   cohesion  consensus       plan  consensus
#> 1834    emotion   cohesion coregulate       <NA>       <NA>       <NA>
#> 1835       plan       plan    monitor  consensus    discuss  synthesis
#> 1836  consensus coregulate       plan    monitor    discuss       plan
#> 1837  consensus    emotion   cohesion       plan  consensus  consensus
#> 1838    discuss  consensus    emotion  consensus coregulate coregulate
#> 1839    monitor    discuss      adapt    discuss coregulate  consensus
#> 1840       plan  consensus coregulate  consensus       plan       plan
#> 1841    discuss  consensus  consensus coregulate       plan    emotion
#> 1842       plan    monitor       plan       plan       plan  consensus
#> 1843   cohesion  consensus  consensus    discuss  consensus       plan
#> 1844    discuss    emotion  consensus    discuss  synthesis  consensus
#> 1845       plan  consensus coregulate    monitor    discuss  synthesis
#> 1846       plan       plan    discuss coregulate    discuss  consensus
#> 1847  consensus    discuss    discuss  synthesis      adapt   cohesion
#> 1848       plan    monitor    discuss  consensus    discuss  synthesis
#> 1849    emotion  consensus coregulate       <NA>       <NA>       <NA>
#> 1850       plan       plan    discuss  consensus    emotion       plan
#> 1851  consensus  synthesis  consensus    discuss       <NA>       <NA>
#> 1852   cohesion coregulate  consensus    discuss  synthesis      adapt
#> 1853 coregulate       plan       plan       plan       plan       plan
#> 1854    discuss    discuss coregulate    discuss coregulate    monitor
#> 1855    discuss coregulate    discuss    discuss    monitor    discuss
#> 1856  consensus coregulate    discuss      adapt   cohesion coregulate
#> 1857      adapt  consensus       plan  consensus       plan  consensus
#> 1858  consensus       plan  consensus    discuss    discuss  synthesis
#> 1859    emotion  consensus       plan  consensus    discuss  synthesis
#> 1860    emotion  consensus       plan       plan       <NA>       <NA>
#> 1861       plan       plan    emotion  consensus       plan  consensus
#> 1862 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1863    emotion  consensus coregulate    discuss    discuss  consensus
#> 1864       plan  consensus       plan  consensus    emotion   cohesion
#> 1865       plan       plan    monitor       plan       plan       plan
#> 1866 coregulate       plan    emotion   cohesion  consensus       plan
#> 1867  consensus       plan       plan  consensus       plan       plan
#> 1868   cohesion coregulate  consensus   cohesion  consensus    discuss
#> 1869    discuss      adapt  consensus coregulate  synthesis    emotion
#> 1870       plan       plan    emotion    emotion   cohesion    emotion
#> 1871       plan  consensus    emotion       plan    monitor    emotion
#> 1872       plan  consensus coregulate    discuss      adapt    discuss
#> 1873  consensus       plan       plan    emotion    emotion  consensus
#> 1874    monitor    monitor       <NA>       <NA>       <NA>       <NA>
#> 1875   cohesion    emotion    discuss    discuss    discuss    discuss
#> 1876    monitor  consensus       plan       plan       plan       plan
#> 1877       plan  consensus coregulate       plan   cohesion  consensus
#> 1878    discuss   cohesion  consensus    discuss   cohesion  consensus
#> 1879  consensus       plan       plan   cohesion  consensus   cohesion
#> 1880       plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 1881       plan  consensus coregulate    discuss  consensus coregulate
#> 1882    discuss    discuss  synthesis coregulate  consensus coregulate
#> 1883       plan       plan       plan    monitor   cohesion       plan
#> 1884  consensus       plan    emotion       plan  consensus    monitor
#> 1885       plan    emotion   cohesion coregulate    discuss      adapt
#> 1886       plan       plan       plan  consensus       plan  consensus
#> 1887    discuss    emotion  consensus       plan  consensus       plan
#> 1888    emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 1889    discuss  synthesis coregulate    discuss    discuss    discuss
#> 1890       plan  consensus       plan  consensus  consensus coregulate
#> 1891  consensus       plan  consensus       plan       <NA>       <NA>
#> 1892    emotion   cohesion  consensus  consensus    monitor    discuss
#> 1893   cohesion    emotion   cohesion    monitor       plan       plan
#> 1894       plan    discuss  synthesis       <NA>       <NA>       <NA>
#> 1895  consensus       plan       plan    monitor    discuss  consensus
#> 1896       plan    emotion   cohesion    emotion   cohesion  consensus
#> 1897       plan  consensus       plan  consensus    monitor    emotion
#> 1898 coregulate    discuss  synthesis  consensus  consensus       plan
#> 1899      adapt  consensus    discuss      adapt    discuss coregulate
#> 1900      adapt   cohesion  consensus    discuss coregulate       plan
#> 1901      adapt   cohesion coregulate       plan       plan       plan
#> 1902   cohesion    discuss    discuss  consensus coregulate    emotion
#> 1903       plan       plan    emotion       <NA>       <NA>       <NA>
#> 1904    emotion   cohesion       plan       plan       plan   cohesion
#> 1905    monitor       plan  consensus coregulate    discuss  consensus
#> 1906    emotion    emotion  consensus    discuss    discuss      adapt
#> 1907    discuss   cohesion    discuss   cohesion    discuss  consensus
#> 1908    emotion  consensus       plan       plan       plan       plan
#> 1909    discuss  consensus       plan  consensus    discuss  consensus
#> 1910  consensus       plan       plan       plan       plan       plan
#> 1911   cohesion       plan  consensus       plan    emotion    monitor
#> 1912  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1913   cohesion  consensus       plan  consensus       plan       plan
#> 1914    discuss    discuss    discuss  consensus       plan       plan
#> 1915       plan       plan  consensus       plan       <NA>       <NA>
#> 1916       plan    monitor    discuss  synthesis      adapt  consensus
#> 1917       plan    discuss    monitor   cohesion  consensus       <NA>
#> 1918       plan    discuss    emotion   cohesion coregulate       plan
#> 1919    discuss    discuss  consensus coregulate    discuss  synthesis
#> 1920    discuss    discuss    monitor       plan  consensus  consensus
#> 1921       plan       plan    discuss  consensus       plan    emotion
#> 1922  synthesis  consensus       plan       plan       plan  consensus
#> 1923      adapt  consensus       plan       plan       plan       plan
#> 1924  consensus       plan    monitor   cohesion    discuss    emotion
#> 1925       plan  consensus  consensus       plan  consensus       plan
#> 1926  consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 1927       plan  consensus       plan  consensus    emotion   cohesion
#> 1928    discuss    discuss coregulate   cohesion    monitor    discuss
#> 1929  consensus coregulate    monitor  consensus coregulate       plan
#> 1930       plan    emotion    discuss  consensus    emotion   cohesion
#> 1931    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1932  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1933       plan    discuss      adapt  consensus    emotion coregulate
#> 1934      adapt  consensus       plan    emotion  consensus    monitor
#> 1935       plan   cohesion    emotion   cohesion  consensus coregulate
#> 1936   cohesion       plan       plan       plan       plan       plan
#> 1937   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1938  consensus  synthesis   cohesion  consensus  consensus       plan
#> 1939    emotion    discuss    discuss       <NA>       <NA>       <NA>
#> 1940 coregulate    discuss    emotion  consensus       plan       <NA>
#> 1941  consensus    discuss    discuss       <NA>       <NA>       <NA>
#> 1942    discuss    discuss  consensus    monitor       <NA>       <NA>
#> 1943   cohesion  consensus       plan       plan  consensus       plan
#> 1944    emotion   cohesion    emotion  consensus       plan       plan
#> 1945 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1946       plan    discuss    discuss  consensus       plan       plan
#> 1947    monitor    discuss  synthesis  consensus       plan  consensus
#> 1948  consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 1949  consensus       plan  consensus       plan  consensus    discuss
#> 1950    monitor   cohesion  consensus    discuss    discuss   cohesion
#> 1951       plan       plan  consensus       <NA>       <NA>       <NA>
#> 1952  consensus  consensus       plan  consensus       <NA>       <NA>
#> 1953       plan    discuss coregulate   cohesion       plan       plan
#> 1954       plan  consensus    discuss      adapt   cohesion coregulate
#> 1955    discuss      adapt   cohesion  consensus    monitor       <NA>
#> 1956  consensus  consensus       plan  consensus coregulate  consensus
#> 1957       plan  consensus       plan       plan  consensus    emotion
#> 1958  consensus coregulate  consensus       plan  consensus coregulate
#> 1959   cohesion  consensus       plan       <NA>       <NA>       <NA>
#> 1960    discuss    emotion  consensus       plan       plan  consensus
#> 1961    emotion   cohesion  consensus       plan       plan       plan
#> 1962    discuss coregulate    emotion coregulate    emotion   cohesion
#> 1963    discuss  synthesis  consensus       plan       <NA>       <NA>
#> 1964 coregulate    discuss    discuss      adapt  consensus       plan
#> 1965    discuss    discuss  synthesis      adapt    discuss  consensus
#> 1966    discuss    discuss coregulate    monitor  consensus       plan
#> 1967    discuss  consensus       plan       plan    monitor    emotion
#> 1968    discuss  synthesis  consensus coregulate  consensus  consensus
#> 1969    emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 1970 coregulate    emotion   cohesion    discuss    discuss    discuss
#> 1971  consensus       plan    monitor  consensus  consensus coregulate
#> 1972  consensus    discuss  synthesis      adapt       plan  consensus
#> 1973       plan    emotion    discuss  synthesis   cohesion       <NA>
#> 1974  consensus coregulate   cohesion   cohesion  consensus       plan
#> 1975   cohesion coregulate    discuss  consensus       plan       plan
#> 1976       plan   cohesion coregulate    emotion       plan    monitor
#> 1977       plan       plan       plan  consensus       plan  consensus
#> 1978  consensus coregulate  consensus       plan  consensus    discuss
#> 1979    discuss    monitor coregulate    monitor    emotion   cohesion
#> 1980 coregulate    discuss    discuss      adapt  consensus coregulate
#> 1981    discuss  synthesis    emotion  consensus coregulate    monitor
#> 1982    discuss coregulate    emotion   cohesion       plan  consensus
#> 1983  synthesis  consensus  consensus       plan  consensus       <NA>
#> 1984  consensus       plan       plan    monitor       <NA>       <NA>
#> 1985  consensus    discuss  synthesis    discuss  synthesis    monitor
#> 1986   cohesion  consensus    monitor  consensus    monitor    discuss
#> 1987    emotion  consensus       plan    monitor    discuss    discuss
#> 1988    monitor    discuss      adapt  consensus coregulate  consensus
#> 1989       plan       plan    emotion   cohesion  consensus coregulate
#> 1990    discuss    discuss    emotion   cohesion  consensus       plan
#> 1991   cohesion  consensus       plan  consensus       plan       plan
#> 1992    emotion    emotion    emotion  consensus       <NA>       <NA>
#> 1993    discuss  consensus       plan    monitor  consensus       plan
#> 1994    emotion   cohesion    monitor    discuss       plan    discuss
#> 1995    discuss  synthesis    discuss coregulate    emotion    emotion
#> 1996       plan    discuss    discuss  synthesis coregulate    discuss
#> 1997  consensus       plan    monitor  consensus    discuss  consensus
#> 1998    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1999       plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 2000       plan       plan    emotion  consensus       plan       plan
#>              T8         T9        T10        T11        T12        T13
#> 1     consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 2          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 3       monitor  consensus coregulate    emotion   cohesion   cohesion
#> 4      cohesion  consensus       plan       plan    emotion       <NA>
#> 5          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 6          plan    emotion       plan       plan  consensus coregulate
#> 7          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 8          plan    monitor    discuss  consensus       plan  consensus
#> 9          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 10         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 11      emotion       plan    emotion   cohesion  consensus    discuss
#> 12    consensus    discuss    discuss   cohesion    monitor  consensus
#> 13         plan    monitor   cohesion  consensus    discuss  synthesis
#> 14         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 15      discuss  synthesis  consensus    discuss  consensus    discuss
#> 16      emotion    discuss   cohesion    emotion   cohesion coregulate
#> 17         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 18    consensus       plan  consensus       plan       plan  consensus
#> 19         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 20         plan       plan       plan       <NA>       <NA>       <NA>
#> 21    consensus       plan  consensus    emotion   cohesion  consensus
#> 22         plan coregulate  consensus coregulate       plan  consensus
#> 23         plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 24         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 25         plan   cohesion  consensus       plan       plan    emotion
#> 26         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 27     cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 28      emotion      adapt  consensus       plan  consensus       plan
#> 29         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 30    synthesis       plan    monitor  consensus    discuss  consensus
#> 31         plan    emotion  consensus coregulate    discuss  consensus
#> 32      emotion    discuss  consensus       plan  consensus    monitor
#> 33    consensus    emotion   cohesion  consensus       plan    monitor
#> 34    consensus    discuss    discuss coregulate    emotion   cohesion
#> 35    synthesis  consensus   cohesion  consensus    discuss  consensus
#> 36     cohesion  consensus  consensus    discuss    emotion  consensus
#> 37     cohesion  consensus    discuss    discuss  consensus       plan
#> 38   coregulate       plan    emotion       plan       <NA>       <NA>
#> 39         plan       plan  consensus       plan       plan coregulate
#> 40         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 41         plan    emotion  consensus coregulate    monitor    discuss
#> 42      monitor  synthesis  consensus coregulate       plan  consensus
#> 43         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 44         plan  consensus    discuss  consensus coregulate       plan
#> 45    consensus       plan       plan    emotion   cohesion       plan
#> 46     cohesion  consensus  consensus       plan  consensus       <NA>
#> 47        adapt   cohesion    emotion  consensus       plan    monitor
#> 48    consensus    discuss  consensus       plan       plan    monitor
#> 49    consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 50    consensus       plan       plan    monitor  consensus coregulate
#> 51         plan       plan    emotion    discuss  consensus    discuss
#> 52      discuss  consensus       plan       plan  consensus coregulate
#> 53    synthesis      adapt  consensus coregulate       plan  consensus
#> 54         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 55    consensus       plan    monitor       plan    discuss coregulate
#> 56         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 57      emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 58      monitor coregulate    monitor    discuss coregulate       plan
#> 59         plan coregulate       plan       plan    discuss    discuss
#> 60         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 61         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 62         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 63      discuss    discuss  synthesis  consensus       plan       plan
#> 64         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 65         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 66    consensus    discuss  consensus       plan       plan    emotion
#> 67    consensus    discuss coregulate    monitor      adapt  consensus
#> 68    consensus       plan    emotion       plan       plan       <NA>
#> 69         plan    emotion    discuss    discuss    discuss       <NA>
#> 70         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 71      discuss    discuss    emotion  consensus coregulate    monitor
#> 72    consensus       plan  consensus       plan       plan    discuss
#> 73         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 74      discuss  synthesis       plan    monitor    discuss    emotion
#> 75     cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 76         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 77    consensus       plan  consensus coregulate    emotion   cohesion
#> 78         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 79    consensus   cohesion coregulate    emotion    discuss    discuss
#> 80         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 81         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 82         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 83    consensus    discuss    discuss  synthesis  consensus    discuss
#> 84         plan       plan       plan    emotion  consensus       plan
#> 85    consensus    emotion  consensus    discuss  consensus       plan
#> 86         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 87      monitor       plan coregulate  consensus    emotion   cohesion
#> 88   coregulate   cohesion    emotion  consensus    discuss    emotion
#> 89         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 90    consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 91      monitor    discuss  consensus  consensus       <NA>       <NA>
#> 92      monitor    discuss   cohesion  consensus coregulate    monitor
#> 93    consensus    emotion   cohesion  consensus    discuss    emotion
#> 94         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 95    consensus       plan       plan    discuss  synthesis  consensus
#> 96      discuss  consensus    discuss    emotion  consensus    discuss
#> 97      emotion  consensus       plan   cohesion  consensus       plan
#> 98    consensus  consensus  consensus    discuss    emotion    discuss
#> 99    consensus       plan    emotion  consensus    emotion   cohesion
#> 100        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 101   consensus  consensus  consensus coregulate  consensus       plan
#> 102     discuss    emotion       plan       <NA>       <NA>       <NA>
#> 103        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 104        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 105   consensus       plan       plan       plan    discuss  consensus
#> 106     emotion   cohesion coregulate    discuss  consensus       plan
#> 107     discuss  consensus    emotion    discuss  consensus       plan
#> 108     discuss  synthesis  consensus       plan  consensus   cohesion
#> 109    cohesion    emotion       plan       plan    emotion    discuss
#> 110     discuss  consensus    monitor  consensus       plan       plan
#> 111   consensus coregulate      adapt       <NA>       <NA>       <NA>
#> 112     emotion   cohesion  synthesis       plan    discuss  consensus
#> 113        plan    emotion    monitor       <NA>       <NA>       <NA>
#> 114   consensus       plan       plan    emotion       plan    emotion
#> 115     discuss  consensus    emotion  consensus  consensus  synthesis
#> 116   consensus       plan       plan  consensus       plan coregulate
#> 117        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 118        plan  consensus       plan    monitor    discuss  consensus
#> 119     emotion   cohesion   cohesion   cohesion  consensus    emotion
#> 120    cohesion  consensus       <NA>       <NA>       <NA>       <NA>
#> 121        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 122   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 123     discuss  consensus    discuss       <NA>       <NA>       <NA>
#> 124        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 125        plan    monitor  consensus       plan  consensus    monitor
#> 126        plan  consensus coregulate    monitor       plan    emotion
#> 127   consensus   cohesion  consensus       plan       plan  consensus
#> 128        plan  consensus       plan       plan coregulate       plan
#> 129        plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 130     emotion  consensus    discuss  consensus    discuss  consensus
#> 131        plan  consensus       plan       <NA>       <NA>       <NA>
#> 132     monitor    discuss coregulate       plan  consensus    emotion
#> 133        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 134   consensus    discuss       plan  consensus       plan       plan
#> 135     emotion  consensus       plan    emotion       plan    discuss
#> 136     discuss  consensus    discuss coregulate    discuss  synthesis
#> 137     emotion   cohesion       plan       plan  consensus       plan
#> 138  coregulate       plan  consensus    discuss coregulate       plan
#> 139   synthesis       <NA>       <NA>       <NA>       <NA>       <NA>
#> 140        plan  consensus    discuss  synthesis    emotion  consensus
#> 141        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 142   consensus    emotion    monitor       plan       plan       plan
#> 143   consensus coregulate    emotion   cohesion    discuss    emotion
#> 144  coregulate       plan coregulate       <NA>       <NA>       <NA>
#> 145    cohesion  consensus  consensus       plan  consensus       plan
#> 146     monitor    discuss  synthesis  consensus       plan       plan
#> 147        plan       plan       plan  consensus coregulate       plan
#> 148     emotion   cohesion  consensus       plan       plan    discuss
#> 149        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 150        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 151        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 152        plan       plan  consensus       plan       plan    monitor
#> 153        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 154        plan   cohesion    emotion    discuss   cohesion  consensus
#> 155        plan       plan  consensus  consensus  synthesis    emotion
#> 156   consensus       plan    monitor       plan       plan       plan
#> 157        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 158     discuss coregulate       <NA>       <NA>       <NA>       <NA>
#> 159    cohesion  consensus       plan  consensus coregulate       plan
#> 160     discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 161     discuss  consensus       plan    emotion  consensus    emotion
#> 162   consensus       plan  consensus       plan  consensus       plan
#> 163        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 164        plan    emotion   cohesion  consensus  consensus    discuss
#> 165        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 166     discuss    discuss    discuss    emotion  consensus    discuss
#> 167     discuss    emotion       plan    emotion    discuss  consensus
#> 168     discuss  consensus    discuss    emotion   cohesion  consensus
#> 169   consensus       plan    emotion    discuss  consensus coregulate
#> 170        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 171   consensus    emotion       plan       plan    emotion   cohesion
#> 172     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 173     discuss    discuss   cohesion  consensus    discuss       <NA>
#> 174        plan       plan  consensus       <NA>       <NA>       <NA>
#> 175     monitor       plan       plan       <NA>       <NA>       <NA>
#> 176     discuss    discuss    discuss  synthesis      adapt  consensus
#> 177        plan  consensus       plan   cohesion  consensus    discuss
#> 178   consensus       plan       plan    monitor coregulate    discuss
#> 179        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 180     emotion  consensus    discuss    discuss    discuss  consensus
#> 181     monitor  consensus       plan       plan    emotion   cohesion
#> 182        plan  consensus    discuss  synthesis  consensus       plan
#> 183        plan  consensus    discuss  consensus       plan       plan
#> 184   synthesis  consensus       <NA>       <NA>       <NA>       <NA>
#> 185     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 186        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 187   synthesis  consensus       plan  consensus  consensus       <NA>
#> 188   consensus       plan    discuss  consensus    discuss  consensus
#> 189        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 190   synthesis      adapt    monitor    discuss  consensus  synthesis
#> 191     discuss   cohesion  consensus       plan   cohesion  consensus
#> 192        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 193        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 194        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 195        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 196        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 197   consensus    discuss  consensus       plan    monitor    emotion
#> 198        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 199  coregulate   cohesion       <NA>       <NA>       <NA>       <NA>
#> 200        plan       plan  consensus       plan    discuss    emotion
#> 201        plan    discuss    emotion       plan       plan    emotion
#> 202   consensus       plan    emotion       <NA>       <NA>       <NA>
#> 203        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 204   consensus    discuss    emotion       plan   cohesion  consensus
#> 205        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 206     emotion    emotion   cohesion  consensus       plan    emotion
#> 207     emotion    emotion    discuss    emotion   cohesion    monitor
#> 208        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 209     emotion  consensus       plan    discuss    discuss  consensus
#> 210        plan    emotion       plan  consensus       <NA>       <NA>
#> 211   consensus       plan  consensus       plan       plan    monitor
#> 212     emotion  consensus    emotion  consensus    discuss    discuss
#> 213        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 214   consensus    discuss    discuss coregulate  consensus   cohesion
#> 215        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 216        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 217  coregulate  consensus coregulate       <NA>       <NA>       <NA>
#> 218     emotion       plan       <NA>       <NA>       <NA>       <NA>
#> 219     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 220        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 221     emotion   cohesion coregulate    discuss   cohesion       <NA>
#> 222        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 223     emotion  consensus    emotion    emotion  consensus       plan
#> 224     emotion    discuss  consensus coregulate    emotion       plan
#> 225     emotion  consensus    monitor       plan       plan    emotion
#> 226     emotion   cohesion  consensus    emotion   cohesion  consensus
#> 227   consensus       plan       plan       plan  consensus    emotion
#> 228    cohesion coregulate       plan    monitor  consensus   cohesion
#> 229        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 230        plan       plan    monitor  consensus coregulate       plan
#> 231   consensus    discuss    monitor  consensus    discuss    emotion
#> 232        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 233   consensus    emotion   cohesion  consensus       plan    discuss
#> 234   consensus    discuss  consensus    discuss  consensus coregulate
#> 235  coregulate  consensus coregulate      adapt    discuss      adapt
#> 236     emotion   cohesion coregulate       plan    emotion    emotion
#> 237        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 238     emotion   cohesion    emotion       plan       plan       plan
#> 239     emotion  consensus    discuss  consensus coregulate       <NA>
#> 240     emotion       plan  consensus    emotion   cohesion  consensus
#> 241   consensus    emotion  consensus coregulate    discuss    emotion
#> 242     discuss    discuss  consensus    emotion   cohesion   cohesion
#> 243        plan  consensus    monitor   cohesion    emotion coregulate
#> 244   consensus  consensus       <NA>       <NA>       <NA>       <NA>
#> 245   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 246        plan  consensus       plan       plan    discuss   cohesion
#> 247   consensus    discuss  synthesis  consensus       plan    discuss
#> 248    cohesion    emotion  synthesis       plan       plan  consensus
#> 249        plan       plan  consensus    monitor    discuss  consensus
#> 250  coregulate       plan       plan  consensus       plan  consensus
#> 251    cohesion    emotion  consensus       <NA>       <NA>       <NA>
#> 252     emotion   cohesion coregulate  consensus       plan    discuss
#> 253        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 254        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 255     discuss  consensus    discuss  consensus       plan       <NA>
#> 256  coregulate    discuss  synthesis  consensus coregulate  consensus
#> 257     discuss    emotion       plan       plan  consensus       <NA>
#> 258     emotion  consensus coregulate   cohesion  consensus       plan
#> 259        plan    emotion   cohesion  consensus       <NA>       <NA>
#> 260  coregulate    emotion       plan       <NA>       <NA>       <NA>
#> 261        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 262  coregulate    discuss       plan    emotion  consensus       plan
#> 263   consensus coregulate       plan       plan       plan       plan
#> 264     emotion   cohesion  consensus  consensus    emotion       plan
#> 265        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 266     discuss coregulate    discuss  consensus    discuss  synthesis
#> 267        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 268     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 269     emotion  consensus    discuss  consensus   cohesion  consensus
#> 270     emotion    discuss    emotion  consensus       plan       plan
#> 271        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 272   consensus    discuss  consensus    emotion    discuss    discuss
#> 273     discuss  consensus    discuss coregulate    emotion  consensus
#> 274        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 275  coregulate    discuss  consensus       plan  consensus       plan
#> 276     emotion    emotion   cohesion  consensus       plan  consensus
#> 277   consensus  consensus  consensus       plan       plan coregulate
#> 278        plan  consensus       plan  consensus       plan  consensus
#> 279     discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 280        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 281     emotion    discuss  consensus       plan       plan    emotion
#> 282  coregulate    emotion   cohesion    emotion  consensus       plan
#> 283   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 284     discuss      adapt  consensus coregulate       plan       plan
#> 285        plan  consensus    monitor       plan       plan       plan
#> 286        plan  consensus       plan       <NA>       <NA>       <NA>
#> 287        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 288        plan       plan       plan       plan  consensus coregulate
#> 289  coregulate    monitor    discuss    discuss coregulate    discuss
#> 290        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 291     emotion    discuss  consensus       plan  consensus coregulate
#> 292   consensus coregulate  consensus       plan  consensus    monitor
#> 293   synthesis      adapt  consensus    discuss  consensus       plan
#> 294        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 295     discuss  consensus       plan       plan       plan    emotion
#> 296       adapt  consensus coregulate  consensus coregulate    discuss
#> 297   consensus    discuss  synthesis       plan       <NA>       <NA>
#> 298        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 299     discuss    discuss       plan  consensus       plan       plan
#> 300     discuss  consensus    emotion coregulate    monitor  consensus
#> 301     emotion   cohesion  consensus    discuss       <NA>       <NA>
#> 302    cohesion  consensus    discuss coregulate       <NA>       <NA>
#> 303  coregulate       plan       <NA>       <NA>       <NA>       <NA>
#> 304        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 305   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 306     discuss   cohesion  consensus    emotion    discuss coregulate
#> 307   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 308        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 309   consensus    discuss  consensus  consensus    monitor    discuss
#> 310        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 311        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 312   consensus       plan    discuss coregulate    discuss coregulate
#> 313     emotion    discuss  consensus       <NA>       <NA>       <NA>
#> 314        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 315        plan       plan  consensus    monitor  synthesis       plan
#> 316        plan       plan  consensus       plan    emotion    discuss
#> 317   consensus       plan    emotion       <NA>       <NA>       <NA>
#> 318        plan  consensus    discuss  synthesis  consensus    discuss
#> 319        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 320   consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 321        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 322        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 323     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 324        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 325   consensus coregulate    discuss  consensus       plan  consensus
#> 326        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 327        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 328     emotion  consensus       plan  consensus    discuss  synthesis
#> 329  coregulate    discuss  consensus       plan  consensus coregulate
#> 330        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 331        plan       plan    emotion   cohesion       plan       <NA>
#> 332   consensus coregulate       plan  consensus    discuss  synthesis
#> 333        plan       plan       plan    emotion       <NA>       <NA>
#> 334        plan       plan       plan    emotion    emotion   cohesion
#> 335   consensus coregulate      adapt  consensus       <NA>       <NA>
#> 336        plan    emotion    emotion   cohesion       plan       <NA>
#> 337        plan    emotion   cohesion  consensus       plan       plan
#> 338     emotion   cohesion    emotion   cohesion coregulate   cohesion
#> 339  coregulate    monitor       plan  consensus coregulate    discuss
#> 340  coregulate       plan       plan       plan       plan       plan
#> 341        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 342   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 343   consensus coregulate  synthesis  consensus       plan    emotion
#> 344     monitor   cohesion  consensus       plan       plan   cohesion
#> 345     discuss  consensus       plan    emotion   cohesion  consensus
#> 346     discuss coregulate       plan  consensus    monitor    discuss
#> 347    cohesion  consensus       plan  consensus       plan  consensus
#> 348        plan  consensus       plan    monitor    discuss   cohesion
#> 349        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 350        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 351     discuss   cohesion  consensus       plan  consensus       plan
#> 352        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 353   consensus    emotion    monitor    emotion  consensus  consensus
#> 354        plan    emotion  consensus    discuss    discuss  consensus
#> 355        plan       plan       plan    emotion  consensus    discuss
#> 356        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 357     discuss      adapt    emotion    discuss    emotion  consensus
#> 358     discuss coregulate    emotion  consensus       <NA>       <NA>
#> 359        plan coregulate    monitor    discuss  consensus       <NA>
#> 360     emotion       plan       <NA>       <NA>       <NA>       <NA>
#> 361     discuss    discuss    emotion  consensus       plan       plan
#> 362     emotion   cohesion    discuss    discuss  consensus       plan
#> 363     discuss    emotion   cohesion  consensus coregulate      adapt
#> 364    cohesion    discuss   cohesion  consensus       plan       plan
#> 365   consensus       plan  consensus    discuss    discuss  consensus
#> 366        plan  consensus       plan  consensus    emotion       plan
#> 367        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 368        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 369        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 370   consensus  consensus coregulate       plan    emotion   cohesion
#> 371   consensus coregulate       plan coregulate  consensus       plan
#> 372        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 373     discuss  consensus    discuss   cohesion    monitor  consensus
#> 374   consensus coregulate       plan    monitor  consensus  consensus
#> 375   consensus       plan    emotion   cohesion  consensus    discuss
#> 376        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 377  coregulate    discuss  synthesis      adapt    emotion  consensus
#> 378        plan  consensus       plan  consensus coregulate   cohesion
#> 379        plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 380        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 381     discuss    emotion    discuss    emotion  consensus       plan
#> 382    cohesion  consensus    discuss    emotion   cohesion  consensus
#> 383        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 384   consensus       plan  consensus    discuss  consensus       plan
#> 385        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 386        plan    emotion  consensus       plan    emotion   cohesion
#> 387  coregulate    monitor  consensus    discuss coregulate  consensus
#> 388        plan    monitor  consensus       plan coregulate       plan
#> 389        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 390        plan    discuss  consensus       plan    emotion    emotion
#> 391     discuss  consensus    discuss  consensus       plan  consensus
#> 392        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 393     discuss    emotion   cohesion       plan coregulate  consensus
#> 394        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 395     discuss  synthesis    emotion  consensus       <NA>       <NA>
#> 396   consensus    discuss coregulate    discuss coregulate  synthesis
#> 397        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 398     emotion  consensus coregulate  consensus  consensus    discuss
#> 399   consensus coregulate    discuss  consensus       plan       plan
#> 400     emotion    emotion  consensus       plan coregulate    discuss
#> 401        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 402    cohesion  synthesis      adapt    monitor    discuss  consensus
#> 403        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 404   consensus       plan    discuss  consensus    monitor    emotion
#> 405     discuss    discuss  consensus    discuss    monitor    discuss
#> 406        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 407     discuss    emotion   cohesion       plan    discuss coregulate
#> 408     discuss    emotion   cohesion  consensus       plan       <NA>
#> 409        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 410  coregulate       plan       plan    monitor    discuss    discuss
#> 411   consensus       plan       plan    monitor       plan       plan
#> 412        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 413        plan    emotion   cohesion    emotion    discuss    emotion
#> 414    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 415   consensus       plan       plan  consensus    discuss       <NA>
#> 416     discuss  consensus    discuss  consensus       plan  consensus
#> 417   consensus    discuss  consensus       <NA>       <NA>       <NA>
#> 418   consensus  synthesis       <NA>       <NA>       <NA>       <NA>
#> 419        plan    monitor      adapt    emotion  consensus    discuss
#> 420     monitor    discuss  synthesis  consensus coregulate    monitor
#> 421     monitor    monitor    discuss  synthesis  consensus       plan
#> 422     monitor    discuss   cohesion  consensus       plan  consensus
#> 423        plan  consensus    discuss    discuss       plan       plan
#> 424   consensus    emotion  consensus coregulate    discuss  consensus
#> 425        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 426   consensus    discuss      adapt    emotion  consensus  consensus
#> 427     monitor   cohesion coregulate    discuss  consensus coregulate
#> 428        plan  consensus       plan    emotion   cohesion  consensus
#> 429     discuss  synthesis  consensus coregulate      adapt   cohesion
#> 430   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 431   consensus       plan       plan  consensus    discuss   cohesion
#> 432     discuss    discuss  consensus    discuss  consensus coregulate
#> 433     emotion   cohesion  consensus       plan    discuss  consensus
#> 434        plan    discuss    discuss  synthesis  consensus coregulate
#> 435     discuss  synthesis      adapt   cohesion    emotion  consensus
#> 436        plan       plan  consensus    emotion    emotion   cohesion
#> 437     discuss      adapt    emotion  consensus coregulate    emotion
#> 438   consensus coregulate    emotion   cohesion    emotion       <NA>
#> 439        plan       plan  consensus    emotion   cohesion coregulate
#> 440     monitor       <NA>       <NA>       <NA>       <NA>       <NA>
#> 441        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 442        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 443     emotion    monitor   cohesion  consensus       plan       plan
#> 444        plan    emotion    discuss  consensus coregulate    discuss
#> 445        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 446        plan       plan       <NA>       <NA>       <NA>       <NA>
#> 447    cohesion  consensus  consensus coregulate       plan  consensus
#> 448        plan  consensus coregulate  consensus  consensus    discuss
#> 449        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 450     discuss      adapt    discuss  consensus coregulate       plan
#> 451        plan       plan    emotion  consensus    emotion    monitor
#> 452  coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 453        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 454        plan  consensus       plan    emotion   cohesion  consensus
#> 455     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 456    cohesion  consensus       plan       plan       plan       <NA>
#> 457        plan       plan       plan  consensus    discuss  consensus
#> 458        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 459  coregulate       plan  consensus    discuss    discuss   cohesion
#> 460        plan       plan    discuss    discuss    emotion  consensus
#> 461   consensus    emotion    emotion   cohesion  consensus       plan
#> 462     emotion  consensus       plan       plan    emotion    discuss
#> 463     emotion  consensus    discuss  consensus coregulate    discuss
#> 464        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 465        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 466        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 467   synthesis  consensus       plan  consensus    monitor       plan
#> 468   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 469  coregulate       plan    discuss    emotion   cohesion    emotion
#> 470        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 471        plan       plan       plan       plan       plan    monitor
#> 472     discuss  consensus coregulate    monitor    discuss   cohesion
#> 473        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 474     emotion       plan   cohesion   cohesion    emotion       plan
#> 475        plan       plan       plan    emotion       plan    emotion
#> 476        plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 477        plan    emotion  consensus    emotion  consensus coregulate
#> 478   consensus    discuss    emotion  consensus coregulate    discuss
#> 479     emotion    monitor       plan    emotion    discuss  synthesis
#> 480     discuss  consensus    discuss coregulate    emotion    emotion
#> 481        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 482     emotion   cohesion  consensus    emotion    monitor coregulate
#> 483        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 484        plan    emotion  consensus  consensus  consensus  consensus
#> 485        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 486        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 487     emotion   cohesion coregulate    emotion   cohesion    monitor
#> 488        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 489        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 490   consensus coregulate       plan coregulate       plan  consensus
#> 491        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 492     emotion  consensus       plan       plan    discuss  synthesis
#> 493     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 494  coregulate    emotion   cohesion  consensus   cohesion       plan
#> 495     emotion    emotion  consensus    discuss  consensus   cohesion
#> 496        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 497        plan       plan  consensus       plan    emotion   cohesion
#> 498     discuss  synthesis       <NA>       <NA>       <NA>       <NA>
#> 499        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 500   consensus  consensus coregulate    emotion    discuss  consensus
#> 501   consensus       plan    discuss  consensus       plan  consensus
#> 502  coregulate    discuss  consensus       plan    monitor       plan
#> 503   consensus       plan    emotion  consensus    discuss   cohesion
#> 504   consensus    emotion   cohesion  consensus    discuss  consensus
#> 505        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 506    cohesion  consensus coregulate    monitor       plan   cohesion
#> 507   consensus  consensus coregulate       plan    monitor    discuss
#> 508     discuss  consensus coregulate  synthesis    emotion       <NA>
#> 509     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 510        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 511     emotion    discuss    discuss  consensus    monitor       plan
#> 512   consensus       plan       plan       plan       plan       plan
#> 513        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 514        plan  consensus coregulate    discuss    discuss  consensus
#> 515        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 516   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 517  coregulate       plan  consensus    emotion  consensus    emotion
#> 518   consensus       plan  consensus       plan   cohesion       plan
#> 519        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 520   consensus    discuss    discuss    discuss    discuss    emotion
#> 521        plan    monitor    discuss  consensus  consensus       plan
#> 522        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 523    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 524        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 525        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 526        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 527    cohesion       plan       <NA>       <NA>       <NA>       <NA>
#> 528        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 529        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 530        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 531        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 532        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 533   consensus    emotion       <NA>       <NA>       <NA>       <NA>
#> 534        plan       plan       plan       plan  consensus    discuss
#> 535   consensus    emotion   cohesion  consensus coregulate       <NA>
#> 536        plan    emotion   cohesion  consensus coregulate  consensus
#> 537        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 538     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 539        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 540   consensus    discuss coregulate    emotion   cohesion  consensus
#> 541     emotion  consensus coregulate       plan  consensus       plan
#> 542        plan    emotion  consensus       plan  consensus       <NA>
#> 543     discuss    discuss  consensus       plan       plan  consensus
#> 544     discuss  consensus       plan    discuss coregulate  consensus
#> 545     emotion    discuss       plan    monitor    emotion coregulate
#> 546   consensus    discuss    emotion  consensus       plan       <NA>
#> 547        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 548        plan  consensus  consensus coregulate    emotion  consensus
#> 549     discuss coregulate    emotion  consensus       plan       plan
#> 550   consensus       plan  consensus       plan    emotion  consensus
#> 551     monitor    emotion    discuss    emotion    discuss  consensus
#> 552    cohesion       plan    emotion       <NA>       <NA>       <NA>
#> 553        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 554     discuss  consensus    emotion   cohesion  consensus    discuss
#> 555   synthesis       plan  consensus       plan       plan       plan
#> 556     discuss    discuss    monitor    discuss  consensus      adapt
#> 557     discuss    discuss  consensus    discuss      adapt  consensus
#> 558   consensus       plan    emotion    emotion       <NA>       <NA>
#> 559   consensus coregulate    discuss  synthesis       plan       <NA>
#> 560  coregulate       plan       plan  consensus       plan       plan
#> 561     emotion  consensus       plan    emotion  consensus    monitor
#> 562   consensus coregulate    discuss  consensus       plan    monitor
#> 563   consensus    discuss coregulate       plan  consensus       <NA>
#> 564        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 565        plan       plan  consensus       plan       plan  consensus
#> 566        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 567   consensus    discuss  consensus coregulate    emotion  consensus
#> 568     discuss  consensus       plan  consensus       plan       plan
#> 569        plan  consensus coregulate       plan    emotion   cohesion
#> 570        plan  consensus  consensus       plan    emotion  consensus
#> 571   consensus       plan       plan    discuss  synthesis      adapt
#> 572   consensus coregulate    discuss  consensus  consensus  consensus
#> 573     discuss  consensus coregulate    discuss    discuss    discuss
#> 574     discuss   cohesion  consensus  consensus  consensus       plan
#> 575        plan  consensus    discuss  consensus       plan    emotion
#> 576     discuss  consensus       plan    monitor    monitor   cohesion
#> 577     emotion    discuss    discuss   cohesion       <NA>       <NA>
#> 578        plan    emotion    emotion    discuss   cohesion  consensus
#> 579        plan  consensus    emotion   cohesion coregulate    monitor
#> 580   consensus  consensus       plan       plan       plan  consensus
#> 581        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 582        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 583        plan       plan    monitor  consensus coregulate    discuss
#> 584        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 585        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 586   consensus       plan  consensus    emotion       plan  consensus
#> 587  coregulate    discuss  synthesis  consensus       plan       plan
#> 588        plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 589  coregulate    emotion   cohesion       plan       <NA>       <NA>
#> 590    cohesion  consensus       plan    emotion   cohesion  consensus
#> 591     discuss    monitor       plan  consensus    emotion       <NA>
#> 592     discuss   cohesion  consensus coregulate    emotion  consensus
#> 593     emotion    discuss    monitor    emotion  consensus coregulate
#> 594     emotion    discuss    discuss    discuss  consensus    discuss
#> 595        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 596        plan       plan       <NA>       <NA>       <NA>       <NA>
#> 597        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 598     emotion  consensus       plan  consensus       plan    emotion
#> 599        plan       plan       plan       plan  consensus coregulate
#> 600        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 601        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 602     emotion  consensus coregulate    discuss  consensus coregulate
#> 603   synthesis      adapt  consensus       plan       plan       <NA>
#> 604     emotion   cohesion    emotion   cohesion  consensus    discuss
#> 605        plan  consensus    emotion       <NA>       <NA>       <NA>
#> 606     emotion   cohesion       plan       plan  consensus coregulate
#> 607   consensus coregulate    discuss  consensus  consensus    emotion
#> 608     emotion  consensus       plan  consensus    emotion    monitor
#> 609        plan    emotion  consensus    emotion  consensus       plan
#> 610        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 611     discuss    emotion   cohesion  consensus    emotion   cohesion
#> 612        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 613     emotion   cohesion  consensus    discuss  consensus       plan
#> 614        plan    monitor       plan       plan       <NA>       <NA>
#> 615   consensus       plan    emotion       plan    emotion    emotion
#> 616        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 617        plan       plan       plan  consensus       plan    monitor
#> 618        plan       plan    emotion       <NA>       <NA>       <NA>
#> 619        plan       plan    emotion   cohesion    discuss    emotion
#> 620   consensus coregulate       plan  consensus       plan       <NA>
#> 621        plan    discuss    monitor    discuss  consensus coregulate
#> 622   consensus    discuss  consensus       plan    monitor    emotion
#> 623     discuss    monitor    discuss  consensus       plan    emotion
#> 624     monitor   cohesion  consensus    emotion  consensus coregulate
#> 625        plan       plan  consensus       plan       plan    emotion
#> 626        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 627        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 628        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 629        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 630        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 631   consensus    emotion   cohesion  consensus    emotion  consensus
#> 632     discuss  consensus       plan  consensus coregulate    emotion
#> 633        plan  consensus       plan    emotion   cohesion  consensus
#> 634     discuss  synthesis  consensus coregulate       plan       plan
#> 635  coregulate      adapt   cohesion       plan  consensus       plan
#> 636        plan       plan    discuss  consensus    monitor    discuss
#> 637        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 638     emotion  consensus    discuss    emotion  consensus       plan
#> 639        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 640        plan    discuss    discuss coregulate  consensus   cohesion
#> 641        plan    discuss    emotion    monitor    discuss    emotion
#> 642   consensus       plan  consensus       plan  consensus    discuss
#> 643     discuss  consensus    discuss    emotion    emotion       <NA>
#> 644     monitor    emotion       plan  consensus    discuss  synthesis
#> 645   synthesis      adapt   cohesion  consensus       plan  consensus
#> 646     emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 647        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 648     discuss  consensus    discuss    monitor    discuss       plan
#> 649        plan       plan    emotion   cohesion  consensus  consensus
#> 650        plan       plan       <NA>       <NA>       <NA>       <NA>
#> 651        plan  consensus    discuss    discuss  consensus       plan
#> 652        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 653   consensus       plan       plan  consensus       plan    discuss
#> 654  coregulate    discuss       <NA>       <NA>       <NA>       <NA>
#> 655  coregulate       plan    monitor       plan    monitor    discuss
#> 656   consensus    discuss  consensus  consensus       <NA>       <NA>
#> 657     monitor    discuss  synthesis    emotion  consensus coregulate
#> 658     discuss coregulate    monitor    discuss   cohesion       plan
#> 659     discuss    discuss  consensus       plan       plan   cohesion
#> 660   consensus    discuss    discuss    emotion   cohesion    emotion
#> 661        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 662       adapt  consensus coregulate    discuss  consensus       plan
#> 663        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 664        plan       plan    emotion  consensus       plan  consensus
#> 665        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 666        plan       plan  consensus       <NA>       <NA>       <NA>
#> 667        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 668    cohesion  consensus    discuss coregulate    emotion    emotion
#> 669     discuss    emotion       plan       plan       plan  consensus
#> 670        plan       plan  consensus       <NA>       <NA>       <NA>
#> 671   consensus       plan  consensus       plan coregulate    discuss
#> 672        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 673     emotion       plan       plan    emotion  consensus    discuss
#> 674     discuss  consensus    discuss    emotion    discuss  synthesis
#> 675        plan  consensus    emotion   cohesion  consensus       plan
#> 676     discuss coregulate    emotion   cohesion    emotion   cohesion
#> 677   consensus       plan       plan       plan  consensus       plan
#> 678     monitor  consensus coregulate    emotion   cohesion  consensus
#> 679     discuss   cohesion    discuss  consensus       plan    emotion
#> 680     monitor    discuss    discuss    discuss  synthesis  consensus
#> 681     monitor    discuss  consensus       <NA>       <NA>       <NA>
#> 682        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 683        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 684     emotion    emotion    emotion  consensus       plan    emotion
#> 685        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 686        plan  consensus       plan    monitor    discuss    monitor
#> 687        plan coregulate    monitor       plan    discuss    discuss
#> 688   consensus       plan       plan  consensus    discuss    emotion
#> 689    cohesion  consensus    monitor    discuss  consensus    emotion
#> 690  coregulate       plan  consensus    monitor   cohesion  consensus
#> 691     discuss  consensus    discuss  consensus       plan       plan
#> 692     emotion   cohesion    discuss  consensus coregulate    discuss
#> 693        plan    emotion  consensus coregulate  consensus  consensus
#> 694     emotion  consensus coregulate       plan       plan       plan
#> 695     discuss  consensus       plan    discuss       <NA>       <NA>
#> 696     emotion    discuss  consensus    emotion   cohesion  consensus
#> 697        plan  consensus    emotion  consensus coregulate       plan
#> 698  coregulate       plan  consensus       <NA>       <NA>       <NA>
#> 699        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 700     discuss coregulate    emotion   cohesion   cohesion coregulate
#> 701        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 702     monitor    discuss coregulate  consensus    discuss  consensus
#> 703        plan    discuss  consensus    discuss  consensus    monitor
#> 704        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 705        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 706        plan       plan    emotion       <NA>       <NA>       <NA>
#> 707   consensus coregulate    discuss       <NA>       <NA>       <NA>
#> 708        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 709     emotion    monitor       plan       plan    emotion  consensus
#> 710        plan    emotion   cohesion  consensus       plan  consensus
#> 711   consensus   cohesion       <NA>       <NA>       <NA>       <NA>
#> 712        plan       plan  consensus coregulate  consensus  consensus
#> 713        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 714     discuss  consensus       plan       plan       plan  consensus
#> 715     emotion    monitor    emotion       plan       plan  consensus
#> 716   consensus coregulate    discuss coregulate  consensus coregulate
#> 717   consensus    emotion   cohesion       plan  consensus       plan
#> 718        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 719        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 720  coregulate    emotion    emotion   cohesion       plan    discuss
#> 721        plan  consensus       plan       plan    emotion  consensus
#> 722     emotion    emotion  consensus       <NA>       <NA>       <NA>
#> 723        plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 724        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 725   consensus       plan  consensus       plan       plan       plan
#> 726        plan    emotion  synthesis  consensus       <NA>       <NA>
#> 727        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 728   consensus       plan       plan  consensus       plan    emotion
#> 729        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 730        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 731     emotion       plan    monitor       plan       plan  consensus
#> 732        plan    monitor   cohesion  consensus    emotion    monitor
#> 733     discuss coregulate       plan    monitor  consensus       plan
#> 734     discuss  synthesis  consensus coregulate       plan       plan
#> 735     emotion   cohesion  consensus       plan       plan    emotion
#> 736        plan    monitor    discuss    emotion    discuss  consensus
#> 737     emotion   cohesion coregulate       plan       plan  consensus
#> 738   consensus       plan  consensus       plan  consensus    discuss
#> 739        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 740  coregulate   cohesion  consensus    discuss  consensus coregulate
#> 741        plan    emotion   cohesion  consensus       plan    emotion
#> 742     discuss  consensus coregulate       plan coregulate    discuss
#> 743        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 744   synthesis  consensus    emotion  consensus       plan  consensus
#> 745     discuss    discuss  consensus       plan       plan       plan
#> 746   consensus       plan  consensus    discuss      adapt   cohesion
#> 747        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 748        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 749     discuss    discuss  consensus coregulate  consensus    emotion
#> 750       adapt  consensus    discuss  consensus       plan       plan
#> 751   consensus  consensus       plan    discuss  consensus       plan
#> 752        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 753        plan  consensus    emotion   cohesion    discuss  synthesis
#> 754        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 755   synthesis  consensus       plan  consensus  consensus    emotion
#> 756        plan    emotion    emotion    discuss      adapt   cohesion
#> 757     emotion    discuss  consensus coregulate       plan       plan
#> 758   consensus coregulate  consensus       plan    monitor    discuss
#> 759        plan    monitor coregulate       plan coregulate       plan
#> 760   consensus    discuss    emotion       <NA>       <NA>       <NA>
#> 761        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 762        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 763        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 764        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 765        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 766  coregulate coregulate    discuss    discuss       <NA>       <NA>
#> 767    cohesion   cohesion  consensus       plan    emotion   cohesion
#> 768    cohesion  consensus    discuss    discuss  synthesis       plan
#> 769        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 770  coregulate    monitor  consensus       plan    emotion  consensus
#> 771   synthesis       <NA>       <NA>       <NA>       <NA>       <NA>
#> 772        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 773     discuss  consensus       plan  consensus       <NA>       <NA>
#> 774  coregulate    discuss    monitor  consensus coregulate  consensus
#> 775        plan    emotion       plan       plan       plan  consensus
#> 776        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 777        plan    discuss  consensus coregulate       <NA>       <NA>
#> 778   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 779   consensus    discuss    monitor    emotion       <NA>       <NA>
#> 780        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 781    cohesion    discuss  consensus coregulate       <NA>       <NA>
#> 782     discuss    discuss  consensus coregulate    emotion       plan
#> 783     monitor       plan       plan    emotion    discuss       <NA>
#> 784  coregulate    monitor       plan       plan  consensus coregulate
#> 785   consensus    discuss    discuss  consensus       plan       plan
#> 786        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 787        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 788        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 789     monitor    discuss    emotion  consensus  consensus       plan
#> 790        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 791   consensus       plan    discuss       <NA>       <NA>       <NA>
#> 792     monitor coregulate    discuss  consensus    emotion  consensus
#> 793   consensus    emotion  consensus       plan       plan       plan
#> 794        plan       plan       plan    monitor    discuss    discuss
#> 795        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 796        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 797        plan  consensus       plan       plan       plan  consensus
#> 798        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 799   consensus  consensus       plan       plan  consensus  consensus
#> 800   consensus       plan  consensus    discuss    emotion  consensus
#> 801        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 802        plan  consensus       plan  consensus       plan    emotion
#> 803        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 804   synthesis       plan    emotion   cohesion  consensus       plan
#> 805   consensus coregulate  consensus       plan       plan       plan
#> 806  coregulate    emotion   cohesion  consensus coregulate       plan
#> 807     emotion   cohesion       plan       <NA>       <NA>       <NA>
#> 808    cohesion coregulate       <NA>       <NA>       <NA>       <NA>
#> 809        plan       plan  consensus    discuss  synthesis       plan
#> 810        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 811        plan    emotion  consensus       plan    emotion    emotion
#> 812    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 813     emotion  consensus    discuss    monitor       plan       plan
#> 814        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 815     discuss   cohesion coregulate       plan       <NA>       <NA>
#> 816        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 817   consensus       plan    emotion    discuss  consensus coregulate
#> 818  coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 819     monitor       plan  consensus       plan    emotion   cohesion
#> 820     discuss    discuss       <NA>       <NA>       <NA>       <NA>
#> 821     emotion    monitor    discuss  synthesis       plan    monitor
#> 822        plan   cohesion    emotion    discuss   cohesion  consensus
#> 823    cohesion  consensus       plan  consensus coregulate       plan
#> 824    cohesion    discuss  consensus       plan  consensus       plan
#> 825  coregulate    monitor    discuss    discuss coregulate       plan
#> 826        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 827     discuss  consensus  consensus    discuss  consensus       plan
#> 828     discuss  consensus    discuss  consensus  consensus   cohesion
#> 829     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 830   consensus coregulate    discuss  consensus       plan       plan
#> 831     discuss    discuss  synthesis       plan    emotion  consensus
#> 832     emotion  consensus       plan  consensus    discuss  synthesis
#> 833    cohesion  synthesis  consensus       plan    monitor    discuss
#> 834   consensus    monitor       plan       plan    emotion  consensus
#> 835     discuss    emotion   cohesion  consensus       plan coregulate
#> 836        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 837     emotion  consensus  consensus  synthesis    emotion   cohesion
#> 838        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 839        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 840  coregulate    discuss    discuss  consensus    discuss       <NA>
#> 841        plan    monitor       plan       <NA>       <NA>       <NA>
#> 842        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 843        plan  consensus coregulate    discuss  consensus    discuss
#> 844        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 845        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 846     discuss  consensus    discuss    discuss  consensus    discuss
#> 847        plan  consensus    emotion   cohesion    emotion  consensus
#> 848     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 849   consensus    discuss    discuss  consensus       plan       plan
#> 850     monitor  consensus    emotion    emotion   cohesion  consensus
#> 851        plan   cohesion  consensus       plan    emotion       plan
#> 852        plan    discuss  consensus       plan  consensus       plan
#> 853   consensus       plan    emotion   cohesion    emotion    discuss
#> 854        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 855        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 856   consensus       plan  consensus       <NA>       <NA>       <NA>
#> 857   consensus    monitor       plan       plan    monitor    discuss
#> 858  coregulate  consensus    emotion  consensus    discuss       <NA>
#> 859        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 860        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 861        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 862     emotion       plan    discuss  consensus       plan    discuss
#> 863   consensus       plan       plan    discuss    emotion   cohesion
#> 864        plan       plan       plan       <NA>       <NA>       <NA>
#> 865   consensus    discuss  consensus  consensus coregulate   cohesion
#> 866        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 867   consensus coregulate    discuss      adapt  consensus       plan
#> 868        plan       plan  consensus    monitor  consensus    monitor
#> 869     discuss    discuss  consensus       plan  consensus    emotion
#> 870    cohesion       plan    emotion  consensus coregulate  consensus
#> 871        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 872        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 873   consensus       plan  consensus    emotion    discuss  synthesis
#> 874        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 875        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 876     emotion    emotion    emotion  consensus       plan    monitor
#> 877   consensus       plan  consensus coregulate       plan  consensus
#> 878    cohesion    emotion  consensus coregulate    emotion    discuss
#> 879     discuss   cohesion coregulate       <NA>       <NA>       <NA>
#> 880        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 881   consensus       plan       plan       plan  consensus       <NA>
#> 882     discuss    emotion   cohesion       <NA>       <NA>       <NA>
#> 883        plan  consensus    discuss coregulate    discuss  consensus
#> 884        plan    monitor    discuss  consensus       plan       plan
#> 885  coregulate    emotion   cohesion  consensus    discuss  synthesis
#> 886   consensus       plan  consensus coregulate    discuss  consensus
#> 887        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 888        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 889        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 890        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 891        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 892   consensus       plan  consensus    monitor       plan   cohesion
#> 893   consensus    discuss    emotion       <NA>       <NA>       <NA>
#> 894        plan    monitor       plan       plan       plan       plan
#> 895        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 896        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 897        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 898     emotion       plan    discuss  consensus    discuss  consensus
#> 899        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 900   consensus    emotion  consensus       plan  consensus       plan
#> 901        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 902        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 903    cohesion  consensus    emotion  consensus    discuss  consensus
#> 904        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 905     monitor    discuss       <NA>       <NA>       <NA>       <NA>
#> 906   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 907        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 908        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 909   consensus       plan  consensus       plan    emotion  consensus
#> 910     discuss  consensus       plan coregulate    emotion       plan
#> 911     emotion    discuss  synthesis  consensus    discuss  consensus
#> 912  coregulate  consensus    emotion       <NA>       <NA>       <NA>
#> 913        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 914     emotion    discuss    discuss  consensus    discuss      adapt
#> 915        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 916     discuss  consensus    monitor  consensus       plan  consensus
#> 917     monitor    discuss  consensus       plan    monitor    emotion
#> 918        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 919    cohesion  consensus       plan       plan  consensus coregulate
#> 920        plan       plan  consensus    emotion       plan       plan
#> 921  coregulate    discuss  consensus    emotion       <NA>       <NA>
#> 922     discuss  consensus  consensus coregulate    emotion  consensus
#> 923   consensus coregulate    emotion  consensus coregulate       plan
#> 924        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 925   consensus    emotion   cohesion  consensus coregulate    discuss
#> 926     emotion   cohesion    discuss  consensus       plan       plan
#> 927        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 928        plan       plan    emotion  consensus       plan    monitor
#> 929     discuss  consensus    emotion  consensus    emotion coregulate
#> 930        plan  consensus  consensus coregulate       plan  consensus
#> 931     discuss    emotion   cohesion       plan       plan  consensus
#> 932       adapt   cohesion  consensus       <NA>       <NA>       <NA>
#> 933   consensus  consensus       plan       plan    monitor       plan
#> 934   consensus coregulate    discuss  consensus   cohesion       <NA>
#> 935  coregulate   cohesion       plan    emotion   cohesion  consensus
#> 936     emotion  consensus coregulate       plan    monitor coregulate
#> 937        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 938        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 939     discuss  consensus    monitor    discuss  consensus coregulate
#> 940     emotion  consensus       plan  consensus  consensus    emotion
#> 941     discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 942        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 943    cohesion  consensus       <NA>       <NA>       <NA>       <NA>
#> 944  coregulate  synthesis       <NA>       <NA>       <NA>       <NA>
#> 945  coregulate    emotion   cohesion       <NA>       <NA>       <NA>
#> 946        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 947        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 948     monitor   cohesion coregulate       plan    emotion   cohesion
#> 949     discuss    discuss    discuss    discuss    emotion    discuss
#> 950        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 951     emotion   cohesion  consensus coregulate  consensus    discuss
#> 952        plan  consensus       plan    monitor  consensus    discuss
#> 953   consensus    discuss    discuss coregulate  consensus       plan
#> 954   consensus  consensus       plan  consensus    discuss  consensus
#> 955   consensus    discuss coregulate       plan       plan    emotion
#> 956     discuss  consensus       plan  consensus       plan    discuss
#> 957        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 958        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 959   consensus       plan  consensus       <NA>       <NA>       <NA>
#> 960        plan    discuss  consensus       plan  consensus coregulate
#> 961        plan  consensus coregulate       plan    monitor    emotion
#> 962  coregulate    monitor   cohesion    emotion   cohesion       plan
#> 963   consensus       plan    emotion    monitor       plan    emotion
#> 964        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 965        plan  consensus       plan       plan       <NA>       <NA>
#> 966     emotion   cohesion    monitor       plan  synthesis  consensus
#> 967    cohesion    emotion  consensus       plan  consensus    discuss
#> 968        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 969     discuss coregulate       <NA>       <NA>       <NA>       <NA>
#> 970   consensus    monitor  consensus    discuss    discuss    emotion
#> 971     emotion  consensus coregulate    discuss  synthesis coregulate
#> 972        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 973  coregulate       plan    discuss      adapt   cohesion  consensus
#> 974     monitor  consensus    discuss  consensus    emotion  consensus
#> 975        plan  consensus       plan       plan  consensus       plan
#> 976     monitor  consensus coregulate    monitor       plan       plan
#> 977   consensus  consensus    discuss  consensus   cohesion  consensus
#> 978        plan  consensus coregulate    discuss  consensus       plan
#> 979        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 980        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 981        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 982        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 983     emotion       plan  consensus    discuss    discuss  consensus
#> 984        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 985     discuss    discuss  consensus  consensus       plan  consensus
#> 986     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 987        plan       plan       plan    monitor    emotion   cohesion
#> 988     emotion  consensus    discuss    discuss  consensus coregulate
#> 989   consensus  synthesis  consensus  consensus coregulate    emotion
#> 990        plan       plan  consensus    emotion  consensus coregulate
#> 991        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 992   consensus    discuss  consensus       plan       <NA>       <NA>
#> 993        plan  consensus coregulate    emotion    monitor    discuss
#> 994     monitor  consensus    discuss    discuss       <NA>       <NA>
#> 995        plan    emotion    emotion  consensus       plan    monitor
#> 996   consensus       plan    monitor    discuss    discuss    discuss
#> 997   consensus       plan       plan    emotion    emotion   cohesion
#> 998        plan       plan       plan       plan    emotion    discuss
#> 999        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1000  consensus       plan       plan       <NA>       <NA>       <NA>
#> 1001    discuss  consensus       plan       plan   cohesion  consensus
#> 1002 coregulate    discuss    discuss      adapt       <NA>       <NA>
#> 1003       plan  consensus coregulate  consensus    discuss    discuss
#> 1004  synthesis  consensus    discuss  synthesis      adapt  consensus
#> 1005  consensus    monitor    discuss  synthesis coregulate  consensus
#> 1006    emotion       plan       <NA>       <NA>       <NA>       <NA>
#> 1007       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1008       plan       plan    discuss    discuss      adapt  consensus
#> 1009    discuss    discuss      adapt    discuss    discuss    emotion
#> 1010       plan  consensus       plan   cohesion coregulate    discuss
#> 1011       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1012       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1013       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1014  consensus    discuss    monitor       plan  consensus coregulate
#> 1015       plan       plan       plan  consensus    discuss    emotion
#> 1016  consensus       plan       plan       plan  consensus coregulate
#> 1017       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1018       plan    emotion   cohesion    discuss       <NA>       <NA>
#> 1019    discuss      adapt    emotion    emotion       <NA>       <NA>
#> 1020       plan    emotion   cohesion  consensus  consensus coregulate
#> 1021   cohesion  consensus       plan       plan       <NA>       <NA>
#> 1022       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1023    monitor    discuss  synthesis coregulate       plan    emotion
#> 1024  consensus       plan    discuss  consensus coregulate    emotion
#> 1025       plan    discuss       <NA>       <NA>       <NA>       <NA>
#> 1026       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1027  consensus coregulate    discuss      adapt  consensus       plan
#> 1028    discuss   cohesion    discuss       <NA>       <NA>       <NA>
#> 1029       plan       plan  consensus coregulate       <NA>       <NA>
#> 1030  consensus coregulate       plan  consensus coregulate       plan
#> 1031    discuss    discuss    discuss    emotion  consensus    monitor
#> 1032    monitor    discuss    discuss    discuss  synthesis coregulate
#> 1033       plan       plan       plan       plan    discuss  synthesis
#> 1034    emotion   cohesion       plan  consensus       plan    emotion
#> 1035    discuss      adapt       plan       plan    monitor    emotion
#> 1036       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1037       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1038    discuss      adapt    discuss  synthesis  consensus       plan
#> 1039    emotion    discuss    discuss  synthesis      adapt   cohesion
#> 1040       plan    discuss       <NA>       <NA>       <NA>       <NA>
#> 1041    emotion   cohesion  consensus    emotion       <NA>       <NA>
#> 1042  consensus       plan       plan       plan       plan    discuss
#> 1043  consensus coregulate    discuss  consensus       <NA>       <NA>
#> 1044    emotion    monitor    emotion  consensus    monitor  consensus
#> 1045       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1046    monitor       plan       plan       plan       plan    discuss
#> 1047       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1048       plan       plan    monitor  consensus       <NA>       <NA>
#> 1049  consensus       plan  consensus coregulate       <NA>       <NA>
#> 1050  synthesis  consensus    monitor    discuss       <NA>       <NA>
#> 1051       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1052    emotion   cohesion    monitor    discuss  consensus    monitor
#> 1053       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1054       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1055       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1056    monitor    discuss  synthesis      adapt  consensus    monitor
#> 1057       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1058    discuss    discuss       <NA>       <NA>       <NA>       <NA>
#> 1059   cohesion coregulate coregulate  consensus    emotion  consensus
#> 1060       plan    emotion   cohesion    emotion   cohesion  consensus
#> 1061       plan  consensus       plan       <NA>       <NA>       <NA>
#> 1062       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1063       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1064    discuss    emotion       <NA>       <NA>       <NA>       <NA>
#> 1065  consensus coregulate    discuss      adapt       <NA>       <NA>
#> 1066       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1067       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1068       plan    discuss      adapt  consensus  consensus coregulate
#> 1069       plan       plan  consensus    emotion    discuss    emotion
#> 1070       plan  consensus  consensus  consensus       <NA>       <NA>
#> 1071       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1072  synthesis  consensus       plan  consensus    monitor       plan
#> 1073    discuss    discuss    discuss  consensus       plan       plan
#> 1074       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1075  consensus    monitor    discuss  consensus       plan  consensus
#> 1076  consensus    emotion  consensus       <NA>       <NA>       <NA>
#> 1077   cohesion    monitor    discuss    discuss    discuss  synthesis
#> 1078    discuss  consensus       plan       <NA>       <NA>       <NA>
#> 1079  consensus       plan    monitor       plan    emotion   cohesion
#> 1080       plan       plan       plan       plan    monitor  consensus
#> 1081  consensus  consensus coregulate    discuss  consensus       plan
#> 1082  consensus       plan    emotion    emotion   cohesion  consensus
#> 1083       plan    emotion  consensus       plan  consensus   cohesion
#> 1084       plan    discuss   cohesion   cohesion  consensus       plan
#> 1085       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1086    emotion    emotion coregulate       plan    monitor    discuss
#> 1087       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1088    discuss    emotion   cohesion    monitor    emotion    monitor
#> 1089    emotion       plan       plan       plan    monitor  consensus
#> 1090    discuss    discuss coregulate    emotion   cohesion  consensus
#> 1091    discuss  consensus   cohesion  consensus coregulate  consensus
#> 1092 coregulate    discuss coregulate       plan    discuss      adapt
#> 1093  consensus       plan    discuss    discuss coregulate    monitor
#> 1094  consensus  consensus   cohesion coregulate  consensus       plan
#> 1095  consensus       plan    monitor  consensus       plan    emotion
#> 1096       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1097       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1098       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1099       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1100   cohesion    discuss      adapt       plan  consensus       plan
#> 1101       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1102       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1103   cohesion  consensus    emotion   cohesion  consensus       plan
#> 1104  consensus       plan  consensus    discuss  synthesis      adapt
#> 1105       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1106  consensus       plan       plan  consensus    discuss  consensus
#> 1107    discuss    discuss   cohesion  consensus       plan    emotion
#> 1108  consensus       plan       plan       plan       <NA>       <NA>
#> 1109    monitor coregulate    discuss  synthesis  consensus    discuss
#> 1110       plan       plan       plan       plan  consensus       plan
#> 1111    discuss      adapt  consensus       plan  consensus coregulate
#> 1112       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1113    emotion   cohesion coregulate       <NA>       <NA>       <NA>
#> 1114       plan    emotion    discuss    emotion  consensus       <NA>
#> 1115    emotion  consensus  consensus       plan  consensus    monitor
#> 1116       plan  consensus       plan    emotion   cohesion coregulate
#> 1117  consensus      adapt  consensus       plan    monitor coregulate
#> 1118       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1119  consensus       plan       plan       <NA>       <NA>       <NA>
#> 1120    monitor    emotion    emotion  consensus    monitor    discuss
#> 1121      adapt   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1122 coregulate    monitor    emotion       <NA>       <NA>       <NA>
#> 1123       plan   cohesion       plan       plan  consensus       plan
#> 1124  consensus       plan       plan       plan    discuss    emotion
#> 1125       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1126      adapt  consensus    monitor    discuss    discuss coregulate
#> 1127    discuss coregulate  consensus    monitor       <NA>       <NA>
#> 1128 coregulate    discuss    monitor coregulate       plan       plan
#> 1129  synthesis      adapt    emotion  consensus       plan  consensus
#> 1130  consensus       plan  consensus coregulate  consensus coregulate
#> 1131       plan       plan  consensus    discuss    emotion    discuss
#> 1132  consensus       plan       plan    monitor  consensus       <NA>
#> 1133  consensus    monitor       plan       plan    discuss    discuss
#> 1134       plan    emotion   cohesion  consensus  consensus       plan
#> 1135  consensus    emotion    emotion    emotion   cohesion    emotion
#> 1136  consensus       plan       plan   cohesion  consensus       plan
#> 1137  consensus    emotion    monitor    discuss  consensus       plan
#> 1138  consensus    discuss    discuss coregulate  consensus coregulate
#> 1139 coregulate    monitor    discuss  consensus    discuss  synthesis
#> 1140 coregulate    discuss      adapt    emotion   cohesion       plan
#> 1141       plan    emotion       plan    emotion  consensus coregulate
#> 1142  synthesis       plan       plan       plan       plan  consensus
#> 1143       plan       plan       plan       plan       plan    discuss
#> 1144    discuss    discuss       plan    monitor       plan       plan
#> 1145    emotion   cohesion       plan    monitor       <NA>       <NA>
#> 1146       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1147  consensus    discuss coregulate       plan       plan  consensus
#> 1148    emotion   cohesion coregulate    emotion   cohesion  consensus
#> 1149  consensus coregulate    discuss  consensus       plan    monitor
#> 1150    emotion   cohesion  consensus coregulate    monitor  consensus
#> 1151       plan coregulate   cohesion  consensus       plan    emotion
#> 1152  consensus       plan  consensus    emotion  consensus  consensus
#> 1153 coregulate  synthesis      adapt   cohesion       plan  consensus
#> 1154    emotion   cohesion  consensus       plan       plan  consensus
#> 1155      adapt  consensus       plan       plan    monitor    discuss
#> 1156       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1157       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1158    discuss      adapt    emotion  consensus       plan  consensus
#> 1159       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1160       plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 1161       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1162 coregulate  consensus       plan       plan       <NA>       <NA>
#> 1163    discuss    discuss  consensus       plan  consensus       plan
#> 1164       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1165    emotion   cohesion  consensus coregulate    discuss    discuss
#> 1166  consensus    discuss    discuss    emotion  consensus      adapt
#> 1167    discuss      adapt  consensus       plan       plan    emotion
#> 1168    discuss coregulate   cohesion       <NA>       <NA>       <NA>
#> 1169    emotion   cohesion  consensus       <NA>       <NA>       <NA>
#> 1170   cohesion  consensus       plan    emotion   cohesion       plan
#> 1171    discuss      adapt   cohesion       plan       plan    discuss
#> 1172       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1173    discuss  synthesis      adapt   cohesion  consensus       plan
#> 1174       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1175    discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 1176       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1177       plan       plan    emotion    discuss   cohesion    discuss
#> 1178       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1179    discuss  synthesis    discuss  consensus coregulate    discuss
#> 1180  consensus coregulate  synthesis  consensus       plan       plan
#> 1181   cohesion  consensus       plan   cohesion  consensus       plan
#> 1182       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1183       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1184    emotion  consensus coregulate    emotion       plan       plan
#> 1185  synthesis      adapt   cohesion coregulate       plan    monitor
#> 1186   cohesion  consensus coregulate    discuss coregulate    emotion
#> 1187 coregulate       plan       plan    emotion  consensus       plan
#> 1188       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1189    emotion  consensus       plan    emotion   cohesion  consensus
#> 1190    monitor    discuss    emotion    discuss       <NA>       <NA>
#> 1191  consensus    monitor       plan  consensus       plan  consensus
#> 1192  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1193  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1194       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1195  consensus    discuss  synthesis  consensus       plan  consensus
#> 1196       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1197  consensus coregulate       plan       plan    emotion  consensus
#> 1198       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1199    monitor    discuss    emotion  consensus coregulate    discuss
#> 1200       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1201       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1202       plan   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1203    discuss    discuss  consensus       plan    emotion  consensus
#> 1204   cohesion    emotion   cohesion  consensus       plan    discuss
#> 1205 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1206 coregulate       plan       <NA>       <NA>       <NA>       <NA>
#> 1207    monitor    discuss  synthesis  consensus coregulate  consensus
#> 1208   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1209 coregulate       plan  consensus    emotion    emotion    emotion
#> 1210    emotion  consensus       plan       plan       plan  consensus
#> 1211    monitor    emotion       <NA>       <NA>       <NA>       <NA>
#> 1212    discuss  consensus       plan       plan       <NA>       <NA>
#> 1213    emotion   cohesion coregulate    discuss    discuss    emotion
#> 1214       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1215  consensus       plan  consensus    discuss  synthesis  consensus
#> 1216       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1217       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1218 coregulate       plan  consensus coregulate       plan       plan
#> 1219  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1220       plan       plan       plan       plan       plan       plan
#> 1221  consensus    discuss  consensus       <NA>       <NA>       <NA>
#> 1222       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1223       plan       plan  consensus       plan    monitor       plan
#> 1224    discuss    discuss   cohesion       <NA>       <NA>       <NA>
#> 1225       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1226  consensus coregulate    monitor       plan  consensus    monitor
#> 1227       plan    discuss  consensus       plan       plan       plan
#> 1228       plan       plan  consensus       plan  consensus    monitor
#> 1229 coregulate       plan       plan       plan  consensus    monitor
#> 1230       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1231       plan    emotion   cohesion  consensus       plan    emotion
#> 1232       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1233    monitor    emotion       plan coregulate       plan  consensus
#> 1234       plan       plan       plan    emotion    monitor    emotion
#> 1235      adapt    emotion    emotion   cohesion       plan    discuss
#> 1236 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1237       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1238   cohesion coregulate    discuss    discuss  consensus    emotion
#> 1239  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1240  consensus       plan       plan       plan  consensus       plan
#> 1241       plan       plan    emotion    discuss  synthesis  consensus
#> 1242  consensus       plan    monitor       plan    monitor    monitor
#> 1243    emotion   cohesion    emotion  consensus    discuss      adapt
#> 1244   cohesion coregulate    emotion coregulate    discuss   cohesion
#> 1245      adapt   cohesion  consensus      adapt  consensus       plan
#> 1246  consensus       plan       plan  consensus  consensus  consensus
#> 1247   cohesion  consensus coregulate coregulate       plan    monitor
#> 1248    emotion   cohesion    emotion       plan  consensus       plan
#> 1249       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1250    monitor    discuss    discuss  synthesis       <NA>       <NA>
#> 1251  consensus       plan       plan    monitor  consensus coregulate
#> 1252 coregulate       plan       plan       plan  consensus       plan
#> 1253    monitor  consensus coregulate       <NA>       <NA>       <NA>
#> 1254       plan       plan       plan    emotion       plan  consensus
#> 1255    discuss    emotion  consensus       plan  consensus coregulate
#> 1256       plan       plan    discuss      adapt  consensus       plan
#> 1257  synthesis  consensus    discuss  consensus       plan       <NA>
#> 1258       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1259       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1260    discuss  synthesis  consensus       plan       plan       plan
#> 1261    emotion   cohesion coregulate       <NA>       <NA>       <NA>
#> 1262       plan    monitor    discuss       <NA>       <NA>       <NA>
#> 1263       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1264  consensus       plan    emotion  consensus       plan    discuss
#> 1265       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1266       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1267       plan  consensus coregulate  consensus coregulate   cohesion
#> 1268    monitor       plan       plan  consensus coregulate  consensus
#> 1269  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1270 coregulate    discuss coregulate       plan  consensus       plan
#> 1271  consensus coregulate    emotion  consensus       plan  consensus
#> 1272  synthesis      adapt   cohesion  consensus    emotion       plan
#> 1273       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1274  consensus coregulate  consensus       plan       plan  consensus
#> 1275       plan       plan    emotion  consensus       plan    emotion
#> 1276    discuss    emotion    monitor    discuss    discuss  consensus
#> 1277       plan       plan  consensus       <NA>       <NA>       <NA>
#> 1278       plan  consensus       plan    emotion  consensus    emotion
#> 1279       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1280       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1281    discuss coregulate       plan  consensus    discuss    emotion
#> 1282       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1283    discuss    discuss      adapt       <NA>       <NA>       <NA>
#> 1284    emotion       plan  consensus   cohesion       <NA>       <NA>
#> 1285       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1286  consensus       plan       plan  consensus coregulate       plan
#> 1287 coregulate    discuss  synthesis   cohesion  consensus       plan
#> 1288    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1289      adapt   cohesion       plan    emotion    emotion  consensus
#> 1290       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1291       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1292       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1293      adapt  consensus coregulate    emotion  consensus coregulate
#> 1294    discuss  synthesis  consensus    monitor    discuss    discuss
#> 1295  consensus       plan    discuss  synthesis    discuss    discuss
#> 1296       plan       plan       plan  consensus       plan  consensus
#> 1297  consensus  consensus       plan  consensus  consensus    discuss
#> 1298       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1299 coregulate    discuss  synthesis   cohesion coregulate    emotion
#> 1300   cohesion  consensus coregulate    discuss  consensus coregulate
#> 1301       plan       plan       plan       <NA>       <NA>       <NA>
#> 1302  consensus       plan coregulate       plan  consensus       plan
#> 1303    emotion   cohesion   cohesion  consensus       plan       <NA>
#> 1304  consensus    emotion coregulate       plan  consensus coregulate
#> 1305  consensus  consensus       plan    monitor       plan  consensus
#> 1306  consensus    discuss  consensus       plan    discuss coregulate
#> 1307       plan  consensus       plan       plan       plan       plan
#> 1308 coregulate    discuss       <NA>       <NA>       <NA>       <NA>
#> 1309       plan    discuss    discuss       <NA>       <NA>       <NA>
#> 1310    monitor       plan    discuss    discuss    discuss coregulate
#> 1311   cohesion    monitor    emotion  consensus       <NA>       <NA>
#> 1312       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1313       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1314       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1315    monitor    discuss    monitor    discuss  synthesis  consensus
#> 1316       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1317       plan    monitor coregulate    emotion       <NA>       <NA>
#> 1318       plan  consensus coregulate       <NA>       <NA>       <NA>
#> 1319 coregulate  consensus       plan  consensus  consensus coregulate
#> 1320 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1321   cohesion  consensus       plan    monitor  consensus       plan
#> 1322  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1323       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1324  consensus    discuss   cohesion  consensus coregulate    discuss
#> 1325       plan  consensus    discuss  consensus       <NA>       <NA>
#> 1326       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1327       plan coregulate  consensus  consensus       plan    discuss
#> 1328  consensus  synthesis  consensus  consensus       plan  consensus
#> 1329 coregulate       plan       plan  consensus  consensus       plan
#> 1330       plan   cohesion       plan  consensus coregulate    emotion
#> 1331       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1332  consensus    monitor       plan    discuss  consensus       <NA>
#> 1333    emotion       plan    discuss    emotion       <NA>       <NA>
#> 1334  consensus       plan       plan    emotion    monitor coregulate
#> 1335       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1336  synthesis  consensus    monitor  consensus       plan  consensus
#> 1337       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1338 coregulate  consensus    emotion  consensus       <NA>       <NA>
#> 1339 coregulate    emotion       <NA>       <NA>       <NA>       <NA>
#> 1340       plan       plan  consensus  synthesis      adapt    monitor
#> 1341 coregulate coregulate       plan       plan  consensus       <NA>
#> 1342       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1343    discuss  synthesis      adapt   cohesion  consensus coregulate
#> 1344    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1345  consensus       plan       plan  consensus       plan  consensus
#> 1346       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1347       plan  consensus       plan    monitor coregulate    discuss
#> 1348      adapt    emotion       plan    monitor    discuss coregulate
#> 1349       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1350       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1351    emotion  consensus coregulate      adapt    emotion       plan
#> 1352       plan       plan       plan  consensus    monitor    discuss
#> 1353       plan    emotion   cohesion  consensus       <NA>       <NA>
#> 1354    monitor    emotion    emotion       plan  consensus    discuss
#> 1355       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1356    discuss  synthesis      adapt  consensus coregulate    discuss
#> 1357  consensus       plan  consensus    discuss    discuss  consensus
#> 1358       plan       plan       plan       plan       plan  consensus
#> 1359    monitor   cohesion    monitor coregulate    emotion    emotion
#> 1360       plan       plan       plan  consensus       plan  consensus
#> 1361  consensus    monitor       plan       plan       plan       plan
#> 1362    emotion   cohesion  consensus       plan       plan    discuss
#> 1363   cohesion    discuss    emotion  consensus coregulate  consensus
#> 1364       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1365       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1366       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1367 coregulate  consensus    emotion  consensus    discuss  consensus
#> 1368       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1369    discuss      adapt  consensus  consensus  consensus coregulate
#> 1370       plan    discuss  consensus       plan       plan       plan
#> 1371 coregulate    monitor   cohesion  consensus    discuss  consensus
#> 1372       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1373       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1374       plan  consensus       plan       plan       plan       plan
#> 1375       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1376    discuss  synthesis  consensus       plan       plan  consensus
#> 1377       plan  consensus       plan coregulate       plan  consensus
#> 1378       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1379       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1380  consensus       plan  consensus coregulate    discuss      adapt
#> 1381       plan    emotion    emotion  consensus       plan  consensus
#> 1382       plan       plan  consensus       plan       plan  consensus
#> 1383    discuss      adapt  consensus       <NA>       <NA>       <NA>
#> 1384  consensus    discuss    discuss    emotion   cohesion    emotion
#> 1385  consensus       plan       plan   cohesion coregulate    emotion
#> 1386       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1387       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1388       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1389       plan  consensus coregulate    emotion       plan    discuss
#> 1390       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1391 coregulate  consensus  consensus    emotion  consensus    emotion
#> 1392       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1393  consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 1394    discuss  consensus  consensus       plan       plan  consensus
#> 1395    discuss  consensus coregulate    discuss    discuss  synthesis
#> 1396       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1397       plan    emotion coregulate    discuss  consensus       <NA>
#> 1398    emotion   cohesion    discuss    emotion   cohesion  consensus
#> 1399    emotion    emotion  consensus       plan    emotion   cohesion
#> 1400       plan       plan       plan  consensus    discuss    discuss
#> 1401 coregulate  consensus  consensus       plan  consensus       plan
#> 1402       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1403       plan  consensus    discuss  synthesis      adapt    emotion
#> 1404  consensus       plan  consensus       plan       plan       plan
#> 1405       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1406       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1407  consensus coregulate       plan       plan       plan  consensus
#> 1408       plan  consensus coregulate  consensus    emotion  consensus
#> 1409    emotion  consensus coregulate    emotion   cohesion  consensus
#> 1410       plan       plan       plan  consensus       plan       plan
#> 1411       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1412       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1413  consensus    discuss  consensus       plan       plan  consensus
#> 1414       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1415    monitor    emotion   cohesion  consensus      adapt   cohesion
#> 1416       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1417  consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 1418    monitor       plan       <NA>       <NA>       <NA>       <NA>
#> 1419    monitor    discuss coregulate    monitor    discuss    emotion
#> 1420       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1421       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1422       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1423   cohesion  consensus coregulate       plan    discuss    discuss
#> 1424  consensus       plan       plan       plan       plan       plan
#> 1425 coregulate       plan       plan    monitor       plan  consensus
#> 1426    discuss coregulate    discuss    discuss    discuss  synthesis
#> 1427       plan       plan       plan  consensus       plan       plan
#> 1428       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1429  consensus       plan       plan  consensus       plan  consensus
#> 1430       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1431       plan       plan  consensus       <NA>       <NA>       <NA>
#> 1432    emotion    discuss  synthesis  consensus       plan       plan
#> 1433    monitor coregulate    monitor coregulate  consensus  consensus
#> 1434  consensus       plan  consensus    emotion       <NA>       <NA>
#> 1435      adapt    emotion    monitor      adapt    emotion       <NA>
#> 1436       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1437   cohesion  consensus    discuss  consensus       plan  consensus
#> 1438       plan       plan       plan    monitor    discuss  consensus
#> 1439       plan    discuss  synthesis coregulate   cohesion coregulate
#> 1440    emotion   cohesion  consensus       plan       plan       plan
#> 1441       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1442       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1443    discuss    discuss  consensus    emotion   cohesion coregulate
#> 1444       plan  consensus       plan       plan  consensus       plan
#> 1445       plan  consensus coregulate    discuss      adapt   cohesion
#> 1446    emotion  consensus    discuss    emotion  consensus coregulate
#> 1447 coregulate  consensus  consensus       plan    discuss  consensus
#> 1448       plan       plan   cohesion    monitor  consensus    monitor
#> 1449       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1450       plan       plan    discuss  consensus coregulate       plan
#> 1451 coregulate  consensus       plan    discuss    discuss    emotion
#> 1452       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1453  synthesis      adapt    discuss    discuss    emotion   cohesion
#> 1454   cohesion coregulate  consensus coregulate       plan   cohesion
#> 1455    discuss  consensus coregulate       plan       plan       plan
#> 1456       plan  consensus       plan       plan    emotion coregulate
#> 1457       plan       plan       plan       <NA>       <NA>       <NA>
#> 1458       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1459 coregulate       plan       plan  consensus       plan       <NA>
#> 1460       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1461       plan       plan       plan   cohesion       plan    monitor
#> 1462    discuss  consensus       plan       plan    monitor    emotion
#> 1463       plan       plan       plan    discuss   cohesion    discuss
#> 1464       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1465      adapt  consensus    monitor    discuss  consensus  consensus
#> 1466 coregulate    discuss   cohesion coregulate    discuss    monitor
#> 1467    discuss      adapt    discuss    discuss    emotion  consensus
#> 1468   cohesion  consensus  consensus       plan  consensus    discuss
#> 1469       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1470  consensus coregulate    discuss  consensus coregulate       plan
#> 1471  consensus    discuss  synthesis  consensus       plan       <NA>
#> 1472 coregulate  consensus       plan    emotion  consensus       plan
#> 1473  consensus    emotion   cohesion  consensus       plan  consensus
#> 1474  consensus       plan    emotion   cohesion  consensus    discuss
#> 1475  consensus       plan  consensus    monitor    discuss  consensus
#> 1476    emotion  consensus       plan  consensus    discuss       <NA>
#> 1477    emotion  consensus    discuss      adapt  consensus       plan
#> 1478       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1479    discuss  synthesis  consensus       <NA>       <NA>       <NA>
#> 1480  synthesis  consensus  consensus       <NA>       <NA>       <NA>
#> 1481       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1482       plan       plan       plan  consensus    discuss  synthesis
#> 1483       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1484       plan       plan    emotion  consensus coregulate       <NA>
#> 1485    emotion   cohesion  consensus   cohesion    emotion  consensus
#> 1486 coregulate    monitor    discuss  consensus       plan       plan
#> 1487    monitor    discuss       <NA>       <NA>       <NA>       <NA>
#> 1488    emotion   cohesion coregulate    monitor    emotion  consensus
#> 1489       plan       plan  consensus    monitor   cohesion       plan
#> 1490 coregulate       plan    monitor    discuss coregulate  consensus
#> 1491       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1492      adapt  consensus      adapt   cohesion  consensus       plan
#> 1493       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1494  consensus    discuss  consensus   cohesion coregulate    emotion
#> 1495    monitor  consensus    discuss      adapt    emotion   cohesion
#> 1496    discuss   cohesion coregulate       <NA>       <NA>       <NA>
#> 1497  consensus       plan    monitor    discuss coregulate    monitor
#> 1498       plan  consensus    monitor    discuss      adapt  consensus
#> 1499  consensus coregulate    emotion  consensus    monitor    discuss
#> 1500       plan       plan       plan       plan    discuss coregulate
#> 1501       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1502       plan    monitor    discuss    monitor       plan  consensus
#> 1503   cohesion    emotion    discuss coregulate       plan       plan
#> 1504  synthesis coregulate  synthesis    discuss      adapt    monitor
#> 1505    emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 1506       plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 1507       plan       plan       plan       plan    monitor       plan
#> 1508       plan  consensus    emotion    emotion    emotion    emotion
#> 1509 coregulate  consensus       plan       plan       plan    monitor
#> 1510  consensus coregulate    discuss  synthesis    emotion   cohesion
#> 1511       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1512  consensus       plan  consensus       plan       plan       plan
#> 1513       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1514       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1515  consensus  consensus coregulate    discuss    discuss    discuss
#> 1516  consensus       plan  consensus  consensus    discuss  synthesis
#> 1517       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1518    emotion    discuss coregulate  consensus    monitor  consensus
#> 1519       plan       plan  consensus coregulate    discuss      adapt
#> 1520       plan    discuss  consensus    monitor    monitor    discuss
#> 1521       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1522       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1523       plan    emotion   cohesion    discuss  synthesis  consensus
#> 1524       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1525       plan       plan  consensus coregulate    monitor   cohesion
#> 1526 coregulate   cohesion coregulate    monitor  consensus coregulate
#> 1527       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1528       plan    emotion coregulate    discuss  consensus  consensus
#> 1529       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1530       plan  consensus  consensus       plan  consensus       plan
#> 1531       plan       plan       plan  consensus  consensus  consensus
#> 1532  consensus  consensus       plan  consensus    emotion    monitor
#> 1533    monitor       plan  consensus       plan       <NA>       <NA>
#> 1534    discuss  consensus    discuss      adapt  consensus coregulate
#> 1535  consensus coregulate       plan    discuss    discuss  consensus
#> 1536       plan    discuss coregulate       plan  consensus coregulate
#> 1537  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1538       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1539  consensus coregulate    discuss  synthesis      adapt  consensus
#> 1540    discuss      adapt   cohesion       plan       plan       plan
#> 1541       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1542       plan  consensus       plan    discuss  synthesis    discuss
#> 1543       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1544  consensus       plan       plan    emotion    discuss      adapt
#> 1545    discuss  consensus coregulate  consensus       plan       <NA>
#> 1546       plan       plan coregulate    discuss    discuss coregulate
#> 1547       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1548       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1549    discuss  synthesis       <NA>       <NA>       <NA>       <NA>
#> 1550  consensus    emotion   cohesion coregulate    discuss  synthesis
#> 1551       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1552       plan coregulate    discuss  consensus       <NA>       <NA>
#> 1553       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1554      adapt  consensus coregulate    monitor    discuss   cohesion
#> 1555       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1556       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1557       plan       plan  consensus       plan  consensus    discuss
#> 1558       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1559  consensus       plan       plan    discuss    discuss  synthesis
#> 1560       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1561       plan    emotion  consensus       plan  consensus    emotion
#> 1562       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1563       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1564       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1565       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1566  synthesis  consensus       plan    monitor       <NA>       <NA>
#> 1567       plan  consensus coregulate    emotion  consensus       plan
#> 1568  consensus    discuss coregulate    discuss    discuss  consensus
#> 1569  consensus       plan       plan  consensus    discuss coregulate
#> 1570  consensus       plan       plan       plan       plan    discuss
#> 1571       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1572       plan       plan    emotion   cohesion  consensus       <NA>
#> 1573       plan  consensus       plan  consensus   cohesion    discuss
#> 1574       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1575  synthesis    discuss  synthesis      adapt   cohesion  consensus
#> 1576       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1577  consensus       plan  consensus coregulate coregulate coregulate
#> 1578       plan    monitor  consensus       plan       plan       plan
#> 1579       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1580       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1581  consensus       plan    monitor       plan       plan    emotion
#> 1582 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1583       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1584    discuss coregulate    discuss    emotion       plan       plan
#> 1585       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1586       plan    emotion  consensus coregulate    emotion    discuss
#> 1587   cohesion       plan       <NA>       <NA>       <NA>       <NA>
#> 1588       plan    monitor  consensus    discuss  synthesis  consensus
#> 1589    discuss  synthesis      adapt       <NA>       <NA>       <NA>
#> 1590    discuss coregulate       plan       plan    discuss coregulate
#> 1591       plan    emotion   cohesion  consensus coregulate  consensus
#> 1592       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1593    emotion    discuss  consensus       plan  consensus       plan
#> 1594       plan       plan       plan       plan    discuss coregulate
#> 1595 coregulate    discuss  consensus       plan  consensus       plan
#> 1596       plan       plan       plan    discuss    discuss    discuss
#> 1597    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1598       plan       plan       plan       <NA>       <NA>       <NA>
#> 1599      adapt  consensus       plan  consensus coregulate  consensus
#> 1600       plan       plan    monitor   cohesion    discuss  consensus
#> 1601  consensus    monitor       plan  consensus  consensus       plan
#> 1602   cohesion       plan       plan  consensus coregulate       plan
#> 1603  consensus coregulate    discuss       <NA>       <NA>       <NA>
#> 1604  consensus    emotion  consensus       plan       plan       plan
#> 1605       plan       plan       plan       plan       plan       plan
#> 1606    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1607  consensus    emotion  consensus       plan       plan       plan
#> 1608 coregulate    emotion       plan       plan       plan  consensus
#> 1609 coregulate  consensus  consensus    monitor   cohesion  consensus
#> 1610  consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 1611       plan  consensus coregulate    emotion  consensus coregulate
#> 1612       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1613    emotion  consensus       plan       plan       <NA>       <NA>
#> 1614       plan       plan       plan    discuss  consensus       <NA>
#> 1615       plan       plan    discuss  synthesis  consensus coregulate
#> 1616    emotion  consensus       plan       plan  consensus  consensus
#> 1617       plan       plan       plan       plan       <NA>       <NA>
#> 1618   cohesion    discuss    emotion       <NA>       <NA>       <NA>
#> 1619       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1620 coregulate       plan    emotion  consensus       plan  consensus
#> 1621       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1622       plan  consensus coregulate    emotion   cohesion    discuss
#> 1623       plan    emotion  consensus       <NA>       <NA>       <NA>
#> 1624       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1625   cohesion       plan       plan       <NA>       <NA>       <NA>
#> 1626       plan    emotion   cohesion    emotion   cohesion       <NA>
#> 1627    emotion  consensus coregulate    discuss  consensus       <NA>
#> 1628   cohesion  consensus coregulate       <NA>       <NA>       <NA>
#> 1629    discuss  consensus    monitor    discuss      adapt  consensus
#> 1630  consensus    emotion  consensus       <NA>       <NA>       <NA>
#> 1631  consensus       plan    emotion  consensus    discuss      adapt
#> 1632      adapt   cohesion coregulate  consensus  consensus    emotion
#> 1633    discuss    discuss   cohesion  consensus coregulate  consensus
#> 1634       plan    monitor    monitor    discuss    monitor    discuss
#> 1635 coregulate    emotion    emotion   cohesion    emotion  consensus
#> 1636  synthesis      adapt    discuss      adapt   cohesion  consensus
#> 1637 coregulate    discuss  synthesis    monitor  synthesis    discuss
#> 1638    monitor       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1639    monitor       plan       plan  consensus       plan  consensus
#> 1640 coregulate    discuss coregulate    monitor       plan    emotion
#> 1641    discuss  synthesis      adapt   cohesion  consensus  consensus
#> 1642 coregulate    discuss coregulate       <NA>       <NA>       <NA>
#> 1643    discuss      adapt   cohesion       <NA>       <NA>       <NA>
#> 1644   cohesion  consensus       plan  consensus    discuss  consensus
#> 1645       plan       plan       plan       plan       plan       plan
#> 1646      adapt   cohesion    monitor    discuss    discuss  consensus
#> 1647  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1648    emotion    discuss  consensus       plan       plan       plan
#> 1649       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1650       plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 1651       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1652  consensus    discuss  consensus       plan       plan    monitor
#> 1653   cohesion    emotion    monitor  consensus       plan       plan
#> 1654    discuss    discuss    emotion  consensus       plan  consensus
#> 1655       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1656       plan       plan    monitor    discuss       <NA>       <NA>
#> 1657    discuss coregulate    discuss  consensus       plan       plan
#> 1658      adapt  consensus coregulate    monitor       plan       <NA>
#> 1659    discuss    discuss  synthesis    emotion   cohesion    monitor
#> 1660   cohesion  consensus    discuss  consensus coregulate   cohesion
#> 1661  synthesis coregulate  consensus coregulate       <NA>       <NA>
#> 1662       plan    discuss  synthesis  consensus  consensus  synthesis
#> 1663   cohesion coregulate    emotion    discuss      adapt  consensus
#> 1664  consensus coregulate    monitor  consensus    monitor    emotion
#> 1665       plan    monitor    discuss  consensus coregulate       <NA>
#> 1666  consensus       plan  consensus  consensus    discuss    discuss
#> 1667    discuss  synthesis      adapt       <NA>       <NA>       <NA>
#> 1668       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1669       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1670       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1671       plan    emotion  consensus    discuss       plan    monitor
#> 1672  consensus    monitor   cohesion coregulate       plan  consensus
#> 1673  consensus       plan  consensus    discuss       <NA>       <NA>
#> 1674       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1675   cohesion  consensus    discuss    discuss  synthesis  consensus
#> 1676       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1677       plan  consensus  synthesis       <NA>       <NA>       <NA>
#> 1678  synthesis    emotion       <NA>       <NA>       <NA>       <NA>
#> 1679       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1680       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1681   cohesion  consensus    emotion   cohesion    emotion coregulate
#> 1682       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1683  consensus    monitor    emotion  consensus  consensus    discuss
#> 1684       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1685       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1686  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1687   cohesion  consensus coregulate    discuss  consensus    monitor
#> 1688       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1689       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1690       plan    emotion  consensus       <NA>       <NA>       <NA>
#> 1691       plan       plan    emotion    monitor  synthesis      adapt
#> 1692       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1693       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1694  consensus    monitor    emotion   cohesion coregulate  synthesis
#> 1695    discuss    discuss    discuss      adapt   cohesion       <NA>
#> 1696  consensus  consensus    discuss  synthesis      adapt  consensus
#> 1697       plan       plan       plan  consensus coregulate coregulate
#> 1698    discuss  consensus    discuss  synthesis       <NA>       <NA>
#> 1699   cohesion coregulate  consensus  consensus       plan       plan
#> 1700  consensus coregulate    monitor coregulate    monitor    discuss
#> 1701  consensus    emotion       plan    emotion    monitor       plan
#> 1702       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1703    emotion  consensus    discuss      adapt coregulate   cohesion
#> 1704       plan    emotion   cohesion       plan    emotion   cohesion
#> 1705    emotion  consensus       plan  consensus       plan       plan
#> 1706       plan  consensus coregulate  consensus  consensus    discuss
#> 1707    discuss    discuss    discuss    monitor  consensus coregulate
#> 1708       plan       plan coregulate    discuss  consensus    monitor
#> 1709       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1710       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1711  consensus    emotion    emotion   cohesion coregulate    discuss
#> 1712  consensus       plan       plan    monitor       <NA>       <NA>
#> 1713    monitor    discuss   cohesion    emotion       plan  consensus
#> 1714       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1715    discuss    discuss  consensus       plan       plan  consensus
#> 1716    discuss  synthesis  consensus coregulate    discuss  synthesis
#> 1717 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1718    emotion coregulate    discuss  consensus       plan    emotion
#> 1719  consensus       plan  consensus coregulate  synthesis  consensus
#> 1720 coregulate    emotion  consensus       plan       plan       plan
#> 1721       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1722  consensus       plan       plan       plan       plan    emotion
#> 1723       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1724    discuss  consensus       plan       plan       plan  consensus
#> 1725       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1726       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1727    monitor coregulate   cohesion  consensus       plan    discuss
#> 1728    emotion  consensus    discuss       <NA>       <NA>       <NA>
#> 1729    discuss  consensus coregulate       <NA>       <NA>       <NA>
#> 1730    monitor    discuss coregulate       <NA>       <NA>       <NA>
#> 1731       plan       plan    discuss       <NA>       <NA>       <NA>
#> 1732      adapt coregulate    discuss  consensus       <NA>       <NA>
#> 1733       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1734 coregulate    discuss    discuss  synthesis       plan    monitor
#> 1735    discuss  consensus    emotion  consensus    discuss    discuss
#> 1736    discuss    emotion   cohesion    emotion   cohesion  consensus
#> 1737  synthesis       plan       plan       plan  consensus       plan
#> 1738  consensus    monitor coregulate    discuss    discuss  consensus
#> 1739       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1740    discuss    discuss    discuss      adapt   cohesion       plan
#> 1741       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1742  synthesis  consensus coregulate  synthesis    discuss  consensus
#> 1743    discuss    discuss    emotion  consensus       plan       plan
#> 1744       plan       plan       plan    emotion       <NA>       <NA>
#> 1745       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1746       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1747    emotion  consensus coregulate    monitor  consensus       plan
#> 1748 coregulate    discuss    emotion  consensus       plan  consensus
#> 1749    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1750    emotion    discuss  consensus coregulate       plan  consensus
#> 1751  consensus  consensus       plan       plan    discuss       <NA>
#> 1752    emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 1753    discuss coregulate       plan  consensus   cohesion    emotion
#> 1754    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1755    discuss      adapt   cohesion  consensus  consensus       plan
#> 1756       plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 1757    discuss coregulate  consensus       plan  consensus       plan
#> 1758    emotion   cohesion  consensus       plan       plan  consensus
#> 1759       plan       plan  consensus    discuss  synthesis      adapt
#> 1760       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1761    emotion       plan       plan  consensus       plan       plan
#> 1762    monitor  consensus    discuss    discuss   cohesion    monitor
#> 1763  consensus coregulate    discuss      adapt   cohesion  consensus
#> 1764       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1765       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1766       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1767  consensus    emotion   cohesion coregulate       plan       plan
#> 1768    emotion  consensus    discuss  consensus       plan       <NA>
#> 1769  consensus       plan       plan       plan coregulate       <NA>
#> 1770 coregulate       plan  consensus    emotion    emotion  consensus
#> 1771  consensus    discuss  consensus       plan       plan  consensus
#> 1772       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1773    emotion  consensus       plan coregulate    discuss    discuss
#> 1774 coregulate       plan    discuss  synthesis  consensus       plan
#> 1775      adapt  consensus    emotion    monitor    discuss  consensus
#> 1776    discuss  consensus    monitor    discuss    emotion    emotion
#> 1777  synthesis       plan  consensus       plan  consensus coregulate
#> 1778       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1779  consensus   cohesion  consensus coregulate    discuss    discuss
#> 1780  synthesis      adapt  consensus    emotion coregulate       plan
#> 1781    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1782    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1783  consensus    monitor    monitor    discuss  consensus       plan
#> 1784  consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 1785       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1786    discuss    discuss    discuss  synthesis    discuss  consensus
#> 1787       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1788       plan       plan    discuss  consensus coregulate  consensus
#> 1789       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1790    discuss  synthesis   cohesion  consensus coregulate       <NA>
#> 1791       plan    emotion       plan    monitor coregulate       plan
#> 1792       plan       plan    monitor       <NA>       <NA>       <NA>
#> 1793       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1794  consensus  consensus       plan       plan  consensus    monitor
#> 1795      adapt   cohesion    emotion  consensus       plan       <NA>
#> 1796   cohesion  consensus       plan  consensus       plan       plan
#> 1797       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1798    emotion   cohesion  consensus coregulate  consensus coregulate
#> 1799  consensus    monitor       <NA>       <NA>       <NA>       <NA>
#> 1800  consensus    discuss  synthesis  consensus coregulate    discuss
#> 1801    emotion    discuss       <NA>       <NA>       <NA>       <NA>
#> 1802 coregulate    emotion       <NA>       <NA>       <NA>       <NA>
#> 1803    discuss  synthesis  consensus       plan       plan       <NA>
#> 1804      adapt    discuss      adapt  consensus coregulate       plan
#> 1805  consensus       plan       plan  consensus       plan       plan
#> 1806    discuss   cohesion    discuss  synthesis  consensus coregulate
#> 1807 coregulate    discuss       <NA>       <NA>       <NA>       <NA>
#> 1808       plan       plan       plan  consensus       <NA>       <NA>
#> 1809 coregulate    emotion   cohesion    emotion   cohesion  consensus
#> 1810    discuss       plan  consensus       <NA>       <NA>       <NA>
#> 1811   cohesion    monitor coregulate  consensus  consensus       plan
#> 1812   cohesion    emotion    emotion    monitor  consensus       plan
#> 1813       plan       plan    monitor    discuss      adapt  consensus
#> 1814       plan       plan  consensus coregulate      adapt   cohesion
#> 1815    emotion  consensus       plan  consensus       plan    discuss
#> 1816  consensus    monitor coregulate       plan       plan       plan
#> 1817       plan  consensus       plan  consensus       plan       plan
#> 1818    discuss  consensus       plan    monitor  consensus       plan
#> 1819       plan  consensus    discuss      adapt   cohesion  consensus
#> 1820 coregulate    emotion  consensus  synthesis       <NA>       <NA>
#> 1821       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1822       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1823       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1824    discuss    discuss    discuss    discuss  consensus       plan
#> 1825       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1826       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1827   cohesion       plan       plan  consensus       plan    discuss
#> 1828    emotion   cohesion  consensus       plan       plan    emotion
#> 1829  consensus  consensus       plan       plan    emotion    emotion
#> 1830    discuss  synthesis  consensus    monitor       plan       plan
#> 1831 coregulate    discuss       <NA>       <NA>       <NA>       <NA>
#> 1832  consensus       plan    monitor    discuss    discuss  synthesis
#> 1833    monitor    discuss  consensus       plan    discuss  synthesis
#> 1834       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1835    discuss   cohesion  consensus       <NA>       <NA>       <NA>
#> 1836  consensus    discuss coregulate    discuss    discuss  consensus
#> 1837       plan       plan    discuss    emotion    emotion    emotion
#> 1838    emotion  consensus    monitor    discuss coregulate    discuss
#> 1839    emotion       plan       plan       plan    monitor  consensus
#> 1840       plan coregulate    discuss  synthesis      adapt  consensus
#> 1841   cohesion  consensus  consensus coregulate       plan    monitor
#> 1842       plan       plan  consensus coregulate       plan   cohesion
#> 1843       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1844       plan       plan       plan       plan    emotion   cohesion
#> 1845    emotion   cohesion  consensus coregulate    emotion  consensus
#> 1846       plan  consensus       plan       plan    monitor    discuss
#> 1847  consensus coregulate    emotion       <NA>       <NA>       <NA>
#> 1848  consensus    emotion coregulate       plan  consensus coregulate
#> 1849       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1850       plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 1851       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1852 coregulate  consensus coregulate      adapt  consensus  consensus
#> 1853       plan  consensus  consensus coregulate    discuss  consensus
#> 1854  consensus    discuss  consensus       plan       plan    monitor
#> 1855    emotion       plan       plan  consensus       plan       plan
#> 1856  consensus       plan    emotion  consensus coregulate  consensus
#> 1857    discuss    discuss    emotion  consensus       plan  consensus
#> 1858  consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 1859  consensus       plan       plan       <NA>       <NA>       <NA>
#> 1860       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1861       plan  consensus coregulate  synthesis      adapt  consensus
#> 1862       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1863    emotion   cohesion  consensus    discuss      adapt  consensus
#> 1864    monitor    discuss    discuss    discuss  consensus       plan
#> 1865  consensus       plan       plan    discuss  consensus       plan
#> 1866  consensus    monitor  consensus       plan  consensus       plan
#> 1867  consensus coregulate       plan  consensus coregulate    emotion
#> 1868 coregulate  consensus    emotion    discuss      adapt   cohesion
#> 1869  consensus    emotion  consensus  consensus coregulate    monitor
#> 1870    emotion  consensus       plan       plan       plan       plan
#> 1871    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1872 coregulate   cohesion  consensus       plan    monitor    discuss
#> 1873       plan  consensus    discuss       <NA>       <NA>       <NA>
#> 1874       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1875 coregulate  consensus       plan    monitor    emotion   cohesion
#> 1876  consensus    monitor coregulate    discuss    discuss    discuss
#> 1877       plan    discuss  consensus       plan    monitor    emotion
#> 1878    monitor    discuss    monitor       plan  consensus       plan
#> 1879       plan       plan       plan  consensus  consensus  consensus
#> 1880       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1881    discuss  synthesis      adapt   cohesion    emotion    emotion
#> 1882    discuss  synthesis  consensus       plan    emotion  consensus
#> 1883       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1884    emotion   cohesion coregulate    monitor  consensus       plan
#> 1885  consensus    monitor    discuss    discuss    discuss    discuss
#> 1886       plan    monitor coregulate  consensus       plan       plan
#> 1887       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1888       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1889       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1890      adapt  consensus       plan       plan       plan       plan
#> 1891       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1892  consensus    discuss  synthesis  consensus       plan  consensus
#> 1893    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1894       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1895    discuss    emotion  consensus coregulate    discuss    discuss
#> 1896       plan       plan  consensus       <NA>       <NA>       <NA>
#> 1897       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1898       plan       plan  consensus    emotion   cohesion    discuss
#> 1899    discuss  consensus coregulate    monitor    discuss    emotion
#> 1900    monitor  consensus       plan  consensus       plan  consensus
#> 1901    emotion   cohesion coregulate    emotion  consensus       plan
#> 1902   cohesion  consensus       plan       plan    monitor    discuss
#> 1903       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1904 coregulate    emotion      adapt   cohesion    discuss      adapt
#> 1905       plan    discuss  consensus coregulate coregulate    monitor
#> 1906  consensus  consensus       plan    discuss    discuss coregulate
#> 1907 coregulate    discuss coregulate  consensus       plan       plan
#> 1908       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1909    monitor       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1910    emotion    discuss  synthesis    discuss  synthesis    discuss
#> 1911    discuss  synthesis      adapt  consensus coregulate    emotion
#> 1912       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1913  consensus       plan       plan       plan       plan       plan
#> 1914    monitor       plan    monitor   cohesion  consensus coregulate
#> 1915       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1916   cohesion  consensus       plan       <NA>       <NA>       <NA>
#> 1917       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1918       plan  consensus    discuss       <NA>       <NA>       <NA>
#> 1919       plan       plan  consensus  consensus       plan  consensus
#> 1920       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1921  consensus    discuss  consensus       plan    emotion   cohesion
#> 1922       plan    discuss      adapt  consensus       plan       plan
#> 1923       plan  consensus       plan    emotion    emotion       plan
#> 1924   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1925  consensus       plan    emotion coregulate  consensus    discuss
#> 1926       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1927       plan  consensus       plan       plan  consensus       <NA>
#> 1928  consensus       plan       plan       plan    discuss  consensus
#> 1929       plan    discuss       <NA>       <NA>       <NA>       <NA>
#> 1930  consensus       plan    discuss coregulate    emotion   cohesion
#> 1931       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1932       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1933    discuss  synthesis    emotion  consensus    monitor       <NA>
#> 1934   cohesion       plan       plan  consensus       <NA>       <NA>
#> 1935    discuss  consensus       plan  consensus    discuss coregulate
#> 1936  consensus       plan       plan  consensus    discuss    emotion
#> 1937       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1938    discuss  consensus coregulate  consensus       plan  consensus
#> 1939       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1940       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1941       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1942       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1943       plan    monitor       plan  consensus       plan       plan
#> 1944  consensus       plan       plan  consensus       plan  consensus
#> 1945       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1946       plan  consensus       plan    monitor       plan  consensus
#> 1947       plan       plan    discuss  synthesis      adapt   cohesion
#> 1948       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1949    discuss      adapt       <NA>       <NA>       <NA>       <NA>
#> 1950  consensus    discuss    discuss  consensus       plan       plan
#> 1951       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1952       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1953   cohesion    discuss coregulate    discuss    emotion       <NA>
#> 1954    discuss  synthesis  consensus       plan    discuss       <NA>
#> 1955       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1956 coregulate  consensus coregulate    discuss    discuss  consensus
#> 1957    emotion    emotion  consensus  consensus    discuss  consensus
#> 1958    discuss  consensus       plan  consensus       plan       plan
#> 1959       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1960       plan  consensus       plan  consensus    discuss  synthesis
#> 1961       plan       plan       plan       <NA>       <NA>       <NA>
#> 1962       plan  consensus    emotion  consensus    discuss  synthesis
#> 1963       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1964       plan  consensus    discuss  consensus       plan  consensus
#> 1965    discuss      adapt   cohesion       plan       plan  consensus
#> 1966       plan    emotion  consensus       plan  consensus coregulate
#> 1967   cohesion       plan    emotion       <NA>       <NA>       <NA>
#> 1968    monitor    discuss       <NA>       <NA>       <NA>       <NA>
#> 1969       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1970      adapt       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1971    monitor       plan  consensus       plan       <NA>       <NA>
#> 1972 coregulate    discuss    discuss    emotion  consensus       plan
#> 1973       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1974  consensus       plan  consensus    discuss coregulate       plan
#> 1975       plan    discuss    emotion  consensus       plan    discuss
#> 1976  consensus    discuss    discuss coregulate    discuss  consensus
#> 1977       plan  consensus      adapt  consensus  consensus       <NA>
#> 1978  consensus    emotion    emotion   cohesion coregulate       plan
#> 1979  consensus    monitor    discuss  consensus coregulate  consensus
#> 1980    emotion coregulate    monitor    emotion  consensus    emotion
#> 1981    discuss   cohesion  consensus coregulate  consensus coregulate
#> 1982  consensus       plan       plan  consensus  consensus    discuss
#> 1983       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1984       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1985       plan       plan    monitor    discuss  consensus coregulate
#> 1986  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1987    discuss  synthesis  consensus       plan       <NA>       <NA>
#> 1988  consensus  consensus       plan    emotion   cohesion  consensus
#> 1989   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1990       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1991    discuss    monitor  consensus  consensus       plan    monitor
#> 1992       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1993       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1994  synthesis      adapt    monitor    discuss    discuss  consensus
#> 1995  consensus    discuss coregulate    discuss    discuss   cohesion
#> 1996    emotion  consensus    emotion    discuss  consensus       plan
#> 1997 coregulate       plan       plan  consensus       plan   cohesion
#> 1998       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1999       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 2000  consensus       plan       plan    discuss       <NA>       <NA>
#>             T14        T15        T16        T17        T18        T19
#> 1          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 2          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 3     consensus       plan    emotion       plan       plan  consensus
#> 4          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 5          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 6         adapt   cohesion  consensus       plan  consensus       <NA>
#> 7          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 8    coregulate  consensus  consensus       plan       plan  consensus
#> 9          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 10         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 11    consensus       plan    emotion   cohesion       plan  consensus
#> 12   coregulate    monitor       plan  consensus    emotion  consensus
#> 13         plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 14         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 15   coregulate    emotion   cohesion       <NA>       <NA>       <NA>
#> 16         plan    emotion   cohesion  consensus coregulate    emotion
#> 17         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 18   coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 19         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 20         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 21      discuss   cohesion  consensus  consensus    emotion       plan
#> 22   coregulate    emotion  consensus       <NA>       <NA>       <NA>
#> 23         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 24         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 25     cohesion  consensus    discuss  synthesis  consensus       plan
#> 26         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 27         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 28      emotion   cohesion       plan    discuss       <NA>       <NA>
#> 29         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 30      discuss    discuss    discuss  consensus    emotion    emotion
#> 31   coregulate   cohesion coregulate    discuss  consensus    discuss
#> 32      discuss  consensus       plan coregulate       plan       plan
#> 33      discuss    emotion   cohesion  consensus    discuss  consensus
#> 34   coregulate  synthesis    emotion  consensus       plan    monitor
#> 35      emotion  consensus       plan  consensus       plan    discuss
#> 36    consensus    emotion       <NA>       <NA>       <NA>       <NA>
#> 37      emotion   cohesion  consensus       plan  consensus       plan
#> 38         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 39         plan    emotion  consensus    emotion   cohesion       plan
#> 40         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 41   coregulate       plan  consensus    discuss  consensus coregulate
#> 42         plan    emotion coregulate    emotion  consensus    discuss
#> 43         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 44    consensus coregulate    discuss       <NA>       <NA>       <NA>
#> 45    consensus       plan    emotion   cohesion       plan    emotion
#> 46         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 47      discuss  consensus       plan       <NA>       <NA>       <NA>
#> 48      discuss   cohesion  consensus    discuss  synthesis  consensus
#> 49         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 50   coregulate      adapt       <NA>       <NA>       <NA>       <NA>
#> 51    consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 52         plan       plan   cohesion coregulate    discuss    emotion
#> 53    consensus  consensus    discuss  consensus       plan   cohesion
#> 54         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 55      emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 56         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 57         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 58     cohesion  consensus    discuss   cohesion    emotion   cohesion
#> 59      monitor    discuss    discuss  consensus       plan    emotion
#> 60         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 61         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 62         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 63      discuss       plan    monitor  consensus coregulate       plan
#> 64         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 65         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 66    consensus       plan       plan       plan  consensus       plan
#> 67     cohesion    discuss  consensus       plan  consensus    discuss
#> 68         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 69         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 70         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 71    synthesis      adapt  consensus    discuss  synthesis      adapt
#> 72     cohesion  consensus    discuss    emotion       <NA>       <NA>
#> 73         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 74     cohesion  consensus       plan    monitor       <NA>       <NA>
#> 75         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 76         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 77      emotion   cohesion coregulate       plan       plan       plan
#> 78         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 79        adapt   cohesion       plan       plan  consensus       plan
#> 80         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 81         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 82         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 83      emotion  consensus    discuss  synthesis       plan       plan
#> 84      emotion  consensus       plan  consensus       plan    emotion
#> 85      monitor    discuss  consensus       plan  consensus    discuss
#> 86         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 87    consensus      adapt   cohesion  consensus       plan  consensus
#> 88    consensus coregulate       plan       <NA>       <NA>       <NA>
#> 89         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 90         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 91         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 92    consensus    discuss    emotion  consensus       plan    emotion
#> 93      discuss    discuss       plan    emotion       <NA>       <NA>
#> 94         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 95      discuss  synthesis  consensus    emotion       plan    emotion
#> 96   coregulate  consensus       plan  consensus coregulate    monitor
#> 97         plan  consensus       plan coregulate    emotion       <NA>
#> 98    consensus       plan  consensus       plan  consensus    discuss
#> 99    consensus  consensus  consensus    discuss  consensus    emotion
#> 100        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 101        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 102        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 103        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 104        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 105        plan       plan    monitor    discuss  consensus       plan
#> 106  coregulate    monitor       plan       plan  consensus       plan
#> 107        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 108   consensus    discuss  consensus       <NA>       <NA>       <NA>
#> 109   consensus    discuss coregulate    discuss       <NA>       <NA>
#> 110   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 111        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 112   consensus    discuss  consensus    emotion  consensus       plan
#> 113        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 114    cohesion       plan  synthesis       plan       <NA>       <NA>
#> 115  coregulate    discuss  consensus    monitor    emotion       plan
#> 116     emotion   cohesion  consensus    discuss       <NA>       <NA>
#> 117        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 118     discuss  consensus  consensus       plan    emotion    emotion
#> 119    cohesion  consensus coregulate       <NA>       <NA>       <NA>
#> 120        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 121        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 122        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 123        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 124        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 125        plan       plan       plan  consensus coregulate  consensus
#> 126        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 127   consensus    discuss  consensus coregulate       plan    emotion
#> 128   consensus    discuss    emotion       plan  consensus    monitor
#> 129        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 130        plan   cohesion       plan    monitor  consensus    discuss
#> 131        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 132        plan    emotion   cohesion       plan       plan  consensus
#> 133        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 134     emotion    discuss    emotion    discuss    emotion       plan
#> 135   consensus       plan  consensus    discuss  consensus   cohesion
#> 136       adapt  consensus    discuss  consensus       <NA>       <NA>
#> 137     emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 138        plan  consensus    discuss  consensus       plan    emotion
#> 139        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 140   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 141        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 142        plan  consensus coregulate       plan       plan  consensus
#> 143    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 144        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 145   consensus    discuss    discuss       <NA>       <NA>       <NA>
#> 146     emotion       plan    emotion  consensus       plan coregulate
#> 147  coregulate       plan    emotion   cohesion  consensus       <NA>
#> 148     discuss  consensus   cohesion  consensus coregulate       plan
#> 149        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 150        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 151        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 152     emotion   cohesion  consensus    emotion    emotion   cohesion
#> 153        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 154  coregulate  consensus    discuss  consensus       plan       plan
#> 155   consensus    emotion  consensus       plan  consensus coregulate
#> 156   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 157        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 158        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 159        plan  consensus    discuss       <NA>       <NA>       <NA>
#> 160        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 161   consensus       plan       plan       plan  consensus coregulate
#> 162   consensus  consensus   cohesion  consensus    discuss       <NA>
#> 163        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 164   synthesis  consensus       plan       plan       plan       plan
#> 165        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 166   synthesis  consensus coregulate    discuss   cohesion       <NA>
#> 167   consensus       plan       plan  consensus       plan       plan
#> 168   consensus  consensus    emotion       <NA>       <NA>       <NA>
#> 169   consensus    emotion    discuss       <NA>       <NA>       <NA>
#> 170        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 171     emotion  consensus  consensus   cohesion       <NA>       <NA>
#> 172        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 173        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 174        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 175        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 176     discuss  consensus    discuss    discuss    monitor  synthesis
#> 177       adapt  consensus    discuss coregulate  consensus       <NA>
#> 178     discuss    discuss  consensus    discuss  consensus coregulate
#> 179        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 180     discuss    discuss    emotion  consensus    discuss  synthesis
#> 181        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 182        plan  consensus       plan       <NA>       <NA>       <NA>
#> 183        plan  consensus       plan   cohesion  consensus    emotion
#> 184        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 185        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 186        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 187        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 188        plan  consensus    discuss    discuss  consensus    discuss
#> 189        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 190        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 191        plan  consensus  consensus       <NA>       <NA>       <NA>
#> 192        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 193        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 194        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 195        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 196        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 197    cohesion  consensus    discuss    discuss  consensus    discuss
#> 198        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 199        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 200        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 201   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 202        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 203        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 204     discuss coregulate    discuss    discuss  consensus       plan
#> 205        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 206        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 207     discuss   cohesion       <NA>       <NA>       <NA>       <NA>
#> 208        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 209        plan  consensus    discuss    discuss  synthesis  consensus
#> 210        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 211   consensus       plan  consensus       plan       plan       plan
#> 212     emotion  synthesis  consensus    discuss      adapt  consensus
#> 213        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 214  coregulate       plan    discuss  consensus    discuss   cohesion
#> 215        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 216        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 217        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 218        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 219        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 220        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 221        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 222        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 223        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 224        plan       plan  consensus    discuss  consensus  consensus
#> 225   consensus coregulate       plan    emotion   cohesion  consensus
#> 226  coregulate    discuss       <NA>       <NA>       <NA>       <NA>
#> 227    cohesion  consensus coregulate    discuss  consensus coregulate
#> 228        plan   cohesion  consensus    monitor       plan  consensus
#> 229        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 230  coregulate  consensus    monitor    discuss  consensus  consensus
#> 231   consensus  consensus       plan  consensus  consensus       plan
#> 232        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 233  coregulate    discuss  consensus    discuss  consensus coregulate
#> 234        plan       plan       plan  consensus       plan    emotion
#> 235    cohesion  synthesis  consensus  consensus    emotion    discuss
#> 236    cohesion  consensus       plan  consensus       plan       plan
#> 237        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 238        plan       plan    emotion       <NA>       <NA>       <NA>
#> 239        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 240     monitor    emotion       plan    monitor    discuss    emotion
#> 241   consensus       plan    emotion   cohesion    discuss   cohesion
#> 242    cohesion  consensus    emotion   cohesion  consensus       plan
#> 243     discuss  consensus       plan  consensus  synthesis       <NA>
#> 244        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 245        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 246   consensus    discuss  consensus    emotion    emotion       <NA>
#> 247     discuss  synthesis      adapt   cohesion  consensus    discuss
#> 248        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 249     emotion    discuss       <NA>       <NA>       <NA>       <NA>
#> 250        plan       plan  consensus       plan    monitor    emotion
#> 251        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 252   consensus  consensus coregulate       plan    discuss  consensus
#> 253        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 254        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 255        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 256        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 257        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 258   consensus       plan  consensus coregulate       <NA>       <NA>
#> 259        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 260        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 261        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 262        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 263   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 264     emotion   cohesion coregulate    emotion  consensus       plan
#> 265        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 266   consensus       plan       plan    emotion  consensus       plan
#> 267        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 268        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 269     discuss    discuss coregulate  synthesis    emotion   cohesion
#> 270     discuss  consensus       plan  consensus       <NA>       <NA>
#> 271        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 272   consensus       plan  consensus       plan       plan       plan
#> 273     discuss  synthesis       <NA>       <NA>       <NA>       <NA>
#> 274        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 275   consensus    emotion  consensus coregulate    discuss coregulate
#> 276        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 277        plan  consensus    discuss coregulate       plan    monitor
#> 278   consensus       plan       plan  consensus  consensus  synthesis
#> 279        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 280        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 281        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 282   consensus coregulate coregulate  consensus    discuss       <NA>
#> 283        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 284        plan       plan       plan       plan    emotion  consensus
#> 285   consensus    discuss  consensus       <NA>       <NA>       <NA>
#> 286        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 287        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 288     emotion   cohesion  consensus       plan  consensus coregulate
#> 289  coregulate    discuss  consensus coregulate    discuss       <NA>
#> 290        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 291     emotion   cohesion       plan  consensus       plan    monitor
#> 292     emotion  consensus coregulate       plan   cohesion  consensus
#> 293     emotion  consensus    discuss    emotion  consensus    discuss
#> 294        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 295     discuss  consensus    emotion  consensus    discuss    emotion
#> 296   consensus    discuss coregulate    discuss       <NA>       <NA>
#> 297        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 298        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 299        plan    monitor    emotion  consensus    discuss      adapt
#> 300     monitor    discuss       plan  consensus coregulate  consensus
#> 301        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 302        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 303        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 304        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 305        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 306     discuss coregulate    discuss  consensus       <NA>       <NA>
#> 307        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 308        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 309        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 310        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 311        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 312     discuss coregulate    monitor  consensus    discuss  consensus
#> 313        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 314        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 315     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 316   consensus       plan  consensus coregulate coregulate       plan
#> 317        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 318   synthesis  consensus       plan       plan coregulate    discuss
#> 319        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 320        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 321        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 322        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 323        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 324        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 325    cohesion    emotion    monitor       <NA>       <NA>       <NA>
#> 326        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 327        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 328   consensus coregulate       plan       plan    monitor    emotion
#> 329     emotion   cohesion  consensus       plan    discuss  consensus
#> 330        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 331        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 332   consensus    emotion   cohesion       plan  consensus    discuss
#> 333        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 334  coregulate    discuss  consensus       plan    emotion   cohesion
#> 335        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 336        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 337   consensus coregulate    discuss    emotion   cohesion    emotion
#> 338     emotion   cohesion  consensus    discuss coregulate  consensus
#> 339   consensus    discuss  consensus coregulate    discuss    emotion
#> 340   consensus       plan       plan    monitor  consensus    discuss
#> 341        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 342        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 343   consensus coregulate  consensus    discuss    discuss    discuss
#> 344   consensus       plan  consensus       <NA>       <NA>       <NA>
#> 345     discuss    emotion  consensus       <NA>       <NA>       <NA>
#> 346    cohesion       plan   cohesion  consensus       <NA>       <NA>
#> 347  coregulate   cohesion       plan    emotion       plan       plan
#> 348   consensus    discuss    discuss  synthesis  consensus    emotion
#> 349        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 350        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 351   consensus       plan  consensus  consensus    discuss  consensus
#> 352        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 353     monitor coregulate  consensus       <NA>       <NA>       <NA>
#> 354     emotion  synthesis       plan       plan  consensus coregulate
#> 355       adapt  consensus       <NA>       <NA>       <NA>       <NA>
#> 356        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 357        plan       plan       plan  consensus       plan    emotion
#> 358        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 359        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 360        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 361        plan       plan       plan  consensus       plan  consensus
#> 362        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 363   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 364        plan coregulate    discuss       <NA>       <NA>       <NA>
#> 365  coregulate    emotion    discuss       <NA>       <NA>       <NA>
#> 366  coregulate  consensus       <NA>       <NA>       <NA>       <NA>
#> 367        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 368        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 369        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 370     emotion   cohesion  consensus coregulate    discuss coregulate
#> 371        plan       plan  consensus  synthesis       <NA>       <NA>
#> 372        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 373        plan  consensus    discuss  consensus       plan       plan
#> 374     discuss  synthesis    emotion    monitor   cohesion  consensus
#> 375   consensus    discuss  consensus    discuss  consensus    emotion
#> 376        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 377        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 378        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 379        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 380        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 381   consensus    monitor       plan       plan       plan  consensus
#> 382     discuss  consensus    discuss  consensus      adapt   cohesion
#> 383        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 384        plan       plan  consensus    discuss    emotion  consensus
#> 385        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 386     discuss  synthesis  consensus    discuss  consensus       <NA>
#> 387        plan    discuss       <NA>       <NA>       <NA>       <NA>
#> 388        plan    monitor    discuss  consensus  consensus       plan
#> 389        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 390   consensus    discuss  consensus  consensus       plan  consensus
#> 391  coregulate    emotion    monitor  consensus    discuss       <NA>
#> 392        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 393     discuss  consensus       plan       plan coregulate    monitor
#> 394        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 395        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 396   consensus    discuss    emotion       <NA>       <NA>       <NA>
#> 397        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 398     emotion       plan    discuss  synthesis  consensus       plan
#> 399        plan    monitor       plan  consensus    discuss  consensus
#> 400        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 401        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 402     discuss  consensus       plan    monitor   cohesion  consensus
#> 403        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 404   consensus       plan    discuss  consensus  consensus       plan
#> 405   synthesis    emotion   cohesion       plan    emotion       <NA>
#> 406        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 407        plan       plan       plan       plan       plan       plan
#> 408        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 409        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 410  coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 411        plan       plan    emotion       plan  consensus  consensus
#> 412        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 413        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 414        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 415        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 416  coregulate    monitor    emotion  consensus  consensus      adapt
#> 417        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 418        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 419     discuss  consensus  synthesis  consensus    monitor    emotion
#> 420   consensus  consensus    discuss   cohesion       plan       plan
#> 421     discuss coregulate       plan       plan       plan       <NA>
#> 422   consensus       plan    emotion       plan       plan    monitor
#> 423     emotion    discuss   cohesion  consensus    discuss       <NA>
#> 424        plan  consensus   cohesion    emotion  consensus    discuss
#> 425        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 426  coregulate    monitor       plan       <NA>       <NA>       <NA>
#> 427     discuss    emotion   cohesion       plan       plan       plan
#> 428     discuss  consensus       plan       plan   cohesion  consensus
#> 429     emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 430        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 431     emotion   cohesion    emotion       plan       plan    emotion
#> 432     monitor       plan    emotion   cohesion    emotion  consensus
#> 433     discuss  synthesis  consensus    discuss  synthesis  consensus
#> 434     monitor       plan       plan  consensus    discuss    emotion
#> 435     discuss    emotion   cohesion  consensus       <NA>       <NA>
#> 436   consensus    monitor    discuss  consensus coregulate    emotion
#> 437        plan  consensus    monitor       plan       plan    emotion
#> 438        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 439     emotion       plan       <NA>       <NA>       <NA>       <NA>
#> 440        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 441        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 442        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 443     emotion   cohesion    emotion       plan       plan    emotion
#> 444   consensus    discuss  consensus    discuss    discuss    discuss
#> 445        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 446        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 447        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 448        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 449        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 450        plan  consensus    discuss    emotion       <NA>       <NA>
#> 451     discuss  consensus coregulate    monitor       <NA>       <NA>
#> 452        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 453        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 454     discuss  consensus coregulate       <NA>       <NA>       <NA>
#> 455        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 456        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 457        plan       plan coregulate  consensus  consensus   cohesion
#> 458        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 459        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 460        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 461     emotion   cohesion    discuss  synthesis  consensus       plan
#> 462   synthesis  consensus  synthesis   cohesion  consensus       <NA>
#> 463   consensus    discuss      adapt  consensus    discuss  consensus
#> 464        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 465        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 466        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 467   consensus       plan       plan       <NA>       <NA>       <NA>
#> 468        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 469   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 470        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 471    cohesion  consensus coregulate    emotion  consensus    discuss
#> 472        plan  consensus coregulate       plan    emotion    monitor
#> 473        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 474   consensus    emotion    discuss  consensus       plan coregulate
#> 475   consensus       plan       plan    monitor  synthesis  consensus
#> 476        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 477     discuss coregulate       plan    discuss coregulate       plan
#> 478     discuss  consensus       plan  consensus    monitor    emotion
#> 479   consensus       plan  consensus    discuss    emotion    monitor
#> 480    cohesion   cohesion  consensus  consensus       plan  consensus
#> 481        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 482     emotion  consensus       plan       plan   cohesion  consensus
#> 483        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 484        plan    emotion  consensus    discuss  synthesis       plan
#> 485        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 486        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 487     discuss  synthesis  consensus   cohesion  consensus    discuss
#> 488        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 489        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 490        plan       plan       plan    emotion  consensus    emotion
#> 491        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 492   consensus    discuss  consensus coregulate    discuss  consensus
#> 493        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 494        plan       plan       plan    discuss  synthesis  consensus
#> 495        plan       plan    emotion    discuss    discuss       <NA>
#> 496        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 497   consensus       plan       plan  consensus    discuss  consensus
#> 498        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 499        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 500     discuss  consensus    discuss  consensus    discuss       <NA>
#> 501     discuss    emotion       plan       plan       <NA>       <NA>
#> 502     emotion coregulate       plan       plan       <NA>       <NA>
#> 503   consensus coregulate    monitor    discuss       <NA>       <NA>
#> 504        plan    emotion  consensus       plan       plan       plan
#> 505        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 506        plan  consensus    discuss    discuss  consensus       <NA>
#> 507   consensus       plan    emotion  consensus  consensus       plan
#> 508        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 509        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 510        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 511        plan  consensus coregulate       plan   cohesion  consensus
#> 512     emotion  consensus       plan    monitor    discuss  consensus
#> 513        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 514        plan    monitor       plan  consensus       plan    emotion
#> 515        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 516        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 517   consensus    discuss  consensus    discuss  consensus    discuss
#> 518   consensus       plan    emotion    discuss    discuss  consensus
#> 519        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 520    cohesion  consensus    discuss    discuss       <NA>       <NA>
#> 521     emotion   cohesion    discuss      adapt  consensus    emotion
#> 522        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 523        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 524        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 525        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 526        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 527        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 528        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 529        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 530        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 531        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 532        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 533        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 534   consensus coregulate    emotion   cohesion  consensus       plan
#> 535        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 536  coregulate    discuss coregulate    discuss  consensus    discuss
#> 537        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 538        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 539        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 540  coregulate    emotion    emotion   cohesion  consensus       plan
#> 541   consensus       plan       plan    emotion  consensus       <NA>
#> 542        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 543   consensus    discuss coregulate    discuss  consensus       <NA>
#> 544     emotion    discuss    discuss    monitor       <NA>       <NA>
#> 545     monitor    discuss  consensus    discuss       <NA>       <NA>
#> 546        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 547        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 548     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 549   consensus       plan    emotion   cohesion   cohesion    emotion
#> 550     discuss  consensus       plan  consensus    emotion   cohesion
#> 551        plan    emotion       plan  consensus    discuss    emotion
#> 552        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 553        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 554     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 555        plan    emotion   cohesion  consensus       plan       plan
#> 556     emotion    monitor    discuss coregulate       plan       plan
#> 557     discuss  consensus    discuss coregulate       plan    discuss
#> 558        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 559        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 560        plan    monitor    discuss  synthesis  consensus       plan
#> 561        plan       plan  consensus       plan    emotion    discuss
#> 562  coregulate  consensus       plan       <NA>       <NA>       <NA>
#> 563        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 564        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 565        plan  consensus       plan       plan       plan    emotion
#> 566        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 567     monitor       <NA>       <NA>       <NA>       <NA>       <NA>
#> 568     monitor  synthesis   cohesion  consensus       plan       plan
#> 569   consensus       plan  consensus    emotion  consensus coregulate
#> 570  coregulate    discuss    discuss       <NA>       <NA>       <NA>
#> 571     emotion       plan    emotion    discuss    discuss    discuss
#> 572        plan  consensus       plan    monitor    discuss  synthesis
#> 573    cohesion  consensus       plan  consensus       <NA>       <NA>
#> 574    cohesion       plan    emotion       plan       plan    monitor
#> 575    cohesion       plan    emotion    discuss  synthesis  consensus
#> 576        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 577        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 578        plan    emotion       plan  consensus    discuss    discuss
#> 579    cohesion  consensus    emotion       <NA>       <NA>       <NA>
#> 580        plan  consensus   cohesion coregulate    discuss    discuss
#> 581        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 582        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 583   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 584        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 585        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 586  coregulate    monitor    discuss    discuss  consensus       plan
#> 587        plan   cohesion    emotion    emotion   cohesion  consensus
#> 588        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 589        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 590     discuss  consensus    discuss  consensus    emotion       <NA>
#> 591        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 592  coregulate       plan       plan    emotion  consensus   cohesion
#> 593     discuss coregulate       plan    emotion   cohesion  consensus
#> 594   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 595        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 596        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 597        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 598        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 599     discuss  consensus      adapt  consensus coregulate    discuss
#> 600        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 601        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 602        plan       plan  consensus  consensus       plan  consensus
#> 603        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 604        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 605        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 606        plan  consensus    emotion  consensus    discuss coregulate
#> 607  coregulate       plan    discuss  consensus       <NA>       <NA>
#> 608        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 609     emotion   cohesion       plan  consensus coregulate   cohesion
#> 610        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 611   consensus    emotion  consensus    discuss  consensus    discuss
#> 612        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 613   consensus       plan  consensus    discuss    discuss  consensus
#> 614        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 615    cohesion    emotion  consensus       plan    emotion       <NA>
#> 616        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 617        plan       plan  consensus    discuss    discuss  consensus
#> 618        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 619     monitor       plan    monitor    discuss    discuss   cohesion
#> 620        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 621        plan    monitor coregulate       plan       plan    monitor
#> 622     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 623     emotion   cohesion    emotion    discuss    emotion   cohesion
#> 624     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 625        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 626        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 627        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 628        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 629        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 630        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 631        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 632     emotion    discuss   cohesion  consensus coregulate    monitor
#> 633     discuss coregulate    discuss  consensus       plan  consensus
#> 634     discuss  consensus    discuss  consensus coregulate    discuss
#> 635     discuss  consensus    discuss   cohesion  consensus       plan
#> 636   synthesis  consensus coregulate  consensus    discuss  consensus
#> 637        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 638        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 639        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 640   consensus    discuss    discuss       <NA>       <NA>       <NA>
#> 641   consensus    emotion   cohesion   cohesion  consensus    discuss
#> 642    cohesion  consensus       plan       plan    emotion    emotion
#> 643        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 644   consensus coregulate    emotion       plan       plan    emotion
#> 645     discuss    discuss  consensus       <NA>       <NA>       <NA>
#> 646        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 647        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 648   consensus       plan       plan       plan    monitor      adapt
#> 649        plan       plan    discuss  consensus  consensus    discuss
#> 650        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 651     emotion    monitor   cohesion       <NA>       <NA>       <NA>
#> 652        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 653     discuss  consensus    emotion  consensus       plan  consensus
#> 654        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 655     emotion    discuss       <NA>       <NA>       <NA>       <NA>
#> 656        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 657     emotion  consensus    discuss    discuss  consensus coregulate
#> 658        plan    discuss  consensus       plan  consensus coregulate
#> 659    cohesion       plan       plan    monitor  consensus       plan
#> 660   consensus    emotion    emotion   cohesion       plan    emotion
#> 661        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 662        plan  consensus       plan    monitor       plan       plan
#> 663        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 664       adapt  consensus       plan       plan       plan  consensus
#> 665        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 666        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 667        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 668   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 669   synthesis  consensus       <NA>       <NA>       <NA>       <NA>
#> 670        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 671   consensus       plan    emotion  consensus       plan    emotion
#> 672        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 673   consensus       plan       plan  consensus coregulate    discuss
#> 674   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 675   consensus coregulate    monitor    emotion    emotion    emotion
#> 676        plan       plan       <NA>       <NA>       <NA>       <NA>
#> 677    cohesion       plan  consensus  consensus       plan  consensus
#> 678     emotion   cohesion  consensus coregulate coregulate    discuss
#> 679     emotion      adapt  consensus       <NA>       <NA>       <NA>
#> 680        plan  consensus       plan    emotion   cohesion    emotion
#> 681        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 682        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 683        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 684   consensus coregulate       plan       <NA>       <NA>       <NA>
#> 685        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 686     discuss  consensus    discuss  consensus coregulate       plan
#> 687   consensus       plan    discuss  consensus       plan       plan
#> 688    cohesion  consensus       plan   cohesion  consensus       <NA>
#> 689    cohesion  consensus       plan  consensus    discuss    discuss
#> 690     discuss    emotion   cohesion coregulate  consensus    discuss
#> 691        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 692     discuss  consensus    monitor    discuss       <NA>       <NA>
#> 693  coregulate       plan  consensus    monitor  consensus coregulate
#> 694     emotion       plan       <NA>       <NA>       <NA>       <NA>
#> 695        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 696  coregulate    monitor    discuss  consensus    discuss       <NA>
#> 697        plan       plan       plan    discuss coregulate    monitor
#> 698        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 699        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 700     emotion       plan  consensus coregulate       <NA>       <NA>
#> 701        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 702     emotion    discuss    discuss  consensus       <NA>       <NA>
#> 703     discuss  consensus       plan       plan       plan       plan
#> 704        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 705        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 706        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 707        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 708        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 709        plan    monitor  consensus    discuss  consensus  consensus
#> 710     emotion  consensus       plan       plan  consensus    discuss
#> 711        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 712     emotion   cohesion    emotion coregulate    discuss  synthesis
#> 713        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 714        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 715     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 716     discuss  consensus coregulate       <NA>       <NA>       <NA>
#> 717        plan    monitor    emotion   cohesion  consensus       plan
#> 718        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 719        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 720   consensus       plan       plan  consensus coregulate  consensus
#> 721        plan       plan       plan    emotion  consensus coregulate
#> 722        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 723        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 724        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 725        plan    emotion   cohesion    emotion  consensus    monitor
#> 726        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 727        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 728   consensus    discuss    emotion   cohesion coregulate  consensus
#> 729        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 730        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 731     discuss  consensus       plan  consensus    discuss      adapt
#> 732        plan    emotion   cohesion       plan    emotion   cohesion
#> 733        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 734        plan coregulate  consensus       plan  consensus       <NA>
#> 735   consensus    emotion   cohesion  consensus    emotion  consensus
#> 736     discuss    discuss    discuss    discuss    emotion coregulate
#> 737  coregulate  consensus       <NA>       <NA>       <NA>       <NA>
#> 738   consensus    discuss    discuss  consensus       <NA>       <NA>
#> 739        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 740     discuss  synthesis  consensus  consensus    emotion  consensus
#> 741     discuss  consensus  consensus       <NA>       <NA>       <NA>
#> 742    cohesion coregulate    discuss       <NA>       <NA>       <NA>
#> 743        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 744        plan    discuss  consensus    discuss    discuss    emotion
#> 745        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 746   consensus    emotion   cohesion  consensus       plan coregulate
#> 747        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 748        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 749     monitor   cohesion  consensus    discuss  consensus    emotion
#> 750     monitor  consensus coregulate   cohesion  consensus  consensus
#> 751     monitor    emotion   cohesion coregulate    emotion    emotion
#> 752        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 753        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 754        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 755        plan  consensus       plan       plan       plan  consensus
#> 756   consensus       plan       plan  consensus    emotion       <NA>
#> 757        plan    emotion    discuss    emotion  consensus    monitor
#> 758   synthesis       <NA>       <NA>       <NA>       <NA>       <NA>
#> 759     monitor    discuss      adapt  consensus       <NA>       <NA>
#> 760        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 761        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 762        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 763        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 764        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 765        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 766        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 767    cohesion  consensus    discuss  consensus  consensus       plan
#> 768   consensus    emotion    discuss    discuss  consensus       plan
#> 769        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 770        plan    emotion    discuss coregulate    discuss    emotion
#> 771        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 772        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 773        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 774   consensus    discuss  consensus coregulate    emotion       plan
#> 775     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 776        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 777        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 778        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 779        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 780        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 781        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 782   consensus       plan  synthesis   cohesion  consensus       plan
#> 783        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 784     monitor      adapt   cohesion    discuss  consensus       plan
#> 785     emotion   cohesion       plan       plan       <NA>       <NA>
#> 786        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 787        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 788        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 789   consensus       plan  consensus coregulate       plan    emotion
#> 790        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 791        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 792    cohesion       plan       plan       plan    emotion       plan
#> 793        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 794     emotion    emotion  consensus       plan    monitor       <NA>
#> 795        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 796        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 797        plan  consensus    discuss  consensus       plan       plan
#> 798        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 799        plan  consensus       plan       plan    discuss    discuss
#> 800        plan  consensus       plan  consensus       <NA>       <NA>
#> 801        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 802    cohesion  consensus    discuss       <NA>       <NA>       <NA>
#> 803        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 804    cohesion  consensus       plan coregulate       plan    emotion
#> 805   consensus    emotion       plan  consensus coregulate    monitor
#> 806   consensus    discuss  consensus    discuss  consensus coregulate
#> 807        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 808        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 809   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 810        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 811     discuss  consensus  consensus    emotion  consensus       <NA>
#> 812        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 813        plan    emotion coregulate    emotion       <NA>       <NA>
#> 814        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 815        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 816        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 817     emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 818        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 819   consensus       plan coregulate    emotion    discuss   cohesion
#> 820        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 821        plan    discuss       <NA>       <NA>       <NA>       <NA>
#> 822        plan  consensus    discuss    discuss  consensus       plan
#> 823   consensus    emotion    discuss  synthesis       <NA>       <NA>
#> 824     monitor   cohesion       plan    emotion    discuss  consensus
#> 825   consensus coregulate    discuss  consensus    monitor    discuss
#> 826        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 827   consensus       plan    monitor       <NA>       <NA>       <NA>
#> 828     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 829        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 830     emotion  consensus    discuss    emotion       plan    emotion
#> 831        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 832   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 833   consensus    discuss  consensus coregulate    discuss  consensus
#> 834  coregulate    monitor    emotion   cohesion    emotion    discuss
#> 835        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 836        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 837     emotion  consensus  consensus coregulate   cohesion  consensus
#> 838        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 839        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 840        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 841        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 842        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 843     discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 844        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 845        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 846        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 847  coregulate   cohesion       <NA>       <NA>       <NA>       <NA>
#> 848        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 849        plan  consensus       plan    monitor       <NA>       <NA>
#> 850  coregulate    discuss    discuss  consensus    discuss   cohesion
#> 851        plan       plan       plan       plan    emotion   cohesion
#> 852        plan       plan       <NA>       <NA>       <NA>       <NA>
#> 853     discuss    discuss  synthesis      adapt   cohesion  consensus
#> 854        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 855        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 856        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 857   consensus       plan       plan   cohesion       <NA>       <NA>
#> 858        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 859        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 860        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 861        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 862     emotion  consensus  consensus coregulate    emotion  consensus
#> 863   consensus coregulate       plan       plan    emotion   cohesion
#> 864        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 865        plan    monitor       plan    monitor       plan       plan
#> 866        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 867     emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 868     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 869   consensus coregulate    emotion  consensus       plan       plan
#> 870        plan    emotion    discuss  synthesis      adapt    discuss
#> 871        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 872        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 873     emotion   cohesion       plan    emotion  consensus       plan
#> 874        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 875        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 876  coregulate    emotion coregulate  consensus    emotion    discuss
#> 877     emotion   cohesion       plan       plan       plan  consensus
#> 878   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 879        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 880        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 881        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 882        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 883   consensus       plan  consensus   cohesion    monitor    discuss
#> 884        plan    emotion   cohesion coregulate       plan       <NA>
#> 885     emotion  consensus coregulate      adapt    emotion  consensus
#> 886     discuss    emotion       <NA>       <NA>       <NA>       <NA>
#> 887        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 888        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 889        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 890        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 891        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 892   consensus       plan    emotion   cohesion       <NA>       <NA>
#> 893        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 894        plan    emotion    discuss  consensus    monitor       <NA>
#> 895        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 896        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 897        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 898  coregulate    monitor    discuss coregulate    discuss  synthesis
#> 899        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 900        plan    emotion   cohesion  consensus    emotion    discuss
#> 901        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 902        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 903        plan       plan    emotion    emotion  consensus  synthesis
#> 904        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 905        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 906        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 907        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 908        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 909     discuss  consensus coregulate    discuss    emotion    discuss
#> 910   consensus  consensus    emotion    emotion       plan    emotion
#> 911        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 912        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 913        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 914        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 915        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 916        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 917     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 918        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 919   consensus    emotion  consensus       plan    emotion   cohesion
#> 920    cohesion  consensus coregulate       plan    emotion   cohesion
#> 921        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 922        plan  consensus coregulate       <NA>       <NA>       <NA>
#> 923        plan       plan    monitor    discuss       <NA>       <NA>
#> 924        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 925     monitor    emotion  consensus       plan       plan       plan
#> 926        plan    emotion   cohesion coregulate       <NA>       <NA>
#> 927        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 928     discuss   cohesion  consensus  consensus       <NA>       <NA>
#> 929        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 930        plan    discuss coregulate      adapt       <NA>       <NA>
#> 931        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 932        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 933   consensus       plan    emotion  synthesis      adapt  consensus
#> 934        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 935     discuss    discuss       <NA>       <NA>       <NA>       <NA>
#> 936        plan  consensus coregulate       <NA>       <NA>       <NA>
#> 937        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 938        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 939   consensus       plan       plan coregulate    discuss  consensus
#> 940    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 941        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 942        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 943        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 944        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 945        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 946        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 947        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 948   consensus    discuss    discuss    emotion   cohesion  consensus
#> 949   consensus       plan       plan  consensus    discuss coregulate
#> 950        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 951     monitor    discuss  consensus       plan  consensus coregulate
#> 952     discuss    discuss      adapt  consensus    discuss    emotion
#> 953   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 954  coregulate    discuss   cohesion    emotion   cohesion  consensus
#> 955    cohesion    emotion       plan       <NA>       <NA>       <NA>
#> 956   consensus       plan       plan  consensus       plan    emotion
#> 957        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 958        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 959        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 960     discuss  consensus       plan       plan  consensus coregulate
#> 961     monitor    discuss    emotion  consensus coregulate    discuss
#> 962   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 963     emotion    discuss       <NA>       <NA>       <NA>       <NA>
#> 964        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 965        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 966        plan    emotion   cohesion       plan       plan  consensus
#> 967   consensus       plan       plan    monitor    discuss    discuss
#> 968        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 969        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 970     discuss    discuss  consensus    emotion    emotion    discuss
#> 971     monitor       plan  consensus       plan       plan    emotion
#> 972        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 973        plan    monitor       plan    monitor       <NA>       <NA>
#> 974        plan    monitor    emotion    monitor       plan coregulate
#> 975   consensus coregulate       plan    emotion       plan  consensus
#> 976        plan       plan       <NA>       <NA>       <NA>       <NA>
#> 977        plan  consensus       plan    discuss  consensus       plan
#> 978        plan  consensus  consensus coregulate    emotion  consensus
#> 979        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 980        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 981        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 982        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 983  coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 984        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 985        plan       plan       plan    emotion       <NA>       <NA>
#> 986        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 987        plan  consensus    discuss  consensus   cohesion  consensus
#> 988     emotion    discuss       <NA>       <NA>       <NA>       <NA>
#> 989    cohesion  consensus    discuss    discuss  consensus    emotion
#> 990     emotion    discuss  consensus coregulate    emotion    discuss
#> 991        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 992        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 993   consensus    monitor    emotion       <NA>       <NA>       <NA>
#> 994        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 995        plan    emotion  consensus    discuss    emotion  consensus
#> 996     emotion   cohesion    emotion    emotion   cohesion    discuss
#> 997   consensus    discuss  synthesis       plan       plan  consensus
#> 998     discuss  consensus    discuss    discuss    discuss      adapt
#> 999        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1000       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1001       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1002       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1003  synthesis    discuss coregulate       <NA>       <NA>       <NA>
#> 1004    emotion    emotion    emotion   cohesion coregulate  consensus
#> 1005       plan  consensus    emotion       <NA>       <NA>       <NA>
#> 1006       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1007       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1008    emotion  consensus       plan    discuss  consensus coregulate
#> 1009       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1010  synthesis  consensus    discuss  consensus       plan       plan
#> 1011       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1012       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1013       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1014    monitor coregulate  consensus       plan  consensus    discuss
#> 1015   cohesion    monitor   cohesion coregulate       plan  consensus
#> 1016    monitor  consensus coregulate    emotion  consensus    discuss
#> 1017       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1018       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1019       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1020       plan       plan  consensus       plan       plan       plan
#> 1021       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1022       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1023  consensus  consensus  consensus       plan       plan       plan
#> 1024    emotion  consensus       plan       plan       plan  consensus
#> 1025       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1026       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1027       plan    discuss      adapt  consensus       plan       plan
#> 1028       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1029       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1030  consensus    discuss  consensus       plan       plan    discuss
#> 1031    discuss    discuss    emotion   cohesion  consensus  consensus
#> 1032    emotion       plan       plan  consensus    discuss  synthesis
#> 1033      adapt   cohesion    monitor    discuss coregulate       plan
#> 1034  consensus coregulate    discuss    discuss    discuss       <NA>
#> 1035  consensus coregulate coregulate    emotion    emotion coregulate
#> 1036       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1037       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1038  consensus       plan  consensus coregulate    emotion       <NA>
#> 1039  consensus  consensus       <NA>       <NA>       <NA>       <NA>
#> 1040       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1041       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1042   cohesion    emotion       plan       plan       <NA>       <NA>
#> 1043       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1044       plan       plan       plan    discuss    monitor       <NA>
#> 1045       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1046    emotion    monitor    discuss  synthesis      adapt   cohesion
#> 1047       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1048       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1049       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1050       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1051       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1052    discuss    discuss    discuss       <NA>       <NA>       <NA>
#> 1053       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1054       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1055       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1056   cohesion  consensus       plan       plan       <NA>       <NA>
#> 1057       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1058       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1059       plan       plan       plan       plan  consensus       plan
#> 1060       plan  consensus coregulate    monitor  consensus       plan
#> 1061       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1062       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1063       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1064       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1065       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1066       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1067       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1068    emotion  consensus coregulate    discuss  synthesis       <NA>
#> 1069       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1070       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1071       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1072  consensus coregulate  consensus       <NA>       <NA>       <NA>
#> 1073       plan       plan       plan  consensus       <NA>       <NA>
#> 1074       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1075       plan  consensus    discuss    emotion       plan       plan
#> 1076       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1077       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1078       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1079  consensus    emotion    emotion    emotion   cohesion    monitor
#> 1080  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1081       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1082    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1083    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1084       plan    discuss   cohesion    emotion       <NA>       <NA>
#> 1085       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1086    discuss    discuss    discuss      adapt  consensus coregulate
#> 1087       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1088    discuss      adapt    discuss  synthesis       <NA>       <NA>
#> 1089       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1090       plan  consensus       plan       plan       plan  consensus
#> 1091       plan  consensus    discuss    emotion       <NA>       <NA>
#> 1092  consensus       plan       plan  consensus coregulate    discuss
#> 1093  synthesis      adapt    discuss    emotion  consensus    monitor
#> 1094  consensus  consensus    monitor       plan       <NA>       <NA>
#> 1095       plan    emotion  consensus coregulate coregulate  consensus
#> 1096       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1097       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1098       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1099       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1100       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1101       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1102       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1103       plan    emotion   cohesion  consensus       plan  consensus
#> 1104   cohesion  consensus       plan    emotion   cohesion  consensus
#> 1105       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1106 coregulate    discuss    discuss  consensus       plan       plan
#> 1107   cohesion  consensus       plan  consensus       plan       plan
#> 1108       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1109    discuss      adapt    emotion    discuss coregulate    discuss
#> 1110    discuss    discuss    discuss  synthesis       <NA>       <NA>
#> 1111    discuss  synthesis  consensus  consensus       plan       plan
#> 1112       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1113       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1114       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1115    discuss  synthesis  consensus       <NA>       <NA>       <NA>
#> 1116    emotion   cohesion  consensus       plan       plan  consensus
#> 1117 coregulate       plan       plan  consensus    discuss  synthesis
#> 1118       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1119       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1120  synthesis  consensus    discuss  synthesis  consensus       plan
#> 1121       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1122       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1123    monitor    discuss    discuss    emotion       plan  consensus
#> 1124  consensus       plan       plan  consensus       plan   cohesion
#> 1125       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1126    discuss coregulate    discuss coregulate    discuss  consensus
#> 1127       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1128       plan       plan    emotion  consensus coregulate    emotion
#> 1129       plan       plan       plan  consensus    discuss    discuss
#> 1130   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1131      adapt  consensus  consensus       plan    emotion coregulate
#> 1132       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1133       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1134    discuss    monitor    discuss    discuss    discuss  consensus
#> 1135    discuss  synthesis coregulate    discuss   cohesion       plan
#> 1136       plan  consensus       plan       plan   cohesion    discuss
#> 1137 coregulate    discuss  synthesis    emotion       plan       plan
#> 1138  consensus    monitor    discuss      adapt       <NA>       <NA>
#> 1139      adapt coregulate      adapt  consensus       plan       plan
#> 1140  consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 1141    discuss  consensus    monitor  consensus coregulate    discuss
#> 1142    discuss   cohesion  consensus    discuss  consensus    monitor
#> 1143    discuss      adapt coregulate    emotion       plan  consensus
#> 1144       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1145       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1146       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1147    emotion    discuss  consensus       plan       plan       plan
#> 1148       plan  consensus       plan  consensus       plan coregulate
#> 1149  consensus       plan       plan  consensus       plan       plan
#> 1150       plan coregulate       plan       plan       plan    discuss
#> 1151  consensus       plan       plan       plan    emotion  consensus
#> 1152       plan   cohesion coregulate       <NA>       <NA>       <NA>
#> 1153 coregulate    monitor       plan  consensus coregulate    monitor
#> 1154 coregulate  consensus coregulate  synthesis      adapt  consensus
#> 1155    emotion  consensus  consensus    emotion       <NA>       <NA>
#> 1156       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1157       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1158    discuss coregulate       plan       plan  consensus    monitor
#> 1159       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1160       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1161       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1162       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1163  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1164       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1165  consensus       plan  consensus       <NA>       <NA>       <NA>
#> 1166   cohesion  consensus    monitor    emotion  consensus       <NA>
#> 1167  consensus       plan       plan    discuss  consensus       plan
#> 1168       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1169       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1170  consensus  synthesis      adapt   cohesion    discuss       <NA>
#> 1171    discuss      adapt    discuss    monitor    discuss  synthesis
#> 1172       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1173       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1174       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1175       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1176       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1177       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1178       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1179      adapt   cohesion  consensus coregulate    discuss    discuss
#> 1180       plan       plan       plan       plan       plan       plan
#> 1181  consensus       plan    monitor       <NA>       <NA>       <NA>
#> 1182       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1183       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1184    emotion  consensus       plan  consensus    emotion  consensus
#> 1185       plan    emotion   cohesion  consensus    monitor    discuss
#> 1186    discuss    discuss  synthesis    monitor    discuss    emotion
#> 1187  consensus  synthesis  consensus       plan       plan       plan
#> 1188       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1189 coregulate    emotion coregulate    discuss       <NA>       <NA>
#> 1190       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1191       plan  consensus       plan       plan       plan  consensus
#> 1192       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1193       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1194       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1195       plan    discuss    discuss       <NA>       <NA>       <NA>
#> 1196       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1197    monitor    discuss coregulate    monitor       plan       plan
#> 1198       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1199  synthesis    discuss  consensus       plan  consensus  consensus
#> 1200       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1201       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1202       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1203  synthesis   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1204 coregulate    emotion    monitor  consensus  synthesis      adapt
#> 1205       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1206       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1207    discuss    emotion  consensus       plan    emotion  consensus
#> 1208       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1209  consensus coregulate       plan      adapt  consensus coregulate
#> 1210    monitor    discuss    discuss       plan       plan  consensus
#> 1211       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1212       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1213   cohesion    discuss  consensus  consensus       plan       plan
#> 1214       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1215    discuss  synthesis  consensus       plan  consensus      adapt
#> 1216       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1217       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1218  consensus    monitor       <NA>       <NA>       <NA>       <NA>
#> 1219       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1220       plan    monitor    discuss    discuss   cohesion       plan
#> 1221       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1222       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1223       plan    monitor    emotion  consensus       plan       plan
#> 1224       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1225       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1226       plan       plan   cohesion    emotion  consensus coregulate
#> 1227  consensus coregulate    emotion  consensus       plan  consensus
#> 1228   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1229    discuss  consensus    discuss    monitor    monitor  consensus
#> 1230       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1231   cohesion coregulate       <NA>       <NA>       <NA>       <NA>
#> 1232       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1233    emotion   cohesion  consensus       plan       <NA>       <NA>
#> 1234  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1235  consensus    monitor  consensus       plan  consensus coregulate
#> 1236       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1237       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1238    discuss    discuss  consensus    monitor   cohesion   cohesion
#> 1239       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1240       plan       plan    monitor       plan    emotion    emotion
#> 1241       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1242    discuss       plan  consensus    emotion    discuss      adapt
#> 1243   cohesion    emotion   cohesion  consensus    monitor   cohesion
#> 1244       plan       plan    monitor    discuss    discuss    emotion
#> 1245  consensus       plan       plan       <NA>       <NA>       <NA>
#> 1246 coregulate  consensus       plan       plan    emotion  consensus
#> 1247    discuss coregulate    monitor    discuss  consensus       plan
#> 1248    emotion       plan  consensus       plan    emotion coregulate
#> 1249       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1250       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1251       plan  consensus       plan       plan  consensus coregulate
#> 1252       plan    emotion  consensus    monitor       <NA>       <NA>
#> 1253       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1254       plan       plan  consensus       plan       plan    monitor
#> 1255       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1256  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1257       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1258       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1259       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1260    emotion  consensus    discuss    discuss  synthesis       <NA>
#> 1261       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1262       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1263       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1264  synthesis      adapt       <NA>       <NA>       <NA>       <NA>
#> 1265       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1266       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1267       plan  consensus  synthesis  consensus    discuss      adapt
#> 1268    monitor       plan    monitor    discuss    discuss      adapt
#> 1269       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1270       plan       plan  consensus       plan   cohesion       plan
#> 1271 coregulate    emotion   cohesion    emotion       <NA>       <NA>
#> 1272    emotion   cohesion       plan  consensus coregulate    discuss
#> 1273       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1274 coregulate       plan    discuss  synthesis  consensus coregulate
#> 1275   cohesion  consensus coregulate       <NA>       <NA>       <NA>
#> 1276       plan    discuss       plan       plan  consensus       plan
#> 1277       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1278  consensus    discuss    monitor       plan       plan  consensus
#> 1279       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1280       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1281   cohesion  consensus    discuss    discuss  consensus coregulate
#> 1282       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1283       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1284       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1285       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1286    emotion    discuss    monitor    discuss    emotion   cohesion
#> 1287       plan       plan  consensus       plan  consensus       <NA>
#> 1288       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1289    discuss  synthesis  consensus    emotion  consensus       <NA>
#> 1290       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1291       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1292       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1293    emotion    discuss    emotion       <NA>       <NA>       <NA>
#> 1294  consensus    discuss    discuss    monitor    discuss coregulate
#> 1295    discuss  synthesis    discuss  synthesis   cohesion    discuss
#> 1296       plan  consensus  consensus       plan    discuss  synthesis
#> 1297  synthesis      adapt  consensus       plan  consensus coregulate
#> 1298       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1299  consensus       plan    emotion  consensus    discuss  consensus
#> 1300  consensus       plan    emotion    emotion   cohesion  consensus
#> 1301       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1302       plan    discuss    monitor    discuss    emotion  consensus
#> 1303       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1304       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1305    emotion    discuss    discuss    discuss coregulate       plan
#> 1306    discuss      adapt  consensus    monitor    discuss    emotion
#> 1307       plan    discuss coregulate       <NA>       <NA>       <NA>
#> 1308       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1309       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1310  synthesis  consensus       plan       <NA>       <NA>       <NA>
#> 1311       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1312       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1313       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1314       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1315    discuss  consensus       plan       plan  consensus coregulate
#> 1316       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1317       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1318       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1319  consensus coregulate    discuss coregulate    discuss       <NA>
#> 1320       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1321    monitor    discuss      adapt       <NA>       <NA>       <NA>
#> 1322       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1323       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1324  synthesis       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1325       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1326       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1327       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1328       plan  consensus       plan    discuss      adapt       <NA>
#> 1329  consensus       plan       plan       plan       plan    discuss
#> 1330    discuss    emotion   cohesion coregulate  consensus       plan
#> 1331       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1332       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1333       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1334  consensus    monitor    discuss    discuss      adapt coregulate
#> 1335       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1336    emotion   cohesion       plan coregulate    discuss coregulate
#> 1337       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1338       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1339       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1340       plan    discuss coregulate    monitor       <NA>       <NA>
#> 1341       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1342       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1343  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1344       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1345       plan  consensus       plan    monitor    discuss    discuss
#> 1346       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1347  synthesis  consensus       plan    discuss    emotion    discuss
#> 1348  consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 1349       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1350       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1351    monitor    discuss  synthesis    emotion  consensus    monitor
#> 1352  synthesis    emotion   cohesion coregulate       plan    emotion
#> 1353       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1354    discuss       plan    monitor    discuss    monitor    discuss
#> 1355       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1356   cohesion       plan    monitor    discuss    emotion    discuss
#> 1357      adapt  consensus       plan  consensus       plan       plan
#> 1358       plan       plan  consensus coregulate       plan       plan
#> 1359   cohesion    emotion    emotion  consensus coregulate       plan
#> 1360    emotion    discuss coregulate    emotion    emotion       <NA>
#> 1361    emotion    discuss    discuss    discuss       <NA>       <NA>
#> 1362  synthesis      adapt  consensus       plan  consensus    discuss
#> 1363    discuss      adapt  consensus coregulate    discuss  synthesis
#> 1364       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1365       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1366       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1367       plan    emotion  consensus    monitor  consensus  consensus
#> 1368       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1369  consensus coregulate    discuss  consensus    emotion    discuss
#> 1370       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1371 coregulate       plan    discuss coregulate    discuss      adapt
#> 1372       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1373       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1374       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1375       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1376    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1377       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1378       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1379       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1380   cohesion  consensus       plan       plan   cohesion  consensus
#> 1381 coregulate    monitor    discuss    monitor       <NA>       <NA>
#> 1382    emotion       plan    emotion       plan    discuss      adapt
#> 1383       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1384    discuss      adapt   cohesion  consensus coregulate       plan
#> 1385   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1386       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1387       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1388       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1389    discuss  consensus  consensus       plan  consensus    discuss
#> 1390       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1391       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1392       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1393       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1394       plan  consensus    discuss       <NA>       <NA>       <NA>
#> 1395       plan  consensus    monitor       <NA>       <NA>       <NA>
#> 1396       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1397       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1398       plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 1399 coregulate    discuss  synthesis  consensus       <NA>       <NA>
#> 1400  synthesis    emotion   cohesion    monitor       plan    emotion
#> 1401  consensus  consensus coregulate       plan  consensus    emotion
#> 1402       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1403  consensus       plan       plan    monitor    discuss      adapt
#> 1404 coregulate  consensus       <NA>       <NA>       <NA>       <NA>
#> 1405       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1406       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1407       plan       plan  consensus    emotion  consensus       plan
#> 1408       plan       plan    discuss    discuss    discuss  synthesis
#> 1409       plan       plan  consensus  consensus       plan       plan
#> 1410       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1411       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1412       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1413  consensus       plan       plan  consensus       plan       <NA>
#> 1414       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1415  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1416       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1417       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1418       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1419    discuss    discuss coregulate       <NA>       <NA>       <NA>
#> 1420       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1421       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1422       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1423       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1424  consensus coregulate       plan       plan       plan    discuss
#> 1425       plan    monitor       plan  consensus       plan    emotion
#> 1426  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1427       plan       plan   cohesion       plan       plan       <NA>
#> 1428       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1429       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1430       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1431       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1432  consensus coregulate    monitor    discuss  synthesis  consensus
#> 1433 coregulate  consensus    discuss      adapt    emotion   cohesion
#> 1434       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1435       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1436       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1437    monitor    discuss  synthesis       <NA>       <NA>       <NA>
#> 1438       plan       plan  consensus    monitor  consensus  consensus
#> 1439    emotion coregulate    emotion    discuss    emotion       <NA>
#> 1440    emotion       plan    monitor       <NA>       <NA>       <NA>
#> 1441       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1442       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1443    discuss coregulate       <NA>       <NA>       <NA>       <NA>
#> 1444       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1445    emotion   cohesion coregulate       <NA>       <NA>       <NA>
#> 1446    discuss  synthesis coregulate    emotion    discuss    discuss
#> 1447 coregulate    emotion       plan       plan  consensus       plan
#> 1448    discuss coregulate    discuss  synthesis      adapt  consensus
#> 1449       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1450  consensus       plan  consensus    emotion       <NA>       <NA>
#> 1451  consensus       plan    emotion    emotion   cohesion       plan
#> 1452       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1453  consensus  consensus    monitor    discuss  consensus       plan
#> 1454 coregulate    discuss  consensus       plan       plan       plan
#> 1455       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1456    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1457       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1458       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1459       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1460       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1461       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1462       plan       plan       plan    monitor  consensus    discuss
#> 1463       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1464       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1465       plan    emotion   cohesion       plan       plan       plan
#> 1466       plan  consensus coregulate    discuss  consensus    monitor
#> 1467    monitor  consensus    monitor coregulate       plan  consensus
#> 1468  consensus    monitor  synthesis    emotion  consensus       plan
#> 1469       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1470  consensus coregulate  consensus       plan    emotion       <NA>
#> 1471       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1472  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1473 coregulate    monitor    emotion  consensus coregulate    discuss
#> 1474    discuss    discuss      adapt  consensus coregulate       plan
#> 1475    emotion       plan       plan       plan  consensus       plan
#> 1476       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1477    monitor    discuss    discuss  synthesis      adapt   cohesion
#> 1478       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1479       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1480       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1481       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1482  consensus       plan    emotion coregulate    emotion   cohesion
#> 1483       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1484       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1485       plan    emotion  consensus       plan  consensus coregulate
#> 1486       plan  consensus  consensus    discuss    discuss  consensus
#> 1487       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1488       plan  consensus    emotion  consensus    emotion    discuss
#> 1489  consensus       plan    discuss   cohesion  consensus    discuss
#> 1490       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1491       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1492       plan  consensus coregulate    discuss  synthesis      adapt
#> 1493       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1494   cohesion       plan       <NA>       <NA>       <NA>       <NA>
#> 1495    emotion  consensus       plan       plan  consensus       plan
#> 1496       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1497    discuss  consensus       plan  consensus       plan    emotion
#> 1498       plan       plan    emotion    emotion    emotion    emotion
#> 1499  synthesis  consensus       plan       <NA>       <NA>       <NA>
#> 1500    discuss    discuss  consensus    emotion    discuss      adapt
#> 1501       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1502       plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 1503  consensus  consensus coregulate    discuss  consensus       plan
#> 1504       plan    emotion    emotion    monitor       plan    emotion
#> 1505       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1506       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1507   cohesion  consensus       plan  consensus coregulate    discuss
#> 1508       plan    emotion    emotion    emotion   cohesion  consensus
#> 1509    discuss    discuss    emotion   cohesion coregulate  consensus
#> 1510       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1511       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1512  consensus    emotion  consensus       plan       plan       <NA>
#> 1513       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1514       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1515    discuss    monitor       plan    emotion       plan    monitor
#> 1516    emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1517       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1518 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1519  consensus    emotion   cohesion       plan  consensus    discuss
#> 1520    discuss    discuss      adapt       <NA>       <NA>       <NA>
#> 1521       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1522       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1523  consensus       plan  consensus    emotion   cohesion  consensus
#> 1524       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1525  consensus       plan       plan  consensus    monitor       plan
#> 1526    monitor  consensus       plan   cohesion    emotion    emotion
#> 1527       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1528       plan  consensus    discuss  synthesis  consensus  consensus
#> 1529       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1530    discuss    discuss  consensus   cohesion       <NA>       <NA>
#> 1531       plan       plan    emotion  consensus    discuss    emotion
#> 1532    emotion    discuss  consensus       plan       plan       plan
#> 1533       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1534    discuss  synthesis      adapt  consensus coregulate       plan
#> 1535       plan       plan       plan  consensus       plan   cohesion
#> 1536    emotion   cohesion  consensus       plan    emotion   cohesion
#> 1537       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1538       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1539       plan       plan  consensus    discuss  consensus    monitor
#> 1540       plan       plan    emotion  consensus       plan       <NA>
#> 1541       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1542  consensus coregulate       plan  consensus    discuss  synthesis
#> 1543       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1544  consensus       plan       plan       <NA>       <NA>       <NA>
#> 1545       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1546    discuss    emotion   cohesion  consensus       <NA>       <NA>
#> 1547       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1548       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1549       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1550  consensus       plan  consensus       plan   cohesion  consensus
#> 1551       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1552       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1553       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1554       plan  consensus       plan       plan  consensus coregulate
#> 1555       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1556       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1557 coregulate       plan    monitor       plan  consensus       plan
#> 1558       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1559    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1560       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1561  consensus    discuss    discuss   cohesion  consensus       plan
#> 1562       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1563       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1564       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1565       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1566       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1567    emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1568  consensus    emotion  consensus       plan    emotion  consensus
#> 1569    emotion    discuss  consensus coregulate    discuss    discuss
#> 1570 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1571       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1572       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1573   cohesion coregulate       <NA>       <NA>       <NA>       <NA>
#> 1574       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1575  consensus       plan    emotion   cohesion  consensus       plan
#> 1576       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1577  consensus       plan       plan       plan    discuss      adapt
#> 1578    discuss    discuss    discuss    discuss    discuss    discuss
#> 1579       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1580       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1581   cohesion    emotion  consensus    discuss    emotion   cohesion
#> 1582       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1583       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1584       plan  consensus  consensus    monitor    discuss       plan
#> 1585       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1586 coregulate    monitor       plan       <NA>       <NA>       <NA>
#> 1587       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1588       plan    discuss  synthesis      adapt  consensus coregulate
#> 1589       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1590  consensus       plan       plan       plan    emotion    discuss
#> 1591       plan       plan    discuss coregulate    discuss  consensus
#> 1592       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1593       plan       plan  consensus    discuss  consensus    monitor
#> 1594    emotion  consensus    monitor       <NA>       <NA>       <NA>
#> 1595    discuss      adapt    emotion       <NA>       <NA>       <NA>
#> 1596      adapt  consensus       <NA>       <NA>       <NA>       <NA>
#> 1597       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1598       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1599       plan       plan       plan    emotion   cohesion       plan
#> 1600       plan       plan  consensus  consensus       plan       plan
#> 1601       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1602    emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 1603       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1604  consensus       plan  consensus       <NA>       <NA>       <NA>
#> 1605    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1606       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1607  consensus coregulate  consensus coregulate    discuss  consensus
#> 1608       plan    discuss coregulate    discuss       <NA>       <NA>
#> 1609       plan    monitor    discuss  synthesis    discuss       plan
#> 1610       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1611    discuss  synthesis coregulate       plan       plan       plan
#> 1612       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1613       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1614       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1615       plan  consensus    monitor    discuss  consensus       plan
#> 1616       plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 1617       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1618       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1619       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1620       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1621       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1622  consensus       plan       plan  consensus       plan    monitor
#> 1623       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1624       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1625       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1626       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1627       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1628       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1629    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1630       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1631   cohesion       plan  consensus  consensus    discuss       <NA>
#> 1632  consensus    emotion    emotion   cohesion  consensus       plan
#> 1633       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1634      adapt       plan    discuss  consensus       plan       plan
#> 1635    discuss  consensus coregulate  consensus       plan    emotion
#> 1636       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1637    emotion       plan  consensus       plan  consensus       plan
#> 1638       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1639       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1640    emotion       plan  consensus   cohesion coregulate   cohesion
#> 1641       plan       plan       plan    emotion    monitor    emotion
#> 1642       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1643       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1644 coregulate       plan       plan       plan       <NA>       <NA>
#> 1645  consensus       plan    emotion    monitor    monitor    discuss
#> 1646       plan       plan       plan  consensus    discuss    discuss
#> 1647       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1648    discuss      adapt    discuss coregulate       plan       plan
#> 1649       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1650       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1651       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1652  consensus    monitor    discuss  synthesis       <NA>       <NA>
#> 1653       plan       plan       plan    emotion    emotion   cohesion
#> 1654    discuss  consensus       plan       plan  consensus coregulate
#> 1655       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1656       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1657    emotion       plan  consensus       <NA>       <NA>       <NA>
#> 1658       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1659       plan       plan       plan       plan    monitor    discuss
#> 1660       plan  consensus    discuss  synthesis      adapt    monitor
#> 1661       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1662    discuss    discuss    discuss      adapt    emotion    emotion
#> 1663       plan    discuss    discuss  consensus coregulate    discuss
#> 1664   cohesion    discuss coregulate    monitor   cohesion    emotion
#> 1665       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1666       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1667       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1668       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1669       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1670       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1671  consensus       plan       plan  consensus       plan  consensus
#> 1672       plan       plan    emotion       plan  consensus       plan
#> 1673       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1674       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1675    discuss  synthesis      adapt   cohesion  consensus    emotion
#> 1676       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1677       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1678       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1679       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1680       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1681  consensus    discuss    monitor    discuss coregulate    emotion
#> 1682       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1683  synthesis       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1684       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1685       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1686       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1687       plan  consensus    discuss coregulate       plan    discuss
#> 1688       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1689       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1690       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1691  consensus       plan       plan    monitor       plan       plan
#> 1692       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1693       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1694  consensus coregulate  consensus       plan       plan  consensus
#> 1695       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1696 coregulate    emotion       plan    emotion  consensus coregulate
#> 1697    emotion  consensus       plan  consensus       plan    monitor
#> 1698       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1699    monitor    emotion       plan   cohesion  consensus       plan
#> 1700  synthesis  consensus    discuss    emotion   cohesion    emotion
#> 1701       plan  consensus       plan       <NA>       <NA>       <NA>
#> 1702       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1703 coregulate       plan       plan    discuss       <NA>       <NA>
#> 1704       plan    discuss    discuss  synthesis  consensus       <NA>
#> 1705       plan    emotion   cohesion       <NA>       <NA>       <NA>
#> 1706    discuss  synthesis       <NA>       <NA>       <NA>       <NA>
#> 1707       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1708       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1709       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1710       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1711  consensus  consensus       plan       <NA>       <NA>       <NA>
#> 1712       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1713   cohesion  consensus coregulate  consensus       plan  consensus
#> 1714       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1715       plan       plan  consensus    monitor    discuss  consensus
#> 1716  consensus       plan       plan       <NA>       <NA>       <NA>
#> 1717       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1718   cohesion  consensus  consensus      adapt  consensus       plan
#> 1719 coregulate       plan    emotion  consensus       plan       plan
#> 1720       plan    discuss  consensus       plan    emotion  consensus
#> 1721       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1722   cohesion    emotion   cohesion  consensus       plan    emotion
#> 1723       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1724       plan  consensus coregulate  consensus       plan       plan
#> 1725       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1726       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1727    emotion  consensus  synthesis  consensus    discuss    discuss
#> 1728       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1729       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1730       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1731       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1732       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1733       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1734    emotion   cohesion  consensus       plan    monitor  consensus
#> 1735 coregulate    discuss      adapt  consensus       plan   cohesion
#> 1736       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1737       plan       plan    emotion  consensus       <NA>       <NA>
#> 1738       plan  consensus coregulate       <NA>       <NA>       <NA>
#> 1739       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1740       plan       plan       plan       plan  consensus coregulate
#> 1741       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1742 coregulate    emotion   cohesion       <NA>       <NA>       <NA>
#> 1743       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1744       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1745       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1746       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1747  consensus coregulate  consensus       <NA>       <NA>       <NA>
#> 1748 coregulate    emotion  consensus       plan  consensus coregulate
#> 1749       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1750       plan       plan       plan  consensus    monitor    discuss
#> 1751       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1752       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1753       plan    emotion   cohesion  consensus  consensus    emotion
#> 1754       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1755       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1756       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1757       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1758       plan coregulate    discuss    discuss coregulate       plan
#> 1759    discuss      adapt       <NA>       <NA>       <NA>       <NA>
#> 1760       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1761    monitor  consensus       plan  consensus coregulate    discuss
#> 1762       plan  consensus    monitor       <NA>       <NA>       <NA>
#> 1763  consensus       plan  consensus       plan  consensus       plan
#> 1764       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1765       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1766       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1767 coregulate    emotion   cohesion coregulate    discuss    emotion
#> 1768       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1769       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1770       plan    discuss  consensus       plan    emotion   cohesion
#> 1771    discuss  synthesis  consensus coregulate    discuss    discuss
#> 1772       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1773  synthesis  consensus       plan       plan       plan       plan
#> 1774       plan    discuss    discuss  synthesis  consensus    discuss
#> 1775       plan    emotion coregulate  consensus  consensus       plan
#> 1776  consensus       plan       plan  consensus coregulate coregulate
#> 1777       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1778       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1779  consensus       plan       plan       plan       plan  consensus
#> 1780       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1781       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1782       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1783  consensus coregulate       plan       <NA>       <NA>       <NA>
#> 1784       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1785       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1786 coregulate  synthesis  consensus       plan       plan       plan
#> 1787       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1788 coregulate    discuss  synthesis      adapt       <NA>       <NA>
#> 1789       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1790       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1791    monitor    discuss coregulate    discuss  consensus       plan
#> 1792       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1793       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1794       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1795       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1796       plan    monitor  consensus coregulate       plan       plan
#> 1797       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1798    discuss    discuss  consensus coregulate    emotion    discuss
#> 1799       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1800  consensus       plan       plan       plan       plan  consensus
#> 1801       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1802       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1803       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1804       plan  consensus coregulate coregulate       plan       plan
#> 1805       plan  consensus       plan       plan    discuss  synthesis
#> 1806    discuss    discuss    emotion   cohesion  consensus    discuss
#> 1807       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1808       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1809  consensus  consensus       <NA>       <NA>       <NA>       <NA>
#> 1810       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1811       plan       plan       plan       plan    discuss    discuss
#> 1812    monitor       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1813    monitor    emotion   cohesion    discuss coregulate    discuss
#> 1814  consensus coregulate    discuss  consensus       plan  consensus
#> 1815  synthesis    discuss    emotion coregulate       <NA>       <NA>
#> 1816       plan       plan       plan    discuss  consensus    discuss
#> 1817  consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 1818       plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 1819       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1820       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1821       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1822       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1823       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1824  consensus coregulate    emotion    monitor  consensus    discuss
#> 1825       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1826       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1827    discuss  consensus  consensus       plan    emotion   cohesion
#> 1828   cohesion    discuss coregulate    discuss  consensus coregulate
#> 1829   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1830       plan  consensus coregulate    discuss      adapt    discuss
#> 1831       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1832  consensus       plan  consensus       <NA>       <NA>       <NA>
#> 1833       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1834       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1835       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1836    emotion   cohesion    emotion   cohesion       plan  consensus
#> 1837 coregulate    discuss    monitor   cohesion       plan  consensus
#> 1838  synthesis  consensus      adapt       <NA>       <NA>       <NA>
#> 1839       plan  consensus  consensus       <NA>       <NA>       <NA>
#> 1840 coregulate  consensus coregulate       plan  consensus       plan
#> 1841    discuss   cohesion coregulate    emotion       <NA>       <NA>
#> 1842 coregulate    discuss  synthesis      adapt    emotion  consensus
#> 1843       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1844    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1845  synthesis      adapt  consensus    discuss    discuss       <NA>
#> 1846      adapt  consensus       plan  consensus       plan       plan
#> 1847       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1848  consensus    discuss    discuss  synthesis      adapt  consensus
#> 1849       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1850       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1851       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1852 coregulate  consensus coregulate    emotion  consensus    discuss
#> 1853       plan    emotion       plan       plan       plan    emotion
#> 1854       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1855    emotion   cohesion  consensus       plan  consensus coregulate
#> 1856       plan       plan  consensus       plan    emotion  consensus
#> 1857       plan       plan    emotion  consensus coregulate  consensus
#> 1858       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1859       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1860       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1861       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1862       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1863 coregulate    emotion       <NA>       <NA>       <NA>       <NA>
#> 1864  consensus    emotion    emotion  consensus       plan  consensus
#> 1865       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1866  consensus    discuss coregulate    discuss    discuss    discuss
#> 1867    emotion    discuss  consensus coregulate    monitor    discuss
#> 1868  consensus    discuss    discuss  synthesis       <NA>       <NA>
#> 1869       plan  consensus       plan       <NA>       <NA>       <NA>
#> 1870       plan  consensus coregulate    discuss  synthesis      adapt
#> 1871       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1872    discuss    emotion   cohesion    emotion      adapt       <NA>
#> 1873       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1874       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1875    discuss  consensus       plan       plan  consensus       plan
#> 1876      adapt  consensus  consensus       plan       <NA>       <NA>
#> 1877   cohesion  consensus    discuss      adapt       <NA>       <NA>
#> 1878       plan coregulate       plan       plan    emotion    emotion
#> 1879  consensus       plan    discuss    discuss  synthesis  consensus
#> 1880       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1881 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1882 coregulate   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1883       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1884       plan  consensus coregulate    emotion       plan       <NA>
#> 1885    discuss  consensus       plan       plan       <NA>       <NA>
#> 1886    emotion       plan       plan       <NA>       <NA>       <NA>
#> 1887       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1888       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1889       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1890  consensus    discuss  synthesis    emotion  consensus    discuss
#> 1891       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1892       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1893       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1894       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1895    discuss    discuss  consensus    discuss       <NA>       <NA>
#> 1896       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1897       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1898  consensus coregulate    monitor       plan  consensus    monitor
#> 1899   cohesion    monitor    discuss  consensus coregulate       plan
#> 1900 coregulate       plan  consensus coregulate  consensus       plan
#> 1901       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1902    discuss  synthesis  consensus       plan       plan       plan
#> 1903       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1904   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1905 coregulate       plan       <NA>       <NA>       <NA>       <NA>
#> 1906       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1907  consensus    emotion  consensus    discuss  synthesis    discuss
#> 1908       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1909       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1910  synthesis       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1911  consensus       plan    emotion    discuss   cohesion  consensus
#> 1912       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1913       plan  consensus    emotion   cohesion       plan       <NA>
#> 1914    discuss    monitor       plan  consensus coregulate    discuss
#> 1915       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1916       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1917       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1918       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1919       plan       plan    emotion   cohesion    emotion    monitor
#> 1920       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1921  consensus       plan    emotion    emotion  consensus       plan
#> 1922       plan       plan    emotion   cohesion    emotion   cohesion
#> 1923       plan coregulate    discuss    emotion   cohesion    monitor
#> 1924       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1925  consensus coregulate       plan       plan  consensus       plan
#> 1926       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1927       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1928       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1929       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1930    discuss    emotion   cohesion  consensus coregulate  consensus
#> 1931       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1932       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1933       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1934       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1935    discuss coregulate  consensus    monitor    discuss  consensus
#> 1936   cohesion  consensus    emotion  consensus       plan    discuss
#> 1937       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1938       plan  consensus coregulate    discuss  consensus coregulate
#> 1939       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1940       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1941       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1942       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1943       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1944 coregulate       plan       plan    monitor  synthesis    discuss
#> 1945       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1946       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1947       plan  consensus       plan       plan  consensus       plan
#> 1948       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1949       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1950       plan    discuss       <NA>       <NA>       <NA>       <NA>
#> 1951       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1952       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1953       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1954       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1955       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1956  consensus       plan  consensus       plan       plan    emotion
#> 1957    emotion       plan    discuss      adapt    emotion    emotion
#> 1958  consensus       plan    discuss  consensus       plan  consensus
#> 1959       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1960       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1961       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1962  consensus       plan    monitor    emotion   cohesion  consensus
#> 1963       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1964 coregulate    discuss coregulate       <NA>       <NA>       <NA>
#> 1965    emotion    emotion    emotion  consensus    monitor    discuss
#> 1966    discuss      adapt    discuss  consensus coregulate       <NA>
#> 1967       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1968       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1969       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1970       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1971       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1972       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1973       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1974  consensus       plan       plan coregulate    discuss  consensus
#> 1975    discuss    discuss  synthesis      adapt  consensus coregulate
#> 1976    emotion    emotion   cohesion    discuss    discuss  consensus
#> 1977       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1978       plan    monitor   cohesion       <NA>       <NA>       <NA>
#> 1979  consensus coregulate    emotion    discuss    emotion       plan
#> 1980    discuss      adapt       plan    discuss       <NA>       <NA>
#> 1981       plan  consensus    monitor    discuss  synthesis       <NA>
#> 1982       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1983       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1984       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1985       plan       plan    monitor       <NA>       <NA>       <NA>
#> 1986       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1987       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1988    discuss       plan    emotion    discuss  consensus coregulate
#> 1989       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1990       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1991    discuss  consensus       plan       plan       <NA>       <NA>
#> 1992       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1993       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1994  consensus    emotion   cohesion  consensus       plan  consensus
#> 1995 coregulate  consensus       <NA>       <NA>       <NA>       <NA>
#> 1996  consensus coregulate    discuss    discuss    emotion  consensus
#> 1997    discuss      adapt   cohesion coregulate       plan       <NA>
#> 1998       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1999       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 2000       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#>             T20        T21        T22        T23        T24        T25
#> 1          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 2          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 3       monitor  consensus       plan    emotion       plan  consensus
#> 4          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 5          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 6          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 7          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 8       emotion   cohesion  consensus    emotion       plan       plan
#> 9          <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 10         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 11         plan       plan    monitor       <NA>       <NA>       <NA>
#> 12   coregulate    emotion   cohesion       plan  consensus       <NA>
#> 13         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 14         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 15         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 16    consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 17         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 18         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 19         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 20         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 21    consensus       plan       plan  consensus       plan       <NA>
#> 22         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 23         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 24         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 25         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 26         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 27         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 28         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 29         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 30     cohesion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 31   coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 32      emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 33         plan    emotion   cohesion  consensus       <NA>       <NA>
#> 34         plan       plan  consensus    discuss  consensus       <NA>
#> 35   coregulate   cohesion  consensus coregulate       <NA>       <NA>
#> 36         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 37      emotion  consensus       plan       <NA>       <NA>       <NA>
#> 38         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 39      discuss    discuss    discuss       <NA>       <NA>       <NA>
#> 40         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 41      emotion    discuss  consensus coregulate       plan    monitor
#> 42    consensus    discuss  consensus       <NA>       <NA>       <NA>
#> 43         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 44         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 45     cohesion  consensus    emotion    emotion       <NA>       <NA>
#> 46         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 47         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 48      emotion  consensus   cohesion    emotion   cohesion  consensus
#> 49         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 50         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 51         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 52         plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 53    consensus    emotion  consensus       <NA>       <NA>       <NA>
#> 54         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 55         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 56         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 57         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 58         plan  consensus       plan       plan  consensus       <NA>
#> 59      discuss  synthesis  consensus    discuss  consensus       <NA>
#> 60         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 61         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 62         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 63      emotion    discuss       <NA>       <NA>       <NA>       <NA>
#> 64         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 65         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 66         plan    discuss       <NA>       <NA>       <NA>       <NA>
#> 67      discuss  consensus   cohesion  consensus       plan       plan
#> 68         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 69         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 70         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 71         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 72         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 73         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 74         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 75         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 76         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 77         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 78         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 79         plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 80         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 81         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 82         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 83      emotion       plan  consensus       plan  consensus       plan
#> 84     cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 85   coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 86         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 87         plan       plan       <NA>       <NA>       <NA>       <NA>
#> 88         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 89         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 90         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 91         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 92      discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 93         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 94         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 95    consensus coregulate       plan       <NA>       <NA>       <NA>
#> 96      discuss  synthesis      adapt  consensus coregulate    emotion
#> 97         <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 98      monitor   cohesion       <NA>       <NA>       <NA>       <NA>
#> 99     cohesion       plan  consensus    discuss    discuss    emotion
#> 100        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 101        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 102        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 103        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 104        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 105     emotion  consensus    monitor       plan    emotion       <NA>
#> 106        plan       plan      adapt  consensus       plan       <NA>
#> 107        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 108        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 109        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 110        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 111        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 112   consensus coregulate    discuss       <NA>       <NA>       <NA>
#> 113        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 114        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 115     discuss  synthesis  consensus    emotion       plan  consensus
#> 116        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 117        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 118        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 119        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 120        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 121        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 122        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 123        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 124        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 125     discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 126        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 127     emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 128     discuss  consensus       plan   cohesion  consensus    discuss
#> 129        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 130   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 131        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 132     discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 133        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 134   consensus       plan  consensus   cohesion    emotion   cohesion
#> 135        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 136        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 137        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 138    cohesion    emotion    discuss  consensus    discuss  consensus
#> 139        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 140        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 141        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 142  coregulate       plan    emotion   cohesion  consensus       <NA>
#> 143        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 144        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 145        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 146     monitor    discuss  consensus coregulate    emotion       <NA>
#> 147        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 148        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 149        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 150        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 151        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 152     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 153        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 154   consensus    discuss  consensus       <NA>       <NA>       <NA>
#> 155        plan       plan    monitor    discuss    emotion       <NA>
#> 156        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 157        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 158        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 159        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 160        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 161     emotion    discuss    discuss    emotion  consensus       plan
#> 162        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 163        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 164        plan   cohesion coregulate    discuss   cohesion  consensus
#> 165        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 166        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 167        plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 168        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 169        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 170        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 171        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 172        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 173        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 174        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 175        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 176   consensus coregulate    emotion       plan    monitor    discuss
#> 177        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 178        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 179        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 180        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 181        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 182        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 183        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 184        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 185        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 186        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 187        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 188     emotion  consensus       plan  consensus       plan    emotion
#> 189        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 190        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 191        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 192        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 193        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 194        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 195        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 196        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 197     emotion  consensus       plan       plan       <NA>       <NA>
#> 198        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 199        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 200        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 201        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 202        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 203        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 204   consensus    discuss  synthesis  consensus       <NA>       <NA>
#> 205        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 206        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 207        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 208        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 209        plan    emotion   cohesion       <NA>       <NA>       <NA>
#> 210        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 211     emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 212  coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 213        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 214     monitor    discuss   cohesion       <NA>       <NA>       <NA>
#> 215        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 216        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 217        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 218        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 219        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 220        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 221        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 222        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 223        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 224    cohesion  consensus       <NA>       <NA>       <NA>       <NA>
#> 225     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 226        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 227     discuss  consensus    discuss    discuss    emotion  consensus
#> 228   consensus    discuss    discuss  synthesis       <NA>       <NA>
#> 229        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 230        plan  consensus    emotion       <NA>       <NA>       <NA>
#> 231     emotion  consensus       plan       plan       <NA>       <NA>
#> 232        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 233     discuss  consensus    discuss    emotion  consensus       plan
#> 234    cohesion  consensus    discuss  synthesis       plan    emotion
#> 235   synthesis       plan    monitor    discuss  synthesis       <NA>
#> 236        plan    emotion       plan  consensus       plan    emotion
#> 237        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 238        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 239        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 240    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 241   consensus   cohesion       <NA>       <NA>       <NA>       <NA>
#> 242     emotion  consensus coregulate    emotion  consensus       plan
#> 243        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 244        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 245        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 246        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 247     discuss    emotion       plan  consensus    discuss  consensus
#> 248        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 249        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 250   consensus   cohesion       <NA>       <NA>       <NA>       <NA>
#> 251        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 252        plan       plan       <NA>       <NA>       <NA>       <NA>
#> 253        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 254        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 255        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 256        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 257        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 258        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 259        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 260        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 261        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 262        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 263        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 264        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 265        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 266   consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 267        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 268        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 269   consensus       plan  consensus    emotion       <NA>       <NA>
#> 270        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 271        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 272  coregulate       plan  consensus       <NA>       <NA>       <NA>
#> 273        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 274        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 275     emotion    discuss   cohesion  consensus       plan    discuss
#> 276        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 277   consensus    monitor       <NA>       <NA>       <NA>       <NA>
#> 278        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 279        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 280        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 281        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 282        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 283        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 284     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 285        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 286        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 287        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 288     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 289        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 290        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 291        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 292    cohesion  consensus coregulate    monitor   cohesion       plan
#> 293   consensus    discuss    emotion  consensus coregulate    discuss
#> 294        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 295   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 296        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 297        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 298        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 299    cohesion  consensus       plan  consensus       plan       <NA>
#> 300        plan       plan    emotion  consensus       plan       <NA>
#> 301        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 302        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 303        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 304        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 305        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 306        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 307        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 308        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 309        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 310        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 311        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 312  coregulate       plan       plan  consensus       <NA>       <NA>
#> 313        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 314        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 315        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 316        plan  consensus    discuss    discuss       <NA>       <NA>
#> 317        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 318   consensus       plan       plan  consensus       plan       plan
#> 319        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 320        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 321        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 322        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 323        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 324        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 325        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 326        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 327        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 328   consensus  consensus coregulate       plan       <NA>       <NA>
#> 329     discuss  consensus    monitor       <NA>       <NA>       <NA>
#> 330        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 331        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 332        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 333        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 334   consensus       plan  consensus       <NA>       <NA>       <NA>
#> 335        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 336        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 337     monitor    discuss    discuss  consensus coregulate       <NA>
#> 338     discuss    emotion       plan  consensus    emotion       <NA>
#> 339     monitor    emotion       <NA>       <NA>       <NA>       <NA>
#> 340   consensus       plan  consensus       plan       <NA>       <NA>
#> 341        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 342        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 343   consensus       plan       plan       plan       <NA>       <NA>
#> 344        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 345        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 346        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 347        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 348        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 349        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 350        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 351     discuss    discuss       plan       plan    monitor       <NA>
#> 352        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 353        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 354     emotion  consensus       plan       <NA>       <NA>       <NA>
#> 355        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 356        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 357   consensus       plan       plan  consensus       plan       <NA>
#> 358        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 359        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 360        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 361     emotion       plan       plan  consensus coregulate    discuss
#> 362        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 363        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 364        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 365        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 366        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 367        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 368        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 369        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 370        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 371        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 372        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 373        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 374  coregulate       plan    discuss  synthesis       <NA>       <NA>
#> 375   consensus    monitor       plan       <NA>       <NA>       <NA>
#> 376        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 377        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 378        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 379        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 380        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 381  coregulate    monitor       plan       plan  consensus       plan
#> 382   consensus  consensus       plan       plan   cohesion       <NA>
#> 383        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 384   consensus       plan       plan       plan       <NA>       <NA>
#> 385        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 386        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 387        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 388   consensus    discuss  consensus coregulate      adapt       <NA>
#> 389        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 390     discuss    emotion   cohesion  consensus       plan       <NA>
#> 391        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 392        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 393     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 394        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 395        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 396        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 397        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 398        plan    discuss coregulate  consensus coregulate       <NA>
#> 399        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 400        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 401        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 402     emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 403        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 404     emotion  consensus    discuss    discuss       <NA>       <NA>
#> 405        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 406        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 407   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 408        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 409        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 410        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 411        plan  consensus       plan coregulate       <NA>       <NA>
#> 412        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 413        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 414        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 415        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 416    cohesion    emotion  consensus   cohesion  consensus    discuss
#> 417        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 418        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 419     emotion  consensus    emotion  consensus coregulate       <NA>
#> 420     emotion   cohesion       plan    monitor       <NA>       <NA>
#> 421        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 422   consensus  consensus    monitor    discuss       <NA>       <NA>
#> 423        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 424   consensus    discuss  synthesis       <NA>       <NA>       <NA>
#> 425        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 426        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 427   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 428        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 429        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 430        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 431    cohesion  consensus  consensus       plan  consensus       <NA>
#> 432        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 433        plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 434        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 435        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 436     discuss  consensus    discuss    discuss       <NA>       <NA>
#> 437   consensus coregulate    emotion    emotion       <NA>       <NA>
#> 438        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 439        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 440        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 441        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 442        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 443   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 444     emotion    discuss  consensus    discuss  synthesis  consensus
#> 445        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 446        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 447        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 448        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 449        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 450        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 451        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 452        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 453        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 454        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 455        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 456        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 457   consensus       plan       plan       plan       <NA>       <NA>
#> 458        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 459        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 460        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 461        plan       plan  consensus       plan  consensus       plan
#> 462        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 463        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 464        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 465        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 466        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 467        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 468        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 469        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 470        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 471  coregulate    emotion    emotion   cohesion  consensus       plan
#> 472   synthesis  consensus       plan   cohesion       plan    monitor
#> 473        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 474     emotion   cohesion coregulate    discuss    emotion    discuss
#> 475   consensus       plan  consensus    discuss    discuss coregulate
#> 476        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 477        plan       plan   cohesion  consensus  consensus  consensus
#> 478     monitor    discuss    discuss  consensus       <NA>       <NA>
#> 479     discuss    monitor  consensus    discuss  synthesis  consensus
#> 480   consensus    emotion       <NA>       <NA>       <NA>       <NA>
#> 481        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 482  coregulate       plan  consensus       plan    monitor       <NA>
#> 483        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 484        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 485        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 486        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 487   consensus       plan       plan       plan       plan    emotion
#> 488        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 489        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 490        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 491        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 492   consensus       plan    emotion    discuss       <NA>       <NA>
#> 493        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 494     discuss  consensus    monitor    discuss       <NA>       <NA>
#> 495        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 496        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 497        plan  consensus    monitor       <NA>       <NA>       <NA>
#> 498        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 499        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 500        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 501        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 502        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 503        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 504        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 505        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 506        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 507   consensus       plan    monitor    discuss   cohesion       <NA>
#> 508        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 509        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 510        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 511     discuss  consensus    discuss    emotion       plan    emotion
#> 512  coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 513        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 514        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 515        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 516        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 517        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 518     emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 519        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 520        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 521   consensus    discuss  consensus    monitor       plan       <NA>
#> 522        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 523        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 524        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 525        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 526        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 527        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 528        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 529        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 530        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 531        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 532        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 533        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 534     discuss    discuss  synthesis       <NA>       <NA>       <NA>
#> 535        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 536   synthesis      adapt  consensus       <NA>       <NA>       <NA>
#> 537        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 538        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 539        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 540        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 541        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 542        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 543        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 544        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 545        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 546        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 547        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 548        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 549    cohesion       plan    emotion       <NA>       <NA>       <NA>
#> 550  coregulate  synthesis  consensus coregulate       <NA>       <NA>
#> 551    cohesion  consensus       plan       <NA>       <NA>       <NA>
#> 552        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 553        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 554        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 555   consensus    emotion  consensus       <NA>       <NA>       <NA>
#> 556   consensus coregulate  consensus       plan       plan       <NA>
#> 557   synthesis       plan       plan       plan       <NA>       <NA>
#> 558        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 559        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 560     discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 561   consensus    emotion  consensus       <NA>       <NA>       <NA>
#> 562        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 563        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 564        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 565   consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 566        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 567        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 568        plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 569        plan       plan       plan       <NA>       <NA>       <NA>
#> 570        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 571     discuss  consensus    emotion    monitor       <NA>       <NA>
#> 572   consensus    discuss  synthesis      adapt       <NA>       <NA>
#> 573        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 574   consensus       plan       plan       <NA>       <NA>       <NA>
#> 575        plan    emotion   cohesion    emotion   cohesion    discuss
#> 576        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 577        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 578     discuss       plan       <NA>       <NA>       <NA>       <NA>
#> 579        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 580       adapt  consensus       plan coregulate   cohesion    discuss
#> 581        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 582        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 583        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 584        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 585        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 586        plan       plan    monitor   cohesion  consensus    discuss
#> 587        plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 588        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 589        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 590        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 591        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 592   consensus coregulate    discuss       <NA>       <NA>       <NA>
#> 593     emotion   cohesion coregulate  consensus coregulate       <NA>
#> 594        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 595        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 596        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 597        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 598        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 599   synthesis  consensus       plan  consensus coregulate       plan
#> 600        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 601        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 602     emotion   cohesion   cohesion  consensus    monitor    discuss
#> 603        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 604        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 605        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 606        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 607        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 608        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 609        plan    monitor       plan       plan       plan       plan
#> 610        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 611     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 612        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 613        plan       plan       plan       plan    monitor    discuss
#> 614        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 615        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 616        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 617     discuss  synthesis  consensus    discuss       <NA>       <NA>
#> 618        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 619   consensus    discuss    discuss       <NA>       <NA>       <NA>
#> 620        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 621     discuss    emotion       <NA>       <NA>       <NA>       <NA>
#> 622        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 623   consensus    emotion       plan       plan       <NA>       <NA>
#> 624        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 625        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 626        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 627        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 628        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 629        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 630        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 631        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 632   consensus       plan  consensus       plan    discuss  consensus
#> 633     emotion   cohesion  consensus       plan  consensus    discuss
#> 634       adapt   cohesion  consensus coregulate    emotion       plan
#> 635        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 636        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 637        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 638        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 639        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 640        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 641        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 642    cohesion  consensus    discuss       <NA>       <NA>       <NA>
#> 643        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 644   consensus    discuss       <NA>       <NA>       <NA>       <NA>
#> 645        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 646        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 647        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 648   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 649     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 650        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 651        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 652        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 653     discuss    emotion   cohesion  consensus coregulate      adapt
#> 654        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 655        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 656        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 657   consensus       plan    monitor       plan       plan  consensus
#> 658     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 659        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 660   consensus    discuss  consensus       <NA>       <NA>       <NA>
#> 661        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 662     discuss  consensus coregulate    discuss coregulate    emotion
#> 663        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 664  coregulate       plan       plan       <NA>       <NA>       <NA>
#> 665        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 666        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 667        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 668        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 669        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 670        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 671        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 672        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 673     emotion   cohesion  consensus       plan       plan       <NA>
#> 674        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 675   consensus       plan  consensus    discuss  synthesis  consensus
#> 676        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 677  coregulate    discuss    emotion   cohesion       plan       plan
#> 678     discuss    emotion       <NA>       <NA>       <NA>       <NA>
#> 679        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 680        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 681        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 682        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 683        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 684        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 685        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 686    cohesion    discuss  synthesis  consensus       <NA>       <NA>
#> 687     emotion  consensus       plan       <NA>       <NA>       <NA>
#> 688        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 689        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 690  coregulate       plan  consensus       <NA>       <NA>       <NA>
#> 691        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 692        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 693     monitor  consensus       <NA>       <NA>       <NA>       <NA>
#> 694        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 695        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 696        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 697        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 698        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 699        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 700        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 701        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 702        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 703        plan  consensus  consensus  consensus coregulate       <NA>
#> 704        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 705        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 706        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 707        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 708        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 709     monitor       plan  consensus       <NA>       <NA>       <NA>
#> 710     discuss    monitor    discuss  consensus    discuss   cohesion
#> 711        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 712   consensus       plan   cohesion coregulate       plan  consensus
#> 713        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 714        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 715        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 716        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 717     discuss  consensus       plan       plan       plan    emotion
#> 718        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 719        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 720        plan  consensus       plan  consensus coregulate    emotion
#> 721     monitor  consensus coregulate    emotion       <NA>       <NA>
#> 722        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 723        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 724        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 725        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 726        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 727        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 728   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 729        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 730        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 731    cohesion  consensus       <NA>       <NA>       <NA>       <NA>
#> 732     emotion   cohesion    monitor    discuss  consensus  consensus
#> 733        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 734        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 735     monitor    discuss    monitor    discuss       <NA>       <NA>
#> 736     monitor       plan  consensus    discuss  synthesis  consensus
#> 737        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 738        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 739        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 740        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 741        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 742        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 743        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 744    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 745        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 746     discuss  consensus       plan  consensus coregulate       plan
#> 747        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 748        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 749        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 750     discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 751    cohesion  consensus    discuss       <NA>       <NA>       <NA>
#> 752        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 753        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 754        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 755   consensus       plan    emotion  consensus    discuss  consensus
#> 756        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 757        plan coregulate    discuss    discuss  consensus    emotion
#> 758        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 759        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 760        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 761        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 762        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 763        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 764        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 765        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 766        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 767   consensus       plan  consensus coregulate    emotion  consensus
#> 768        plan    monitor       plan       plan  consensus coregulate
#> 769        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 770        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 771        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 772        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 773        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 774     monitor   cohesion  consensus       plan    discuss    discuss
#> 775        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 776        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 777        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 778        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 779        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 780        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 781        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 782     monitor       plan   cohesion       plan       plan  consensus
#> 783        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 784        plan    discuss  synthesis    discuss       <NA>       <NA>
#> 785        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 786        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 787        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 788        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 789        plan    discuss coregulate  consensus       plan       plan
#> 790        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 791        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 792   consensus    discuss  consensus       plan       plan       plan
#> 793        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 794        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 795        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 796        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 797        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 798        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 799   consensus    discuss  consensus       <NA>       <NA>       <NA>
#> 800        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 801        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 802        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 803        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 804   consensus       plan    emotion  consensus       <NA>       <NA>
#> 805     discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 806        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 807        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 808        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 809        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 810        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 811        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 812        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 813        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 814        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 815        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 816        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 817        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 818        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 819  coregulate       plan    discuss  consensus    discuss       <NA>
#> 820        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 821        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 822   consensus       plan       plan       plan       <NA>       <NA>
#> 823        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 824    cohesion       plan       plan  consensus       plan       plan
#> 825   consensus    discuss  consensus       plan       plan  consensus
#> 826        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 827        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 828        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 829        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 830    cohesion  consensus       <NA>       <NA>       <NA>       <NA>
#> 831        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 832        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 833     discuss  consensus coregulate       plan  consensus       plan
#> 834   consensus    monitor  consensus       <NA>       <NA>       <NA>
#> 835        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 836        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 837  coregulate    emotion       plan    monitor    emotion  synthesis
#> 838        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 839        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 840        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 841        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 842        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 843        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 844        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 845        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 846        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 847        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 848        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 849        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 850        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 851   consensus    emotion   cohesion  consensus  consensus    discuss
#> 852        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 853        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 854        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 855        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 856        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 857        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 858        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 859        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 860        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 861        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 862    cohesion  consensus    emotion  consensus coregulate    emotion
#> 863   consensus    emotion       <NA>       <NA>       <NA>       <NA>
#> 864        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 865     emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 866        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 867        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 868        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 869     emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 870   consensus    discuss  consensus       <NA>       <NA>       <NA>
#> 871        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 872        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 873        plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 874        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 875        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 876   consensus       plan   cohesion  consensus       <NA>       <NA>
#> 877        plan       plan  consensus       plan    emotion       <NA>
#> 878        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 879        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 880        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 881        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 882        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 883     emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 884        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 885  coregulate    discuss   cohesion  consensus    emotion       <NA>
#> 886        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 887        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 888        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 889        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 890        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 891        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 892        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 893        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 894        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 895        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 896        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 897        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 898        plan    emotion    discuss    discuss       <NA>       <NA>
#> 899        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 900     discuss  consensus       plan       plan       <NA>       <NA>
#> 901        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 902        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 903   consensus       plan    monitor  consensus    emotion   cohesion
#> 904        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 905        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 906        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 907        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 908        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 909   consensus    discuss    emotion       <NA>       <NA>       <NA>
#> 910    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 911        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 912        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 913        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 914        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 915        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 916        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 917        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 918        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 919     discuss  consensus       plan  consensus       <NA>       <NA>
#> 920        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 921        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 922        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 923        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 924        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 925        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 926        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 927        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 928        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 929        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 930        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 931        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 932        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 933     discuss    discuss  synthesis       <NA>       <NA>       <NA>
#> 934        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 935        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 936        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 937        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 938        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 939        plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 940        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 941        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 942        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 943        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 944        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 945        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 946        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 947        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 948  coregulate    discuss       <NA>       <NA>       <NA>       <NA>
#> 949   consensus       plan       plan       plan       <NA>       <NA>
#> 950        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 951     discuss    discuss  synthesis       <NA>       <NA>       <NA>
#> 952        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 953        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 954   consensus  consensus       plan       <NA>       <NA>       <NA>
#> 955        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 956     emotion  consensus coregulate    emotion       plan       plan
#> 957        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 958        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 959        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 960     discuss      adapt  consensus       plan       <NA>       <NA>
#> 961        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 962        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 963        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 964        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 965        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 966     discuss  consensus       plan       plan    emotion   cohesion
#> 967     emotion   cohesion  consensus  consensus       <NA>       <NA>
#> 968        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 969        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 970     discuss  consensus coregulate    discuss       plan       <NA>
#> 971   consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 972        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 973        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 974        plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 975  coregulate       plan  consensus       <NA>       <NA>       <NA>
#> 976        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 977        plan    emotion coregulate       plan  consensus       plan
#> 978        plan    monitor    discuss    discuss    emotion       <NA>
#> 979        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 980        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 981        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 982        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 983        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 984        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 985        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 986        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 987        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 988        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 989     discuss  consensus    discuss  synthesis       <NA>       <NA>
#> 990  coregulate      adapt  consensus       plan       <NA>       <NA>
#> 991        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 992        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 993        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 994        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 995        plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 996        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 997        plan       plan    monitor       <NA>       <NA>       <NA>
#> 998    cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 999        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1000       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1001       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1002       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1003       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1004  consensus coregulate    discuss       <NA>       <NA>       <NA>
#> 1005       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1006       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1007       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1008    monitor       plan       plan  consensus    discuss       <NA>
#> 1009       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1010    emotion   cohesion  consensus  consensus    discuss      adapt
#> 1011       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1012       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1013       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1014    discuss    discuss  consensus       plan       plan  consensus
#> 1015    emotion   cohesion coregulate    discuss      adapt   cohesion
#> 1016      adapt  consensus coregulate    discuss  consensus       <NA>
#> 1017       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1018       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1019       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1020       plan  consensus coregulate       plan       plan coregulate
#> 1021       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1022       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1023       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1024  consensus    monitor    discuss  synthesis  consensus coregulate
#> 1025       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1026       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1027  consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 1028       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1029       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1030  consensus coregulate   cohesion    emotion   cohesion  consensus
#> 1031       plan    emotion coregulate    monitor  consensus       <NA>
#> 1032 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1033    emotion    monitor    discuss coregulate       <NA>       <NA>
#> 1034       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1035       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1036       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1037       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1038       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1039       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1040       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1041       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1042       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1043       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1044       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1045       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1046  consensus    discuss    discuss    discuss    discuss      adapt
#> 1047       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1048       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1049       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1050       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1051       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1052       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1053       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1054       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1055       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1056       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1057       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1058       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1059    emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1060       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1061       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1062       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1063       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1064       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1065       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1066       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1067       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1068       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1069       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1070       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1071       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1072       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1073       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1074       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1075  consensus       plan   cohesion  consensus       <NA>       <NA>
#> 1076       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1077       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1078       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1079    discuss      adapt  consensus       plan       plan       <NA>
#> 1080       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1081       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1082       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1083       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1084       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1085       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1086    discuss coregulate       <NA>       <NA>       <NA>       <NA>
#> 1087       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1088       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1089       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1090  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1091       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1092    emotion    emotion    discuss coregulate    discuss  synthesis
#> 1093  consensus  consensus       plan    discuss  synthesis       <NA>
#> 1094       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1095       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1096       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1097       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1098       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1099       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1100       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1101       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1102       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1103 coregulate  consensus coregulate    discuss       <NA>       <NA>
#> 1104    discuss  consensus    emotion    discuss       <NA>       <NA>
#> 1105       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1106  consensus    discuss      adapt   cohesion    emotion    discuss
#> 1107       plan  consensus       plan  consensus       plan       plan
#> 1108       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1109       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1110       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1111  consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 1112       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1113       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1114       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1115       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1116    monitor   cohesion    emotion       plan       plan       <NA>
#> 1117      adapt   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1118       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1119       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1120       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1121       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1122       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1123    monitor    monitor    emotion  consensus    discuss  synthesis
#> 1124    monitor       plan    discuss  consensus    discuss  synthesis
#> 1125       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1126 coregulate    emotion coregulate coregulate    monitor    discuss
#> 1127       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1128  consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 1129    discuss    emotion  consensus  consensus       plan       <NA>
#> 1130       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1131    monitor    monitor       <NA>       <NA>       <NA>       <NA>
#> 1132       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1133       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1134       plan  consensus    monitor    discuss    discuss       <NA>
#> 1135       plan       plan    discuss coregulate    discuss       <NA>
#> 1136  synthesis coregulate    discuss  synthesis       <NA>       <NA>
#> 1137    monitor       plan       plan       plan       <NA>       <NA>
#> 1138       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1139  consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 1140       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1141      adapt  consensus       plan       plan coregulate    discuss
#> 1142    emotion   cohesion  consensus       plan       <NA>       <NA>
#> 1143       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1144       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1145       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1146       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1147  consensus coregulate       plan       plan       <NA>       <NA>
#> 1148    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1149  consensus       plan  consensus   cohesion       <NA>       <NA>
#> 1150 coregulate       plan       plan       <NA>       <NA>       <NA>
#> 1151       plan coregulate  consensus       plan       plan       plan
#> 1152       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1153    discuss    monitor       plan  consensus       plan  consensus
#> 1154    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1155       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1156       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1157       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1158 coregulate    discuss       <NA>       <NA>       <NA>       <NA>
#> 1159       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1160       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1161       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1162       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1163       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1164       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1165       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1166       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1167       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1168       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1169       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1170       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1171  consensus coregulate    emotion    emotion       <NA>       <NA>
#> 1172       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1173       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1174       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1175       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1176       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1177       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1178       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1179  synthesis       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1180       plan       plan       plan    emotion  consensus       plan
#> 1181       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1182       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1183       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1184    emotion   cohesion coregulate  consensus       plan       plan
#> 1185       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1186   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1187    emotion    emotion       plan       plan       plan       <NA>
#> 1188       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1189       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1190       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1191    monitor    discuss    monitor       <NA>       <NA>       <NA>
#> 1192       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1193       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1194       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1195       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1196       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1197       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1198       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1199       plan       plan  consensus    discuss       <NA>       <NA>
#> 1200       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1201       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1202       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1203       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1204  consensus coregulate    discuss       <NA>       <NA>       <NA>
#> 1205       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1206       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1207    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1208       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1209  consensus coregulate    emotion       plan  consensus    monitor
#> 1210       plan    monitor       plan       <NA>       <NA>       <NA>
#> 1211       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1212       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1213   cohesion coregulate  consensus       plan       <NA>       <NA>
#> 1214       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1215   cohesion  consensus coregulate  consensus    emotion   cohesion
#> 1216       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1217       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1218       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1219       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1220       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1221       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1222       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1223       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1224       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1225       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1226    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1227    monitor    discuss      adapt       <NA>       <NA>       <NA>
#> 1228       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1229 coregulate    monitor       <NA>       <NA>       <NA>       <NA>
#> 1230       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1231       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1232       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1233       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1234       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1235       plan  consensus       plan       <NA>       <NA>       <NA>
#> 1236       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1237       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1238 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1239       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1240  consensus       plan       plan    discuss    discuss       <NA>
#> 1241       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1242       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1243       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1244       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1245       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1246 coregulate  consensus       plan    emotion       <NA>       <NA>
#> 1247       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1248       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1249       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1250       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1251       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1252       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1253       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1254    discuss  consensus       plan    emotion       <NA>       <NA>
#> 1255       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1256       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1257       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1258       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1259       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1260       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1261       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1262       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1263       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1264       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1265       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1266       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1267  consensus       plan       plan    discuss coregulate    discuss
#> 1268  consensus  consensus       plan    discuss  synthesis       <NA>
#> 1269       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1270    discuss  synthesis      adapt  consensus coregulate       <NA>
#> 1271       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1272  synthesis  consensus coregulate coregulate    discuss  synthesis
#> 1273       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1274       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1275       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1276       plan       plan    discuss  consensus coregulate    emotion
#> 1277       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1278       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1279       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1280       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1281  consensus       plan  consensus       <NA>       <NA>       <NA>
#> 1282       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1283       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1284       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1285       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1286  consensus    discuss      adapt  consensus       <NA>       <NA>
#> 1287       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1288       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1289       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1290       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1291       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1292       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1293       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1294       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1295      adapt       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1296  consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 1297  consensus       plan       plan   cohesion coregulate    emotion
#> 1298       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1299       plan    discuss    discuss  synthesis       plan       plan
#> 1300 coregulate    discuss       <NA>       <NA>       <NA>       <NA>
#> 1301       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1302    discuss    discuss      adapt  consensus    monitor       plan
#> 1303       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1304       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1305   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1306    emotion   cohesion  consensus       <NA>       <NA>       <NA>
#> 1307       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1308       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1309       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1310       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1311       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1312       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1313       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1314       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1315       plan    emotion  consensus       plan       <NA>       <NA>
#> 1316       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1317       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1318       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1319       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1320       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1321       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1322       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1323       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1324       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1325       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1326       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1327       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1328       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1329  synthesis  consensus    monitor       <NA>       <NA>       <NA>
#> 1330       plan  consensus       plan       <NA>       <NA>       <NA>
#> 1331       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1332       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1333       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1334       plan    emotion  consensus    discuss  consensus    discuss
#> 1335       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1336       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1337       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1338       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1339       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1340       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1341       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1342       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1343       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1344       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1345  consensus       plan  consensus coregulate    emotion    discuss
#> 1346       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1347       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1348       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1349       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1350       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1351  consensus       plan       plan       plan       <NA>       <NA>
#> 1352    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1353       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1354    discuss  consensus       <NA>       <NA>       <NA>       <NA>
#> 1355       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1356  synthesis       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1357    discuss  synthesis coregulate       <NA>       <NA>       <NA>
#> 1358    emotion       plan       plan       <NA>       <NA>       <NA>
#> 1359    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1360       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1361       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1362    discuss    discuss      adapt  consensus coregulate       plan
#> 1363       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1364       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1365       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1366       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1367    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1368       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1369  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1370       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1371   cohesion coregulate       plan       <NA>       <NA>       <NA>
#> 1372       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1373       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1374       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1375       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1376       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1377       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1378       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1379       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1380       plan       plan       plan       <NA>       <NA>       <NA>
#> 1381       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1382   cohesion  consensus       plan    monitor       plan    monitor
#> 1383       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1384    emotion  consensus    monitor    emotion  consensus coregulate
#> 1385       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1386       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1387       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1388       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1389  synthesis  consensus       plan       plan       <NA>       <NA>
#> 1390       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1391       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1392       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1393       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1394       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1395       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1396       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1397       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1398       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1399       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1400   cohesion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1401   cohesion  consensus    discuss    discuss    emotion  consensus
#> 1402       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1403    discuss      adapt   cohesion  consensus       plan    discuss
#> 1404       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1405       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1406       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1407    emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 1408    emotion       plan    discuss  synthesis      adapt  consensus
#> 1409    discuss  synthesis      adapt  consensus    discuss  synthesis
#> 1410       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1411       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1412       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1413       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1414       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1415       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1416       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1417       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1418       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1419       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1420       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1421       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1422       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1423       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1424       plan   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1425  consensus  consensus    emotion   cohesion       <NA>       <NA>
#> 1426       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1427       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1428       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1429       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1430       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1431       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1432 coregulate    emotion coregulate       plan       plan    emotion
#> 1433    emotion   cohesion  consensus    discuss  consensus coregulate
#> 1434       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1435       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1436       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1437       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1438       plan    emotion    monitor    discuss   cohesion    discuss
#> 1439       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1440       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1441       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1442       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1443       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1444       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1445       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1446    emotion   cohesion  consensus coregulate       plan       plan
#> 1447       plan    emotion  consensus    discuss       <NA>       <NA>
#> 1448 coregulate       plan  consensus    discuss      adapt   cohesion
#> 1449       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1450       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1451       plan       plan       plan       plan       plan       <NA>
#> 1452       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1453  consensus    discuss    discuss    emotion   cohesion       <NA>
#> 1454    emotion    emotion   cohesion       <NA>       <NA>       <NA>
#> 1455       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1456       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1457       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1458       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1459       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1460       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1461       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1462      adapt       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1463       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1464       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1465    monitor       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1466       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1467       plan    emotion       plan       plan       <NA>       <NA>
#> 1468    emotion   cohesion  consensus coregulate       plan       <NA>
#> 1469       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1470       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1471       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1472       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1473 coregulate       plan       plan       <NA>       <NA>       <NA>
#> 1474  consensus       plan       plan coregulate    emotion  consensus
#> 1475       plan       plan    emotion   cohesion  consensus       <NA>
#> 1476       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1477 coregulate coregulate       plan       plan       <NA>       <NA>
#> 1478       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1479       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1480       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1481       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1482       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1483       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1484       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1485  consensus coregulate    discuss       <NA>       <NA>       <NA>
#> 1486       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1487       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1488  synthesis  consensus       <NA>       <NA>       <NA>       <NA>
#> 1489  consensus    monitor       plan       <NA>       <NA>       <NA>
#> 1490       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1491       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1492  consensus       plan    emotion  consensus       plan       plan
#> 1493       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1494       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1495       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1496       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1497  consensus    monitor    discuss    discuss    discuss    discuss
#> 1498   cohesion coregulate  synthesis       <NA>       <NA>       <NA>
#> 1499       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1500    emotion   cohesion  consensus    discuss coregulate    discuss
#> 1501       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1502       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1503  consensus    monitor  consensus coregulate       <NA>       <NA>
#> 1504 coregulate    discuss  consensus       plan       plan    emotion
#> 1505       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1506       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1507       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1508 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1509       plan  consensus  consensus coregulate       <NA>       <NA>
#> 1510       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1511       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1512       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1513       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1514       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1515 coregulate    discuss    discuss      adapt       <NA>       <NA>
#> 1516       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1517       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1518       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1519  synthesis      adapt    discuss  synthesis       <NA>       <NA>
#> 1520       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1521       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1522       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1523 coregulate       plan       <NA>       <NA>       <NA>       <NA>
#> 1524       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1525       plan       plan       plan       plan       plan    monitor
#> 1526    emotion       plan       <NA>       <NA>       <NA>       <NA>
#> 1527       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1528  consensus       plan    emotion coregulate    emotion   cohesion
#> 1529       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1530       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1531   cohesion  consensus    discuss    discuss    emotion       <NA>
#> 1532  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1533       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1534       plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 1535    emotion  consensus       plan       <NA>       <NA>       <NA>
#> 1536  consensus       plan       plan       <NA>       <NA>       <NA>
#> 1537       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1538       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1539  consensus coregulate       plan       plan    monitor    discuss
#> 1540       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1541       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1542  consensus       plan  consensus       plan       <NA>       <NA>
#> 1543       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1544       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1545       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1546       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1547       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1548       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1549       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1550       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1551       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1552       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1553       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1554       plan       plan       <NA>       <NA>       <NA>       <NA>
#> 1555       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1556       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1557       plan    discuss       <NA>       <NA>       <NA>       <NA>
#> 1558       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1559       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1560       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1561    emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1562       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1563       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1564       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1565       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1566       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1567       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1568 coregulate    emotion   cohesion       <NA>       <NA>       <NA>
#> 1569    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1570       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1571       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1572       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1573       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1574       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1575    emotion   cohesion coregulate  consensus    discuss  consensus
#> 1576       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1577    emotion   cohesion  consensus       <NA>       <NA>       <NA>
#> 1578  consensus       plan    discuss    discuss       <NA>       <NA>
#> 1579       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1580       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1581  consensus       plan   cohesion coregulate    discuss       <NA>
#> 1582       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1583       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1584    monitor    emotion   cohesion    emotion       <NA>       <NA>
#> 1585       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1586       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1587       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1588    emotion    monitor  synthesis      adapt  consensus       <NA>
#> 1589       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1590  consensus       plan       plan       plan       plan  consensus
#> 1591  consensus coregulate       plan       <NA>       <NA>       <NA>
#> 1592       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1593    discuss    discuss    discuss  consensus       <NA>       <NA>
#> 1594       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1595       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1596       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1597       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1598       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1599    emotion       plan       plan       plan    discuss       <NA>
#> 1600    monitor   cohesion  consensus    discuss       <NA>       <NA>
#> 1601       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1602       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1603       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1604       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1605       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1606       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1607  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1608       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1609       plan  consensus    emotion   cohesion    discuss    discuss
#> 1610       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1611  consensus    emotion   cohesion       plan  consensus    emotion
#> 1612       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1613       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1614       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1615       plan       plan    emotion   cohesion    emotion       <NA>
#> 1616       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1617       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1618       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1619       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1620       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1621       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1622  consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 1623       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1624       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1625       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1626       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1627       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1628       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1629       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1630       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1631       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1632       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1633       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1634  consensus  consensus    emotion  consensus       plan       plan
#> 1635   cohesion    monitor    discuss      adapt       <NA>       <NA>
#> 1636       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1637  consensus       plan  consensus       <NA>       <NA>       <NA>
#> 1638       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1639       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1640  consensus       plan       plan       plan       plan       plan
#> 1641 coregulate       plan    emotion       plan       plan  consensus
#> 1642       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1643       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1644       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1645    discuss    emotion   cohesion    emotion   cohesion       <NA>
#> 1646       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1647       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1648    discuss    discuss       <NA>       <NA>       <NA>       <NA>
#> 1649       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1650       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1651       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1652       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1653    emotion  consensus  consensus    discuss       <NA>       <NA>
#> 1654    emotion    discuss       <NA>       <NA>       <NA>       <NA>
#> 1655       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1656       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1657       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1658       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1659    discuss    discuss    discuss  consensus  consensus    discuss
#> 1660    discuss    discuss    emotion coregulate  consensus       plan
#> 1661       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1662       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1663    emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1664       plan    discuss    discuss   cohesion  consensus       plan
#> 1665       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1666       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1667       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1668       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1669       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1670       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1671    discuss    emotion       <NA>       <NA>       <NA>       <NA>
#> 1672       plan       plan       plan    monitor       plan       plan
#> 1673       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1674       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1675       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1676       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1677       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1678       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1679       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1680       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1681    monitor    discuss       <NA>       <NA>       <NA>       <NA>
#> 1682       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1683       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1684       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1685       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1686       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1687    emotion coregulate       plan       plan       <NA>       <NA>
#> 1688       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1689       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1690       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1691  consensus       plan       <NA>       <NA>       <NA>       <NA>
#> 1692       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1693       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1694  consensus    discuss      adapt   cohesion coregulate   cohesion
#> 1695       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1696 coregulate       plan       <NA>       <NA>       <NA>       <NA>
#> 1697       plan       plan  consensus       plan    emotion       plan
#> 1698       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1699  consensus   cohesion    monitor       <NA>       <NA>       <NA>
#> 1700    emotion  consensus       plan    discuss  synthesis coregulate
#> 1701       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1702       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1703       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1704       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1705       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1706       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1707       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1708       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1709       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1710       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1711       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1712       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1713       plan       plan    monitor  consensus       plan       plan
#> 1714       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1715 coregulate  consensus    discuss    monitor       plan       plan
#> 1716       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1717       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1718       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1719       plan    monitor  synthesis  consensus    discuss  synthesis
#> 1720       plan       plan  consensus       plan  consensus       <NA>
#> 1721       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1722   cohesion  consensus    emotion  consensus       plan    monitor
#> 1723       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1724       plan       plan  consensus       plan    emotion  consensus
#> 1725       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1726       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1727    monitor   cohesion  consensus coregulate       <NA>       <NA>
#> 1728       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1729       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1730       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1731       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1732       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1733       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1734       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1735  consensus    discuss coregulate    discuss    emotion    discuss
#> 1736       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1737       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1738       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1739       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1740  consensus       plan       plan    emotion    discuss    discuss
#> 1741       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1742       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1743       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1744       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1745       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1746       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1747       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1748    monitor   cohesion  consensus coregulate       <NA>       <NA>
#> 1749       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1750   cohesion coregulate  consensus    monitor   cohesion  consensus
#> 1751       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1752       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1753       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1754       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1755       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1756       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1757       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1758       plan       plan       plan    monitor       <NA>       <NA>
#> 1759       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1760       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1761       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1762       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1763    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1764       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1765       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1766       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1767  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1768       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1769       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1770 coregulate       plan       plan    monitor       <NA>       <NA>
#> 1771  synthesis      adapt       <NA>       <NA>       <NA>       <NA>
#> 1772       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1773       plan  consensus coregulate       plan  consensus       plan
#> 1774  synthesis    monitor   cohesion  consensus       plan       <NA>
#> 1775    emotion  consensus       <NA>       <NA>       <NA>       <NA>
#> 1776  synthesis  consensus       plan       plan       plan       <NA>
#> 1777       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1778       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1779       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1780       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1781       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1782       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1783       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1784       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1785       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1786       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1787       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1788       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1789       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1790       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1791    emotion    monitor    discuss  consensus       <NA>       <NA>
#> 1792       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1793       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1794       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1795       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1796       plan  consensus coregulate coregulate    discuss       <NA>
#> 1797       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1798  consensus    discuss    emotion  consensus       <NA>       <NA>
#> 1799       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1800 coregulate    discuss       <NA>       <NA>       <NA>       <NA>
#> 1801       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1802       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1803       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1804    monitor      adapt  consensus       <NA>       <NA>       <NA>
#> 1805 coregulate       plan       plan       plan  consensus       plan
#> 1806 coregulate       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1807       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1808       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1809       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1810       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1811 coregulate    discuss  consensus coregulate       <NA>       <NA>
#> 1812       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1813  consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 1814 coregulate       plan       plan  consensus       plan       <NA>
#> 1815       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1816  consensus       plan       plan       <NA>       <NA>       <NA>
#> 1817       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1818       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1819       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1820       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1821       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1822       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1823       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1824 coregulate    monitor       <NA>       <NA>       <NA>       <NA>
#> 1825       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1826       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1827  consensus coregulate coregulate    emotion   cohesion       <NA>
#> 1828    monitor    emotion  consensus       plan  consensus       <NA>
#> 1829       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1830      adapt    emotion    emotion  consensus       plan       <NA>
#> 1831       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1832       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1833       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1834       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1835       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1836       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1837 coregulate       plan    discuss    emotion   cohesion       <NA>
#> 1838       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1839       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1840       plan       plan       plan       plan  consensus       plan
#> 1841       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1842    monitor  consensus       plan       <NA>       <NA>       <NA>
#> 1843       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1844       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1845       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1846    discuss   cohesion  consensus    emotion       <NA>       <NA>
#> 1847       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1848       plan    monitor  synthesis  consensus       plan  consensus
#> 1849       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1850       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1851       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1852    emotion  consensus       plan       <NA>       <NA>       <NA>
#> 1853       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1854       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1855    discuss      adapt  consensus       plan       plan  consensus
#> 1856       plan    monitor       <NA>       <NA>       <NA>       <NA>
#> 1857       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1858       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1859       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1860       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1861       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1862       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1863       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1864       plan  consensus       <NA>       <NA>       <NA>       <NA>
#> 1865       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1866      adapt    emotion coregulate       plan       <NA>       <NA>
#> 1867    discuss    emotion       plan       plan       plan       <NA>
#> 1868       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1869       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1870    emotion       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1871       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1872       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1873       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1874       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1875      adapt  consensus       <NA>       <NA>       <NA>       <NA>
#> 1876       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1877       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1878    emotion  consensus       plan  consensus coregulate    emotion
#> 1879       plan       plan       plan   cohesion coregulate       plan
#> 1880       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1881       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1882       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1883       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1884       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1885       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1886       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1887       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1888       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1889       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1890  consensus    discuss      adapt  consensus    emotion   cohesion
#> 1891       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1892       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1893       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1894       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1895       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1896       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1897       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1898   cohesion  consensus  consensus       plan  consensus       plan
#> 1899       plan       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1900    emotion   cohesion    monitor    discuss    emotion   cohesion
#> 1901       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1902       plan    emotion       <NA>       <NA>       <NA>       <NA>
#> 1903       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1904       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1905       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1906       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1907       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1908       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1909       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1910       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1911       plan    emotion  consensus       plan    discuss  synthesis
#> 1912       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1913       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1914      adapt  consensus coregulate    emotion       <NA>       <NA>
#> 1915       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1916       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1917       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1918       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1919       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1920       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1921       plan  consensus    monitor       plan       <NA>       <NA>
#> 1922  consensus       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1923 coregulate    monitor    emotion  consensus       <NA>       <NA>
#> 1924       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1925  consensus       plan       plan    emotion  consensus       <NA>
#> 1926       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1927       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1928       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1929       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1930    emotion   cohesion    emotion coregulate    monitor    discuss
#> 1931       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1932       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1933       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1934       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1935    emotion   cohesion       <NA>       <NA>       <NA>       <NA>
#> 1936 coregulate    emotion   cohesion    emotion    discuss    monitor
#> 1937       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1938    discuss       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1939       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1940       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1941       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1942       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1943       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1944  synthesis    emotion  consensus    monitor    discuss       <NA>
#> 1945       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1946       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1947  consensus       plan    discuss  consensus       plan       <NA>
#> 1948       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1949       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1950       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1951       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1952       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1953       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1954       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1955       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1956 coregulate    monitor    discuss      adapt    discuss  consensus
#> 1957   cohesion coregulate       plan       plan    monitor coregulate
#> 1958       plan    discuss  consensus  consensus    monitor    emotion
#> 1959       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1960       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1961       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1962       plan       plan  consensus    discuss    discuss      adapt
#> 1963       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1964       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1965    discuss  consensus    monitor       <NA>       <NA>       <NA>
#> 1966       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1967       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1968       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1969       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1970       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1971       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1972       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1973       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1974       plan    emotion    emotion  consensus       plan       plan
#> 1975    emotion   cohesion    monitor    discuss       <NA>       <NA>
#> 1976       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1977       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1978       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1979  consensus coregulate       <NA>       <NA>       <NA>       <NA>
#> 1980       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1981       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1982       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1983       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1984       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1985       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1986       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1987       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1988    discuss  consensus  consensus    discuss coregulate       plan
#> 1989       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1990       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1991       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1992       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1993       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1994    discuss      adapt       <NA>       <NA>       <NA>       <NA>
#> 1995       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1996 coregulate    discuss  consensus    discuss       <NA>       <NA>
#> 1997       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1998       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 1999       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 2000       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#>             T26
#> 1          <NA>
#> 2          <NA>
#> 3       discuss
#> 4          <NA>
#> 5          <NA>
#> 6          <NA>
#> 7          <NA>
#> 8          <NA>
#> 9          <NA>
#> 10         <NA>
#> 11         <NA>
#> 12         <NA>
#> 13         <NA>
#> 14         <NA>
#> 15         <NA>
#> 16         <NA>
#> 17         <NA>
#> 18         <NA>
#> 19         <NA>
#> 20         <NA>
#> 21         <NA>
#> 22         <NA>
#> 23         <NA>
#> 24         <NA>
#> 25         <NA>
#> 26         <NA>
#> 27         <NA>
#> 28         <NA>
#> 29         <NA>
#> 30         <NA>
#> 31         <NA>
#> 32         <NA>
#> 33         <NA>
#> 34         <NA>
#> 35         <NA>
#> 36         <NA>
#> 37         <NA>
#> 38         <NA>
#> 39         <NA>
#> 40         <NA>
#> 41      emotion
#> 42         <NA>
#> 43         <NA>
#> 44         <NA>
#> 45         <NA>
#> 46         <NA>
#> 47         <NA>
#> 48         <NA>
#> 49         <NA>
#> 50         <NA>
#> 51         <NA>
#> 52         <NA>
#> 53         <NA>
#> 54         <NA>
#> 55         <NA>
#> 56         <NA>
#> 57         <NA>
#> 58         <NA>
#> 59         <NA>
#> 60         <NA>
#> 61         <NA>
#> 62         <NA>
#> 63         <NA>
#> 64         <NA>
#> 65         <NA>
#> 66         <NA>
#> 67      emotion
#> 68         <NA>
#> 69         <NA>
#> 70         <NA>
#> 71         <NA>
#> 72         <NA>
#> 73         <NA>
#> 74         <NA>
#> 75         <NA>
#> 76         <NA>
#> 77         <NA>
#> 78         <NA>
#> 79         <NA>
#> 80         <NA>
#> 81         <NA>
#> 82         <NA>
#> 83         plan
#> 84         <NA>
#> 85         <NA>
#> 86         <NA>
#> 87         <NA>
#> 88         <NA>
#> 89         <NA>
#> 90         <NA>
#> 91         <NA>
#> 92         <NA>
#> 93         <NA>
#> 94         <NA>
#> 95         <NA>
#> 96     cohesion
#> 97         <NA>
#> 98         <NA>
#> 99         <NA>
#> 100        <NA>
#> 101        <NA>
#> 102        <NA>
#> 103        <NA>
#> 104        <NA>
#> 105        <NA>
#> 106        <NA>
#> 107        <NA>
#> 108        <NA>
#> 109        <NA>
#> 110        <NA>
#> 111        <NA>
#> 112        <NA>
#> 113        <NA>
#> 114        <NA>
#> 115        <NA>
#> 116        <NA>
#> 117        <NA>
#> 118        <NA>
#> 119        <NA>
#> 120        <NA>
#> 121        <NA>
#> 122        <NA>
#> 123        <NA>
#> 124        <NA>
#> 125        <NA>
#> 126        <NA>
#> 127        <NA>
#> 128        <NA>
#> 129        <NA>
#> 130        <NA>
#> 131        <NA>
#> 132        <NA>
#> 133        <NA>
#> 134  coregulate
#> 135        <NA>
#> 136        <NA>
#> 137        <NA>
#> 138        <NA>
#> 139        <NA>
#> 140        <NA>
#> 141        <NA>
#> 142        <NA>
#> 143        <NA>
#> 144        <NA>
#> 145        <NA>
#> 146        <NA>
#> 147        <NA>
#> 148        <NA>
#> 149        <NA>
#> 150        <NA>
#> 151        <NA>
#> 152        <NA>
#> 153        <NA>
#> 154        <NA>
#> 155        <NA>
#> 156        <NA>
#> 157        <NA>
#> 158        <NA>
#> 159        <NA>
#> 160        <NA>
#> 161        <NA>
#> 162        <NA>
#> 163        <NA>
#> 164        <NA>
#> 165        <NA>
#> 166        <NA>
#> 167        <NA>
#> 168        <NA>
#> 169        <NA>
#> 170        <NA>
#> 171        <NA>
#> 172        <NA>
#> 173        <NA>
#> 174        <NA>
#> 175        <NA>
#> 176        <NA>
#> 177        <NA>
#> 178        <NA>
#> 179        <NA>
#> 180        <NA>
#> 181        <NA>
#> 182        <NA>
#> 183        <NA>
#> 184        <NA>
#> 185        <NA>
#> 186        <NA>
#> 187        <NA>
#> 188     monitor
#> 189        <NA>
#> 190        <NA>
#> 191        <NA>
#> 192        <NA>
#> 193        <NA>
#> 194        <NA>
#> 195        <NA>
#> 196        <NA>
#> 197        <NA>
#> 198        <NA>
#> 199        <NA>
#> 200        <NA>
#> 201        <NA>
#> 202        <NA>
#> 203        <NA>
#> 204        <NA>
#> 205        <NA>
#> 206        <NA>
#> 207        <NA>
#> 208        <NA>
#> 209        <NA>
#> 210        <NA>
#> 211        <NA>
#> 212        <NA>
#> 213        <NA>
#> 214        <NA>
#> 215        <NA>
#> 216        <NA>
#> 217        <NA>
#> 218        <NA>
#> 219        <NA>
#> 220        <NA>
#> 221        <NA>
#> 222        <NA>
#> 223        <NA>
#> 224        <NA>
#> 225        <NA>
#> 226        <NA>
#> 227        <NA>
#> 228        <NA>
#> 229        <NA>
#> 230        <NA>
#> 231        <NA>
#> 232        <NA>
#> 233        <NA>
#> 234        <NA>
#> 235        <NA>
#> 236        <NA>
#> 237        <NA>
#> 238        <NA>
#> 239        <NA>
#> 240        <NA>
#> 241        <NA>
#> 242        <NA>
#> 243        <NA>
#> 244        <NA>
#> 245        <NA>
#> 246        <NA>
#> 247        <NA>
#> 248        <NA>
#> 249        <NA>
#> 250        <NA>
#> 251        <NA>
#> 252        <NA>
#> 253        <NA>
#> 254        <NA>
#> 255        <NA>
#> 256        <NA>
#> 257        <NA>
#> 258        <NA>
#> 259        <NA>
#> 260        <NA>
#> 261        <NA>
#> 262        <NA>
#> 263        <NA>
#> 264        <NA>
#> 265        <NA>
#> 266        <NA>
#> 267        <NA>
#> 268        <NA>
#> 269        <NA>
#> 270        <NA>
#> 271        <NA>
#> 272        <NA>
#> 273        <NA>
#> 274        <NA>
#> 275     discuss
#> 276        <NA>
#> 277        <NA>
#> 278        <NA>
#> 279        <NA>
#> 280        <NA>
#> 281        <NA>
#> 282        <NA>
#> 283        <NA>
#> 284        <NA>
#> 285        <NA>
#> 286        <NA>
#> 287        <NA>
#> 288        <NA>
#> 289        <NA>
#> 290        <NA>
#> 291        <NA>
#> 292        <NA>
#> 293        <NA>
#> 294        <NA>
#> 295        <NA>
#> 296        <NA>
#> 297        <NA>
#> 298        <NA>
#> 299        <NA>
#> 300        <NA>
#> 301        <NA>
#> 302        <NA>
#> 303        <NA>
#> 304        <NA>
#> 305        <NA>
#> 306        <NA>
#> 307        <NA>
#> 308        <NA>
#> 309        <NA>
#> 310        <NA>
#> 311        <NA>
#> 312        <NA>
#> 313        <NA>
#> 314        <NA>
#> 315        <NA>
#> 316        <NA>
#> 317        <NA>
#> 318        <NA>
#> 319        <NA>
#> 320        <NA>
#> 321        <NA>
#> 322        <NA>
#> 323        <NA>
#> 324        <NA>
#> 325        <NA>
#> 326        <NA>
#> 327        <NA>
#> 328        <NA>
#> 329        <NA>
#> 330        <NA>
#> 331        <NA>
#> 332        <NA>
#> 333        <NA>
#> 334        <NA>
#> 335        <NA>
#> 336        <NA>
#> 337        <NA>
#> 338        <NA>
#> 339        <NA>
#> 340        <NA>
#> 341        <NA>
#> 342        <NA>
#> 343        <NA>
#> 344        <NA>
#> 345        <NA>
#> 346        <NA>
#> 347        <NA>
#> 348        <NA>
#> 349        <NA>
#> 350        <NA>
#> 351        <NA>
#> 352        <NA>
#> 353        <NA>
#> 354        <NA>
#> 355        <NA>
#> 356        <NA>
#> 357        <NA>
#> 358        <NA>
#> 359        <NA>
#> 360        <NA>
#> 361        <NA>
#> 362        <NA>
#> 363        <NA>
#> 364        <NA>
#> 365        <NA>
#> 366        <NA>
#> 367        <NA>
#> 368        <NA>
#> 369        <NA>
#> 370        <NA>
#> 371        <NA>
#> 372        <NA>
#> 373        <NA>
#> 374        <NA>
#> 375        <NA>
#> 376        <NA>
#> 377        <NA>
#> 378        <NA>
#> 379        <NA>
#> 380        <NA>
#> 381        <NA>
#> 382        <NA>
#> 383        <NA>
#> 384        <NA>
#> 385        <NA>
#> 386        <NA>
#> 387        <NA>
#> 388        <NA>
#> 389        <NA>
#> 390        <NA>
#> 391        <NA>
#> 392        <NA>
#> 393        <NA>
#> 394        <NA>
#> 395        <NA>
#> 396        <NA>
#> 397        <NA>
#> 398        <NA>
#> 399        <NA>
#> 400        <NA>
#> 401        <NA>
#> 402        <NA>
#> 403        <NA>
#> 404        <NA>
#> 405        <NA>
#> 406        <NA>
#> 407        <NA>
#> 408        <NA>
#> 409        <NA>
#> 410        <NA>
#> 411        <NA>
#> 412        <NA>
#> 413        <NA>
#> 414        <NA>
#> 415        <NA>
#> 416     emotion
#> 417        <NA>
#> 418        <NA>
#> 419        <NA>
#> 420        <NA>
#> 421        <NA>
#> 422        <NA>
#> 423        <NA>
#> 424        <NA>
#> 425        <NA>
#> 426        <NA>
#> 427        <NA>
#> 428        <NA>
#> 429        <NA>
#> 430        <NA>
#> 431        <NA>
#> 432        <NA>
#> 433        <NA>
#> 434        <NA>
#> 435        <NA>
#> 436        <NA>
#> 437        <NA>
#> 438        <NA>
#> 439        <NA>
#> 440        <NA>
#> 441        <NA>
#> 442        <NA>
#> 443        <NA>
#> 444     discuss
#> 445        <NA>
#> 446        <NA>
#> 447        <NA>
#> 448        <NA>
#> 449        <NA>
#> 450        <NA>
#> 451        <NA>
#> 452        <NA>
#> 453        <NA>
#> 454        <NA>
#> 455        <NA>
#> 456        <NA>
#> 457        <NA>
#> 458        <NA>
#> 459        <NA>
#> 460        <NA>
#> 461     emotion
#> 462        <NA>
#> 463        <NA>
#> 464        <NA>
#> 465        <NA>
#> 466        <NA>
#> 467        <NA>
#> 468        <NA>
#> 469        <NA>
#> 470        <NA>
#> 471   consensus
#> 472   synthesis
#> 473        <NA>
#> 474   consensus
#> 475        <NA>
#> 476        <NA>
#> 477     emotion
#> 478        <NA>
#> 479        <NA>
#> 480        <NA>
#> 481        <NA>
#> 482        <NA>
#> 483        <NA>
#> 484        <NA>
#> 485        <NA>
#> 486        <NA>
#> 487   consensus
#> 488        <NA>
#> 489        <NA>
#> 490        <NA>
#> 491        <NA>
#> 492        <NA>
#> 493        <NA>
#> 494        <NA>
#> 495        <NA>
#> 496        <NA>
#> 497        <NA>
#> 498        <NA>
#> 499        <NA>
#> 500        <NA>
#> 501        <NA>
#> 502        <NA>
#> 503        <NA>
#> 504        <NA>
#> 505        <NA>
#> 506        <NA>
#> 507        <NA>
#> 508        <NA>
#> 509        <NA>
#> 510        <NA>
#> 511    cohesion
#> 512        <NA>
#> 513        <NA>
#> 514        <NA>
#> 515        <NA>
#> 516        <NA>
#> 517        <NA>
#> 518        <NA>
#> 519        <NA>
#> 520        <NA>
#> 521        <NA>
#> 522        <NA>
#> 523        <NA>
#> 524        <NA>
#> 525        <NA>
#> 526        <NA>
#> 527        <NA>
#> 528        <NA>
#> 529        <NA>
#> 530        <NA>
#> 531        <NA>
#> 532        <NA>
#> 533        <NA>
#> 534        <NA>
#> 535        <NA>
#> 536        <NA>
#> 537        <NA>
#> 538        <NA>
#> 539        <NA>
#> 540        <NA>
#> 541        <NA>
#> 542        <NA>
#> 543        <NA>
#> 544        <NA>
#> 545        <NA>
#> 546        <NA>
#> 547        <NA>
#> 548        <NA>
#> 549        <NA>
#> 550        <NA>
#> 551        <NA>
#> 552        <NA>
#> 553        <NA>
#> 554        <NA>
#> 555        <NA>
#> 556        <NA>
#> 557        <NA>
#> 558        <NA>
#> 559        <NA>
#> 560        <NA>
#> 561        <NA>
#> 562        <NA>
#> 563        <NA>
#> 564        <NA>
#> 565        <NA>
#> 566        <NA>
#> 567        <NA>
#> 568        <NA>
#> 569        <NA>
#> 570        <NA>
#> 571        <NA>
#> 572        <NA>
#> 573        <NA>
#> 574        <NA>
#> 575        <NA>
#> 576        <NA>
#> 577        <NA>
#> 578        <NA>
#> 579        <NA>
#> 580   consensus
#> 581        <NA>
#> 582        <NA>
#> 583        <NA>
#> 584        <NA>
#> 585        <NA>
#> 586   consensus
#> 587        <NA>
#> 588        <NA>
#> 589        <NA>
#> 590        <NA>
#> 591        <NA>
#> 592        <NA>
#> 593        <NA>
#> 594        <NA>
#> 595        <NA>
#> 596        <NA>
#> 597        <NA>
#> 598        <NA>
#> 599        <NA>
#> 600        <NA>
#> 601        <NA>
#> 602        <NA>
#> 603        <NA>
#> 604        <NA>
#> 605        <NA>
#> 606        <NA>
#> 607        <NA>
#> 608        <NA>
#> 609        <NA>
#> 610        <NA>
#> 611        <NA>
#> 612        <NA>
#> 613     discuss
#> 614        <NA>
#> 615        <NA>
#> 616        <NA>
#> 617        <NA>
#> 618        <NA>
#> 619        <NA>
#> 620        <NA>
#> 621        <NA>
#> 622        <NA>
#> 623        <NA>
#> 624        <NA>
#> 625        <NA>
#> 626        <NA>
#> 627        <NA>
#> 628        <NA>
#> 629        <NA>
#> 630        <NA>
#> 631        <NA>
#> 632        plan
#> 633   consensus
#> 634        <NA>
#> 635        <NA>
#> 636        <NA>
#> 637        <NA>
#> 638        <NA>
#> 639        <NA>
#> 640        <NA>
#> 641        <NA>
#> 642        <NA>
#> 643        <NA>
#> 644        <NA>
#> 645        <NA>
#> 646        <NA>
#> 647        <NA>
#> 648        <NA>
#> 649        <NA>
#> 650        <NA>
#> 651        <NA>
#> 652        <NA>
#> 653    cohesion
#> 654        <NA>
#> 655        <NA>
#> 656        <NA>
#> 657        <NA>
#> 658        <NA>
#> 659        <NA>
#> 660        <NA>
#> 661        <NA>
#> 662   consensus
#> 663        <NA>
#> 664        <NA>
#> 665        <NA>
#> 666        <NA>
#> 667        <NA>
#> 668        <NA>
#> 669        <NA>
#> 670        <NA>
#> 671        <NA>
#> 672        <NA>
#> 673        <NA>
#> 674        <NA>
#> 675        <NA>
#> 676        <NA>
#> 677        <NA>
#> 678        <NA>
#> 679        <NA>
#> 680        <NA>
#> 681        <NA>
#> 682        <NA>
#> 683        <NA>
#> 684        <NA>
#> 685        <NA>
#> 686        <NA>
#> 687        <NA>
#> 688        <NA>
#> 689        <NA>
#> 690        <NA>
#> 691        <NA>
#> 692        <NA>
#> 693        <NA>
#> 694        <NA>
#> 695        <NA>
#> 696        <NA>
#> 697        <NA>
#> 698        <NA>
#> 699        <NA>
#> 700        <NA>
#> 701        <NA>
#> 702        <NA>
#> 703        <NA>
#> 704        <NA>
#> 705        <NA>
#> 706        <NA>
#> 707        <NA>
#> 708        <NA>
#> 709        <NA>
#> 710        <NA>
#> 711        <NA>
#> 712        plan
#> 713        <NA>
#> 714        <NA>
#> 715        <NA>
#> 716        <NA>
#> 717   consensus
#> 718        <NA>
#> 719        <NA>
#> 720        <NA>
#> 721        <NA>
#> 722        <NA>
#> 723        <NA>
#> 724        <NA>
#> 725        <NA>
#> 726        <NA>
#> 727        <NA>
#> 728        <NA>
#> 729        <NA>
#> 730        <NA>
#> 731        <NA>
#> 732        <NA>
#> 733        <NA>
#> 734        <NA>
#> 735        <NA>
#> 736        <NA>
#> 737        <NA>
#> 738        <NA>
#> 739        <NA>
#> 740        <NA>
#> 741        <NA>
#> 742        <NA>
#> 743        <NA>
#> 744        <NA>
#> 745        <NA>
#> 746   consensus
#> 747        <NA>
#> 748        <NA>
#> 749        <NA>
#> 750        <NA>
#> 751        <NA>
#> 752        <NA>
#> 753        <NA>
#> 754        <NA>
#> 755        <NA>
#> 756        <NA>
#> 757   consensus
#> 758        <NA>
#> 759        <NA>
#> 760        <NA>
#> 761        <NA>
#> 762        <NA>
#> 763        <NA>
#> 764        <NA>
#> 765        <NA>
#> 766        <NA>
#> 767        <NA>
#> 768     emotion
#> 769        <NA>
#> 770        <NA>
#> 771        <NA>
#> 772        <NA>
#> 773        <NA>
#> 774        <NA>
#> 775        <NA>
#> 776        <NA>
#> 777        <NA>
#> 778        <NA>
#> 779        <NA>
#> 780        <NA>
#> 781        <NA>
#> 782   consensus
#> 783        <NA>
#> 784        <NA>
#> 785        <NA>
#> 786        <NA>
#> 787        <NA>
#> 788        <NA>
#> 789        plan
#> 790        <NA>
#> 791        <NA>
#> 792   consensus
#> 793        <NA>
#> 794        <NA>
#> 795        <NA>
#> 796        <NA>
#> 797        <NA>
#> 798        <NA>
#> 799        <NA>
#> 800        <NA>
#> 801        <NA>
#> 802        <NA>
#> 803        <NA>
#> 804        <NA>
#> 805        <NA>
#> 806        <NA>
#> 807        <NA>
#> 808        <NA>
#> 809        <NA>
#> 810        <NA>
#> 811        <NA>
#> 812        <NA>
#> 813        <NA>
#> 814        <NA>
#> 815        <NA>
#> 816        <NA>
#> 817        <NA>
#> 818        <NA>
#> 819        <NA>
#> 820        <NA>
#> 821        <NA>
#> 822        <NA>
#> 823        <NA>
#> 824        <NA>
#> 825        <NA>
#> 826        <NA>
#> 827        <NA>
#> 828        <NA>
#> 829        <NA>
#> 830        <NA>
#> 831        <NA>
#> 832        <NA>
#> 833        plan
#> 834        <NA>
#> 835        <NA>
#> 836        <NA>
#> 837        <NA>
#> 838        <NA>
#> 839        <NA>
#> 840        <NA>
#> 841        <NA>
#> 842        <NA>
#> 843        <NA>
#> 844        <NA>
#> 845        <NA>
#> 846        <NA>
#> 847        <NA>
#> 848        <NA>
#> 849        <NA>
#> 850        <NA>
#> 851   synthesis
#> 852        <NA>
#> 853        <NA>
#> 854        <NA>
#> 855        <NA>
#> 856        <NA>
#> 857        <NA>
#> 858        <NA>
#> 859        <NA>
#> 860        <NA>
#> 861        <NA>
#> 862        <NA>
#> 863        <NA>
#> 864        <NA>
#> 865        <NA>
#> 866        <NA>
#> 867        <NA>
#> 868        <NA>
#> 869        <NA>
#> 870        <NA>
#> 871        <NA>
#> 872        <NA>
#> 873        <NA>
#> 874        <NA>
#> 875        <NA>
#> 876        <NA>
#> 877        <NA>
#> 878        <NA>
#> 879        <NA>
#> 880        <NA>
#> 881        <NA>
#> 882        <NA>
#> 883        <NA>
#> 884        <NA>
#> 885        <NA>
#> 886        <NA>
#> 887        <NA>
#> 888        <NA>
#> 889        <NA>
#> 890        <NA>
#> 891        <NA>
#> 892        <NA>
#> 893        <NA>
#> 894        <NA>
#> 895        <NA>
#> 896        <NA>
#> 897        <NA>
#> 898        <NA>
#> 899        <NA>
#> 900        <NA>
#> 901        <NA>
#> 902        <NA>
#> 903        <NA>
#> 904        <NA>
#> 905        <NA>
#> 906        <NA>
#> 907        <NA>
#> 908        <NA>
#> 909        <NA>
#> 910        <NA>
#> 911        <NA>
#> 912        <NA>
#> 913        <NA>
#> 914        <NA>
#> 915        <NA>
#> 916        <NA>
#> 917        <NA>
#> 918        <NA>
#> 919        <NA>
#> 920        <NA>
#> 921        <NA>
#> 922        <NA>
#> 923        <NA>
#> 924        <NA>
#> 925        <NA>
#> 926        <NA>
#> 927        <NA>
#> 928        <NA>
#> 929        <NA>
#> 930        <NA>
#> 931        <NA>
#> 932        <NA>
#> 933        <NA>
#> 934        <NA>
#> 935        <NA>
#> 936        <NA>
#> 937        <NA>
#> 938        <NA>
#> 939        <NA>
#> 940        <NA>
#> 941        <NA>
#> 942        <NA>
#> 943        <NA>
#> 944        <NA>
#> 945        <NA>
#> 946        <NA>
#> 947        <NA>
#> 948        <NA>
#> 949        <NA>
#> 950        <NA>
#> 951        <NA>
#> 952        <NA>
#> 953        <NA>
#> 954        <NA>
#> 955        <NA>
#> 956        plan
#> 957        <NA>
#> 958        <NA>
#> 959        <NA>
#> 960        <NA>
#> 961        <NA>
#> 962        <NA>
#> 963        <NA>
#> 964        <NA>
#> 965        <NA>
#> 966        <NA>
#> 967        <NA>
#> 968        <NA>
#> 969        <NA>
#> 970        <NA>
#> 971        <NA>
#> 972        <NA>
#> 973        <NA>
#> 974        <NA>
#> 975        <NA>
#> 976        <NA>
#> 977        <NA>
#> 978        <NA>
#> 979        <NA>
#> 980        <NA>
#> 981        <NA>
#> 982        <NA>
#> 983        <NA>
#> 984        <NA>
#> 985        <NA>
#> 986        <NA>
#> 987        <NA>
#> 988        <NA>
#> 989        <NA>
#> 990        <NA>
#> 991        <NA>
#> 992        <NA>
#> 993        <NA>
#> 994        <NA>
#> 995        <NA>
#> 996        <NA>
#> 997        <NA>
#> 998        <NA>
#> 999        <NA>
#> 1000       <NA>
#> 1001       <NA>
#> 1002       <NA>
#> 1003       <NA>
#> 1004       <NA>
#> 1005       <NA>
#> 1006       <NA>
#> 1007       <NA>
#> 1008       <NA>
#> 1009       <NA>
#> 1010  consensus
#> 1011       <NA>
#> 1012       <NA>
#> 1013       <NA>
#> 1014       <NA>
#> 1015  consensus
#> 1016       <NA>
#> 1017       <NA>
#> 1018       <NA>
#> 1019       <NA>
#> 1020       <NA>
#> 1021       <NA>
#> 1022       <NA>
#> 1023       <NA>
#> 1024       <NA>
#> 1025       <NA>
#> 1026       <NA>
#> 1027       <NA>
#> 1028       <NA>
#> 1029       <NA>
#> 1030    discuss
#> 1031       <NA>
#> 1032       <NA>
#> 1033       <NA>
#> 1034       <NA>
#> 1035       <NA>
#> 1036       <NA>
#> 1037       <NA>
#> 1038       <NA>
#> 1039       <NA>
#> 1040       <NA>
#> 1041       <NA>
#> 1042       <NA>
#> 1043       <NA>
#> 1044       <NA>
#> 1045       <NA>
#> 1046       <NA>
#> 1047       <NA>
#> 1048       <NA>
#> 1049       <NA>
#> 1050       <NA>
#> 1051       <NA>
#> 1052       <NA>
#> 1053       <NA>
#> 1054       <NA>
#> 1055       <NA>
#> 1056       <NA>
#> 1057       <NA>
#> 1058       <NA>
#> 1059       <NA>
#> 1060       <NA>
#> 1061       <NA>
#> 1062       <NA>
#> 1063       <NA>
#> 1064       <NA>
#> 1065       <NA>
#> 1066       <NA>
#> 1067       <NA>
#> 1068       <NA>
#> 1069       <NA>
#> 1070       <NA>
#> 1071       <NA>
#> 1072       <NA>
#> 1073       <NA>
#> 1074       <NA>
#> 1075       <NA>
#> 1076       <NA>
#> 1077       <NA>
#> 1078       <NA>
#> 1079       <NA>
#> 1080       <NA>
#> 1081       <NA>
#> 1082       <NA>
#> 1083       <NA>
#> 1084       <NA>
#> 1085       <NA>
#> 1086       <NA>
#> 1087       <NA>
#> 1088       <NA>
#> 1089       <NA>
#> 1090       <NA>
#> 1091       <NA>
#> 1092  consensus
#> 1093       <NA>
#> 1094       <NA>
#> 1095       <NA>
#> 1096       <NA>
#> 1097       <NA>
#> 1098       <NA>
#> 1099       <NA>
#> 1100       <NA>
#> 1101       <NA>
#> 1102       <NA>
#> 1103       <NA>
#> 1104       <NA>
#> 1105       <NA>
#> 1106    discuss
#> 1107  consensus
#> 1108       <NA>
#> 1109       <NA>
#> 1110       <NA>
#> 1111       <NA>
#> 1112       <NA>
#> 1113       <NA>
#> 1114       <NA>
#> 1115       <NA>
#> 1116       <NA>
#> 1117       <NA>
#> 1118       <NA>
#> 1119       <NA>
#> 1120       <NA>
#> 1121       <NA>
#> 1122       <NA>
#> 1123       <NA>
#> 1124       <NA>
#> 1125       <NA>
#> 1126       <NA>
#> 1127       <NA>
#> 1128       <NA>
#> 1129       <NA>
#> 1130       <NA>
#> 1131       <NA>
#> 1132       <NA>
#> 1133       <NA>
#> 1134       <NA>
#> 1135       <NA>
#> 1136       <NA>
#> 1137       <NA>
#> 1138       <NA>
#> 1139       <NA>
#> 1140       <NA>
#> 1141       <NA>
#> 1142       <NA>
#> 1143       <NA>
#> 1144       <NA>
#> 1145       <NA>
#> 1146       <NA>
#> 1147       <NA>
#> 1148       <NA>
#> 1149       <NA>
#> 1150       <NA>
#> 1151       <NA>
#> 1152       <NA>
#> 1153       plan
#> 1154       <NA>
#> 1155       <NA>
#> 1156       <NA>
#> 1157       <NA>
#> 1158       <NA>
#> 1159       <NA>
#> 1160       <NA>
#> 1161       <NA>
#> 1162       <NA>
#> 1163       <NA>
#> 1164       <NA>
#> 1165       <NA>
#> 1166       <NA>
#> 1167       <NA>
#> 1168       <NA>
#> 1169       <NA>
#> 1170       <NA>
#> 1171       <NA>
#> 1172       <NA>
#> 1173       <NA>
#> 1174       <NA>
#> 1175       <NA>
#> 1176       <NA>
#> 1177       <NA>
#> 1178       <NA>
#> 1179       <NA>
#> 1180       <NA>
#> 1181       <NA>
#> 1182       <NA>
#> 1183       <NA>
#> 1184       plan
#> 1185       <NA>
#> 1186       <NA>
#> 1187       <NA>
#> 1188       <NA>
#> 1189       <NA>
#> 1190       <NA>
#> 1191       <NA>
#> 1192       <NA>
#> 1193       <NA>
#> 1194       <NA>
#> 1195       <NA>
#> 1196       <NA>
#> 1197       <NA>
#> 1198       <NA>
#> 1199       <NA>
#> 1200       <NA>
#> 1201       <NA>
#> 1202       <NA>
#> 1203       <NA>
#> 1204       <NA>
#> 1205       <NA>
#> 1206       <NA>
#> 1207       <NA>
#> 1208       <NA>
#> 1209 coregulate
#> 1210       <NA>
#> 1211       <NA>
#> 1212       <NA>
#> 1213       <NA>
#> 1214       <NA>
#> 1215       <NA>
#> 1216       <NA>
#> 1217       <NA>
#> 1218       <NA>
#> 1219       <NA>
#> 1220       <NA>
#> 1221       <NA>
#> 1222       <NA>
#> 1223       <NA>
#> 1224       <NA>
#> 1225       <NA>
#> 1226       <NA>
#> 1227       <NA>
#> 1228       <NA>
#> 1229       <NA>
#> 1230       <NA>
#> 1231       <NA>
#> 1232       <NA>
#> 1233       <NA>
#> 1234       <NA>
#> 1235       <NA>
#> 1236       <NA>
#> 1237       <NA>
#> 1238       <NA>
#> 1239       <NA>
#> 1240       <NA>
#> 1241       <NA>
#> 1242       <NA>
#> 1243       <NA>
#> 1244       <NA>
#> 1245       <NA>
#> 1246       <NA>
#> 1247       <NA>
#> 1248       <NA>
#> 1249       <NA>
#> 1250       <NA>
#> 1251       <NA>
#> 1252       <NA>
#> 1253       <NA>
#> 1254       <NA>
#> 1255       <NA>
#> 1256       <NA>
#> 1257       <NA>
#> 1258       <NA>
#> 1259       <NA>
#> 1260       <NA>
#> 1261       <NA>
#> 1262       <NA>
#> 1263       <NA>
#> 1264       <NA>
#> 1265       <NA>
#> 1266       <NA>
#> 1267 coregulate
#> 1268       <NA>
#> 1269       <NA>
#> 1270       <NA>
#> 1271       <NA>
#> 1272  consensus
#> 1273       <NA>
#> 1274       <NA>
#> 1275       <NA>
#> 1276       <NA>
#> 1277       <NA>
#> 1278       <NA>
#> 1279       <NA>
#> 1280       <NA>
#> 1281       <NA>
#> 1282       <NA>
#> 1283       <NA>
#> 1284       <NA>
#> 1285       <NA>
#> 1286       <NA>
#> 1287       <NA>
#> 1288       <NA>
#> 1289       <NA>
#> 1290       <NA>
#> 1291       <NA>
#> 1292       <NA>
#> 1293       <NA>
#> 1294       <NA>
#> 1295       <NA>
#> 1296       <NA>
#> 1297       <NA>
#> 1298       <NA>
#> 1299       <NA>
#> 1300       <NA>
#> 1301       <NA>
#> 1302       plan
#> 1303       <NA>
#> 1304       <NA>
#> 1305       <NA>
#> 1306       <NA>
#> 1307       <NA>
#> 1308       <NA>
#> 1309       <NA>
#> 1310       <NA>
#> 1311       <NA>
#> 1312       <NA>
#> 1313       <NA>
#> 1314       <NA>
#> 1315       <NA>
#> 1316       <NA>
#> 1317       <NA>
#> 1318       <NA>
#> 1319       <NA>
#> 1320       <NA>
#> 1321       <NA>
#> 1322       <NA>
#> 1323       <NA>
#> 1324       <NA>
#> 1325       <NA>
#> 1326       <NA>
#> 1327       <NA>
#> 1328       <NA>
#> 1329       <NA>
#> 1330       <NA>
#> 1331       <NA>
#> 1332       <NA>
#> 1333       <NA>
#> 1334   cohesion
#> 1335       <NA>
#> 1336       <NA>
#> 1337       <NA>
#> 1338       <NA>
#> 1339       <NA>
#> 1340       <NA>
#> 1341       <NA>
#> 1342       <NA>
#> 1343       <NA>
#> 1344       <NA>
#> 1345  synthesis
#> 1346       <NA>
#> 1347       <NA>
#> 1348       <NA>
#> 1349       <NA>
#> 1350       <NA>
#> 1351       <NA>
#> 1352       <NA>
#> 1353       <NA>
#> 1354       <NA>
#> 1355       <NA>
#> 1356       <NA>
#> 1357       <NA>
#> 1358       <NA>
#> 1359       <NA>
#> 1360       <NA>
#> 1361       <NA>
#> 1362  consensus
#> 1363       <NA>
#> 1364       <NA>
#> 1365       <NA>
#> 1366       <NA>
#> 1367       <NA>
#> 1368       <NA>
#> 1369       <NA>
#> 1370       <NA>
#> 1371       <NA>
#> 1372       <NA>
#> 1373       <NA>
#> 1374       <NA>
#> 1375       <NA>
#> 1376       <NA>
#> 1377       <NA>
#> 1378       <NA>
#> 1379       <NA>
#> 1380       <NA>
#> 1381       <NA>
#> 1382       <NA>
#> 1383       <NA>
#> 1384       plan
#> 1385       <NA>
#> 1386       <NA>
#> 1387       <NA>
#> 1388       <NA>
#> 1389       <NA>
#> 1390       <NA>
#> 1391       <NA>
#> 1392       <NA>
#> 1393       <NA>
#> 1394       <NA>
#> 1395       <NA>
#> 1396       <NA>
#> 1397       <NA>
#> 1398       <NA>
#> 1399       <NA>
#> 1400       <NA>
#> 1401 coregulate
#> 1402       <NA>
#> 1403       <NA>
#> 1404       <NA>
#> 1405       <NA>
#> 1406       <NA>
#> 1407       <NA>
#> 1408       <NA>
#> 1409      adapt
#> 1410       <NA>
#> 1411       <NA>
#> 1412       <NA>
#> 1413       <NA>
#> 1414       <NA>
#> 1415       <NA>
#> 1416       <NA>
#> 1417       <NA>
#> 1418       <NA>
#> 1419       <NA>
#> 1420       <NA>
#> 1421       <NA>
#> 1422       <NA>
#> 1423       <NA>
#> 1424       <NA>
#> 1425       <NA>
#> 1426       <NA>
#> 1427       <NA>
#> 1428       <NA>
#> 1429       <NA>
#> 1430       <NA>
#> 1431       <NA>
#> 1432       <NA>
#> 1433       <NA>
#> 1434       <NA>
#> 1435       <NA>
#> 1436       <NA>
#> 1437       <NA>
#> 1438       <NA>
#> 1439       <NA>
#> 1440       <NA>
#> 1441       <NA>
#> 1442       <NA>
#> 1443       <NA>
#> 1444       <NA>
#> 1445       <NA>
#> 1446    discuss
#> 1447       <NA>
#> 1448    emotion
#> 1449       <NA>
#> 1450       <NA>
#> 1451       <NA>
#> 1452       <NA>
#> 1453       <NA>
#> 1454       <NA>
#> 1455       <NA>
#> 1456       <NA>
#> 1457       <NA>
#> 1458       <NA>
#> 1459       <NA>
#> 1460       <NA>
#> 1461       <NA>
#> 1462       <NA>
#> 1463       <NA>
#> 1464       <NA>
#> 1465       <NA>
#> 1466       <NA>
#> 1467       <NA>
#> 1468       <NA>
#> 1469       <NA>
#> 1470       <NA>
#> 1471       <NA>
#> 1472       <NA>
#> 1473       <NA>
#> 1474       <NA>
#> 1475       <NA>
#> 1476       <NA>
#> 1477       <NA>
#> 1478       <NA>
#> 1479       <NA>
#> 1480       <NA>
#> 1481       <NA>
#> 1482       <NA>
#> 1483       <NA>
#> 1484       <NA>
#> 1485       <NA>
#> 1486       <NA>
#> 1487       <NA>
#> 1488       <NA>
#> 1489       <NA>
#> 1490       <NA>
#> 1491       <NA>
#> 1492       <NA>
#> 1493       <NA>
#> 1494       <NA>
#> 1495       <NA>
#> 1496       <NA>
#> 1497    emotion
#> 1498       <NA>
#> 1499       <NA>
#> 1500       <NA>
#> 1501       <NA>
#> 1502       <NA>
#> 1503       <NA>
#> 1504       <NA>
#> 1505       <NA>
#> 1506       <NA>
#> 1507       <NA>
#> 1508       <NA>
#> 1509       <NA>
#> 1510       <NA>
#> 1511       <NA>
#> 1512       <NA>
#> 1513       <NA>
#> 1514       <NA>
#> 1515       <NA>
#> 1516       <NA>
#> 1517       <NA>
#> 1518       <NA>
#> 1519       <NA>
#> 1520       <NA>
#> 1521       <NA>
#> 1522       <NA>
#> 1523       <NA>
#> 1524       <NA>
#> 1525       plan
#> 1526       <NA>
#> 1527       <NA>
#> 1528    discuss
#> 1529       <NA>
#> 1530       <NA>
#> 1531       <NA>
#> 1532       <NA>
#> 1533       <NA>
#> 1534       <NA>
#> 1535       <NA>
#> 1536       <NA>
#> 1537       <NA>
#> 1538       <NA>
#> 1539  synthesis
#> 1540       <NA>
#> 1541       <NA>
#> 1542       <NA>
#> 1543       <NA>
#> 1544       <NA>
#> 1545       <NA>
#> 1546       <NA>
#> 1547       <NA>
#> 1548       <NA>
#> 1549       <NA>
#> 1550       <NA>
#> 1551       <NA>
#> 1552       <NA>
#> 1553       <NA>
#> 1554       <NA>
#> 1555       <NA>
#> 1556       <NA>
#> 1557       <NA>
#> 1558       <NA>
#> 1559       <NA>
#> 1560       <NA>
#> 1561       <NA>
#> 1562       <NA>
#> 1563       <NA>
#> 1564       <NA>
#> 1565       <NA>
#> 1566       <NA>
#> 1567       <NA>
#> 1568       <NA>
#> 1569       <NA>
#> 1570       <NA>
#> 1571       <NA>
#> 1572       <NA>
#> 1573       <NA>
#> 1574       <NA>
#> 1575       <NA>
#> 1576       <NA>
#> 1577       <NA>
#> 1578       <NA>
#> 1579       <NA>
#> 1580       <NA>
#> 1581       <NA>
#> 1582       <NA>
#> 1583       <NA>
#> 1584       <NA>
#> 1585       <NA>
#> 1586       <NA>
#> 1587       <NA>
#> 1588       <NA>
#> 1589       <NA>
#> 1590 coregulate
#> 1591       <NA>
#> 1592       <NA>
#> 1593       <NA>
#> 1594       <NA>
#> 1595       <NA>
#> 1596       <NA>
#> 1597       <NA>
#> 1598       <NA>
#> 1599       <NA>
#> 1600       <NA>
#> 1601       <NA>
#> 1602       <NA>
#> 1603       <NA>
#> 1604       <NA>
#> 1605       <NA>
#> 1606       <NA>
#> 1607       <NA>
#> 1608       <NA>
#> 1609       <NA>
#> 1610       <NA>
#> 1611   cohesion
#> 1612       <NA>
#> 1613       <NA>
#> 1614       <NA>
#> 1615       <NA>
#> 1616       <NA>
#> 1617       <NA>
#> 1618       <NA>
#> 1619       <NA>
#> 1620       <NA>
#> 1621       <NA>
#> 1622       <NA>
#> 1623       <NA>
#> 1624       <NA>
#> 1625       <NA>
#> 1626       <NA>
#> 1627       <NA>
#> 1628       <NA>
#> 1629       <NA>
#> 1630       <NA>
#> 1631       <NA>
#> 1632       <NA>
#> 1633       <NA>
#> 1634  consensus
#> 1635       <NA>
#> 1636       <NA>
#> 1637       <NA>
#> 1638       <NA>
#> 1639       <NA>
#> 1640       <NA>
#> 1641       <NA>
#> 1642       <NA>
#> 1643       <NA>
#> 1644       <NA>
#> 1645       <NA>
#> 1646       <NA>
#> 1647       <NA>
#> 1648       <NA>
#> 1649       <NA>
#> 1650       <NA>
#> 1651       <NA>
#> 1652       <NA>
#> 1653       <NA>
#> 1654       <NA>
#> 1655       <NA>
#> 1656       <NA>
#> 1657       <NA>
#> 1658       <NA>
#> 1659       <NA>
#> 1660    discuss
#> 1661       <NA>
#> 1662       <NA>
#> 1663       <NA>
#> 1664       <NA>
#> 1665       <NA>
#> 1666       <NA>
#> 1667       <NA>
#> 1668       <NA>
#> 1669       <NA>
#> 1670       <NA>
#> 1671       <NA>
#> 1672  consensus
#> 1673       <NA>
#> 1674       <NA>
#> 1675       <NA>
#> 1676       <NA>
#> 1677       <NA>
#> 1678       <NA>
#> 1679       <NA>
#> 1680       <NA>
#> 1681       <NA>
#> 1682       <NA>
#> 1683       <NA>
#> 1684       <NA>
#> 1685       <NA>
#> 1686       <NA>
#> 1687       <NA>
#> 1688       <NA>
#> 1689       <NA>
#> 1690       <NA>
#> 1691       <NA>
#> 1692       <NA>
#> 1693       <NA>
#> 1694       <NA>
#> 1695       <NA>
#> 1696       <NA>
#> 1697       plan
#> 1698       <NA>
#> 1699       <NA>
#> 1700       <NA>
#> 1701       <NA>
#> 1702       <NA>
#> 1703       <NA>
#> 1704       <NA>
#> 1705       <NA>
#> 1706       <NA>
#> 1707       <NA>
#> 1708       <NA>
#> 1709       <NA>
#> 1710       <NA>
#> 1711       <NA>
#> 1712       <NA>
#> 1713       <NA>
#> 1714       <NA>
#> 1715       <NA>
#> 1716       <NA>
#> 1717       <NA>
#> 1718       <NA>
#> 1719  consensus
#> 1720       <NA>
#> 1721       <NA>
#> 1722  consensus
#> 1723       <NA>
#> 1724       <NA>
#> 1725       <NA>
#> 1726       <NA>
#> 1727       <NA>
#> 1728       <NA>
#> 1729       <NA>
#> 1730       <NA>
#> 1731       <NA>
#> 1732       <NA>
#> 1733       <NA>
#> 1734       <NA>
#> 1735  synthesis
#> 1736       <NA>
#> 1737       <NA>
#> 1738       <NA>
#> 1739       <NA>
#> 1740       <NA>
#> 1741       <NA>
#> 1742       <NA>
#> 1743       <NA>
#> 1744       <NA>
#> 1745       <NA>
#> 1746       <NA>
#> 1747       <NA>
#> 1748       <NA>
#> 1749       <NA>
#> 1750       <NA>
#> 1751       <NA>
#> 1752       <NA>
#> 1753       <NA>
#> 1754       <NA>
#> 1755       <NA>
#> 1756       <NA>
#> 1757       <NA>
#> 1758       <NA>
#> 1759       <NA>
#> 1760       <NA>
#> 1761       <NA>
#> 1762       <NA>
#> 1763       <NA>
#> 1764       <NA>
#> 1765       <NA>
#> 1766       <NA>
#> 1767       <NA>
#> 1768       <NA>
#> 1769       <NA>
#> 1770       <NA>
#> 1771       <NA>
#> 1772       <NA>
#> 1773    emotion
#> 1774       <NA>
#> 1775       <NA>
#> 1776       <NA>
#> 1777       <NA>
#> 1778       <NA>
#> 1779       <NA>
#> 1780       <NA>
#> 1781       <NA>
#> 1782       <NA>
#> 1783       <NA>
#> 1784       <NA>
#> 1785       <NA>
#> 1786       <NA>
#> 1787       <NA>
#> 1788       <NA>
#> 1789       <NA>
#> 1790       <NA>
#> 1791       <NA>
#> 1792       <NA>
#> 1793       <NA>
#> 1794       <NA>
#> 1795       <NA>
#> 1796       <NA>
#> 1797       <NA>
#> 1798       <NA>
#> 1799       <NA>
#> 1800       <NA>
#> 1801       <NA>
#> 1802       <NA>
#> 1803       <NA>
#> 1804       <NA>
#> 1805       <NA>
#> 1806       <NA>
#> 1807       <NA>
#> 1808       <NA>
#> 1809       <NA>
#> 1810       <NA>
#> 1811       <NA>
#> 1812       <NA>
#> 1813       <NA>
#> 1814       <NA>
#> 1815       <NA>
#> 1816       <NA>
#> 1817       <NA>
#> 1818       <NA>
#> 1819       <NA>
#> 1820       <NA>
#> 1821       <NA>
#> 1822       <NA>
#> 1823       <NA>
#> 1824       <NA>
#> 1825       <NA>
#> 1826       <NA>
#> 1827       <NA>
#> 1828       <NA>
#> 1829       <NA>
#> 1830       <NA>
#> 1831       <NA>
#> 1832       <NA>
#> 1833       <NA>
#> 1834       <NA>
#> 1835       <NA>
#> 1836       <NA>
#> 1837       <NA>
#> 1838       <NA>
#> 1839       <NA>
#> 1840       plan
#> 1841       <NA>
#> 1842       <NA>
#> 1843       <NA>
#> 1844       <NA>
#> 1845       <NA>
#> 1846       <NA>
#> 1847       <NA>
#> 1848       <NA>
#> 1849       <NA>
#> 1850       <NA>
#> 1851       <NA>
#> 1852       <NA>
#> 1853       <NA>
#> 1854       <NA>
#> 1855       <NA>
#> 1856       <NA>
#> 1857       <NA>
#> 1858       <NA>
#> 1859       <NA>
#> 1860       <NA>
#> 1861       <NA>
#> 1862       <NA>
#> 1863       <NA>
#> 1864       <NA>
#> 1865       <NA>
#> 1866       <NA>
#> 1867       <NA>
#> 1868       <NA>
#> 1869       <NA>
#> 1870       <NA>
#> 1871       <NA>
#> 1872       <NA>
#> 1873       <NA>
#> 1874       <NA>
#> 1875       <NA>
#> 1876       <NA>
#> 1877       <NA>
#> 1878    emotion
#> 1879  consensus
#> 1880       <NA>
#> 1881       <NA>
#> 1882       <NA>
#> 1883       <NA>
#> 1884       <NA>
#> 1885       <NA>
#> 1886       <NA>
#> 1887       <NA>
#> 1888       <NA>
#> 1889       <NA>
#> 1890       <NA>
#> 1891       <NA>
#> 1892       <NA>
#> 1893       <NA>
#> 1894       <NA>
#> 1895       <NA>
#> 1896       <NA>
#> 1897       <NA>
#> 1898       <NA>
#> 1899       <NA>
#> 1900       <NA>
#> 1901       <NA>
#> 1902       <NA>
#> 1903       <NA>
#> 1904       <NA>
#> 1905       <NA>
#> 1906       <NA>
#> 1907       <NA>
#> 1908       <NA>
#> 1909       <NA>
#> 1910       <NA>
#> 1911    discuss
#> 1912       <NA>
#> 1913       <NA>
#> 1914       <NA>
#> 1915       <NA>
#> 1916       <NA>
#> 1917       <NA>
#> 1918       <NA>
#> 1919       <NA>
#> 1920       <NA>
#> 1921       <NA>
#> 1922       <NA>
#> 1923       <NA>
#> 1924       <NA>
#> 1925       <NA>
#> 1926       <NA>
#> 1927       <NA>
#> 1928       <NA>
#> 1929       <NA>
#> 1930       <NA>
#> 1931       <NA>
#> 1932       <NA>
#> 1933       <NA>
#> 1934       <NA>
#> 1935       <NA>
#> 1936    discuss
#> 1937       <NA>
#> 1938       <NA>
#> 1939       <NA>
#> 1940       <NA>
#> 1941       <NA>
#> 1942       <NA>
#> 1943       <NA>
#> 1944       <NA>
#> 1945       <NA>
#> 1946       <NA>
#> 1947       <NA>
#> 1948       <NA>
#> 1949       <NA>
#> 1950       <NA>
#> 1951       <NA>
#> 1952       <NA>
#> 1953       <NA>
#> 1954       <NA>
#> 1955       <NA>
#> 1956       <NA>
#> 1957       <NA>
#> 1958       <NA>
#> 1959       <NA>
#> 1960       <NA>
#> 1961       <NA>
#> 1962  consensus
#> 1963       <NA>
#> 1964       <NA>
#> 1965       <NA>
#> 1966       <NA>
#> 1967       <NA>
#> 1968       <NA>
#> 1969       <NA>
#> 1970       <NA>
#> 1971       <NA>
#> 1972       <NA>
#> 1973       <NA>
#> 1974       <NA>
#> 1975       <NA>
#> 1976       <NA>
#> 1977       <NA>
#> 1978       <NA>
#> 1979       <NA>
#> 1980       <NA>
#> 1981       <NA>
#> 1982       <NA>
#> 1983       <NA>
#> 1984       <NA>
#> 1985       <NA>
#> 1986       <NA>
#> 1987       <NA>
#> 1988    emotion
#> 1989       <NA>
#> 1990       <NA>
#> 1991       <NA>
#> 1992       <NA>
#> 1993       <NA>
#> 1994       <NA>
#> 1995       <NA>
#> 1996       <NA>
#> 1997       <NA>
#> 1998       <NA>
#> 1999       <NA>
#> 2000       <NA>

net <- build_tna(trajectories)
extract_sequences(net)
#>              1          2          3          4          5          6
#> 1       Active Disengaged Disengaged Disengaged     Active     Active
#> 2      Average    Average    Average    Average    Average    Average
#> 3      Average     Active     Active     Active     Active     Active
#> 4       Active     Active     Active     Active     Active    Average
#> 5       Active     Active     Active    Average     Active    Average
#> 6      Average    Average Disengaged    Average Disengaged Disengaged
#> 7       Active    Average    Average    Average    Average    Average
#> 8      Average Disengaged    Average Disengaged Disengaged Disengaged
#> 9       Active     Active     Active     Active     Active     Active
#> 10     Average    Average    Average    Average    Average    Average
#> 11     Average    Average    Average     Active    Average     Active
#> 12     Average    Average    Average    Average    Average Disengaged
#> 13      Active     Active     Active     Active     Active     Active
#> 14      Active     Active     Active     Active     Active     Active
#> 15     Average    Average Disengaged Disengaged    Average Disengaged
#> 16     Average     Active     Active    Average    Average    Average
#> 17      Active Disengaged     Active Disengaged Disengaged Disengaged
#> 18     Average    Average    Average    Average    Average     Active
#> 19      Active    Average Disengaged    Average Disengaged    Average
#> 20     Average    Average    Average    Average    Average    Average
#> 21      Active Disengaged    Average Disengaged    Average Disengaged
#> 22     Average    Average     Active     Active       <NA>       <NA>
#> 23     Average Disengaged Disengaged     Active    Average    Average
#> 24      Active     Active     Active    Average    Average    Average
#> 25      Active     Active     Active     Active     Active     Active
#> 26      Active     Active     Active     Active     Active     Active
#> 27  Disengaged    Average     Active    Average     Active    Average
#> 28  Disengaged    Average    Average    Average Disengaged    Average
#> 29  Disengaged    Average    Average    Average    Average     Active
#> 30  Disengaged Disengaged       <NA>       <NA>       <NA>       <NA>
#> 31  Disengaged    Average     Active    Average    Average    Average
#> 32      Active     Active     Active     Active     Active     Active
#> 33  Disengaged       <NA>       <NA>       <NA>       <NA>       <NA>
#> 34      Active     Active     Active     Active     Active    Average
#> 35  Disengaged Disengaged Disengaged    Average    Average    Average
#> 36     Average    Average    Average    Average Disengaged    Average
#> 37      Active     Active    Average     Active    Average    Average
#> 38     Average    Average Disengaged    Average Disengaged     Active
#> 39      Active    Average     Active     Active     Active     Active
#> 40     Average     Active     Active     Active    Average     Active
#> 41     Average Disengaged Disengaged Disengaged Disengaged     Active
#> 42      Active Disengaged     Active     Active     Active     Active
#> 43      Active    Average    Average    Average    Average    Average
#> 44     Average    Average    Average    Average    Average    Average
#> 45     Average Disengaged Disengaged Disengaged Disengaged Disengaged
#> 46     Average Disengaged Disengaged Disengaged Disengaged Disengaged
#> 47      Active     Active     Active     Active     Active    Average
#> 48  Disengaged Disengaged Disengaged    Average Disengaged Disengaged
#> 49     Average     Active    Average    Average    Average     Active
#> 50      Active     Active     Active     Active     Active     Active
#> 51     Average Disengaged Disengaged Disengaged Disengaged Disengaged
#> 52  Disengaged    Average    Average    Average    Average    Average
#> 53     Average    Average Disengaged Disengaged Disengaged     Active
#> 54  Disengaged Disengaged Disengaged Disengaged Disengaged    Average
#> 55      Active    Average    Average     Active     Active    Average
#> 56     Average    Average    Average    Average    Average    Average
#> 57     Average Disengaged Disengaged    Average    Average Disengaged
#> 58     Average    Average    Average    Average    Average    Average
#> 59  Disengaged    Average    Average    Average Disengaged    Average
#> 60     Average Disengaged Disengaged Disengaged    Average Disengaged
#> 61     Average    Average Disengaged    Average    Average    Average
#> 62  Disengaged     Active    Average    Average    Average    Average
#> 63  Disengaged    Average    Average    Average Disengaged    Average
#> 64      Active Disengaged Disengaged Disengaged       <NA>       <NA>
#> 65     Average    Average    Average    Average    Average Disengaged
#> 66      Active     Active     Active     Active     Active     Active
#> 67     Average    Average    Average    Average    Average    Average
#> 68      Active     Active     Active     Active     Active     Active
#> 69     Average    Average Disengaged    Average    Average    Average
#> 70  Disengaged    Average Disengaged Disengaged Disengaged Disengaged
#> 71  Disengaged    Average    Average     Active     Active    Average
#> 72  Disengaged    Average Disengaged Disengaged Disengaged Disengaged
#> 73      Active    Average    Average     Active     Active    Average
#> 74  Disengaged    Average    Average    Average    Average    Average
#> 75      Active     Active     Active     Active     Active     Active
#> 76      Active     Active    Average    Average     Active    Average
#> 77        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 78      Active     Active     Active    Average    Average     Active
#> 79      Active     Active     Active    Average     Active     Active
#> 80     Average     Active     Active     Active     Active     Active
#> 81     Average     Active    Average     Active    Average     Active
#> 82      Active     Active     Active     Active     Active     Active
#> 83     Average    Average    Average    Average    Average Disengaged
#> 84      Active    Average    Average    Average     Active    Average
#> 85     Average     Active     Active     Active    Average     Active
#> 86      Active     Active     Active     Active     Active     Active
#> 87      Active     Active     Active    Average     Active     Active
#> 88     Average     Active    Average     Active     Active Disengaged
#> 89      Active     Active     Active     Active     Active     Active
#> 90     Average     Active     Active     Active    Average     Active
#> 91     Average    Average    Average     Active    Average    Average
#> 92      Active     Active     Active     Active     Active     Active
#> 93     Average    Average    Average Disengaged    Average    Average
#> 94      Active     Active     Active     Active     Active     Active
#> 95  Disengaged    Average Disengaged     Active    Average    Average
#> 96  Disengaged    Average    Average    Average Disengaged    Average
#> 97  Disengaged       <NA>       <NA>       <NA>       <NA>       <NA>
#> 98  Disengaged       <NA>       <NA>       <NA>       <NA>       <NA>
#> 99     Average Disengaged Disengaged Disengaged       <NA>       <NA>
#> 100 Disengaged    Average     Active    Average    Average    Average
#> 101 Disengaged Disengaged Disengaged Disengaged Disengaged Disengaged
#> 102 Disengaged       <NA>       <NA>       <NA>       <NA>       <NA>
#> 103 Disengaged Disengaged       <NA>       <NA>       <NA>       <NA>
#> 104       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 105 Disengaged Disengaged Disengaged    Average       <NA>       <NA>
#> 106 Disengaged       <NA>       <NA>       <NA>       <NA>       <NA>
#> 107     Active    Average    Average     Active     Active    Average
#> 108     Active    Average    Average     Active     Active    Average
#> 109 Disengaged    Average    Average    Average Disengaged    Average
#> 110    Average    Average    Average    Average    Average Disengaged
#> 111     Active    Average Disengaged    Average Disengaged    Average
#> 112     Active    Average    Average     Active     Active    Average
#> 113    Average    Average Disengaged Disengaged Disengaged     Active
#> 114 Disengaged    Average     Active    Average     Active    Average
#> 115    Average     Active     Active     Active     Active     Active
#> 116     Active    Average    Average    Average    Average    Average
#> 117     Active Disengaged    Average Disengaged    Average Disengaged
#> 118 Disengaged    Average    Average    Average    Average     Active
#> 119     Active     Active     Active     Active     Active    Average
#> 120     Active     Active     Active     Active     Active    Average
#> 121 Disengaged    Average     Active    Average    Average    Average
#> 122     Active Disengaged     Active     Active     Active     Active
#> 123    Average    Average Disengaged    Average Disengaged Disengaged
#> 124     Active     Active     Active     Active     Active     Active
#> 125    Average     Active     Active     Active    Average     Active
#> 126     Active     Active     Active     Active     Active    Average
#> 127    Average    Average    Average    Average    Average Disengaged
#> 128 Disengaged Disengaged Disengaged    Average Disengaged Disengaged
#> 129 Disengaged    Average     Active    Average     Active    Average
#> 130    Average    Average    Average    Average Disengaged    Average
#> 131    Average     Active    Average    Average    Average     Active
#> 132    Average     Active     Active    Average    Average    Average
#> 133     Active     Active     Active     Active     Active     Active
#> 134    Average    Average    Average    Average    Average Disengaged
#> 135     Active     Active     Active     Active     Active     Active
#> 136     Active Disengaged     Active     Active     Active     Active
#> 137     Active     Active    Average    Average     Active    Average
#> 138 Disengaged    Average     Active    Average    Average    Average
#>              7          8          9         10         11         12
#> 1       Active     Active    Average     Active    Average     Active
#> 2      Average    Average    Average     Active    Average    Average
#> 3       Active     Active     Active    Average     Active     Active
#> 4      Average     Active     Active     Active    Average     Active
#> 5       Active     Active     Active Disengaged     Active     Active
#> 6      Average    Average    Average Disengaged    Average Disengaged
#> 7      Average    Average    Average Disengaged    Average    Average
#> 8       Active     Active     Active     Active     Active     Active
#> 9      Average     Active    Average    Average     Active     Active
#> 10      Active    Average    Average Disengaged    Average    Average
#> 11      Active     Active     Active    Average    Average     Active
#> 12     Average Disengaged    Average    Average Disengaged Disengaged
#> 13     Average    Average     Active     Active     Active     Active
#> 14      Active     Active     Active     Active     Active     Active
#> 15     Average Disengaged    Average    Average Disengaged Disengaged
#> 16     Average    Average    Average     Active    Average    Average
#> 17  Disengaged     Active Disengaged    Average    Average     Active
#> 18      Active    Average    Average    Average Disengaged    Average
#> 19      Active     Active    Average    Average    Average    Average
#> 20      Active    Average    Average    Average Disengaged    Average
#> 21     Average    Average    Average    Average    Average    Average
#> 22        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 23     Average     Active    Average     Active    Average    Average
#> 24     Average    Average    Average       <NA>       <NA>       <NA>
#> 25     Average     Active     Active    Average    Average    Average
#> 26      Active     Active Disengaged     Active     Active     Active
#> 27      Active    Average Disengaged     Active     Active    Average
#> 28  Disengaged    Average Disengaged    Average    Average    Average
#> 29  Disengaged    Average Disengaged    Average    Average     Active
#> 30        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 31  Disengaged    Average Disengaged    Average    Average    Average
#> 32      Active     Active     Active     Active     Active     Active
#> 33        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 34  Disengaged    Average Disengaged    Average     Active    Average
#> 35     Average       <NA>       <NA>       <NA>       <NA>       <NA>
#> 36  Disengaged    Average Disengaged     Active    Average    Average
#> 37  Disengaged     Active     Active     Active     Active     Active
#> 38     Average     Active    Average     Active    Average    Average
#> 39      Active     Active     Active    Average    Average     Active
#> 40     Average     Active     Active    Average    Average     Active
#> 41      Active Disengaged    Average    Average     Active    Average
#> 42      Active     Active     Active     Active     Active     Active
#> 43     Average    Average     Active     Active     Active     Active
#> 44     Average    Average    Average     Active    Average    Average
#> 45     Average Disengaged    Average Disengaged    Average Disengaged
#> 46     Average Disengaged    Average Disengaged Disengaged Disengaged
#> 47     Average    Average    Average     Active     Active    Average
#> 48  Disengaged Disengaged Disengaged Disengaged Disengaged Disengaged
#> 49      Active Disengaged Disengaged     Active    Average Disengaged
#> 50      Active Disengaged Disengaged    Average    Average     Active
#> 51  Disengaged Disengaged Disengaged Disengaged Disengaged Disengaged
#> 52     Average Disengaged Disengaged    Average    Average     Active
#> 53     Average Disengaged Disengaged    Average    Average Disengaged
#> 54  Disengaged     Active Disengaged Disengaged Disengaged    Average
#> 55     Average Disengaged Disengaged    Average    Average Disengaged
#> 56     Average Disengaged Disengaged Disengaged    Average Disengaged
#> 57  Disengaged Disengaged Disengaged Disengaged Disengaged    Average
#> 58     Average    Average Disengaged    Average Disengaged     Active
#> 59     Average     Active    Average    Average Disengaged    Average
#> 60  Disengaged    Average    Average Disengaged Disengaged    Average
#> 61  Disengaged     Active     Active    Average Disengaged     Active
#> 62      Active     Active    Average    Average    Average    Average
#> 63     Average    Average    Average    Average    Average    Average
#> 64        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 65     Average    Average    Average Disengaged Disengaged    Average
#> 66     Average     Active    Average     Active    Average     Active
#> 67  Disengaged    Average    Average Disengaged    Average    Average
#> 68      Active     Active    Average    Average    Average     Active
#> 69     Average    Average    Average    Average Disengaged    Average
#> 70  Disengaged    Average    Average Disengaged Disengaged Disengaged
#> 71     Average    Average    Average    Average    Average    Average
#> 72  Disengaged    Average    Average Disengaged Disengaged    Average
#> 73     Average    Average    Average     Active     Active     Active
#> 74     Average     Active     Active     Active    Average     Active
#> 75      Active     Active     Active     Active     Active     Active
#> 76      Active     Active    Average     Active     Active     Active
#> 77        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 78      Active    Average     Active     Active     Active     Active
#> 79      Active     Active     Active     Active     Active     Active
#> 80      Active     Active     Active     Active     Active     Active
#> 81      Active    Average     Active     Active     Active     Active
#> 82      Active     Active     Active     Active     Active     Active
#> 83     Average     Active    Average    Average    Average     Active
#> 84     Average    Average     Active    Average    Average     Active
#> 85      Active    Average    Average     Active     Active     Active
#> 86      Active     Active     Active     Active     Active     Active
#> 87      Active     Active     Active     Active     Active     Active
#> 88     Average    Average    Average    Average    Average       <NA>
#> 89      Active     Active     Active     Active     Active     Active
#> 90      Active    Average     Active    Average     Active     Active
#> 91     Average    Average    Average    Average     Active     Active
#> 92      Active     Active     Active     Active     Active     Active
#> 93      Active     Active    Average    Average    Average     Active
#> 94      Active     Active     Active     Active     Active     Active
#> 95      Active    Average    Average     Active    Average    Average
#> 96  Disengaged    Average Disengaged    Average Disengaged    Average
#> 97        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 98        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 99        <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 100    Average    Average    Average    Average    Average    Average
#> 101 Disengaged    Average       <NA>       <NA>       <NA>       <NA>
#> 102       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 103       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 104       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 105       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 106       <NA>       <NA>       <NA>       <NA>       <NA>       <NA>
#> 107    Average Disengaged Disengaged    Average    Average Disengaged
#> 108    Average    Average    Average     Active     Active     Active
#> 109 Disengaged    Average Disengaged    Average    Average    Average
#> 110    Average     Active    Average    Average    Average     Active
#> 111     Active     Active    Average    Average    Average    Average
#> 112    Average Disengaged Disengaged    Average    Average Disengaged
#> 113    Average Disengaged Disengaged    Average    Average Disengaged
#> 114     Active    Average Disengaged     Active     Active    Average
#> 115     Active     Active     Active     Active     Active     Active
#> 116    Average    Average     Active     Active     Active     Active
#> 117    Average    Average    Average    Average    Average    Average
#> 118 Disengaged    Average Disengaged    Average    Average     Active
#> 119    Average     Active     Active     Active    Average     Active
#> 120    Average    Average    Average     Active     Active    Average
#> 121 Disengaged    Average Disengaged    Average    Average    Average
#> 122     Active     Active     Active     Active     Active     Active
#> 123    Average    Average    Average Disengaged    Average Disengaged
#> 124     Active     Active     Active     Active     Active     Active
#> 125    Average     Active     Active    Average    Average     Active
#> 126    Average     Active     Active     Active    Average     Active
#> 127    Average    Average    Average Disengaged Disengaged    Average
#> 128 Disengaged Disengaged Disengaged Disengaged Disengaged Disengaged
#> 129     Active    Average Disengaged     Active     Active    Average
#> 130 Disengaged    Average Disengaged     Active    Average    Average
#> 131     Active Disengaged Disengaged     Active    Average Disengaged
#> 132    Average    Average    Average     Active    Average    Average
#> 133    Average     Active     Active    Average    Average    Average
#> 134    Average    Average    Average Disengaged Disengaged    Average
#> 135    Average     Active    Average    Average     Active     Active
#> 136     Active     Active     Active     Active     Active     Active
#> 137     Active     Active    Average     Active     Active     Active
#> 138    Average    Average    Average    Average    Average    Average
#>             13         14         15
#> 1         <NA>       <NA>       <NA>
#> 2      Average    Average    Average
#> 3       Active     Active     Active
#> 4      Average     Active    Average
#> 5       Active    Average    Average
#> 6      Average    Average    Average
#> 7      Average    Average    Average
#> 8       Active     Active     Active
#> 9       Active     Active     Active
#> 10     Average    Average    Average
#> 11      Active     Active    Average
#> 12     Average    Average    Average
#> 13      Active     Active     Active
#> 14      Active     Active    Average
#> 15  Disengaged Disengaged    Average
#> 16     Average    Average    Average
#> 17      Active     Active     Active
#> 18     Average Disengaged     Active
#> 19      Active    Average     Active
#> 20  Disengaged     Active    Average
#> 21     Average    Average    Average
#> 22        <NA>       <NA>       <NA>
#> 23     Average    Average    Average
#> 24        <NA>       <NA>       <NA>
#> 25      Active     Active    Average
#> 26      Active    Average     Active
#> 27     Average Disengaged Disengaged
#> 28     Average    Average Disengaged
#> 29      Active    Average     Active
#> 30        <NA>       <NA>       <NA>
#> 31     Average    Average Disengaged
#> 32      Active     Active     Active
#> 33        <NA>       <NA>       <NA>
#> 34     Average    Average Disengaged
#> 35        <NA>       <NA>       <NA>
#> 36     Average Disengaged     Active
#> 37      Active     Active    Average
#> 38     Average    Average    Average
#> 39     Average    Average     Active
#> 40     Average     Active    Average
#> 41      Active Disengaged     Active
#> 42      Active     Active     Active
#> 43      Active     Active    Average
#> 44     Average    Average    Average
#> 45  Disengaged Disengaged     Active
#> 46  Disengaged Disengaged    Average
#> 47     Average    Average    Average
#> 48  Disengaged Disengaged Disengaged
#> 49  Disengaged Disengaged Disengaged
#> 50      Active     Active     Active
#> 51  Disengaged Disengaged Disengaged
#> 52      Active     Active     Active
#> 53  Disengaged Disengaged Disengaged
#> 54     Average    Average    Average
#> 55  Disengaged Disengaged Disengaged
#> 56  Disengaged Disengaged Disengaged
#> 57     Average    Average    Average
#> 58      Active     Active     Active
#> 59     Average    Average    Average
#> 60     Average    Average    Average
#> 61      Active     Active     Active
#> 62     Average    Average    Average
#> 63     Average    Average    Average
#> 64        <NA>       <NA>       <NA>
#> 65     Average    Average    Average
#> 66      Active     Active     Active
#> 67     Average    Average    Average
#> 68      Active     Active     Active
#> 69     Average    Average    Average
#> 70     Average Disengaged       <NA>
#> 71     Average    Average    Average
#> 72     Average    Average    Average
#> 73     Average Disengaged     Active
#> 74     Average Disengaged    Average
#> 75     Average    Average    Average
#> 76      Active Disengaged     Active
#> 77        <NA>       <NA>       <NA>
#> 78     Average    Average     Active
#> 79      Active     Active     Active
#> 80      Active     Active     Active
#> 81      Active     Active     Active
#> 82      Active     Active     Active
#> 83     Average Disengaged    Average
#> 84     Average    Average    Average
#> 85     Average     Active     Active
#> 86      Active     Active     Active
#> 87      Active     Active    Average
#> 88        <NA>       <NA>       <NA>
#> 89      Active     Active    Average
#> 90      Active    Average     Active
#> 91      Active    Average    Average
#> 92      Active    Average    Average
#> 93     Average Disengaged     Active
#> 94      Active     Active     Active
#> 95      Active     Active    Average
#> 96     Average    Average       <NA>
#> 97        <NA>       <NA>       <NA>
#> 98        <NA>       <NA>       <NA>
#> 99        <NA>       <NA>       <NA>
#> 100    Average    Average    Average
#> 101       <NA>       <NA>       <NA>
#> 102       <NA>       <NA>       <NA>
#> 103       <NA>       <NA>       <NA>
#> 104       <NA>       <NA>       <NA>
#> 105       <NA>       <NA>       <NA>
#> 106       <NA>       <NA>       <NA>
#> 107 Disengaged Disengaged Disengaged
#> 108    Average Disengaged     Active
#> 109    Average    Average Disengaged
#> 110    Average Disengaged    Average
#> 111     Active    Average     Active
#> 112 Disengaged Disengaged Disengaged
#> 113 Disengaged Disengaged Disengaged
#> 114    Average Disengaged Disengaged
#> 115     Active     Active     Active
#> 116     Active     Active    Average
#> 117    Average    Average    Average
#> 118     Active    Average     Active
#> 119    Average     Active    Average
#> 120    Average    Average    Average
#> 121    Average    Average Disengaged
#> 122     Active     Active     Active
#> 123    Average    Average    Average
#> 124     Active     Active     Active
#> 125    Average     Active    Average
#> 126    Average     Active    Average
#> 127    Average    Average    Average
#> 128 Disengaged Disengaged Disengaged
#> 129    Average Disengaged Disengaged
#> 130    Average Disengaged     Active
#> 131 Disengaged Disengaged Disengaged
#> 132    Average    Average    Average
#> 133     Active     Active    Average
#> 134    Average    Average    Average
#> 135     Active     Active     Active
#> 136     Active     Active     Active
#> 137     Active Disengaged     Active
#> 138    Average    Average    Average
```
