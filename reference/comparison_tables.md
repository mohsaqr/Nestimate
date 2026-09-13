# Tables of a network comparison

Named accessors for the tables inside a `net_network_comparison` object
(from
[`compare_networks()`](https://saqr.me/Nestimate/reference/compare_networks.md)).
Each returns a plain data frame with one row per unit and no row names.
Numeric columns are rounded to `digits` decimals; p-value columns are
never rounded.

## Usage

``` r
# S3 method for class 'net_network_comparison'
summary(object, pair = NULL, digits = 2L, ...)

edge_differences(x, pair = NULL, digits = 2L)

node_differences(x, pair = NULL, measure = NULL, digits = 2L)

global_differences(x, pair = NULL, digits = 2L)

network_metrics(x, digits = 2L)

# S3 method for class 'net_table'
print(x, digits = attr(x, "digits") %||% 2L, ...)
```

## Arguments

- pair:

  Optional selection of comparisons. Any of: the pair name(s) as printed
  (`"A vs B"`); the two network names (`c("A", "B")`, either order); one
  network name (`"A"`, every pair it takes part in); or index/indices
  into the pair table (`1`, `c(1, 3)`). `NULL` keeps all.

- digits:

  Decimals kept in numeric columns. Default `2`.

- ...:

  Ignored.

- x, object:

  A `net_network_comparison` object (for `print.net_table()`, a table
  returned by one of these verbs).

- measure:

  Optional centrality measure name(s) to keep.

## Value

- [`summary()`](https://rdrr.io/r/base/summary.html): one row per pair –
  `pair`, `network_a`, `network_b`, `n_cells`, `n_differing`,
  `share_higher_a`, `share_higher_b`, `mean_abs_diff`, `max_abs_diff`,
  `pearson`, `spearman`, `cosine`, `jaccard`, `top_edge`,
  `top_edge_higher`, and with inference `n_sig_edges`, `n_sig_nodes`,
  `m_stat`, `m_p`, `s_stat`, `s_p`. The per-edge, per-node, per-metric
  and per-network tables are the four verbs below.

- `edge_differences()`: one row per pair x transition: `pair`,
  `network_a`, `network_b`, `from`, `to`, `weight_a`, `weight_b`,
  `diff`, `abs_diff`, `rel_diff`, `ratio`, `log_ratio`, `rank_a`,
  `rank_b`, `rank_diff`, `percentile_diff`, `status`, `higher`, and the
  inference columns when `test` was used.

- `node_differences()`: one row per pair x state x measure: `pair`,
  `network_a`, `network_b`, `node`, `measure`, `value_a`, `value_b`,
  `diff`, `abs_diff`, `rank_a`, `rank_b`, `higher`, plus inference
  columns. Errors (class `nestimate_compare_no_nodes`) when
  [`compare_networks()`](https://saqr.me/Nestimate/reference/compare_networks.md)
  was called with no `measures`.

- `global_differences()`: one row per pair x metric: `pair`,
  `network_a`, `network_b`, `category`, `metric`, `key`, `value`, plus
  inference rows and columns.

- `network_metrics()`: one row per network x structural metric:
  `network`, `metric`, `value`.

Every table is a data frame of class `net_table` whose
[`print()`](https://rdrr.io/r/base/print.html) shows whole numbers
without decimals, exact zeros as `0`, other values with `digits`
decimals, and p-values with three decimals.

## Examples

``` r
achievers <- build_network(group_regulation_long, method = "relative",
                           actor = "Actor", action = "Action",
                           time = "Time", group = "Achiever")
cmp <- compare_networks(achievers)
summary(cmp)
#>         pair network_a network_b n_cells n_differing share_higher_a
#>  High vs Low      High       Low      81          78           0.51
#>  share_higher_b mean_abs_diff max_abs_diff pearson spearman cosine jaccard
#>            0.46          0.03         0.21    0.92     0.92   0.95    0.75
#>              top_edge top_edge_higher
#>  discuss -> consensus            High
edge_differences(cmp)
#>         pair network_a network_b       from         to weight_a weight_b  diff
#>  High vs Low      High       Low      adapt      adapt        0        0     0
#>  High vs Low      High       Low      adapt   cohesion     0.26     0.28 -0.01
#>  High vs Low      High       Low      adapt  consensus     0.52     0.46  0.06
#>  High vs Low      High       Low      adapt coregulate        0     0.03 -0.03
#>  High vs Low      High       Low      adapt    discuss     0.04     0.07 -0.03
#>  High vs Low      High       Low      adapt    emotion     0.14     0.11  0.03
#>  High vs Low      High       Low      adapt    monitor     0.03     0.04 -0.01
#>  High vs Low      High       Low      adapt       plan     0.01     0.02     0
#>  High vs Low      High       Low      adapt  synthesis        0        0     0
#>  High vs Low      High       Low   cohesion      adapt     0.01        0  0.01
#>  High vs Low      High       Low   cohesion   cohesion     0.04     0.01  0.04
#>  High vs Low      High       Low   cohesion  consensus     0.54     0.45  0.09
#>  High vs Low      High       Low   cohesion coregulate     0.08     0.17 -0.09
#>  High vs Low      High       Low   cohesion    discuss     0.04     0.08 -0.04
#>  High vs Low      High       Low   cohesion    emotion     0.12     0.11  0.01
#>  High vs Low      High       Low   cohesion    monitor     0.02     0.05 -0.04
#>  High vs Low      High       Low   cohesion       plan     0.15     0.13  0.02
#>  High vs Low      High       Low   cohesion  synthesis     0.01        0  0.01
#>  High vs Low      High       Low  consensus      adapt        0     0.01     0
#>  High vs Low      High       Low  consensus   cohesion     0.02     0.01  0.01
#>  High vs Low      High       Low  consensus  consensus     0.08     0.08     0
#>  High vs Low      High       Low  consensus coregulate     0.17     0.21 -0.04
#>  High vs Low      High       Low  consensus    discuss     0.23     0.14  0.10
#>  High vs Low      High       Low  consensus    emotion     0.08     0.06  0.02
#>  High vs Low      High       Low  consensus    monitor     0.04     0.06 -0.02
#>  High vs Low      High       Low  consensus       plan     0.36     0.43 -0.07
#>  High vs Low      High       Low  consensus  synthesis     0.01     0.01     0
#>  High vs Low      High       Low coregulate      adapt     0.02     0.01  0.01
#>  High vs Low      High       Low coregulate   cohesion     0.04     0.04     0
#>  High vs Low      High       Low coregulate  consensus     0.11     0.16 -0.05
#>  High vs Low      High       Low coregulate coregulate     0.01     0.03 -0.02
#>  High vs Low      High       Low coregulate    discuss     0.23     0.31 -0.07
#>  High vs Low      High       Low coregulate    emotion     0.20     0.15  0.06
#>  High vs Low      High       Low coregulate    monitor     0.10     0.08  0.02
#>  High vs Low      High       Low coregulate       plan     0.27     0.22  0.05
#>  High vs Low      High       Low coregulate  synthesis     0.02     0.02     0
#>  High vs Low      High       Low    discuss      adapt     0.02     0.12 -0.10
#>  High vs Low      High       Low    discuss   cohesion     0.06     0.03  0.03
#>  High vs Low      High       Low    discuss  consensus     0.42     0.21  0.21
#>  High vs Low      High       Low    discuss coregulate     0.07     0.10 -0.02
#>  High vs Low      High       Low    discuss    discuss     0.17     0.22 -0.05
#>  High vs Low      High       Low    discuss    emotion     0.11     0.10  0.01
#>  High vs Low      High       Low    discuss    monitor     0.02     0.03 -0.01
#>  High vs Low      High       Low    discuss       plan     0.01     0.01     0
#>  High vs Low      High       Low    discuss  synthesis     0.11     0.18 -0.07
#>  High vs Low      High       Low    emotion      adapt        0        0     0
#>  High vs Low      High       Low    emotion   cohesion     0.33     0.32     0
#>  High vs Low      High       Low    emotion  consensus     0.34     0.30  0.03
#>  High vs Low      High       Low    emotion coregulate     0.02     0.05 -0.02
#>  High vs Low      High       Low    emotion    discuss     0.12     0.08  0.04
#>  High vs Low      High       Low    emotion    emotion     0.06     0.09 -0.03
#>  High vs Low      High       Low    emotion    monitor     0.03     0.04 -0.01
#>  High vs Low      High       Low    emotion       plan     0.09     0.11 -0.02
#>  High vs Low      High       Low    emotion  synthesis     0.01        0  0.01
#>  High vs Low      High       Low    monitor      adapt     0.01     0.01     0
#>  High vs Low      High       Low    monitor   cohesion     0.05     0.06 -0.01
#>  High vs Low      High       Low    monitor  consensus     0.16     0.16     0
#>  High vs Low      High       Low    monitor coregulate     0.05     0.06 -0.01
#>  High vs Low      High       Low    monitor    discuss     0.37     0.38 -0.01
#>  High vs Low      High       Low    monitor    emotion     0.10     0.09  0.01
#>  High vs Low      High       Low    monitor    monitor     0.02     0.02     0
#>  High vs Low      High       Low    monitor       plan     0.23     0.21  0.02
#>  High vs Low      High       Low    monitor  synthesis     0.02     0.01  0.01
#>  High vs Low      High       Low       plan      adapt        0        0     0
#>  High vs Low      High       Low       plan   cohesion     0.03     0.02  0.01
#>  High vs Low      High       Low       plan  consensus     0.29     0.29  0.01
#>  High vs Low      High       Low       plan coregulate     0.02     0.01  0.01
#>  High vs Low      High       Low       plan    discuss     0.06     0.07 -0.01
#>  High vs Low      High       Low       plan    emotion     0.18     0.12  0.07
#>  High vs Low      High       Low       plan    monitor     0.08     0.08     0
#>  High vs Low      High       Low       plan       plan     0.33     0.42 -0.09
#>  High vs Low      High       Low       plan  synthesis        0        0     0
#>  High vs Low      High       Low  synthesis      adapt     0.14     0.30 -0.16
#>  High vs Low      High       Low  synthesis   cohesion     0.03     0.04 -0.01
#>  High vs Low      High       Low  synthesis  consensus     0.58     0.39  0.19
#>  High vs Low      High       Low  synthesis coregulate     0.01     0.07 -0.05
#>  High vs Low      High       Low  synthesis    discuss     0.03     0.09 -0.06
#>  High vs Low      High       Low  synthesis    emotion     0.06     0.07 -0.01
#>  High vs Low      High       Low  synthesis    monitor        0     0.02 -0.02
#>  High vs Low      High       Low  synthesis       plan     0.14     0.02  0.12
#>  High vs Low      High       Low  synthesis  synthesis        0        0     0
#>  abs_diff rel_diff ratio log_ratio rank_a rank_b rank_diff percentile_diff
#>         0       NA    NA         0      3   3.50     -0.50           -0.01
#>      0.01     0.03  0.95     -0.01     70     70         0               0
#>      0.06     0.06  1.12      0.04     79     81        -2           -0.02
#>      0.03        1     0     -0.03      3     26       -23           -0.26
#>      0.03     0.31  0.52     -0.03     35     40        -5           -0.06
#>      0.03     0.12  1.27      0.03     58     54         4            0.05
#>      0.01     0.11  0.80     -0.01     29     29         0               0
#>         0     0.07  0.87         0     17     19        -2           -0.02
#>         0       NA    NA         0      3   3.50     -0.50           -0.01
#>      0.01        1    NA      0.01     11   3.50      7.50            0.06
#>      0.04     0.74  6.62      0.04     38     11        27            0.33
#>      0.09     0.09  1.19      0.06     80     80         0               0
#>      0.09     0.35  0.49     -0.08     47     63       -16           -0.20
#>      0.04     0.35  0.49     -0.04     37     47       -10           -0.12
#>      0.01     0.03  1.05      0.01     56     55         1            0.01
#>      0.04     0.51  0.32     -0.03     20     34       -14           -0.17
#>      0.02     0.08  1.18      0.02     61     58         3            0.04
#>      0.01        1    NA      0.01     12   3.50      8.50            0.07
#>         0     0.14  0.76         0      9     10        -1           -0.01
#>      0.01     0.36  2.15      0.01     24     13        11            0.14
#>         0     0.02  1.04         0     49     46         3            0.04
#>      0.04     0.10  0.83     -0.03     64     65        -1           -0.01
#>      0.10     0.26  1.70      0.08     68     59         9            0.11
#>      0.02     0.13  1.30      0.02     48     37        11            0.14
#>      0.02     0.25  0.59     -0.02     34     35        -1           -0.01
#>      0.07     0.08  0.84     -0.05     76     79        -3           -0.04
#>         0     0.05  1.11         0     13     12         1            0.01
#>      0.01     0.33  2.01      0.01     25     15        10            0.12
#>         0     0.01  0.99         0     36     30         6            0.07
#>      0.05     0.18  0.69     -0.04     54     61        -7           -0.09
#>      0.02     0.40  0.42     -0.02     16     27       -11           -0.14
#>      0.07     0.13  0.77     -0.06     69     74        -5           -0.06
#>      0.06     0.17  1.40      0.05     66     60         6            0.07
#>      0.02     0.10  1.23      0.02     51     45         6            0.07
#>      0.05     0.10  1.23      0.04     71     68         3            0.04
#>         0     0.01  1.02         0     23     21         2            0.02
#>      0.10     0.67  0.20     -0.09     28     57       -29           -0.36
#>      0.03     0.31  1.88      0.03     42     28        14            0.17
#>      0.21     0.33  1.98      0.16     78     67        11            0.14
#>      0.02     0.14  0.75     -0.02     45     51        -6           -0.07
#>      0.05     0.13  0.76     -0.04     63     69        -6           -0.07
#>      0.01     0.06  1.13      0.01     55     52         3            0.04
#>      0.01     0.26  0.58     -0.01     19     25        -6           -0.07
#>         0     0.07  1.16         0     15     14         1            0.01
#>      0.07     0.25  0.60     -0.06     53     64       -11           -0.14
#>         0     0.35  2.08         0      7      9        -2           -0.02
#>         0        0     1         0     73     75        -2           -0.02
#>      0.03     0.05  1.11      0.03     75     72         3            0.04
#>      0.02     0.34  0.49     -0.02     26     33        -7           -0.09
#>      0.04     0.22  1.57      0.04     57     44        13            0.16
#>      0.03     0.20  0.67     -0.03     43     50        -7           -0.09
#>      0.01     0.14  0.75     -0.01     33     32         1            0.01
#>      0.02     0.10  0.81     -0.02     50     53        -3           -0.04
#>      0.01        1    NA      0.01     10   3.50      6.50            0.05
#>         0     0.01  0.98         0     14     16        -2           -0.02
#>      0.01     0.11  0.80     -0.01     39     36         3            0.04
#>         0        0  1.01         0     62     62         0               0
#>      0.01     0.12  0.79     -0.01     40     38         2            0.02
#>      0.01     0.01  0.97     -0.01     77     76         1            0.01
#>      0.01     0.06  1.12      0.01     52     48         4            0.05
#>         0     0.04  1.08         0  21.50     20      1.50            0.02
#>      0.02     0.04  1.09      0.02     67     66         1            0.01
#>      0.01     0.16  1.38      0.01  21.50     18      3.50            0.05
#>         0     0.39  2.26         0      6      8        -2           -0.02
#>      0.01     0.23  1.61      0.01     32     22        10            0.12
#>      0.01     0.01  1.02      0.01     72     71         1            0.01
#>      0.01     0.36  2.11      0.01     27     17        10            0.12
#>      0.01     0.11  0.81     -0.01     41     41         0               0
#>      0.07     0.22  1.58      0.06     65     56         9            0.11
#>         0        0  1.01         0     46     43         3            0.04
#>      0.09     0.12  0.79     -0.06     74     78        -4           -0.05
#>         0     0.84 11.29         0      8      7         1            0.01
#>      0.16     0.35  0.48     -0.13  59.50     73    -13.50           -0.16
#>      0.01     0.13  0.77     -0.01  30.50     31     -0.50               0
#>      0.19     0.20  1.49      0.13     81     77         4            0.05
#>      0.05     0.65  0.22     -0.05     18     39       -21           -0.26
#>      0.06     0.51  0.33     -0.06  30.50     49    -18.50           -0.22
#>      0.01     0.07  0.86     -0.01     44     42         2            0.02
#>      0.02        1     0     -0.02      3     23       -20           -0.22
#>      0.12     0.71  5.98      0.11  59.50     24     35.50            0.44
#>         0       NA    NA         0      3   3.50     -0.50           -0.01
#>   status higher
#>  neither  equal
#>     both    Low
#>     both   High
#>   only_b    Low
#>     both    Low
#>     both   High
#>     both    Low
#>     both    Low
#>  neither  equal
#>   only_a   High
#>     both   High
#>     both   High
#>     both    Low
#>     both    Low
#>     both   High
#>     both    Low
#>     both   High
#>   only_a   High
#>     both    Low
#>     both   High
#>     both   High
#>     both    Low
#>     both   High
#>     both   High
#>     both    Low
#>     both    Low
#>     both   High
#>     both   High
#>     both    Low
#>     both    Low
#>     both    Low
#>     both    Low
#>     both   High
#>     both   High
#>     both   High
#>     both   High
#>     both    Low
#>     both   High
#>     both   High
#>     both    Low
#>     both    Low
#>     both   High
#>     both    Low
#>     both   High
#>     both    Low
#>     both   High
#>     both   High
#>     both   High
#>     both    Low
#>     both   High
#>     both    Low
#>     both    Low
#>     both    Low
#>   only_a   High
#>     both    Low
#>     both    Low
#>     both   High
#>     both    Low
#>     both    Low
#>     both   High
#>     both   High
#>     both   High
#>     both   High
#>     both   High
#>     both   High
#>     both   High
#>     both   High
#>     both    Low
#>     both   High
#>     both   High
#>     both    Low
#>     both   High
#>     both    Low
#>     both    Low
#>     both   High
#>     both    Low
#>     both    Low
#>     both    Low
#>   only_b    Low
#>     both   High
#>  neither  equal
edge_differences(cmp, digits = 4)
#>         pair network_a network_b       from         to weight_a weight_b
#>  High vs Low      High       Low      adapt      adapt        0        0
#>  High vs Low      High       Low      adapt   cohesion   0.2624   0.2772
#>  High vs Low      High       Low      adapt  consensus   0.5177   0.4620
#>  High vs Low      High       Low      adapt coregulate        0   0.0299
#>  High vs Low      High       Low      adapt    discuss   0.0355   0.0679
#>  High vs Low      High       Low      adapt    emotion   0.1418   0.1114
#>  High vs Low      High       Low      adapt    monitor   0.0284   0.0353
#>  High vs Low      High       Low      adapt       plan   0.0142   0.0163
#>  High vs Low      High       Low      adapt  synthesis        0        0
#>  High vs Low      High       Low   cohesion      adapt   0.0053        0
#>  High vs Low      High       Low   cohesion   cohesion   0.0437   0.0066
#>  High vs Low      High       Low   cohesion  consensus   0.5362   0.4505
#>  High vs Low      High       Low   cohesion coregulate   0.0810   0.1664
#>  High vs Low      High       Low   cohesion    discuss   0.0405   0.0832
#>  High vs Low      High       Low   cohesion    emotion   0.1183   0.1123
#>  High vs Low      High       Low   cohesion    monitor   0.0171   0.0528
#>  High vs Low      High       Low   cohesion       plan   0.1514   0.1281
#>  High vs Low      High       Low   cohesion  synthesis   0.0064        0
#>  High vs Low      High       Low  consensus      adapt   0.0041   0.0054
#>  High vs Low      High       Low  consensus   cohesion   0.0198   0.0092
#>  High vs Low      High       Low  consensus  consensus   0.0834   0.0804
#>  High vs Low      High       Low  consensus coregulate   0.1710   0.2070
#>  High vs Low      High       Low  consensus    discuss   0.2326   0.1365
#>  High vs Low      High       Low  consensus    emotion   0.0814   0.0626
#>  High vs Low      High       Low  consensus    monitor   0.0354   0.0596
#>  High vs Low      High       Low  consensus       plan   0.3644   0.4321
#>  High vs Low      High       Low  consensus  synthesis   0.0080   0.0072
#>  High vs Low      High       Low coregulate      adapt   0.0224   0.0112
#>  High vs Low      High       Low coregulate   cohesion   0.0358   0.0362
#>  High vs Low      High       Low coregulate  consensus   0.1085   0.1561
#>  High vs Low      High       Low coregulate coregulate   0.0134   0.0316
#>  High vs Low      High       Low coregulate    discuss   0.2349   0.3058
#>  High vs Low      High       Low coregulate    emotion   0.2036   0.1459
#>  High vs Low      High       Low coregulate    monitor   0.0962   0.0781
#>  High vs Low      High       Low coregulate       plan   0.2662   0.2165
#>  High vs Low      High       Low coregulate  synthesis   0.0190   0.0186
#>  High vs Low      High       Low    discuss      adapt   0.0240   0.1201
#>  High vs Low      High       Low    discuss   cohesion   0.0619   0.0329
#>  High vs Low      High       Low    discuss  consensus   0.4249   0.2146
#>  High vs Low      High       Low    discuss coregulate   0.0724   0.0965
#>  High vs Low      High       Low    discuss    discuss   0.1692   0.2213
#>  High vs Low      High       Low    discuss    emotion   0.1123   0.0991
#>  High vs Low      High       Low    discuss    monitor   0.0165   0.0282
#>  High vs Low      High       Low    discuss       plan   0.0125   0.0108
#>  High vs Low      High       Low    discuss  synthesis   0.1063   0.1766
#>  High vs Low      High       Low    emotion      adapt   0.0032   0.0016
#>  High vs Low      High       Low    emotion   cohesion   0.3258   0.3248
#>  High vs Low      High       Low    emotion  consensus   0.3361   0.3015
#>  High vs Low      High       Low    emotion coregulate   0.0232   0.0474
#>  High vs Low      High       Low    emotion    discuss   0.1219   0.0777
#>  High vs Low      High       Low    emotion    emotion   0.0626   0.0940
#>  High vs Low      High       Low    emotion    monitor   0.0316   0.0420
#>  High vs Low      High       Low    emotion       plan   0.0903   0.1111
#>  High vs Low      High       Low    emotion  synthesis   0.0052        0
#>  High vs Low      High       Low    monitor      adapt   0.0111   0.0112
#>  High vs Low      High       Low    monitor   cohesion   0.0490   0.0612
#>  High vs Low      High       Low    monitor  consensus   0.1596   0.1588
#>  High vs Low      High       Low    monitor coregulate   0.0506   0.0638
#>  High vs Low      High       Low    monitor    discuss   0.3697   0.3800
#>  High vs Low      High       Low    monitor    emotion   0.0964   0.0862
#>  High vs Low      High       Low    monitor    monitor   0.0190   0.0175
#>  High vs Low      High       Low    monitor       plan   0.2259   0.2075
#>  High vs Low      High       Low    monitor  synthesis   0.0190   0.0138
#>  High vs Low      High       Low       plan      adapt   0.0014   0.0006
#>  High vs Low      High       Low       plan   cohesion   0.0315   0.0196
#>  High vs Low      High       Low       plan  consensus   0.2939   0.2873
#>  High vs Low      High       Low       plan coregulate   0.0239   0.0113
#>  High vs Low      High       Low       plan    discuss   0.0602   0.0747
#>  High vs Low      High       Low       plan    emotion   0.1822   0.1155
#>  High vs Low      High       Low       plan    monitor   0.0757   0.0753
#>  High vs Low      High       Low       plan       plan   0.3278   0.4153
#>  High vs Low      High       Low       plan  synthesis   0.0035   0.0003
#>  High vs Low      High       Low  synthesis      adapt   0.1439   0.3021
#>  High vs Low      High       Low  synthesis   cohesion   0.0288   0.0374
#>  High vs Low      High       Low  synthesis  consensus   0.5755   0.3850
#>  High vs Low      High       Low  synthesis coregulate   0.0144   0.0668
#>  High vs Low      High       Low  synthesis    discuss   0.0288   0.0882
#>  High vs Low      High       Low  synthesis    emotion   0.0647   0.0749
#>  High vs Low      High       Low  synthesis    monitor        0   0.0214
#>  High vs Low      High       Low  synthesis       plan   0.1439   0.0241
#>  High vs Low      High       Low  synthesis  synthesis        0        0
#>     diff abs_diff rel_diff   ratio log_ratio  rank_a rank_b rank_diff
#>        0        0       NA      NA         0       3 3.5000   -0.5000
#>  -0.0148   0.0148   0.0274  0.9467   -0.0116      70     70         0
#>   0.0558   0.0558   0.0569  1.1207    0.0374      79     81        -2
#>  -0.0299   0.0299        1       0   -0.0295       3     26       -23
#>  -0.0325   0.0325   0.3141  0.5220   -0.0309      35     40        -5
#>   0.0304   0.0304   0.1202  1.2731    0.0270      58     54         4
#>  -0.0070   0.0070   0.1092  0.8031   -0.0067      29     29         0
#>  -0.0021   0.0021   0.0695  0.8700   -0.0021      17     19        -2
#>        0        0       NA      NA         0       3 3.5000   -0.5000
#>   0.0053   0.0053        1      NA    0.0053      11 3.5000    7.5000
#>   0.0371   0.0371   0.7375  6.6177    0.0362      38     11        27
#>   0.0858   0.0858   0.0869  1.1904    0.0575      80     80         0
#>  -0.0854   0.0854   0.3452  0.4868   -0.0761      47     63       -16
#>  -0.0427   0.0427   0.3452  0.4868   -0.0402      37     47       -10
#>   0.0061   0.0061   0.0262  1.0539    0.0054      56     55         1
#>  -0.0358   0.0358   0.5119  0.3228   -0.0346      20     34       -14
#>   0.0232   0.0232   0.0832  1.1814    0.0204      61     58         3
#>   0.0064   0.0064        1      NA    0.0064      12 3.5000    8.5000
#>  -0.0013   0.0013   0.1379  0.7576   -0.0013       9     10        -1
#>   0.0106   0.0106   0.3648  2.1486    0.0104      24     13        11
#>   0.0031   0.0031   0.0188  1.0383    0.0028      49     46         3
#>  -0.0360   0.0360   0.0953  0.8260   -0.0303      64     65        -1
#>   0.0961   0.0961   0.2603  1.7037    0.0811      68     59         9
#>   0.0187   0.0187   0.1300  1.2988    0.0175      48     37        11
#>  -0.0242   0.0242   0.2549  0.5937   -0.0231      34     35        -1
#>  -0.0677   0.0677   0.0850  0.8433   -0.0484      76     79        -3
#>   0.0008   0.0008   0.0536  1.1132    0.0008      13     12         1
#>   0.0112   0.0112   0.3347  2.0060    0.0110      25     15        10
#>  -0.0005   0.0005   0.0063  0.9876   -0.0004      36     30         6
#>  -0.0476   0.0476   0.1800  0.6949   -0.0421      54     61        -7
#>  -0.0182   0.0182   0.4037  0.4248   -0.0178      16     27       -11
#>  -0.0709   0.0709   0.1311  0.7682   -0.0558      69     74        -5
#>   0.0577   0.0577   0.1650  1.3952    0.0491      66     60         6
#>   0.0181   0.0181   0.1040  1.2322    0.0167      51     45         6
#>   0.0497   0.0497   0.1029  1.2294    0.0400      71     68         3
#>   0.0004   0.0004   0.0114  1.0230    0.0004      23     21         2
#>  -0.0962   0.0962   0.6674  0.1995   -0.0898      28     57       -29
#>   0.0291   0.0291   0.3066  1.8843    0.0277      42     28        14
#>   0.2103   0.2103   0.3289  1.9800    0.1597      78     67        11
#>  -0.0241   0.0241   0.1428  0.7501   -0.0222      45     51        -6
#>  -0.0520   0.0520   0.1332  0.7649   -0.0435      63     69        -6
#>   0.0133   0.0133   0.0627  1.1338    0.0120      55     52         3
#>  -0.0118   0.0118   0.2630  0.5835   -0.0115      19     25        -6
#>   0.0017   0.0017   0.0731  1.1578    0.0017      15     14         1
#>  -0.0703   0.0703   0.2483  0.6022   -0.0616      53     64       -11
#>   0.0017   0.0017   0.3498  2.0758    0.0017       7      9        -2
#>   0.0010   0.0010   0.0016  1.0031    0.0008      73     75        -2
#>   0.0347   0.0347   0.0543  1.1149    0.0263      75     72         3
#>  -0.0242   0.0242   0.3423  0.4900   -0.0233      26     33        -7
#>   0.0442   0.0442   0.2216  1.5693    0.0402      57     44        13
#>  -0.0314   0.0314   0.2007  0.6656   -0.0292      43     50        -7
#>  -0.0103   0.0103   0.1406  0.7534   -0.0100      33     32         1
#>  -0.0208   0.0208   0.1032  0.8129   -0.0189      50     53        -3
#>   0.0052   0.0052        1      NA    0.0051      10 3.5000    6.5000
#>  -0.0002   0.0002   0.0086  0.9830   -0.0002      14     16        -2
#>  -0.0123   0.0123   0.1114  0.7996   -0.0116      39     36         3
#>   0.0008   0.0008   0.0025  1.0051    0.0007      62     62         0
#>  -0.0132   0.0132   0.1155  0.7930   -0.0125      40     38         2
#>  -0.0103   0.0103   0.0138  0.9728   -0.0075      77     76         1
#>   0.0101   0.0101   0.0554  1.1173    0.0093      52     48         4
#>   0.0015   0.0015   0.0400  1.0833    0.0014 21.5000     20    1.5000
#>   0.0184   0.0184   0.0425  1.0887    0.0151      67     66         1
#>   0.0052   0.0052   0.1592  1.3787    0.0051 21.5000     18    3.5000
#>   0.0008   0.0008   0.3861  2.2580    0.0008       6      8        -2
#>   0.0119   0.0119   0.2323  1.6053    0.0116      32     22        10
#>   0.0066   0.0066   0.0114  1.0231    0.0051      72     71         1
#>   0.0125   0.0125   0.3560  2.1054    0.0123      27     17        10
#>  -0.0146   0.0146   0.1080  0.8051   -0.0136      41     41         0
#>   0.0668   0.0668   0.2243  1.5782    0.0581      65     56         9
#>   0.0004   0.0004   0.0025  1.0051    0.0004      46     43         3
#>  -0.0875   0.0875   0.1178  0.7893   -0.0638      74     78        -4
#>   0.0032   0.0032   0.8373 11.2898    0.0031       8      7         1
#>  -0.1583   0.1583   0.3548  0.4762   -0.1296 59.5000     73  -13.5000
#>  -0.0087   0.0087   0.1307  0.7688   -0.0084 30.5000     31   -0.5000
#>   0.1905   0.1905   0.1983  1.4948    0.1289      81     77         4
#>  -0.0525   0.0525   0.6457  0.2153   -0.0504      18     39       -21
#>  -0.0595   0.0595   0.5081  0.3261   -0.0562 30.5000     49  -18.5000
#>  -0.0101   0.0101   0.0725  0.8649   -0.0095      44     42         2
#>  -0.0214   0.0214        1       0   -0.0212       3     23       -20
#>   0.1198   0.1198   0.7134  5.9792    0.1107 59.5000     24   35.5000
#>        0        0       NA      NA         0       3 3.5000   -0.5000
#>  percentile_diff  status higher
#>          -0.0123 neither  equal
#>                0    both    Low
#>          -0.0247    both   High
#>          -0.2593  only_b    Low
#>          -0.0617    both    Low
#>           0.0494    both   High
#>                0    both    Low
#>          -0.0247    both    Low
#>          -0.0123 neither  equal
#>           0.0617  only_a   High
#>           0.3333    both   High
#>                0    both   High
#>          -0.1975    both    Low
#>          -0.1235    both    Low
#>           0.0123    both   High
#>          -0.1728    both    Low
#>           0.0370    both   High
#>           0.0741  only_a   High
#>          -0.0123    both    Low
#>           0.1358    both   High
#>           0.0370    both   High
#>          -0.0123    both    Low
#>           0.1111    both   High
#>           0.1358    both   High
#>          -0.0123    both    Low
#>          -0.0370    both    Low
#>           0.0123    both   High
#>           0.1235    both   High
#>           0.0741    both    Low
#>          -0.0864    both    Low
#>          -0.1358    both    Low
#>          -0.0617    both    Low
#>           0.0741    both   High
#>           0.0741    both   High
#>           0.0370    both   High
#>           0.0247    both   High
#>          -0.3580    both    Low
#>           0.1728    both   High
#>           0.1358    both   High
#>          -0.0741    both    Low
#>          -0.0741    both    Low
#>           0.0370    both   High
#>          -0.0741    both    Low
#>           0.0123    both   High
#>          -0.1358    both    Low
#>          -0.0247    both   High
#>          -0.0247    both   High
#>           0.0370    both   High
#>          -0.0864    both    Low
#>           0.1605    both   High
#>          -0.0864    both    Low
#>           0.0123    both    Low
#>          -0.0370    both    Low
#>           0.0494  only_a   High
#>          -0.0247    both    Low
#>           0.0370    both    Low
#>                0    both   High
#>           0.0247    both    Low
#>           0.0123    both    Low
#>           0.0494    both   High
#>           0.0247    both   High
#>           0.0123    both   High
#>           0.0494    both   High
#>          -0.0247    both   High
#>           0.1235    both   High
#>           0.0123    both   High
#>           0.1235    both   High
#>                0    both    Low
#>           0.1111    both   High
#>           0.0370    both   High
#>          -0.0494    both    Low
#>           0.0123    both   High
#>          -0.1605    both    Low
#>                0    both    Low
#>           0.0494    both   High
#>          -0.2593    both    Low
#>          -0.2222    both    Low
#>           0.0247    both    Low
#>          -0.2222  only_b    Low
#>           0.4444    both   High
#>          -0.0123 neither  equal
node_differences(cmp, measure = "InStrength")
#>         pair network_a network_b       node    measure value_a value_b  diff
#>  High vs Low      High       Low      adapt InStrength    0.22    0.45 -0.24
#>  High vs Low      High       Low   cohesion InStrength    0.81    0.80  0.02
#>  High vs Low      High       Low  consensus InStrength    2.95    2.42  0.54
#>  High vs Low      High       Low coregulate InStrength    0.44    0.69 -0.25
#>  High vs Low      High       Low    discuss InStrength    1.12    1.21 -0.09
#>  High vs Low      High       Low    emotion InStrength       1    0.81  0.19
#>  High vs Low      High       Low    monitor InStrength    0.30    0.39 -0.09
#>  High vs Low      High       Low       plan InStrength    1.27    1.15  0.12
#>  High vs Low      High       Low  synthesis InStrength    0.17    0.22 -0.05
#>  abs_diff rank_a rank_b higher
#>      0.24      2      3    Low
#>      0.02      5      5   High
#>      0.54      9      9   High
#>      0.25      4      4    Low
#>      0.09      7      8    Low
#>      0.19      6      6   High
#>      0.09      3      2    Low
#>      0.12      8      7   High
#>      0.05      1      1    Low
global_differences(cmp)
#>         pair network_a network_b             category               metric
#>  High vs Low      High       Low    Weight Deviations      Mean Abs. Diff.
#>  High vs Low      High       Low    Weight Deviations    Median Abs. Diff.
#>  High vs Low      High       Low    Weight Deviations            RMS Diff.
#>  High vs Low      High       Low    Weight Deviations       Max Abs. Diff.
#>  High vs Low      High       Low    Weight Deviations Rel. Mean Abs. Diff.
#>  High vs Low      High       Low    Weight Deviations             CV Ratio
#>  High vs Low      High       Low         Correlations              Pearson
#>  High vs Low      High       Low         Correlations             Spearman
#>  High vs Low      High       Low         Correlations              Kendall
#>  High vs Low      High       Low         Correlations             Distance
#>  High vs Low      High       Low      Dissimilarities            Euclidean
#>  High vs Low      High       Low      Dissimilarities            Manhattan
#>  High vs Low      High       Low      Dissimilarities             Canberra
#>  High vs Low      High       Low      Dissimilarities          Bray-Curtis
#>  High vs Low      High       Low      Dissimilarities            Frobenius
#>  High vs Low      High       Low         Similarities               Cosine
#>  High vs Low      High       Low         Similarities              Jaccard
#>  High vs Low      High       Low         Similarities                 Dice
#>  High vs Low      High       Low         Similarities              Overlap
#>  High vs Low      High       Low         Similarities                   RV
#>  High vs Low      High       Low Pattern Similarities       Rank Agreement
#>  High vs Low      High       Low Pattern Similarities       Sign Agreement
#>              key value
#>    mean_abs_diff  0.03
#>  median_abs_diff  0.02
#>         rms_diff  0.05
#>     max_abs_diff  0.21
#>     rel_mean_abs  0.29
#>         cv_ratio  1.10
#>          pearson  0.92
#>         spearman  0.92
#>          kendall  0.77
#>     distance_cor  0.84
#>        euclidean  0.47
#>        manhattan  2.61
#>         canberra 14.76
#>      bray_curtis  0.15
#>        frobenius  0.22
#>           cosine  0.95
#>          jaccard  0.75
#>             dice  0.85
#>          overlap  0.85
#>               rv  0.90
#>   rank_agreement  0.82
#>   sign_agreement  0.94
network_metrics(cmp)
#>  network                      metric value
#>     High                  Node Count     9
#>     High                  Edge Count    76
#>     High             Network Density     1
#>     High               Mean Distance  0.04
#>     High           Mean Out-Strength     1
#>     High             SD Out-Strength  0.91
#>     High            Mean In-Strength     1
#>     High              SD In-Strength     0
#>     High             Mean Out-Degree  8.44
#>     High               SD Out-Degree  1.13
#>     High Centralization (Out-Degree)  0.05
#>     High  Centralization (In-Degree)  0.05
#>     High                 Reciprocity  0.96
#>      Low                  Node Count     9
#>      Low                  Edge Count    75
#>      Low             Network Density     1
#>      Low               Mean Distance  0.06
#>      Low           Mean Out-Strength     1
#>      Low             SD Out-Strength  0.72
#>      Low            Mean In-Strength     1
#>      Low              SD In-Strength     0
#>      Low             Mean Out-Degree  8.33
#>      Low               SD Out-Degree  0.87
#>      Low Centralization (Out-Degree)  0.06
#>      Low  Centralization (In-Degree)  0.06
#>      Low                 Reciprocity  0.94
```
