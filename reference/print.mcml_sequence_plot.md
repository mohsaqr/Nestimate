# Draw a stacked multichannel mcml sequence plot

Print method for the figure
[`sequence_plot`](https://saqr.me/Nestimate/reference/sequence_plot.md)
returns for an `mcml` with more than one channel: one panel per channel
(the macro `Summary` and one per cluster), each with its own legend.

## Usage

``` r
# S3 method for class 'mcml_sequence_plot'
print(x, ...)
```

## Arguments

- x:

  An `mcml_sequence_plot` (a `gtable`).

- ...:

  Ignored.

## Value

`x`, invisibly. Called for the side effect of drawing it on a new page
of the current graphics device.
