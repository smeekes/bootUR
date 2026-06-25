# Plot Missing Values

Plots missing values of different types for a time series dataset.

## Usage

``` r
plot_missing_values(y, show_names = FALSE, show_legend = TRUE,
  axis_text_size = NULL, legend_size = NULL, cols = NULL)
```

## Arguments

- y:

  A (\\T\\x\\N\\)-matrix of \\N\\ time series with \\T\\ observations.
  Data may also be in a time series format (e.g. `ts`, `zoo` or `xts`)
  or data frame.

- show_names:

  Show the time series' names on the plot (`TRUE`) or not (`FALSE`).
  Default is `TRUE`.

- show_legend:

  Logical indicator whether a legend should be displayed. Default is
  `TRUE`.

- axis_text_size:

  Size of the text on the axis. Default takes `ggplot2` defaults.

- legend_size:

  Size of the text in the legend if `show_legend = TRUE`. Default takes
  `ggplot2` defaults.

- cols:

  Vector with colours for displaying the different types of data. If the
  default is overwritten, four colours must be supplied.

## Value

A `ggplot2` object containing the missing values plot.

## Details

The function distinguish four types of data: observed data (non-missing)
and three missing types. Type `"Balanced NA"` indicates where entire
rows are missing (`NA`). These do not cause unbalancedness as the
missing rows can simply be deleted. Type `"Unbalanced NA"` are missing
values on the beginning or end of the sample, which cause
unbalancedness. These affect some (but not all) bootstrap methods, see
e.g.~[`boot_fdr`](https://smeekes.github.io/bootUR/reference/boot_fdr.md).
Type `"Internal NA"` are missing values inside the sample, which need to
be removed before the bootstrap unit root tests can be used.

This function requires the package `ggplot2` to be installed. If the
package is not found, plotting is aborted.

## References

Smeekes, S. and Wilms, I. (2023). bootUR: An R Package for Bootstrap
Unit Root Tests. *Journal of Statistical Software*, 106(12), 1-39.
