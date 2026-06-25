# Plot Orders of Integration

Plots a vector with orders of integration of time series.

## Usage

``` r
plot_order_integration(orders, show_names = TRUE, show_legend = TRUE,
  names_size = NULL, legend_size = NULL, cols = NULL)
```

## Arguments

- orders:

  A `"bootUR", "order_integration"` object obtained from the function
  `order_integration`, or a vector with found orders of integration.

- show_names:

  Show the time series' names on the plot (`TRUE`) or not (`FALSE`).
  Default is `TRUE`.

- show_legend:

  Logical indicator whether a legend should be displayed. Default is
  `TRUE`.

- names_size:

  Size of the time series' names if `show_names = TRUE`. Default takes
  `ggplot2` defaults.

- legend_size:

  Size of the text in the legend if `show_legend = TRUE`. Default takes
  `ggplot2` defaults.

- cols:

  Vector with colours for displaying the orders of integration. At least
  as many colours as orders of integration need to be supplied. Default
  supplies 4 colours for displaying up to \\I(3)\\ series.

## Value

A `ggplot2` object containing the plot of the orders of integration.

## Details

This function requires the package `ggplot2` to be installed. If the
package is not found, plotting is aborted.

## References

Smeekes, S. and Wilms, I. (2023). bootUR: An R Package for Bootstrap
Unit Root Tests. *Journal of Statistical Software*, 106(12), 1-39.

## See also

[`order_integration`](https://smeekes.github.io/bootUR/reference/order_integration.md)
