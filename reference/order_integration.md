# Determine Order of Integration

Determines the order of integration for each time series in a dataset
via a sequence of unit root tests, and differences the data accordingly
to eliminate stochastic trends.

## Usage

``` r
order_integration(data, max_order = 2, method = "boot_ur", level = 0.05,
  plot_orders = FALSE, data_name = NULL, ...)
```

## Arguments

- data:

  A (\\T\\x\\N\\)-matrix of \\N\\ time series with \\T\\ observations.
  Data may also be in a time series format (e.g. `ts`, `zoo` or `xts`)
  or data frame.

- max_order:

  The maximum order of integration of the time series. Default is 2.

- method:

  The unit root tests to be used in the procedure. For multiple time
  series the options are "boot_ur", "boot_sqt" and "boot_fdr", with
  "boot_ur" the default. For single time series the options are "adf",
  boot_adf", "boot_union" and "boot_ur", with the latter the default.

- level:

  Desired significance level of the unit root test. Default is 0.05.

- plot_orders:

  Logical indicator whether the resulting orders of integration should
  be plotted. Default is `FALSE`.

- data_name:

  Optional name for the data, to be used in the output. The default uses
  the name of the 'data' argument.

- ...:

  Optional arguments passed to the chosen unit root test function.

## Value

An object of class `"bootUR", "order_integration"` with the following
components

- `order_int`:

  A vector with the found orders of integration of each time series.

- `diff_data`:

  The appropriately differenced data according to `order_int` in the
  same format as the original data.

## Details

The function follows the approach laid out in Smeekes and Wijler (2020),
where all series is differenced \\d-1\\ times, where \\d\\ is the
specified maximum order, and these differenced series are tested for
unit roots. The series for which the unit root null is not rejected, are
classified as \\I(d)\\ and removed from consideration. The remaining
series are integrated, and tested for unit roots again, leading to a
classification of \\I(d-1)\\ series if the null is not rejected. This is
continued until a non-rejection is observed for all time series, or the
series are integrated back to their original level. The series for which
the null hypothesis is rejected in the final stage are classified as
\\I(0)\\.

Care must be taken when using
[`boot_sqt`](https://smeekes.github.io/bootUR/reference/boot_sqt.md)
when the argument `steps` is given as a sequence of integers. As at each
step series are removed, one may end up with fewer series to test than
indicated in `steps`. While integers larger than the number of series
will automatically be removed - along with a warning - by the test, it
is recommend to set `steps` in the form of quantiles.

Plotting the orders of integration requires the `ggplot2` package to be
installed; plot will be skipped and a warning is given if not. For plots
the function
[`plot_order_integration`](https://smeekes.github.io/bootUR/reference/plot_order_integration.md)
is called. The user may prefer to set `plot_orders = FALSE` and call
this function directly using the returned value of `order_int` in order
to have more control over plot settings and save the plot object.

## References

Smeekes, S. and Wilms, I. (2023). bootUR: An R Package for Bootstrap
Unit Root Tests. *Journal of Statistical Software*, 106(12), 1-39.

Smeekes, S. and Wijler, E. (2020). Unit roots and cointegration. In P.
Fuleky (Ed.) *Macroeconomic Forecasting in the Era of Big Data*, Chapter
17, pp. 541-584. *Advanced Studies in Theoretical and Applied
Econometrics*, vol. 52. Springer.

## Examples

``` r
# Use "boot_ur" to determine the order of GDP_BE and GDP_DE
orders_tseries <- order_integration(MacroTS[, 1:2], method = "boot_ur", B = 199,
do_parallel = FALSE, show_progress = FALSE)
```
