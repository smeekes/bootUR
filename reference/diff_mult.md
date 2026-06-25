# Differences of Multiple Time Series

Performs differencing of multiple time series, with possibly different
orders for each time series.

## Usage

``` r
diff_mult(data, d, keep_NAs = TRUE)
```

## Arguments

- data:

  A (\\T\\x\\N\\)-matrix of \\N\\ time series with \\T\\ observations.
  Data may also be in a time series format (e.g. `ts`, `zoo` or `xts`)
  or data frame.

- d:

  An \\N\\-dimensional vector containing the orders

- keep_NAs:

  Logical indicator whether or not to keep the `NA` values resulting
  from differencing at the beginning of the sample. Default is `TRUE`.
  If `FALSE`, the entire row containing the `NA` values is removed.

## Value

The appropriately differenced data in the same format as the original
data.

## References

Smeekes, S. and Wilms, I. (2023). bootUR: An R Package for Bootstrap
Unit Root Tests. *Journal of Statistical Software*, 106(12), 1-39.
