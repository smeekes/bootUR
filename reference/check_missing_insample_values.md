# Check Missing Values in Sample

Check Missing Values in Sample

## Usage

``` r
check_missing_insample_values(X)
```

## Arguments

- X:

  A (\\T\\x\\N\\)-matrix of \\N\\ time series with \\T\\ observations.
  Data may also be in a time series format (e.g. `ts`, `zoo` or `xts`)
  or data frame.

## Value

\\N\\-dimensional vector, for each series whether missing values are
present (`TRUE`) or not (`FALSE`)

## References

Smeekes, S. and Wilms, I. (2023). bootUR: An R Package for Bootstrap
Unit Root Tests. *Journal of Statistical Software*, 106(12), 1-39.
