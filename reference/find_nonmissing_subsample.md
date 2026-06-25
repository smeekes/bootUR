# Find Non-Missing Subsamples

Find Non-Missing Subsamples

## Usage

``` r
find_nonmissing_subsample(X)
```

## Arguments

- X:

  A (\\T\\x\\N\\)-matrix of \\N\\ time series with \\T\\ observations.
  Data may also be in a time series format (e.g. `ts`, `zoo` or `xts`)
  or data frame. Assumes a prior check on missing values in-sample has
  been done.

## Value

A list with the following components

- `range`:

  (2x\\N\\)-dimensional matrix containing the first and last non-missing
  observation in each column of X.

- `all_equal`:

  Logical value indicating whether all series have the same non-missing
  indices.

## References

Smeekes, S. and Wilms, I. (2023). bootUR: An R Package for Bootstrap
Unit Root Tests. *Journal of Statistical Software*, 106(12), 1-39.
