# Augmented Dickey-Fuller Unit Root Test

This function performs a standard augmented Dickey-Fuller unit root test
on a single time series.

## Usage

``` r
adf(data, data_name = NULL, deterministics = "intercept", min_lag = 0,
  max_lag = NULL, criterion = "MAIC", criterion_scale = TRUE,
  two_step = TRUE)
```

## Arguments

- data:

  A \\T\\-dimensional vector to be tested for unit roots. Data may also
  be in a time series format (e.g. `ts`, `zoo` or `xts`), or a data
  frame.

- data_name:

  Optional name for the data, to be used in the output. The default uses
  the name of the 'data' argument.

- deterministics:

  String indicating the deterministic specification. Only relevant if
  `union = FALSE`. Options are

  `"none":` no deterministics;

  `"intercept":` intercept only;

  `"trend":` intercept and trend.

  If `union = FALSE`, the default is adding an intercept (a warning is
  given).

- min_lag:

  Minimum lag length in the augmented Dickey-Fuller regression. Default
  is 0.

- max_lag:

  Maximum lag length in the augmented Dickey-Fuller regression. Default
  uses the sample size-based rule \\12(T/100)^{1/4}\\.

- criterion:

  String for information criterion used to select the lag length in the
  augmented Dickey-Fuller regression. Options are: `"AIC"`, `"BIC"`,
  `"MAIC"`, `"MBIC"`. Default is `"MAIC"` (Ng and Perron, 2001).

- criterion_scale:

  Logical indicator whether or not to use the rescaled information
  criteria of Cavaliere et al. (2015) (`TRUE`) or not (`FALSE`). Default
  is `TRUE`.

- two_step:

  Logical indicator whether to use one-step (`two_step = FALSE`) or
  two-step (`two_step = TRUE`) detrending. The default is two-step
  detrending.

## Value

An object of class `"bootUR"`, `"htest"` with the following components:

- `method`:

  The name of the hypothesis test method;

- `data.name`:

  The name of the variable on which the method is performed;

- `null.value`:

  The value of the (gamma) parameter of the lagged dependent variable in
  the ADF regression under the null hypothesis. Under the null, the
  series has a unit root. Testing the null of a unit root then boils
  down to testing the significance of the gamma parameter;

- `alternative`:

  A character string specifying the direction of the alternative
  hypothesis relative to the null value. The alternative postulates that
  the series is stationary;

- `estimate`:

  The estimated value of the (gamma) parameter of the lagged dependent
  variable in the ADF regression;

- `statistic`:

  The value of the test statistic of the ADF unit root test;

- `p.value`:

  The p-value of the ADF unit root test.

- `specifications`:

  The specifications used in the test.

## Details

The function encompasses the standard augmented Dickey-Fuller test. The
reported p-values are MacKinnon's unit root p-values taken from the
package urca.

Lag length selection is done automatically in the ADF regression with
the specified information criterion. If one of the modified criteria of
Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is
applied. For very short time series (fewer than 50 time points) the
maximum lag length is adjusted downward to avoid potential
multicollinearity issues in the bootstrap. To overwrite data-driven lag
length selection with a pre-specified lag length, simply set both the
minimum \`min_lag\` and maximum lag length \`max_lag\` for the selection
algorithm equal to the desired lag length.

## Errors and warnings

- `Error: Multiple time series not allowed. Switch to a multivariate method such as boot_ur, or change argument data to a univariate time series.`:

  The function provides a standard ADF test with asymptotic p-value. It
  does not support multiple time series

## References

Smeekes, S. and Wilms, I. (2023). bootUR: An R Package for Bootstrap
Unit Root Tests. *Journal of Statistical Software*, 106(12), 1-39.

## Examples

``` r
# standard ADF test on GDP_BE
GDP_BE_adf <- adf(MacroTS[, 1], deterministics = "trend")
```
