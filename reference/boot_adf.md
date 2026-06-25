# Bootstrap augmented Dickey-Fuller Unit Root Test

This function performs a standard augmented Dickey-Fuller bootstrap unit
root test on a single time series.

## Usage

``` r
boot_adf(data, data_name = NULL, bootstrap = "AWB", B = 1999,
  block_length = NULL, ar_AWB = NULL, deterministics = "intercept",
  min_lag = 0, max_lag = NULL, criterion = "MAIC", detrend = "OLS",
  criterion_scale = TRUE, show_progress = TRUE, do_parallel = TRUE,
  cores = NULL)
```

## Arguments

- data:

  A \\T\\-dimensional vector to be tested for unit roots. Data may also
  be in a time series format (e.g. `ts`, `zoo` or `xts`), or a data
  frame.

- data_name:

  Optional name for the data, to be used in the output. The default uses
  the name of the 'data' argument.

- bootstrap:

  String for bootstrap method to be used. Options are

  `"MBB"`

  :   Moving blocks bootstrap (Paparoditis and Politis, 2003);

  `"BWB"`

  :   Block wild bootstrap (Shao, 2011);

  `"DWB"`

  :   Dependent wild bootstrap (Shao, 2010; Rho and Shao, 2019);

  `"AWB"`

  :   Autoregressive wild bootstrap (Smeekes and Urbain, 2014a;
      Friedrich, Smeekes and Urbain, 2020), this is the default;

  `"SB"`

  :   Sieve bootstrap (Chang and Park, 2003; Palm, Smeekes and Urbain,
      2008; Smeekes, 2013);

  `"SWB"`

  :   Sieve wild bootstrap (Cavaliere and Taylor, 2009; Smeekes and
      Taylor, 2012).

- B:

  Number of bootstrap replications. Default is 1999.

- block_length:

  Desired 'block length' in the bootstrap. For the MBB, BWB and DWB
  bootstrap, this is a genuine block length. For the AWB bootstrap, the
  block length is transformed into an autoregressive parameter via the
  formula \\0.01^(1/block_length)\\ as in Smeekes and Urbain (2014a);
  this can be overwritten by setting `ar_AWB` directly. Default sets the
  block length as a function of the time series length T, via the rule
  \\block_length = 1.75 T^(1/3)\\ of Palm, Smeekes and Urbain (2011).

- ar_AWB:

  Autoregressive parameter used in the AWB bootstrap method
  (`bootstrap = "AWB"`). Can be used to set the parameter directly
  rather than via the default link to the block length.

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

- detrend:

  String indicating the type of detrending to be performed. Only
  relevant if `union = FALSE`. Options are: `"OLS"` or `"QD"` (typically
  also called GLS, see Elliott, Rothenberg and Stock, 1996). The default
  is `"OLS"`.

- criterion_scale:

  Logical indicator whether or not to use the rescaled information
  criteria of Cavaliere et al. (2015) (`TRUE`) or not (`FALSE`). Default
  is `TRUE`.

- show_progress:

  Logical indicator whether a bootstrap progress update should be
  printed to the console. Default is FALSE.

- do_parallel:

  Logical indicator whether bootstrap loop should be executed in
  parallel. Default is TRUE.

- cores:

  The number of cores to be used in the parallel loops. Default is to
  use all but one.

## Value

An object of class `"bootUR"`, `"htest"` with the following components:

- `method`:

  The name of the hypothesis test method;

- `data.name`:

  The name of the data on which the method is performed;

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
  variable in the ADF regression.;

- `statistic`:

  The value of the test statistic of the unit root test;

- `p.value`:

  The p-value of the unit root test;

- `details`:

  A list containing the detailed outcomes of the performed test, such as
  selected lags, individual estimates and p-values.

- `specifications`:

  The specifications used in the test.

## Details

The options encompass many test proposed in the literature.
`detrend = "OLS"` gives the standard augmented Dickey-Fuller test, while
`detrend = "QD"` provides the DF-GLS test of Elliott, Rothenberg and
Stock (1996). The bootstrap algorithm is always based on a residual
bootstrap (under the alternative) to obtain residuals rather than a
difference-based bootstrap (under the null), see e.g. Palm, Smeekes and
Urbain (2008).

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

  The function is a simple wrapper around
  [`boot_ur`](https://smeekes.github.io/bootUR/reference/boot_ur.md) to
  facilitate use for single time series. It does not support multiple
  time series, as
  [`boot_ur`](https://smeekes.github.io/bootUR/reference/boot_ur.md) is
  specifically suited for that.

## References

Smeekes, S. and Wilms, I. (2023). bootUR: An R Package for Bootstrap
Unit Root Tests. *Journal of Statistical Software*, 106(12), 1-39.

Chang, Y. and Park, J. (2003). A sieve bootstrap for the test of a unit
root. *Journal of Time Series Analysis*, 24(4), 379-400.

Cavaliere, G. and Taylor, A.M.R (2009). Heteroskedastic time series with
a unit root. *Econometric Theory*, 25, 1228–1276.

Cavaliere, G., Phillips, P.C.B., Smeekes, S., and Taylor, A.M.R. (2015).
Lag length selection for unit root tests in the presence of
nonstationary volatility. *Econometric Reviews*, 34(4), 512-536.

Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996). Efficient tests
for an autoregressive unit root. *Econometrica*, 64(4), 813-836.

Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild
bootstrap inference for nonparametric trends. *Journal of Econometrics*,
214(1), 81-109.

Ng, S. and Perron, P. (2001). Lag Length Selection and the Construction
of Unit Root Tests with Good Size and Power. *Econometrica*, 69(6),
1519-1554,

Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008). Bootstrap unit root
tests: Comparison and extensions. *Journal of Time Series Analysis*,
29(1), 371-401.

Paparoditis, E. and Politis, D.N. (2003). Residual-based block bootstrap
for unit root testing. *Econometrica*, 71(3), 813-855.

Perron, P. and Qu, Z. (2008). A simple modification to improve the
finite sample properties of Ng and Perron's unit root tests. *Economic
Letters*, 94(1), 12-19.

Rho, Y. and Shao, X. (2019). Bootstrap-assisted unit root testing with
piecewise locally stationary errors. *Econometric Theory*, 35(1),
142-166.

Smeekes, S. (2013). Detrending bootstrap unit root tests. *Econometric
Reviews*, 32(8), 869-891.

Shao, X. (2010). The dependent wild bootstrap. *Journal of the American
Statistical Association*, 105(489), 218-235.

Shao, X. (2011). A bootstrap-assisted spectral test of white noise under
unknown dependence. *Journal of Econometrics*, 162, 213-224.

Smeekes, S. and Taylor, A.M.R. (2012). Bootstrap union tests for unit
roots in the presence of nonstationary volatility. *Econometric Theory*,
28(2), 422-456.

Smeekes, S. and Urbain, J.-P. (2014a). A multivariate invariance
principle for modified wild bootstrap methods with an application to
unit root testing. GSBE Research Memorandum No. RM/14/008, Maastricht
University

## See also

[`boot_ur`](https://smeekes.github.io/bootUR/reference/boot_ur.md)

## Examples

``` r
# boot_adf on GDP_BE
GDP_BE_adf <- boot_adf(MacroTS[, 1], B = 199, deterministics = "trend", detrend = "OLS",
                       do_parallel = FALSE, show_progress = FALSE)
print(GDP_BE_adf)
#> 
#>  AWB bootstrap OLS test (with intercept and trend) on a single time
#>  series
#> 
#> data: MacroTS[, 1]
#> null hypothesis: Series has a unit root
#> alternative hypothesis: Series is stationary
#> 
#>              estimate largest root statistic p-value
#> MacroTS[, 1]                0.9304    -2.792  0.1106
#> 
```
