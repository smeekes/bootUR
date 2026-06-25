# Individual Unit Root Tests without multiple testing control

This function performs bootstrap unit root tests on each time series
individually.

## Usage

``` r
boot_ur(data, data_name = NULL, bootstrap = "AWB", B = 1999,
  block_length = NULL, ar_AWB = NULL, level = NULL, union = TRUE,
  union_quantile = 0.05, deterministics = NULL, detrend = NULL,
  min_lag = 0, max_lag = NULL, criterion = "MAIC",
  criterion_scale = TRUE, show_progress = TRUE, do_parallel = TRUE,
  cores = NULL)
```

## Arguments

- data:

  A \\T\\-dimensional vector or a (\\T\\ x \\N\\)-matrix of \\N\\ time
  series with \\T\\ observations to be tested for unit roots. Data may
  also be in a time series format (e.g. `ts`, `zoo` or `xts`), or a data
  frame, as long as each column represents a single time series.

- data_name:

  Optional name for the data, to be used in the output. The default uses
  the name of the 'data' argument.

- bootstrap:

  String for bootstrap method to be used. Options are

  `"MBB"`

  :   Moving block bootstrap (Paparoditis and Politis, 2003; Palm,
      Smeekes and Urbain, 2011);

  `"BWB"`

  :   Block wild bootstrap (Shao, 2011; Smeekes and Urbain, 2014a);

  `"DWB"`

  :   Dependent wild bootstrap (Shao, 2010; Smeekes and Urbain, 2014a;
      Rho and Shao, 2019);

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

- level:

  The desired significance level of the test (optional). This is only
  used for multivariate series to be able to provide a boolean vector
  with rejections of the null hypothesis or not for easy
  post-processing. Default is `NULL`, in which case no such vector is
  given.

- union:

  Logical indicator whether or not to use bootstrap union tests (`TRUE`)
  or not (`FALSE`), see Smeekes and Taylor (2012). Default is `TRUE`.

- union_quantile:

  The quantile of the bootstrap distribution used for scaling the
  individual statistics in the union. Ideally this should equal the
  desired significance level of the test. Default is 0.05. This
  parameter is overwritten when a significance level is provided in the
  argument `level`.

- deterministics:

  String indicating the deterministic specification. Only relevant if
  `union = FALSE`. Options are

  `"none":` no deterministics;

  `"intercept":` intercept only;

  `"trend":` intercept and trend.

  If `union = FALSE`, the default is adding an intercept (a warning is
  given).

- detrend:

  String indicating the type of detrending to be performed. Only
  relevant if `union = FALSE`. Options are: `"OLS"` or `"QD"` (typically
  also called GLS, see Elliott, Rothenberg and Stock, 1996). The default
  is `"OLS"`.

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

An object of class `"bootUR"`, `"\*"`, where `"\*"` is `"mult_htest"`
for multiple time series or `"htest"` for single time series, with the
following components:

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

  The estimated value(s) of the (gamma) parameter of the lagged
  dependent variable in the ADF regressions. Note that for the union
  test (`union = TRUE`), this estimate is not defined, hence `NA` is
  returned;

- `statistic`:

  The value(s) of the test statistic of the unit root test(s);

- `p.value`:

  The p-value(s) of the unit root test(s);

- `rejections`:

  For `"mult_htest"` only. A vector with logical indicators for each
  time series whether the null hypothesis of a unit root is rejected
  (`TRUE`) or not (`FALSE`). This is only supplied when an optional
  significance level is given, otherwise `NULL` is returned;

- `details`:

  A list containing the detailed outcomes of the performed tests, such
  as selected lags, individual estimates and p-values.

- `series.names`:

  For `"mult_htest"` only. The names of the series that the tests are
  performed on;

- `specifications`:

  The specifications used in the test(s).

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

## Warnings

The function may give the following warnings.

- `Warning: Missing values cause resampling bootstrap to be executed for each time series individually.`:

  If the time series in `data` have different starting and end points
  (and thus some series contain `NA` values at the beginning and/or end
  of the sample, the resampling-based moving block bootstrap (MBB) and
  sieve bootstrap (SB) cannot be used directly, as they create holes
  (internal missings) in the bootstrap samples. These bootstrap methods
  are therefore not applied jointly as usual, but individually to each
  series.

- `Warning: Deterministic specification in argument deterministics is ignored, as union test is applied.`:

  The union test calculates the union of all four combinations of
  deterministic components (intercept or intercept and trend) and
  detrending methods (OLS or QD). Setting deterministic components
  manually therefore has no effect.

- `Warning: Detrending method in argument detrend is ignored, as union test is applied.`:

  The union test calculates the union of all four combinations of
  deterministic components (intercept or intercept and trend) and
  detrending methods (OLS or QD). Setting detrending methods manually
  therefore has no effect.

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

Palm, F. C., Smeekes, S., and Urbain, J.-.P. (2011). Cross-sectional
dependence robust block bootstrap panel unit root tests. *Journal of
Econometrics*, 163(1), 85-104.

Paparoditis, E. and Politis, D.N. (2003). Residual-based block bootstrap
for unit root testing. *Econometrica*, 71(3), 813-855.

Perron, P. and Qu, Z. (2008). A simple modification to improve the
finite sample properties of Ng and Perron's unit root tests. *Economic
Letters*, 94(1), 12-19.

Rho, Y. and Shao, X. (2019). Bootstrap-assisted unit root testing with
piecewise locally stationary errors. *Econometric Theory*, 35(1),
142-166.

Shao, X. (2010). The dependent wild bootstrap. *Journal of the American
Statistical Association*, 105(489), 218-235.

Shao, X. (2011). A bootstrap-assisted spectral test of white noise under
unknown dependence. *Journal of Econometrics*, 162, 213-224.

Smeekes, S. (2013). Detrending bootstrap unit root tests. *Econometric
Reviews*, 32(8), 869-891.

Smeekes, S. and Taylor, A.M.R. (2012). Bootstrap union tests for unit
roots in the presence of nonstationary volatility. *Econometric Theory*,
28(2), 422-456.

Smeekes, S. and Urbain, J.-P. (2014a). A multivariate invariance
principle for modified wild bootstrap methods with an application to
unit root testing. GSBE Research Memorandum No. RM/14/008, Maastricht
University

## Examples

``` r
# boot_ur on GDP_BE and GDP_DE
two_series_boot_ur <- boot_ur(MacroTS[, 1:2], bootstrap = "MBB", B = 199,
                              do_parallel = FALSE, show_progress = FALSE)
print(two_series_boot_ur)
#> 
#>  MBB bootstrap union test on each individual series (no multiple testing
#>  correction)
#> 
#> data: MacroTS[, 1:2]
#> null hypothesis: Series has a unit root
#> alternative hypothesis: Series is stationary
#> 
#> Tests performed on each series: 
#>        estimate largest root statistic p-value
#> GDP_BE                    NA   -0.8123 0.34673
#> GDP_DE                    NA   -1.0958 0.08543
#> 
```
