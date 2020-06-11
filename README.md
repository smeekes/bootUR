
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bootUR

`bootUR` implements several bootstrap tests for unit roots, both for
single time series and for (potentially) large systems of time series.

## Installation

`bootUR` can be installed using

``` r
devtools::install_github("smeekes/bootUR")
```

## Inspect data for missing values

The bootstrap tests in `bootUR` do not work with missing data, although
multivariate time series with different start and end dates (unbalanced
panels) are allowed. `bootUR` provides a simple function to check if
your data contain missing values. We will illustrate this on the
`MacroTS` dataset of macreconomic time series that comes with the
package.

``` r
library(bootUR)
data("MacroTS")
check_missing_insample_values(MacroTS)
#>  GDP_BE  GDP_DE  GDP_FR  GDP_NL  GDP_UK CONS_BE CONS_DE CONS_FR CONS_NL CONS_UK 
#>   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
#> HICP_BE HICP_DE HICP_FR HICP_NL HICP_UK   UR_BE   UR_DE   UR_FR   UR_NL   UR_UK 
#>   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE
```

## Checking start and end points of time series

If your time series have different starting and end points (and thus
some series contain `NA`s at the beginning and/or end of your sample,
the resampling-based moving block bootstrap (MBB) and sieve bootstrap
(SB) cannot be used. `bootUR` lets you check the start and end points as
follows:

``` r
sample_check <- find_nonmissing_subsample(MacroTS)
sample_check$range  # Provides the number of the first and last non-missing observation for each series
#>       GDP_BE GDP_DE GDP_FR GDP_NL GDP_UK CONS_BE CONS_DE CONS_FR CONS_NL
#> first      1      1      1      5      1       1       1       1       5
#> last     100    100    100    100    100     100     100     100     100
#>       CONS_UK HICP_BE HICP_DE HICP_FR HICP_NL HICP_UK UR_BE UR_DE UR_FR UR_NL
#> first       1       9       9       9       9       9     1     1     1     1
#> last      100     100     100     100     100     100   100   100   100   100
#>       UR_UK
#> first     1
#> last    100
sample_check$all_equal  # Gives TRUE if the time series all start and end at the same observation
#> [1] FALSE
```

## Univariate Dickey-Fuller unit root test

To perform a standard augmented Dickey-Fuller (ADF) bootstrap unit root
test on a single time series, use the `boot_df()` function. The function
allows to set many options, including the bootstrap method used (option
`boot`), the deterministic components included (option `dc`) and the
type of detrending used. While `dc = "OLS"` gives the standard ADF test,
`dc = "QD"` provides the powerful DF-GLS test of Elliott, Rothenberg and
Stock (1996). Here we use the terminology Quasi-Differencing (QD) rather
than GLS as this conveys the meaning less ambiguously and is the same
terminology used by Smeekes and Taylor (2012) and Smeekes (2013).

### Lag selection

Lag length selection is done automatically in the ADF regression; the
default is by the modified Akaike information criterion (MAIC) proposed
by Ng and Perron (2001) with the correction of Perron and Qu (2008).
Other options include the regular Akaike information criterion (AIC), as
well as the Bayesian information criterion and its modified variant. In
addition, the rescaling suggested by Calaiere et al. (2015) is
implemented to improve the power of the test under hetereskedasticity;
this can be turned off by setting `rescale = FALSE`. To overwrite
data-driven lag length selection with a pre-specified lag length, simply
set both the minimum `p.min` and maximum lag length `p.max` for the
selection algorithm equal to the desired lag length.

### Implementation

We illustrate the bootstrap ADF test here on Dutch GDP, with the sieve
bootstrap (`boot = SB`) to as the specification used by Palm, Smeekes
and Urbain (2008) and Smeekes (2013). We set only 399 bootstrap
replications (`B = 399`) to prevent the code from running too long. We
add an intercept and a trend (`dc = 2`), and compare OLS with QD (GLS)
detrending. The option `verbose = TRUE` prints easy to read output on
the console. As random number generation is required to draw bootstrap
samples, we first set the seed of the random number generator to obtain
replicable results.

``` r
set.seed(155776)
GDP_NL <- MacroTS[, 4]
adf_out <- boot_df(GDP_NL, B = 399, boot = "SB", dc = 2, detr = c("OLS", "QD"), verbose = TRUE)
#> Bootstrap DF Test with SB bootstrap method.
#> ----------------------------------------
#> Type of unit root test performed: detr = OLS, dc = intercept and trend
#> test statistic        p-value 
#>     -2.5152854      0.1353383 
#> ----------------------------------------
#> Type of unit root test performed: detr = QD, dc = intercept and trend
#> test statistic        p-value 
#>     -1.5965001      0.4185464
```

## Union of Rejections Test

Use `boot_union()` for a test based on the union of rejections of 4
tests with different number of deterministic components and different
type of detrending (Smeekes and Taylor, 2012). The advantage of the
union test is that you don’t have to specify these (rather influential)
specification tests. This makes the union test a safe option for quick
or automatic unit root testing where careful manual specification setup
is not viable. Here we illustrate it with the sieve wild bootstrap as
proposed by Smeekes and Taylor (2012).

``` r
union_out <- boot_union(GDP_NL, B = 399, boot = "SWB", verbose = TRUE)
#> Bootstrap Test with SWB bootstrap method.
#> Bootstrap Union Test:
#> The null hypothesis of a unit root is not rejected at a significance level of 0.05.
#> test statistic        p-value 
#>     -0.6722611      0.6090226
```

## Panel Unit Root Test

The function `panel_test` performs a test on a multivariate (panel) time
series by testing the null hypothesis that all series have a unit root.
A rejection is typically interpreted as evidence that a ‘significant
proportion’ of the series is stationary, although how large that
proportion is - or which series are stationary - is not given by the
test. The test is based on averaging the individual test statistics,
also called the Group-Mean (GM) test in Palm, Smeekes and Urbain (2011).

Palm, Smeekes and Urbain (2011) introduced this test with the moving
block bootstrap (`boot = "MBB"`), which is the standard option. However,
this resampling-based method cannot handle unbalancedness, and will
therefore gives an error when applied to `MacroTS`. Therefore, you
should switch to one of the wild bootstrap methods. Here we illustrate
it with the dependent wild bootstrap (DWB) of Shao (2010) and Rho and
Shao (2019).

By default the union test is used for each series (`union = TRUE`), if
this is set to `FALSE` the deterministic components and detrending
methods can be specified as in the univariate Dickey-Fuller test.

Although the sieve bootstrap method `"SB"` and `"SWB"` can be used
(historically they have been popular among practitioners), Smeekes and
Urbain (2014b) show that these are not suited to capture general forms
of dependence across units. The code will give a warning to recommend
using a different bootstrap method.

``` r
# This will give an error!
# panel_out <- paneltest(MacroTS, boot = "MBB", B = 399, verbose = TRUE)
# This works!
panel_out <- paneltest(MacroTS, boot = "DWB", B = 399, verbose = TRUE)
#> Panel Bootstrap Group-Mean Union Test
#> The null hypothesis that all series have a unit root, is rejected at a significance level of 0.05.
#>      test statistic    p-value
#> [1,]     -0.8537718 0.03508772
```

## Individual ADF Tests for Multiple Time Series

To perform individual ADF tests on multiple time series simultaneously,
the function `iADFtest()` can be used. As the bootstrap is performed for
all series simultaneously, resampling-based bootstrap methods `"MBB"`
and `"SB"` cannot be used directly in case of unbalanced panels. If they
are used anyway, the function will revert to splitting the bootstrap up
and perfoming it individually per time series. In this case a warning is
given to alert the user. The other options are the same as for
`paneltest`.

``` r
iADF_out <- iADFtest(MacroTS[, 1:5], boot = "MBB", B = 399, verbose = TRUE, union = FALSE, dc = 2, detr = "OLS")
#> Warning in generate_inputs(y = y, BSQT_test = FALSE, iADF_test = TRUE, level =
#> level, : Missing values cause resampling bootstrap to be executed for each time
#> series individually.
#> ----------------------------------------
#> Type of unit root test performed: detr = OLS, dc = intercept and trend
#> There are 0 stationary time series
#>        test statistic   p-value
#> GDP_BE      -2.792169 0.2330827
#> GDP_DE      -2.774320 0.0877193
#> GDP_FR      -2.048760 0.4837093
#> GDP_NL      -2.515285 0.1954887
#> GDP_UK      -2.449065 0.2857143
```

Note that `iADFtest` (intentionally) does not provide a correction for
multiple testing; of course, if we perform each test with a signifiance
level of 5%, the probability of making a mistake in all these tests
becomes (much, if `N` is large) more than 5%. To explicitly account for
multiple testing, use the functions `BSQTtest()` or `bFDRtest()`.

## Sequential Tests for Multiple Time Series

The function `BSQTtest()` performs the Bootstrap Sequential Quantile
Test (BSQT) proposed by Smeekes (2015). Here we split the series in
groups which are consecutively tested for unit roots, starting with the
group most likely to be stationary (having the smallest ADF statistics).
If the unit root hypothesis cannot be rejected for the first group, the
algorithm stops; if there is a rejection, the second group is tested,
and so on.

Most options are the same as for `paneltest`. The most important new
parameter to set here is the group sizes. These can either be set in
units, or in fractions of the total number of series (i.e. quantiles,
hence the name) via the parameter `q`. If we have `N` time series,
setting `q = 0:N` means each unit should be tested sequentially. To
split the series in four equally sized groups (regardless of many series
there are), use `q = 0:4 / 4`. By convention and in accordance with
notation in Smeekes (2015), the first entry of the vector should be
equal to zero, while the second entry indicates the start of the second
group, and so on. However, if the initial zero is accidentally omitted,
it is automatically added by the function. Similarly, if the final value
is not equal to `1` (in case of quantiles) or `N` to end the last group,
this is added by the function.

Testing individual series consecutively is easiest for interpretation,
but is only meaningful if `N` is small. In this case the method is
equivalent to the bootstrap StepM method of Romano and Wolf (2005),
which controls the familywise error rate, that is the probability of
making at least one false rejection. This can get very conservative if
`N` is large, and you would typically end up not rejecting any null
hypothesis. The method is illustrated with the autoregressive wild
bootstrap of Smeekes and Urbain (2014a) and Friedrich, Smeekes and
Urbain (2020).

``` r
N <- ncol(MacroTS)
# Test each unit sequentially
BSQT_out1 <- BSQTtest(MacroTS, q = 0:N, boot = "AWB", B = 399, verbose = TRUE)
#> There are 2 stationary time series, namely: HICP_BE HICP_DE.
#> Details of the BSQT ssquential tests:
#>        Unit H0 Unit H1 Test statistic    p-value
#> Step 1       0       1      -1.619140 0.01754386
#> Step 2       1       2      -1.550431 0.03007519
#> Step 3       2       3      -1.340742 0.18546366
# Split in four equally sized groups (motivated by the 4 series per country)
BSQT_out2 <- BSQTtest(MacroTS, q = 0:4 / 4, boot = "AWB", B = 399, verbose = TRUE)
#> There are 0 stationary time series.
#> Details of the BSQT ssquential tests:
#>        Unit H0 Unit H1 Test statistic    p-value
#> Step 1       0       5       -1.00877 0.05513784
```

## Bootstrap FDR controlling tests

The function `bFDRtest()` controls for multiple testing by controlling
the false discovery rate (FDR), which is defined as the expected
proportion of false rejections relative to the total number of
rejections. This scales with the total number of tests, making it more
suitable for large `N` than the familywise error rate.

The bootstrap method for controlling FDR was introduced by Romano,
Shaikh and Wolf (2008), who showed that, unlike the classical way to
control FDR, the bootstrap is appropriate under general forms of
dependence between series. Moon and Perron (2012) applied this method to
unit root testing; it is essentially their method which is implemented
in `bFDRtest()` though again with the option to change the bootstrap
used (their suggestion was MBB). The arguments to be set are the same as
for the other multivariate unit root tests, though the meaning of
`level` changes from regular significance level to FDR level. As BSQT,
the method only report those tests until no rejection occurs.

We illustrate it here with the final available bootstrap method, the
block wild bootstrap of Shao (2011) and Smeekes and Urbain (2014a).

``` r
N <- ncol(MacroTS)
bFDR_out <- bFDRtest(MacroTS[, 1:10], level = 0.1, boot = "BWB", B = 399, verbose = TRUE)
#> There are 0 stationary time series
#> Details of the FDR sequential tests:
#>        test statistic critical value
#> GDP_DE      -1.060895      -1.577345
```

## References

  - Cavaliere, G., Phillips, P.C.B., Smeekes, S., and Taylor, A.M.R.
    (2015). Lag length selection for unit root tests in the presence of
    nonstationary volatility. *Econometric Reviews*, 34(4), 512-536.
  - Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996). Efficient
    tests for an autoregressive unit root. *Econometrica*, 64(4),
    813-836.
  - Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive
    wild bootstrap inference for nonparametric trends. *Journal of
    Econometrics*, 214(1), 81-109.
  - Moon, H.R. and Perron, B. (2012). Beyond panel unit root tests:
    Using multiple testing to determine the non stationarity properties
    of individual series in a panel. *Journal of Econometrics*, 169(1),
    29-33.
  - Ng, S. and Perron, P. (2001). Lag Length Selection and the
    Construction of Unit Root Tests with Good Size and Power.
    *Econometrica*, 69(6), 1519-1554,
  - Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008). Bootstrap unit
    root tests: Comparison and extensions. *Journal of Time Series
    Analysis*, 29(1), 371-401.
  - Palm, F. C., Smeekes, S., and Urbain, J.-.P. (2011). Cross-sectional
    dependence robust block bootstrap panel unit root tests. *Journal of
    Econometrics*, 163(1), 85-104.
  - Paparoditis, E. and Politis, D.N. (2003). Residual‐based block
    bootstrap for unit root testing. *Econometrica*, 71(3), 813-855.
  - Perron, P. and Qu, Z. (2008). A simple modification to improve the
    finite sample properties of Ng and Perron’s unit root tests.
    *Economic Letters*, 94(1), 12-19.
  - Rho, Y. and Shao, X. (2019). Bootstrap-assisted unit root testing
    with piecewise locally stationary errors. *Econometric Theory*,
    35(1), 142-166.
  - Romano, J.P., Shaikh, A.M., and Wolf, M. (2008). Control of the
    false discovery rate under dependence using the bootstrap and
    subsampling. *Test*, 17(3), 417.
  - Romano, J. P. and Wolf, M. (2005). Stepwise multiple testing as
    formalized data snooping. *Econometrica*, 73(4), 1237-1282.
  - Shao, X. (2010). The dependent wild bootstrap. *Journal of the
    American Statistical Association*, 105(489), 218-235.
  - Shao, X. (2011). A bootstrap-assisted spectral test of white noise
    under unknown dependence. *Journal of Econometrics*, 162, 213-224.
  - Smeekes (2013)
  - Smeekes, S. (2015). Bootstrap sequential tests to determine the
    order of integration of individual units in a time series panel.
    *Journal of Time Series Analysis*, 36(3), 398-415.
  - Smeekes, S. and Taylor, A.M.R. (2012). Bootstrap union tests for
    unit roots in the presence of nonstationary volatility. *Econometric
    Theory*, 28(2), 422-456.
  - Smeekes, S. and Urbain, J.-P. (2014a). A multivariate invariance
    principle for modified wild bootstrap methods with an application to
    unit root testing. GSBE Research Memorandum No. RM/14/008,
    Maastricht University.
  - Smeekes, S. and Urbain, J.-P. (2014b). On the applicability of the
    sieve bootstrap in time series panels. *Oxford Bulletin of Economics
    and Statistics*, 76(1), 139-151.
