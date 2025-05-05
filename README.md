
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bootUR: Bootstrap Unit Root Tests

<!-- badges: start -->

[![CRAN_Version_Badge](https://www.r-pkg.org/badges/version/bootUR)](https://cran.r-project.org/package=bootUR)
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/grand-total/bootUR)](https://cran.r-project.org/package=bootUR)
[![License_GPLv2_Badge](https://img.shields.io/badge/License-GPLv2-yellow.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
[![License_GPLv3_Badge](https://img.shields.io/badge/License-GPLv3-yellow.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![R-CMD-check](https://github.com/smeekes/bootUR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/smeekes/bootUR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The R package `bootUR` implements several bootstrap tests for unit
roots, both for single time series and for (potentially) large systems
of time series.

## Installation and Loading

### Installation

The package can be installed from CRAN using

``` r
install.packages("bootUR")
```

The development version of the `bootUR` package can be installed from
GitHub using

``` r
# install.packages("devtools")
devtools::install_github("smeekes/bootUR")
```

When installing from GitHub, in order to build the package from source,
you need to have the appropriate R development tools installed
([Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, or
[these tools](https://mac.r-project.org/tools/) on Mac).

If you want the vignette to appear in your package when installing from
GitHub, use

``` r
# install.packages("devtools")
devtools::install_github("smeekes/bootUR", build_vignettes = TRUE, dependencies = TRUE)
```

instead. As building the vignette may take a bit of time (all bootstrap
code below is run), package installation will be slower this way.

### Load Package

After installation, the package can be loaded in the standard way:

``` r
library(bootUR)
```

## Preliminary Analysis: Missing Values

`bootUR` provides a few simple tools to check if your data are suitable
to be bootstrapped.

### Inspect Data for Missing Values

The bootstrap tests in `bootUR` do not work with missing data, although
multivariate time series with different start and end dates (unbalanced
panels) are allowed. `bootUR` provides a simple function to check if
your data contain missing values. We will illustrate this on the
`MacroTS` dataset of macroeconomic time series that comes with the
package.

``` r
data("MacroTS")
check_missing_insample_values(MacroTS)
#>  GDP_BE  GDP_DE  GDP_FR  GDP_NL  GDP_UK CONS_BE CONS_DE CONS_FR CONS_NL CONS_UK 
#>   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE 
#> HICP_BE HICP_DE HICP_FR HICP_NL HICP_UK   UR_BE   UR_DE   UR_FR   UR_NL   UR_UK 
#>   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE   FALSE
```

### Checking Start and End Points of Time Series

If your time series have different starting and end points (and thus
some series contain `NA`s at the beginning and/or end of your sample,
the resampling-based moving block bootstrap (MBB) and sieve bootstrap
(SB) cannot be used. `bootUR` lets you check the start and end points as
follows:

``` r
sample_check <- find_nonmissing_subsample(MacroTS)
# Provides the number of the first and last non-missing observation for each series:
sample_check$range 
#>       GDP_BE GDP_DE GDP_FR GDP_NL GDP_UK CONS_BE CONS_DE CONS_FR CONS_NL
#> first      1      1      1      5      1       1       1       1       5
#> last     100    100    100    100    100     100     100     100     100
#>       CONS_UK HICP_BE HICP_DE HICP_FR HICP_NL HICP_UK UR_BE UR_DE UR_FR UR_NL
#> first       1       9       9       9       9       9     1     1     1     1
#> last      100     100     100     100     100     100   100   100   100   100
#>       UR_UK
#> first     1
#> last    100
# Gives TRUE if the time series all start and end at the same observation:
sample_check$all_equal
#> [1] FALSE
```

### Visualizing Missing Data

If you have `ggplot2` installed, you can also plot the missing data
patterns in your series to get a quick overview. You may need to
manipulate some arguments to get the plot properly sized (therefore it
is not run here automatically).

``` r
plot_missing_values(MacroTS, show_names = TRUE, axis_text_size = 5, legend_size = 6)
```

## Augmented Dickey-Fuller Test

As the standard test for unit roots, `bootUR` also has an implementation
of the standard, non-bootstrap, augmented Dickey-Fuller (ADF) test
(though its use is not recommended if sample sizes are small). For this
purpose the `adf()` function can be used. The function allows to set
many options. First, one can choose between the classical single-step
procedure (`two_step = FALSE`), in which deterministic components are
directly included in the test regression, and the more flexible and
modern two-step procedure (`two_step = TRUE`) where deterministic
components are first removed before applying the unit root test to
detrended data. For the standard ADF test, the two specifications
generally yield nearly identical results.

**Lag selection**

Lag length selection is done automatically in the ADF regression; the
default is by the modified Akaike information criterion (MAIC) proposed
by Ng and Perron (2001) with the correction of Perron and Qu (2008).
Other options include the regular Akaike information criterion (AIC), as
well as the Bayesian information criterion and its modified variant. In
addition, the rescaling suggested by Cavaliere et al. (2015) is
implemented to improve the power of the test under heteroskedasticity;
this can be turned off by setting `criterion_scale = FALSE`. To
overwrite data-driven lag length selection with a pre-specified lag
length, simply set both the minimum `min_lag` and maximum lag length
`max_lag` for the selection algorithm equal to the desired lag length.

**Implementation**

We illustrate the ADF test here on Dutch GDP for the two-step
specification, including a linear trend in the specification.

``` r
GDP_NL <- MacroTS[, 4]
adf(GDP_NL, deterministics = "trend")
#> 
#>  Two-step ADF test (with trend) on a single time series
#> 
#> data: GDP_NL
#> null hypothesis: Series has a unit root
#> alternative hypothesis: Series is stationary
#> 
#>        estimate largest root statistic p-value
#> GDP_NL                0.9471    -2.515  0.3202
```

## Univariate Bootstrap Unit Root Tests

### Augmented Dickey-Fuller Test

To perform a bootstrap version of the ADF unit root test on a single
time series, use the `boot_adf()` function. The function allows to set
many options, including the bootstrap method used (option `bootstrap`),
the deterministic components included (option `deterministics`) and the
type of detrending used (option `detrend`). While `detrend = "OLS"`
gives the standard ADF test, `detrend = "QD"` provides the powerful
DF-GLS test of Elliott, Rothenberg and Stock (1996). Here we use the
terminology Quasi-Differencing (QD) rather than GLS as this conveys the
meaning less ambiguously and is the same terminology used by Smeekes and
Taylor (2012) and Smeekes (2013). In all cases, two-step detrending is
used.

**Implementation**

We illustrate the bootstrap ADF test here on Dutch GDP, with the sieve
bootstrap (`bootstrap = SB`) as in the specification used by Palm,
Smeekes and Urbain (2008) and Smeekes (2013). To get the well-known test
proposed by Paparoditis and Politis (2003), set `bootstrap = "MBB"`. We
set only 399 bootstrap replications (`B = 399`) to prevent the code from
running too long. We add an intercept and a trend
(`deterministics = "trend"`) and OLS detrending. The console gives you
live updates on the bootstrap progress. To turn these off, set
`show_progress = FALSE`. The bootstrap loop can be run in parallel by
setting `do_parallel = TRUE` (the default).

As random number generation is required to draw bootstrap samples, we
first set the seed of the random number generator to obtain replicable
results.

``` r
set.seed(155776)
boot_adf(GDP_NL, B = 399, bootstrap = "SB", deterministics = "trend", 
                    detrend = "OLS", do_parallel = FALSE)
#> Progress: |------------------| 
#>           ********************
#> 
#>  SB bootstrap OLS test (with intercept and trend) on a single time
#>  series
#> 
#> data: GDP_NL
#> null hypothesis: Series has a unit root
#> alternative hypothesis: Series is stationary
#> 
#>        estimate largest root statistic p-value
#> GDP_NL                0.9471    -2.515  0.1454
```

### Union of Rejections Test

Use `boot_union()` for a test based on the union of rejections of 4
tests with different number of deterministic components and different
type of detrending (Smeekes and Taylor, 2012). The advantage of the
union test is that you don’t have to specify these (rather influential)
specification tests. This makes the union test a safe option for quick
or automatic unit root testing where careful manual specification setup
is not viable. Here we illustrate it with the sieve wild bootstrap as
proposed by Smeekes and Taylor (2012).

``` r
boot_union(GDP_NL, B = 399, bootstrap = "SWB", do_parallel = FALSE)
#> Progress: |------------------| 
#>           ********************
#> 
#>  SWB bootstrap union test on a single time series
#> 
#> data: GDP_NL
#> null hypothesis: Series has a unit root
#> alternative hypothesis: Series is stationary
#> 
#>        estimate largest root statistic p-value
#> GDP_NL                    NA   -0.7115   0.614
```

## Panel Unit Root Test

The function `boot_panel` performs a test on a multivariate (panel) time
series by testing the null hypothesis that all series have a unit root.
A rejection is typically interpreted as evidence that a ‘significant
proportion’ of the series is stationary, although how large that
proportion is - or which series are stationary - is not given by the
test. The test is based on averaging the individual test statistics,
also called the Group-Mean (GM) test in Palm, Smeekes and Urbain (2011).

Palm, Smeekes and Urbain (2011) introduced this test with the moving
block bootstrap (`bootstrap = "MBB"`). However, this resampling-based
method cannot handle unbalancedness, and will therefore give an error
when applied to `MacroTS`:

``` r
boot_panel(MacroTS, bootstrap = "MBB", B = 399, do_parallel = FALSE)
#> Error in check_inputs(data = data, boot_sqt_test = boot_sqt_test, boot_ur_test = boot_ur_test, : Resampling-based bootstraps MBB and SB cannot handle unbalanced series.
```

Therefore, you should switch to one of the wild bootstrap methods. Here
we illustrate it with the dependent wild bootstrap (DWB) of Shao (2010)
and Rho and Shao (2019).

By default the union test is used for each series (`union = TRUE`), if
this is set to `FALSE` the deterministic components and detrending
methods can be specified as in the univariate Dickey-Fuller test.

Although the sieve bootstrap method `"SB"` and `"SWB"` can be used
(historically they have been popular among practitioners), Smeekes and
Urbain (2014b) show that these are not suited to capture general forms
of dependence across units. The code will give a warning to recommend
using a different bootstrap method.

``` r
boot_panel(MacroTS, bootstrap = "DWB", B = 399, do_parallel = FALSE)
#> Progress: |------------------| 
#>           ********************
#> 
#>  Panel DWB bootstrap group-mean union test
#> 
#> data: MacroTS
#> null hypothesis: All series have a unit root
#> alternative hypothesis: Some series are stationary
#> 
#>         estimate largest root statistic p-value
#> MacroTS                    NA   -0.8621  0.1103
```

## Tests for Multiple Time Series

### Individual ADF Tests

To perform individual ADF tests on multiple time series simultaneously,
the function `boot_ur()` can still be used. As the bootstrap is
performed for all series simultaneously, resampling-based bootstrap
methods `"MBB"` and `"SB"` cannot be used directly in case of unbalanced
panels. If they are used anyway, the function will revert to splitting
the bootstrap up and performing it individually per time series. In this
case a warning is given to alert the user.

``` r
ADFtests_out <- boot_ur(MacroTS[, 1:5], bootstrap = "MBB", B = 399, union = FALSE, 
                        deterministics = "trend", detrend = "OLS", do_parallel = FALSE)
#> Warning in check_inputs(data = data, boot_sqt_test = boot_sqt_test,
#> boot_ur_test = boot_ur_test, : Missing values cause resampling bootstrap to be
#> executed for each time series individually.
#> Progress: |------------------| 
#>           ********************
print(ADFtests_out)
#> 
#>  MBB bootstrap ADF test (with intercept and trend) on each individual
#>  series (no multiple testing correction)
#> 
#> data: MacroTS[, 1:5]
#> null hypothesis: Series has a unit root
#> alternative hypothesis: Series is stationary
#> 
#> Tests performed on each series: 
#>        estimate largest root statistic p-value
#> GDP_BE                0.9304    -2.792  0.1955
#> GDP_DE                0.8911    -2.774  0.1003
#> GDP_FR                0.9655    -2.049  0.5113
#> GDP_NL                0.9471    -2.515  0.1930
#> GDP_UK                0.9600    -2.449  0.2882
```

Note that `boot_ur` (intentionally) does not provide a correction for
multiple testing; of course, if we perform each test with a significance
level of 5%, the probability of making a mistake in all these tests
becomes (much, if `N` is large) more than 5%. To explicitly account for
multiple testing, use the functions `boot_sqt()` or `boot_fdr()`.

### Bootstrap Sequential Tests

The function `boot_sqt()` performs the Bootstrap Sequential Quantile
Test (BSQT) proposed by Smeekes (2015). Here we split the series in
groups which are consecutively tested for unit roots, starting with the
group most likely to be stationary (having the smallest ADF statistics).
If the unit root hypothesis cannot be rejected for the first group, the
algorithm stops; if there is a rejection, the second group is tested,
and so on.

Most options are the same as for `boot_panel`. The parameter `SQT_level`
controls the significance level of the individual tests performed in the
sequence, with a default value of 0.05. The other important new
parameter to set here is the group sizes. These can either be set in
units, or in fractions of the total number of series (i.e. quantiles,
hence the name) via the parameter `steps`. If we have `N` time series,
setting `steps = 0:N` means each unit should be tested sequentially. To
split the series in four equally sized groups (regardless of many series
there are), use `steps = 0:4 / 4`. By convention and in accordance with
notation in Smeekes (2015), the first entry of the vector should be
equal to zero, while the second entry indicates the end of the first
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
boot_sqt(MacroTS, steps = 0:N, bootstrap = "AWB", B = 399, do_parallel = FALSE)
#> Progress: |------------------| 
#>           ********************
#> 
#>  AWB bootstrap sequential quantile union test
#> 
#> data: MacroTS
#> null hypothesis: Series has a unit root
#> alternative hypothesis: Series is stationary
#> 
#> Sequence of tests: 
#>        H0: # I(0) H1: # I(0)  tstat p-value
#> Step 1          0          1 -1.661 0.02256
#> Step 2          1          2 -1.413 0.12281
# Split in four equally sized groups (motivated by the 4 series per country)
boot_sqt(MacroTS, steps = 0:4 / 4, bootstrap = "AWB", B = 399, do_parallel = FALSE)
#> Progress: |------------------| 
#>           ********************
#> 
#>  AWB bootstrap sequential quantile union test
#> 
#> data: MacroTS
#> null hypothesis: Series has a unit root
#> alternative hypothesis: Series is stationary
#> 
#> Sequence of tests: 
#>        H0: # I(0) H1: # I(0)  tstat p-value
#> Step 1          0          5 -1.052 0.05013
```

### Bootstrap FDR Controlling Tests

The function `boot_fdr()` controls for multiple testing by controlling
the false discovery rate (FDR), which is defined as the expected
proportion of false rejections relative to the total number of
rejections. This scales with the total number of tests, making it more
suitable for large `N` than the familywise error rate.

The bootstrap method for controlling FDR was introduced by Romano,
Shaikh and Wolf (2008), who showed that, unlike the classical way to
control FDR, the bootstrap is appropriate under general forms of
dependence between series. Moon and Perron (2012) applied this method to
unit root testing; it is essentially their method which is implemented
in `boot_fdr()` though again with the option to change the bootstrap
used (their suggestion was MBB). The arguments to be set are the same as
for the other multivariate unit root tests, with the exception of
`FDR_level` wihch controls the FDR level. As BSQT, the method only
report those tests until no rejection occurs.

We illustrate it here with the final available bootstrap method, the
block wild bootstrap of Shao (2011) and Smeekes and Urbain (2014a).

``` r
N <- ncol(MacroTS)
boot_fdr(MacroTS[, 1:10], FDR_level = 0.1, bootstrap = "BWB", B = 399, do_parallel = FALSE)
#> Progress: |------------------| 
#>           ********************
#> 
#>  BWB bootstrap union test with false discovery rate control
#> 
#> data: MacroTS[, 1:10]
#> null hypothesis: Series has a unit root
#> alternative hypothesis: Series is stationary
#> 
#> Sequence of tests: 
#>         tstat critical value
#> GDP_DE -1.077         -1.581
```

## Determining Order of Integration

Generally the unit root tests above would only be used as a single step
in a larger algorithm to determine the orders of integration of the time
series in the dataset. In particular, many economic datasets contain
variables that have order of integration 2, and would so need to be
differenced twice to eliminate all trends. A standard unit root test
cannot determine this however. For this purpose, we add the function
`order_integration()` which performs a sequence of unit root tests to
determine the orders of each time series.

### How does it work

Starting from a maximum order (by default equal to 2), it differences
the data time until there can be at most one unit root. If the test is
not rejected for a particular series, we know this series if of order .
The series for which we do reject are integrated once (such that they
are differenced times from their original level), and the test is
repeated. By doing so until we have classified all series, we obtain a
full specification of the orders of all time series.

### Implementation

The function allows us to choose which unit root test we want to use.
Here we take `boot_fdr`. We don’t only get the orders out, but also the
appropriately differenced data.

``` r
out_orders <- order_integration(MacroTS[, 11:15], method = "boot_fdr", B = 399, 
                                do_parallel = FALSE)
#> Progress: |------------------| 
#>           ********************
#> Progress: |------------------| 
#>           ********************
# Orders
out_orders$order_int
#> HICP_BE HICP_DE HICP_FR HICP_NL HICP_UK 
#>       0       0       1       0       1
# Differenced data
stationary_data <- out_orders$diff_data
```

To achieve the differencing, `order_integration()` uses the function
`diff_mult()` which is also available as stand-alone function in the
package. Finally, a function is provided to plot the found orders (not
run):

``` r
plot_order_integration(out_orders)
```

## References

- Cavaliere, G., Phillips, P.C.B., Smeekes, S., and Taylor, A.M.R.
  (2015). Lag length selection for unit root tests in the presence of
  nonstationary volatility. *Econometric Reviews*, 34(4), 512-536.
- Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996). Efficient tests
  for an autoregressive unit root. *Econometrica*, 64(4), 813-836.
- Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive
  wild bootstrap inference for nonparametric trends. *Journal of
  Econometrics*, 214(1), 81-109.
- Moon, H.R. and Perron, B. (2012). Beyond panel unit root tests: Using
  multiple testing to determine the non stationarity properties of
  individual series in a panel. *Journal of Econometrics*, 169(1),
  29-33.
- Ng, S. and Perron, P. (2001). Lag Length Selection and the
  Construction of Unit Root Tests with Good Size and Power.
  *Econometrica*, 69(6), 1519-1554,
- Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008). Bootstrap unit root
  tests: Comparison and extensions. *Journal of Time Series Analysis*,
  29(1), 371-401.
- Palm, F. C., Smeekes, S., and Urbain, J.-.P. (2011). Cross-sectional
  dependence robust block bootstrap panel unit root tests. *Journal of
  Econometrics*, 163(1), 85-104.
- Paparoditis, E. and Politis, D.N. (2003). Residual‐based block
  bootstrap for unit root testing. *Econometrica*, 71(3), 813-855.
- Perron, P. and Qu, Z. (2008). A simple modification to improve the
  finite sample properties of Ng and Perron’s unit root tests. *Economic
  Letters*, 94(1), 12-19.
- Rho, Y. and Shao, X. (2019). Bootstrap-assisted unit root testing with
  piecewise locally stationary errors. *Econometric Theory*, 35(1),
  142-166.
- Romano, J.P., Shaikh, A.M., and Wolf, M. (2008). Control of the false
  discovery rate under dependence using the bootstrap and subsampling.
  *Test*, 17(3), 417.
- Romano, J. P. and Wolf, M. (2005). Stepwise multiple testing as
  formalized data snooping. *Econometrica*, 73(4), 1237-1282.
- Shao, X. (2010). The dependent wild bootstrap. *Journal of the
  American Statistical Association*, 105(489), 218-235.
- Shao, X. (2011). A bootstrap-assisted spectral test of white noise
  under unknown dependence. *Journal of Econometrics*, 162, 213-224.
- Smeekes (2013). Detrending bootstrap unit root tests. *Econometric
  Reviews*, 32(8), 869-891.
- Smeekes, S. (2015). Bootstrap sequential tests to determine the order
  of integration of individual units in a time series panel. *Journal of
  Time Series Analysis*, 36(3), 398-415.
- Smeekes, S. and Taylor, A.M.R. (2012). Bootstrap union tests for unit
  roots in the presence of nonstationary volatility. *Econometric
  Theory*, 28(2), 422-456.
- Smeekes, S. and Urbain, J.-P. (2014a). A multivariate invariance
  principle for modified wild bootstrap methods with an application to
  unit root testing. GSBE Research Memorandum No. RM/14/008, Maastricht
  University.
- Smeekes, S. and Urbain, J.-P. (2014b). On the applicability of the
  sieve bootstrap in time series panels. *Oxford Bulletin of Economics
  and Statistics*, 76(1), 139-151.
