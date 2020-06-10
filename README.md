
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bootUR

`bootUR` implements several bootstrap tests for unit, both for single
time series and for (potentially) large systems of time series.

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

We illustrate the bootstrap ADF test here on Dutch GDP, with the sieve
bootstrap (`boot = SB`) to as the specification used by Palm, Smeekes
and Urbain (2008) and Smeekes (2013). We set only 399 bootstrap
replications (`B = 399`) to prevent the code from running too long. We
add an intercept and a trend (`dc = 2`), and compare OLS with QD(GLS)
detrending. The option `verbose = TRUE` prints easy to read output on
the console.

``` r
GDP_NL <- MacroTS[, 4]
adf_out <- boot_df(GDP_NL, B = 399, boot = "SB", dc = 2, detr = c("OLS", "QD"), verbose = TRUE)
#> Bootstrap DF Test with SB bootstrap method.
#> ----------------------------------------
#> Type of unit root test performed: detr = OLS, dc = intercept and trend
#> test statistic        p-value 
#>     -2.5152854      0.1553885 
#> ----------------------------------------
#> Type of unit root test performed: detr = QD, dc = intercept and trend
#> test statistic        p-value 
#>     -1.5965001      0.4461153
```

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
#>     -0.6923184      0.6140351
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
it with the dependent wild bootstrap (DWB) of Shao (2010).

By default the union test is used for each series (`union = TRUE`), if
this is set to `FALSE` the deterministic components and detrending
methods can be specified as in the univariate Dickey-Fuller test.

Although the sieve bootstrap method `"SB"` and `"SWB"` can be used
(historically they have been popular among practitioners), Smeekes and
Urbain (2014) show that these are not suited to capture general forms of
dependence across units. The code will give a warning to recommend using
a different bootstrap method.

``` r
# This will give an error!
# panel_out <- paneltest(MacroTS, boot = "MBB", B = 399, verbose = TRUE)
panel_out <- paneltest(MacroTS, boot = "DWB", B = 399, verbose = TRUE)
#> Panel Bootstrap Group-Mean Union Test
#> The null hypothesis that all series have a unit root, is rejected at a significance level of 0.05.
#>      test statistic    p-value
#> [1,]     -0.8525506 0.03258145
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
#>        test statistic    p-value
#> GDP_BE      -2.792169 0.22556391
#> GDP_DE      -2.774320 0.09273183
#> GDP_FR      -2.048760 0.52882206
#> GDP_NL      -2.515285 0.17543860
#> GDP_UK      -2.449065 0.28320802
```

## Sequential Tests for Multiple Time Series

Smeekes (2015). Bootstrap Sequential Quantile Test (BSQT). The other
options are the same as for
`paneltest`.

``` r
BSQT_out <- BSQTtest(MacroTS[, 16:20], boot = "AWB", B = 399, verbose = TRUE)
#> There are 0 stationary time series.
#> Details of the BSQT ssquential tests:
#>        Unit H0 Unit H1 Test statistic   p-value
#> Step 1       0       1      -0.924728 0.6641604
```

## References

  - Elliott, Rothenberg and Stock (1996)
  - Palm, Smeekes and Urbain (2011)
  - Shao (2010)
  - Smeekes (2013)
  - Smeekes and Urbain (2014)
  - Smeekes and Taylor (2012)
