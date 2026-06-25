# bootUR: Bootstrap Unit Root Tests

The R package `bootUR` implements several bootstrap tests for unit
roots, both for single time series and for (potentially) large systems
of time series.

## Installation and Loading

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

After installation, the package can be loaded in the standard way:

``` r

library(bootUR)
```

## Functionality

A quick overview of the package functionality is provided in the
vignette
[`vignette("bootUR")`](https://smeekes.github.io/bootUR/articles/bootUR.md).

A further investigation of the functionalities is provided in the
*Journal of Statistical Software* article [bootUR: An R Package for
Bootstrap Unit Root Tests](https://doi.org/10.18637/jss.v106.i12).
