## Version 1.0.4

### Miscellaneous
* Added RcppParallel cxxflags to makevars, needed for compatibility with Windows ARM64
* GitHub Actions added
* New package website

## Version 1.0.3

### Bug Fixes
* Removed the unused argument 'h_rs' from the documentation of the function chech_inputs() (CRAN note).

## Version 1.0.2

### Bug Fixes
* Fixed the CRAN check Additional Issues LTO warnings by avoiding the use of arma::cube in C++ code.

## Version 1.0.0

### New Functionality
* Updated citation and references to incorporate Journal of Statistical Software publication.

### Bug Fixes
* Fixed the error in boot_panel() when not performing union.

## Version 0.5.0

### New Functionality
* Output slots added that store the details of the tests such as individual test 
statistics and selected lag lengths.
* Removed the deprecated test names.

### Bug Fixes
* Fixed the error in order_of_integration() with boot_ur() test applied to single time series.

## Version 0.4.2

### New Functionality
* Output slot added that stores te specifications of the tests.

### Bug Fixes
* Visually differentiated test print output by making the name reflect the most 
important specifications.
* Forced normalisation of the signs of the eigenvectors calculated for DWB to
ensure reproducability of results across OSes.

## Version 0.4.1

### New Functionality
* Parallel computing via the RcppParallel package replaces OpenMP.

### Bug Fixes
* Fixed order_integration() to allow for all possible tests.

## Version 0.3.0

### New Functionality
* Added new functions  boot_adf(), boot_ur(), boot_frd(), boot_sqt() and 
boot_panel() that replace the old functions boot_df(),  iADFtest(), bFDRtest(),
BQSTtest() and paneltest(). The old functions are still available in the package
as deprecated functions.
* For additional clarity, new arguments names are used in the new functions.
* New adf() function added which implements the standard, non-bootstrap, augmented Dickey-Fuller (ADF) test 


## Version 0.2.0

### New Functionality
* Added parallel loops using OpenMP.
* Added functions order_integration() to determine the order of integration of 
each time series in a dataset. In addition mult_diff() is added to difference
the time series accordingly to eliminate stochastic trends, and 
plot_order_integration() is added to plot the found orders of integration.
* Added function plot_missing_values() to give a visual representation of the
pattern of missing values in the data.

### Bug Fixes
* Fixed an incorrect warning in check_inputs() function which was given
whenever NAs were present and S(W)B was used, even if the dataset remained
balanced.
* Fixed bug in check_inputs() function that also set bootstrap to be
performed individually in the presence of NAs even if S(W)B was not used.
* Fixed bug/incorrect warning in BSQTtest() that incorrectly added another 0
to argument q even if one was present already.
* Fixed possibility of taking a too fine grid for argument q in BSQTtest()
leading to duplicates. Duplicates are now removed with a warning given.
* Fixed bug in SB where the last l observations were excluded from resampling.
* Fixed bug in S(W)B where too short vector with ADF residuals was used 
(supplemented with zeros) in bootstrap.

## Version 0.1.0
First release; package available on CRAN.
