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
