## Version 0.1.0.901

### Bug Fixes
* Fixed an incorrect warning in check_inputs() function which was given
whenever NAs were present and S(W)B was used, even if the dataset remained
balanced.
* Fixed bug in check_inputs() function that also set bootstrap to be
performed individually in the presence of NAs even if S(W)B was not used.
* Fixed bug/incorrect version in BSQTtest() that incorrectly added another 0
to argument q even if one was present already.

## Version 0.1.0
First release; package available on CRAN.
