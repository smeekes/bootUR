# Auxiliary Function (not accessible to users) to create all bootstrap statistics used to perform the unit root tests.

Auxiliary Function (not accessible to users) to create all bootstrap
statistics used to perform the unit root tests.

## Usage

``` r
do_tests_and_bootstrap(data, boot_sqt_test, boot_ur_test, level, bootstrap, B,
  block_length, ar_AWB, union, min_lag, max_lag, criterion, deterministics,
  detrend, criterion_scale, steps, h_rs, show_progress, do_parallel, cores,
  data_name)
```

## Arguments

- data:

  A \\T\\-dimensional vector or a (\\T\\ x \\N\\)-matrix of \\N\\ time
  series with \\T\\ observations to be tested for unit roots. Data may
  also be in a time series format (e.g. `ts`, `zoo` or `xts`), or a data
  frame, as long as each column represents a single time series.

- boot_sqt_test:

  Logical indicator whether or not to perform the Bootstrap Quantile
  Test (`TRUE`) or not (`FALSE`).

- boot_ur_test:

  Logical indicator whether or not to perform the individual ADF Tests
  (`TRUE`) or not (`FALSE`).

- level:

  Desired significance level of the unit root test.

- bootstrap:

  String for bootstrap method to be used. Options are

  `"MBB"`

  :   Moving blocks bootstrap;

  `"BWB"`

  :   Block wild bootstrap;

  `"DWB"`

  :   Dependent wild bootstrap;

  `"AWB"`

  :   Autoregressive wild bootstrap;

  `"SB"`

  :   Sieve bootstrap;

  `"SWB"`

  :   Sieve wild bootstrap.

- B:

  Number of bootstrap replications. Default is 1999.

- block_length:

  Desired 'block length' in the bootstrap. For the MBB, BWB and DWB
  bootstrap, this is a genuine block length. For the AWB bootstrap, the
  block length is transformed into an autoregressive parameter via the
  formula \\0.01^(1/block_length)\\; this can be overwritten by setting
  `ar_AWB` directly. If NULL, sets the block length as a function of the
  time series length T, via the rule \\block_length = 1.75 T^(1/3)\\.

- ar_AWB:

  Autoregressive parameter used in the AWB bootstrap method
  (`bootstrap = "AWB"`). Can be used to set the parameter directly
  rather than via the default link to the block length.

- union:

  Logical indicator whether or not to use bootstrap union tests (`TRUE`)
  or not (`FALSE`).

- min_lag:

  Minimum lag length in the augmented Dickey-Fuller regression.

- max_lag:

  Maximum lag length in the augmented Dickey-Fuller regression.

- criterion:

  String for information criterion used to select the lag length in the
  augmented Dickey-Fuller regression. Options are: `"AIC"`, `"BIC"`,
  `"MAIC"`, `"MBIC`.

- deterministics:

  Numeric vector indicating the deterministic specification. Only
  relevant if `union = FALSE`. Options are (combinations of): `0` - no
  deterministics; `1` - intercept only; `2` - intercept and trend. If
  `union = FALSE` and `deterministics = NULL`, an intercept is added.

- detrend:

  String vector indicating the type of detrending to be performed. Only
  relevant if `union = FALSE`. Options are: `"OLS"` and/or `"QD"`
  (typically also called GLS). If NULL, set to `"OLS"`.

- criterion_scale:

  Logical indicator whether or not to use the rescaled information
  criteria (`TRUE`) or not (`FALSE`).

- steps:

  Numeric vector of quantiles to be tested. Default is to test each unit
  sequentially.

- h_rs:

  Bandwidth used in rescaled information criteria.

- show_progress:

  Logical indicator whether a bootstrap progress update should be
  printed to the console. Default is FALSE.

- do_parallel:

  Logical indicator whether bootstrap loop should be executed in
  parallel. Default is TRUE.

- cores:

  The number of cores to be used in the parallel loops. Default is to
  use all but one.

- data_name:

  Optional name for the data, to be used in the output. The default uses
  the name of the 'data' argument.

## See also

[`boot_ur`](https://smeekes.github.io/bootUR/reference/boot_ur.md),
[`boot_sqt`](https://smeekes.github.io/bootUR/reference/boot_sqt.md),
[`boot_fdr`](https://smeekes.github.io/bootUR/reference/boot_fdr.md)
