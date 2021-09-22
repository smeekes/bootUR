#' Individual Unit Root Tests without multiple testing control
#' @description This function performs bootstrap unit root tests on each time series individually.
#' @param data A \eqn{T}-dimensional vector or a (\eqn{T} x \eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations to be tested for unit roots. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}), or a data frame, as long as each column represents a single time series.
#' @param level Desired significance level of the unit root test. Default is 0.05.
#' @param bootstrap String for bootstrap method to be used. Options are
#' \describe{
#' \item{\code{"MBB"}}{Moving blocks bootstrap (Paparoditis and Politis, 2003; Palm, Smeekes and Urbain, 2011);}
#' \item{\code{"BWB"}}{Block wild bootstrap (Shao, 2011; Smeekes and Urbain, 2014a);}
#' \item{\code{"DWB"}}{Dependent wild bootstrap (Shao, 2010; Smeekes and Urbain, 2014a; Rho and Shao, 2019);}
#' \item{\code{"AWB"}}{Autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2020), this is the default;}
#' \item{\code{"SB"}}{Sieve bootstrap (Chang and Park, 2003; Palm, Smeekes and Urbain, 2008; Smeekes, 2013);}
#' \item{\code{"SWB"}}{Sieve wild boostrap (Cavaliere and Taylor, 2009; Smeekes and Taylor, 2012).}
#' }
#' @param B Number of bootstrap replications. Default is 1999.
#' @param block_length Desired 'block length' in the bootstrap. For the MBB, BWB and DWB boostrap, this is a genuine block length. For the AWB boostrap, the block length is transformed into an autoregressive parameter via the formula \eqn{0.01^(1/block_length)} as in Smeekes and Urbain (2014a); this can be overwritten by setting \code{ar_AWB} directly. Default sets the block length as a function of the time series length T, via the rule \eqn{block_length = 1.75 T^(1/3)} of Palm, Smeekes and Urbain (2011).
#' @param ar_AWB Autoregressive parameter used in the AWB bootstrap method (\code{bootstrap = "AWB"}). Can be used to set the parameter directly rather than via the default link to the block length.
#' @param union Logical indicator whether or not to use bootstrap union tests (\code{TRUE}) or not (\code{FALSE}), see Smeekes and Taylor (2012). Default is \code{TRUE}.
#' @param min_lag Minimum lag length in the augmented Dickey-Fuller regression. Default is 0.
#' @param max_lag Maximum lag length in the augmented Dickey-Fuller regression. Default uses the sample size-based rule \eqn{12(T/100)^{1/4}}.
#' @param criterion String for information criterion used to select the lag length in the augmented Dickey-Fuller regression. Options are: \code{"AIC"}, \code{"BIC"}, \code{"MAIC"}, \code{"MBIC"}. Default is \code{"MAIC"} (Ng and Perron, 2001).
#' @param deterministics String indicating the deterministic specification. Only relevant if \code{union = FALSE}. Options are
#'
#' \verb{"none":} no deterministics;
#'
#' \verb{"intercept":} intercept only;
#'
#' \verb{"trend":} intercept and trend.
#'
#' If \code{union = FALSE}, the default is adding an intercept (a warning is given).
#' @param detrend String indicating the type of detrending to be performed. Only relevant if \code{union = FALSE}. Options are: \code{"OLS"} or \code{"QD"} (typically also called GLS, see Elliott, Rothenberg and Stock, 1996). The default is \code{"OLS"}.
#' @param criterion_scale Logical indicator whether or not to use the rescaled information criteria of Cavaliere et al. (2015) (\code{TRUE}) or not (\code{FALSE}). Default is \code{TRUE}.
#' @param show_progress Logical indicator whether a bootstrap progress update should be printed to the console. Default is FALSE.
#' @param do_parallel Logical indicator whether bootstrap loop should be executed in parallel. Parallel computing is only available if OpenMP can be used, if not this option is ignored. Default is FALSE.
#' @param cores The number of cores to be used in the parallel loops. Default is to use all but one.
#' @details The options encompass many test proposed in the literature. \code{detrend = "OLS"} gives the standard augmented Dickey-Fuller test, while \code{detrend = "QD"} provides the DF-GLS test of Elliott, Rothenberg and Stock (1996). The bootstrap algorithm is always based on a residual bootstrap (under the alternative) to obtain residuals rather than a difference-based bootstrap (under the null), see e.g. Palm, Smeekes and Urbain (2008).
#'
#' Lag length selection is done automatically in the ADF regression with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. For very short time series (fewer than 50 time points) the maximum lag length is adjusted downward to avoid potential multicollinearity issues in the bootstrap. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `min_lag` and maximum lag length `max_lag` for the selection algorithm equal to the desired lag length.
#'
#' @return A list with N components, one for each variable, where each element of the list returns an object of class \code{htest} containing
#' \item{\code{method}}{The name of the hypothesis test. For boot_ur these are ADF tests (for \code{union = FALSE}) or Union tests (for \code{union = TRUE}) on each individual series (no multiple testing correction).;}
#' \item{\code{data.name}}{The name of the variable on which the ADF test is performed.;}
#' \item{\code{null.value}}{The value of the (gamma) parameter of the lagged dependent variable in the ADF regression under the null hypothesis. Under the null, the series has a unit root. Testing the null of a unit root then boils down to testing the significance of the gamma parameter;}
#' \item{\code{alternative}}{A character string specifying the direction of the alternative hypothesis relative to the null value. The alternative postulates that the series is stationary;}
#' \item{\code{estimate}}{The estimated value of the (gamma) parameter of the lagged dependent variable in the ADF regression. Note that for the union test (\code{union = TRUE}), this estimate is not defined, hence NA is returned.;}
#' \item{\code{statistic}}{The value of the test statistic of the unit root test.;}
#' \item{\code{p.value}}{P-value of the unit root test.}
#' @section Warnings:
#' The function may give the following warnings.
#' \describe{
#' \item{\code{Warning: Missing values cause resampling bootstrap to be executed for each time series individually.}}{If the time series in \code{data} have different starting and end points (and thus some series contain \code{NA} values at the beginning and/or end of the sample, the resampling-based moving block bootstrap (MBB) and sieve bootstrap (SB) cannot be used directly, as they create holes (internal missings) in the bootstrap samples. These bootstrap methods are therefore not applied jointly as usual, but individually to each series.}
#' \item{\code{Warning: Deterministic specification in argument deterministics is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting deterministic components manually therefore has no effect.}
#' \item{\code{Warning: Detrending method in argument detrend is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting detrending methods manually therefore has no effect.}
#' }
#' @references Chang, Y. and Park, J. (2003). A sieve bootstrap for the test of a unit root. \emph{Journal of Time Series Analysis}, 24(4), 379-400.
#' @references Cavaliere, G. and Taylor, A.M.R (2009). Heteroskedastic time series with a unit root. \emph{Econometric Theory}, 25, 1228–1276.
#' @references Cavaliere, G., Phillips, P.C.B., Smeekes, S., and Taylor, A.M.R.
#' (2015). Lag length selection for unit root tests in the presence of nonstationary volatility. \emph{Econometric Reviews}, 34(4), 512-536.
#' @references Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996). Efficient tests for an autoregressive unit root. \emph{Econometrica}, 64(4), 813-836.
#' @references Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild bootstrap inference for nonparametric trends. \emph{Journal of Econometrics}, 214(1), 81-109.
#' @references Ng, S. and Perron, P. (2001). Lag Length Selection and the Construction of Unit Root Tests with Good Size and Power. \emph{Econometrica}, 69(6), 1519-1554,
#' @references Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008). Bootstrap unit root tests: Comparison and extensions. \emph{Journal of Time Series Analysis}, 29(1), 371-401.
#' @references Palm, F. C., Smeekes, S., and Urbain, J.-.P. (2011). Cross-sectional dependence robust block bootstrap panel unit root tests. \emph{Journal of Econometrics}, 163(1), 85-104.
#' @references Paparoditis, E. and Politis, D.N. (2003). Residual-based block bootstrap for unit root testing. \emph{Econometrica}, 71(3), 813-855.
#' @references Perron, P. and Qu, Z. (2008). A simple modification to improve the finite sample properties of Ng and Perron's unit root tests. \emph{Economic Letters}, 94(1), 12-19.
#' @references Rho, Y. and Shao, X. (2019). Bootstrap-assisted unit root testing with piecewise locally stationary errors. \emph{Econometric Theory}, 35(1), 142-166.
#' @references Shao, X. (2010). The dependent wild bootstrap. \emph{Journal of the American Statistical Association}, 105(489), 218-235.
#' @references Shao, X. (2011). A bootstrap-assisted spectral test of white noise under unknown dependence. \emph{Journal of Econometrics}, 162, 213-224.
#' @references Smeekes, S. (2013). Detrending bootstrap unit root tests. \emph{Econometric Reviews}, 32(8), 869-891.
#' @references Smeekes, S. and Taylor, A.M.R. (2012). Bootstrap union tests for unit roots in the presence of nonstationary volatility. \emph{Econometric Theory}, 28(2), 422-456.
#' @references Smeekes, S. and Urbain, J.-P. (2014a). A multivariate invariance principle for modified wild bootstrap methods with an application to unit root testing. GSBE Research Memorandum No. RM/14/008, Maastricht University
#' @examples
#' # boot_ur on GDP_BE and GDP_DE
#' two_series_boot_ur <- boot_ur(MacroTS[, 1:2], bootstrap = "MBB", B = 399)
#' print(two_series_boot_ur)
#' @export
boot_ur <- function(data, level = 0.05, bootstrap = "AWB", B = 1999, block_length = NULL,
                    ar_AWB = NULL, union = TRUE, min_lag = 0, max_lag = NULL,
                    criterion = "MAIC", deterministics = NULL, detrend = NULL, criterion_scale = TRUE,
                    show_progress = FALSE, do_parallel = FALSE, cores = NULL){
  
  inputs <- do_tests_and_bootstrap(data = data, boot_sqt_test = FALSE, boot_ur_test = TRUE,
                                   level = level, bootstrap = bootstrap, B = B,
                                   block_length = block_length, ar_AWB = ar_AWB,
                                   union = union, min_lag = min_lag, max_lag = max_lag,
                                   criterion = criterion, deterministics = deterministics,
                                   detrend = detrend, criterion_scale = criterion_scale,
                                   steps = NULL, h_rs = 0.1, show_progress = show_progress,
                                   do_parallel = do_parallel, cores = cores)

  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(data))
  }

  # Results
  if (union) { # Union test
    iADFout <- iADF_cpp(test_i = inputs$test_stats, t_star = inputs$test_stats_star, level = inputs$level)
    iADFout <- cbind(rep(NA, nrow(iADFout)), iADFout)
    if (NCOL(data) > 1) {
      colnames(iADFout) <- c("gamma", "tstat", "pvalue")
    }
    # Parameter estimates, tstats and p-values. Note: Parameter Estimates not defined for union test

    if (NCOL(data) == 1){ # boot_union
      method_name <- "Bootstrap Union test on a single time series"
    } else {
      method_name <- "Bootstrap Union tests on each individual series (no multiple testing correction)"
    }

  } else { # No union test
    iADFout <- iADF_cpp(test_i = matrix(inputs$tests_i[1, ], nrow = 1),
                        t_star = matrix(inputs$t_star[ , 1, ], nrow = B), level = inputs$level)
    iADFout <- cbind(t(inputs$param_i), iADFout) # Parameter estimates, tstats and p-values

    if (NCOL(data) == 1){ # boot_adf
      method_name <- "Bootstrap ADF test on a single time series"
    } else {
      method_name <- "Bootstrap ADF tests on each individual series (no multiple testing correction)"
    }
  }

  rej_H0 <- (iADFout[, 3] < level)
  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(data))
  }

  if (NCOL(data) > 1) {
    boot_ur_output <- list(method = method_name, data.name = "data", details = iADFout,
                       alternative = "less", null.value =  c("gamma" = 0), rejections = rej_H0,
                       estimate = iADFout[, 1], statistic = iADFout[, 2], p.value = iADFout[, 3])
    class(boot_ur_output) <- "mult_htest"
  } else {
    iADFtstat <- iADFout[1, 2]
    attr(iADFtstat, "names") <- "tstat"
    boot_ur_output <- list(method = method_name, data.name = var_names, null.value = c("gamma" = 0),
                         alternative = "less", estimate = iADFout[1, 1], statistic = iADFtstat, p.value = iADFout[1, 3])
    class(boot_ur_output) <- "htest"

  }
  
  return(boot_ur_output)
}


#' Bootstrap augmented Dickey-Fuller Unit Root Test
#' @description This function performs a standard augmented Dickey-Fuller bootstrap unit root test on a single time series.
#' @inheritParams boot_ur
#' @param data A \eqn{T}-dimensional vector to be tested for unit roots. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}), or a data frame.
#' @param bootstrap String for bootstrap method to be used. Options are
#' \describe{
#' \item{\code{"MBB"}}{Moving blocks bootstrap (Paparoditis and Politis, 2003);}
#' \item{\code{"BWB"}}{Block wild bootstrap (Shao, 2011);}
#' \item{\code{"DWB"}}{Dependent wild bootstrap (Shao, 2010; Rho and Shao, 2019);}
#' \item{\code{"AWB"}}{Autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2020), this is the default;}
#' \item{\code{"SB"}}{Sieve bootstrap (Chang and Park, 2003; Palm, Smeekes and Urbain, 2008; Smeekes, 2013);}
#' \item{\code{"SWB"}}{Sieve wild boostrap (Cavaliere and Taylor, 2009; Smeekes and Taylor, 2012).}
#' }
#' @details The options encompass many test proposed in the literature. \code{detrend = "OLS"} gives the standard augmented Dickey-Fuller test, while \code{detrend = "QD"} provides the DF-GLS test of Elliott, Rothenberg and Stock (1996). The bootstrap algorithm is always based on a residual bootstrap (under the alternative) to obtain residuals rather than a difference-based bootstrap (under the null), see e.g. Palm, Smeekes and Urbain (2008).
#'
#' Lag length selection is done automatically in the ADF regression with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. For very short time series (fewer than 50 time points) the maximum lag length is adjusted downward to avoid potential multicollinearity issues in the bootstrap. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `min_lag` and maximum lag length `max_lag` for the selection algorithm equal to the desired lag length.
#' @export
#' @return An object of class \code{htest} containing
#' \item{\code{method}}{The name of the hypothesis test. For boot_adf this is the ADF test on a single time series.;}
#' \item{\code{data.name}}{The name of the variable on which the ADF test is performed.;}
#' \item{\code{null.value}}{The value of the (gamma) parameter of the lagged dependent variable in the ADF regression under the null hypothesis. Under the null, the series has a unit root. Testing the null of a unit root then boils down to testing the significance of the gamma parameter;}
#' \item{\code{alternative}}{A character string specifying the direction of the alternative hypothesis relative to the null value. The alternative postulates that the series is stationary;}
#' \item{\code{estimate}}{The estimated value of the (gamma) parameter of the lagged dependent variable in the ADF regression.;}
#' \item{\code{statistic}}{The value of the test statistic of the unit root test.;}
#' \item{\code{p.value}}{P-value of the unit root test.}
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Multiple time series not allowed. Switch to a multivariate method such as boot_ur, or change argument data to a univariate time series.}}{The function is a simple wrapper around \code{\link{boot_ur}} to facilitate use for single time series. It does not support multiple time series, as \code{\link{boot_ur}} is specifically suited for that.}
#' }
#' @references Chang, Y. and Park, J. (2003). A sieve bootstrap for the test of a unit root. \emph{Journal of Time Series Analysis}, 24(4), 379-400.
#' @references Cavaliere, G. and Taylor, A.M.R (2009). Heteroskedastic time series with a unit root. \emph{Econometric Theory}, 25, 1228–1276.
#' @references Cavaliere, G., Phillips, P.C.B., Smeekes, S., and Taylor, A.M.R. (2015). Lag length selection for unit root tests in the presence of nonstationary volatility. \emph{Econometric Reviews}, 34(4), 512-536.
#' @references Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996). Efficient tests for an autoregressive unit root. \emph{Econometrica}, 64(4), 813-836.
#' @references Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild bootstrap inference for nonparametric trends. \emph{Journal of Econometrics}, 214(1), 81-109.
#' @references Ng, S. and Perron, P. (2001). Lag Length Selection and the Construction of Unit Root Tests with Good Size and Power. \emph{Econometrica}, 69(6), 1519-1554,
#' @references Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008). Bootstrap unit root tests: Comparison and extensions. \emph{Journal of Time Series Analysis}, 29(1), 371-401.
#' @references Paparoditis, E. and Politis, D.N. (2003). Residual-based block bootstrap for unit root testing. \emph{Econometrica}, 71(3), 813-855.
#' @references Perron, P. and Qu, Z. (2008). A simple modification to improve the finite sample properties of Ng and Perron's unit root tests. \emph{Economic Letters}, 94(1), 12-19.
#' @references Rho, Y. and Shao, X. (2019). Bootstrap-assisted unit root testing with piecewise locally stationary errors. \emph{Econometric Theory}, 35(1), 142-166.
#' @references Smeekes, S. (2013). Detrending bootstrap unit root tests. \emph{Econometric Reviews}, 32(8), 869-891.
#' @references Shao, X. (2010). The dependent wild bootstrap. \emph{Journal of the American Statistical Association}, 105(489), 218-235.
#' @references Shao, X. (2011). A bootstrap-assisted spectral test of white noise under unknown dependence. \emph{Journal of Econometrics}, 162, 213-224.
#' @references Smeekes, S. and Taylor, A.M.R. (2012). Bootstrap union tests for unit roots in the presence of nonstationary volatility. \emph{Econometric Theory}, 28(2), 422-456.
#' @references Smeekes, S. and Urbain, J.-P. (2014a). A multivariate invariance principle for modified wild bootstrap methods with an application to unit root testing. GSBE Research Memorandum No. RM/14/008, Maastricht University
#' @seealso \code{\link{boot_ur}}
#' @examples
#' # boot_adf on GDP_BE
#' GDP_BE_adf <- boot_adf(MacroTS[, 1], B = 399, deterministics = "trend",
#' detrend = "OLS")
#' print(GDP_BE_adf)
boot_adf <- function(data, level = 0.05, bootstrap = "AWB", B = 1999, block_length = NULL, ar_AWB = NULL,
                     min_lag = 0, max_lag = NULL, criterion = "MAIC", deterministics = "intercept",
                     detrend = "OLS", criterion_scale = TRUE, show_progress = FALSE,
                     do_parallel = FALSE, cores = NULL){

  if (NCOL(data) > 1) {
    stop("Multiple time series not allowed. Switch to a multivariate method such as boot_ur,
         or change argument data to a univariate time series.")
  }

  out <- boot_ur(data = matrix(data, ncol = 1), level = level, bootstrap = bootstrap, B = B,
                 block_length = block_length, ar_AWB = ar_AWB, union = FALSE, min_lag = min_lag,
                 max_lag = max_lag, criterion = criterion, deterministics = deterministics,
                 detrend = detrend, criterion_scale = criterion_scale,
                 show_progress = show_progress, do_parallel = do_parallel, cores = cores)

  return(out)
}


#' Bootstrap Union Test for Unit Roots
#' @description Performs bootstrap unit root test based on the union of rejections of 4 tests with different number of deterministic components and different type of detrending (Harvey, Leybourne and Taylor, 2012; Smeekes and Taylor, 2012).
#' @inheritParams boot_ur
#' @param data A \eqn{T}-dimensional vector to be tested for unit roots. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}), or a data frame.
#' @param bootstrap String for bootstrap method to be used. Options are
#' \describe{
#' \item{\code{"MBB"}}{Moving blocks bootstrap (Paparoditis and Politis, 2003);}
#' \item{\code{"BWB"}}{Block wild bootstrap (Shao, 2011);}
#' \item{\code{"DWB"}}{Dependent wild bootstrap (Shao, 2010; Rho and Shao, 2019);}
#' \item{\code{"AWB"}}{Autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2020), this is the default;}
#' \item{\code{"SB" }}{Sieve bootstrap (Palm, Smeekes and Urbain, 2008);}
#' \item{\code{"SWB"}}{Sieve wild boostrap (Cavaliere and Taylor, 2009; Smeekes and Taylor, 2012).}
#' }
#' @details The union is taken over the combination of tests with intercept only and intercept plus trend, coupled with OLS detrending and QD detrending, as in Harvey, Leybourne and Taylor (2012) and Smeekes an Taylor (2012). The bootstrap algorithm is always based on a residual bootstrap (under the alternative) to obtain residuals rather than a difference-based bootstrap (under the null), see e.g. Palm, Smeekes and Urbain (2008).
#'
#' Lag length selection is done automatically in the ADF regressions with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `min_lag` and maximum lag length `max_lag` for the selection algorithm equal to the desired lag length.
#' @export
#' @return An object of class \code{htest} containing
#' \item{\code{method}}{The name of the hypothesis test. For boot_adf this is the Union test on a single time series.;}
#' \item{\code{data.name}}{The name of the variable on which the union test is performed.;}
#' \item{\code{null.value}}{The value of the (gamma) parameter of the lagged dependent variable in the ADF regression under the null hypothesis. Under the null, the series has a unit root. Testing the null of a unit root then boils down to testing the significance of the gamma parameter;}
#' \item{\code{alternative}}{A character string specifying the direction of the alternative hypothesis relative to the null value. The alternative postulates that the series is stationary;}
#' \item{\code{estimate}}{For the union test, the estimated value of the (gamma) parameter of the lagged dependent variable in the ADF regression is not defined.;}
#' \item{\code{statistic}}{The value of the test statistic of the unit root test.;}
#' \item{\code{p.value}}{P-value of the unit root test.}
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Multiple time series not allowed. Switch to a multivariate method such as boot_ur, or change argument data to a univariate time series.}}{The function is a simple wrapper around \code{\link{boot_ur}} to facilitate use for single time series. It does not support multiple time series, as \code{\link{boot_ur}} is specifically suited for that.}
#' }
#' @references Chang, Y. and Park, J. (2003). A sieve bootstrap for the test of a unit root. \emph{Journal of Time Series Analysis}, 24(4), 379-400.
#' @references Cavaliere, G. and Taylor, A.M.R (2009). Heteroskedastic time series with a unit root. \emph{Econometric Theory}, 25, 1228–1276.
#' @references Cavaliere, G., Phillips, P.C.B., Smeekes, S., and Taylor, A.M.R.
#' (2015). Lag length selection for unit root tests in the presence of nonstationary volatility. \emph{Econometric Reviews}, 34(4), 512-536.
#' @references Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild bootstrap inference for nonparametric trends. \emph{Journal of Econometrics}, 214(1), 81-109.
#' @references Harvey, D.I., Leybourne, S.J., and Taylor, A.M.R. (2012). Testing for unit roots in the presence of uncertainty over both the trend and initial condition. \emph{Journal of Econometrics}, 169(2), 188-195.
#' @references Ng, S. and Perron, P. (2001). Lag Length Selection and the Construction of Unit Root Tests with Good Size and Power. \emph{Econometrica}, 69(6), 1519-1554,
#' @references Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008). Bootstrap unit root tests: Comparison and extensions. \emph{Journal of Time Series Analysis}, 29(1), 371-401.
#' @references Paparoditis, E. and Politis, D.N. (2003). Residual-based block bootstrap for unit root testing. \emph{Econometrica}, 71(3), 813-855.
#' @references Perron, P. and Qu, Z. (2008). A simple modification to improve the finite sample properties of Ng and Perron's unit root tests. \emph{Economic Letters}, 94(1), 12-19.
#' @references Rho, Y. and Shao, X. (2019). Bootstrap-assisted unit root testing with piecewise locally stationary errors. \emph{Econometric Theory}, 35(1), 142-166.
#' @references Shao, X. (2010). The dependent wild bootstrap. \emph{Journal of the American Statistical Association}, 105(489), 218-235.
#' @references Shao, X. (2011). A bootstrap-assisted spectral test of white noise under unknown dependence. \emph{Journal of Econometrics}, 162, 213-224.
#' @references Smeekes, S. (2013). Detrending bootstrap unit root tests. \emph{Econometric Reviews}, 32(8), 869-891.
#' @references Smeekes, S. and Taylor, A.M.R. (2012). Bootstrap union tests for unit roots in the presence of nonstationary volatility. \emph{Econometric Theory}, 28(2), 422-456.
#' @references Smeekes, S. and Urbain, J.-P. (2014a). A multivariate invariance principle for modified wild bootstrap methods with an application to unit root testing. GSBE Research Memorandum No. RM/14/008, Maastricht University
#' @seealso \code{\link{boot_ur}}
#' @examples
#' # boot_union on GDP_BE
#' GDP_BE_df <- boot_union(MacroTS[, 1], B = 399)
#' print(GDP_BE_df)
boot_union <- function(data, level = 0.05, bootstrap = "AWB", B = 1999, block_length = NULL,
                       ar_AWB = NULL, min_lag = 0, max_lag = NULL, criterion = "MAIC",
                       criterion_scale = TRUE, show_progress = FALSE, do_parallel = FALSE,
                       cores = NULL){

  if (NCOL(data) > 1) {
    stop("Multiple time series not allowed. Switch to a multivariate method such as boot_ur,
         or change argument data to a univariate time series.")
  }
  out <- boot_ur(data = data, level = level, bootstrap = bootstrap, B = B,
                 block_length = block_length, ar_AWB = ar_AWB, union = TRUE, min_lag = min_lag,
                 max_lag = max_lag, criterion = criterion, criterion_scale = criterion_scale,
                 show_progress = show_progress, do_parallel = do_parallel,
                 cores = cores)

  return(out)
}

#' Bootstrap Unit Root Tests with False Discovery Rate control
#' @description Controls for multiple testing by controlling the false discovery rate (FDR), see Moon and Perron (2012) and Romano, Shaikh and Wolf (2008).
#' @inheritParams boot_ur
#' @param level Desired False Discovery Rate level of the unit root tests. Default is 0.05.
#' @param union Logical indicator whether or not to use bootstrap union tests (\code{TRUE}) or not (\code{FALSE}), see Smeekes and Taylor (2012). Default is \code{TRUE}.
#' @details The false discovery rate FDR is defined as the expected proportion of false rejections relative to the total number of rejections.
#'
#' See \code{\link{boot_ur}} for details on the bootstrap algorithm and lag selection.
#' @return A list with the following components
#' \item{\code{rej_H0}}{Logical indicator whether the null hypothesis of a unit root is rejected (\code{TRUE}) or not (\code{FALSE});}
#' \item{\code{FDR_sequence}}{Details on the unit root tests: value of the test statistics and critical values.}
#' For the union test (\code{union = TRUE}), the output is arranged per time series. If \code{union = FALSE}, the output is arranged per time series, type of deterministic component (\code{deterministics}) and detrending method (\code{detrend}).
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Resampling-based bootstraps MBB and SB cannot handle missing values.}}{If the time series in \code{data} have different starting and end points (and thus some series contain \code{NA} values at the beginning and/or end of the sample, the resampling-based moving block bootstrap (MBB) and sieve bootstrap (SB) cannot be used, as they create holes (internal missings) in the bootstrap samples. Switch to another bootstrap method or truncate your sample to eliminate \code{NA} values.}
#' \item{\code{Warning: SB and SWB bootstrap only recommended for boot_ur; see help for details.}}{Although the sieve bootstrap methods \code{"SB"} and \code{"SWB"} can be used, Smeekes and Urbain (2014b) show that these are not suited to capture general forms of dependence across units, and using them for joint or multiple testing is not valid. This warning thereofre serves to recommend the user to consider a different bootstrap method.}
#' \item{\code{Warning: Deterministic specification in argument deterministics is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting deterministic components manually therefore has no effect.}
#' \item{\code{Warning: Detrending method in argument detrend is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting detrending methods manually therefore has no effect.}
#' }
#' @references Chang, Y. and Park, J. (2003). A sieve bootstrap for the test of a unit root. \emph{Journal of Time Series Analysis}, 24(4), 379-400.
#' @references Cavaliere, G. and Taylor, A.M.R (2009). Heteroskedastic time series with a unit root. \emph{Econometric Theory}, 25, 1228–1276.
#' @references Cavaliere, G., Phillips, P.C.B., Smeekes, S., and Taylor, A.M.R.
#' (2015). Lag length selection for unit root tests in the presence of nonstationary volatility. \emph{Econometric Reviews}, 34(4), 512-536.
#' @references Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996). Efficient tests for an autoregressive unit root. \emph{Econometrica}, 64(4), 813-836.
#' @references Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild bootstrap inference for nonparametric trends. \emph{Journal of Econometrics}, 214(1), 81-109.
#' @references Moon, H.R. and Perron, B. (2012). Beyond panel unit root tests: Using multiple testing to determine the non stationarity properties of individual series in a panel. Journal of Econometrics, 169(1), 29-33.
#' @references Ng, S. and Perron, P. (2001). Lag Length Selection and the Construction of Unit Root Tests with Good Size and Power. \emph{Econometrica}, 69(6), 1519-1554,
#' @references Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008). Bootstrap unit root tests: Comparison and extensions. \emph{Journal of Time Series Analysis}, 29(1), 371-401.
#' @references Palm, F. C., Smeekes, S., and Urbain, J.-.P. (2011). Cross-sectional dependence robust block bootstrap panel unit root tests. \emph{Journal of Econometrics}, 163(1), 85-104.
#' @references Paparoditis, E. and Politis, D.N. (2003). Residual-based block bootstrap for unit root testing. \emph{Econometrica}, 71(3), 813-855.
#' @references Perron, P. and Qu, Z. (2008). A simple modification to improve the finite sample properties of Ng and Perron's unit root tests. \emph{Economic Letters}, 94(1), 12-19.
#' @references Rho, Y. and Shao, X. (2019). Bootstrap-assisted unit root testing with piecewise locally stationary errors. \emph{Econometric Theory}, 35(1), 142-166.
#' @references Romano, J.P., Shaikh, A.M., and Wolf, M. (2008). Control of the false discovery rate under dependence using the bootstrap and subsampling. \emph{Test}, 17(3), 417.
#' @references Shao, X. (2010). The dependent wild bootstrap. \emph{Journal of the American Statistical Association}, 105(489), 218-235.
#' @references Shao, X. (2011). A bootstrap-assisted spectral test of white noise under unknown dependence. \emph{Journal of Econometrics}, 162, 213-224.
#' @references Smeekes, S. (2013). Detrending bootstrap unit root tests. \emph{Econometric Reviews}, 32(8), 869-891.
#' @references Smeekes, S. and Taylor, A.M.R. (2012). Bootstrap union tests for unit roots in the presence of nonstationary volatility. \emph{Econometric Theory}, 28(2), 422-456.
#' @references Smeekes, S. and Urbain, J.-P. (2014a). A multivariate invariance principle for modified wild bootstrap methods with an application to unit root testing. GSBE Research Memorandum No. RM/14/008, Maastricht University
#' @references Smeekes, S. and Urbain, J.-P. (2014b). On the applicability of the sieve bootstrap in time series panels. \emph{Oxford Bulletin of Economics and Statistics}, 76(1), 139-151.
#' @seealso \code{\link{boot_ur}}
#' @examples
#' # boot_fdr on GDP_BE and GDP_DE
#' two_series_boot_fdr <- boot_fdr(MacroTS[, 1:2], bootstrap = "MBB", B = 399)
#' print(two_series_boot_fdr)
#' @export
boot_fdr <- function(data, level = 0.05,  bootstrap = "AWB", B = 1999, block_length = NULL,
                     ar_AWB = NULL, union = TRUE, min_lag = 0, max_lag = NULL, criterion = "MAIC",
                     deterministics = NULL, detrend = NULL, criterion_scale = TRUE,
                     show_progress = FALSE, do_parallel = FALSE, cores = NULL){

  inputs <- do_tests_and_bootstrap(data = data, boot_sqt_test = FALSE, boot_ur_test = FALSE,
                                   level = level, bootstrap = bootstrap, B = B,
                                   block_length = block_length, ar_AWB = ar_AWB, union = union,
                                   min_lag = min_lag, max_lag = max_lag, criterion = criterion,
                                   deterministics = deterministics, detrend = detrend,
                                   criterion_scale = criterion_scale, steps = NULL, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel,
                                   cores = cores)

  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(data))
  }

  if (union) { # Union Tests
    bFDRout <- FDR_cpp(test_i = inputs$test_stats, t_star = inputs$test_stats_star, level = inputs$level)
    rej_H0 <- matrix(bFDRout$rej_H0 == 1, nrow = NCOL(data))
    rownames(rej_H0) <- var_names
    colnames(rej_H0) <- "Reject null"
    FDR_seq <- bFDRout$FDR_Tests[, -1, drop = FALSE]
    rownames(FDR_seq) <- var_names[bFDRout$FDR_Tests[, 1, drop = FALSE]]
    colnames(FDR_seq) <- c("tstat", "critical value")

    method_name <- "Bootstrap Union Tests with False Discovery Rate control"
  } else { # No Union Tests
      bFDRout <- FDR_cpp(test_i = matrix(inputs$tests_i[1, ], nrow = 1), t_star = inputs$t_star[ , 1,],
                         level = inputs$level)
      rej_H0 <- matrix(bFDRout$rej_H0 == 1, nrow = NCOL(data))
      rownames(rej_H0) <- var_names
      colnames(rej_H0) <- "Reject null"
      FDR_seq <- bFDRout$FDR_Tests[, -1, drop = FALSE]
      rownames(FDR_seq) <- var_names[bFDRout$FDR_Tests[, 1, drop = FALSE]]
      colnames(FDR_seq) <- c("tstat", "critical value")
      method_name <- "Bootstrap ADF Tests with False Discovery Rate control"
  }

  fdr_output <- list(method = method_name, data.name = "data", details = FDR_seq,
                     alternative = "less", null.value =  c("gamma" = 0), rejections = rej_H0,
                     estimate = NULL, statistic = NULL, p.value = NULL)
  class(fdr_output) <- "mult_htest"
  
  return(fdr_output)
}

#' Bootstrap Sequential Quantile Test
#' @description Performs the Bootstrap Sequential Quantile Test (BSQT) proposed by Smeekes (2015).
#' @inheritParams boot_ur
#' @param steps Numeric vector of quantiles or units to be tested. Default is to test each unit sequentially.
#' @param union Logical indicator whether or not to use bootstrap union tests (\code{TRUE}) or not (\code{FALSE}), see Smeekes and Taylor (2012). Default is \code{TRUE}.
#' @details The parameter \code{steps} can either be set as an increasing sequence of integers smaller or equal to the number of series \code{N}, or fractions of the total number of series (quantiles). For \code{N} time series, setting \code{steps = 0:N} means each unit should be tested sequentially. In this case the method is equivalent to the StepM method of Romano and Wolf (2005), and therefore controls the familywise error rate. To split the series in \code{K} equally sized groups, use \code{steps = 0:K / K}.
#'
#' By convention and in accordance with notation in Smeekes (2015), the first entry of the vector should be equal to zero, while the second entry indicates the end of the first group, and so on. If the initial \code{0} or final value (\code{1} or \code{N}) are omitted, they are automatically added by the function.
#'
#' See \code{\link{boot_ur}} for details on the bootstrap algorithm and lag selection.
#' @return A list with the following components
#' \item{\code{rej_H0}}{Logical indicator whether the null hypothesis of a unit root is rejected (\code{TRUE}) or not (\code{FALSE});}
#' \item{\code{BSQT_sequence}}{Details on the unit root tests: outcome of the sequential steps, value of the test statistics and p-values.}
#' For the union test (\code{union = TRUE}), the output is arranged per time series. If \code{union = FALSE}, the output is arranged per time series, type of deterministic component (\code{deterministics}) and detrending method (\code{detrend}).
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Resampling-based bootstraps MBB and SB cannot handle missing values.}}{If the time series in \code{data} have different starting and end points (and thus some series contain \code{NA} values at the beginning and/or end of the sample, the resampling-based moving block bootstrap (MBB) and sieve bootstrap (SB) cannot be used, as they create holes (internal missings) in the bootstrap samples. Switch to another bootstrap method or truncate your sample to eliminate \code{NA} values.}
#' \item{\code{Error: Invalid input values for steps: must be quantiles or positive integers.}}{Construction of \code{steps} does not satisfy the criteria listed under 'Details'.}
#' \item{\code{Warning: SB and SWB bootstrap only recommended for boot_ur; see help for details.}}{Although the sieve bootstrap methods \code{"SB"} and \code{"SWB"} can be used, Smeekes and Urbain (2014b) show that these are not suited to capture general forms of dependence across units, and using them for joint or multiple testing is not valid. This warning thereofre serves to recommend the user to consider a different bootstrap method.}
#' \item{\code{Warning: Deterministic specification in argument deterministics is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting deterministic components manually therefore has no effect.}
#' \item{\code{Warning: Detrending method in argument detrend is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting detrending methods manually therefore has no effect.}
#' }
#' @references Chang, Y. and Park, J. (2003). A sieve bootstrap for the test of a unit root. \emph{Journal of Time Series Analysis}, 24(4), 379-400.
#' @references Cavaliere, G. and Taylor, A.M.R (2009). Heteroskedastic time series with a unit root. \emph{Econometric Theory}, 25, 1228–1276.
#' @references Cavaliere, G., Phillips, P.C.B., Smeekes, S., and Taylor, A.M.R.
#' (2015). Lag length selection for unit root tests in the presence of nonstationary volatility. \emph{Econometric Reviews}, 34(4), 512-536.
#' @references Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996). Efficient tests for an autoregressive unit root. \emph{Econometrica}, 64(4), 813-836.
#' @references Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild bootstrap inference for nonparametric trends. \emph{Journal of Econometrics}, 214(1), 81-109.
#' @references Ng, S. and Perron, P. (2001). Lag Length Selection and the Construction of Unit Root Tests with Good Size and Power. \emph{Econometrica}, 69(6), 1519-1554,
#' @references Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008). Bootstrap unit root tests: Comparison and extensions. \emph{Journal of Time Series Analysis}, 29(1), 371-401.
#' @references Palm, F. C., Smeekes, S., and Urbain, J.-.P. (2011). Cross-sectional dependence robust block bootstrap panel unit root tests. \emph{Journal of Econometrics}, 163(1), 85-104.
#' @references Paparoditis, E. and Politis, D.N. (2003). Residual-based block bootstrap for unit root testing. \emph{Econometrica}, 71(3), 813-855.
#' @references Perron, P. and Qu, Z. (2008). A simple modification to improve the finite sample properties of Ng and Perron's unit root tests. \emph{Economic Letters}, 94(1), 12-19.
#' @references Rho, Y. and Shao, X. (2019). Bootstrap-assisted unit root testing with piecewise locally stationary errors. \emph{Econometric Theory}, 35(1), 142-166.
#' @references Romano, J. P. and Wolf, M. (2005). Stepwise multiple testing as formalized data snooping. \emph{Econometrica}, 73(4), 1237-1282.
#' #' @references Shao, X. (2010). The dependent wild bootstrap. \emph{Journal of the American Statistical Association}, 105(489), 218-235.
#' @references Shao, X. (2011). A bootstrap-assisted spectral test of white noise under unknown dependence. \emph{Journal of Econometrics}, 162, 213-224.
#' @references Smeekes, S. (2013). Detrending bootstrap unit root tests. \emph{Econometric Reviews}, 32(8), 869-891.
#' @references Smeekes, S. (2015). Bootstrap sequential tests to determine the order of integration of individual units in a time series panel. \emph{Journal of Time Series Analysis}, 36(3), 398-415.
#' @references Smeekes, S. and Taylor, A.M.R. (2012). Bootstrap union tests for unit roots in the presence of nonstationary volatility. \emph{Econometric Theory}, 28(2), 422-456.
#' @references Smeekes, S. and Urbain, J.-P. (2014a). A multivariate invariance principle for modified wild bootstrap methods with an application to unit root testing. GSBE Research Memorandum No. RM/14/008, Maastricht University
#' @references Smeekes, S. and Urbain, J.-P. (2014b). On the applicability of the sieve bootstrap in time series panels. \emph{Oxford Bulletin of Economics and Statistics}, 76(1), 139-151.
#' @seealso \code{\link{boot_ur}}
#' @examples
#' # boot_sqt on GDP_BE and GDP_DE
#' two_series_boot_sqt <- boot_sqt(MacroTS[, 1:2], bootstrap = "AWB", B = 399)
#' print(two_series_boot_fdr)
#' @export
boot_sqt <- function(data, steps = 0:NCOL(data), level = 0.05,  bootstrap = "AWB",
                     B = 1999, block_length = NULL, ar_AWB = NULL, union = TRUE,
                     min_lag = 0, max_lag = NULL, criterion = "MAIC", deterministics = NULL,
                     detrend = NULL, criterion_scale = TRUE, show_progress = FALSE,
                     do_parallel = FALSE, cores = NULL){

  inputs <- do_tests_and_bootstrap(data = data, boot_sqt_test = TRUE, boot_ur_test = FALSE,
                                   level = level, bootstrap = bootstrap, B = B,
                                   block_length = block_length, ar_AWB = ar_AWB, union = union,
                                   min_lag = min_lag, max_lag = max_lag, criterion = criterion,
                                   deterministics = deterministics, detrend = detrend,
                                   criterion_scale = criterion_scale, steps = steps, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel,
                                   cores = cores)

  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(data))
  }

  if (union) { # Union Tests
    BSQTout <- BSQT_cpp(pvec = inputs$p_vec, test_i = inputs$test_stats,
                        t_star = inputs$test_stats_star, level = inputs$level)
    rej_H0 <- matrix(BSQTout$rej_H0 == 1, nrow = NCOL(data))
    rownames(rej_H0) <- var_names
    colnames(rej_H0) <- "Reject H0"
    BSQT_seq <- BSQTout$BSQT_steps[, -3, drop = FALSE]
    rownames(BSQT_seq) <- paste("Step", 1:nrow(BSQT_seq))
    colnames(BSQT_seq) <- c("Unit H0", "Unit H1", "tstat", "p-value")

    method_name <- "Bootstrap Sequential Quantile Union Test"
  } else { # No Union Tests
    BSQTout <- BSQT_cpp(pvec = inputs$p_vec, test_i = matrix(inputs$tests_i[1, ], nrow = 1),
                        t_star = inputs$t_star[ , 1,], level = inputs$level)
    rej_H0 <- matrix(BSQTout$rej_H0 == 1, nrow = NCOL(data))
    rownames(rej_H0) <- var_names
    colnames(rej_H0) <- "Reject null"
    BSQT_seq <- BSQTout$BSQT_steps[, -3, drop = FALSE]
    rownames(BSQT_seq) <- paste("Step", 1:nrow(BSQT_seq))
    colnames(BSQT_seq) <- c("Unit H0", "Unit H1", "tstat", "p-value")
  }

  sqt_output <- list(method = method_name, data.name = "data", details = BSQT_seq,
                     alternative = "less", null.value =  c("gamma" = 0), rejections = rej_H0,
                     estimate = NULL, statistic = NULL, p.value = NULL)
  class(sqt_output) <- "mult_htest"
  return(sqt_output)
}

#' Panel Unit Root Test
#' @description Performs a test on a multivariate (panel) time series by testing the null hypothesis that all series have a unit root. The test is based on averaging the individual test statistics, also called the Group-Mean (GM) test in Palm, Smeekes and Urbain (2011).
#' @inheritParams boot_ur
#' @param union Logical indicator whether or not to use bootstrap union tests (\code{TRUE}) or not (\code{FALSE}), see Smeekes and Taylor (2012). Default is \code{TRUE}.
#' @details See \code{\link{boot_ur}} for details on the bootstrap algorithm and lag selection.
#'
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Resampling-based bootstraps MBB and SB cannot handle missing values.}}{If the time series in \code{data} have different starting and end points (and thus some series contain \code{NA} values at the beginning and/or end of the sample, the resampling-based moving block bootstrap (MBB) and sieve bootstrap (SB) cannot be used, as they create holes (internal missings) in the bootstrap samples. Switch to another bootstrap method or truncate your sample to eliminate \code{NA} values.}
#' \item{\code{Warning: SB and SWB bootstrap only recommended for boot_ur; see help for details.}}{Although the sieve bootstrap methods \code{"SB"} and \code{"SWB"} can be used, Smeekes and Urbain (2014b) show that these are not suited to capture general forms of dependence across units, and using them for joint or multiple testing is not valid. This warning thereofre serves to recommend the user to consider a different bootstrap method.}
#' \item{\code{Warning: Deterministic specification in argument deterministics is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting deterministic components manually therefore has no effect.}
#' \item{\code{Warning: Detrending method in argument detrend is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting detrending methods manually therefore has no effect.}
#' }
#' @return An object of class \code{htest} containing
#' \item{\code{method}}{The name of the hypothesis test. For boot_panel this is the Panel Bootstrap Group-Mean Union Test (for \code{union = TRUE}) or the "Panel Bootstrap Group-Mean Test" (for \code{union = FALSE});}
#' \item{\code{data.name}}{The name of the variable on which the test is performed.;}
#' \item{\code{null.value}}{The value of the (gamma) parameter of the lagged dependent variable in the ADF regression under the null hypothesis. Under the null, the series has a unit root. Testing the null of a unit root then boils down to testing the significance of the gamma parameter;}
#' \item{\code{alternative}}{A character string specifying the direction of the alternative hypothesis relative to the null value. The alternative postulates that the series is stationary;}
#' \item{\code{estimate}}{For the panel test, the estimated value of the (gamma) parameter of the lagged dependent variable in the ADF regression is not defined.;}
#' \item{\code{statistic}}{The value of the test statistic of the panel unit root test.;}
#' \item{\code{p.value}}{P-value of the panel unit root test.}
#' @references Chang, Y. and Park, J. (2003). A sieve bootstrap for the test of a unit root. \emph{Journal of Time Series Analysis}, 24(4), 379-400.
#' @references Cavaliere, G. and Taylor, A.M.R (2009). Heteroskedastic time series with a unit root. \emph{Econometric Theory}, 25, 1228–1276.
#' @references Cavaliere, G., Phillips, P.C.B., Smeekes, S., and Taylor, A.M.R.
#' (2015). Lag length selection for unit root tests in the presence of nonstationary volatility. \emph{Econometric Reviews}, 34(4), 512-536.
#' @references Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996). Efficient tests for an autoregressive unit root. \emph{Econometrica}, 64(4), 813-836.
#' @references Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild bootstrap inference for nonparametric trends. \emph{Journal of Econometrics}, 214(1), 81-109.
#' @references Ng, S. and Perron, P. (2001). Lag Length Selection and the Construction of Unit Root Tests with Good Size and Power. \emph{Econometrica}, 69(6), 1519-1554,
#' @references Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008). Bootstrap unit root tests: Comparison and extensions. \emph{Journal of Time Series Analysis}, 29(1), 371-401.
#' @references Palm, F. C., Smeekes, S., and Urbain, J.-.P. (2011). Cross-sectional dependence robust block bootstrap panel unit root tests. \emph{Journal of Econometrics}, 163(1), 85-104.
#' @references Paparoditis, E. and Politis, D.N. (2003). Residual-based block bootstrap for unit root testing. \emph{Econometrica}, 71(3), 813-855.
#' @references Perron, P. and Qu, Z. (2008). A simple modification to improve the finite sample properties of Ng and Perron's unit root tests. \emph{Economic Letters}, 94(1), 12-19.
#' @references Rho, Y. and Shao, X. (2019). Bootstrap-assisted unit root testing with piecewise locally stationary errors. \emph{Econometric Theory}, 35(1), 142-166.
#' @references Shao, X. (2010). The dependent wild bootstrap. \emph{Journal of the American Statistical Association}, 105(489), 218-235.
#' @references Shao, X. (2011). A bootstrap-assisted spectral test of white noise under unknown dependence. \emph{Journal of Econometrics}, 162, 213-224.
#' @references Smeekes, S. (2013). Detrending bootstrap unit root tests. \emph{Econometric Reviews}, 32(8), 869-891.
#' @references Smeekes, S. and Taylor, A.M.R. (2012). Bootstrap union tests for unit roots in the presence of nonstationary volatility. \emph{Econometric Theory}, 28(2), 422-456.
#' @references Smeekes, S. and Urbain, J.-P. (2014a). A multivariate invariance principle for modified wild bootstrap methods with an application to unit root testing. GSBE Research Memorandum No. RM/14/008, Maastricht University
#' @references Smeekes, S. and Urbain, J.-P. (2014b). On the applicability of the sieve bootstrap in time series panels. \emph{Oxford Bulletin of Economics and Statistics}, 76(1), 139-151.
#' @seealso \code{\link{boot_ur}}
#' @examples
#' # boot_panel on GDP_BE and GDP_DE
#' two_series_boot_panel <- boot_panel(MacroTS[, 1:2], bootstrap = "AWB", B = 399)
#' print(two_series_boot_panel)
#' @export
boot_panel <- function(data, level = 0.05,  bootstrap = "AWB", B = 1999, block_length = NULL,
                       ar_AWB = NULL, union = TRUE, min_lag = 0, max_lag = NULL, criterion = "MAIC",
                       deterministics = NULL, detrend = NULL, criterion_scale = TRUE,
                       show_progress = FALSE, do_parallel = FALSE, cores = NULL){

  inputs <- do_tests_and_bootstrap(data = data, boot_sqt_test = FALSE, boot_ur_test = FALSE,
                                   level = level, bootstrap = bootstrap, B = B,
                                   block_length = block_length, ar_AWB = ar_AWB, union = union,
                                   min_lag = min_lag, max_lag = max_lag, criterion = criterion,
                                   deterministics = deterministics, detrend = detrend,
                                   criterion_scale = criterion_scale, steps = NULL, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel,
                                   cores = cores)


  if (union) { # Union Test
    GM_test <- mean(inputs$test_stats)
    t_star <- rowMeans(inputs$test_stats_star)
    p_val <- mean(t_star < GM_test)
    method_name <- "Panel Bootstrap Group-Mean Union Test"
  } else { # No Union Test
    GM_test <- rowMeans(inputs$tests_i)
    t_star <- apply(inputs$t_star, 1:2, mean)
    p_val <- sapply(1, function(i){mean(t_star[, i] < GM_test[i])})
    method_name <- "Panel Bootstrap Group-Mean Test"
  }

  attr(GM_test, "names") <- "tstat"
  gamma_hat <- NA
  attr(gamma_hat, "names") <- "gamma"
  panel_output <- list(method = method_name, data.name = "panel", null.value = c("gamma" = 0),
                       alternative = "less", estimate = gamma_hat, statistic = GM_test, p.value = p_val)
  class(panel_output) <- "htest"
  return(panel_output)
}