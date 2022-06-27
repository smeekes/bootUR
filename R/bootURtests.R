#' Individual Unit Root Tests without multiple testing control
#' @description This function performs bootstrap unit root tests on each time series individually.
#' @param data A \eqn{T}-dimensional vector or a (\eqn{T} x \eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations to be tested for unit roots. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}), or a data frame, as long as each column represents a single time series.
#' @param data_name Optional name for the data, to be used in the output. The default uses the name of the 'data' argument.
#' @param bootstrap String for bootstrap method to be used. Options are
#' \describe{
#' \item{\code{"MBB"}}{Moving block bootstrap (Paparoditis and Politis, 2003; Palm, Smeekes and Urbain, 2011);}
#' \item{\code{"BWB"}}{Block wild bootstrap (Shao, 2011; Smeekes and Urbain, 2014a);}
#' \item{\code{"DWB"}}{Dependent wild bootstrap (Shao, 2010; Smeekes and Urbain, 2014a; Rho and Shao, 2019);}
#' \item{\code{"AWB"}}{Autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2020), this is the default;}
#' \item{\code{"SB"}}{Sieve bootstrap (Chang and Park, 2003; Palm, Smeekes and Urbain, 2008; Smeekes, 2013);}
#' \item{\code{"SWB"}}{Sieve wild bootstrap (Cavaliere and Taylor, 2009; Smeekes and Taylor, 2012).}
#' }
#' @param B Number of bootstrap replications. Default is 1999.
#' @param block_length Desired 'block length' in the bootstrap. For the MBB, BWB and DWB bootstrap, this is a genuine block length. For the AWB bootstrap, the block length is transformed into an autoregressive parameter via the formula \eqn{0.01^(1/block_length)} as in Smeekes and Urbain (2014a); this can be overwritten by setting \code{ar_AWB} directly. Default sets the block length as a function of the time series length T, via the rule \eqn{block_length = 1.75 T^(1/3)} of Palm, Smeekes and Urbain (2011).
#' @param ar_AWB Autoregressive parameter used in the AWB bootstrap method (\code{bootstrap = "AWB"}). Can be used to set the parameter directly rather than via the default link to the block length.
#' @param union_quantile The quantile of the bootstrap distribution used for scaling the individual statistics in the union. Ideally this should equal the desired significance level of the test. Default is 0.05. This parameter is overwritten when a significance level is provided in the argument \code{level}.
#' @param level The desired significance level of the test (optional). This is only used for multivariate series to be able to provide a boolean vector with rejections of the null hypothesis or not for easy post-processing. Default is \code{NULL}, in which case no such vector is given.
#' @param union Logical indicator whether or not to use bootstrap union tests (\code{TRUE}) or not (\code{FALSE}), see Smeekes and Taylor (2012). Default is \code{TRUE}.
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
#' @param min_lag Minimum lag length in the augmented Dickey-Fuller regression. Default is 0.
#' @param max_lag Maximum lag length in the augmented Dickey-Fuller regression. Default uses the sample size-based rule \eqn{12(T/100)^{1/4}}.
#' @param criterion String for information criterion used to select the lag length in the augmented Dickey-Fuller regression. Options are: \code{"AIC"}, \code{"BIC"}, \code{"MAIC"}, \code{"MBIC"}. Default is \code{"MAIC"} (Ng and Perron, 2001).
#' @param criterion_scale Logical indicator whether or not to use the rescaled information criteria of Cavaliere et al. (2015) (\code{TRUE}) or not (\code{FALSE}). Default is \code{TRUE}.
#' @param show_progress Logical indicator whether a bootstrap progress update should be printed to the console. Default is FALSE.
#' @param do_parallel Logical indicator whether bootstrap loop should be executed in parallel. Default is TRUE.
#' @param cores The number of cores to be used in the parallel loops. Default is to use all but one.
#' @details The options encompass many test proposed in the literature. \code{detrend = "OLS"} gives the standard augmented Dickey-Fuller test, while \code{detrend = "QD"} provides the DF-GLS test of Elliott, Rothenberg and Stock (1996). The bootstrap algorithm is always based on a residual bootstrap (under the alternative) to obtain residuals rather than a difference-based bootstrap (under the null), see e.g. Palm, Smeekes and Urbain (2008).
#'
#' Lag length selection is done automatically in the ADF regression with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. For very short time series (fewer than 50 time points) the maximum lag length is adjusted downward to avoid potential multicollinearity issues in the bootstrap. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `min_lag` and maximum lag length `max_lag` for the selection algorithm equal to the desired lag length.
#'
#' @return An object of class \code{"bootUR"}, \code{"\*"}, where \code{"\*"} is \code{"mult_htest"} for multiple time series or \code{"htest"} for single time series, with the following components:
#' \item{\code{method}}{The name of the hypothesis test method;}
#' \item{\code{data.name}}{The name of the data on which the method is performed;}
#' \item{\code{null.value}}{The value of the (gamma) parameter of the lagged dependent variable in the ADF regression under the null hypothesis. Under the null, the series has a unit root. Testing the null of a unit root then boils down to testing the significance of the gamma parameter;}
#' \item{\code{alternative}}{A character string specifying the direction of the alternative hypothesis relative to the null value. The alternative postulates that the series is stationary;}
#' \item{\code{estimate}}{The estimated value(s) of the (gamma) parameter of the lagged dependent variable in the ADF regressions. Note that for the union test (\code{union = TRUE}), this estimate is not defined, hence \code{NA} is returned;}
#' \item{\code{statistic}}{The value(s) of the test statistic of the unit root test(s);}
#' \item{\code{p.value}}{The p-value(s) of the unit root test(s);}
#' \item{\code{rejections}}{For \code{"mult_htest"} only. A vector with logical indicators for each time series whether the null hypothesis of a unit root is rejected (\code{TRUE}) or not (\code{FALSE}). This is only supplied when an optional significance level is given, otherwise \code{NULL} is returned;}
#' \item{\code{details}}{A list containing the detailed outcomes of the performed tests, such as selected lags, individual estimates and p-values.}
#' \item{\code{series.names}}{For \code{"mult_htest"} only. The names of the series that the tests are performed on;}
#' \item{\code{specifications}}{The specifications used in the test(s).}
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
#' two_series_boot_ur <- boot_ur(MacroTS[, 1:2], bootstrap = "MBB", B = 199,
#'                               do_parallel = FALSE, show_progress = FALSE)
#' print(two_series_boot_ur)
#' @export
boot_ur <- function(data, data_name = NULL, bootstrap = "AWB", B = 1999, block_length = NULL,
                    ar_AWB = NULL, level = NULL, union = TRUE, union_quantile = 0.05,
                    deterministics = NULL, detrend = NULL, min_lag = 0, max_lag = NULL,
                    criterion = "MAIC", criterion_scale = TRUE, show_progress = TRUE,
                    do_parallel = TRUE, cores = NULL){

  if (!is.null(level)) {
    union_quantile <- level
  }
  inputs <- do_tests_and_bootstrap(data = data, boot_sqt_test = FALSE, boot_ur_test = TRUE,
                                   level = union_quantile, bootstrap = bootstrap, B = B,
                                   block_length = block_length, ar_AWB = ar_AWB,
                                   union = union, min_lag = min_lag, max_lag = max_lag,
                                   criterion = criterion, deterministics = deterministics,
                                   detrend = detrend, criterion_scale = criterion_scale,
                                   steps = NULL, h_rs = 0.1, show_progress = show_progress,
                                   do_parallel = do_parallel, cores = cores,
                                   data_name = data_name)

  if (is.null(data_name)) {
    data_name <- deparse(substitute(data))
  }
  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    if (NCOL(data) > 1) {
      var_names <- paste0(data_name, " var", 1:NCOL(data))
    } else {
      var_names <- data_name
    }
  }

  spec <- list("bootstrap" = bootstrap, "B" = B, "block_length" = inputs$inputs$l,
               "ar_AWB" = inputs$inputs$ar_AWB, "level" = level,
               "union" = union, "union_quantile" = inputs$inputs$union_quantile,
               "deterministics" = inputs$inputs$deterministics,
               "detrend" = inputs$inputs$detrend, "min_lag" = min_lag,
               "max_lag" = inputs$inputs$p_max, "criterion" = inputs$inputs$criterion,
               "criterion_scale" = inputs$inputs$criterion_scale, "mult_test_ctrl" = "none")

  # Results
  if (union) { # Union test
    iADFout <- iADF_cpp(test_i = inputs$test_stats, t_star = inputs$test_stats_star)
    iADFout <- cbind(rep(NA, nrow(iADFout)), t(inputs$test_stats), iADFout)
    if (NCOL(data) > 1) {
      rownames(iADFout) <- var_names
      colnames(iADFout) <- c("gamma", "tstat", "p-value")
    }
    # Parameter estimates, tstats and p-values. Note: Parameter Estimates not defined for union test

    details <- list("individual estimates" = t(inputs$indiv_par_est),
                    "individual statistics" = t(inputs$indiv_test_stats),
                    "individual p-values" = inputs$indiv_pval,
                    "selected lags" = t(inputs$indiv_lags),
                    "txt_null" = "Series has a unit root",
                    "txt_alternative" = "Series is stationary")
    rownames(details$"individual estimates") <- var_names
    colnames(details$"individual estimates") <- c(t(outer(c("OLS", "QD"),
                                                  c("intercept", "intercept and trend"),
                                                  function(x,y){paste0(x, "/", y)})))
    rownames(details$"individual statistics") <- var_names
    colnames(details$"individual statistics") <- c(t(outer(c("OLS", "QD"),
                                                   c("intercept", "intercept and trend"),
                                                   function(x,y){paste0(x, "/", y)})))
    rownames(details$"individual p-values") <- var_names
    colnames(details$"individual p-values") <- c(t(outer(c("OLS", "QD"),
                                                 c("intercept", "intercept and trend"),
                                                 function(x,y){paste0(x, "/", y)})))
    rownames(details$"selected lags") <- var_names
    colnames(details$"selected lags") <- c(t(outer(c("OLS", "QD"),
                                           c("intercept", "intercept and trend"),
                                           function(x,y){paste0(x, "/", y)})))
  } else { # No union test
    iADFout <- iADF_cpp(test_i = matrix(inputs$indiv_test_stats[1, ], nrow = 1),
                        t_star = matrix(inputs$t_star[ , 1, ], nrow = B))
    iADFout <- cbind(t(inputs$indiv_par_est), t(inputs$indiv_test_stats), iADFout)
    # Parameter estimates, tstats and p-values

    details <- list("individual estimates" = t(inputs$indiv_par_est),
                    "individual statistics" = t(inputs$indiv_test_stats),
                    "individual p-values" = inputs$indiv_pval,
                    "selected lags" = t(inputs$indiv_lags),
                    "txt_null" = "Series has a unit root",
                    "txt_alternative" = "Series is stationary")
    rownames(details$"individual estimates") <- var_names
    colnames(details$"individual estimates") <- paste0(inputs$inputs$detrend, "/",
                                                       inputs$inputs$deterministics)
    rownames(details$"individual statistics") <- var_names
    colnames(details$"individual statistics") <- paste0(inputs$inputs$detrend, "/",
                                                        inputs$inputs$deterministics)
    rownames(details$"individual p-values") <- var_names
    colnames(details$"individual p-values") <- paste0(inputs$inputs$detrend, "/",
                                                      inputs$inputs$deterministics)
    rownames(details$"selected lags") <- var_names
    colnames(details$"selected lags") <- paste0(inputs$inputs$detrend, "/",
                                                inputs$inputs$deterministics)

    if (NCOL(data) > 1) {
      rownames(iADFout) <- var_names
      colnames(iADFout) <- c("gamma", "statistic", "p.value")
      colnames(iADFout) <- c("gamma", "statistic", "p.value")
    }
  }

  if (!is.null(level)) {
    rej_H0 <- (iADFout[, 3] < level)
  } else {
    rej_H0 <- NULL
  }

  if (NCOL(data) > 1) {
    if (union) {
      method_name <- paste0(bootstrap,
                           " bootstrap union test on each individual series (no multiple testing correction)")
    } else {
      method_name <- paste0(bootstrap, " bootstrap ", inputs$inputs$name,
                            " test (with " , inputs$inputs$deterministics,
                            ") on each individual series (no multiple testing correction)")
    }
    boot_ur_output <- list(method = method_name, data.name = data_name,
                           null.value =  c("gamma" = 0), alternative = "less",
                           estimate = iADFout[, 1], statistic = iADFout[, 2],
                           p.value = iADFout[, 3], rejections = rej_H0,
                           details = details, series.names = var_names, specifications = spec)
    class(boot_ur_output) <- c("bootUR", "mult_htest")


  } else {
    param <- drop(iADFout[1, 1])
    attr(param, "names") <- "gamma"
    iADFtstat <- iADFout[1, 2]
    attr(iADFtstat, "names") <- "tstat"
    p_val <- drop(iADFout[1, 3])
    attr(p_val, "names") <- "p-value"

    if (union) {
      method_name <- paste0(bootstrap, " bootstrap union test on a single time series")
    } else {
      method_name <- paste0(bootstrap, " bootstrap ", inputs$inputs$detrend,
                           " test (with " , inputs$inputs$deterministics,
                           ") on a single time series")
    }
    boot_ur_output <- list(method = method_name, data.name = var_names,
                           null.value = c("gamma" = 0), alternative = "less",
                           estimate = param, statistic = iADFtstat, p.value = p_val,
                           details = details, series.names = var_names, specifications = spec)
    class(boot_ur_output) <- c("bootUR", "htest")
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
#' \item{\code{"SWB"}}{Sieve wild bootstrap (Cavaliere and Taylor, 2009; Smeekes and Taylor, 2012).}
#' }
#' @details The options encompass many test proposed in the literature. \code{detrend = "OLS"} gives the standard augmented Dickey-Fuller test, while \code{detrend = "QD"} provides the DF-GLS test of Elliott, Rothenberg and Stock (1996). The bootstrap algorithm is always based on a residual bootstrap (under the alternative) to obtain residuals rather than a difference-based bootstrap (under the null), see e.g. Palm, Smeekes and Urbain (2008).
#'
#' Lag length selection is done automatically in the ADF regression with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. For very short time series (fewer than 50 time points) the maximum lag length is adjusted downward to avoid potential multicollinearity issues in the bootstrap. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `min_lag` and maximum lag length `max_lag` for the selection algorithm equal to the desired lag length.
#' @export
#' @return An object of class \code{"bootUR"}, \code{"htest"} with the following components:
#' \item{\code{method}}{The name of the hypothesis test method;}
#' \item{\code{data.name}}{The name of the data on which the method is performed;}
#' \item{\code{null.value}}{The value of the (gamma) parameter of the lagged dependent variable in the ADF regression under the null hypothesis. Under the null, the series has a unit root. Testing the null of a unit root then boils down to testing the significance of the gamma parameter;}
#' \item{\code{alternative}}{A character string specifying the direction of the alternative hypothesis relative to the null value. The alternative postulates that the series is stationary;}
#' \item{\code{estimate}}{The estimated value of the (gamma) parameter of the lagged dependent variable in the ADF regression.;}
#' \item{\code{statistic}}{The value of the test statistic of the unit root test;}
#' \item{\code{p.value}}{The p-value of the unit root test;}
#' \item{\code{details}}{A list containing the detailed outcomes of the performed test, such as selected lags, individual estimates and p-values.}
#' \item{\code{specifications}}{The specifications used in the test.}
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
#' GDP_BE_adf <- boot_adf(MacroTS[, 1], B = 199, deterministics = "trend", detrend = "OLS",
#'                        do_parallel = FALSE, show_progress = FALSE)
#' print(GDP_BE_adf)
boot_adf <- function(data, data_name = NULL, bootstrap = "AWB", B = 1999,
                     block_length = NULL, ar_AWB = NULL, deterministics = "intercept",
                     min_lag = 0, max_lag = NULL, criterion = "MAIC",
                     detrend = "OLS", criterion_scale = TRUE, show_progress = TRUE,
                     do_parallel = TRUE, cores = NULL){

  if (NCOL(data) > 1) {
    stop("Multiple time series not allowed. Switch to a multivariate method such as boot_ur,
         or change argument data to a univariate time series.")
  }
  if (is.null(data_name)) {
    if (is.null(colnames(data))) {
      data_name <- deparse(substitute(data))
    } else {
      data_name <- colnames(data)[1]
    }
  }
  out <- boot_ur(data = data, bootstrap = bootstrap, B = B,
                 block_length = block_length, ar_AWB = ar_AWB, union = FALSE,
                 min_lag = min_lag, max_lag = max_lag, criterion = criterion,
                 deterministics = deterministics, detrend = detrend,
                 criterion_scale = criterion_scale, show_progress = show_progress,
                 do_parallel = do_parallel, cores = cores, data_name = data_name)

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
#' \item{\code{"SWB"}}{Sieve wild bootstrap (Cavaliere and Taylor, 2009; Smeekes and Taylor, 2012).}
#' }
#' @details The union is taken over the combination of tests with intercept only and intercept plus trend, coupled with OLS detrending and QD detrending, as in Harvey, Leybourne and Taylor (2012) and Smeekes an Taylor (2012). The bootstrap algorithm is always based on a residual bootstrap (under the alternative) to obtain residuals rather than a difference-based bootstrap (under the null), see e.g. Palm, Smeekes and Urbain (2008).
#'
#' Lag length selection is done automatically in the ADF regressions with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `min_lag` and maximum lag length `max_lag` for the selection algorithm equal to the desired lag length.
#' @export
#' @return An object of class \code{"bootUR"}, \code{"htest"} with the following components:
#' \item{\code{method}}{The name of the hypothesis test method;}
#' \item{\code{data.name}}{The name of the variable on which the method is performed;}
#' \item{\code{null.value}}{The value of the (gamma) parameter of the lagged dependent variable in the ADF regression under the null hypothesis. Under the null, the series has a unit root. Testing the null of a unit root then boils down to testing the significance of the gamma parameter;}
#' \item{\code{alternative}}{A character string specifying the direction of the alternative hypothesis relative to the null value. The alternative postulates that the series is stationary;}
#' \item{\code{estimate}}{For the union test, the estimated value of the (gamma) parameter of the lagged dependent variable in the ADF regression is not defined, hence NA is given;}
#' \item{\code{statistic}}{The value of the test statistic of the unit root test;}
#' \item{\code{p.value}}{The p-value of the unit root test;}
#' \item{\code{details}}{A list containing the detailed outcomes of the performed tests, such as selected lags, individual estimates and p-values.}
#' \item{\code{specifications}}{The specifications used in the test.}
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
#' GDP_BE_df <- boot_union(MacroTS[, 1], B = 199, do_parallel = FALSE, show_progress = FALSE)
#' print(GDP_BE_df)
boot_union <- function(data, data_name = NULL, bootstrap = "AWB", B = 1999, block_length = NULL,
                       ar_AWB = NULL, min_lag = 0, max_lag = NULL, criterion = "MAIC",
                       criterion_scale = TRUE, union_quantile = 0.05, show_progress = TRUE,
                       do_parallel = TRUE, cores = NULL){

  if (NCOL(data) > 1) {
    stop("Multiple time series not allowed. Switch to a multivariate method such as boot_ur,
         or change argument data to a univariate time series.")
  }
  if (is.null(data_name)) {
    if (is.null(colnames(data))) {
      data_name <- deparse(substitute(data))
    } else {
      data_name <- colnames(data)[1]
    }
  }
  out <- boot_ur(data = data, union_quantile = union_quantile, bootstrap = bootstrap, B = B,
                 block_length = block_length, ar_AWB = ar_AWB, union = TRUE, min_lag = min_lag,
                 max_lag = max_lag, criterion = criterion, criterion_scale = criterion_scale,
                 show_progress = show_progress, do_parallel = do_parallel,
                 cores = cores, data_name = data_name)

  return(out)
}

#' Bootstrap Unit Root Tests with False Discovery Rate control
#' @description Controls for multiple testing by controlling the false discovery rate (FDR), see Moon and Perron (2012) and Romano, Shaikh and Wolf (2008).
#' @inheritParams boot_ur
#' @param FDR_level Desired False Discovery Rate level of the unit root tests. Default is 0.05.
#' @param union Logical indicator whether or not to use bootstrap union tests (\code{TRUE}) or not (\code{FALSE}), see Smeekes and Taylor (2012). Default is \code{TRUE}.
#' @details The false discovery rate FDR is defined as the expected proportion of false rejections relative to the total number of rejections.
#'
#' See \code{\link{boot_ur}} for details on the bootstrap algorithm and lag selection.
#' @return An object of class \code{"bootUR"}, \code{"mult_htest"} with the following components:
#' \item{\code{method}}{The name of the hypothesis test method;}
#' \item{\code{data.name}}{The name of the data on which the method is performed;}
#' \item{\code{null.value}}{The value of the (gamma) parameter of the lagged dependent variable in the ADF regression under the null hypothesis. Under the null, the series has a unit root. Testing the null of a unit root then boils down to testing the significance of the gamma parameter;}
#' \item{\code{alternative}}{A character string specifying the direction of the alternative hypothesis relative to the null value. The alternative postulates that the series is stationary;}
#' \item{\code{estimate}}{The estimated values of the (gamma) parameter of the lagged dependent variable in the ADF regressions. Note that for the union test (\code{union = TRUE}), this estimate is not defined, hence NA is returned;}
#' \item{\code{statistic}}{The value of the test statistic of the unit root tests;}
#' \item{\code{p.value}}{A vector with \code{NA} values, as p-values are not available for the FDR method;}
#' \item{\code{rejections}}{A vector with logical indicators for each time series whether the null hypothesis of a unit root is rejected (\code{TRUE}) or not (\code{FALSE});}
#' \item{\code{details}}{A list containing the detailed outcomes of the performed tests, such as selected lags, individual estimates and p-values. In addtion, the slot \code{FDR} contains a matrix with for each step the test statistics and critical value, up to non-rejection.}
#' \item{\code{series.names}}{The names of the series that the tests are performed on;}
#' \item{\code{specifications}}{The specifications used in the test(s).}
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
#' two_series_boot_fdr <- boot_fdr(MacroTS[, 1:2], bootstrap = "MBB", B = 199,
#'                                 do_parallel = FALSE, show_progress = FALSE)
#' print(two_series_boot_fdr)
#' @export
boot_fdr <- function(data, data_name = NULL, bootstrap = "AWB", B = 1999, block_length = NULL,
                     ar_AWB = NULL, FDR_level = 0.05, union = TRUE, deterministics = NULL,
                     detrend = NULL, min_lag = 0, max_lag = NULL, criterion = "MAIC",
                     criterion_scale = TRUE, show_progress = TRUE, do_parallel = TRUE,
                     cores = NULL){

  inputs <- do_tests_and_bootstrap(data = data, boot_sqt_test = FALSE, boot_ur_test = FALSE,
                                   level = FDR_level, bootstrap = bootstrap, B = B,
                                   block_length = block_length, ar_AWB = ar_AWB, union = union,
                                   min_lag = min_lag, max_lag = max_lag, criterion = criterion,
                                   deterministics = deterministics, detrend = detrend,
                                   criterion_scale = criterion_scale, steps = NULL, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel,
                                   cores = cores, data_name = data_name)

  if (is.null(data_name)) {
    data_name <- deparse(substitute(data))
  }
  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    if (NCOL(data) > 1) {
      var_names <- paste0(data_name, " var", 1:NCOL(data))
    } else {
      var_names <- data_name
    }
  }

  spec <- list("bootstrap" = bootstrap, "B" = B, "block_length" = inputs$inputs$l,
               "ar_AWB" = inputs$inputs$ar_AWB, "FDR_level" = FDR_level, "union" = union,
               "deterministics" = inputs$inputs$deterministics,
               "detrend" = inputs$inputs$detrend, "min_lag" = min_lag,
               "max_lag" = inputs$inputs$p_max, "criterion" = inputs$inputs$criterion,
               "criterion_scale" = inputs$inputs$criterion_scale, "mult_test_ctrl" = "FDR")

  if (union) { # Union Tests
    bFDRout <- FDR_cpp(test_i = inputs$test_stats, t_star = inputs$test_stats_star,
                       level = inputs$level)
    estimates <- rep(NA, NCOL(data))
    tstats <- drop(inputs$test_stats)
    method_name <- paste0(bootstrap, " bootstrap union test with false discovery rate control")
    details <- list("individual estimates" = t(inputs$indiv_par_est),
                    "individual statistics" = t(inputs$indiv_test_stats),
                    "individual p-values" = inputs$indiv_pval,
                    "selected lags" = t(inputs$indiv_lags),
                    "txt_null" = "Series has a unit root",
                    "txt_alternative" = "Series is stationary")
    rownames(details$"individual estimates") <- var_names
    colnames(details$"individual estimates") <- c(t(outer(c("OLS", "QD"),
                                                  c("intercept", "intercept and trend"),
                                                  function(x,y){paste0(x, "/", y)})))
    rownames(details$"individual statistics") <- var_names
    colnames(details$"individual statistics") <- c(t(outer(c("OLS", "QD"),
                                                   c("intercept", "intercept and trend"),
                                                   function(x,y){paste0(x, "/", y)})))
    rownames(details$"individual p-values") <- var_names
    colnames(details$"individual p-values") <- c(t(outer(c("OLS", "QD"),
                                                 c("intercept", "intercept and trend"),
                                                 function(x,y){paste0(x, "/", y)})))
    rownames(details$"selected lags") <- var_names
    colnames(details$"selected lags") <- c(t(outer(c("OLS", "QD"),
                                           c("intercept", "intercept and trend"),
                                           function(x,y){paste0(x, "/", y)})))
  } else { # No Union Tests
    bFDRout <- FDR_cpp(test_i = matrix(inputs$indiv_test_stats[1, ], nrow = 1),
                       t_star = inputs$t_star[ , 1,], level = inputs$level)
    estimates <- t(inputs$indiv_par_est)
    tstats <- drop(inputs$tests_i[1, ])
    method_name <- paste0(bootstrap, " bootstrap ", inputs$inputs$name, " tests (with " ,
                          inputs$inputs$deterministics, ") with false discovery rate control")
    details <- list("individual estimates" = t(inputs$indiv_par_est),
                    "individual statistics" = t(inputs$indiv_test_stats),
                    "individual p-values" = inputs$indiv_pval,
                    "selected lags" = t(inputs$indiv_lags),
                    "txt_null" = "Series has a unit root",
                    "txt_alternative" = "Series is stationary")
    rownames(details$"individual estimates") <- var_names
    colnames(details$"individual estimates") <- paste0(inputs$inputs$detrend, "/",
                                                       inputs$inputs$deterministics)
    rownames(details$"individual statistics") <- var_names
    colnames(details$"individual statistics") <- paste0(inputs$inputs$detrend, "/",
                                                        inputs$inputs$deterministics)
    rownames(details$"individual p-values") <- var_names
    colnames(details$"individual p-values") <- paste0(inputs$inputs$detrend, "/",
                                                      inputs$inputs$deterministics)
    rownames(details$"selected lags") <- var_names
    colnames(details$"selected lags") <- paste0(inputs$inputs$detrend, "/",
                                                inputs$inputs$deterministics)
  }
  rej_H0 <- matrix(bFDRout$rej_H0 == 1, nrow = NCOL(data))
  rownames(rej_H0) <- var_names
  colnames(rej_H0) <- "Reject null"
  FDR_seq <- bFDRout$FDR_Tests[, -1, drop = FALSE]
  rownames(FDR_seq) <- var_names[bFDRout$FDR_Tests[, 1, drop = FALSE]]
  colnames(FDR_seq) <- c("tstat", "critical value")
  p_vals <- rep(NA, NCOL(data))
  names(estimates) <- names(tstats) <- names(p_vals) <- var_names
  details$FDR <- FDR_seq

  fdr_output <- list(method = method_name, data.name = data_name,
                     null.value =  c("gamma" = 0), alternative = "less",
                     estimate = estimates, statistic = tstats, p.value = p_vals,
                     rejections = rej_H0, details = details, series.names = var_names,
                     specifications = spec)
  class(fdr_output) <- c("bootUR", "mult_htest")

  return(fdr_output)
}

#' Bootstrap Sequential Quantile Test
#' @description Performs the Bootstrap Sequential Quantile Test (BSQT) proposed by Smeekes (2015).
#' @inheritParams boot_ur
#' @param SQT_level Desired significance level of the sequential tests performed. Default is 0.05.
#' @param steps Numeric vector of quantiles or units to be tested. Default is to test each unit sequentially.
#' @param union Logical indicator whether or not to use bootstrap union tests (\code{TRUE}) or not (\code{FALSE}), see Smeekes and Taylor (2012). Default is \code{TRUE}.
#' @details The parameter \code{steps} can either be set as an increasing sequence of integers smaller or equal to the number of series \code{N}, or fractions of the total number of series (quantiles). For \code{N} time series, setting \code{steps = 0:N} means each unit should be tested sequentially. In this case the method is equivalent to the StepM method of Romano and Wolf (2005), and therefore controls the familywise error rate. To split the series in \code{K} equally sized groups, use \code{steps = 0:K / K}.
#'
#' By convention and in accordance with notation in Smeekes (2015), the first entry of the vector should be equal to zero, while the second entry indicates the end of the first group, and so on. If the initial \code{0} or final value (\code{1} or \code{N}) are omitted, they are automatically added by the function.
#'
#' See \code{\link{boot_ur}} for details on the bootstrap algorithm and lag selection.
#' @return An object of class \code{"bootUR"}, \code{"mult_htest"} with the following components:
#' \item{\code{method}}{The name of the hypothesis test method;}
#' \item{\code{data.name}}{The name of the data on which the method is performed;}
#' \item{\code{null.value}}{The value of the (gamma) parameter of the lagged dependent variable in the ADF regression under the null hypothesis. Under the null, the series has a unit root. Testing the null of a unit root then boils down to testing the significance of the gamma parameter;}
#' \item{\code{alternative}}{A character string specifying the direction of the alternative hypothesis relative to the null value. The alternative postulates that the series is stationary;}
#' \item{\code{estimate}}{The estimated values of the (gamma) parameter of the lagged dependent variable in the ADF regressions. Note that for the union test (\code{union = TRUE}), this estimate is not defined, hence NA is returned;}
#' \item{\code{statistic}}{The value of the test statistic of the unit root tests;}
#' \item{\code{p.value}}{A vector with \code{NA} values, as p-values per inidividual series are not available.The p-value for each test in the sequence can be found in \code{details};}
#' \item{\code{rejections}}{A vector with logical indicators for each time series whether the null hypothesis of a unit root is rejected (\code{TRUE}) or not (\code{FALSE});}
#' \item{\code{details}}{A list containing the detailed outcomes of the performed tests, such as selected lags, individual estimates and p-values. In addtion, the slot \code{FDR} contains a matrix with for each step the stationary units under the null and alternative hypothesis, the test statistic and the p-value;}
#' \item{\code{series.names}}{The names of the series that the tests are performed on;}
#' \item{\code{specifications}}{The specifications used in the tests.}
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
#' two_series_boot_sqt <- boot_sqt(MacroTS[, 1:2], bootstrap = "AWB", B = 199,
#'                                 do_parallel = FALSE, show_progress = FALSE)
#' print(two_series_boot_sqt)
#' @export
boot_sqt <- function(data, data_name = NULL, steps = 0:NCOL(data), bootstrap = "AWB",
                     B = 1999, block_length = NULL, ar_AWB = NULL, SQT_level = 0.05, union = TRUE,
                     deterministics = NULL, detrend = NULL, min_lag = 0, max_lag = NULL,
                     criterion = "MAIC", criterion_scale = TRUE, show_progress = TRUE,
                     do_parallel = TRUE, cores = NULL){

  inputs <- do_tests_and_bootstrap(data = data, boot_sqt_test = TRUE, boot_ur_test = FALSE,
                                   level = SQT_level, bootstrap = bootstrap, B = B,
                                   block_length = block_length, ar_AWB = ar_AWB, union = union,
                                   min_lag = min_lag, max_lag = max_lag, criterion = criterion,
                                   deterministics = deterministics, detrend = detrend,
                                   criterion_scale = criterion_scale, steps = steps, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel,
                                   cores = cores, data_name = data_name)

  if (is.null(data_name)) {
    data_name <- deparse(substitute(data))
  }
  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    if (NCOL(data) > 1) {
      var_names <- paste0(data_name, " var", 1:NCOL(data))
    } else {
      var_names <- data_name
    }
  }

  spec <- list("steps" = steps, "bootstrap" = bootstrap, "B" = B, "block_length" = inputs$inputs$l,
               "ar_AWB" = inputs$inputs$ar_AWB, "SQT_level" = SQT_level, "union" = union,
               "deterministics" = inputs$inputs$deterministics,
               "detrend" = inputs$inputs$detrend, "min_lag" = min_lag,
               "max_lag" = inputs$inputs$p_max, "criterion" = inputs$inputs$criterion,
               "criterion_scale" = inputs$inputs$criterion_scale, "mult_test_ctrl" = "SQT")

  if (union) { # Union Tests
    BSQTout <- BSQT_cpp(pvec = inputs$p_vec, test_i = inputs$test_stats,
                        t_star = inputs$test_stats_star, level = inputs$level)
    estimates <- rep(NA, NCOL(data))
    tstats <- drop(inputs$test_stats)
    method_name <- paste0(bootstrap, " bootstrap sequential quantile union test")

    details <- list("individual estimates" = t(inputs$indiv_par_est),
                    "individual statistics" = t(inputs$indiv_test_stats),
                    "individual p-values" = inputs$indiv_pval,
                    "selected lags" = t(inputs$indiv_lags),
                    "txt_null" = "Series has a unit root",
                    "txt_alternative" = "Series is stationary")
    rownames(details$"individual estimates") <- var_names
    colnames(details$"individual estimates") <- c(t(outer(c("OLS", "QD"),
                                                  c("intercept", "intercept and trend"),
                                                  function(x,y){paste0(x, "/", y)})))
    rownames(details$"individual statistics") <- var_names
    colnames(details$"individual statistics") <- c(t(outer(c("OLS", "QD"),
                                                   c("intercept", "intercept and trend"),
                                                   function(x,y){paste0(x, "/", y)})))
    rownames(details$"individual p-values") <- var_names
    colnames(details$"individual p-values") <- c(t(outer(c("OLS", "QD"),
                                                 c("intercept", "intercept and trend"),
                                                 function(x,y){paste0(x, "/", y)})))
    rownames(details$"selected lags") <- var_names
    colnames(details$"selected lags") <- c(t(outer(c("OLS", "QD"),
                                           c("intercept", "intercept and trend"),
                                           function(x,y){paste0(x, "/", y)})))
  } else { # No Union Tests
    BSQTout <- BSQT_cpp(pvec = inputs$p_vec,
                        test_i = matrix(inputs$indiv_test_stats[1, ], nrow = 1),
                        t_star = inputs$t_star[ , 1,], level = inputs$level)
    estimates <- t(inputs$indiv_par_est)
    tstats <- drop(inputs$indiv_test_stats)
    method_name <- paste0(bootstrap, " bootstrap sequential quantile ",
                          inputs$inputs$name, " test (with " ,
                          inputs$inputs$deterministics,")")

    details <- list("individual estimates" = t(inputs$indiv_par_est),
                    "individual statistics" = t(inputs$indiv_test_stats),
                    "individual p-values" = inputs$indiv_pval,
                    "selected lags" = t(inputs$indiv_lags),
                    "txt_null" = "Series has a unit root",
                    "txt_alternative" = "Series is stationary")
    rownames(details$"individual estimates") <- var_names
    colnames(details$"individual estimates") <- paste0(inputs$inputs$detrend, "/",
                                                       inputs$inputs$deterministics)
    rownames(details$"individual statistics") <- var_names
    colnames(details$"individual statistics") <- paste0(inputs$inputs$detrend, "/",
                                                        inputs$inputs$deterministics)
    rownames(details$"individual p-values") <- var_names
    colnames(details$"individual p-values") <- paste0(inputs$inputs$detrend, "/",
                                                      inputs$inputs$deterministics)
    rownames(details$"selected lags") <- var_names
    colnames(details$"selected lags") <- paste0(inputs$inputs$detrend, "/",
                                                inputs$inputs$deterministics)
  }
  rej_H0 <- matrix(BSQTout$rej_H0 == 1, nrow = NCOL(data))
  rownames(rej_H0) <- var_names
  colnames(rej_H0) <- "Reject null"
  BSQT_seq <- BSQTout$BSQT_steps[, -3, drop = FALSE]
  rownames(BSQT_seq) <- paste("Step", 1:nrow(BSQT_seq))
  colnames(BSQT_seq) <- c("H0: # I(0)", "H1: # I(0)", "tstat", "p-value")
  p_vals <- rep(NA, NCOL(data))
  names(estimates) <- names(tstats) <- names(p_vals) <- var_names
  details$SQT <- BSQT_seq

  sqt_output <- list(method = method_name, data.name = data_name,
                     null.value =  c("gamma" = 0), alternative = "less",
                     estimate = estimates, statistic = tstats, p.value = p_vals,
                     rejections = rej_H0, details = details, series.names = var_names,
                     specifications = spec)
  class(sqt_output) <- c("bootUR", "mult_htest")
  return(sqt_output)
}

#' Panel Unit Root Test
#' @description Performs a test on a multivariate (panel) time series by testing the null hypothesis that all series have a unit root. The test is based on averaging the individual test statistics, also called the Group-Mean (GM) test in Palm, Smeekes and Urbain (2011).
#' @inheritParams boot_ur
#' @param union Logical indicator whether or not to use bootstrap union tests (\code{TRUE}) or not (\code{FALSE}), see Smeekes and Taylor (2012). Default is \code{TRUE}.
#' @details See \code{\link{boot_ur}} for details on the bootstrap algorithm and lag selection.
#' @return An object of class \code{"bootUR"}, \code{"htest"} with the following components:
#' \item{\code{method}}{The name of the hypothesis test method;}
#' \item{\code{data.name}}{The name of the variable on which the method is performed;}
#' \item{\code{null.value}}{The value of the (gamma) parameter of the lagged dependent variable in the ADF regression under the null hypothesis. Under the null, the series has a unit root. Testing the null of a unit root then boils down to testing the significance of the gamma parameter;}
#' \item{\code{alternative}}{A character string specifying the direction of the alternative hypothesis relative to the null value. The alternative postulates that the series is stationary;}
#' \item{\code{estimate}}{For the union test, the estimated value of the (gamma) parameter of the lagged dependent variable in the ADF regression is not defined, hence NA is given;}
#' \item{\code{statistic}}{The value of the test statistic of the unit root test;}
#' \item{\code{p.value}}{The p-value of the unit root test;}
#' \item{\code{details}}{A list containing the detailed outcomes of the performed tests, such as selected lags, individual estimates and p-values.}
#' \item{\code{specifications}}{The specifications used in the test.}
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
#' two_series_boot_panel <- boot_panel(MacroTS[, 1:2], bootstrap = "AWB", B = 199,
#'                                     do_parallel = FALSE, show_progress = FALSE)
#' print(two_series_boot_panel)
#' @export
boot_panel <- function(data, data_name = NULL, bootstrap = "AWB", B = 1999,
                       block_length = NULL, ar_AWB = NULL, union = TRUE, union_quantile = 0.05,
                       deterministics = NULL, detrend = NULL, min_lag = 0, max_lag = NULL,
                       criterion = "MAIC", criterion_scale = TRUE, show_progress = TRUE,
                       do_parallel = TRUE, cores = NULL){

  inputs <- do_tests_and_bootstrap(data = data, boot_sqt_test = FALSE, boot_ur_test = FALSE,
                                   level = union_quantile, bootstrap = bootstrap, B = B,
                                   block_length = block_length, ar_AWB = ar_AWB, union = union,
                                   min_lag = min_lag, max_lag = max_lag, criterion = criterion,
                                   deterministics = deterministics, detrend = detrend,
                                   criterion_scale = criterion_scale, steps = NULL, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel,
                                   cores = cores, data_name = data_name)

  if (is.null(data_name)) {
    data_name <- deparse(substitute(data))
  }
  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    if (NCOL(data) > 1) {
      var_names <- paste0(data_name, " var", 1:NCOL(data))
    } else {
      var_names <- data_name
    }
  }

  spec <- list("bootstrap" = bootstrap, "B" = B, "block_length" = inputs$inputs$l,
               "ar_AWB" = inputs$inputs$ar_AWB, "union" = union,
               "union_quantile" = inputs$inputs$union_quantile,
               "deterministics" = inputs$inputs$deterministics,
               "detrend" = inputs$inputs$detrend, "min_lag" = min_lag,
               "max_lag" = inputs$inputs$p_max, "criterion" = inputs$inputs$criterion,
               "criterion_scale" = inputs$inputs$criterion_scale, "mult_test_ctrl" = "none")

  if (union) { # Union Test
    GM_test <- mean(inputs$test_stats)
    t_star <- rowMeans(inputs$test_stats_star)
    p_val <- mean(t_star < GM_test)
    method_name <- paste0("Panel ", bootstrap, " bootstrap group-mean union test")

    details <- list("individual estimates" = t(inputs$indiv_par_est),
                    "individual statistics" = t(inputs$indiv_test_stats),
                    "individual p-values" = inputs$indiv_pval,
                    "selected lags" = t(inputs$indiv_lags),
                    "txt_null" = "All series have a unit root",
                    "txt_alternative" = "Some series are stationary")
    rownames(details$"individual estimates") <- var_names
    colnames(details$"individual estimates") <- c(t(outer(c("OLS", "QD"),
                                                  c("intercept", "intercept and trend"),
                                                  function(x,y){paste0(x, "/", y)})))
    rownames(details$"individual statistics") <- var_names
    colnames(details$"individual statistics") <- c(t(outer(c("OLS", "QD"),
                                                   c("intercept", "intercept and trend"),
                                                   function(x,y){paste0(x, "/", y)})))
    rownames(details$"individual p-values") <- var_names
    colnames(details$"individual p-values") <- c(t(outer(c("OLS", "QD"),
                                                 c("intercept", "intercept and trend"),
                                                 function(x,y){paste0(x, "/", y)})))
    rownames(details$"selected lags") <- var_names
    colnames(details$"selected lags") <- c(t(outer(c("OLS", "QD"),
                                           c("intercept", "intercept and trend"),
                                           function(x,y){paste0(x, "/", y)})))
  } else { # No Union Test
    GM_test <- rowMeans(inputs$tests_i)
    t_star <- apply(inputs$t_star, 1:2, mean)
    p_val <- sapply(1, function(i){mean(t_star[, i] < GM_test[i])})
    method_name <- paste0("Panel" , bootstrap, " bootstrap group-mean ", inputs$inputs$name,
                          " test (with " , inputs$inputs$deterministics,")")

    details <- list("individual estimates" = t(inputs$indiv_par_est),
                    "individual statistics" = t(inputs$indiv_test_stats),
                    "individual p-values" = inputs$indiv_pval,
                    "selected lags" = t(inputs$indiv_lags),
                    "txt_null" = "All series have a unit root",
                    "txt_alternative" = "Some series are stationary")
    rownames(details$"individual estimates") <- var_names
    colnames(details$"individual estimates") <- paste0(inputs$inputs$detrend, "/",
                                                inputs$inputs$deterministics)
    rownames(details$"individual statistics") <- var_names
    colnames(details$"individual statistics") <- paste0(inputs$inputs$detrend, "/",
                                                 inputs$inputs$deterministics)
    rownames(details$"individual p-values") <- var_names
    colnames(details$"individual p-values") <- paste0(inputs$inputs$detrend, "/",
                                           inputs$inputs$deterministics)
    rownames(details$"selected lags") <- var_names
    colnames(details$"selected lags") <- paste0(inputs$inputs$detrend, "/",
                                           inputs$inputs$deterministics)
  }

  attr(GM_test, "names") <- "tstat"
  gamma_hat <- NA
  attr(gamma_hat, "names") <- "gamma"
  attr(p_val, "names") <- "p-value"
  panel_output <- list(method = method_name, data.name = data_name,
                       null.value = c("gamma" = 0), alternative = "less",
                       estimate = gamma_hat, statistic = GM_test, p.value = p_val,
                       details = details, series.names = data_name, specifications = spec)
  class(panel_output) <- c("bootUR", "htest")
  return(panel_output)
}
