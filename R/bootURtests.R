#' Individual Unit Root Tests without multiple testing control
#' @description This function performs bootstrap unit root tests on each time series individually.
#' @param y A \eqn{T}-dimensional vector or a (\eqn{T} x \eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations to be tested for unit roots. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}), or a data frame, as long as each column represents a single time series.
#' @param level Desired significance level of the unit root test. Default is 0.05.
#' @param boot String for bootstrap method to be used. Options are
#' \describe{
#' \item{\code{"MBB"}}{Moving blocks bootstrap (Paparoditis and Politis, 2003; Palm, Smeekes and Urbain, 2011);}
#' \item{\code{"BWB"}}{Block wild bootstrap (Shao, 2011; Smeekes and Urbain, 2014a);}
#' \item{\code{"DWB"}}{Dependent wild bootstrap (Shao, 2010; Smeekes and Urbain, 2014a; Rho and Shao, 2019);}
#' \item{\code{"AWB"}}{Autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2020), this is the default;}
#' \item{\code{"SB"}}{Sieve bootstrap (Chang and Park, 2003; Palm, Smeekes and Urbain, 2008; Smeekes, 2013);}
#' \item{\code{"SWB"}}{Sieve wild boostrap (Cavaliere and Taylor, 2009; Smeekes and Taylor, 2012).}
#' }
#' @param B Number of bootstrap replications. Default is 1999.
#' @param l Desired 'block length' in the bootstrap. For the MBB, BWB and DWB boostrap, this is a genuine block length. For the AWB boostrap, the block length is transformed into an autoregressive parameter via the formula \eqn{0.01^(1/l)} as in Smeekes and Urbain (2014a); this can be overwritten by setting \code{ar_AWB} directly. Default sets the block length as a function of the time series length T, via the rule \eqn{l = 1.75 T^(1/3)} of Palm, Smeekes and Urbain (2011).
#' @param ar_AWB Autoregressive parameter used in the AWB bootstrap method (\code{boot = "AWB"}). Can be used to set the parameter directly rather than via the default link to the block length l.
#' @param union Logical indicator whether or not to use bootstrap union tests (\code{TRUE}) or not (\code{FALSE}), see Smeekes and Taylor (2012). Default is \code{TRUE}.
#' @param p_min Minimum lag length in the augmented Dickey-Fuller regression. Default is 0.
#' @param p_max Maximum lag length in the augmented Dickey-Fuller regression. Default uses the sample size-based rule \eqn{12(T/100)^{1/4}}.
#' @param ic String for information criterion used to select the lag length in the augmented Dickey-Fuller regression. Options are: \code{"AIC"}, \code{"BIC"}, \code{"MAIC"}, \code{"MBIC"}. Default is \code{"MAIC"} (Ng and Perron, 2001).
#' @param dc Numeric vector indicating the deterministic specification. Only relevant if \code{union = FALSE}. Options are (combinations of)
#'
#' \verb{0 } no deterministics;
#'
#' \verb{1 } intercept only;
#'
#' \verb{2 } intercept and trend.
#'
#' If \code{union = FALSE}, the default is adding an intercept (a warning is given).
#' @param detr String vector indicating the type of detrending to be performed. Only relevant if \code{union = FALSE}. Options are: \code{"OLS"} and/or \code{"QD"} (typically also called GLS, see Elliott, Rothenberg and Stock, 1996). The default is \code{"OLS"}.
#' @param ic_scale Logical indicator whether or not to use the rescaled information criteria of Cavaliere et al. (2015) (\code{TRUE}) or not (\code{FALSE}). Default is \code{TRUE}.
#' @param verbose Logical indicator whether or not information on the outcome of the unit root test needs to be printed to the console. Default is \code{FALSE}.
#' @param show_progress Logical indicator whether a bootstrap progress update should be printed to the console. Default is FALSE.
#' @param do_parallel Logical indicator whether bootstrap loop should be executed in parallel. Parallel computing is only available if OpenMP can be used, if not this option is ignored. Default is FALSE.
#' @param nc The number of cores to be used in the parallel loops. Default is to use all but one.
#' @details The options encompass many test proposed in the literature. \code{dc = "OLS"} gives the standard augmented Dickey-Fuller test, while \code{dc = "QD"} provides the DF-GLS test of Elliott, Rothenberg and Stock (1996). The bootstrap algorithm is always based on a residual bootstrap (under the alternative) to obtain residuals rather than a difference-based bootstrap (under the null), see e.g. Palm, Smeekes and Urbain (2008).
#'
#' Lag length selection is done automatically in the ADF regression with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. For very short time series (fewer than 50 time points) the maximum lag length is adjusted downward to avoid potential multicollinearity issues in the bootstrap. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `p_min` and maximum lag length `p_max` for the selection algorithm equal to the desired lag length.
#'
#' @export
#' @return A list with the following components
#' \item{\code{rej_H0}}{Logical indicator whether the null hypothesis of a unit root is rejected (\code{TRUE}) or not (\code{FALSE});}
#' \item{\code{ADF_tests}}{Details on the unit root tests: value of the test statistics and p-values.}
#' For the union test (\code{union = TRUE}), the output is arranged per time series. If \code{union = FALSE}, the output is arranged per time series, type of deterministic component (\code{dc}) and detrending method (\code{detr}).
#' @section Warnings:
#' The function may give the following warnings.
#' \describe{
#' \item{\code{Warning: Missing values cause resampling bootstrap to be executed for each time series individually.}}{If the time series in \code{y} have different starting and end points (and thus some series contain \code{NA} values at the beginning and/or end of the sample, the resampling-based moving block bootstrap (MBB) and sieve bootstrap (SB) cannot be used directly, as they create holes (internal missings) in the bootstrap samples. These bootstrap methods are therefore not applied jointly as usual, but individually to each series.}
#' \item{\code{Warning: Deterministic specification in argument dc is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting deterministic components manually therefore has no effect.}
#' \item{\code{Warning: Detrending method in argument detr is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting detrending methods manually therefore has no effect.}
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
#' # iADFtest on GDP_BE and GDP_DE
#' two_series_iADFtest <- iADFtest(MacroTS[, 1:2], boot = "MBB", B = 399,
#' verbose = TRUE)
iADFtest <- function(y, level = 0.05, boot = "AWB", B = 1999, l = NULL,
                     ar_AWB = NULL, union = TRUE, p_min = 0, p_max = NULL,
                     ic = "MAIC", dc = NULL, detr = NULL, ic_scale = TRUE,
                     verbose = FALSE, show_progress = FALSE,
                     do_parallel = FALSE, nc = NULL){

  inputs <- do_tests_and_bootstrap(y = y, BSQT_test = FALSE, iADF_test = TRUE, level = level,
                                   boot = boot, B = B, l = l, ar_AWB = ar_AWB, union = union,
                                   p_min = p_min, p_max = p_max, ic = ic, dc = dc, detr = detr,
                                   ic_scale = ic_scale, q = NULL, h_rs = 0.1,
                                   show_progress = show_progress,
                                   do_parallel = do_parallel, nc = nc)

  if (!is.null(colnames(y))) {
    var_names <- colnames(y)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(y))
  }

  if (union) { # Union Tests
    iADFout <- iADF_cpp(test_i = inputs$test_stats, t_star = inputs$test_stats_star,
                        level = inputs$level)
    rej_H0 <- (iADFout[, 2] < level)
    colnames(iADFout) <- c("test statistic", "p-value")
    rownames(iADFout) <- var_names

    if (verbose) {
      p_hat <- sum(rej_H0)
      if (NCOL(y) > 1) {
        if (p_hat > 1) {
          cat(paste("There are ", p_hat, " stationary time series, namely: ",
                    paste(var_names[rej_H0], collapse = " "), "\n", sep = ""))
        } else {
          if (p_hat == 1) {
            cat(paste("There is ", p_hat, " stationary time series, namely: ",
                      paste(var_names[rej_H0], collapse = " "), "\n", sep = ""))
          } else {
            cat(paste("There are ", p_hat, " stationary time series.\n", sep = ""))
          }
        }
        print(iADFout)
      } else {
        cat("Bootstrap Union Test:\n")
        if (rej_H0) {
          cat(paste("The null hypothesis of a unit root is rejected at a significance
                    level of ", level, ".\n", sep = ""))
        } else {
          cat(paste("The null hypothesis of a unit root is not rejected at a significance
                    level of ", level, ".\n", sep = ""))
        }
        print(iADFout[1, ])
      }
    }
  } else { # No Union Tests
    detr_names <- rep(NA, length(inputs$detr))
    detr_names[inputs$detr=="OLS"] <- c("detr = OLS")
    detr_names[inputs$detr=="QD"] <- c("detr = QD")

    dc_names <- rep(NA, length(inputs$dc))
    dc_names[inputs$dc==0] <- c("dc = none")
    dc_names[inputs$dc==1] <- c("dc = intercept")
    dc_names[inputs$dc==2] <- c("dc = intercept and trend")

    detr_dc_names <- paste(rep(detr_names, each = length(inputs$dc)), ", ",
                           rep(dc_names, length(inputs$detr)), sep = "")
    iADFout <- array(NA, c(NCOL(y), 2, length(detr_dc_names)), dimnames =
                      list(var_names, c("test statistic", "p-value"), detr_dc_names))
    rej_H0 <- array(NA, c(NCOL(y), length(inputs$dc)*length(inputs$detr)), dimnames =
                      list(var_names, detr_dc_names))
    for(i in 1:nrow(inputs$tests_i)){
      iADFout[, , i] <- iADF_cpp(test_i = matrix(inputs$tests_i[i, ], nrow = 1),
                                 t_star = matrix(inputs$t_star[ , i, ], nrow = B),
                                 level = inputs$level)
      rej_H0[, i] <- (iADFout[, 2, i] < level)

      if (verbose) {
        cat(paste(c(rep("-", 40), "\n"), sep = "", collapse = ""))
        cat(paste("Type of unit root test performed: ", detr_dc_names[i], "\n", sep = ""))
        if (NCOL(y) > 1) {
          p_hat <- sum(rej_H0[, i])
          if (p_hat > 1) {
            cat(paste("There are ", p_hat, " stationary time series, namely: ",
                      paste(var_names[rej_H0[, i]], collapse = " "), ".\n", sep = ""))
          } else {
            if (p_hat == 1) {
              cat(paste("There is ", p_hat, " stationary time series, namely: ",
                        paste(var_names[rej_H0[, i]], collapse = " "), ".\n", sep = ""))
            } else {
              cat(paste("There are ", p_hat, " stationary time series", "\n", sep = ""))
            }
          }
          print(iADFout[, , i])
        } else {
          print(iADFout[, , i])
        }
      }
    }
  }
  return(list(rej_H0 = rej_H0, ADF_tests = iADFout))
}

#' Bootstrap augmented Dickey-Fuller Unit Root Test
#' @description This function performs a standard augmented Dickey-Fuller bootstrap unit root test on a single time series.
#' @inheritParams iADFtest
#' @param y A \eqn{T}-dimensional vector to be tested for unit roots. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}), or a data frame.
#' @param boot String for bootstrap method to be used. Options are
#' \describe{
#' \item{\code{"MBB"}}{Moving blocks bootstrap (Paparoditis and Politis, 2003);}
#' \item{\code{"BWB"}}{Block wild bootstrap (Shao, 2011);}
#' \item{\code{"DWB"}}{Dependent wild bootstrap (Shao, 2010; Rho and Shao, 2019);}
#' \item{\code{"AWB"}}{Autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2020), this is the default;}
#' \item{\code{"SB"}}{Sieve bootstrap (Chang and Park, 2003; Palm, Smeekes and Urbain, 2008; Smeekes, 2013);}
#' \item{\code{"SWB"}}{Sieve wild boostrap (Cavaliere and Taylor, 2009; Smeekes and Taylor, 2012).}
#' }
#' @details The options encompass many test proposed in the literature. \code{dc = "OLS"} gives the standard augmented Dickey-Fuller test, while \code{dc = "QD"} provides the DF-GLS test of Elliott, Rothenberg and Stock (1996). The bootstrap algorithm is always based on a residual bootstrap (under the alternative) to obtain residuals rather than a difference-based bootstrap (under the null), see e.g. Palm, Smeekes and Urbain (2008).
#'
#' Lag length selection is done automatically in the ADF regression with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. For very short time series (fewer than 50 time points) the maximum lag length is adjusted downward to avoid potential multicollinearity issues in the bootstrap. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `p_min` and maximum lag length `p_max` for the selection algorithm equal to the desired lag length.
#' @export
#' @return Values of the Dickey-Fuller test statistics and corresponding bootstrap p-values.
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Multiple time series not allowed. Switch to a multivariate method such as iADFtest, or change argument y to a univariate time series.}}{The function is a simple wrapper around \code{\link{iADFtest}} to facilitate use for single time series. It does not support multiple time series, as \code{\link{iADFtest}} is specifically suited for that.}
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
#' @seealso \code{\link{iADFtest}}
#' @examples
#' # boot_df on GDP_BE
#' GDP_BE_df <- boot_df(MacroTS[, 1], B = 399, dc = 2, detr = "OLS", verbose = TRUE)
boot_df <- function(y, level = 0.05, boot = "AWB", B = 1999, l = NULL, ar_AWB = NULL,
                    p_min = 0, p_max = NULL, ic = "MAIC", dc = 1, detr = "OLS",
                    ic_scale = TRUE, verbose = FALSE, show_progress = FALSE,
                    do_parallel = FALSE, nc = NULL){
  if (NCOL(y) > 1) {
    stop("Multiple time series not allowed. Switch to a multivariate method such as iADFtest,
         or change argument y to a univariate time series.")
  }

  if (verbose) {
    cat("Bootstrap DF Test with", boot, "bootstrap method.\n")
  }
  out <- iADFtest(y, level = level, boot = boot, B = B, l = l, ar_AWB = ar_AWB,
                  union = FALSE, p_min = p_min, p_max = p_max, ic = ic, dc = dc,
                  detr = detr, ic_scale = ic_scale, verbose = verbose,
                  show_progress = show_progress, do_parallel = do_parallel, nc = nc)
  return(aperm(out$ADF_tests, 3:1)[, , 1])
}

#' Bootstrap Union Test for Unit Roots
#' @description Performs bootstrap unit root test based on the union of rejections of 4 tests with different number of deterministic components and different type of detrending (Harvey, Leybourne and Taylor, 2012; Smeekes and Taylor, 2012).
#' @inheritParams iADFtest
#' @param y A \eqn{T}-dimensional vector to be tested for unit roots. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}), or a data frame.
#' @param boot String for bootstrap method to be used. Options are
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
#' Lag length selection is done automatically in the ADF regressions with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `p_min` and maximum lag length `p_max` for the selection algorithm equal to the desired lag length.
#' @export
#' @return Value of the union test statistic and the bootstrap p-values.
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Multiple time series not allowed. Switch to a multivariate method such as iADFtest, or change argument y to a univariate time series.}}{The function is a simple wrapper around \code{\link{iADFtest}} to facilitate use for single time series. It does not support multiple time series, as \code{\link{iADFtest}} is specifically suited for that.}
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
#' @seealso \code{\link{iADFtest}}
#' @examples
#' # boot_union on GDP_BE
#' GDP_BE_df <- boot_union(MacroTS[, 1], B = 399, verbose = TRUE)
boot_union <- function(y, level = 0.05, boot = "AWB", B = 1999, l = NULL, ar_AWB = NULL,
                       p_min = 0, p_max = NULL, ic = "MAIC", ic_scale = TRUE, verbose = FALSE,
                       show_progress = FALSE, do_parallel = FALSE, nc = NULL){

  if (verbose) {
    cat("Bootstrap Test with", boot, "bootstrap method.\n")
  }
  if (NCOL(y) > 1) {
    stop("Multiple time series not allowed. Switch to a multivariate method such as iADFtest,
         or change argument y to a univariate time series.")
  }
  out <- iADFtest(y, level = level, boot = boot, B = B, l = l, ar_AWB = ar_AWB,
                  union = TRUE, p_min = p_min, p_max = p_max, ic = ic, ic_scale = ic_scale,
                  verbose = verbose, show_progress = show_progress,
                  do_parallel = do_parallel, nc = nc)
  return(out$ADF_tests[1, ])
}

#' Bootstrap Unit Root Tests with False Discovery Rate control
#' @description Controls for multiple testing by controlling the false discovery rate (FDR), see Moon and Perron (2012) and Romano, Shaikh and Wolf (2008).
#' @inheritParams iADFtest
#' @param level Desired False Discovery Rate level of the unit root tests. Default is 0.05.
#' @details The false discovery rate FDR is defined as the expected proportion of false rejections relative to the total number of rejections.
#'
#' See \code{\link{iADFtest}} for details on the bootstrap algorithm and lag selection.
#' @export
#' @return A list with the following components
#' \item{\code{rej_H0}}{Logical indicator whether the null hypothesis of a unit root is rejected (\code{TRUE}) or not (\code{FALSE});}
#' \item{\code{FDR_sequence}}{Details on the unit root tests: value of the test statistics and critical values.}
#' For the union test (\code{union = TRUE}), the output is arranged per time series. If \code{union = FALSE}, the output is arranged per time series, type of deterministic component (\code{dc}) and detrending method (\code{detr}).
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Resampling-based bootstraps MBB and SB cannot handle missing values.}}{If the time series in \code{y} have different starting and end points (and thus some series contain \code{NA} values at the beginning and/or end of the sample, the resampling-based moving block bootstrap (MBB) and sieve bootstrap (SB) cannot be used, as they create holes (internal missings) in the bootstrap samples. Switch to another bootstrap method or truncate your sample to eliminate \code{NA} values.}
#' \item{\code{Warning: SB and SWB bootstrap only recommended for iADFtest; see help for details.}}{Although the sieve bootstrap methods \code{"SB"} and \code{"SWB"} can be used, Smeekes and Urbain (2014b) show that these are not suited to capture general forms of dependence across units, and using them for joint or multiple testing is not valid. This warning thereofre serves to recommend the user to consider a different bootstrap method.}
#' \item{\code{Warning: Deterministic specification in argument dc is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting deterministic components manually therefore has no effect.}
#' \item{\code{Warning: Detrending method in argument detr is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting detrending methods manually therefore has no effect.}
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
#' @seealso \code{\link{iADFtest}}
#' @examples
#' # bFDRtest on GDP_BE and GDP_DE
#' two_series_bFDRtest <- bFDRtest(MacroTS[, 1:2], boot = "MBB", B = 399,  verbose = TRUE)
bFDRtest <- function(y, level = 0.05,  boot = "AWB", B = 1999, l = NULL, ar_AWB = NULL,
                     union = TRUE, p_min = 0, p_max = NULL, ic = "MAIC", dc = NULL,
                     detr = NULL, ic_scale = TRUE, verbose = FALSE, show_progress = FALSE,
                     do_parallel = FALSE, nc = NULL){

  inputs <- do_tests_and_bootstrap(y = y, BSQT_test = FALSE, iADF_test = FALSE, level = level,
                                   boot = boot, B = B, l = l, ar_AWB = ar_AWB, union = union,
                                   p_min = p_min, p_max = p_max, ic = ic, dc = dc, detr = detr,
                                   ic_scale = ic_scale, q = NULL, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel, nc = nc)

  if (!is.null(colnames(y))) {
    var_names <- colnames(y)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(y))
  }

  if (union) { # Union Tests
    bFDRout <- FDR_cpp(test_i = inputs$test_stats, t_star = inputs$test_stats_star,
                       level = inputs$level)
    rej_H0 <- matrix(bFDRout$rej_H0 == 1, nrow = NCOL(y))
    FDR_seq <- bFDRout$FDR_Tests[, -1, drop = FALSE]

    rownames(rej_H0) <- var_names
    rownames(FDR_seq) <- var_names[bFDRout$FDR_Tests[, 1, drop = FALSE]]

    colnames(rej_H0) <- "Reject H0"
    colnames(FDR_seq) <- c("test statistic", "critical value")

    if (verbose) {
      p_hat <- sum(rej_H0)
      if (p_hat > 1){
        cat(paste("There are ", p_hat, " stationary time series, namely: ",
                  paste(var_names[rej_H0], collapse = " "), "\n", sep = ""))
      } else if (p_hat == 1) {
          cat(paste("There is ", p_hat, " stationary time series, namely: ",
                    paste(var_names[rej_H0], collapse = " "), "\n", sep = ""))
      } else {
          cat(paste("There are ", p_hat, " stationary time series", "\n", sep = ""))
      }
      cat("Details of the FDR sequential tests:\n")
      print(FDR_seq)
    }
  } else { # No Union Tests
    detr_names <- rep(NA, length(inputs$detr))
    detr_names[inputs$detr=="OLS"] <- c("detr = OLS")
    detr_names[inputs$detr=="QD"] <- c("detr = QD")

    dc_names <- rep(NA, length(inputs$dc))
    dc_names[inputs$dc==0] <- c("dc = none")
    dc_names[inputs$dc==1] <- c("dc = intercept")
    dc_names[inputs$dc==2] <- c("dc = intercept and trend")

    detr_dc_names <- paste(rep(detr_names, each = length(inputs$dc)), ", ",
                           rep(dc_names, length(inputs$detr)), sep = "")
    rej_H0 <- matrix(nrow = NCOL(y), ncol = length(inputs$dc)*length(inputs$detr))
    rownames(rej_H0) <- var_names
    colnames(rej_H0) <- detr_dc_names
    FDR_seq <- vector("list", length(inputs$dc)*length(inputs$detr))
    names(FDR_seq) <- detr_dc_names
    for (i in 1:nrow(inputs$tests_i)) {
      bFDRout <- FDR_cpp(test_i = matrix(inputs$tests_i[i, ], nrow = 1),
                         t_star = inputs$t_star[ , i,], level = inputs$level)
      rej_H0[, i] <- bFDRout$rej_H0 == 1
      FDR_seq[[i]] <- bFDRout$FDR_Tests[, -1, drop = FALSE]
      rownames(FDR_seq[[i]]) <- var_names[bFDRout$FDR_Tests[, 1, drop = FALSE]]
      colnames(FDR_seq[[i]]) <- c("test statistic", "critical value")

      if (verbose) {
        p_hat <- sum(rej_H0[, i])
        cat(paste(c(rep("-", 40), "\n"), sep = "", collapse = ""))
        cat(paste("Type of unit root test performed: ", detr_dc_names[i], "\n", sep = ""))
        if (p_hat > 1){
          cat(paste("There are ", p_hat, " stationary time series, namely: ",
                    paste(var_names[rej_H0[, i]], collapse = " "), ".\n", sep = ""))
        } else if (p_hat == 1) {
          cat(paste("There is ", p_hat, " stationary time series, namely: ",
                    paste(var_names[rej_H0[, i]], collapse = " "), ".\n", sep = ""))
        } else {
          cat(paste("There are ", p_hat, " stationary time series.", "\n", sep = ""))
        }
        cat("Details of the FDR sequential tests:\n")
        print(FDR_seq[[i]])
      }
    }
  }
  return(list(rej_H0 = rej_H0, FDR_sequence = FDR_seq))
}


#' Bootstrap Sequential Quantile Test
#' @description Performs the Bootstrap Sequential Quantile Test (BSQT) proposed by Smeekes (2015).
#' @inheritParams iADFtest
#' @param q Numeric vector of quantiles or units to be tested. Default is to test each unit sequentially.
#' @details The parameter \code{q} can either be set as an increasing sequence of integers smaller or equal to the number of series \code{N}, or fractions of the total number of series (quantiles). For \code{N} time series, setting \code{q = 0:N} means each unit should be tested sequentially. In this case the method is equivalent to the StepM method of Romano and Wolf (2005), and therefore controls the familywise error rate. To split the series in \code{K} equally sized groups, use \code{q = 0:K / K}.
#'
#' By convention and in accordance with notation in Smeekes (2015), the first entry of the vector should be equal to zero, while the second entry indicates the end of the first group, and so on. If the initial \code{0} or final value (\code{1} or \code{N}) are omitted, they are automatically added by the function.
#'
#' See \code{\link{iADFtest}} for details on the bootstrap algorithm and lag selection.
#' @export
#' @return A list with the following components
#' \item{\code{rej_H0}}{Logical indicator whether the null hypothesis of a unit root is rejected (\code{TRUE}) or not (\code{FALSE});}
#' \item{\code{BSQT_sequence}}{Details on the unit root tests: outcome of the sequential steps, value of the test statistics and p-values.}
#' For the union test (\code{union = TRUE}), the output is arranged per time series. If \code{union = FALSE}, the output is arranged per time series, type of deterministic component (\code{dc}) and detrending method (\code{detr}).
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Resampling-based bootstraps MBB and SB cannot handle missing values.}}{If the time series in \code{y} have different starting and end points (and thus some series contain \code{NA} values at the beginning and/or end of the sample, the resampling-based moving block bootstrap (MBB) and sieve bootstrap (SB) cannot be used, as they create holes (internal missings) in the bootstrap samples. Switch to another bootstrap method or truncate your sample to eliminate \code{NA} values.}
#' \item{\code{Error: Invalid input values for q: must be quantiles or positive integers.}}{Construction of \code{q} does not satisfy the criteria listed under 'Details'.}
#' \item{\code{Warning: SB and SWB bootstrap only recommended for iADFtest; see help for details.}}{Although the sieve bootstrap methods \code{"SB"} and \code{"SWB"} can be used, Smeekes and Urbain (2014b) show that these are not suited to capture general forms of dependence across units, and using them for joint or multiple testing is not valid. This warning thereofre serves to recommend the user to consider a different bootstrap method.}
#' \item{\code{Warning: Deterministic specification in argument dc is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting deterministic components manually therefore has no effect.}
#' \item{\code{Warning: Detrending method in argument detr is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting detrending methods manually therefore has no effect.}
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
#' @seealso \code{\link{iADFtest}}
#' @examples
#' # BSQTtest on GDP_BE and GDP_DE
#' two_series_BSQTtest <- BSQTtest(MacroTS[, 1:2], boot = "AWB", B = 399,  verbose = TRUE)
BSQTtest <- function(y, q = 0:NCOL(y), level = 0.05,  boot = "AWB", B = 1999, l = NULL,
                     ar_AWB = NULL, union = TRUE, p_min = 0, p_max = NULL, ic = "MAIC", dc = NULL,
                     detr = NULL, ic_scale = TRUE, verbose = FALSE, show_progress = FALSE,
                     do_parallel = FALSE, nc = NULL){

  inputs <- do_tests_and_bootstrap(y = y, BSQT_test = TRUE, iADF_test = FALSE, level = level,
                                   boot = boot, B = B, l = l, ar_AWB = ar_AWB, union = union,
                                   p_min = p_min, p_max = p_max, ic = ic, dc = dc, detr = detr,
                                   ic_scale = ic_scale, q = q, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel, nc = nc)

  if (!is.null(colnames(y))) {
    var_names <- colnames(y)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(y))
  }

  if (union) { # Union Tests
    BSQTout <- BSQT_cpp(pvec = inputs$p_vec, test_i = inputs$test_stats, t_star =
                          inputs$test_stats_star, level = inputs$level)
    rej_H0 <- matrix(BSQTout$rej_H0 == 1, nrow = NCOL(y))
    BSQT_seq <- BSQTout$BSQT_steps[, -3, drop = FALSE]

    rownames(rej_H0) <- var_names
    rownames(BSQT_seq) <- paste("Step", 1:nrow(BSQT_seq))
    colnames(rej_H0) <- "Reject H0"
    colnames(BSQT_seq) <- c("Unit H0", "Unit H1", "Test statistic", "p-value")

    if (verbose) {
      p_hat <- sum(rej_H0)
      if (p_hat > 1) {
        cat(paste("There are ", p_hat, " stationary time series, namely: ",
                  paste(var_names[rej_H0], collapse = " "), ".\n", sep = ""))
      } else {
        if (p_hat == 1) {
          cat(paste("There is ", p_hat, " stationary time series, namely: ",
                    paste(var_names[rej_H0], collapse = " "), ".\n", sep = ""))
        } else {
          cat(paste("There are ", p_hat, " stationary time series.\n", sep = ""))
        }
      }
      cat("Details of the BSQT sequential tests:\n")
      print(BSQT_seq)
    }
  } else { # No Union Tests
    detr_names <- rep(NA, length(inputs$detr))
    detr_names[inputs$detr=="OLS"] <- c("detr = OLS")
    detr_names[inputs$detr=="QD"] <- c("detr = QD")

    dc_names <- rep(NA, length(inputs$dc))
    dc_names[inputs$dc==0] <- c("dc = none")
    dc_names[inputs$dc==1] <- c("dc = intercept")
    dc_names[inputs$dc==2] <- c("dc = intercept and trend")

    detr_dc_names <- paste(rep(detr_names, each = length(inputs$dc)), ", ",
                           rep(dc_names, length(inputs$detr)), sep = "")
    rej_H0 <- matrix(nrow = NCOL(y), ncol = length(inputs$dc)*length(inputs$detr))
    rownames(rej_H0) <- var_names
    colnames(rej_H0) <- detr_dc_names
    BSQT_seq <- vector("list", length(inputs$dc)*length(inputs$detr))
    names(BSQT_seq) <- detr_dc_names
    for (i in 1:nrow(inputs$tests_i)) {
      BSQTout <- BSQT_cpp(pvec = inputs$p_vec, test_i = matrix(inputs$tests_i[i, ], nrow = 1),
                          t_star = inputs$t_star[ , i,], level = inputs$level)
      rej_H0[, i] <- BSQTout$rej_H0 == 1
      BSQT_seq[[i]] <- BSQTout$BSQT_steps[, -3, drop = FALSE]
      rownames(BSQT_seq[[i]]) <- paste("Step", 1:nrow(BSQT_seq[[i]]))
      colnames(BSQT_seq[[i]]) <- c("Unit H0", "Unit H1", "Test statistic", "p-value")

      if (verbose) {
        p_hat <- sum(rej_H0[, i])
        cat(paste(c(rep("-", 40), "\n"), sep = "", collapse = ""))
        cat(paste("Type of unit root test performed:", detr_dc_names[i], "\n", sep = ""))
        if (p_hat > 1){
          cat(paste("There are ", p_hat, "stationary time series, namely:",
                    paste(var_names[rej_H0[, i]], collapse = " "), ".\n", sep = ""))
        } else if (p_hat == 1) {
          cat(paste("There is ", p_hat, "stationary time series, namely:",
                    paste(var_names[rej_H0[, i]], collapse = " "), ".\n", sep = ""))
        } else {
          cat(paste("There are ", p_hat, "stationary time series.\n", sep = ""))
        }
        cat("Details of the BSQT sequential tests:\n")
        print(BSQT_seq[[i]])
      }
    }
  }
  return(list(rej_H0 = rej_H0, BSQT_sequence = BSQT_seq))
}

#' Panel Unit Root Test
#' @description Performs a test on a multivariate (panel) time series by testing the null hypothesis that all series have a unit root. The test is based on averaging the individual test statistics, also called the Group-Mean (GM) test in Palm, Smeekes and Urbain (2011).
#' @inheritParams iADFtest
#' @export
#' @details See \code{\link{iADFtest}} for details on the bootstrap algorithm and lag selection.
#'
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Resampling-based bootstraps MBB and SB cannot handle missing values.}}{If the time series in \code{y} have different starting and end points (and thus some series contain \code{NA} values at the beginning and/or end of the sample, the resampling-based moving block bootstrap (MBB) and sieve bootstrap (SB) cannot be used, as they create holes (internal missings) in the bootstrap samples. Switch to another bootstrap method or truncate your sample to eliminate \code{NA} values.}
#' \item{\code{Warning: SB and SWB bootstrap only recommended for iADFtest; see help for details.}}{Although the sieve bootstrap methods \code{"SB"} and \code{"SWB"} can be used, Smeekes and Urbain (2014b) show that these are not suited to capture general forms of dependence across units, and using them for joint or multiple testing is not valid. This warning thereofre serves to recommend the user to consider a different bootstrap method.}
#' \item{\code{Warning: Deterministic specification in argument dc is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting deterministic components manually therefore has no effect.}
#' \item{\code{Warning: Detrending method in argument detr is ignored, as union test is applied.}}{The union test calculates the union of all four combinations of deterministic components (intercept or intercept and trend) and detrending methods (OLS or QD). Setting detrending methods manually therefore has no effect.}
#' }
#' @return For the union test (\code{union = TRUE}), the test statistic and p-value are returned. If \code{union = FALSE}, the test statistics and p-values are reported per type of deterministic component (\code{dc}) and detrending method (\code{detr}).
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
#' @seealso \code{\link{iADFtest}}
#' @examples
#' # paneltest on GDP_BE and GDP_DE
#' two_series_paneltest <- paneltest(MacroTS[, 1:2], boot = "AWB", B = 399,  verbose = TRUE)
paneltest <- function(y, level = 0.05,  boot = "AWB", B = 1999, l = NULL, ar_AWB = NULL,
                      union = TRUE, p_min = 0, p_max = NULL, ic = "MAIC", dc = NULL, detr = NULL,
                      ic_scale = TRUE, verbose = FALSE, show_progress = FALSE,
                      do_parallel = FALSE, nc = NULL){

  inputs <- do_tests_and_bootstrap(y = y, BSQT_test = FALSE, iADF_test = FALSE, level = level,
                                   boot = boot, B = B, l = l, ar_AWB = ar_AWB, union = union,
                                   p_min = p_min, p_max = p_max, ic = ic, dc = dc, detr = detr,
                                   ic_scale = ic_scale, q = NULL, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel, nc = nc)

  if (union) { # Union Tests
    GM_test <- mean(inputs$test_stats)
    t_star <- rowMeans(inputs$test_stats_star)
    p_val <- mean(t_star < GM_test)
    out <- cbind("test statistic" = GM_test, "p-value" = p_val)

    if (verbose) {
      cat("Panel Bootstrap Group-Mean Union Test\n")
      if (p_val < level){
        cat(paste("The null hypothesis that all series have a unit root, is
                  rejected at a significance level of ", level, ".\n", sep = ""))
      } else {
        cat(paste("The null hypothesis that all series have a unit root, is not
                  rejected at a significance level of ", level, ".\n", sep = ""))
      }
      print(out)
    }
  } else { # No Union Tests
    detr_names <- rep(NA, length(inputs$detr))
    detr_names[inputs$detr=="OLS"] <- c("detr = OLS")
    detr_names[inputs$detr=="QD"] <- c("detr = QD")

    dc_names <- rep(NA, length(inputs$dc))
    dc_names[inputs$dc==0] <- c("dc = none")
    dc_names[inputs$dc==1] <- c("dc = intercept")
    dc_names[inputs$dc==2] <- c("dc = intercept and trend")

    detr_dc_names <- paste(rep(detr_names, each = length(inputs$dc)), ", ",
                           rep(dc_names, length(inputs$detr)), sep = "")
    GM_test <- rowMeans(inputs$tests_i)
    t_star <- apply(inputs$t_star, 1:2, mean)
    p_val <- sapply(1:length(detr_dc_names), function(i){mean(t_star[, i] < GM_test[i])})
    out <- cbind("test statistic" = GM_test, "p-value" = p_val)
    rownames(out) <- detr_dc_names

    if (verbose) {
      if (nrow(out) == 1){
        cat("Panel Bootstrap Group-Mean Test\n")
      } else {
        cat("Panel Bootstrap Group-Mean Tests\n")
      }
      for (i in 1:nrow(out)) {
        if (out[i, 2] < level) {
          cat(paste(rownames(out)[i], ": The null hypothesis that all series have a unit root,
                    is rejected at a significance level of ", level, ".\n", sep = ""))
        } else {
          cat(paste(rownames(out)[i], ": The null hypothesis that all series have a unit root,
                    is not rejected at a significance level of ", level, ".\n", sep = ""))
        }
      }
      print(out)
    }
  }
  return(out)
}

#' Differences of Multiple Time Series
#' @description Performs differencing of multiple time series, with possibly different orders for each time series.
#' @param y A (\eqn{T}x\eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}) or data frame.
#' @param d An \eqn{N}-dimensional vector containing the orders
#' @param keep_NAs Logical indicator whether or not to keep the \code{NA} values resulting from differencing at the beginning of the sample. Default is \code{TRUE}. If \code{FALSE}, the entire row containing the \code{NA} values is removed.
#' @export
#' @return The appropriately differenced data in the same format as the original data.
diff_mult <- function(y, d, keep_NAs = TRUE) {
  x <- as.matrix(y)
  diffed_x <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  if (length(c(d)) != ncol(x)) {
    stop("Argument d should have length equal to columns of y.")
  }
  if (any(d != round(d) | d < 0 )) {
    stop("Argument d may only contain integer larger or equal to 0.")
  }
  for (i in 1:ncol(x)) {
    if (d[i] >= 1) {
      diffed_x[(1 + d[i]):nrow(x), i] <- diff(x[, i], differences = d[i])
    } else if (d[i] == 0) {
      diffed_x[, i] <- x[, i]
    }
  }
  diffed_y <- y
  diffed_y[] <- diffed_x
  if (!keep_NAs & (max(d) > 0)) {
    diffed_y <- diffed_y[-(1:max(d)), ]
  }
  return(diffed_y)
}

#' Determine Order of Integration
#' @description Determines the order of integration for each time series in a dataset via a sequence of unit root tests, and differences the data accordingly to eliminate stochastic trends.
#' @param y A (\eqn{T}x\eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}) or data frame.
#' @param max_order The maximum order of integration of the time series. Default is 2.
#' @param test The unit root tests to be used in the procedure. For multiple time series the options are "iADFtest", "BSQTtest" and "bFDRtest", with "iADFtest" the default. For single time series the options are "boot_df" and "boot_union", with the latter the default.
#' @param plot_orders Logical indicator whether the resulting orders of integration should be plotted. Default is \code{FALSE}.
#' @param ... Optional arguments passed to the chosen unit root test function.
#' @details The function follows the approach laid out in Smeekes and Wijler (2020), where all series is differenced \eqn{d-1} times, where \eqn{d} is the specified maximum order, and these differenced series are tested for unit roots. The series foe which the unit root null is not rejected, are classified as \eqn{I(d)} and removed from consideration. The remaining series are integrated, and tested for unit roots again, leading to a classification of \eqn{I(d-1)} series if the null is not rejected. This is continued until a non-rejection is observed for all time series, or the series are integrated back to their original level. The series for which the null hypothesis is rejected in the final stage are classified as \eqn{I(0)}.
#'
#' Care must be taken when using \code{\link{BSQTtest}} when the argument \code{q} is given as a sequence of integers. As at each step series are removed, one may end up with fewer series to test than indicated in \code{q}. While integers larger than the number of series will automatically be removed - along with a warning - by the test, it is recommend to set \code{q} in the form of quantiles.
#'
#' Plotting the orders of integration requires the \code{ggplot2} package to be installed; plot will be skipped and a warning is given if not. For plots the function \code{\link{plot_order_integration}} is called. The user may prefer to set \code{plot_orders = FALSE} and call this function directly using the returned value of \code{order_int} in order to have more control over plot settings and save the plot object.
#' @export
#' @return A list with the following components
#' \item{\code{order_int}}{A vector with the found orders of integration of each time series.}
#' \item{\code{diff_data}}{The appropriately differenced data according to \code{order_int} in the same format as the original data.}
#' @references Smeekes, S. and Wijler, E. (2020). Unit roots and cointegration. In P. Fuleky (Ed.) \emph{Macroeconomic Forecasting in the Era of Big Data}, Chapter 17, pp. 541-584. \emph{Advanced Studies in Theoretical and Applied Econometrics}, vol. 52. Springer.
#' @examples
#' # Use iADFtest to determine the order of GDP_BE and GDP_DE
#' orders_iADFtest <- order_integration(MacroTS[, 1:2], test = "iADFtest", B = 199)
order_integration <- function(y, max_order = 2, test = NULL, plot_orders = FALSE, ...) {
  N <- NCOL(y)
  if (is.null(test)) {
    if (N == 1) {
      test = "boot_union"
    } else {
      test = "iADFtest"
    }
  }
  if (!is.null(colnames(y))) {
    var_names <- colnames(y)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(y))
  }
  d <- rep(NA, N)
  names(d) <- var_names
  yd <- y
  i_in_yd <- 1:N
  for (d_i in (max_order - 1):0) {
    yd <- diff_mult(y[, i_in_yd], rep(d_i, length(i_in_yd)), keep_NAs = FALSE)
    if (test == "iADFtest") {
      out <- iADFtest(yd, ...)
    } else if (test == "bFDRtest" & N > 1) {
      out <- bFDRtest(yd, ...)
    } else if (test == "BSQTtest" & N > 1) {
      out <- BSQTtest(yd, ...)
    } else if (test == "boot_df" & N == 1) {
      out <- boot_df(yd, ...)
    } else if (test == "boot_union" & N == 1) {
      out <- boot_union(yd, ...)
    } else {
      stop("Invalid test argument.")
    }
    d[i_in_yd[!out$rej_H0]] <- d_i + 1
    if (any(out$rej_H0)) {
      if (d_i == 0) {
        d[i_in_yd[out$rej_H0]] <- 0
      } else {
        i_in_yd <- i_in_yd[out$rej_H0]
      }
    } else {
      break
    }
  }
  if (plot_orders) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Cannot plot orders of integration as package ggplot2 not installed.")
    } else {
      g <- plot_order_integration(d)
      print(g)
    }
  }
  return(list(order_int = d, diff_data = diff_mult(y, d)))
}

#' Plot Orders of Integration
#' @description Plots a vector with orders of integration of time series.
#' @param d T\eqn{N}-dimensional vector with time series' orders of integration. Elements should be named after the respective time series to ensure easy interpretation of the plot.
#' @param show_names Show the time series' names on the plot (\code{TRUE}) or not (\code{FALSE}). Default is \code{TRUE}.
#' @param show_legend Logical indicator whether a legend should be displayed. Default is \code{TRUE}.
#' @param names_size Size of the time series' names if \code{show_names = TRUE}. Default takes \code{ggplot2} defaults.
#' @param legend_size Size of the text in the legend if \code{show_legend = TRUE}. Default takes \code{ggplot2} defaults.
#' @param cols Vector with colours for displaying the orders of integration. At least as many colours as orders of integration need to be supplied. Default supplies 4 colours for displaying up to \eqn{I(3)} series.
#' @export
#' @details This function requires the package \code{ggplot2} to be installed. If the package is not found, plotting is aborted.
#' @return A \code{ggplot2} object containing the plot of the orders of integration.
plot_order_integration <- function(d, show_names = TRUE, show_legend = TRUE,
                                   names_size = NULL, legend_size = NULL, cols = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Cannot plot orders of integration as package ggplot2 not installed.")
  } else {
    if (is.null(cols)) {
      cols <- c("#1B9E77", "#7570B3", "#D95F02", "#E7298A")
    }
    if (max(d) > length(cols)) {
      stop("Insufficient number of colours supplied to display all orders of integration.")
    }
    cols <- cols[unique(d) + 1]
    n_g <- floor(max(1, length(d) - 5)^(1/3))
    group <- rep(paste0("g", 1:n_g), each = ceiling(length(d) / n_g))[1:length(d)]
    df <- data.frame(names = factor(names(d), levels = rev(names(d))),
                     order = paste0("I(", d, ")"), group = group)
    if (show_names) {
      g <- ggplot2::ggplot(df) +
        ggplot2::geom_col(ggplot2::aes(x = names, y = order, fill = order),
                          show.legend = show_legend) +
        ggplot2::labs(y = "Order of Integration", x = "Variables") +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::theme_minimal() +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.title = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(size = legend_size),
                       axis.text.x = ggplot2::element_text(size = names_size))
    } else {
      g <- ggplot2::ggplot(df) +
        ggplot2::geom_col(ggplot2::aes(x = names, y = order, fill = order),
                          show.legend = show_legend) +
        ggplot2::labs(y = "Order of Integration", x = "Variables") +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::theme_minimal() +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.title = ggplot2::element_blank(),
                       legend.text = ggplot2::element_text(size = legend_size),
                       axis.text.x = ggplot2::element_blank())
    }

    if (n_g > 1) {
      g <- g + ggplot2::facet_wrap(ggplot2::vars(df$group), nrow = 1, scales = "free_y") +
        ggplot2::theme(strip.text = ggplot2::element_blank())
    }
    return(g)
  }
}

#' Check Missing Values in Sample
#' @param X A (\eqn{T}x\eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}) or data frame.
#' @export
#' @return \eqn{N}-dimensional vector, for each series whether missing values are present (\code{TRUE}) or not (\code{FALSE})
check_missing_insample_values <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  check <- apply(X, 2, function(x){any(diff((1:n)[!is.na(x)]) != 1)})
  return(check)
}

#' Find Non-Missing Subsamples
#' @param X A (\eqn{T}x\eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}) or data frame. Assumes a prior check on missing values in-sample has been done.
#' @export
#' @return A list with the following components
#' \item{\code{range}}{(2x\eqn{N})-dimensional matrix containing the first and last non-missing observation in each column of X.}
#' \item{\code{all_equal}}{Logical value indicating whether all series have the same non-missing indices.}
find_nonmissing_subsample <- function(X) {
  names <- colnames(X)
  X <- as.matrix(X)
  n <- NROW(X)
  N <- NCOL(X)
  ind <- matrix(1:n, nrow = n, ncol = N)
  ind[is.na(X)] <- NA
  range_nonmiss <- matrix(apply(ind, 2, function(x){c(min(stats::na.omit(x)),
                                                      max(stats::na.omit(x)))}), ncol = N)
  colnames(range_nonmiss) <- names
  rownames(range_nonmiss) <- c("first", "last")
  all_equal <- (max(range_nonmiss[1, ]) - min(range_nonmiss[1, ]) == 0) &
    (max(range_nonmiss[2, ]) - min(range_nonmiss[2, ]) == 0)
  return(list(range = range_nonmiss, all_equal = all_equal))
}

#' Plot Missing Values
#' @description Plots missing values of different types for a time series dataset.
#' @param y A (\eqn{T}x\eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}) or data frame.
#' @param show_names Show the time series' names on the plot (\code{TRUE}) or not (\code{FALSE}). Default is \code{TRUE}.
#' @param show_legend Logical indicator whether a legend should be displayed. Default is \code{TRUE}.
#' @param axis_text_size Size of the text on the axis. Default takes \code{ggplot2} defaults.
#' @param legend_size Size of the text in the legend if \code{show_legend = TRUE}. Default takes \code{ggplot2} defaults.
#' @param cols Vector with colours for displaying the different types of data. If the default is overwriten, four colours must be supplied.
#' @export
#' @details The function distinguish four types of data: observed data (non-missing) and three missing types. Type \code{"Balanced NA"} indicates where entire rows are missing (\code{NA}). These do not cause unbalancedness as the missing rows can simply be deleted.  Type \code{"Unalanced NA"} are missing values on the beginning or end of the sample, which cause unbalancedness. These affect some (but not all) bootstrap methods, see e.g.~\code{\link{bFDRtest}}. Type \code{"Internal NA"} are missing values inside the sample, which need to be removed before the bootstrap unit root tests can be used.
#'
#' This function requires the package \code{ggplot2} to be installed. If the package is not found, plotting is aborted.
#' @return A \code{ggplot2} object containing the missing values plot.
plot_missing_values <- function(y, show_names = FALSE, show_legend = TRUE,
                                axis_text_size = NULL, legend_size = NULL,
                                cols = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Cannot plot missing values as package ggplot2 not installed.")
  } else {
    y <- as.matrix(y)
    if (!is.null(colnames(y))) {
      var_names <- colnames(y)
    } else {
      var_names <- paste0("Variable ", 1:NCOL(y))
    }
    check_missing <- check_missing_insample_values(y)
    fns <- find_nonmissing_subsample(y)
    range <- fns$range
    all_equal <- fns$all_equal
    missing_type <- c("Observed", "Balanced NA", "Unbalanced NA", "Internal NA")
    if (is.null(cols)){
      cols <- c("#4DAF4A", "#377EB8", "#984EA3", "#E41A1C")
    }
    if (length(cols) < 4) {
      stop("Insufficient colors supplied.")
    }
    names(cols) <- missing_type
    x <- array(missing_type[1], dim = dim(y))
    min_index <- min(range[1, ])
    max_index <- max(range[2, ])
    if (min_index > 1) {
      x[1:(min_index - 1), ] <- missing_type[2]
    }
    if (max_index < nrow(x)) {
      x[(max_index + 1):nrow(x), ] <- missing_type[2]
    }
    if (!all_equal) {
      unb_miss <- is.na(y) * ((1:nrow(x)) >= min_index) * ((1:nrow(x)) <= max_index)
      x[unb_miss == 1] <- missing_type[3]
    }
    if (any(check_missing)) {
      for (i in (1:ncol(x))[check_missing]) {
        int_miss_i <- is.na(y[, i]) * ((1:nrow(x)) > range[1, i]) * ((1:nrow(x)) < range[2, i])
        x[int_miss_i == 1, i] <- missing_type[4]
      }
    }
    col_sel <- cols[sapply(1:4, function(i){any(x == missing_type[i])})]
    obs <- 1:nrow(x)
    df <- data.frame(var_names = factor(rep(var_names, each = nrow(x)), levels = var_names),
                     obs = factor(rep(obs, ncol(x)), levels = rev(obs)),
                     missing_type = factor(c(x), levels = missing_type))

    if (show_names) {
      g <- ggplot2::ggplot(data = df, ggplot2::aes(x = var_names, y = obs,
                                                   fill = missing_type)) +
        ggplot2::geom_tile(height = 0.8, width = 0.8, show.legend = show_legend) +
        ggplot2::labs(x = "Variables", y = "Observations", fill = "Missing Type") +
        ggplot2::scale_fill_manual(values = col_sel) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = axis_text_size),
                       axis.text.y = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(size = legend_size),
                       legend.text = ggplot2::element_text(size = legend_size),
                       panel.border = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(size = axis_text_size))
    } else {
      g <- ggplot2::ggplot(data = df, ggplot2::aes(x = var_names, y = obs,
                                                   fill = missing_type)) +
        ggplot2::geom_tile(height = 0.8, width = 0.8, show.legend = show_legend) +
        ggplot2::labs(x = "Variables", y = "Observations", fill = "Missing Type") +
        ggplot2::scale_fill_manual(values = col_sel) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.title = ggplot2::element_text(size = legend_size),
                       legend.text = ggplot2::element_text(size = legend_size),
                       axis.text.x = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(size = axis_text_size))
    }
    return(g)
  }
}
