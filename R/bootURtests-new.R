#' Bootstrap augmented Dickey-Fuller Unit Root Tests without multiple testing control
#' @description This function performs a standard augmented Dickey-Fuller bootstrap unit root test on a single time series or on each time series individually in case multiple time series are provided.
#' @param data A \eqn{T}-dimensional vector or a (\eqn{T} x \eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations to be tested for unit roots. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}), or a data frame, as long as each column represents a single time series.
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
#' @param p_min Minimum lag length in the augmented Dickey-Fuller regression. Default is 0.
#' @param p_max Maximum lag length in the augmented Dickey-Fuller regression. Default uses the sample size-based rule \eqn{12(T/100)^{1/4}}.
#' @param ic String for information criterion used to select the lag length in the augmented Dickey-Fuller regression. Options are: \code{"AIC"}, \code{"BIC"}, \code{"MAIC"}, \code{"MBIC"}. Default is \code{"MAIC"} (Ng and Perron, 2001).
#' @param dc Numeric vector indicating the deterministic specification. Options are (combinations of)
#'
#' \verb{0 } no deterministics;
#'
#' \verb{1 } intercept only;
#'
#' \verb{2 } intercept and trend.
#'
#' The default is adding an intercept.
#' @param detr String vector indicating the type of detrending to be performed. Options are: \code{"OLS"} and/or \code{"QD"} (typically also called GLS, see Elliott, Rothenberg and Stock, 1996). The default is \code{"OLS"}.
#' @param ic_scale Logical indicator whether or not to use the rescaled information criteria of Cavaliere et al. (2015) (\code{TRUE}) or not (\code{FALSE}). Default is \code{TRUE}.
#' @param verbose Logical indicator whether or not information on the outcome of the unit root test needs to be printed to the console. Default is \code{FALSE}.
#' @param show_progress Logical indicator whether a bootstrap progress update should be printed to the console. Default is FALSE.
#' @param do_parallel Logical indicator whether bootstrap loop should be executed in parallel. Parallel computing is only available if OpenMP can be used, if not this option is ignored. Default is FALSE.
#' @param nc The number of cores to be used in the parallel loops. Default is to use all but one.
#' @details The options encompass many test proposed in the literature. \code{dc = "OLS"} gives the standard augmented Dickey-Fuller test, while \code{dc = "QD"} provides the DF-GLS test of Elliott, Rothenberg and Stock (1996). The bootstrap algorithm is always based on a residual bootstrap (under the alternative) to obtain residuals rather than a difference-based bootstrap (under the null), see e.g. Palm, Smeekes and Urbain (2008).
#'
#' Lag length selection is done automatically in the ADF regression with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. For very short time series (fewer than 50 time points) the maximum lag length is adjusted downward to avoid potential multicollinearity issues in the bootstrap. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `p_min` and maximum lag length `p_max` for the selection algorithm equal to the desired lag length.
#'
#' @return A list with the following components
#' \item{\code{rej_H0}}{Logical indicator whether the null hypothesis of a unit root is rejected (\code{TRUE}) or not (\code{FALSE});}
#' \item{\code{ADF_tests}}{Details on the unit root tests: value of the test statistics and p-values.}
#' The output is arranged per time series, type of deterministic component (\code{dc}) and detrending method (\code{detr}).
#' @section Warnings:
#' The function may give the following warnings.
#' \describe{
#' \item{\code{Warning: Missing values cause resampling bootstrap to be executed for each time series individually.}}{If the time series in \code{data} have different starting and end points (and thus some series contain \code{NA} values at the beginning and/or end of the sample, the resampling-based moving block bootstrap (MBB) and sieve bootstrap (SB) cannot be used directly, as they create holes (internal missings) in the bootstrap samples. These bootstrap methods are therefore not applied jointly as usual, but individually to each series.}
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
#' # boot_adf on GDP_BE and GDP_DE
#' two_series_boot_adf <- boot_adf(MacroTS[, 1:2], boot = "MBB", B = 399,
#' verbose = TRUE)
#' @export
boot_adf <- function(data, level = 0.05, boot = "AWB", B = 1999, l = NULL,
                     ar_AWB = NULL, p_min = 0, p_max = NULL,
                     ic = "MAIC", dc = 1, detr = "OLS", ic_scale = TRUE,
                     verbose = FALSE, show_progress = FALSE,
                     do_parallel = FALSE, nc = NULL){
  
  inputs <- do_tests_and_bootstrap(y = data, BSQT_test = FALSE, iADF_test = TRUE, level = level,
                                   boot = boot, B = B, l = l, ar_AWB = ar_AWB, union = FALSE,
                                   p_min = p_min, p_max = p_max, ic = ic, dc = dc, detr = detr,
                                   ic_scale = ic_scale, q = NULL, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel, nc = nc)
  
  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(data))
  }
  
    # No Union Tests
    detr_names <- rep(NA, length(inputs$detr))
    detr_names[inputs$detr=="OLS"] <- c("detr = OLS")
    detr_names[inputs$detr=="QD"] <- c("detr = QD")
    
    dc_names <- rep(NA, length(inputs$dc))
    dc_names[inputs$dc==0] <- c("dc = none")
    dc_names[inputs$dc==1] <- c("dc = intercept")
    dc_names[inputs$dc==2] <- c("dc = intercept and trend")
    
    detr_dc_names <- paste(rep(detr_names, each = length(inputs$dc)), ", ",
                           rep(dc_names, length(inputs$detr)), sep = "")
    iADFout <- array(NA, c(NCOL(data), 2, length(detr_dc_names)), dimnames =
                       list(var_names, c("test statistic", "p-value"), detr_dc_names))
    rej_H0 <- array(NA, c(NCOL(data), length(inputs$dc)*length(inputs$detr)), dimnames =
                      list(var_names, detr_dc_names))
    for(i in 1:nrow(inputs$tests_i)){
      iADFout[, , i] <- iADF_cpp(test_i = matrix(inputs$tests_i[i, ], nrow = 1),
                                 t_star = matrix(inputs$t_star[ , i, ], nrow = B),
                                 level = inputs$level)
      rej_H0[, i] <- (iADFout[, 2, i] < level)
      
      if (verbose) {
        cat(paste(c(rep("-", 40), "\n"), sep = "", collapse = ""))
        cat(paste("Type of unit root test performed: ", detr_dc_names[i], "\n", sep = ""))
        if (NCOL(data) > 1) {
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

  return(list(rej_H0 = rej_H0, ADF_tests = iADFout))
}


#' Bootstrap Union Test for Unit Roots
#' @description Performs bootstrap unit root tests based on the union of rejections of 4 tests with different number of deterministic components and different type of detrending (Harvey, Leybourne and Taylor, 2012; Smeekes and Taylor, 2012).
#' @inheritParams boot_adf
#' @param boot String for bootstrap method to be used. Options are
#' \describe{
#' \item{\code{"MBB"}}{Moving blocks bootstrap (Paparoditis and Politis, 2003; Palm, Smeekes and Urbain, 2011);}
#' \item{\code{"BWB"}}{Block wild bootstrap (Shao, 2011; Smeekes and Urbain, 2014a);}
#' \item{\code{"DWB"}}{Dependent wild bootstrap (Shao, 2010; Smeekes and Urbain, 2014a; Rho and Shao, 2019);}
#' \item{\code{"AWB"}}{Autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2020), this is the default;}
#' \item{\code{"SB"}}{Sieve bootstrap (Chang and Park, 2003; Palm, Smeekes and Urbain, 2008; Smeekes, 2013);}
#' \item{\code{"SWB"}}{Sieve wild boostrap (Cavaliere and Taylor, 2009; Smeekes and Taylor, 2012).}
#' }
#' @details The union is taken over the combination of tests with intercept only and intercept plus trend, coupled with OLS detrending and QD detrending, as in Harvey, Leybourne and Taylor (2012) and Smeekes an Taylor (2012). The bootstrap algorithm is always based on a residual bootstrap (under the alternative) to obtain residuals rather than a difference-based bootstrap (under the null), see e.g. Palm, Smeekes and Urbain (2008).
#'
#' Lag length selection is done automatically in the ADF regressions with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `p_min` and maximum lag length `p_max` for the selection algorithm equal to the desired lag length.
#' @export
#' @return Value of the union test statistic and the bootstrap p-values.
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
#' @examples
#' # boot_union on GDP_BE
#' GDP_BE_df <- boot_union(MacroTS[, 1], B = 399, verbose = TRUE)
boot_union <- function(data, level = 0.05, boot = "AWB", B = 1999, l = NULL, ar_AWB = NULL,
                       p_min = 0, p_max = NULL, ic = "MAIC", ic_scale = TRUE, verbose = FALSE,
                       show_progress = FALSE, do_parallel = FALSE, nc = NULL){
  
  if (verbose) {
    cat("Bootstrap Test with", boot, "bootstrap method.\n")
  }
  
  inputs <- do_tests_and_bootstrap(y = data, BSQT_test = FALSE, iADF_test = TRUE, level = level,
                                   boot = boot, B = B, l = l, ar_AWB = ar_AWB, union = TRUE,
                                   p_min = p_min, p_max = p_max, ic = ic, dc = NULL, detr = NULL,
                                   ic_scale = ic_scale, q = NULL, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel, nc = nc)
  
  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(data))
  }
  
  # Union Tests
    iADFout <- iADF_cpp(test_i = inputs$test_stats, t_star = inputs$test_stats_star,
                        level = inputs$level)
    rej_H0 <- (iADFout[, 2] < level)
    colnames(iADFout) <- c("test statistic", "p-value")
    rownames(iADFout) <- var_names
    
    if (verbose) {
      p_hat <- sum(rej_H0)
      if (NCOL(data) > 1) {
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
  
  return(list(rej_H0 = rej_H0, ADF_tests = iADFout))
}

#' Bootstrap Unit Root Tests with False Discovery Rate control
#' @description Controls for multiple testing by controlling the false discovery rate (FDR), see Moon and Perron (2012) and Romano, Shaikh and Wolf (2008).
#' @inheritParams boot_adf
#' @param level Desired False Discovery Rate level of the unit root tests. Default is 0.05.
#' @param union Logical indicator whether or not to use bootstrap union tests (\code{TRUE}) or not (\code{FALSE}), see Smeekes and Taylor (2012). Default is \code{TRUE}.
#' @details The false discovery rate FDR is defined as the expected proportion of false rejections relative to the total number of rejections.
#'
#' See \code{\link{iADFtest}} for details on the bootstrap algorithm and lag selection.
#' @return A list with the following components
#' \item{\code{rej_H0}}{Logical indicator whether the null hypothesis of a unit root is rejected (\code{TRUE}) or not (\code{FALSE});}
#' \item{\code{FDR_sequence}}{Details on the unit root tests: value of the test statistics and critical values.}
#' For the union test (\code{union = TRUE}), the output is arranged per time series. If \code{union = FALSE}, the output is arranged per time series, type of deterministic component (\code{dc}) and detrending method (\code{detr}).
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Resampling-based bootstraps MBB and SB cannot handle missing values.}}{If the time series in \code{data} have different starting and end points (and thus some series contain \code{NA} values at the beginning and/or end of the sample, the resampling-based moving block bootstrap (MBB) and sieve bootstrap (SB) cannot be used, as they create holes (internal missings) in the bootstrap samples. Switch to another bootstrap method or truncate your sample to eliminate \code{NA} values.}
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
#' # boot_fdr on GDP_BE and GDP_DE
#' two_series_boot_fdr <- boot_fdr(MacroTS[, 1:2], boot = "MBB", B = 399,  verbose = TRUE)
#' @export
boot_fdr <- function(data, level = 0.05,  boot = "AWB", B = 1999, l = NULL, ar_AWB = NULL,
                     union = TRUE, p_min = 0, p_max = NULL, ic = "MAIC", dc = NULL,
                     detr = NULL, ic_scale = TRUE, verbose = FALSE, show_progress = FALSE,
                     do_parallel = FALSE, nc = NULL){
  
  inputs <- do_tests_and_bootstrap(y = data, BSQT_test = FALSE, iADF_test = FALSE, level = level,
                                   boot = boot, B = B, l = l, ar_AWB = ar_AWB, union = union,
                                   p_min = p_min, p_max = p_max, ic = ic, dc = dc, detr = detr,
                                   ic_scale = ic_scale, q = NULL, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel, nc = nc)
  
  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(data))
  }
  
  if (union) { # Union Tests
    bFDRout <- FDR_cpp(test_i = inputs$test_stats, t_star = inputs$test_stats_star,
                       level = inputs$level)
    rej_H0 <- matrix(bFDRout$rej_H0 == 1, nrow = NCOL(data))
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
    rej_H0 <- matrix(nrow = NCOL(data), ncol = length(inputs$dc)*length(inputs$detr))
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
#' @inheritParams boot_adf
#' @param q Numeric vector of quantiles or units to be tested. Default is to test each unit sequentially.
#' @param union Logical indicator whether or not to use bootstrap union tests (\code{TRUE}) or not (\code{FALSE}), see Smeekes and Taylor (2012). Default is \code{TRUE}.
#' @details The parameter \code{q} can either be set as an increasing sequence of integers smaller or equal to the number of series \code{N}, or fractions of the total number of series (quantiles). For \code{N} time series, setting \code{q = 0:N} means each unit should be tested sequentially. In this case the method is equivalent to the StepM method of Romano and Wolf (2005), and therefore controls the familywise error rate. To split the series in \code{K} equally sized groups, use \code{q = 0:K / K}.
#'
#' By convention and in accordance with notation in Smeekes (2015), the first entry of the vector should be equal to zero, while the second entry indicates the end of the first group, and so on. If the initial \code{0} or final value (\code{1} or \code{N}) are omitted, they are automatically added by the function.
#'
#' See \code{\link{iADFtest}} for details on the bootstrap algorithm and lag selection.
#' @return A list with the following components
#' \item{\code{rej_H0}}{Logical indicator whether the null hypothesis of a unit root is rejected (\code{TRUE}) or not (\code{FALSE});}
#' \item{\code{BSQT_sequence}}{Details on the unit root tests: outcome of the sequential steps, value of the test statistics and p-values.}
#' For the union test (\code{union = TRUE}), the output is arranged per time series. If \code{union = FALSE}, the output is arranged per time series, type of deterministic component (\code{dc}) and detrending method (\code{detr}).
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Resampling-based bootstraps MBB and SB cannot handle missing values.}}{If the time series in \code{data} have different starting and end points (and thus some series contain \code{NA} values at the beginning and/or end of the sample, the resampling-based moving block bootstrap (MBB) and sieve bootstrap (SB) cannot be used, as they create holes (internal missings) in the bootstrap samples. Switch to another bootstrap method or truncate your sample to eliminate \code{NA} values.}
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
#' # boot_sqt on GDP_BE and GDP_DE
#' two_series_boot_sqt <- boot_sqt(MacroTS[, 1:2], boot = "AWB", B = 399,  verbose = TRUE)
#' @export
boot_sqt <- function(data, q = 0:NCOL(data), level = 0.05,  boot = "AWB", B = 1999, l = NULL,
                     ar_AWB = NULL, union = TRUE, p_min = 0, p_max = NULL, ic = "MAIC", dc = NULL,
                     detr = NULL, ic_scale = TRUE, verbose = FALSE, show_progress = FALSE,
                     do_parallel = FALSE, nc = NULL){
  
  inputs <- do_tests_and_bootstrap(y = data, BSQT_test = TRUE, iADF_test = FALSE, level = level,
                                   boot = boot, B = B, l = l, ar_AWB = ar_AWB, union = union,
                                   p_min = p_min, p_max = p_max, ic = ic, dc = dc, detr = detr,
                                   ic_scale = ic_scale, q = q, h_rs = 0.1,
                                   show_progress = show_progress, do_parallel = do_parallel, nc = nc)
  
  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(data))
  }
  
  if (union) { # Union Tests
    BSQTout <- BSQT_cpp(pvec = inputs$p_vec, test_i = inputs$test_stats, t_star =
                          inputs$test_stats_star, level = inputs$level)
    rej_H0 <- matrix(BSQTout$rej_H0 == 1, nrow = NCOL(data))
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
    rej_H0 <- matrix(nrow = NCOL(data), ncol = length(inputs$dc)*length(inputs$detr))
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
#' @inheritParams boot_adf
#' @param union Logical indicator whether or not to use bootstrap union tests (\code{TRUE}) or not (\code{FALSE}), see Smeekes and Taylor (2012). Default is \code{TRUE}.
#' @details See \code{\link{iADFtest}} for details on the bootstrap algorithm and lag selection.
#'
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Resampling-based bootstraps MBB and SB cannot handle missing values.}}{If the time series in \code{data} have different starting and end points (and thus some series contain \code{NA} values at the beginning and/or end of the sample, the resampling-based moving block bootstrap (MBB) and sieve bootstrap (SB) cannot be used, as they create holes (internal missings) in the bootstrap samples. Switch to another bootstrap method or truncate your sample to eliminate \code{NA} values.}
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
#' # boot_panel on GDP_BE and GDP_DE
#' two_series_boot_panel <- boot_panel(MacroTS[, 1:2], boot = "AWB", B = 399,  verbose = TRUE)
#' @export
boot_panel <- function(data, level = 0.05,  boot = "AWB", B = 1999, l = NULL, ar_AWB = NULL,
                      union = TRUE, p_min = 0, p_max = NULL, ic = "MAIC", dc = NULL, detr = NULL,
                      ic_scale = TRUE, verbose = FALSE, show_progress = FALSE,
                      do_parallel = FALSE, nc = NULL){
  
  inputs <- do_tests_and_bootstrap(y = data, BSQT_test = FALSE, iADF_test = FALSE, level = level,
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