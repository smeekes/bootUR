#' Individual Unit Root Tests without multiple testing control
#' @description This function performs bootstrap unit root tests on each time series individually.
#' @param y A \eqn{T}-dimensional vector or a (\eqn{T} x \eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations to be tested for unit roots. Data may also be in a time series format (e.g. \verb{ts}, \verb{zoo} or \verb{xts}).
#' @param level Desired significance level of the unit root test. Default is 0.05.
#' @param boot String for bootstrap method to be used. Options are
#' \describe{
#' \item{\verb{"MBB"}}{Moving blocks bootstrap (Paparoditis and Politis, 2003; Palm, Smeekes and Urbain, 2011), this is the default;}
#' \item{\verb{"BWB"}}{Block wild bootstrap (Shao, 2011);}
#' \item{\verb{"DWB"}}{Dependent wild bootstrap (Shao, 2010; Smeekes and Urbain, 2014a; Rho and Shao, 2019);}
#' \item{\verb{"AWB"}}{Autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2020);}
#' \item{\verb{"SB"}}{Sieve bootstrap (Palm, Smeekes and Urbain,2008; Smeekes and Urbain, 2014a);}
#' \item{\verb{"SWB"}}{Sieve wild boostrap (Cavaliere and Taylor, 2008; Smeekes and Taylor, 2014b).}
#' }
#' @param B Number of bootstrap replications. Default is 9999.
#' @param l Desired 'block length' in the bootstrap. For the MBB, BWB and DWB boostrap, this is a genuine block length. For the AWB boostrap, the block length is transformed into an autoregressive parameter via the formula \eqn{0.01^(1/l)} as in Smeekes and Urbain (2014); this can be overwritten by setting \verb{ar_AWB} directly. Default sets the block length as a function of the time series length T, via the rule \eqn{l = 1.75 T^(1/3)}.
#' @param ar_AWB Autoregressive parameter used in the AWB bootstrap method (\verb{boot = "AWB"}). Can be used to set the parameter directly rather than via the default link to the block length l.
#' @param union Logical indicator whether or not to use bootstrap union tests (\verb{TRUE}) or not (\verb{FALSE}), see Harvey, Leybourne and Taylor (2012) and Smeekes and Taylor (2012). Default is \verb{TRUE}.
#' @param p.min Minimum lag length in the augmented Dickey-Fuller regression. Default is 0.
#' @param p.max Maximum lag length in the augmented Dickey-Fuller regression. Default uses the sample size-based rule \eqn{12(T/100)^{1/4}}.
#' @param ic String for information criterion used to select the lag length in the augmented Dickey-Fuller regression. Options are: \verb{"AIC"}, \verb{"BIC"}, \verb{"MAIC"}, \verb{"MBIC}. Default is \verb{"MAIC"} (Ng and Perron, 2001).
#' @param dc Numeric vector indicating the deterministic specification. Only relevant if \verb{union = FALSE}. Options are: 0: no deterministics, 1: intercept only, 2: intercept and trend. 
#' \describe{
#' \item{\emph{0}}{no deterministics;}
#' \item{\emph{1}}{intercept only;}
#' \item{\emph{2}}{intercept and trend.}
#' Combinations thereof are allowed. Default is the union test (\verb{union = TRUE}), in which case this is not relevant. If \verb{union = FALSE}, the default is adding an intercept (a warning is given).
#' }
#' @param detr String vector indicating the type of detrending to be performed. Only relevant if \verb{unionunion = FALSE}. Options are \verb{"OLS"} and/or \verb{"QD"} (typically also called GLS, see Elliott, Rothenberg and Stock, 1996). Default is \verb{"OLS"}.
#' @param ic.scale Logical indicator whether or not to use the rescaled
#' information criteria of Cavaliere et al. (2015) (\verb{TRUE}) or not (\verb{FALSE}).
#' Default is \verb{TRUE}.
#' @param verbose Logical indicator whether or not information on the outcome of
#'  the unit root test needs to be printed to the console. Default is \verb{FALSE}.
#' @details Lag selection (rescaling and Perron-Qu correction). Residual bootstrap.
#' @export
#' @return A list with the following components:
#' \describe{
#' \item{\verb{rej_H0}}{Logical indicator whether the null hypothesis of a unit root is rejected (\verb{TRUE}) or not (\verb{FALSE});}
#' \item{\verb{ADF_tests}}{Details on the unit root tests: value of the test statistics and p-values.}
#' For the union test (\verb{union = TRUE}), the output is arranged per time series. If \verb{union = FALSE}, the output is arranged per time series, type of deterministic component (verb{dc}) and detrending method (verb{detr}).
#' }
#' @references Chang, Y. and Park, J. (2003). A sieve bootstrap for the test of a unit root. \emph{Journal of Time Series Analysis}, 24(4), 379-400.
#' @references Cavaliere, G. and Taylor, A.M.R (2009). Heteroskedastic time
#' series with a unit root. \emph{Econometric Theory}, 25, 1228–1276.
#' @references Cavaliere, G., Phillips, P.C.B., Smeekes, S., and Taylor, A.M.R.
#' (2015). Lag length selection for unit root tests in the presence of nonstationary volatility. \emph{Econometric Reviews}, 34(4), 512-536.
#' @references Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996). Efficient tests for an autoregressive unit root. \emph{Econometrica}, 64(4), 813-836.
#' @references Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild bootstrap inference for nonparametric trends. \emph{Journal of Econometrics}, 214(1), 81-109.
#' @references Harvey, D.I., Leybourne, S.J., and Taylor, A.M.R. (2012). Testing for unit roots in the presence of uncertainty over both the trend and initial condition. \emph{Journal of Econometrics}, 169(2), 188-195.
#' @references Moon, H.R. and Perron, B. (2012). Beyond panel unit root tests: Using multiple testing to determine the non stationarity properties of individual series in a panel. Journal of Econometrics, 169(1), 29-33.
#' @references Ng, S. and Perron, P. (2001). Lag Length Selection and the Construction of Unit Root Tests with Good Size and Power. \emph{Econometrica}, 69(6), 1519-1554,
#' @references Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008). Bootstrap unit root tests: Comparison and extensions. \emph{Journal of Time Series Analysis}, 29(1), 371-401.
#' @references Perron, P. and Qu, Z. (2008). A simple modification to improve the finite sample properties of Ng and Perron's unit root tests. \emph{Economic Letters}, 94(1), 12-19.
#' @references Rho, Y. and Shao, X. (2019). Bootstrap-assisted unit root testing with piecewise locally stationary errors. \emph{Econometric Theory}, 35(1), 142-166.
#' @references Shao, X. (2010). The dependent wild bootstrap. \emph{Journal of the American Statistical Association}, 105(489), 218-235.
#' @references Shao, X. (2011). A bootstrap-assisted spectral test of white noise under unknown dependence. \emph{Journal of Econometrics}, 162, 213-224.
#' @references Smeekes, S. (2015). Bootstrap sequential tests to determine the order of integration of individual units in a time series panel. \emph{Journal of Time Series Analysis}, 36(3), 398-415.
#' @references Smeekes, S. and Taylor, A.M.R. (2012). Bootstrap union tests for unit roots in the presence of nonstationary volatility. \emph{Econometric Theory}, 28(2), 422-456.
#' @references Smeekes, S. and Urbain, J.-P. (2014a). A multivariate invariance principle for modified wild bootstrap methods with an application to unit root testing. GSBE Research Memorandum No. RM/14/008, Maastricht University
#' @references Smeekes, S. and Urbain, J.-P. (2014b). On the applicability of the sieve bootstrap in time series panels. \emph{Oxford Bulletin of Economics and Statistics}, 76(1), 139-151.
#' @examples
#' # iADFtest on GDP_BE and GDP_DE
#' two_series_iADFtest <- iADFtest(eurostat[, 1:2], boot = "MBB", B=399, verbose = TRUE)
iADFtest <- function(y, level = 0.05, boot = "MBB", B = 9999, l = NULL,
                     ar_AWB = NULL, union = TRUE, p.min = 0, p.max = NULL,
                     ic = "MAIC", dc = NULL, detr = NULL, ic.scale = TRUE,
                     verbose = FALSE){

  inputs <- generate_inputs(y = y, BSQT_test = FALSE, iADF_test = TRUE,
                            level = level, boot = boot, B = B, union = union,
                            p.min = p.min, p.max = p.max, ic = ic, dc = dc,
                            detr = detr, q = NULL, l = l, ic.scale = ic.scale,
                            h.rs = 0.1, k_DWB = "k.TBB", ar_AWB = ar_AWB)

  if (!is.null(colnames(y))) {
    var_names <- colnames(y)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(y))
  }

  if (union) { # Union Tests
    iADFout <- iADF_cpp(test_i = inputs$test_stats, t_star = inputs$test_stats_star, level = inputs$level)
    rej_H0 <- (iADFout[, 2] < level)
    colnames(iADFout) <- c("test statistic", "p-value")
    rownames(iADFout) <- var_names

    if (verbose) {
      p_hat <- sum(rej_H0)
      if (NCOL(y) > 1) {
        if (p_hat > 1) {
          cat(paste("There are ", p_hat, " stationary time series, namely: ", paste(var_names[rej_H0], collapse = " "), "\n", sep = ""))
        } else {
          if (p_hat == 1) {
            cat(paste("There is ", p_hat, " stationary time series, namely: ", paste(var_names[rej_H0], collapse = " "), "\n", sep = ""))
          } else {
            cat(paste("There are ", p_hat, " stationary time series.\n", sep = ""))
          }
        }
        print(iADFout)
      } else {
        cat("Bootstrap Union Test:\n")
        if (rej_H0) {
          cat(paste("The null hypothesis of a unit root is rejected at a significance level of ", level, ".\n", sep = ""))
        } else {
          cat(paste("The null hypothesis of a unit root is not rejected at a significance level of ", level, ".\n", sep = ""))
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

    detr_dc_names <- paste(rep(detr_names, each = length(inputs$dc)), ", ", rep(dc_names, length(inputs$detr)), sep = "")
    iADFout <- array(NA, c(NCOL(y), 2, length(detr_dc_names)), dimnames = list(var_names, c("test statistic", "p-value"), detr_dc_names))
    rej_H0 <- array(NA, c(NCOL(y), length(inputs$dc)*length(inputs$detr)), dimnames = list(var_names, detr_dc_names))
    for(i in 1:nrow(inputs$tests_i)){
      iADFout[, , i] <- iADF_cpp(test_i = matrix(inputs$tests_i[i, ], nrow = 1), t_star = matrix(inputs$t_star[ , i, ], nrow = B), level = inputs$level)
      rej_H0[, i] <- (iADFout[, 2, i] < level)

      if (verbose) {
        cat(paste(c(rep("-", 40), "\n"), sep = "", collapse = ""))
        cat(paste("Type of unit root test performed: ", detr_dc_names[i], "\n", sep = ""))
        if (NCOL(y) > 1) {
          p_hat <- sum(rej_H0[, i])
          if (p_hat > 1) {
            cat(paste("There are ", p_hat, " stationary time series, namely: ", paste(var_names[rej_H0[, i]], collapse = " "), ".\n", sep = ""))
          } else {
            if (p_hat == 1) {
              cat(paste("There is ", p_hat, " stationary time series, namely: ", paste(var_names[rej_H0[, i]], collapse = " "), ".\n", sep = ""))
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

#' Bootstrap Dickey-Fuller Unit Root Test
#' @param y Vector with time series data to be tested for unit roots. Data may also be in a time series format (e.g. ts, zoo or xts).
#' @param level Desired significance level of the unit root test. Default is 0.05.
#' @param boot String for bootstrap method to be used. Options are MBB - moving blocks bootstrap (Moon and Perron, 2012; Smeekes, 2015), BWB - block wild bootstrap (Shao, 2011), DWB - dependent wild bootstrap (Shao, 2010; Smeekes and Urbain, 2014a; Rho and Shao, 2019), AWB - autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2018), SB - sieve bootstrap (Palm, Smeekes and Urbain, 2008; Smeekes and Urbain, 2014a), SWB - sieve wild boostrap (Smeekes and Taylor, 2014b). Default is MBB.
#' @param B Number of bootstrap replications. Default is 9999.
#' @param l Desired 'block length' in the bootstrap. For the MBB, BWB and DWB, this is a genuine block length. For AWB, the blcok length is transformed into an autoregressive parameter via the formula 0.01^(1/l) as in Smeekes and Urbain (2014); this can be overwritten by setting ar_AWB directly. Default sets the block length as a function of the time series length T, via the rule l = 1.75*T^(1/3).
#' @param ar_AWB Autoregressive parameter used in the AWB bootstrap method (boot = "AWB"). Can be used to set the parameter directly rather than via the default link to the block length l.
#' @param p.min Minimum lag length in augmented Dickey-Fuller (ADF) regression. Default is 0.
#' @param p.max Maximum lag length in ADF regression. Default NULL uses the sample size-based rule 12*(T/100)^(1/4).
#' @param ic String for information criterion used to select the lag length in the ADF regression. Options are: AIC, BIC, MAIC, MBIC. Default is MAIC.
#' @param dc Numeric vector indicating the deterministic specification. Only relevant if union=FALSE. Options are: 0: no deterministics, 1: intercept only, 2: intercept and trend. Combinations thereof are allowed. Default is the union test, in which case this is not relevant. Note: program gives error if dc not set if union is FALSE. (Or change to constant as default?)
#' @param detr String vector indicating the type of detrending to be performed. Only relevant if union=FALSE. Options are OLS and/or QD (typically also called GLS, see Elliott, Rothenberg and Stock, 1996). Default is OLS.
#' @param ic.scale Logical indicator whether or not to use the rescaled information criteria of Cavaliere et al. (2015) (TRUE) or not (FALSE). Default is TRUE.
#' @param verbose Logical indictator wheter or not information on the outcome of the unit root test needs to be printed to the console. Default is FALSE.
#' @export
#' @return A list with the following components:
#'
#' rej_H0: Logical indicator of whether the null hypothesis of a unit root is rejected (TRUE) or not (FALSE).
#'
#' ADF_tests: Details on the unit root tests: value of the test statistics and p-values.
#'
#' The output components are arranged per time series, type of deterministic component and detrending method.
#' #' @references Cavaliere, G., Phillips, P. C. B., Smeekes, S., and Taylor, A. M. R. (2015), Lag length selection for unit root tests in the presence of nonstationary volatility, Econometric Reviews, 34(4), 512-536.
#' @references Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996), Efficient tests for an autoregressive unit root, Econometrica, 64(4), 813-836.
#' @references Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild bootstrap inference for nonparametric trends. \emph{Journal of Econometrics}, 214(1), 81-109.
#' @references Harvey, D. I., Leybourne, S. J., and Taylor, A. M. R. (2012), Testing for unit roots in the presence of uncertainty over both the trend and initial condition, Journal of Econometrics, 169(2), 188-195.
#' @references Moon, H. R. and Perron, B. (2012), Beyond panel unit root tests: Using multiple testing to determine the non stationarity properties of individual series in a panel. Journal of Econometrics, 169(1), 29-33.
#' @references Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008), Bootstrap unit root tests: Comparison and extensions, Journal of Time Series Analysis, 29(1), 371-401.
#' @references Rho, Y. and Shao, X. (2019), Bootstrap-assisted unit root testing with piecewise locally stationary errors, Econometric Theory, 35(1), 142-166.
#' @references Shao, X. (2010), The dependent wild bootstrap, Journal of the American Statistical Association, 105(489), 218-235.
#' @references Shao, X. (2011), A bootstrap-assisted spectral test of white noise under unknown dependence, Journal of Econometrics, 162, 213-224.
#' @references Smeekes, S. (2015), Bootstrap sequential tests to determine the order of integration of individual units in a time series panel, Journal of Time Series Analysis, 36(3), 398-415.
#' @references Smeekes, S. and Taylor, A.M.R. (2012), Bootstrap union tests for unit roots in the presence of nonstationary volatility, Econometric Theory, 28(2), 422-456.
#' @references Smeekes, S. and Urbain, J.-P. (2014a), A multivariate invariance principle for modified wild bootstrap methods with an application to unit root testing (GSBE Research Memorandum No. RM/14/008), Maastricht University
#' @references Smeekes, S. and Urbain, J.-P. (2014b), On the applicability of the sieve bootstrap in time series panels, Oxford Bulletin of Economics and Statistics, 76(1), 139-151.
#' @examples
#' # boot_df on GDP_BE
#' GDP_BE_df <- boot_df(eurostat[, 1], B = 399, dc = 2, detr = "OLS", verbose = TRUE)
boot_df <- function(y, level = 0.05, boot = "MBB", B = 9999, l = NULL, ar_AWB = NULL, p.min = 0,
                    p.max = NULL, ic = "MAIC", dc = 1, detr = "OLS", ic.scale = TRUE, verbose = FALSE){

  if (verbose) {
    cat("Bootstrap DF Test with", boot, "bootstrap method.\n")
  }
  out <- iADFtest(y, level = level, boot = boot, B = B, l = l, ar_AWB = ar_AWB, union = FALSE, p.min = p.min, p.max = p.max,
                       ic = ic, dc = dc, detr = detr, ic.scale = ic.scale, verbose = verbose)
  return(out)
}

#' Bootstrap Union Test for Unit Roots
#' @param y Vector with time series data to be tested for unit roots. Data may also be in a time series format (e.g. ts, zoo or xts).
#' @param level Desired significance level of the unit root test. Default is 0.05.
#' @param boot String for bootstrap method to be used. Options are MBB - moving blocks bootstrap (Moon and Perron, 2012; Smeekes, 2015), BWB - block wild bootstrap (Shao, 2011), DWB - dependent wild bootstrap (Shao, 2010; Smeekes and Urbain, 2014a; Rho and Shao, 2019), AWB - autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2018), SB - sieve bootstrap (Palm, Smeekes and Urbain, 2008; Smeekes and Urbain, 2014a), SWB - sieve wild boostrap (Smeekes and Taylor, 2014b). Default is MBB.
#' @param B Number of bootstrap replications. Default is 9999.
#' @param l Desired 'block length' in the bootstrap. For the MBB, BWB and DWB, this is a genuine block length. For AWB, the blcok length is transformed into an autoregressive parameter via the formula 0.01^(1/l) as in Smeekes and Urbain (2014); this can be overwritten by setting ar_AWB directly. Default sets the block length as a function of the time series length T, via the rule l = 1.75*T^(1/3).
#' @param ar_AWB Autoregressive parameter used in the AWB bootstrap method (boot = "AWB"). Can be used to set the parameter directly rather than via the default link to the block length l.
#' @param p.min Minimum lag length in augmented Dickey-Fuller (ADF) regression. Default is 0.
#' @param p.max Maximum lag length in ADF regression. Default NULL uses the sample size-based rule 12*(T/100)^(1/4).
#' @param ic String for information criterion used to select the lag length in the ADF regression. Options are: AIC, BIC, MAIC, MBIC. Default is MAIC.
#' @param ic.scale Logical indicator whether or not to use the rescaled information criteria of Cavaliere et al. (2015) (TRUE) or not (FALSE). Default is TRUE.
#' @param verbose Logical indictator wheter or not information on the outcome of the unit root test needs to be printed to the console. Default is FALSE.
#' @export
#' @return A list with the following components:
#'
#' rej_H0: Logical indicator of whether the null hypothesis of a unit root is rejected (TRUE) or not (FALSE).
#'
#' ADF_tests: Details on the unit root tests: value of the test statistics and p-values.
#'
#' The output components are arranged per time series
#' @references Cavaliere, G., Phillips, P. C. B., Smeekes, S., and Taylor, A. M. R. (2015), Lag length selection for unit root tests in the presence of nonstationary volatility, Econometric Reviews, 34(4), 512-536.
#' @references Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996), Efficient tests for an autoregressive unit root, Econometrica, 64(4), 813-836.
#' @references Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild bootstrap inference for nonparametric trends. \emph{Journal of Econometrics}, 214(1), 81-109.
#' @references Harvey, D. I., Leybourne, S. J., and Taylor, A. M. R. (2012), Testing for unit roots in the presence of uncertainty over both the trend and initial condition, Journal of Econometrics, 169(2), 188-195.
#' @references Moon, H. R. and Perron, B. (2012), Beyond panel unit root tests: Using multiple testing to determine the non stationarity properties of individual series in a panel. Journal of Econometrics, 169(1), 29-33.
#' @references Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008), Bootstrap unit root tests: Comparison and extensions, Journal of Time Series Analysis, 29(1), 371-401.
#' @references Rho, Y. and Shao, X. (2019), Bootstrap-assisted unit root testing with piecewise locally stationary errors, Econometric Theory, 35(1), 142-166.
#' @references Shao, X. (2010), The dependent wild bootstrap, Journal of the American Statistical Association, 105(489), 218-235.
#' @references Shao, X. (2011), A bootstrap-assisted spectral test of white noise under unknown dependence, Journal of Econometrics, 162, 213-224.
#' @references Smeekes, S. (2015), Bootstrap sequential tests to determine the order of integration of individual units in a time series panel, Journal of Time Series Analysis, 36(3), 398-415.
#' @references Smeekes, S. and Taylor, A.M.R. (2012), Bootstrap union tests for unit roots in the presence of nonstationary volatility, Econometric Theory, 28(2), 422-456.
#' @references Smeekes, S. and Urbain, J.-P. (2014a), A multivariate invariance principle for modified wild bootstrap methods with an application to unit root testing (GSBE Research Memorandum No. RM/14/008), Maastricht University
#' @references Smeekes, S. and Urbain, J.-P. (2014b), On the applicability of the sieve bootstrap in time series panels, Oxford Bulletin of Economics and Statistics, 76(1), 139-151.
#' @examples
#' # boot_union on GDP_BE
#' GDP_BE_df <- boot_union(eurostat[, 1], B = 399, verbose = TRUE)
boot_union <- function(y, level = 0.05, boot = "MBB", B = 9999, l = NULL, ar_AWB = NULL, p.min = 0,
                       p.max = NULL, ic = "MAIC", ic.scale = TRUE, verbose = FALSE){

  if (verbose) {
    cat("Bootstrap Test with", boot, "bootstrap method.\n")
  }
  if (NCOL(y) > 1) {
    stop("Function takes single time series only. Use one of the functions suited for multivariate time series.")
  }
  out <- iADFtest(y, level = level, boot = boot, B = B, l = l, ar_AWB = ar_AWB, union = TRUE, p.min = p.min, p.max = p.max,
                  ic = ic, ic.scale = ic.scale, verbose = verbose)
  return(out)
}

#' False Discovery Rate controlled Unit Root Tests
#' @param y A (TxN)-matrix of N time series with T observations to be tested for unit roots. Data may also be in a time series format (e.g. ts, zoo or xts).
#' @param level Desired significance level of the unit root test. Default is 0.05.
#' @param boot String for bootstrap method to be used. Options are MBB - moving blocks bootstrap (Moon and Perron, 2012; Smeekes, 2015), BWB - block wild bootstrap (Shao, 2011), DWB - dependent wild bootstrap (Shao, 2010; Smeekes and Urbain, 2014a; Rho and Shao, 2019), AWB - autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2018), SB - sieve bootstrap (Palm, Smeekes and Urbain, 2008; Smeekes and Urbain, 2014a), SWB - sieve wild boostrap (Smeekes and Taylor, 2014b). Default is MBB.
#' @param B Number of bootstrap replications. Default is 9999.
#' @param l Desired 'block length' in the bootstrap. For the MBB, BWB and DWB, this is a genuine block length. For AWB, the blcok length is transformed into an autoregressive parameter via the formula 0.01^(1/l) as in Smeekes and Urbain (2014); this can be overwritten by setting ar_AWB directly. Default sets the block length as a function of the time series length n, via the rule l = 1.75*n^(1/3).
#' @param ar_AWB Autoregressive parameter used in the AWB bootstrap method (boot = "AWB"). Can be used to set the parameter directly rather than via the default link to the block length l.
#' @param union Logical indicator whether or not to use bootstrap union tests (TRUE) or not (FALSE), see Harvey et al. (2012), Smeekes and Taylor (2012). Default is TRUE.
#' @param p.min Minimum lag length in augmented Dickey-Fuller (ADF) regression. Default is 0.
#' @param p.max Maximum lag length in ADF regression. Default NULL uses the sample size-based rule 12*(T/100)^(1/4).
#' @param ic String for information criterion used to select the lag length in the ADF regression. Options are: AIC, BIC, MAIC, MBIC. Default is MAIC.
#' @param dc Numeric vector indicating the deterministic specification. Only relevant if union=FALSE. Options are: 0: no deterministics, 1: intercept only, 2: intercept and trend. Combinations thereof are allowed. Default is 1.
#' @param detr String vector indicating the type of detrending to be performed. Only relevant if union=FALSE. Options are OLS and/or QD (typically also called GLS, see Elliott, Rothenberg and Stock, 1996). Default is OLS.
#' @param ic.scale Logical indicator whether or not to use the rescaled information criteria of Cavaliere et al. (2015) (TRUE) or not (FALSE). Default is TRUE.
#' @param verbose Logical indictator wheter or not information on the outcome of the unit root test needs to be printed to the console. Default is FALSE.
#' @export
#' @return A list with the following components:
#'
#' rej_H0: Logical indicator of whether the null hypothesis of a unit root is rejected (TRUE) or not (FALSE).
#'
#' FDR_sequence: Details on the unit root tests: value of the test statistics and critical values.
#'
#' If a union test is used, the output components are arranged per time series. If no union test is used, the output components are arranged per time series, type of deterministic component and detrending method.
#'
#' @references Cavaliere, G., Phillips, P. C. B., Smeekes, S., and Taylor, A. M. R. (2015), Lag length selection for unit root tests in the presence of nonstationary volatility, Econometric Reviews, 34(4), 512-536.
#' @references Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996), Efficient tests for an autoregressive unit root, Econometrica, 64(4), 813-836.
#' @references Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild bootstrap inference for nonparametric trends. \emph{Journal of Econometrics}, 214(1), 81-109.
#' @references Harvey, D. I., Leybourne, S. J., and Taylor, A. M. R. (2012), Testing for unit roots in the presence of uncertainty over both the trend and initial condition, Journal of Econometrics, 169(2), 188-195.
#' @references Moon, H. R. and Perron, B. (2012), Beyond panel unit root tests: Using multiple testing to determine the non stationarity properties of individual series in a panel. Journal of Econometrics, 169(1), 29-33.
#' @references Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008), Bootstrap unit root tests: Comparison and extensions, Journal of Time Series Analysis, 29(1), 371-401.
#' @references Rho, Y. and Shao, X. (2019), Bootstrap-assisted unit root testing with piecewise locally stationary errors, Econometric Theory, 35(1), 142-166.
#' @references Romano, J. P., Shaikh, A. M. ad Wold, M. (2008), Control of the false discovery rate under dependence using the bootstrap and subsampling, Test, 17(3), 417-442.
#' @references Shao, X. (2010), The dependent wild bootstrap, Journal of the American Statistical Association, 105(489), 218-235.
#' @references Shao, X. (2011), A bootstrap-assisted spectral test of white noise under unknown dependence, Journal of Econometrics, 162, 213-224.
#' @references Smeekes, S. (2015), Bootstrap sequential tests to determine the order of integration of individual units in a time series panel, Journal of Time Series Analysis, 36(3), 398-415.
#' @references Smeekes, S. and Taylor, A.M.R. (2012), Bootstrap union tests for unit roots in the presence of nonstationary volatility, Econometric Theory, 28(2), 422-456.
#' @references Smeekes, S. and Urbain, J.-P. (2014a), A multivariate invariance principle for modified wild bootstrap methods with an application to unit root testing (GSBE Research Memorandum No. RM/14/008), Maastricht University
#' @references Smeekes, S. and Urbain, J.-P. (2014b), On the applicability of the sieve bootstrap in time series panels, Oxford Bulletin of Economics and Statistics, 76(1), 139-151.
#' @examples
#' # bFDRtest on GDPC1 and T5YFFM
#' two_series_bFDRtest <- bFDRtest(eurostat[, 1:2], boot = "MBB", B = 399,  verbose = TRUE)
bFDRtest <- function(y, level = 0.05,  boot = "MBB", l = NULL, ar_AWB = NULL, B = 9999, union = TRUE, p.min = 0, p.max = NULL, ic = "MAIC", dc = NULL, detr = NULL, ic.scale = TRUE, verbose = FALSE){

  inputs <- generate_inputs(y = y, BSQT_test = FALSE, iADF_test = FALSE, level = level, boot = boot, B = B, union = union,
                            p.min = p.min, p.max = p.max, ic = ic, dc = dc, detr = detr, q = NULL, l = l,
                            ic.scale = ic.scale, h.rs = 0.1, k_DWB = "k.TBB", ar_AWB = ar_AWB)

  if (!is.null(colnames(y))) {
    var_names <- colnames(y)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(y))
  }

  if (union) { # Union Tests
    bFDRout <- FDR_cpp(test_i = inputs$test_stats, t_star = inputs$test_stats_star, level = inputs$level)
    rej_H0 <- matrix(bFDRout$rej_H0 == 1, nrow = NCOL(y))
    FDR_seq <- bFDRout$FDR_Tests[, -1, drop = FALSE]

    rownames(rej_H0) <- var_names
    rownames(FDR_seq) <- var_names[bFDRout$FDR_Tests[, 1]]

    colnames(rej_H0) <- "Reject H0"
    colnames(FDR_seq) <- c("test statistic", "critical value")

    if (verbose) {
      p_hat <- sum(rej_H0)
      if (p_hat > 1){
        cat(paste("There are ", p_hat, " stationary time series, namely: ", paste(var_names[rej_H0], collapse = " "), "\n", sep = ""))
      } else if (p_hat == 1) {
          cat(paste("There is ", p_hat, " stationary time series, namely: ", paste(var_names[rej_H0], collapse = " "), "\n", sep = ""))
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

    detr_dc_names <- paste(rep(detr_names, each = length(inputs$dc)), ", ", rep(dc_names, length(inputs$detr)), sep = "")
    rej_H0 <- matrix(nrow = NCOL(y), ncol = length(inputs$dc)*length(inputs$detr))
    rownames(rej_H0) <- var_names
    colnames(rej_H0) <- detr_dc_names
    FDR_seq <- vector("list", length(inputs$dc)*length(inputs$detr))
    names(FDR_seq) <- detr_dc_names
    for (i in 1:nrow(inputs$tests_i)) {
      bFDRout <- FDR_cpp(test_i = matrix(inputs$tests_i[i, ], nrow = 1), t_star = inputs$t_star[ , i,], level = inputs$level)
      rej_H0[, i] <- bFDRout$rej_H0 == 1
      FDR_seq[[i]] <- bFDRout$FDR_Tests[, -1]
      rownames(FDR_seq[[i]]) <- var_names[bFDRout$FDR_Tests[, 1]]
      colnames(FDR_seq[[i]]) <- c("test statistic", "critical value")

      if (verbose) {
        p_hat <- sum(rej_H0[, i])
        cat(paste(c(rep("-", 40), "\n"), sep = "", collapse = ""))
        cat(paste("Type of unit root test performed: ", detr_dc_names[i], "\n", sep = ""))
        if (p_hat > 1){
          cat(paste("There are ", p_hat, " stationary time series, namely: ", paste(var_names[rej_H0[, i]], collapse = " "), ".\n", sep = ""))
        } else if (p_hat == 1) {
          cat(paste("There is ", p_hat, " stationary time series, namely: ", paste(var_names[rej_H0[, i]], collapse = " "), ".\n", sep = ""))
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


#' Bootstrap Sequential Panel Test
#' @param y A (TxN)-matrix of N time series with T observations to be tested for unit roots. Data may also be in a time series format (e.g. ts, zoo or xts).
#' @param level Desired significance level of the unit root test. Default is 0.05.
#' @param boot String for bootstrap method to be used. Options are MBB - moving blocks bootstrap (Moon and Perron, 2012; Smeekes, 2015), BWB - block wild bootstrap (Shao, 2011), DWB - dependent wild bootstrap (Shao, 2010; Smeekes and Urbain, 2014a; Rho and Shao, 2019), AWB - autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2018), SB - sieve bootstrap (Palm, Smeekes and Urbain, 2008; Smeekes and Urbain, 2014a), SWB - sieve wild boostrap (Smeekes and Taylor, 2014b). Default is MBB.
#' @param B Number of bootstrap replications. Default is 9999.
#' @param l Desired 'block length' in the bootstrap. For the MBB, BWB and DWB, this is a genuine block length. For AWB, the blcok length is transformed into an autoregressive parameter via the formula 0.01^(1/l) as in Smeekes and Urbain (2014); this can be overwritten by setting ar_AWB directly. Default sets the block length as a function of the time series length n, via the rule l = 1.75*n^(1/3).
#' @param ar_AWB Autoregressive parameter used in the AWB bootstrap method (boot = "AWB"). Can be used to set the parameter directly rather than via the default link to the block length l.
#' @param union Logical indicator whether or not to use bootstrap union tests (TRUE) or not (FALSE), see Harvey et al. (2012), Smeekes and Taylor (2012). Default is TRUE.
#' @param p.min Minimum lag length in augmented Dickey-Fuller (ADF) regression. Default is 0.
#' @param p.max Maximum lag length in ADF regression. Default NULL uses the sample size-based rule 12*(T/100)^(1/4).
#' @param ic String for information criterion used to select the lag length in the ADF regression. Options are: AIC, BIC, MAIC, MBIC. Default is MAIC.
#' @param dc Numeric vector indicating the deterministic specification. Only relevant if union=FALSE. Options are: 0: no deterministics, 1: intercept only, 2: intercept and trend. Combinations thereof are allowed. Default is 1.
#' @param detr String vector indicating the type of detrending to be performed. Only relevant if union=FALSE. Options are OLS and/or QD (typically also called GLS, see Elliott, Rothenberg and Stock, 1996). Default is OLS.
#' @param q Numeric vector quantiles to be tested Default is to test each unit sequentially.
#' @param ic.scale Logical indicator whether or not to use the rescaled information criteria of Cavaliere et al. (2015) (TRUE) or not (FALSE). Default is TRUE.
#' @param verbose Logical indictator wheter or not information on the outcome of the unit root test needs to be printed to the console. Default is FALSE.
#' @export
#' @return A list with the following components:
#'
#' rej_H0: Logical indicator of whether the null hypothesis of a unit root is rejected (TRUE) or not (FALSE).
#'
#' BSQT_sequence: Details on the unit root tests: outcome of the sequential steps, value of the test statistics and p-values.
#'
#' If a union test is used, the output components are arranged per time series. If no union test is used, the output components are arranged per time series, type of deterministic component and detrending method.
#'
#' @references Cavaliere, G., Phillips, P. C. B., Smeekes, S., and Taylor, A. M. R. (2015), Lag length selection for unit root tests in the presence of nonstationary volatility, Econometric Reviews, 34(4), 512-536.
#' @references Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996), Efficient tests for an autoregressive unit root, Econometrica, 64(4), 813-836.
#' @references Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild bootstrap inference for nonparametric trends. \emph{Journal of Econometrics}, 214(1), 81-109.
#' @references Harvey, D. I., Leybourne, S. J., and Taylor, A. M. R. (2012), Testing for unit roots in the presence of uncertainty over both the trend and initial condition, Journal of Econometrics, 169(2), 188-195.
#' @references Moon, H. R. and Perron, B. (2012), Beyond panel unit root tests: Using multiple testing to determine the non stationarity properties of individual series in a panel. Journal of Econometrics, 169(1), 29-33.
#' @references Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008), Bootstrap unit root tests: Comparison and extensions, Journal of Time Series Analysis, 29(1), 371-401.
#' @references Rho, Y. and Shao, X. (2019), Bootstrap-assisted unit root testing with piecewise locally stationary errors, Econometric Theory, 35(1), 142-166.
#' @references Shao, X. (2010), The dependent wild bootstrap, Journal of the American Statistical Association, 105(489), 218-235.
#' @references Shao, X. (2011), A bootstrap-assisted spectral test of white noise under unknown dependence, Journal of Econometrics, 162, 213-224.
#' @references Smeekes, S. (2015), Bootstrap sequential tests to determine the order of integration of individual units in a time series panel, Journal of Time Series Analysis, 36(3), 398-415.
#' @references Smeekes, S. and Taylor, A.M.R. (2012), Bootstrap union tests for unit roots in the presence of nonstationary volatility, Econometric Theory, 28(2), 422-456.
#' @references Smeekes, S. and Urbain, J.-P. (2014a), A multivariate invariance principle for modified wild bootstrap methods with an application to unit root testing (GSBE Research Memorandum No. RM/14/008), Maastricht University
#' @references Smeekes, S. and Urbain, J.-P. (2014b), On the applicability of the sieve bootstrap in time series panels, Oxford Bulletin of Economics and Statistics, 76(1), 139-151.
#' @examples
#' # BSQTtest on GDP_BE and GDP_DE
#' two_series_BSQTtest <- BSQTtest(eurostat[, 1:2], boot = "MBB", B = 399,  verbose = TRUE)
BSQTtest <- function(y, level = 0.05,  boot = "MBB", B = 9999, l = NULL, ar_AWB = NULL, union = TRUE, p.min = 0, p.max = NULL,
                     ic = "MAIC", dc = NULL, detr = NULL, q = 0:NCOL(y), ic.scale = TRUE, verbose = FALSE){

  inputs <- generate_inputs(y = y, BSQT_test = TRUE, iADF_test = FALSE, level = level, boot = boot, B = B, union = union,
                            p.min = p.min, p.max = p.max, ic = ic, dc = dc, detr = detr, q = q, l = l,
                            ic.scale = ic.scale, h.rs = 0.1, k_DWB = "k.TBB", ar_AWB = ar_AWB)

  if (!is.null(colnames(y))) {
    var_names <- colnames(y)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(y))
  }

  if (union) { # Union Tests
    BSQTout <- BSQT_cpp(pvec = inputs$p_vec, test_i = inputs$test_stats, t_star = inputs$test_stats_star, level = inputs$level)
    rej_H0 <- matrix(BSQTout$rej_H0 == 1, nrow = NCOL(y))
    BSQT_seq <- BSQTout$BSQT_steps[, -3, drop = FALSE]

    rownames(rej_H0) <- var_names
    rownames(BSQT_seq) <- paste("Step", 1:nrow(BSQT_seq))
    colnames(rej_H0) <- "Reject H0"
    colnames(BSQT_seq) <- c("Unit H0", "Unit H1", "Test statistic", "p-value")

    if (verbose) {
      p_hat <- sum(rej_H0)
      if (p_hat > 1) {
        cat(paste("There are ", p_hat, " stationary time series, namely: ", paste(var_names[rej_H0], collapse = " "), ".\n", sep = ""))
      } else {
        if (p_hat == 1) {
          cat(paste("There is ", p_hat, " stationary time series, namely: ", paste(var_names[rej_H0], collapse = " "), ".\n", sep = ""))
        } else {
          cat(paste("There are ", p_hat, " stationary time series.\n", sep = ""))
        }
      }
      cat("Details of the BSQT ssquential tests:\n")
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

    detr_dc_names <- paste(rep(detr_names, each = length(inputs$dc)), ", ", rep(dc_names, length(inputs$detr)), sep = "")
    rej_H0 <- matrix(nrow = NCOL(y), ncol = length(inputs$dc)*length(inputs$detr))
    rownames(rej_H0) <- var_names
    colnames(rej_H0) <- detr_dc_names
    BSQT_seq <- vector("list", length(inputs$dc)*length(inputs$detr))
    names(BSQT_seq) <- detr_dc_names
    for (i in 1:nrow(inputs$tests_i)) {
      BSQTout <- BSQT_cpp(pvec = inputs$p_vec, test_i = matrix(inputs$tests_i[i, ], nrow = 1), t_star = inputs$t_star[ , i,], level = inputs$level)
      rej_H0[, i] <- BSQTout$rej_H0 == 1
      BSQT_seq[[i]] <- BSQTout$BSQT_steps[, -3, drop = FALSE]
      rownames(BSQT_seq[[i]]) <- paste("Step", 1:nrow(BSQT_seq[[i]]))
      colnames(BSQT_seq[[i]]) <- c("Unit H0", "Unit H1", "Test statistic", "p-value")

      if (verbose) {
        p_hat <- sum(rej_H0[, i])
        cat(paste(c(rep("-", 40), "\n"), sep = "", collapse = ""))
        cat(paste("Type of unit root test performed:", detr_dc_names[i], "\n", sep = ""))
        if (p_hat > 1){
          cat(paste("There are ", p_hat, "stationary time series, namely:", paste(var_names[rej_H0[, i]], collapse = " "), ".\n", sep = ""))
        } else if (p_hat == 1) {
          cat(paste("There is ", p_hat, "stationary time series, namely:", paste(var_names[rej_H0[, i]], collapse = " "), ".\n", sep = ""))
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
#' @param y A (TxN)-matrix of N time series with T observations to be tested for unit roots. Data may also be in a time series format (e.g. ts, zoo or xts).
#' @param level Desired significance level of the unit root test. Default is 0.05.
#' @param boot String for bootstrap method to be used. Options are MBB - moving blocks bootstrap (Moon and Perron, 2012; Smeekes, 2015), BWB - block wild bootstrap (Shao, 2011), DWB - dependent wild bootstrap (Shao, 2010; Smeekes and Urbain, 2014a; Rho and Shao, 2019), AWB - autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2018), SB - sieve bootstrap (Palm, Smeekes and Urbain, 2008; Smeekes and Urbain, 2014a), SWB - sieve wild boostrap (Smeekes and Taylor, 2014b). Default is MBB.
#' @param B Number of bootstrap replications. Default is 9999.
#' @param l Desired 'block length' in the bootstrap. For the MBB, BWB and DWB, this is a genuine block length. For AWB, the blcok length is transformed into an autoregressive parameter via the formula 0.01^(1/l) as in Smeekes and Urbain (2014); this can be overwritten by setting ar_AWB directly. Default sets the block length as a function of the time series length n, via the rule l = 1.75*n^(1/3).
#' @param ar_AWB Autoregressive parameter used in the AWB boostrap method (boot = "AWB"). Can be used to set the parameter directly rather than via the default link to the block length l.
#' @param union Logical indicator whether or not to use bootstrap union tests (TRUE) or not (FALSE), see Harvey et al. (2012), Smeekes and Taylor (2012). Default is TRUE.
#' @param p.min Minimum lag length in augmented Dickey-Fuller (ADF) regression. Default is 0.
#' @param p.max Maximum lag length in ADF regression. Default NULL uses the sample size-based rule 12*(T/100)^(1/4).
#' @param ic String for information criterion used to select the lag length in the ADF regression. Options are: AIC, BIC, MAIC, MBIC. Default is MAIC.
#' @param dc Numeric vector indicating the deterministic specification. Only relevant if union=FALSE. Options are: 0: no deterministics, 1: intercept only, 2: intercept and trend. Combinations thereof are allowed. Default is 1.
#' @param detr String vector indicating the type of detrending to be performed. Only relevant if union=FALSE. Options are OLS and/or QD (typically also called GLS, see Elliott, Rothenberg and Stock, 1996). Default is OLS.
#' @param ic.scale Logical indicator whether or not to use the rescaled information criteria of Cavaliere et al. (2015) (TRUE) or not (FALSE). Default is TRUE.
#' @param verbose Logical indictator wheter or not information on the outcome of the unit root test needs to be printed to the console. Default is FALSE.
#' @export
#' @return If a union test is used, the test statistic and p-value are returned. If no union test is used, the test statistics and p-values are reported per type of deterministic component and detrending method.
#' @references Cavaliere, G., Phillips, P. C. B., Smeekes, S., and Taylor, A. M. R. (2015), Lag length selection for unit root tests in the presence of nonstationary volatility, Econometric Reviews, 34(4), 512-536.
#' @references Elliott, G., Rothenberg, T.J., and Stock, J.H. (1996), Efficient tests for an autoregressive unit root, Econometrica, 64(4), 813-836.
#' @references Friedrich, M., Smeekes, S. and Urbain, J.-P. (2020). Autoregressive wild bootstrap inference for nonparametric trends. \emph{Journal of Econometrics}, 214(1), 81-109.
#' @references Harvey, D. I., Leybourne, S. J., and Taylor, A. M. R. (2012), Testing for unit roots in the presence of uncertainty over both the trend and initial condition, Journal of Econometrics, 169(2), 188-195.
#' @references Moon, H. R. and Perron, B. (2012), Beyond panel unit root tests: Using multiple testing to determine the non stationarity properties of individual series in a panel. Journal of Econometrics, 169(1), 29-33.
#' @references Palm, F.C., Smeekes, S. and Urbain, J.-P. (2008), Bootstrap unit root tests: Comparison and extensions, Journal of Time Series Analysis, 29(1), 371-401.
#' @references Palm, F. C., Smeekes, S., and Urbain, J.-.P. (2011), Cross-sectional dependence robust block bootstrap panel unit root tests, Journal of Econometrics, 163(1), 85-104.
#' @references Rho, Y. and Shao, X. (2019), Bootstrap-assisted unit root testing with piecewise locally stationary errors, Econometric Theory, 35(1), 142-166.
#' @references Shao, X. (2010), The dependent wild bootstrap, Journal of the American Statistical Association, 105(489), 218-235.
#' @references Shao, X. (2011), A bootstrap-assisted spectral test of white noise under unknown dependence, Journal of Econometrics, 162, 213-224.
#' @references Smeekes, S. (2015), Bootstrap sequential tests to determine the order of integration of individual units in a time series panel, Journal of Time Series Analysis, 36(3), 398-415.
#' @references Smeekes, S. and Taylor, A.M.R. (2012), Bootstrap union tests for unit roots in the presence of nonstationary volatility, Econometric Theory, 28(2), 422-456.
#' @references Smeekes, S. and Urbain, J.-P. (2014a), A multivariate invariance principle for modified wild bootstrap methods with an application to unit root testing (GSBE Research Memorandum No. RM/14/008), Maastricht University
#' @references Smeekes, S. and Urbain, J.-P. (2014b), On the applicability of the sieve bootstrap in time series panels, Oxford Bulletin of Economics and Statistics, 76(1), 139-151.
#' @examples
#' # paneltest on GDP_BE and GDP_DE
#' two_series_paneltest <- paneltest(eurostat[, 1:2], boot = "MBB", B = 399,  verbose = TRUE)
paneltest <- function(y, level = 0.05,  boot = "MBB", B = 9999, l = NULL, ar_AWB = NULL, union = TRUE, p.min = 0, p.max = NULL,
                      ic = "MAIC", dc = NULL, detr = NULL, ic.scale = TRUE, verbose = FALSE){

  inputs <- generate_inputs(y = y, BSQT_test = FALSE, iADF_test = FALSE, level = level, boot = boot, B = B, union = union,
                            p.min = p.min, p.max = p.max, ic = ic, dc = dc, detr = detr, q = NULL, l = l,
                            ic.scale = ic.scale, h.rs = 0.1, k_DWB = "k.TBB", ar_AWB = ar_AWB)

  if (union) { # Union Tests
    GM_test <- mean(inputs$test_stats)
    t_star <- rowMeans(inputs$test_stats_star)
    p_val <- mean(t_star < GM_test)
    out <- cbind("test statistic" = GM_test, "p-value" = p_val)

    if (verbose) {
      cat("Panel Bootstrap Group-Mean Union Test\n")
      if (p_val < level){
        cat(paste("The null hypothesis that all series have a unit root, is rejected at a significance level of ", level, ".\n", sep = ""))
      } else {
        cat(paste("The null hypothesis that all series have a unit root, is not rejected at a significance level of ", level, ".\n", sep = ""))
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

    detr_dc_names <- paste(rep(detr_names, each = length(inputs$dc)), ", ", rep(dc_names, length(inputs$detr)), sep = "")
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
          cat(paste(rownames(out)[i], ": The null hypothesis that all series have a unit root, is rejected at a significance level of ", level, ".\n", sep = ""))
        } else {
          cat(paste(rownames(out)[i], ": The null hypothesis that all series have a unit root, is not rejected at a significance level of ", level, ".\n", sep = ""))
        }
      }
      print(out)
    }
  }
  return(out)
}

#' Check Missing Values in Sample
#' @param X A (TxN)-matrix of N time series with T observations. Data may also be in a time series format (e.g. ts, zoo or xts).
#' @export
#' @return N-dimensional vector, for each series whether missing values are present (TRUE) or not (FALSE)
check_missing_insample_values <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  check <- apply(X, 2, function(x){any(diff((1:n)[!is.na(x)]) != 1)})
  return(check)
}

#' Find Non-Missing Subsamples
#' @param X A (TxN)-matrix of N time series with T observations. Data may also be in a time series format (e.g. ts, zoo or xts). Assumes a prior check on missing vaues in-sample has been done.
#' @export
#' @return A list with the following components: range: (2xN)-dimensional matrix containing the first and last non-missing observation in each column of X; all_equal logical value indicating whether all series have the same non-missing indices
find_nonmissing_subsample <- function(X) {
  n <- NROW(X)
  N <- NCOL(X)
  ind <- matrix(1:n, nrow = n, ncol = N)
  ind[is.na(X)] <- NA
  range_nonmiss <- matrix(apply(ind, 2, function(x){c(min(stats::na.omit(x)), max(stats::na.omit(x)))}), ncol = N)
  all_equal <- (max(range_nonmiss[1, ]) - min(range_nonmiss[1, ]) == 0) & (max(range_nonmiss[2, ]) - min(range_nonmiss[2, ]) > 0)
  return(list(range = range_nonmiss, all_equal = all_equal))
}
