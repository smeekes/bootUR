#' Auxiliary Function (not accessible to users) to create all inputs for the unit root test functions.
#' @param y A \eqn{T}-dimensional vector or a (\eqn{T} x \eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations to be tested for unit roots. Data may also be in a time series format (e.g. \verb{ts}, \verb{zoo} or \verb{xts}).
#' @param BSQT_test Logical indicator whether or not to perform the Bootstrap Quantile Test (\verb{TRUE}) or not (\verb{FALSE}).
#' @param iADF_test Logical indicator whether or not to perform the individual ADF Tests (\verb{TRUE}) or not (\verb{FALSE}).
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
#' @param dc Numeric vector indicating the deterministic specification. Only relevant if \verb{union = FALSE}. Options are
#' \describe{
#' \item{\emph{0}}{no deterministics;}
#' \item{\emph{1}}{intercept only;}
#' \item{\emph{2}}{intercept and trend.}
#' Combinations thereof are allowed. Default is the union test (\verb{union = TRUE}), in which case this is not relevant. If \verb{union = FALSE}, the default is adding an intercept (a warning is given).
#' }
#' @param detr String vector indicating the type of detrending to be performed. Only relevant if \verb{unionunion = FALSE}. Options are: \verb{"OLS"} and/or \verb{"QD"} (typically also called GLS, see Elliott, Rothenberg and Stock, 1996). Default is \verb{"OLS"}.
#' @param ic.scale Logical indicator whether or not to use the rescaled
#' @param q Numeric vector of quantiles to be tested. Default is to test each unit sequentially.
#' @param h.rs Bandwidth used in rescaled information criteria.
#' @param k_DWB Kernel used in DWB bootstrap method.
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
generate_inputs <- function(y, BSQT_test, iADF_test, level, boot, B, l, ar_AWB, union, p.min, p.max, ic, dc, detr, ic.scale, q, h.rs, k_DWB){
  if (level * (B + 1) < 1) {
    stop("Bootstrap iterations B too low to perform test at desired significance level.")
  }
  y <- as.matrix(y)

  # Dimensions
  n <- nrow(y)
  N <- ncol(y)

  # Check for missing values or unbalanced panels (MBB, SB)
  check_missing <- check_missing_insample_values(y)
  if (any(check_missing)) {
    stop("Missing values detected inside samples. Series may only contain missing values at the beginning or end for wild bootstrap methods.")
  } else if (anyNA(y)) {
    if (boot %in% c("MBB", "SB")) {
      if (!iADF_test) {
        stop("Resampling-based bootstraps MBB and SB cannot handle missing values.")
      } else {
        warning("Missing values cause resampling bootstrap to be executed for each time series individually.")
      }
    }
    check_nonmiss <- find_nonmissing_subsample(y)
    range_nonmiss <- check_nonmiss$range - 1
    joint <- check_nonmiss$all_equal
  } else {
    range_nonmiss <- matrix(c(0, n - 1), nrow = 2, ncol = N)
    joint <- TRUE
  }

  # Checks on inputs
  if (!(boot %in% c("MBB", "BWB", "DWB", "AWB", "SB", "SWB")) | length(boot) > 1) {
    stop("The argument boot should be equal to either MBB, BWB, DWB, AWB, SB or SWB")
  } else if (boot %in% c("SB", "SWB") & !iADF_test) {
    warning("SB and SWB bootstrap only recommended for iADF_test; see help for details.")
  }
  boot <- 1*(boot=="MBB") + 2*(boot=="BWB") + 3*(boot=="DWB") + 4*(boot=="AWB") + 5*(boot=="SB") + 6*(boot=="SWB")

  if (any(!is.element(ic, c("AIC", "BIC", "MAIC", "MBIC"))) | length(ic) > 1) {
    stop("The argument ic should be equal to either AIC, BIC, MAIC, MBIC)")
  }
  ic <- 1*(ic=="AIC") + 2*(ic=="BIC") + 3*(ic=="MAIC") + 4*(ic=="MBIC")

  # Bootstrap Union Tests: Settings
  if (union) {
    if (!is.null(dc)) {
      warning("Deterministic specification in argument dc is ignored, as union test is applied.")
    }
    if (!is.null(detr)) {
      warning("Detrending method in argument detr is ignored, as union test is applied.")
    }
    dc <- c(1,2)
    dc.boot <- 2
    detr_int <- 1:2
  } else {
    if (is.null(dc)) {
      stop("No deterministic specification set. Set dc to 0, 1 or 2, or set union to TRUE for union test.")
    } else if (any(!is.element(dc, 0:2))){
      stop("The argument dc should only contain values 0, 1, and/or 2: (0: no deterministics, 1: intercept, 2: intercept and trend)")
    }
    dc <- sort(dc)
    dc.boot <- max(dc)
    if (is.null(detr)) {
      warning("No detrending specification set. Using OLS detrending.")
      detr <- "OLS"
    } else if(any(!is.element(detr, c("OLS", "QD")))) {
      stop("The argument detr should only contain the strings OLS and/or QD")
    }
    detr_int <- 1*(detr=="OLS") + 2*(detr=="QD")
    detr_int <- sort(detr_int)
  }

  if (is.null(l)) {
    l <- round(1.75*NROW(y)^(1/3))
  }

  if (is.null(ar_AWB)) {
    ar_AWB <- 0.01^(1/l)
  } else if (boot != 4){
    warning("Argument ar_AWB set, but AWB method not used. ar_AWB is ignored.")
  }

  # Defaults
  p_vec <- NULL
  if (BSQT_test) {
    if (is.numeric(q) && all(q == floor(q) & q >= 0)) {
      p_vec <- sort(unique(q))
    } else if (is.numeric(q) && all(q >= 0 & q <= 1)) {
      p_vec <- sort(unique(round(q*N)))
    } else {
      stop("Invalid input values for q: must be quantiles or positive integers")
    }
  }
  if(is.null(p.max)){
    p.max = round(12*(n/100)^(1/4))
  }
  s_DWB <- matrix(0, n, n)
  if (boot == 3) {
    # Self-convolution
    self.conv <- function(t, c) {
      y <- apply(as.matrix(t), c(1,2), FUN = function(x){stats::integrate(f.w, -1, 1, t = x, c = c)$value})
      return(y)
    }
    # Integrand
    f.w <- function(x, t, c) {
      y <- w.trap(x, c) * w.trap(x + abs(t), c)
      return(y)
    }
    # Trapezoid
    w.trap <- function(x, c) {
      y <- (x/c) * (x >= 0) * (x<c) + (x >= c) * (x <= 1-c) + ((1-x) / c) * (x > 1-c) * (x <= 1)
      return(y)
    }
    m <- self.conv(outer(1:n,1:n,"-") / l, 0.43)/ drop(self.conv(0, 0.43))
    ev <- eigen(m)
    e_va <- ev$values
    e_ve <- ev$vectors
    e_va <- diag((e_va > 1e-10) * e_va + (e_va <= 1e-10) * 1e-10)
    s_DWB <- e_ve%*%sqrt(e_va)
  }

  # ADF tests
  panel_est <- adf_panel_bootstrap_dgp_cpp(y = y, pmin = p.min, pmax = p.max, ic = ic, dc = dc.boot, QD = FALSE, trim = FALSE, ic_scale = ic.scale, h_rs = h.rs, range = range_nonmiss)
  u_boot <- panel_est$b_res
  u_boot[is.nan(u_boot)] <- NA
  res <- panel_est$res
  ar_est <- panel_est$par[-1, , drop = FALSE]
  t_star <- bootstrap_cpp(B = B, boot = boot, u = u_boot, e = res, l = l, s = s_DWB, ar = ar_AWB, ar_est = ar_est, y0 = matrix(0, ncol = N), pmin = p.min, pmax = p.max, ic = ic, dc = dc, detr = detr_int, ic_scale = ic.scale, h_rs = h.rs, range_nonmiss, joint)
   tests_i <- adf_tests_panel_cpp(y, pmin = p.min, pmax = p.max, ic = ic, dc = dc, detr = detr_int, ic_scale = ic.scale, h_rs = h.rs, range_nonmiss)

  if (union) {
    scaling <- scaling_factors_cpp(t_star, level)
    if (N > 1) {
      test_stats_star <- union_tests_cpp(t_star, scaling)
      test_stats <- union_tests_cpp(array(tests_i, dim = c(1, length(dc) * length(detr_int), N)), scaling)
    } else {
      test_stats_star <- union_test_cpp(t_star[, , 1], scaling)
      test_stats <- union_test_cpp(array(tests_i, dim = c(1, length(dc) * length(detr_int))), scaling)
    }
   } else {
    test_stats_star <- NULL
    test_stats <- NULL
  }
  out <- list("y" = y, "p_vec" = p_vec, "t_star" = t_star, "test_stats_star" = test_stats_star, "tests_i" = tests_i, "test_stats" = test_stats, "level" = level, "dc" = dc, "detr" = detr)

  return(out)
}
