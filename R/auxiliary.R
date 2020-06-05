#' Auxiliary Function (not accessible to users) to create all inputs for the unit root test functions.
#' @param y A T-dimensional vector or (TxN)-matrix with data to be tested for unit roots
#' @param BSQT_test Logical indicator whether or not to perform the Bootstrap Quantile Test (TRUE) or not (FALSE).
#' @param iADF_test Logical indicator whether or not to perform the individual ADF Tests (TRUE) or not (FALSE).
#' @param level Desired significance level of the unit root test. Default is 0.05.
#' @param boot String for bootstrap method to be used. Options are MBB - moving blocks bootstrap (Moon and Perron, 2012; Smeekes, 2015), BWB - block wild bootstrap (Shao, 2011), DWB - dependent wild bootstrap (Shao, 2010; Smeekes and Urbain, 2014a; Rho and Shao, 2019), AWB - autoregressive wild bootstrap (Smeekes and Urbain, 2014a; Friedrich, Smeekes and Urbain, 2018), SB - sieve bootstrap (Palm, Smeekes and Urbain, 2008; Smeekes and Urbain, 2014a), SWB - sieve wild boostrap (Smeekes and Taylor, 2014b). Default is MBB.
#' @param B Number of bootstrap replications. Default is 9999.
#' @param union Logical indicator whether or not to use bootstrap union tests (TRUE) or not (FALSE), see Harvey et al. (2012), Smeekes and Taylor (2012). Default is TRUE.
#' @param p.min Minimum lag length in augmented Dickey-Fuller (ADF) regression. Default is 0.
#' @param p.max Maximum lag length in ADF regression. Default NULL uses the sample size-based rule 12*(T/100)^(1/4).
#' @param ic String for information criterion used to select the lag length in the ADF regression. Options are: AIC, BIC, MAIC, MBIC. Default is MAIC.
#' @param dc Numeric vector indicating the deterministic specification. Only relevant if union=FALSE. Options are: 0: no deterministics, 1: intercept only, 2: intercept and trend. Combinations thereof are allowed. Default is 1.
#' @param detr String vector indicating the type of detrending to be performed. Only relevant if union=FALSE. Options are OLS and/or QD (typically also called GLS, see Elliott, Rothenberg and Stock, 1996). Default is OLS.
#' @param q Numeric vector quantiles to be tested Default is to test each unit sequentially.
#' @param l Block length used in MBB and BWB.
#' @param ic.scale Logical indicator whether or not to use the rescaled information criteria of Cavaliere et al. (2015) (TRUE) or not (FALSE). Default is TRUE.
#' @param h.rs Bandwidth used in rescaled information criteria
#' @param k_DWB Kernel used in DWB method
#' @param ar_AWB AR parameter used in AWB method.
generate_inputs <- function(y, BSQT_test, iADF_test, level, boot, B, union, p.min, p.max, ic, dc, detr, q, l, ic.scale, h.rs, k_DWB, ar_AWB){
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