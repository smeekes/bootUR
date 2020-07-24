#' Auxiliary Function (not accessible to users) to create all bootstrap statistics used to perform the unit root tests.
#' @inheritParams iADFtest
#' @param BSQT_test Logical indicator whether or not to perform the Bootstrap Quantile Test (\verb{TRUE}) or not (\verb{FALSE}).
#' @param iADF_test Logical indicator whether or not to perform the individual ADF Tests (\verb{TRUE}) or not (\verb{FALSE}).
#' @param level Desired significance level of the unit root test.
#' @param boot String for bootstrap method to be used. Options are
#' \describe{
#' \item{\verb{"MBB"}}{Moving blocks bootstrap;}
#' \item{\verb{"BWB"}}{Block wild bootstrap;}
#' \item{\verb{"DWB"}}{Dependent wild bootstrap;}
#' \item{\verb{"AWB"}}{Autoregressive wild bootstrap;}
#' \item{\verb{"SB"}}{Sieve bootstrap;}
#' \item{\verb{"SWB"}}{Sieve wild boostrap.}
#' }
#' @param l Desired 'block length' in the bootstrap. For the MBB, BWB and DWB boostrap, this is a genuine block length. For the AWB boostrap, the block length is transformed into an autoregressive parameter via the formula \eqn{0.01^(1/l)}; this can be overwritten by setting \verb{ar_AWB} directly. If NULL, sets the block length as a function of the time series length T, via the rule \eqn{l = 1.75 T^(1/3)}.
#' @param ar_AWB Autoregressive parameter used in the AWB bootstrap method (\verb{boot = "AWB"}). Can be used to set the parameter directly rather than via the default link to the block length l.
#' @param union Logical indicator whether or not to use bootstrap union tests (\verb{TRUE}) or not (\verb{FALSE}).
#' @param p_min Minimum lag length in the augmented Dickey-Fuller regression.
#' @param p_max Maximum lag length in the augmented Dickey-Fuller regression.
#' @param ic String for information criterion used to select the lag length in the augmented Dickey-Fuller regression. Options are: \verb{"AIC"}, \verb{"BIC"}, \verb{"MAIC"}, \verb{"MBIC}.
#' @param dc Numeric vector indicating the deterministic specification. Only relevant if \code{union = FALSE}. Options are (combinations of): \code{0} - no deterministics; \code{1} - intercept only; \code{2} - intercept and trend. If \code{union = FALSE} and \code{NULL}, an intercpet is added.
#' @param detr String vector indicating the type of detrending to be performed. Only relevant if \verb{union = FALSE}. Options are: \verb{"OLS"} and/or \verb{"QD"} (typically also called GLS). If NULL, set to \verb{"OLS"}.
#' @param ic_scale Logical indicator whether or not to use the rescaled information criteria (\verb{TRUE}) or not (\verb{FALSE}).
#' @param q Numeric vector of quantiles to be tested. Default is to test each unit sequentially.
#' @param h_rs Bandwidth used in rescaled information criteria.
#' @seealso \code{\link{iADFtest}}, \code{\link{BSQTtest}}, \code{\link{bFDRtest}}
#' @keywords internal
do_tests_and_bootstrap <- function(y, BSQT_test, iADF_test, level, boot, B, l, ar_AWB, union, p_min,
                                   p_max, ic, dc, detr, ic_scale, q, h_rs, show_progress,
                                   do_parallel, nc){
  y <- as.matrix(y)

  # Check correctness arguments and perform initial calculations and transformations
  inputs <- check_inputs(y = y, BSQT_test = BSQT_test, iADF_test = iADF_test, level = level, boot = boot,
               B = B, l = l, ar_AWB = ar_AWB, union = union, p_min = p_min, p_max = p_max, ic = ic,
               dc = dc, detr = detr, q = q, do_parallel = do_parallel, nc = nc)

  boot <- inputs$boot
  l <- inputs$l
  s_DWB <- inputs$s_DWB
  ar_AWB <- inputs$ar_AWB
  dc <- inputs$dc
  dc_boot <- inputs$dc_boot
  detr <- inputs$detr
  detr_int <- inputs$detr_int
  ic <- inputs$ic
  p_max <- inputs$p_max
  p_vec <- inputs$p_vec
  range_nonmiss <- inputs$range_nonmiss
  joint <- inputs$joint
  nc <- inputs$nc

  # Dimensions
  n <- nrow(y)
  N <- ncol(y)

  # ADF tests
  panel_est <- adf_panel_bootstrap_dgp_cpp(y = y, pmin = p_min, pmax = p_max, ic = ic,
                                           dc = dc_boot, QD = FALSE, trim = FALSE,
                                           ic_scale = ic_scale, h_rs = h_rs, range = range_nonmiss)
  u_boot <- panel_est$b_res
  u_boot[is.nan(u_boot)] <- NA
  res <- panel_est$res
  ar_est <- panel_est$par[-1, , drop = FALSE]
  t_star <- bootstrap_cpp(B = B, boot = boot, u = u_boot, e = res, l = l, s = s_DWB, ar = ar_AWB,
                          ar_est = ar_est, y0 = matrix(0, ncol = N), pmin = p_min, pmax = p_max,
                          ic = ic, dc = dc, detr = detr_int, ic_scale = ic_scale, h_rs = h_rs,
                          range = range_nonmiss, joint = joint, show_progress = show_progress,
                          do_parallel = do_parallel, nc = nc)
  tests_i <- adf_tests_panel_cpp(y, pmin = p_min, pmax = p_max, ic = ic, dc = dc, detr = detr_int,
                                  ic_scale = ic_scale, h_rs = h_rs, range = range_nonmiss)

  if (union) {
    scaling <- scaling_factors_cpp(t_star, level)
    if (N > 1) {
      test_stats_star <- union_tests_cpp(t_star, scaling)
      test_stats <- union_tests_cpp(array(tests_i,
                                          dim = c(1, length(dc) * length(detr_int), N)), scaling)
    } else {
      test_stats_star <- union_test_cpp(t_star[, , 1], scaling)
      test_stats <- union_test_cpp(array(tests_i,
                                         dim = c(1, length(dc) * length(detr_int))), scaling)
    }
   } else {
    test_stats_star <- NULL
    test_stats <- NULL
  }
  out <- list("y" = y, "p_vec" = p_vec, "t_star" = t_star, "test_stats_star" = test_stats_star,
              "tests_i" = tests_i, "test_stats" = test_stats, "level" = level, "dc" = dc,
              "detr" = detr)

  return(out)
}

#' Auxiliary Function (not accessible to users) to check if all arguments put in by the user are correct, and to perform some preliminary calculations.
#' @inheritParams iADFtest
#' @param BSQT_test Logical indicator whether or not to perform the Bootstrap Quantile Test (\verb{TRUE}) or not (\verb{FALSE}).
#' @param iADF_test Logical indicator whether or not to perform the individual ADF Tests (\verb{TRUE}) or not (\verb{FALSE}).
#' @param level Desired significance level of the unit root test.
#' @param boot String for bootstrap method to be used. Options are
#' \describe{
#' \item{\verb{"MBB"}}{Moving blocks bootstrap;}
#' \item{\verb{"BWB"}}{Block wild bootstrap;}
#' \item{\verb{"DWB"}}{Dependent wild bootstrap;}
#' \item{\verb{"AWB"}}{Autoregressive wild bootstrap;}
#' \item{\verb{"SB"}}{Sieve bootstrap;}
#' \item{\verb{"SWB"}}{Sieve wild boostrap.}
#' }
#' @param l Desired 'block length' in the bootstrap. For the MBB, BWB and DWB boostrap, this is a genuine block length. For the AWB boostrap, the block length is transformed into an autoregressive parameter via the formula \eqn{0.01^(1/l)}; this can be overwritten by setting \verb{ar_AWB} directly. If NULL, sets the block length as a function of the time series length T, via the rule \eqn{l = 1.75 T^(1/3)}.
#' @param ar_AWB Autoregressive parameter used in the AWB bootstrap method (\verb{boot = "AWB"}). Can be used to set the parameter directly rather than via the default link to the block length l.
#' @param union Logical indicator whether or not to use bootstrap union tests (\verb{TRUE}) or not (\verb{FALSE}).
#' @param p_min Minimum lag length in the augmented Dickey-Fuller regression.
#' @param p_max Maximum lag length in the augmented Dickey-Fuller regression.
#' @param ic String for information criterion used to select the lag length in the augmented Dickey-Fuller regression. Options are: \verb{"AIC"}, \verb{"BIC"}, \verb{"MAIC"}, \verb{"MBIC}.
#' @param dc Numeric vector indicating the deterministic specification. Only relevant if \code{union = FALSE}. Options are (combinations of): \code{0} - no deterministics; \code{1} - intercept only; \code{2} - intercept and trend. If \code{union = FALSE} and \code{NULL}, an intercpet is added.
#' @param detr String vector indicating the type of detrending to be performed. Only relevant if \verb{union = FALSE}. Options are: \verb{"OLS"} and/or \verb{"QD"} (typically also called GLS). If NULL, set to \verb{"OLS"}.
#' @param ic_scale Logical indicator whether or not to use the rescaled information criteria (\verb{TRUE}) or not (\verb{FALSE}).
#' @param q Numeric vector of quantiles to be tested. Default is to test each unit sequentially.
#' @param h_rs Bandwidth used in rescaled information criteria.
#' @seealso \code{\link{iADFtest}}, \code{\link{BSQTtest}}, \code{\link{bFDRtest}}
#' @keywords internal
check_inputs <- function(y, BSQT_test, iADF_test, level, boot, B, l, ar_AWB, union,
                         p_min, p_max, ic, dc, detr, q, do_parallel, nc){

  # Dimensions
  n <- nrow(y)
  N <- ncol(y)

  # Check if sufficient bootstrap replications are done
  if (level * (B + 1) < 1) {
    stop("Bootstrap iterations B too low to perform test at desired significance level.")
  }
  # Set up parallel computing
  if (is.null(nc)) {
    if (do_parallel) {
      nc <- parallel::detectCores() - 1
    } else {
      nc <- 1
    }
  }
  if ((nc != round(nc)) | (nc < 1)) {
    stop("Invalid value for argument nc")
  }

  # Check for missing values or unbalanced panels (MBB, SB)
  check_missing <- check_missing_insample_values(y)
  if (any(check_missing)) {
    stop("Missing values detected inside sample.")
  } else {
    joint <- TRUE
    check_nonmiss <- find_nonmissing_subsample(y)
    range_nonmiss <- check_nonmiss$range - 1
    all_range_equal <- check_nonmiss$all_equal
    if (!all_range_equal) {
      if (boot %in% c("MBB", "SB")) {
        if (!iADF_test) {
          stop("Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
        } else if (N > 1) {
          joint <- FALSE
          warning("Missing values cause resampling bootstrap to be executed for each time series individually.")
        }
      }
    }
  }

  # Checks on inputs
  if (!(boot %in% c("MBB", "BWB", "DWB", "AWB", "SB", "SWB")) | length(boot) > 1) {
    stop("The argument boot should be equal to either MBB, BWB, DWB, AWB, SB or SWB")
  } else if (boot %in% c("SB", "SWB") & !iADF_test) {
    warning("SB and SWB bootstrap only recommended for iADFtest; see help for details.")
  }
  boot <- 1 * (boot == "MBB") + 2*(boot == "BWB") + 3 * (boot == "DWB") +
    4 * (boot == "AWB") + 5 * (boot == "SB") + 6 * (boot == "SWB")

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
    dc_boot <- 2
    detr_int <- 1:2
  } else {
    if (is.null(dc)) {
      stop("No deterministic specification set.
           Set dc to 0, 1 or 2, or set union to TRUE for union test.")
    } else if (any(!is.element(dc, 0:2))){
      stop("The argument dc should only contain values 0, 1, and/or 2:
           (0: no deterministics, 1: intercept, 2: intercept and trend)")
    }
    dc <- sort(dc)
    dc_boot <- max(dc)
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
    l <- round(1.75 * nrow(y)^(1/3))
  }

  if (is.null(ar_AWB)) {
    ar_AWB <- 0.01^(1/l)
  } else if (boot != 4){
    warning("Argument ar_AWB set, but AWB method not used. ar_AWB is ignored.")
  }

  # Defaults
  p_vec <- NULL
  if (BSQT_test) {
    if (is.numeric(q) & all(q == floor(q) & q >= 0) & !anyNA(q)) {
      p_vec <- sort(unique(q[q <= N]))
      if (max(p_vec) < N) {
        p_vec <- c(p_vec, N)
      }
      if (min(p_vec) > 0) {
        p_vec <- c(0, p_vec)
      }
      if (!identical(p_vec, q)) {
        warning(paste0(paste0("Input to argument q transformed to fit sequential test: q = c("),
                       paste0(p_vec, collapse = ", "), ")."))
      }
      if (!identical(unique(p_vec), p_vec)) {
        p_vec <- unique(p_vec)
        warning(paste0(paste0("Input to argument q transformed to remove duplicate groups: q = c("),
                       paste0(p_vec, collapse = ", "), ")."))
      }
    } else if (is.numeric(q) & all(q >= 0 & q <= 1) & !anyNA(q)) {
      q_vec <- sort(unique(q))
      if (max(q_vec) < 1) {
        q_vec <- c(q_vec, 1)
      }
      if (min(q_vec) > 0) {
        q_vec <- c(0, q_vec)
      }
      if (!identical(q_vec, q)) {
        warning(paste0(paste0("Input to argument q transformed to fit sequential test: q = c("),
                       paste0(q_vec, collapse = ", "), ")"))
      }

      p_vec <- round(q_vec * N)
      if (!identical(unique(p_vec), p_vec)) {
        p_vec <- unique(p_vec)
        warning(paste0(paste0("Input to argument q transformed to remove duplicate groups after transformation to integers: q = c("),
                       paste0(p_vec, collapse = ", "), ")."))
      }
    } else {
      stop("Invalid input values for q: must be quantiles or positive integers.")
    }
  }
  if (boot == 1){
    # Obtain the probability that we draw identical blocks for the whole time series:
    # This will cause multicollinearity if p_max is larger than l.
    prob_identical_bl <- 1 - (1 - (1 / (n - l))^(ceiling(n / l) - 1))^B
  }
  if(is.null(p_max)){
    # Correction for small samples as formula doesn't work well for micropanels
    p_max = round(12*(n/100)^(1/4)) - 7*max(1 - n/50, 0)*(n/100)^(1/4)
    if (boot == 1) {
      if (prob_identical_bl > 0.01) {
        # If the probability of obtaining a multicollinear bootstrap sample is too large,
        # force p_max to be smaller than l.
        p_max <- min(p_max, l)
      }
    }
  }
  s_DWB <- matrix(0, n, n)
  if (boot == 3) {
    # Self-convolution
    self.conv <- function(t, c) {
      y <- apply(as.matrix(t), c(1,2), FUN = function(x){
        stats::integrate(f.w, -1, 1, t = x, c = c)$value})
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
    m <- self.conv(outer(1:n, 1:n, "-") / l, 0.43)/ drop(self.conv(0, 0.43))
    ev <- eigen(m)
    e_va <- ev$values
    e_ve <- ev$vectors
    e_va <- diag((e_va > 1e-10) * e_va + (e_va <= 1e-10) * 1e-10)
    s_DWB <- e_ve %*% sqrt(e_va)
  }
  out <- list(boot = boot, l = l, s_DWB = s_DWB, ar_AWB = ar_AWB, dc = dc,
              dc_boot = dc_boot, detr = detr, detr_int = detr_int, ic = ic,
              p_max = p_max, p_vec = p_vec, range_nonmiss = range_nonmiss, joint = joint, nc = nc)

  return(out)
}
