#' Augmented Dickey-Fuller Unit Root Test
#' @description This function performs a standard augmented Dickey-Fuller unit root test on a single time series.
#' @inheritParams boot_adf
#' @param two_step Logical indicator whether to use one-step (\code{two_step = FALSE}) or two-step (\code{two_step = TRUE}) detrending. The default is two-step detrending.
#' @details The function encompasses the standard augmented Dickey-Fuller test. The reported p-values are MacKinnon's unit root p-values taken from the package urca.
#'
#' Lag length selection is done automatically in the ADF regression with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. For very short time series (fewer than 50 time points) the maximum lag length is adjusted downward to avoid potential multicollinearity issues in the bootstrap. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `min_lag` and maximum lag length `max_lag` for the selection algorithm equal to the desired lag length.
#' @export
#' @return An object of class \code{"bootUR"}, \code{"htest"} with the following components:
#' \item{\code{method}}{The name of the hypothesis test method;}
#' \item{\code{data.name}}{The name of the variable on which the method is performed;}
#' \item{\code{null.value}}{The value of the (gamma) parameter of the lagged dependent variable in the ADF regression under the null hypothesis. Under the null, the series has a unit root. Testing the null of a unit root then boils down to testing the significance of the gamma parameter;}
#' \item{\code{alternative}}{A character string specifying the direction of the alternative hypothesis relative to the null value. The alternative postulates that the series is stationary;}
#' \item{\code{estimate}}{The estimated value of the (gamma) parameter of the lagged dependent variable in the ADF regression;}
#' \item{\code{statistic}}{The value of the test statistic of the ADF unit root test;}
#' \item{\code{p.value}}{The p-value of the ADF unit root test.}
#' \item{\code{specifications}}{The specifications used in the test.}
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Multiple time series not allowed. Switch to a multivariate method such as boot_ur, or change argument data to a univariate time series.}}{The function provides a standard ADF test with asymptotic p-value. It does not support multiple time series}
#' }
#' @examples
#' # standard ADF test on GDP_BE
#' GDP_BE_adf <- adf(MacroTS[, 1], deterministics = "trend")
adf <- function(data, data_name = NULL, deterministics = "intercept", min_lag = 0,
                max_lag = NULL, criterion = "MAIC",
                criterion_scale = TRUE, two_step = TRUE){

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

  #### Checking inputs ####
  y <- as.matrix(data)
  n <- nrow(y)

  check_missing <- check_missing_insample_values(y)
  if (any(check_missing)) {
    stop("Missing values detected inside sample.")
  } else {
    joint <- TRUE
    check_nonmiss <- find_nonmissing_subsample(y)
    range_nonmiss <- check_nonmiss$range - 1
    TT <- range_nonmiss[2] - range_nonmiss[1] + 1
  }


  if(is.null(max_lag)){
    # Correction for small samples as formula doesn't work well for micropanels
    max_lag = round(12*(n/100)^(1/4) - 7*max(1 - n/50, 0)*(n/100)^(1/4))
  }

  if (any(!is.element(criterion, c("AIC", "BIC", "MAIC", "MBIC"))) | length(criterion) > 1) {
    stop("The argument criterion should be equal to either AIC, BIC, MAIC, MBIC)")
  }
  ic <- 1*(criterion == "AIC") + 2*(criterion == "BIC") + 3*(criterion == "MAIC")
        + 4*(criterion == "MBIC")

  if (is.null(deterministics)) {
    stop("No deterministic specification set.
         Set deterministics to the strings none, intercept or trend.")
  } else if (any(!is.element(deterministics, c("none", "intercept", "trend"))) |
             length(deterministics) > 1){
    stop("The argument deterministics should be equal to either none, intercept, trend:
           (none: no deterministics, intercept: intercept only, trend: intercept and trend)")
  }
  dc_int <- 0*(deterministics=="none") + 1*(deterministics=="intercept") +
    2*(deterministics=="trend")

  detr <- "OLS"
  detr_int <- 1

  # Get ADF test statistic: two-step detrending
  if (two_step) {
    tests_and_params <- adf_tests_panel_cpp(y, pmin = min_lag, pmax = max_lag, ic = ic,
                                            dc = dc_int, detr = detr_int,
                                            ic_scale = criterion_scale, h_rs = 0.1,
                                            range = range_nonmiss)
  } else {
    tests_and_params <- adf_onestep_tests_panel_cpp(y, pmin = min_lag, pmax = max_lag,
                                                    ic = ic, dc = dc_int,
                                                    ic_scale = criterion_scale,
                                                    h_rs = 0.1, range = range_nonmiss)
  }

  # Collect estimate, test statistic and pvalue ADF test
  tstat <- drop(tests_and_params$tests) # Test statistics
  attr(tstat, "names") <- "tstat"
  param <- c(tests_and_params$par) # Parameter estimates
  attr(param, "names") <- "gamma"

  switch(deterministics,
         "trend" = urtype <- "ct",
         "intercept" = urtype <- "c",
         "none"  =  urtype <- "nc")
  p_val <- drop(urca::punitroot(q = tstat, N = max(TT - max_lag - 1, 20),
                                trend = urtype, statistic = "t"))
  attr(p_val, "names") <- "p-value"

  details <- list("individual estimates" = param,
                  "individual statistics" = tstat,
                  "individual p-values" = p_val,
                  "selected lags" = drop(tests_and_params$lags),
                  "txt_null" = "Series has a unit root",
                  "txt_alternative" = "Series is stationary")

  spec <- list("deterministics" = deterministics, "min_lag" = min_lag, "max_lag" = max_lag,
               "criterion" = criterion, "criterion_scale" = criterion_scale, "two_step" = two_step)
  switch(deterministics,
         "trend" = deterministics <- "intercept and trend",
         "intercept" = deterministics <- "intercept",
         "none"  =  deterministics <- "no deterministics")
  if(two_step){
    method_name <- paste0("Two-step ADF test (with " , deterministics,") on a single time series")
  }else{
    method_name <- paste0("One-step ADF test (with " , deterministics,") on a single time series")
  }

  adf_out <- list(method = method_name, data.name = data_name,
                  null.value = c("gamma" = 0), alternative = "less",
                  estimate = drop(param), statistic = tstat, p.value = p_val,
                  details = details, series.names = data_name, specifications = spec)
  class(adf_out) <- c("bootUR", "htest")

  return(adf_out)
}
