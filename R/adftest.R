#' Augmented Dickey-Fuller Unit Root Test
#' @description This function performs a standard augmented Dickey-Fuller unit root test on a single time series.
#' @inheritParams boot_ur
#' @details The function encompasses the standard augmented Dickey-Fuller test. The reported p-values are MacKinnon's unit root p-values taken from the package urca.
#'
#' Lag length selection is done automatically in the ADF regression with the specified information criterion. If one of the modified criteria of Ng and Perron (2001) is used, the correction of Perron and Qu (2008) is applied. For very short time series (fewer than 50 time points) the maximum lag length is adjusted downward to avoid potential multicollinearity issues in the bootstrap. To overwrite data-driven lag length selection with a pre-specified lag length, simply set both the minimum `min_lag` and maximum lag length `max_lag` for the selection algorithm equal to the desired lag length.
#' @export
#' @return Values of the Dickey-Fuller test statistics and corresponding bootstrap p-values.
#' @section Errors and warnings:
#' \describe{
#' \item{\code{Error: Multiple time series not allowed. Switch to a multivariate method such as boot_ur, or change argument data to a univariate time series.}}{The function provides a standard ADF test with asymptotic p-value. It does not support multiple time series}
#' }
#' @examples
#' # standard ADF test on GDP_BE
#' GDP_BE_adf <- adf(MacroTS[, 1], deterministics = "trend")
adf <- function(data, min_lag = 0, max_lag = NULL, criterion = "MAIC", deterministics = "intercept", criterion_scale = TRUE){
  #### !!! For now only with detr = "OLS" & 2-step detrending !!! ####
  
  if (NCOL(data) > 1) {
    stop("Multiple time series not allowed. Switch to a multivariate method such as boot_ur,
         or change argument data to a univariate time series.")
  }
  
  # Checking inputs 
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
  
  if (any(!is.element(criterion, c("AIC", "BIC", "MAIC", "MBIC"))) | length(criterion) > 1) {
    stop("The argument criterion should be equal to either AIC, BIC, MAIC, MBIC)")
  }
  ic <- 1*(criterion=="AIC") + 2*(criterion=="BIC") + 3*(criterion=="MAIC") + 4*(criterion=="MBIC")
  
  if (is.null(deterministics)) {
    stop("No deterministic specification set.
         Set deterministics to the strings none, intercept or trend.")
  } else if (any(!is.element(deterministics, c("none", "intercept", "trend")))){
    stop("The argument deterministics should only contain the strings none, intercept and/or trend:
           (none: no deterministics, intercept: intercept only, trend: intercept and trend)")
  }
  dc_int <- 0*(deterministics=="none") + 1*(deterministics=="intercept") + 2*(deterministics=="trend")
  dc_int <- sort(dc_int)
  
  # if (is.null(detr)) {
  #   warning("No detrending specification set. Using OLS detrending.")
  #   detr <- "OLS"
  # } else if(any(!is.element(detr, c("OLS", "QD")))) {
  #   stop("The argument detr should only contain the strings OLS and/or QD")
  # }
  detr <- "OLS" # ONLY OLS FOR NOW
  detr_int <- 1*(detr=="OLS") + 2*(detr=="QD")
  detr_int <- sort(detr_int)
  
  if(is.null(max_lag)){
    # Correction for small samples as formula doesn't work well for micropanels
    max_lag = round(12*(n/100)^(1/4)) - 7*max(1 - n/50, 0)*(n/100)^(1/4)
  }
  
  # Get ADF test statistic
  teststats <- adf_tests_panel_cpp(y, pmin = min_lag, pmax = max_lag, ic = ic, dc = dc_int, detr = detr_int,
                                   ic_scale = criterion_scale, h_rs = 0.1, range = range_nonmiss)
  
  # Get p-values from urca package 
  pvalues <- rep(NA, length(deterministics))
  # If deterministics is not a vector, we don't need for-loop and we can use this 
  # switch(deterministics,
  #        "trend" = urtype <- "ct",
  #        "intercept" = urtype <- "c",
  #        "none"  =  urtype <- "nc")
  for(idc in 1:length(dc_int)){
    if(dc_int[idc]==0){
      urtype <- "nc"
    }else if (dc_int[idc]==1){
      urtype <- "c"
    }else{
      urtype <- "ct"
    }
    pvalues[idc] <- urca::punitroot(q = teststats[idc], N = TT - max_lag - 1, trend = urtype, statistic = "t")
  }
  
  
  # Display results for now in matrix form
  results <- cbind(teststats, pvalues)
  colnames(results) <- c("test statistic", "p-value")
  rownames(results)[dc_int==0] <- "none"
  rownames(results)[dc_int==1] <- "intercept"
  rownames(results)[dc_int==2] <- "trend"
  
  return(results)
}
