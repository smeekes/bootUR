#' Differences of Multiple Time Series
#' @description Performs differencing of multiple time series, with possibly different orders for each time series.
#' @param data A (\eqn{T}x\eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}) or data frame.
#' @param d An \eqn{N}-dimensional vector containing the orders
#' @param keep_NAs Logical indicator whether or not to keep the \code{NA} values resulting from differencing at the beginning of the sample. Default is \code{TRUE}. If \code{FALSE}, the entire row containing the \code{NA} values is removed.
#' @export
#' @return The appropriately differenced data in the same format as the original data.
diff_mult <- function(data, d, keep_NAs = TRUE) {
  x <- as.matrix(data)
  diffed_x <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  if (length(c(d)) != ncol(x)) {
    stop("Argument d should have length equal to columns of data.")
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
  diffed_data <- data
  diffed_data[] <- diffed_x
  if (!keep_NAs & (max(d) > 0)) {
    if (NCOL(diffed_data) > 1) {
      diffed_data <- diffed_data[-(1:max(d)), ]
    } else if (NCOL(diffed_data) == 1) {
      diffed_data <- diffed_data[-(1:max(d))]
    }
  }
  return(diffed_data)
}

#' Determine Order of Integration
#' @description Determines the order of integration for each time series in a dataset via a sequence of unit root tests, and differences the data accordingly to eliminate stochastic trends.
#' @param data A (\eqn{T}x\eqn{N})-matrix of \eqn{N} time series with \eqn{T} observations. Data may also be in a time series format (e.g. \code{ts}, \code{zoo} or \code{xts}) or data frame.
#' @param max_order The maximum order of integration of the time series. Default is 2.
#' @param method The unit root tests to be used in the procedure. For multiple time series the options are "boot_ur", "boot_sqt" and "boot_fdr", with "boot_ur" the default. For single time series the options are "adf", boot_adf", "boot_union" and "boot_ur", with the latter the default.
#' @param level Desired significance level of the unit root test. Default is 0.05.
#' @param plot_orders Logical indicator whether the resulting orders of integration should be plotted. Default is \code{FALSE}.
#' @param data_name Optional name for the data, to be used in the output. The default uses the name of the 'data' argument.
#' @param ... Optional arguments passed to the chosen unit root test function.
#' @details The function follows the approach laid out in Smeekes and Wijler (2020), where all series is differenced \eqn{d-1} times, where \eqn{d} is the specified maximum order, and these differenced series are tested for unit roots. The series for which the unit root null is not rejected, are classified as \eqn{I(d)} and removed from consideration. The remaining series are integrated, and tested for unit roots again, leading to a classification of \eqn{I(d-1)} series if the null is not rejected. This is continued until a non-rejection is observed for all time series, or the series are integrated back to their original level. The series for which the null hypothesis is rejected in the final stage are classified as \eqn{I(0)}.
#'
#' Care must be taken when using \code{\link{boot_sqt}} when the argument \code{steps} is given as a sequence of integers. As at each step series are removed, one may end up with fewer series to test than indicated in \code{steps}. While integers larger than the number of series will automatically be removed - along with a warning - by the test, it is recommend to set \code{steps} in the form of quantiles.
#'
#' Plotting the orders of integration requires the \code{ggplot2} package to be installed; plot will be skipped and a warning is given if not. For plots the function \code{\link{plot_order_integration}} is called. The user may prefer to set \code{plot_orders = FALSE} and call this function directly using the returned value of \code{order_int} in order to have more control over plot settings and save the plot object.
#' @export
#' @return A list with the following components
#' \item{\code{order_int}}{A vector with the found orders of integration of each time series.}
#' \item{\code{diff_data}}{The appropriately differenced data according to \code{order_int} in the same format as the original data.}
#' @references Smeekes, S. and Wijler, E. (2020). Unit roots and cointegration. In P. Fuleky (Ed.) \emph{Macroeconomic Forecasting in the Era of Big Data}, Chapter 17, pp. 541-584. \emph{Advanced Studies in Theoretical and Applied Econometrics}, vol. 52. Springer.
#' @examples
#' # Use "boot_ur" to determine the order of GDP_BE and GDP_DE
#' orders_tseries <- order_integration(MacroTS[, 1:2], method = "boot_ur", B = 199)
order_integration <- function(data, max_order = 2, method = "boot_ur", level = 0.05,
                              plot_orders = FALSE, data_name = NULL, ...) {
  N <- NCOL(data)
  if (!is.null(colnames(data))) {
    var_names <- colnames(data)
  } else {
    var_names <- paste0("Variable ", 1:NCOL(data))
  }
  d <- rep(NA, N)
  names(d) <- var_names
  data_mat <- as.matrix(data)
  datad <- data_mat
  i_in_datad <- 1:N
  for (d_i in (max_order - 1):0) {
    datad <- diff_mult(data_mat[, i_in_datad], rep(d_i, length(i_in_datad)), keep_NAs = FALSE)
    if (method == "boot_ur") {
      out <- boot_ur(datad, level = level, ...)
    } else if (method == "boot_fdr" & N > 1) {
      out <- boot_fdr(datad, level = level, ...)
    } else if (method == "boot_sqt" & N > 1) {
      out <- boot_sqt(datad, level = level, ...)
    } else if (method == "boot_adf" & N == 1) {
      test_out <- boot_adf(datad, level = level, ...)
      out <- list("rejections" = test_out$p.value < level)
    } else if (method == "boot_union" & N == 1) {
      test_out <- boot_union(datad, level = level, ...)
      out <- list("rejections" = test_out$p.value < level)
    } else if (method == "adf" & N == 1) {
      test_out <- adf(datad, ...)
      out <- list("rejections" = test_out$p.value < level)
    } else if (method == "adf" & N > 1) {
        rejections <- apply(datad, 2, function(x){
          return(adf(x, ...)$p.value < level)
        })
        out <- list("rejections" = rejections)
      } else if (method == "iADFtest") {
      stop("'iADFtest' is deprecated. Use 'boot_ur' instead.")
#      out <- iADFtest(datad, ...)
    } else if (method == "bFDRtest" & N > 1) {
      stop("'bFDRtest' is deprecated. Use 'boot_fdr' instead.")
#      out <- bFDRtest(datad, ...)
    } else if (method == "BSQTtest" & N > 1) {
      stop("'BSQTtest' is deprecated. Use 'boot_sqt' instead.")
#      out <- BSQTtest(datad, ...)
    } else if (method == "boot_df" & N == 1) {
      stop("'boot_df' is deprecated. Use 'boot_adf' instead.")
#      out <- boot_df(datad, ...)
    } else {
      stop("Invalid 'method' argument.")
    }
    d[i_in_datad[!out$rejections]] <- d_i + 1
    if (any(out$rejections)) {
      if (d_i == 0) {
        d[i_in_datad[out$rejections]] <- 0
      } else {
        i_in_datad <- i_in_datad[out$rejections]
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
  return(list(order_int = d, diff_data = diff_mult(data, d)))
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
#' @param cols Vector with colours for displaying the different types of data. If the default is overwritten, four colours must be supplied.
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
