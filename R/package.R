#' Full-range tail dependence copulas
#'
#' Distribution functions, simulation, dependence measures, and likelihood
#' fitting for the FRA1, GGEE, and PPPP bivariate copula families.
#' @importFrom stats cor dnorm pnorm rgamma rexp runif pbeta qbeta plogis qlogis
#'   optim optimHess nlm uniroot integrate setNames rnorm
#' @importFrom graphics contour plot
#' @importFrom methods is
#' @keywords internal
"_PACKAGE"

#' European stock-index example data
#'
#' Observations for European stock indices supplied with CopulaOne for fitting examples.
#' @format A data frame with 957 rows and 6 columns: date, cac40, dax, ftse,
#'   oseax, and smi. The index columns contain the supplied return observations.
#' @usage data(euro0306)
#' @name euro0306
#' @docType data
#' @keywords datasets
NULL

#' Hourly AUD/USD example data
#'
#' AUD/USD observations supplied with CopulaOne.
#' @format A data frame with 6226 rows and 7 columns: Date, Timestamp, Open,
#'   High, Low, Close, and Volume.
#' @usage data(AUDUSD2015_H1)
#' @name AUDUSD2015_H1
#' @docType data
#' @keywords datasets
NULL

#' Hourly USD/CAD example data
#'
#' USD/CAD observations supplied with CopulaOne.
#' @format A data frame with 6226 rows and 7 columns: Date, Timestamp, Open,
#'   High, Low, Close, and Volume.
#' @usage data(USDCAD2015_H1)
#' @name USDCAD2015_H1
#' @docType data
#' @keywords datasets
NULL
