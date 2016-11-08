# finTools.R

#' Auto select ARMA+GARCH
#'
#' automatically delect orders for ARMA, and distributions for innovations of GARCH
#' @param dat: input of univariate time series data.
#' @param distvec: vector of names of distributions to choose from for innovations
#' @param maxAR: maximum AR order
#' @param maxMA: maximum MA order
#' @param pqGARCH: specify GARCH orders (default: c(1,1))
#' @param windowLength: rolling window length
#' @param sinkflag: whether sink outputs (default: T)
#' @keywords auto ARMA GARCH
#' @export
#' @examples
#' data("USDCAD2015_H1")
#' dat <- USDCAD2015_H1$Close[1:250] 
#' logr <- diff(log(dat))
#' autoARMAGARCH(dat=logr, windowLength=245)

autoARMAGARCH <- function(dat, distvec = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu"), 
                          maxAR = 4, maxMA = 4, pqGARCH = c(1,1), windowLength = 100, sinkflag = T)
{
  library(rugarch)
  if(sinkflag == T){sink(paste("autoARMA_GARCH_output_", format(as.numeric(Sys.time())), ".txt", sep=""), append=T, split = T)}
  
  dat <- ts(dat)

  foreLength <- length(dat) - windowLength
  output <- data.frame(i=integer(foreLength+1),
                       ar=integer(foreLength+1),
                       ma=integer(foreLength+1),
                       p=integer(foreLength+1),
                       q=integer(foreLength+1),
                       dist=character(foreLength+1),
                       aic=numeric(foreLength+1),
                       stringsAsFactors = FALSE)
  
  for (d in 0:foreLength)
  {
    dataOffset <- ts(dat[(1+d):(windowLength+d)])
    
    # Fit the ARIMA model
    final.aic <- Inf
    final.order <- c(0,0,0)
    for (p in 0:maxAR) for (q in 0:maxMA) {
      if ( p == 0 && q == 0) {
        next
      }
      
      arimaFit <- tryCatch( arima(dataOffset, order=c(p, 0, q)),
                           error=function( err ) FALSE,
                           warning=function( err ) FALSE )
      
      if( !is.logical( arimaFit ) ) {
        current.aic <- AIC(arimaFit)
        if (current.aic < final.aic) {
          final.aic <- current.aic
          final.order <- c(p, 0, q)
          final.arima <- arima(dataOffset, order=final.order)
        }
      } else {
        next
      }
    }
    
    # Specify and fit the GARCH model
    final.aic <- Inf
    current.aic <- Inf
    final.dist <- "NA"
    
    # Fit the ARIMA model
    for(dist in distvec)
    {
      spec <- rugarch::ugarchspec(
        variance.model=list(garchOrder=pqGARCH),
        mean.model=list(armaOrder=c(final.order[1], final.order[3]), include.mean=T),
        distribution.model=dist
      )
      
      fit = tryCatch(
        rugarch::ugarchfit(
          spec, dataOffset, solver = 'hybrid'
        ), error=function(err) FALSE, warning=function(err) FALSE
      )
      
      if( !is.logical( fit ) ) {
        current.aic <- rugarch::infocriteria(fit)[1]
        if (current.aic < final.aic) {
          final.aic <- current.aic
          final.dist <- dist
          final.fit <- fit
        }
      } else {
        next
      }
    }
    output$i[d+1] <- d+1
    output$ar[d+1] <- final.order[1]
    output$ma[d+1] <- final.order[3]
    output$p[d+1] <- pqGARCH[1]
    output$q[d+1] <- pqGARCH[2]
    output$dist[d+1] <- final.dist
    output$aic[d+1] <- final.aic
    if(sinkflag == T){cat(d, final.order[1], final.order[3], pqGARCH, final.dist, final.aic, "\n")}
  }
  if(sinkflag == T){sink()} # close sink()
  output
}
