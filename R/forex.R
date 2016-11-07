# forex.R

data("USDCAD2015_H1")
data("AUDUSD2015_H1")

USDCAD2015_H1$Time.format <- chron(dates=as.character(USDCAD2015_H1$Date),times=USDCAD2015_H1$Timestamp,format=c('ymd','h:m:s'))
AUDUSD2015_H1$Time.format <- chron(dates=as.character(AUDUSD2015_H1$Date),times=AUDUSD2015_H1$Timestamp,format=c('ymd','h:m:s'))

alldat <- merge(USDCAD2015_H1, AUDUSD2015_H1, by="Time.format", all.x=FALSE)
datused <- data.frame(alldat$Time.format, 1/alldat$Close.x, alldat$Close.y)
names(datused) <- c("time", "CADUSD.cl", "AUDUSD.cl")
rm(alldat,USDCAD2015_H1,AUDUSD2015_H1)
logr <- apply(datused[,2:3], 2, function(x){diff(log(x))})

# 60-min data
Nmin <- 60
Ndayroll <- 10 #how many days for one rolling window
windowLength <- (60/Nmin)*24*Ndayroll
foreLength <- dim(logr)[1] - windowLength
forecasts <- vector(mode="character", length=foreLength)

#ex <- "CADUSD"
ex <- "AUDUSD"
i_ex <- (ex == "CADUSD")*1 + (ex == "AUDUSD")*2
sink(paste("fitting_", ex, "_", format(Nmin), ".txt", sep=""), append=T, split = T)
for (d in 0:foreLength)
{
  # Obtain the S&P500 rolling window for this day
  logrOffset = ts(logr[(1+d):(windowLength+d),i_ex])   # 1: CADUSD; 2: AUDUSD
  
  # Fit the ARIMA model
  final.aic <- Inf
  final.order <- c(0,0,0)
  for (p in 0:4) for (q in 0:4) {
    if ( p == 0 && q == 0) {
      next
    }
    
    arimaFit = tryCatch( arima(logrOffset, order=c(p, 0, q)),
                         error=function( err ) FALSE,
                         warning=function( err ) FALSE )
    
    if( !is.logical( arimaFit ) ) {
      current.aic <- AIC(arimaFit)
      if (current.aic < final.aic) {
        final.aic <- current.aic
        final.order <- c(p, 0, q)
        final.arima <- arima(logrOffset, order=final.order)
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
  for(dist in c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu"))
  {
    spec = ugarchspec(
      variance.model=list(garchOrder=c(1,1)),
      mean.model=list(armaOrder=c(final.order[1], final.order[3]), include.mean=T),
      distribution.model=dist
    )
    
    fit = tryCatch(
      ugarchfit(
        spec, logrOffset, solver = 'hybrid'
      ), error=function(err) FALSE, warning=function(err) FALSE
    )
    
    if( !is.logical( fit ) ) {
      # cat(dist, infocriteria(fit)[1], "\n")
      current.aic <- infocriteria(fit)[1]
      if (current.aic < final.aic) {
        final.aic <- current.aic
        final.dist <- dist
        final.fit <- fit
      }
    } else {
      next
    }
  }
  cat(d, " final.order = ", final.order, " final.dist = ", final.dist, " final.aic = ", final.aic, "\n")
}
sink()
