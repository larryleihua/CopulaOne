# forex.R

data("USDCAD2015_H1")
data("AUDUSD2015_H1")

USDCAD2015_H1$Time.format <- chron::chron(dates=as.character(USDCAD2015_H1$Date),times=USDCAD2015_H1$Timestamp,format=c('ymd','h:m:s'))
AUDUSD2015_H1$Time.format <- chron::chron(dates=as.character(AUDUSD2015_H1$Date),times=AUDUSD2015_H1$Timestamp,format=c('ymd','h:m:s'))

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

# rst <- autoARMAGARCH(logr[,1], windowLength=windowLength)
