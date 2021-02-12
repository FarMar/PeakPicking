#####################################################################################################

#### Thermalox peak-picking script                                                ###################

#### mark.farrell@csiro.au        +61 8 8303 8664         12/02/2021 ################################

#####################################################################################################

#### Set working directory ####
setwd("/Users/markfarrell/OneDrive - CSIRO/Data/PeakPicking")


#### Install packages as needed ####
install.packages("zoo")
install.packages("DescTools")


#### Load packages ####

library(tidyverse)
library(janitor)
library(zoo)
library(pracma)
library(DescTools)

#### import data ####
n_raw <- read_csv("data/tn.csv")
plot(n_raw, col = 'Gray', type = 'l')


#### smooth ####
n_raw$y.smooth <- loess(n_raw$signal ~ n_raw$time, span=0.005)$fitted
plot(n_raw$time, n_raw$y.smooth, type = 'l')



#x <- n_raw$time
#y <- n_raw$signal


#### old smooth ####
# https://rpubs.com/mengxu/peak_detection?fbclid=IwAR0t8Ni1SvEi86m5yeZL9aQhPE_FxRnR6kFMYVsKX6zQXt3xHVG66g9iPkE
#local regression
n_raw$y.smooth <- loess(n_raw$signal ~ n_raw$time, span=0.008)$fitted
plot(n_raw$time, n_raw$y.smooth, type = 'l')

#resampling
w=30
y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
x.max <- rollapply(zoo(x), 2*w+1, median, align="center")
length(y.max)
length(x.max)
plot(x.max, y.max, col = 'Gray', type='l')
lines(x.max, y.max, col = 'SkyBlue', lwd = 2)

plot(x, y, col = 'Gray', type = 'l')
lines(x, y.smooth, col = 'Black')
lines(x.max, y.max, col = 'SkyBlue', lwd = 2)

legend('topleft', c('1', '2', '3'), cex=0.8, col=c('Gray', 'Black', 'SkyBlue'), lty=c(1,1,1))

#peak detection
n <- length(y)
delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
plot(x.max, delta, type='l')
abline(h = 0, lty='dotted', col = 'red')

argmax <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

test <- function(w, span) {
  peaks <- argmax(x, y, w=w, span=span)
  
  plot(x, y, cex=0.75, col="Gray", main=paste("w = ", w, ", span = ", span, sep=""))
  lines(x, peaks$y.hat,  lwd=2) #$
  y.min <- min(y)
  sapply(peaks$i, function(i) lines(c(x[i],x[i]), c(y.min, peaks$y.hat[i]), col="Red", lty=2))
  points(x[peaks$i], peaks$y.hat[peaks$i], col="Red", pch=19, cex=1.25)
}

test(100, 0.008)

#so the above identifies peaks well, but not obvious how to output them or integrate their areas. Onwards...

#### go pracma ####
#Uses {pracma::findpeaks} to identify peaks

peaks <- data.frame(findpeaks(n_raw$y.smooth, 
                              npeaks=23, 
                              threshold=0, 
                              peakpat = "[+]{1,}[0]*[-]{1,}", 
                              sortstr=TRUE)) #IDs peaks and makes a df that has
#                                            # X1 = height at point of picking
#                                            # X2 = time at point of picking
#                                            # X3 = time at start of peak
#                                            # X4 = time at end of peak

n_raw$n <- seq(1,length(n_raw$y.smooth))
merged <- merge(x=n_raw, 
               y=peaks, 
               by.x="n", 
               by.y="X2", 
               all.x=TRUE, 
               all.y=TRUE)

ggplot(merged, aes(x=time, y=signal)) +
  geom_col(orientation="x", colour = "grey") +
  geom_line(aes(x=time, y=y.smooth))+
  geom_point(aes(x=time, y=X1), colour = "red")

trim <- merged %>% slice(4581:4720)

ggplot(trim, aes(x=time, y=signal)) +
  geom_col(orientation="x", colour = "grey", fill = "grey") +
  geom_line(aes(x=time, y=y.smooth))+
  geom_point(aes(x=time, y=X1), colour = "red", size = 6)

AUC(trim$n, trim$signal, method = "spline")
