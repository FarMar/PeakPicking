#####################################################################################################

#### Thermalox peak-picking script                                                ###################

#### mark.farrell@csiro.au        +61 8 8303 8664         12/02/2021 ################################

#### senani.karunaratne@csiro.au                                     ################################

#####################################################################################################

#### Set working directory ####
setwd("/Users/markfarrell/OneDrive - CSIRO/Data/PeakPicking")


#### Install packages as needed ####
install.packages("zoo")
install.packages("DescTools")
install.packages("ggpmisc")


#### Load packages ####

library(tidyverse)
library(janitor)
library(zoo)
library(pracma)
library(DescTools)
library(ggpmisc)

#### TN ####
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

# sort peaks df
peaks_sort <- peaks %>% 
  arrange(X2)

#generate peak number vector and combine with sorted peaks df
peak_no <- 1:23
peaks_sort <- cbind(peak_no, peaks_sort)

merged <- merge(x=n_raw, 
               y=peaks_sort, 
               by.x="n", 
               by.y="X2", 
               all.x=TRUE, 
               all.y=TRUE)

ggplot(merged, aes(x=time, y=signal)) +
  geom_col(orientation="x", colour = "grey") +
  geom_line(aes(x=time, y=y.smooth))+
  geom_point(aes(x=time, y=X1), colour = "red")


#### individual curve area calculation ####
trim <- merged %>% slice(4581:4720)

ggplot(trima, aes(x=time, y=signal)) +
  geom_col(orientation="x", colour = "grey", fill = "grey") +
  geom_line(aes(x=time, y=y.smooth))+
  geom_point(aes(x=time, y=X1), colour = "red", size = 6)

AUC(trima$n, trima$y.smooth, method = "spline")


#### Glue it all together ####

# Have 2 dfs: one containing the whole trace and also the 23 rows defining the peak, and one with just the peak definitions
# need to work out how to chop the long df into 23 dfs based upon their start (X3) and end (X4) values, +/- 10 secs
# Once done that, need to apply the AUC function to each df, probably using `map` and have it output a vector
# That vector then needs to be bound to the `peaks_sort` df, and saved as a .csv

#slice.peaks = function(target, x, y) {
#  slice([target$X3]:[target$X4])
#}


# Manually, this works. Having sliced the keying df to one row, we then apply the numbers on columns X3 and X4 to slice
# the trace df

# Looks like the step to preceed this then is to slice the key df, in a function, and push the output of `slice(peaks_sorta$X3:peaks_sorta$X4)`
# to new dfs named by the function

# We'd then write another function that punts that list of dfs into the AUC command
#peaks_sorta <- peaks_sort %>% slice(1)

#trima <- merged %>% slice(peaks_sorta$X3:peaks_sorta$X4)

## Create a sequence
index <- seq(1:dim(peaks)[1])

peaks_a <- cbind(index, peaks)
head(peaks_a)
tail(peaks_a)
head(n_raw)

## Empty data frame
results <- NULL 

## Run through the index 
for(i in seq_along(peaks_a$index)){
  index <- i
  min <- peaks_a[i,4]
  max <- peaks_a[i,5]
  temp <- dplyr::filter(n_raw, time >= min & time <= max)
  area <- AUC(temp$time, temp$y.smooth, method = "spline")
  print(area)
  results_temp <- data.frame("index"=index, "Area_f"=area)
  results <- rbind(results,results_temp)
}


## Merge with returned Peak info
conc <- read_csv("data/concs.csv")
final <- dplyr::left_join(peaks_a, results, by='index') %>% 
  arrange(X3) %>% 
  cbind(conc)
head(final)

## plot

my.formula <- y ~ x

ggplot(final, aes(conc, Area_f)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)

## Save results
write_csv(final, "outputs/Results.csv")


#### TC ####
#### import data ####
c_raw <- read_csv("data/tc.csv")
plot(c_raw, col = 'Gray', type = 'l')

#### fix baseline ####
c_raw$b <- c_raw$conc + abs(min(c_raw$conc))


#### smooth ####
c_raw$y.smooth <- loess(c_raw$b ~ c_raw$time, span=0.012)$fitted
plot(c_raw$time, c_raw$y.smooth, type = 'l')


#### go pracma ####
#Uses {pracma::findpeaks} to identify peaks

peaks <- data.frame(findpeaks(c_raw$y.smooth, 
                              npeaks=30, 
                              threshold=0, 
                              peakpat = "[+]{1,}[0]*[-]{1,}", 
                              sortstr=TRUE)) #IDs peaks and makes a df that has
#                                            # X1 = height at point of picking
#                                            # X2 = time at point of picking
#                                            # X3 = time at start of peak
#                                            # X4 = time at end of peak

c_raw$n <- seq(1,length(c_raw$y.smooth))

# sort peaks df
peaks_sort <- peaks %>% 
  arrange(X2)

#generate peak number vector and combine with sorted peaks df
peak_no <- 1:30
peaks_sort <- cbind(peak_no, peaks_sort)

merged <- merge(x=c_raw, 
                y=peaks_sort, 
                by.x="n", 
                by.y="X2", 
                all.x=TRUE, 
                all.y=TRUE)

ggplot(merged, aes(x=time, y=conc)) +
  geom_col(orientation="x", colour = "grey") +
  geom_line(aes(x=time, y=y.smooth))+
  geom_point(aes(x=time, y=X1), colour = "red")


#### cut early washes ####
trim <- merged %>% slice(16731:49621)

ggplot(trim, aes(x=time, y=conc)) +
  geom_col(orientation="x", colour = "grey", fill = "grey") +
  geom_line(aes(x=time, y=y.smooth))+
  geom_point(aes(x=time, y=X1), colour = "red", size = 6)

AUC(trima$n, trima$y.smooth, method = "spline")


#### Glue it all together ####


## Create a sequence
peaks_trim <- peaks_sort %>% slice(8:30)

index <- seq(1:dim(peaks_trim)[1])

peaks_a <- cbind(index, peaks_trim)
head(peaks_a)
tail(peaks_a)
head(c_raw)

## Empty data frame
results <- NULL 

## Run through the index 
for(i in seq_along(peaks_a$index)){
  index <- i
  min <- peaks_a[i,5]
  max <- peaks_a[i,6]
  temp <- dplyr::filter(c_raw, time >= min & time <= max)
  area <- AUC(temp$time, temp$y.smooth, method = "spline")
  print(area)
  results_temp <- data.frame("index"=index, "Area_f"=area)
  results <- rbind(results,results_temp)
}


## Merge with returned Peak info
conc <- read_csv("data/concsc.csv")
final <- dplyr::left_join(peaks_a, results, by='index') %>% 
  arrange(X3) %>% 
  cbind(conc)
head(final)

## plot

my.formula <- y ~ x

ggplot(final, aes(conc, Area_f)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE, formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)

## Save results
write_csv(final, "outputs/CResults.csv")



