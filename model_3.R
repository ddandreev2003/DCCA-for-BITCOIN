library("quantmod")
library('aTSA')
library("forecast")
library("tseries")
library("FinTS")
library("moments")
library("rugarch")
library("lmtest")
library("multicool")
library("rugarch")
library(MASS)
library(stats)
library("ggplot2")
library("e1071")
library('corrplot')
library('vars')
library('stringr')

DCCA_CC=function(x,y,k){
  ## calculate cumulative sum profile of all t
  xx<- cumsum(x - mean(x))  ## Equation 2
  yy<- cumsum(y - mean(y))  ## Equation 2
  
  ## Divide in to overlapping boxes of size k
  
  slide_win_xx = mat_sliding_window(xx,k)
  slide_win_yy = mat_sliding_window(yy,k)
  ## calculate linear fit value in each box 
  x_hat = t(apply(slide_win_xx,1,function(n) (lm(n~seq(1:length(n)))$fitted.values)))
  y_hat = t(apply(slide_win_yy,1,function(n) (lm(n~seq(1:length(n)))$fitted.values)))
  
  ##  Get detrend variance in each box with linear fit value (detrend by local trend).
  F2_dfa_x = c()
  F2_dfa_y = c()
  for(i in 1:nrow(x_hat)){
    ## Equation 4
    F2_dfa_x = c(F2_dfa_x,mean((xx[i:(i+k-1)]-x_hat[i,])^2))
  }
  for(i in 1:nrow(y_hat)){
    ## Equation 4
    F2_dfa_y = c(F2_dfa_y,mean((yy[i:(i+k-1)]-y_hat[i,])^2))
  }
  ## Average detrend variance over all boxes to obtain fluctuation
  F2_dfa_x = mean(F2_dfa_x) ## Equation 3
  F2_dfa_y = mean(F2_dfa_y) ## Equation 3
  
  ## Get detrended covariance of two profile
  F2_dcca = c()
  for(i in 1:nrow(x_hat)){
    ## Equation 5
    F2_dcca = c(F2_dcca,mean((xx[i:(i+k-1)]-x_hat[i,]) * (yy[i:(i+k-1)]-y_hat[i,]) ))
  }
  
  ## Equation 6
  F2_dcca = mean(F2_dcca)
  
  ## Calculate correlation coefficient 
  rho = F2_dcca / sqrt(F2_dfa_x * F2_dfa_y) ## Equation 1
  return(rho)
}
mat_sliding_window = function(xx,k){
  ## Function to generate boxes given dataset(xx) and box size (k)
  slide_mat=c()
  for (i in 1:(length(xx)-k+1)){
    slide_mat = rbind(slide_mat,xx[i:(i+k-1)] )
  }
  return(slide_mat)
}


dt <- read.csv("C:/Users/Dell/Documents/data_table.csv",header = TRUE,sep = ';')
swiss <- read.csv("E:/Games/GBP=X.csv", header = TRUE,sep = ',')
btc <- read.csv('E:/Games/BTC-USD (1).csv', header = TRUE,sep = ',')
nikel <- read.csv('E:/Games/^N225.csv', header = TRUE,sep = ',')
eth <- read.csv('E:/Games/ETH-USD.csv', header = TRUE,sep = ',')
gold <- read.csv('E:/Games/GLD.csv', header = TRUE,sep = ',')
LTC <-  read.csv('E:/Games/LTC-USD.csv', header = TRUE,sep = ',')
swiss <- swiss$Open
btc <- btc$Open
nikel <- nikel$Open
nikel <- as.numeric(nikel)
gold <- gold$Open
str_replace_all(gold, ',','.')
gold <- as.numeric(gold)
LTC <- LTC$Open
swiss <- swiss[1:261]
btc <- btc[1:261]
nikel <- nikel[1:261]
gold <- gold[1:261]
LTC <- LTC[1:261]

btc <- na.omit(btc)
gold <- na.omit(gold)
LTC <- na.omit(LTC)
nikel <- na.omit(nikel)
swiss <- na.omit(swiss)

swiss <- swiss[1:260]
btc <- btc[1:260]
nikel <- nikel[1:260]
gold <- gold[1:260]
LTC <- LTC[1:260]

a <- c(DCCA_CC(btc,btc,6), DCCA_CC(btc,gold,6),DCCA_CC(btc,LTC,6)
       ,DCCA_CC(btc,nikel,6), DCCA_CC(btc,swiss,6))
b <- c(DCCA_CC(gold,btc,6), DCCA_CC(gold,gold,6),DCCA_CC(gold,LTC,6)
       ,DCCA_CC(gold,nikel,6), DCCA_CC(gold,swiss,6))
C <- c(DCCA_CC(LTC,btc,6), DCCA_CC(LTC,gold,6),DCCA_CC(LTC,LTC,6)
       ,DCCA_CC(LTC,nikel,6), DCCA_CC(LTC,swiss,6))
d <- c(DCCA_CC(nikel,btc,6), DCCA_CC(nikel,gold,6),DCCA_CC(nikel,LTC,6)
       ,DCCA_CC(nikel,nikel,6), DCCA_CC(nikel,swiss,6))
e <- c(DCCA_CC(swiss,btc,6), DCCA_CC(swiss,gold,6),DCCA_CC(swiss,LTC,6)
       ,DCCA_CC(swiss,nikel,6), DCCA_CC(swiss,swiss,6))

M_matrix <- rbind(a,b,C,d,e)
rownames(M_matrix) <- c('btc','gold','LTC','nikel','swiss')
colnames(M_matrix) <- c('btc','gold','LTC','nikel','swiss')
corrplot(M_matrix,method = 'pie')
