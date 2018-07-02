rm(list = ls())

library(mcsm)

xdata = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
Time = c(94.32, 15.72, 62.88, 125.76, 5.24, 31.44, 1.05, 1.05, 2.1, 10.48)

alpha <- 1.8
gamma <- 0.01
delta <- 1

nSim <- 50000
nPumps <- 10

beta <- rep(0.1, nSim)
lambda <- matrix(0, nSim, nPumps)

for (i in 2:nSim){
  for (j in 1:nPumps){
    lambda[i,j] <- rgamma(1, sh = xdata[j] + alpha, ra = Time[j] + beta[i-1])
    beta[i] <- rgamma(1, sh = gamma + nPumps * alpha, ra = delta + sum(lambda[i,]))
  }
}

xAxisLimit <-4

par(mfrow = c(5, 2), mar = c(4, 4, 1, 1))

hist(lambda[(nSim/10):nSim, 1], col = "grey", breaks = 25, xlim = c(0,xAxisLimit), xlab = "", 
     main = expression(lambda_1), fre = F)
hist(lambda[(nSim/10):nSim, 2], col = "grey", breaks = 25, xlim = c(0,xAxisLimit), xlab = "", 
     main = expression(lambda_2), fre = F)
hist(lambda[(nSim/10):nSim, 3], col = "grey", breaks = 25, xlim = c(0,xAxisLimit), xlab = "", 
     main = expression(lambda_3), fre = F)
hist(lambda[(nSim/10):nSim, 4], col = "grey", breaks = 25, xlim = c(0,xAxisLimit), xlab = "", 
     main = expression(lambda_4), fre = F)
hist(lambda[(nSim/10):nSim, 5], col = "grey", breaks = 25, xlim = c(0,xAxisLimit), xlab = "", 
     main = expression(lambda_5), fre = F)
hist(lambda[(nSim/10):nSim, 6], col = "grey", breaks = 25, xlim = c(0,xAxisLimit), xlab = "", 
     main = expression(lambda_6), fre = F)
hist(lambda[(nSim/10):nSim, 7], col = "grey", breaks = 25, xlim = c(0,xAxisLimit), xlab = "", 
     main = expression(lambda_7), fre = F)
hist(lambda[(nSim/10):nSim, 8], col = "grey", breaks = 25, xlim = c(0,xAxisLimit), xlab = "", 
     main = expression(lambda_8), fre = F)
hist(lambda[(nSim/10):nSim, 9], col = "grey", breaks = 25, xlim = c(0,xAxisLimit), xlab = "", 
     main = expression(lambda_9), fre = F)
hist(lambda[(nSim/10):nSim, 10], col = "grey", breaks = 25, xlim = c(0,xAxisLimit), xlab = "", 
     main = expression(lambda_10), fre = F)




