rm(list = ls())

library("MultinomialCI")
simVector<-rmultinom(1, size = 12, prob = c(0.6,0.1,0.2,0.8))

trueVector<-matrix(NA,4,1)
trueVector[1]<-125
trueVector[2]<-18
trueVector[3]<-20
trueVector[4]<-34

x<-trueVector


## ------------------- deterministic sequence ------------------
# take the starting point
theta<-0.2

for (i in 1:1000){
  theta<-cbind(theta,NA)
  theta[i+1]<-( (theta[i]/(2+theta[i])) * x[1] + x[4] ) / ((theta[i]/(2+theta[i])) * x[1] + x[2] + x[3] + x[4]) 

}

plot(theta[1,], type='l')

## ----------------- numerical simulation ----------------------

nObservationsBinom<-1000
z<-matrix(NA, nObservationsBinom, 1)
successProbBinom<-matrix(NA, nObservationsBinom, 1)
thetaNum<-0.2

for (i in 1:nObservationsBinom){
  thetaNum<-cbind(thetaNum,NA)
  successProbBinom[i] <- thetaNum[i] / (2 + thetaNum[i]) 
  z[i,1]<-rbinom(1, x[1], successProbBinom[i])
  thetaNum[i+1]<-( mean(z[1:i,1]) + x[4] ) / ( mean(z[1:i,1]) + x[2] + x[3] + x[4]) 
}

lines(thetaNum[1,], col="red")

## ----------------- Gibbs sampler ------------------------------

nSim <- 50000
# Allocate memory  for gibbs samples
thetaGibbs <- zGibbs <- rep(0.5, nSim)

a <- 500
b <- 500

# thetaGibbs[1] <- rbeta(1,a,b)

for (j in 2:nSim){
  thetaGibbs[j] <- rbeta(1, a + zGibbs[j-1] + trueVector[4] + 1, b + trueVector[2] + trueVector[3] + 1)
  zGibbs[j] <- rbinom(1, trueVector[1], thetaGibbs[j] / ( 2 + thetaGibbs[j] ) )
}

hist(thetaGibbs[(nSim/10):nSim], col = "grey", breaks = 250, xlab = "", 
     main = expression(theta), fre = F)

mean(thetaGibbs)
# -----------------------------------------------------------------
  
  





