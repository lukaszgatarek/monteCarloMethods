rm(list = ls())

library(mcsm)
nSim <- 10000

zA <- rep(NA, nSim)
zB <- rep(NA, nSim)

pA <- rep(NA, nSim)
pB <- rep(NA, nSim)
p0 <- rep(NA, nSim)

nA  <- 186 
nB  <- 38
nAB <- 13
n0  <- 284

probabilities <- rdirichlet(1, c(1,1,1)) 

pA[1] <- probabilities[1]
pB[1] <- probabilities[2]
p0[1] <- probabilities[3]



for (i in 2:nSim){
  
  zA[i] <- rbinom(1, nA, pA[i - 1]^2 / (pA[i - 1]^2 + 2 * pA[i - 1] * p0[i - 1]) )
  zB[i] <- rbinom(1, nB, pB^2 / (pB[i - 1]^2 + 2 * pB[i - 1] * p0[i - 1]) )
                                  
  # define the shape parameters
  alphaS <- zA[i] + nA + nAB + 1
  betaS  <- zB[i] + nB + nAB + 1
  gammaS <- nA + nB - zA[i] - zB[i] + 2 * n0 + 1
  
  probabilities <- rdirichlet(1, c(alphaS, betaS, gammaS))  

  pA[i] <- probabilities[1]
  pB[i] <- probabilities[2]
  p0[i] <- probabilities[3]
  
}

zA[1] <- zA[2]
zB[1] <- zB[2]

hist(pA,fre=FALSE,col="grey",nclass=123,xlim=c( quantile(pA, 0.0000001), quantile(pA, 0.999)), main="",xlab=expression(pA))
hist(pB,fre=FALSE,col="grey",nclass=123,xlim=c( quantile(pB, 0.0000001), quantile(pB, 0.999)), main="",xlab=expression(pB))
hist(p0,fre=FALSE,col="grey",nclass=123,xlim=c( quantile(p0, 0.0000001), quantile(p0, 0.999)), main="",xlab=expression(p0))



