################## Lukasz Gatarek, 09-01-2017
################## implementation of gelfant et al 1990 (JASA)
################## posterior simulator for linear hierarchical model

rm(list = ls())
# necessary libraries
library('MASS')
library('MCMCpack')
#

# write a data generating process
eta <- matrix(NA, 2, 1)
eta[1, 1] <- 100
eta[2, 1] <- 15

C <- matrix(0, 2, 2)
C[1,1] <- 40^2
C[2,2] <- 10^2

rho <- 5
R <- matrix(0, 2, 2)
R[1,1] <-10^2
R[2,2] <- .5^2

a <-3
b <- 1/40

# number of individuals
n <- 30
# number of measurements per individual
m <- 5

# dimension of the explanatory variables: constant factor plus variables associated with betas
d <- 2

# draw the true parameters
trueSIGMA <- riwish(v = rho, S = rho * R)
trueTheta0 <- mvrnorm(n = 1, mu = eta, Sigma = C)
trueSigma2 <- rinvgamma(1, a, b)

trueTheta <- array(NA, dim = c(n, d, 1))
for (i in 1:n){
  # draw theta_i = [alpha_i, beta_i] for each individual 
  trueTheta[i,,] <- mvrnorm(n = 1, mu = trueTheta0, Sigma = trueSIGMA) }

trueX <- matrix(NA, n, m)
# gelfand data
for (i in 1:n){
  trueX[i, 1] <- 8
  trueX[i, 2] <- 15
  trueX[i, 3] <- 22
  trueX[i, 4] <- 29
  trueX[i, 5] <- 36 }

Y <- matrix(NA, n, m)
trueAlpha <- matrix(NA, n, 1)
trueBeta <- matrix(NA, n, 1)

for (i in 1:n){
  trueAlpha[i,1] <- trueTheta[i,1,1]
  trueBeta[i,1] <- trueTheta[i,2,1]
  for (j in 1:m){
    Y[i,j] <- rnorm(1, trueAlpha[i,1] + trueBeta[i,1] * trueX[i,j], sqrt(trueSigma2))
  }
}

#################
# for the MCMC algorithm we assume the number of simulations
nSim <- 1000
# and prepare the space for the variables
SIGMA <- array(data = NA, dim = c(nSim, d, d))
theta0 <- array(data = NA, dim = c(nSim, d, 1))
theta <- array(data = NA, dim = c(nSim, d, n))
sigma2 <- matrix(data = NA, nSim, 1)
# average of thetas
thetaBar <- array(data = NA, dim = c(nSim, d, 1))
squares <- rep(NA, n)

thetaSquares <- vector(mode="list", length = n)
#################


# simulate the first Sigma given rho and R from the inverted Wishart distribution
# wishartSigma <- matrix(rho * R, d, d)
# c <- matrix(rWishart(1, rho, ginv(wishartSigma)), d, d) 
# sigmas[1,,] <- ginv(c) # the inverse function does not work extremely well here 

# instead of drawing precision form wishart, we draw sigmas directly frm the inverted wishart
SIGMA[1,,] <- riwish(v = rho, S = rho * R)
theta0[1,,] <- mvrnorm(n = 1, mu = eta, Sigma = C)
sigma2[1] <- rinvgamma(1, a, b)

for (t in 2:nSim){
  for (i in 1:n){
    # complete posterior conditional for theta_i
    X_i <- cbind(rep(1, m), trueX[i,])
    # page 175
    D_theta_i <- ginv( t(X_i) %*% X_i / sigma2[t-1] + ginv(SIGMA[t-1,,]) )
    d_theta_i <- t(X_i) %*% matrix(Y[i,]) / sigma2[t-1] + ginv(SIGMA[t-1,,]) %*% matrix(theta0[t-1,,])
    
    theta[t,,i] <- mvrnorm(n = 1, mu = D_theta_i %*% d_theta_i, Sigma = D_theta_i)
  }
  # posterior of theta_0
  D_theta_0 <- ginv(n * ginv(SIGMA[t-1,,]) + ginv(C))
  # compute the average of thetas in draw t
  thetaBar[t,,] = apply(theta[t,,], 1, mean)
  d_theta_0 <- (n * ginv(SIGMA[t-1,,]) %*% thetaBar[t,,] + ginv(C) %*% eta)
  theta0[t,,1] <- mvrnorm(n = 1, mu = D_theta_0 %*% d_theta_0, Sigma = D_theta_0)
  
  for (i in 1:n){
    X_i <- cbind(rep(1, m), trueX[i,])  
    squares[i] <- t(matrix(Y[i,]) - X_i %*% matrix(theta[t,,i])) %*% (matrix(Y[i,]) - X_i %*% matrix(theta[t,,i]) )
    }
  sigma2[t] <- rinvgamma(1, (m * n)/2 + a, ginv( 1/2 * sum(squares) + b^(-1)  ) )

  # posterior of Sigma

  for (i in 1:n){
    thetaSquares[[i]] <- matrix(theta[t,,i] - theta0[t,,1]) %*% t(matrix(theta[t,,i] - theta0[t,,1]))
  }
  
  SIGMA[t,,] <- riwish(v = n + rho, S = Reduce("+", thetaSquares )  + rho * R )
  
}

#par(mfrow = c(2, 1), mar = c(2, 1, 1, 1))
hist(sigma2[(nSim/10):nSim], col = "grey", breaks = 25, xlab = "", 
     main = expression(sigma^2), fre = F)
abline(v = trueSigma2, lwd = 5, col = 'red')
summary(sigma2)

hist(theta0[(nSim/10):nSim,1,1], col = "grey", breaks = 25, xlab = "", 
     main = expression(alpha_0), fre = F)
abline(v = trueTheta0[1], lwd = 5, col = 'red')
summary(theta0[,1,1])

hist(theta0[(nSim/10):nSim,2,1], col = "grey", breaks = 25, xlab = "", 
     main = expression(beta_0), fre = F)
abline(v = trueTheta0[2], lwd = 5, col = 'red')
summary(theta0[,2,1])

for (i in 1:n){
  hist(theta[(nSim/10):nSim,1,i], col = "grey", breaks = 25, xlab = "", 
       main = paste(expression(alpha_), toString(i)), fre = F)
  abline(v = trueTheta[i,1,1], lwd = 5, col = 'red')
}

for (i in 1:n){
  hist(theta[(nSim/10):nSim,2,i], col = "grey", breaks = 25, xlab = "", 
       main = paste(expression(beta_), toString(i)), fre = F)
  abline(v = trueTheta[i,2,1], lwd = 5, col = 'red')
}

# SIGMA
hist(SIGMA[(nSim/10):nSim,1,1], col = "grey", breaks = 25, xlab = "", 
     main = expression(SIGMA_1_1), fre = F)
abline(v = trueSIGMA[1,1], lwd = 5, col = 'red')

hist(SIGMA[(nSim/10):nSim,1,2], col = "grey", breaks = 25, xlab = "", 
     main = expression(SIGMA_1_2), fre = F)
abline(v = trueSIGMA[1,2], lwd = 5, col = 'red')

hist(SIGMA[(nSim/10):nSim,2,2], col = "grey", breaks = 25, xlab = "", 
     main = expression(SIGMA_2_2), fre = F)
abline(v = trueSIGMA[2,2], lwd = 5, col = 'red')
