################## Lukasz Gatarek, 09-01-2017
################## implementation of gelfant et al 1990 (JASA)
################## posterior simulator for linear hierarchical model

rm(list = ls())
# necessary libraries
library('MASS')
library('MCMCpack')
#

## symulowanie danych zgodnie z procesem generujący dane
# parametry dla rozkładu prior dla parametru theta_0
eta <- matrix(NA, 2, 1)
eta[1, 1] <- 100
eta[2, 1] <- 15

C <- matrix(0, 2, 2)
C[1,1] <- 40^2
C[2,2] <- 10^2

# parametry dla rozkładu prior dla parametru Sigma
rho <- 5
R <- matrix(0, 2, 2)
R[1,1] <-10^2
R[2,2] <- .5^2

# parametry dla rozkładu prior dla parametru sigma2
a <-3
b <- 1/40

# liczba obiektów
n <- 30
# ilość pomiarów dla każdego z obiektów
m <- 5
# wymiar zmiennej objaśniającej, wliczając parametr stałej
d <- 2

# losowanie parametrów procesu generującego dane z odpowiednich rozkładów a priori
trueSIGMA <- riwish(v = rho, S = rho * R)
trueTheta0 <- mvrnorm(n = 1, mu = eta, Sigma = C)
trueSigma2 <- rinvgamma(1, a, b)

# losowanie parametrów theta_i
trueTheta <- array(NA, dim = c(n, d, 1))
for (i in 1:n){
  trueTheta[i,,] <- mvrnorm(n = 1, mu = trueTheta0, Sigma = trueSIGMA) }

# obserwacje X_{i,m}
trueX <- matrix(NA, n, m)
for (i in 1:n){
  for (j in 1:m){  
  trueX[i, j] <- rnorm(1,0,1) } }

#trueX <- matrix(NA, n, m)
# gelfand data
#for (i in 1:n){
#  trueX[i, 1] <- 8
#  trueX[i, 2] <- 15
#  trueX[i, 3] <- 22
#  trueX[i, 4] <- 29
#  trueX[i, 5] <- 36 }


# nazwy obiektów 
objectNames<-'obiekt 1'
for (i in 2:n){
  objectNames <- c(objectNames, paste('obiekt', toString(i)))
}
rownames(trueX)<-objectNames
  
# obserwacje Y zgodnie z procesem generującym dane 
Y <- matrix(NA, n, m)
rownames(Y)<-objectNames

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
# zakładamy ilość symulacji w łańcuchu MCMC
nSim <- 1000
# alokujemy pamięć pod wyniki losowania ...
SIGMA <- array(data = NA, dim = c(nSim, d, d))
theta0 <- array(data = NA, dim = c(nSim, d, 1))
theta <- array(data = NA, dim = c(nSim, d, n))
sigma2 <- matrix(data = NA, nSim, 1)
# ... oraz pod dodatkowe zmienne konieczne do losowania
thetaBar <- array(data = NA, dim = c(nSim, d, 1))
squares <- rep(NA, n)
thetaSquares <- vector(mode="list", length = n)
#################


# simulate the first Sigma given rho and R from the inverted Wishart distribution
# wishartSigma <- matrix(rho * R, d, d)
# c <- matrix(rWishart(1, rho, ginv(wishartSigma)), d, d) 
# sigmas[1,,] <- ginv(c) # the inverse function does not work extremely well here 

# wartości początkowe
SIGMA[1,,] <- riwish(v = rho, S = rho * R)
theta0[1,,] <- mvrnorm(n = 1, mu = eta, Sigma = C)
sigma2[1] <- rinvgamma(1, a, b)

for (t in 2:nSim){
  for (i in 1:n){
    # warunkowy rozkład a posteriori dla theta_i
    X_i <- cbind(rep(1, m), trueX[i,])
    # page 175
    D_theta_i <- ginv( t(X_i) %*% X_i / sigma2[t-1] + ginv(SIGMA[t-1,,]) )
    d_theta_i <- t(X_i) %*% matrix(Y[i,]) / sigma2[t-1] + ginv(SIGMA[t-1,,]) %*% matrix(theta0[t-1,,])
    
    theta[t,,i] <- mvrnorm(n = 1, mu = D_theta_i %*% d_theta_i, Sigma = D_theta_i)
  }
  # warunkowy rozkład a posteriori dla  theta_0
  D_theta_0 <- ginv(n * ginv(SIGMA[t-1,,]) + ginv(C))
  # średnia dla theta_i's
  thetaBar[t,,] = apply(theta[t,,], 1, mean)
  d_theta_0 <- (n * ginv(SIGMA[t-1,,]) %*% thetaBar[t,,] + ginv(C) %*% eta)
  theta0[t,,1] <- mvrnorm(n = 1, mu = D_theta_0 %*% d_theta_0, Sigma = D_theta_0)
  
  # warunkowy rozkład a posteriori dla sigma2
  for (i in 1:n){
    X_i <- cbind(rep(1, m), trueX[i,])  
    squares[i] <- t(matrix(Y[i,]) - X_i %*% matrix(theta[t,,i])) %*% (matrix(Y[i,]) - X_i %*% matrix(theta[t,,i]) )
    }
  sigma2[t] <- rinvgamma(1, (m * n)/2 + a, ginv( 1/2 * sum(squares) + b^(-1)  ) )

  # warunkowy rozkład a posteriori dla  Sigma
  for (i in 1:n){
    thetaSquares[[i]] <- matrix(theta[t,,i] - theta0[t,,1]) %*% t(matrix(theta[t,,i] - theta0[t,,1]))
  }
  SIGMA[t,,] <- riwish(v = n + rho, S = Reduce("+", thetaSquares )  + rho * R )
}

# comparison with OLS estimation
estimatedAlpha <- matrix(NA, n, 1)
estimatedBeta <- matrix(NA, n, 1)
for (i in 1:n){
estimatedAlpha[i,1] <- mean(theta[2:nSim,1,i])
estimatedBeta[i,1]  <- mean(theta[2:nSim,2,i])
}
estimatOLS<-matrix(NA, n, 2)
for (i in 1:n){
    X_i <- cbind(rep(1, m), trueX[i,])  
    
    estimatOLS[i,]<- ginv(t(X_i) %*% X_i) %*% t(X_i) %*% (matrix(Y[i,]))
    
}


plot(estimatedAlpha, type="l")
lines(trueTheta[,1,1], col="red")
lines(estimatOLS[,1], col="green")

