rm(list = ls())
# ---------------- Example 7.5 --------------------------
library("mcsm")

Girls<-c(91,504,557,609,693,727,764,803,857,929,970, 1043,1089,1195,1384,1713)
Boys<-c(457, 645, 736, 790,899,991,1104,1154,1203,1320,1417,1569,1685,1843,2296,2710)

Girls<-log(Girls)
Boys<-log(Boys)

x<-rbind(Girls, Boys)
k<-dim(x)[1]

xbar<-apply(x,1,mean)

n1<-length(Girls)
n2<-length(Boys)

x1bar <- mean(Girls)
x2bar <- mean(Boys)

meanXbar = (n1 * x1bar + n2 * x2bar)/(n1 + n2)
n <- c(n1, n2) # in general case they can differ

a1 = a2 = a3 = 10
b1 = b2 = b3 = 30
mu0 <- meanXbar # assume the same theta0 for each variable in the system

Nsim<-5000

shSig <- (sum(n)/2) + a1
shTau <- (k/2) + a2
shSigMu <- (1/2) + a3

theta = matrix(0, Nsim, k)
sigma2 = tau2 = mu = thetabar = sigmaMu2 = rep(0, Nsim)

#sigma2[1]   <- 1 / rgamma(1, shape = a1, rate = b1)
#tau2[1]     <- 1 / rgamma(1, shape = a2, rate = b2)
#sigmaMu2[1] <- 1 / rgamma(1, shape = a3, rate = b3)

sigma2[1]= tau2[1] = sigmaMu2[1] = (var(Girls) + var(Boys))/2

B <- sigma2[1]/(sigma2[1] + n * tau2)

for (j in 1:k) {
  
  theta[1,j] <- rnorm(1, m = B * mu0 + (1-B) * xbar[j], sd = sqrt(tau2[1] * B)) }

for (i in 2:Nsim){
  for (j in 1:k) {
    
    Bsig <- sigma2[i-1] / (sigma2[i-1] + n[j] * tau2[i-1])
    theta[i,j] <- rnorm(1, m = Bsig * mu[i-1] + (1 - Bsig) * xbar[j], sd = sqrt(tau2[i-1] * Bsig) )
    
    Btau <- tau2[i-1] / (tau2[i-1] + k * sigmaMu2[i-1])
    thetabar[i]<-mean(theta[i,])
    mu[i] <- rnorm(1, m = Btau * mu0 + (1 - Btau) * thetabar[i], sd = sqrt(sigmaMu2[i-1] * Btau) )
  
    raSig <- (1/2) * sum( (apply(x,2,function(x) x - theta[i,] ) )^2 ) + b1
    sigma2[i] <- 1/(rgamma(1, shape = shSig, rate = raSig) )
    
    raTau <- (1/2) * sum( (theta[i,] - mu[i])^2 ) + b2
    tau2[i] <- 1/(rgamma(1, shape = shTau, rate = raSig) )
    
    raSigMu2 <- (1/2) * (mu[i] - mu0)^2 + b3
    sigmaMu2[i] <- 1/(rgamma(1, shape = shSigMu, rate = raSigMu2) )
  }
}

par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))
top = max(mu)
bot = min(mu)
hist(mu[(Nsim/10):Nsim], col = "grey", breaks = 25, xlab = "", 
     main = expression(mu), fre = F)
hist(theta[,1][(Nsim/10):Nsim], col = "grey", breaks = 25, xlab = "", 
     main = expression(theta[1]), fre = F)
hist(theta[,2][(Nsim/10):Nsim], col = "grey", breaks = 25, xlab = "", 
     main = expression(theta[2]), fre = F)
hist(sqrt(sigmaMu2[(Nsim/10):Nsim]), col = "sienna", breaks = 50, 
     xlim = c(0.5, 4), xlab = "", main = expression(sigma[mu]), 
     fre = F)
hist(sqrt(tau2[(Nsim/10):Nsim]), col = "sienna", breaks = 25, 
     xlim = c(0.5, 4), xlab = "", main = expression(tau), 
     fre = F)
hist(sqrt(sigma2[(Nsim/10):Nsim]), col = "sienna", breaks = 25, 
     xlim = c(0.5, 2), xlab = "", main = expression(sigma), 
     fre = F)


# ------------- solution from the book --------------------------

  nsim = 10^3
  a = 10
  b = 30
  

  x1 = log(Girls)
  x2 = log(Boys)
  n1 = length(x1)
  n2 = length(x2)
  n = n1 + n2
  x1bar = mean(x1)
  x2bar = mean(x2)
  xbar = (n1 * x1bar + n2 * x2bar)/n
  k = 2
  a1 = a2 = a3 = a
  b1 = b2 = b3 = b
  mu = rep(xbar, nsim)
  theta1 = rep(x1bar, nsim)
  theta2 = rep(x2bar, nsim)
  sigma2mu = sigma2 = tau2 = rep((var(x1) + var(x2))/2, nsim)
  mu0 = xbar
  for (i in 2:nsim) {
    B1 = sigma2[i - 1]/(sigma2[i - 1] + n1 * tau2[i - 1])
    theta1[i] = rnorm(1, mean = B1 * mu[i - 1] + (1 - B1) * 
                        x1bar, sd = sqrt(tau2[i - 1] * B1))
    B2 = sigma2[i - 1]/(sigma2[i - 1] + n2 * tau2[i - 1])
    theta2[i] = rnorm(1, mean = B2 * mu[i - 1] + (1 - B2) * 
                        x2bar, sd = sqrt(tau2[i - 1] * B2))
    B = tau2[i - 1]/(k * sigma2mu[i - 1] + tau2[i - 1])
    m = B * mu0 + (1 - B) * (n1 * theta1[i] + n2 * theta2[i])/n
    mu[i] = rnorm(1, mean = m, sd = sqrt(sigma2mu[i - 1] * 
      
                                                                                B))
    sh1 = (n/2) + a1
    ra1 = (1/2) * (sum((x1 - theta1[i])^2) + sum((x2 - theta2[i])^2)) + 
      b1
    sigma2[i] = 1/rgamma(1, shape = sh1, rate = ra1)
    sh2 = (k/2) + a2
    ra2 = (1/2) * (sum((theta1[i] - mu[i])^2) + sum((theta2[i] - 
                                                       mu[i])^2)) + b2
    tau2[i] = 1/rgamma(1, shape = sh2, rate = ra2)
    sh3 = (1/2) + a3
    ra3 = (1/2) * (mu[i] - mu0)^2 + b3
    sigma2mu[i] = 1/rgamma(1, shape = sh3, rate = ra3)
  }
  par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))
  top = max(mu)
  bot = min(mu)
  hist(mu[(nsim/10):nsim], col = "grey", breaks = 25, xlab = "", 
       main = expression(mu), fre = F)
  hist(theta1[(nsim/10):nsim], col = "grey", breaks = 25, xlab = "", 
       main = expression(theta[1]), fre = F)
  hist(theta2[(nsim/10):nsim], col = "grey", breaks = 25, xlab = "", 
       main = expression(theta[2]), fre = F)
  hist(sqrt(sigma2mu[(nsim/10):nsim]), col = "sienna", breaks = 50, 
       xlim = c(0.5, 4), xlab = "", main = expression(sigma[mu]), 
       fre = F)
  hist(sqrt(tau2[(nsim/10):nsim]), col = "sienna", breaks = 25, 
       xlim = c(0.5, 4), xlab = "", main = expression(tau), 
       fre = F)
  hist(sqrt(sigma2[(nsim/10):nsim]), col = "sienna", breaks = 25, 
       xlim = c(0.5, 2), xlab = "", main = expression(sigma), 
       fre = F)



