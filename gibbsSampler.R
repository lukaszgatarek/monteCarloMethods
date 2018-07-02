# ---------------- Example 7.3 --------------------------
library("datasets")
library("MCMCpack")
library("pscl")
x<-c(91,504,557,609,693,727,764,803,857,929,970, 1043,1089,1195,1384,1713)
xbar<-mean(x)

n <- length(x)
a = b = 3
tau2 <- 10
theta0 <- 5

Nsim<-50000
sh1 <- (n/2) + a

sigma2 = theta = rep(0, Nsim)
sigma2[1] <- 1 / rgamma(1, shape = a, rate = b)
B <- sigma2[1]/(sigma2[1] + n * tau2)

theta[1] <- rnorm(1, m = B * theta0 + (1-B) * xbar, sd = sqrt(tau2 * B))

for (i in 2:Nsim) {
  B <- sigma2[i-1] / (sigma2[i-1] + n * tau2)
  
  theta[i] <- rnorm(1, m = B * theta[i-1] + (1 - B) * xbar, sd = sqrt(tau2 * B) )
  ra1 <- (1/2) * sum((x-theta[i])^2) + b^(-1)
  #sigma2[i] <- 1/rgamma(1, shape = sh1, rate = ra1)
  sigma2[i] <- rigamma(1, sh1, ra1)
  
}

mean(theta)
mean(sqrt(sigma2))

hist(log(theta))
hist(log(sqrt(sigma2)))


# The inverse-Gamma density for x>0 with parameters α>0 and β>0 is

# (beta^alpha)/Gamma(alpha) x^(-alpha-1) exp(-beta/x)

