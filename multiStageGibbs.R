rm(list = ls())

p<-5
nSim<-10000000
rho <- 0.44

x<-matrix(NA, p, nSim)

# initialize the chain
x[,1]<-0

a <- 1 + (p - 2) * rho

for (i in 2:nSim){
  for (j in 1:p){
    # the x that shouldbe conditioned on is a combination of current and previous column
    if (j>1 && j < p){
      currentlyAvailableXs <- cbind(x[ (j+1) :p, i-1], x[1:(j-1), i]) 
     } else if (j == 1) {
      currentlyAvailableXs <-x[ (j+1) :p, (i-1)]
     } else if (j == p) {
      currentlyAvailableXs <-x[ 1: (p-1), i] 
     }
        x[j,i] <- rnorm(1, m =  ( (p-1) * rho) * mean(currentlyAvailableXs) / a, sd = sqrt( (a - (p - 1) * rho^2) / a) )
  }
}

cov(t(x))
cor(t(x))

