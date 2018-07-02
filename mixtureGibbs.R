rm(list = ls())

# exercise --------------------- 7.10 ------------------------
# simulate data from the mixture of two normals with he ratio of probability 1/4 to 3/4 and the corresponding means of 0 and 2.5
x<-cbind(t(rnorm(9*100)), 2.7 + t(rnorm(1*100)) )
hist(x[1,])
n<-length(x)
# reorder the sample
x<-x[sample(1:n, replace = F)]
hist(x)

# number of components in the mixture
k<-2
v<-1
sigma2<-1
# number of simulations of the MCMC
nSim <- 2000
# matrix to keep the sampled mus
mus <- matrix(0, nSim, k)
p <- 0.7
Z <- matrix(0, nSim, n)
nj <- matrix(0,nSim,k)

probZ1<-matrix(0, nSim, n)
probZ2<-matrix(0, nSim, n)
probBinom<-matrix(0, nSim, n)
# to start the sampler draw the Zs with the probability 1-p, then the succes is getting in the second mixture component
Z[1,] <- rbinom(n, 1, 1-p) 

for (i in 2:nSim){
  for (j in 1:k){
    # number of observations in the jth mixtrue component
    nj[i,j] = sum(Z[i-1,]==j-1)
    mus[i,j]<-rnorm(1, m = v^2 * sum(x[Z[i-1,]==j-1]) / (nj[i,j] * v^2 + 1), sd = sigma2 * v^2 / (nj[i,j] * v^2 + 1)) 
  }  
  
  # condition to avoid the labeling problem
  mus[i,] <-sort(mus[i,], decreasing = FALSE)
  
  for (t in 1:n){
    probZ1[i,t] <- p * exp( -(x[t]-mus[i,1])^2 / 2*sigma2 )
    probZ2[i,t] <- (1-p) * exp( -(x[t]-mus[i,2])^2 / 2*sigma2 )
    
    # we assume the probability of success to be assigned to group 2
    probBinom[i,t] <- probZ2[i,t] / (probZ1[i,t] + probZ2[i,t]) 
  
    Z[i,t] <- rbinom(1, 1, probBinom[i,t]) 
  }
  
print(i)  
}


