rm(list = ls())

set.seed(1233)

xm<-rcauchy(500)
f<-function(y){
  -sum(log(1+(x-y)^2))}

mi<-matrix(0,500,1)
f<-function(y){-sin(y*100)^2-sum(log(1+(x-y)^2))}

for (i in 1:500){
  x<-xm[1:i]
  mi[i, 1]<-optimize(f, interval = c(-10,10), maximum = T)$max}

plot(mi, type="l")

# --------------------- example 5.2 ---------------------------

da<-rbind(rnorm(10^2), 2.5+rnorm(3*10^2))

likelihoodFunction<-function(mu){
  -1* sum(log(0.25*dnorm(da-mu[1]) + (0.75*dnorm(da-mu[2]))))
}

mu1Grid<-seq(-2,5,0.1)
mu2Grid<-seq(-1,4,0.1)

lmu1Grid<-length(mu1Grid)
lmu2Grid<-length(mu2Grid)

# variable for storing the evaluations of the likelihood over the grid
likeEval<-matrix(0,lmu1Grid,lmu2Grid)

for (i in 1:lmu1Grid){
  for (j in 1:lmu2Grid){  
    likeEval[i,j]<-likelihoodFunction(c(mu1Grid[i], mu2Grid[j]))      
    }
} 
contour(mu1Grid, mu2Grid, likeEval)



library(MASS)
Sigma <- diag(2)
startingPointsList<-mvrnorm(n = 100, rep(0, 2), Sigma)

mmu<-array(NA,dim=c(100,2,100))

for (j in 1:100){
 startingPoint<-c(-1,-1)
#  startingPoint<-startingPointsList[j,]

  mmu[1,,j]<-startingPoint

  for (i in 1:(nlm(likelihoodFunction, startingPoint)$iterations)){
    mmu[1+i,,j]<-nlm(likelihoodFunction, startingPoint, iterlim=i)$estimate }

}

contour(mu1Grid, mu2Grid, likeEval)
lines(mmu[,,1], pch=19, lwd=4)

#------------------------ example 5.3 --------------------------

h = function(x){ (cos(50*x) + sin(20*x))^2 }
x<-seq(0,1,0.01)
plot(h)

optimum<-optimize(h, interval = c(0,1), maximum = T)
argMax<-optimum$maximum
argMaxVal<-optimum$objective

# variability of a uniform sampler
rangom<-h(matrix(runif(10^5), ncol=10^3)) # sample multiple sequences from uniform
monitor<- t(apply(rangom,1,cummax)) # returns a vector whose elements are the cumulative maxima

plot(monitor[1,], type="l")
for (i in 2:100){
  lines(monitor[i,],  pch=19, lwd=4)
}
polygon(c(1:10^3,10^3:1), c(apply(monitor, 2, max), rev(apply(monitor, 2, min))), col="grey")
abline(h=optimise(h, int=c(0, 1), maximum=T)$ob, lwd=7)


cau<-rcauchy(10^2)
mcau<-median(cau)
rcau<-diff(quantile(cau, c(.25, .75)))
f<-function(x){
  z<-dcauchy(outer(x, cau, FUN = "-"))
  apply(z,1,mean)
}
fcst<-integrate(f, lower=-10, upper=10)

ft<-function(x){f(x) }

g<-function(x){dt((x-mcau)/rcau, df=49) / rcau}

curve(ft, from = -10, to = 10)
curve(g, add = T, type = "b")


# ----------------------------------
unisan<-matrix(f(runif(5*10^4,-5,5)), ncol=500)
dft<-2
causan<-matrix(f(rt(5*10^4, df=dft) * rcau + mcau), ncol=500)
unimax<-apply(unisan,2,cummax)[10:10^2,]
caumax<-apply(causan,2,cummax)[10:10^2,]
plot(caumax[,1], col="white", ylim=c(.8,1)* max(causan))

polygon(c(10:10^2,10^2:10), c(apply(unimax, 1, max), rev(apply(unimax, 1, min))), col="grey")
polygon(c(10:10^2,10^2:10), c(apply(caumax, 1, max), rev(apply(caumax, 1, min))), col="wheat")


#----------------- Example 5.6 ---------------------------
h <- function(x,y){(x*sin(20*y)+y*sin(20*x))^2 + cosh(sin(10*x)*x) + (x*cos(10*y) - y * sin(10*x))^2 * cosh(cos(20*y)*y)}
x<-seq(-3,3,le=435)
y<-x
z<-outer(x,y,h)
par(bg="wheat", mar=c(1,1,1,1))
persp(x,y,z, theta=155, phi=30, col="green4", ltheta=-120, shade=.75, border=NA, box=FALSE)  
  

# stochastisc search algorithm
start<-c(.65, .8)
theta<-matrix(start, ncol=2)
diff=iter=1
alpha<-1
beta<-1
# initialize space for gradients (we dont know the number of iterations a prior, so it should be just defined as an empty matrix)
grad<-matrix(NA, nrow = 0, ncol = 2)
argNumber<-2 # number of argumnets in the function


while (diff>10^-8){

  alpha[iter]<- 1 / ((iter + 1) )
  beta[iter]<- 1 / (((iter + 1))^0.1)
  
  zeta<-runif(argNumber) * 10
  zeta<-zeta / sqrt(t(zeta) %*% zeta) # normalization with the variance

  # compute the gradients in both directions
  
  grad<-rbind(grad, matrix(NA, nrow = 1, ncol = 2)) # extend the matrix by one row
  
  grad[iter, 1] <- alpha[iter] * zeta[1] * (-1) * ( h(theta[iter,1] + beta[iter]*zeta[1], theta[2]) - h(theta[iter,1]-beta[iter]*zeta[1], theta[2]) ) / (2*beta[iter]) 
  grad[iter, 2] <- alpha[iter] * zeta[2] * (-1) * ( h(theta[1], theta[iter,2] + beta[iter]*zeta[2]) - h(theta[1], theta[iter,2]-beta[iter]*zeta[2]) ) / (2*beta[iter]) 

  theta<-rbind(theta, theta[iter,] + grad[iter,])
  
  print(theta[iter, ])
  
  #diff<-sqrt(t(theta[iter, ])%*%theta[iter, ])
  
  scale<-sqrt(t(grad[iter, ])%*%grad[iter, ])

    print(scale)
    
#  while ( (scale>1) || sum(abs(theta[iter+1,])>3) ){
  while (scale>1){
    zeta<-runif(2) * 10
    zeta<-zeta / sqrt(t(zeta) %*% zeta) # normalization with the variance
    grad[iter, 1] <- alpha[iter] * zeta[1] * (-1) * ( h(theta[iter,1] + beta[iter]*zeta[1], theta[iter, 2]) - h(theta[iter,1]-beta[iter]*zeta[1], theta[iter, 2]) ) / (2*beta[iter]) 
    grad[iter, 2] <- alpha[iter] * zeta[2] * (-1) * ( h(theta[iter, 1], theta[iter,2] + beta[iter]*zeta[2]) - h(theta[iter, 1], theta[iter,2]-beta[iter]*zeta[2]) ) / (2*beta[iter]) 
    
    theta[iter+1,] <- theta[iter,] + grad[iter,]
    
    scale<-sqrt(t(grad[iter, ])%*%grad[iter, ])  
    }
  
  diff<-sqrt(t(grad[iter, ])%*%grad[iter, ])  
    
  iter<-iter+1
    
}


#--------------------- stocjastic approximation of h, probit model---------------
library(MASS)
da<-matrix(NA,200,2)
da[,1]<-Pima.tr$type
da[,1]<-da[,1]-1
da[,2]<-Pima.tr$bmi

probitLikelihood<-function(a,b){
  apply(pnorm(-a-outer(X=b, Y=da[,2], FUN="*"), log=T)*(1-da[,1]) + pnorm(a+outer(X=b, Y=da[,2], FUN="*"), log=T)*da[,1], 1, sum)
}

like<-1
for (i in 1:1000){
  like[i]<-probitLikelihood(i,2) }
plot(like)

# marginalizing theta1 to infere about theta0
margap<-function(a){
  b<-rt(10^3, df=5)
  dtb<-dt(b, 5, log=T)
  b = b*.3+.1
  themar<-0
  for (i in 1:10^3){
    themar <- themar + exp( probitLikelihood(a, b[i]) - dtb[i] )
  }
  return(themar)
}

# marginalizing theta0 to infere about theta1
margbp<-function(b){
  a<-rt(10^3, df=5)
  dta<-dt(a, 5, log=T)
  a = a*.3+.1
  themar<-0
  for (i in 1:10^3){
    themar <- themar + exp( probitLikelihood(a[i], b) - dta[i] )
  themar<-themar/10^3
  }
  return(themar)
}

#there are two level of stochasticity. 
# inside of the margap we need 1000 of simulations from the t. 
# then we run the entire procedure 100 times to approximate the distriution in each point a

# --------------------------- approximation of posterior of theta 0 ------------------------------
nSim<-100
seqTheta0<-seq(-4,1,0.1)
marginalPosteriorApproximationTheta0<-matrix(NA, nSim, length(seqTheta0))
index<-1

for (theta0 in seqTheta0){
  for (sim in 1:nSim){
    print(index)
    marginalPosteriorApproximationTheta0[sim, index] <- margap(theta0)
  }
  index<-index + 1
}


#xaxis<--4:1
plot(seqTheta0, marginalPosteriorApproximationTheta0[1,], ylim=c(0, max(marginalPosteriorApproximationTheta0))) 
for (i in 2:nSim){
  lines(seqTheta0, marginalPosteriorApproximationTheta0[i,])
}
polygon( c(seqTheta0, rev(seqTheta0)), c(apply(marginalPosteriorApproximationTheta0, 2, max), rev(apply(marginalPosteriorApproximationTheta0, 2, min))), col="grey")
lines(seqTheta0, apply(marginalPosteriorApproximationTheta0, 2, mean), lwd = 2)



# --------------------------- approximation of posterior of theta 1 ------------------------------
nSim<-100
seqTheta1<-seq(-1,1,0.1)
marginalPosteriorApproximationTheta1<-matrix(NA, nSim, length(seqTheta1))
index<-1

for (theta1 in seqTheta1){
  for (sim in 1:nSim){
    print(index)
    marginalPosteriorApproximationTheta1[sim, index] <- margbp(theta1)
  }
  index<-index + 1
}


#xaxis<--4:1
plot(seqTheta1, marginalPosteriorApproximationTheta1[1,], ylim=c(0, max(marginalPosteriorApproximationTheta1))) 
for (i in 2:nSim){
  lines(seqTheta1, marginalPosteriorApproximationTheta1[i,])
}
polygon(c(seqTheta1, rev(seqTheta1)), c(apply(marginalPosteriorApproximationTheta1, 2, max), rev(apply(marginalPosteriorApproximationTheta1, 2, min))), col="grey")
lines(seqTheta1, apply(marginalPosteriorApproximationTheta1, 2, mean), lwd = 2)


# example --------------------- 5.12 ------------------------
# simulate data from the mixture of two normals with he ratio of probability 1/4 to 3/4 and the corresponding means of 0 and 2.5
dataSimMix<-cbind(t(rnorm(1*10)), 6 + t(rnorm(3*10)) )
dataSimMix<-dataSimMix[sample(1:length(dataSimMix), replace = F)]


mixtureLikelihood<-function(mu1,mu2,z){
  # sum( log( 0.25*dnorm(dataSimMix[z==1]-mu1) ) ) + sum( log(0.75*dnorm(dataSimMix[z==0]-mu2) ) ) 
  sum(log( 0.25*dnorm(dataSimMix[z==1]-mu1))) + sum( log(0.75*dnorm(dataSimMix[z==0]-mu2) ) ) 
  }

nSim<-1000 # replication of likelihood evaluation for each point on the grid

margMixturePosterior<-function(mu1,mu2){
  z<-matrix(NA, 10^3, length(dataSimMix) )
  for (i in 1:10^3){
    z[i,]<-rbinom(length(dataSimMix), 1, 0.25)
  }
  dz<-dbinom(z, 1, 0.25, log=T)
  
  themar<-0
  for (i in 1:10^3){
    themar <- themar + exp( mixtureLikelihood(mu1,mu2,z[i,]) - sum(dz[i,]) )
    themar<-themar/10^3
    }
  return(themar)
}

mu1Grid<-seq(-2,1,1)
mu2Grid<-seq(0,10,1)

lmu1Grid<-length(mu1Grid)
lmu2Grid<-length(mu2Grid)

# variable for storing the evaluations of the likelihood over the grid
marginalPosteriorApproximationMixture<-array(0, dim = c(nSim, lmu1Grid, lmu2Grid) )

index<-1
for (i in 1:lmu1Grid){
  mu1<-mu1Grid[i]
  for (j in 1:lmu2Grid){
    mu2<-mu2Grid[j]
    for (sim in 1:nSim){
      print(index)
      marginalPosteriorApproximationMixture[sim, i, j] <- margMixturePosterior(mu1, mu2)
      print(log(marginalPosteriorApproximationMixture[sim, i, j]))
    }
    index<-index + 1
  }
}

marginalPosteriorMixtureAverage<-array(0, dim = c(lmu1Grid, lmu2Grid) )

for (i in 1:lmu1Grid){
  for (j in 1:lmu2Grid){
    marginalPosteriorMixtureAverage[i,j]<-mean(marginalPosteriorApproximationMixture[,i,j])
  }
}

contour(mu1Grid, mu2Grid, log(marginalPosteriorMixtureAverage) )




# one dimensional analysis. did not work
marginalPosteriorApproximationMixtureOneDim<-matrix(NA, nSim, lmu1Grid)
index=1
for (i in 1:lmu1Grid){
  mu1<-mu1Grid[i]
    mu2<-2.5
    for (sim in 1:nSim){
      print(index)
      marginalPosteriorApproximationMixtureOneDim[sim, i] <- margMixturePosterior(mu1, mu2)
    }
    index<-index + 1
  }


plot(mu1Grid, marginalPosteriorApproximationMixtureOneDim[1,], ylim=c(0, max(marginalPosteriorApproximationMixtureOneDim))) 

for (i in 2:nSim){
  lines(mu1Grid, marginalPosteriorApproximationMixtureOneDim[i,])
}
polygon(c(mu1Grid, rev(mu1Grid)), c(apply(marginalPosteriorApproximationMixtureOneDim, 2, max), rev(apply(marginalPosteriorApproximationMixtureOneDim, 2, min))), col="grey")
lines(mu1Grid, apply(marginalPosteriorApproximationMixtureOneDim, 2, mean), lwd = 2)

