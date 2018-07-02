h = function(x){(cos(50*x)+ sin(20*x))^2}

par(mar=c(2,2,2,1),mfrow=c(2,1))
curve(h, xlab="Function", ylab="", from=1, to = 3, lwd=3)
analyticalIntegral = integrate(h,0,1)

# takes a sample from the uniform and evaluate the h function on it
nSim=10^4
x=h(runif(nSim))
plot(x)
estInt = cumsum(x)/(1:nSim)
# plot(estInt[400:500])
estVar= (cumsum((x-estInt)^2))/(1:nSim)^2
#plot(estVar, type="l", lwd = 1)
#lines(estInt, type="l", lwd = 2, col= "red")  
estErr = sqrt(estVar)
plot(estInt, xlab="Mean and error range", type = "l", lwd=2, ylim =
      + mean(x)+20*c(-estErr[nSim], estErr[nSim]), ylab="")
lines(estInt+2*estErr, col = "gold", lwd=2)
lines(estInt-2*estErr, col = "gold", lwd=2)
lines(analyticalIntegral$value * rep(1, length(estErr)), col = "red")

## Cauchy vs Normal

nSim=10^6
monteCarloIntegrand = function(h, FUN = rnorm, nSim, ifPlot) {
  # first we draw from Normal and evaluate the Cauchy function
  #h = function(x){(1/(1+x^2))}
  par(mar=c(2,2,2,1),mfrow=c(2,1))
#  curve(h, xlab="Function", ylab="", from=0, to = 1, lwd=3)
  analyticalIntegral = integrate(h,0,1)
  
  # takes a sample from the uniform and evaluate the h function on it
  sampleFromPdf <- FUN(nSim)
  #sampleFromPdf<-sampleFromPdf(sampleFromPdf<1 && sampleFromPdf>0)
  
  x<-h(sampleFromPdf)

  estInt = cumsum(x)/(1:nSim)
  # plot(estInt[400:500])
  estVar= (cumsum((x-estInt)^2))/(1:nSim)^2
  #plot(estVar, type="l", lwd = 1)
  #lines(estInt, type="l", lwd = 2, col= "red")  
  estErr = sqrt(estVar)

  if (ifPlot == 1){
    plotCOnvergence(x, estInt, estErr, nSim, analyticalIntegral)
  }
  outputVariables<-list(analyticalIntegral$value, estInt, estErr)
  
  return(outputVariables)
}

# Construct the integral over the product of Cauchy and Normal
h1 = function(x){(x/(1+x^2))}
h2 = function(x){(1/(1+x^2))}

output1 <- monteCarloIntegrand(h1,rnorm, nSim,1) # check, if we change from runif to rnorm 0< <1 then it should still not be equal to analytical integral 
output2 <- monteCarloIntegrand(h2,rnorm, nSim,0)

h3 = function(x){x/exp((x^2)/2)}
h4 = function(x){1/exp((x^2)/2)}

output3 <- monteCarloIntegrand(h3,rcauchy, nSim,0)
output4 <- monteCarloIntegrand(h4,rcauchy, nSim,0)
plot(output1[[2]]/output2[[2]], type="l", lwd=2, ylim =
       + mean(output1[[2]]/output2[[2]])+20*c(-output1[[3]][nSim], output1[[3]][nSim]))

lines(output3[[2]]/output4[[2]], lwd=2, col="red") 


## Importance sampling
plotConvergence = function(x, estInt, estErr, nSim, analyticalIntegral){
  plot(estInt, xlab="Mean and error range", type = "l", lwd=2, ylim =
         + mean(x)+20*c(-estErr[nSim], estErr[nSim]), ylab="")
  lines(estInt+2*estErr, col = "gold", lwd=2)
  lines(estInt-2*estErr, col = "gold", lwd=2)
  lines(analyticalIntegral$value * rep(1, length(estErr)), col = "red")
}


monteCarloIS = function(h, funF = rnorm, funG = rnorm, nSim, ifPlot, paramsFG = list(0,1,0,1)) {
  
# paramsFG contains the parameters of the distribution F and G in a sequence, first F then G
  analyticalIntegral = integrate(h,0,1)
  
  # derive the dprob functions from rprob functions
  funDensityF<-get(paste("d", substring(substitute(funF),2), sep=""))
  funDensityG<-get(paste("d", substring(substitute(funG),2), sep=""))
  
  
  sampleG <- funG(nSim, paramsFG[[3]], paramsFG[[4]])
  # takes a sample from the g and evaluates the integral on it
  x<-h(sampleG) * funDensityF(sampleG, paramsFG[[1]], paramsFG[[2]]) / funDensityG(sampleG, paramsFG[[3]], paramsFG[[4]]) 
  
  estInt = cumsum(x)/(1:nSim)
  # plot(estInt[400:500])
  estVar= (cumsum((x-estInt)^2))/(1:nSim)^2
  #plot(estVar, type="l", lwd = 1)
  #lines(estInt, type="l", lwd = 2, col= "red")  
  estErr = sqrt(estVar)
  
  if (ifPlot == 1){
    plotConvergence(x, estInt, estErr, nSim, analyticalIntegral)
  }
  outputVariables<-list(analyticalIntegral, estInt, estErr)
  
  return(outputVariables)
}



## comparison of standard monte carlo integration versus importance sampling
h = function(x){exp(-(x-3)^2/2) + exp(-(x-6)^2/2)}
nSim = 1000
monteCarloIntegralApproximation = monteCarloIS(h, runif, runif, nSim, 1, list(0,1,0,1))
## Example 3.3.5
# the importance sampling with f=g is the same as standard monte carlo integration
# we define a function that is squared non integrable 
h = function(x){sqrt(abs(x/(1-x)))}
# lets first use the t distribution to generate the sample. In that case the degrees of freedom equal df and the location parameter is 0  
df = 12
mCIAT = monteCarloIS(h, rt, rt, nSim, 1, list(df,0,df,0))
mCIAC = monteCarloIS(h, rt, rcauchy, nSim, 1, list(df,0,0,1))
mCIAN = monteCarloIS(h, rt, rnorm, nSim, 1, list(df,0,0,df/(df-2)))

# save all the results in the list to compare the convergence
listEstint=list(mCIAT[[2]], mCIAC[[2]], mCIAN[[2]])

compareConvergence = function(listEstint){
  plot(listEstint[[1]], type = "l", lwd=2, ylim =
         + mean(listEstint[[1]])+10*c(-sd(listEstint[[1]]), sd(listEstint[[1]])), ylab="")
  for (i in 2:length(listEstint)) {
    lines(listEstint[[i]], lwd=2, col = i) 
    legend("topleft", legend = 1:length(listEstint), col=1:length(listEstint), pch=1)
    } }
compareConvergence(listEstint)

monteCarloIntegralApproximation = monteCarloIS(h, runif, runif, nSim, 1, list(0,1,0,1))






