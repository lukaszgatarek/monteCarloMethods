################## Lukasz Gatarek, 09-01-2017
################## implementation of gelfant et al 1990 (JASA) with extension to multiple explanatory variables
################## posterior simulator for linear hierarchical model

rm(list = ls())
# necessary libraries
library('MASS')
library('MCMCpack')
library('mvtnorm')

require(xlsx)
library(mice)

setwd("C:/Users/lgatarek/Dropbox/Rcodes/monteCarloMethodsWithR")

df <- read.xlsx("Dane_mod.xlsx", sheetName = "Dane", as.data.frame = T)
setCountries <- unique(df$Country) 

countrySpecificDataCorrected <- NULL

# count the countires coming into analysis
numberCountries <- length(setCountries)
#n <- 27 # liczba krajów uwzględnionych w analizie
n <- numberCountries

  for (i in 1: n) {
    print(setCountries[i])
  
    countrySpecificData <-df[df$Country == setCountries[i], ]
    
    tempData <-mice(countrySpecificData[,1:25],m=5,maxit=50,meth='cart',seed=500)
    tempDataDF<-data.frame(complete(tempData))
  
    countrySpecificDataCorrected <- rbind(countrySpecificDataCorrected, tempDataDF)
  }

# uzupełnienie na całości


    
    countrySpecificData <-df
    tempData <-mice(countrySpecificData[,1:25],m=5,maxit=50,meth='cart',seed=500)
    tempDataDF<-data.frame(complete(tempData))
    
    countrySpecificDataCorrected <-tempDataDF


#############
# print how many nas in each column
affectedVariablesCheckUp<-apply(countrySpecificDataCorrected,2,function(x) sum(is.na(x)) )
affectedVariables <- affectedVariablesCheckUp[affectedVariablesCheckUp!=0]




#countrySpecificDataCorrectedCopy
## print how many nas in each column
#affectedVariablesCheckUp<-apply(countrySpecificDataCorrectedCopy,2,function(x) sum(is.na(x)) )
#affectedVariables <- affectedVariablesCheckUp[affectedVariablesCheckUp!=0]

# which countries affected
affectedCountries <-list()
for (av in 1:length(affectedVariables)){
  
affectedCountries[[av]] <- ( unique(countrySpecificDataCorrected[is.na(countrySpecificDataCorrected[,names(affectedVariables[av])]), 'Country']) ) } 

unique(unlist(affectedCountries))

#affectedCountries <- ( unique(countrySpecificDataCorrected[is.na(countrySpecificDataCorrected[,'GINI']), 'Country']) ) 

#affectedCountries <- ( unique(countrySpecificDataCorrectedCopy[is.na(countrySpecificDataCorrectedCopy[,names(affectedVariables)]), 'Country']) ) 


countrySpecificDataCorrectedCopy<- countrySpecificDataCorrected
for (country in 1:(length(affectedCountries) - 1)){
  countrySpecificDataCorrectedCopy<-countrySpecificDataCorrectedCopy[countrySpecificDataCorrectedCopy[,2]!=affectedCountries[country], ] }


matrix(rbind(countrySpecificDataCorrected[,2], is.na(countrySpecificDataCorrected[,c(3)]) ))

countrySpecificDataCorrected[is.na(countrySpecificDataCorrected[,'Itspend_GDP']), 'Country']
#########################33

setExplanatory <- names(countrySpecificDataCorrected) # all exisiting explanatory variables

indexes <- 3:25
indexes <- indexes[-18] # wyrzucic zmienna objasniana

#combinations <- combinations[-15]

k <- 1
explanatoryVariablesCombination <- list()

for (j in 4:length(indexes) ){
  combinations <- combn(indexes, j)
  numCombinations <- dim(combinations)[2]
  
  for (i in 1:numCombinations){
    explanatoryVariablesCombination[[k]] <- setExplanatory[combinations[,i]]  
    k <- k + 1
  }
  print(c(j,i))
}
#explanatoryVariablesCombination[[1]] <- setExplanatory[c(3,4,5)] # wybieramy zmienne wyjasniajace
#explanatoryVariablesCombination[[2]] <- setExplanatory[c(3,4)]
#explanatoryVariablesCombination[[3]] <- setExplanatory[c(4,5)]
#explanatoryVariablesCombination[[4]] <- setExplanatory[c(3,6)]
#explanatoryVariablesCombination[[5]] <- setExplanatory[c(4,6)]
#explanatoryVariablesCombination[[6]] <- setExplanatory[c(5,6)]
#explanatoryVariablesCombination[[7]] <- setExplanatory[c(4,5,6)]
#explanatoryVariablesCombination[[8]] <- setExplanatory[c(3,6,7)]
#explanatoryVariablesCombination[[9]] <- setExplanatory[c(4,5,8)]
#explanatoryVariablesCombination[[10]] <- setExplanatory[c(6,7,8)]

noModels <- length(explanatoryVariablesCombination)
marginalLikelihood <- rep(NA, noModels)

for (model in 1:noModels){
 # for (model in 690:numCombinations){
  print(model)
  explanatoryVariables <- explanatoryVariablesCombination[[model]]
  
  p <- length(explanatoryVariables) # ilosc zmiennych wyjasniajacych (oprócz stałej)

eta <- matrix(NA, p + 1, 1)
eta[1, 1] <- 0.5 # dla parametru stałej
# dla zmiennych X
for (i in 1:p){
  eta[i + 1, 1] <- 0.5 }
#eta[3,1] <- 40

C <- matrix(0, 1 + p, 1 + p)
C[1,1] <- 40^2
for (i in 1:p){
  C[i + 1, i + 1] <- 10^2 }

# parametry dla rozkładu prior dla parametru Sigma
rho <-  (p+1)
R <- matrix(0, p + 1, p + 1)
R[1,1] <-10^2
for (i in 1:p){
  R[i + 1, i + 1] <- .5^2 }

# parametry dla rozkładu prior dla parametru sigma2
a <-3
b <- 1/40

# wymiar zmiennej objaśniającej, wliczając parametr stałej
d <- p + 1

# ilość pomiarów dla każdego kraju
m <- 12

trueX <- array(NA, dim=c(n, m, p)) # array przechowujacy zmienne wyjasniajace
Y <- matrix(NA, n, m) # array przechowujacy zmienna wyjasniana 


#TO DO: zmianna wyjasniana powinna byc uzmieniona
for (k in 1:p) {
  for (i in 1:n){
    for (j in 1:m){  
      trueX[i, j, k] <- (countrySpecificDataCorrected[countrySpecificDataCorrected$Country == setCountries[i], explanatoryVariables])[j,k]    
      Y[i,j] <- countrySpecificDataCorrected[countrySpecificDataCorrected$Country== setCountries[i], 20][j] # EgovIndex jako zmienna wyjasniana
      } } }
rownames(trueX) <- setCountries[1:n]
rownames(Y) <- setCountries[1:n]


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
 # print(t)
  for (i in 1:n){
    # warunkowy rozkład a posteriori dla theta_i
    X_i <- cbind(rep(1, m), trueX[i,,])
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
    X_i <- cbind(rep(1, m), trueX[i,,])  
    squares[i] <- t(matrix(Y[i,]) - X_i %*% matrix(theta[t,,i])) %*% (matrix(Y[i,]) - X_i %*% matrix(theta[t,,i]) )
    }
  sigma2[t] <- rinvgamma(1, (m * n)/2 + a, ginv( 1/2 * sum(squares) + b^(-1)  ) )

  # warunkowy rozkład a posteriori dla  Sigma
  for (i in 1:n){
    thetaSquares[[i]] <- matrix(theta[t,,i] - theta0[t,,1]) %*% t(matrix(theta[t,,i] - theta0[t,,1]))
  }
  SIGMA[t,,] <- riwish(v = n + rho, S = Reduce("+", thetaSquares )  + rho * R )
}

# par(mfrow = c(5, 2), mar = c(2, 1, 1, 1))
# #par(mfrow = c(10, 3), mar = c(2, 1, 1, 1))
# hist(sigma2[(nSim/10):nSim], col = "grey", breaks = 25, xlab = "", 
#      main = expression(sigma^2), fre = F)
# #abline(v = trueSigma2, lwd = 5, col = 'red')
# #summary(sigma2)
# 
# 
# hist(theta0[(nSim/10):nSim,1,1], col = "grey", breaks = 25, xlab = "", 
#      main = expression(alpha_0), fre = F)
# #abline(v = trueTheta0[1], lwd = 5, col = 'red') 
# 
# for (k in 2:(p+1)) {
#   hist(theta0[(nSim/10):nSim,k,1], col = "grey", breaks = 25, xlab = "", 
#        main = paste(expression(beta_0), explanatoryVariables[k - 1] ), fre = F) }
# #abline(v = trueTheta0[k], lwd = 5, col = 'red') }
# 
# for (i in 1:n){
#   hist(theta[(nSim/10):nSim,1,i], col = "grey", breaks = 25, xlab = "", 
#        main = paste(expression(alpha_), setCountries[i]) , fre = F )
# #  abline(v = trueTheta[i,1,1], lwd = 5, col = 'red')
# }
# 
# for (i in 1:n){
#   for (k in 2:(p+1)) { # przesunięcie ze względu na stałą
#     hist(theta[(nSim/10):nSim,k,i], col = "grey", breaks = 25, xlab = "", 
#          main = paste(expression(beta_), setCountries[i], '_', explanatoryVariables[k - 1] ), fre = F)
#    # abline(v = trueTheta[i,k,1], lwd = 5, col = 'red')
#   }
# }
# 
# # SIGMA
# hist(SIGMA[(nSim/10):nSim,1,1], col = "grey", breaks = 25, xlab = "", 
#      main = expression(SIGMA_1_1), fre = F)
# #abline(v = trueSIGMA[1,1], lwd = 5, col = 'red')
# 
# hist(SIGMA[(nSim/10):nSim,1,2], col = "grey", breaks = 25, xlab = "", 
#      main = expression(SIGMA_1_2), fre = F)
# #abline(v = trueSIGMA[1,2], lwd = 5, col = 'red')
# 
# hist(SIGMA[(nSim/10):nSim,2,2], col = "grey", breaks = 25, xlab = "", 
#      main = expression(SIGMA_2_2), fre = F)

# marginal likelihood evaluation
theta[t,,i] <- mvrnorm(n = 1, mu = D_theta_i %*% d_theta_i, Sigma = D_theta_i)


probTheta_i <-rep(NA, n) # vector of probabilities for mean of theta_i conditional on mean of other parameters

for (i in 1:n){
  # warunkowy rozkład a posteriori dla theta_i
  X_i <- cbind(rep(1, m), trueX[i,,])
  # page 175
  D_theta_i <- ginv( t(X_i) %*% X_i / mean(na.omit(sigma2)) + ginv(apply(na.omit(SIGMA),2:3,mean)) )
  d_theta_i <- t(X_i) %*% matrix(Y[i,]) / mean(na.omit(sigma2)) + ginv(apply(na.omit(SIGMA),2:3,mean)) %*% matrix(apply(na.omit(theta0[,,]),2,mean))
  
 # theta[t,,i] <- mvrnorm(n = 1, mu = D_theta_i %*% d_theta_i, Sigma = D_theta_i)
  
  probTheta_i[i] <-dmvnorm( apply(na.omit(theta[,,i]),2,mean), mean = D_theta_i %*% d_theta_i, sigma = D_theta_i, log=FALSE)
}

exp(sum(log(probTheta_i)))


probTheta0 <-rep(NA, nSim)
# warunkowy rozkład a posteriori dla  theta_0
D_theta_0 <- ginv(n * ginv(apply(na.omit(SIGMA),2:3,mean)) + ginv(C))
# średnia dla theta_i's
for (t in 2: nSim){
  thetaBar[t,,] = apply(theta[t,,], 1, mean)
  d_theta_0 <- (n * ginv(apply(na.omit(SIGMA),2:3,mean)) %*% thetaBar[t,,] + ginv(C) %*% eta)
  probTheta0[t] <- dmvnorm( apply(na.omit(theta0[,,]),2,mean), mean = D_theta_0 %*% d_theta_0, sigma = D_theta_0)
}

probSigma2 <-rep(NA, nSim)
# warunkowy rozkład a posteriori dla sigma2
for (t in 2: nSim){
  for (i in 1:n){
    X_i <- cbind(rep(1, m), trueX[i,,])  
      squares[i] <- t(matrix(Y[i,]) - X_i %*% matrix(theta[t,,i])) %*% (matrix(Y[i,]) - X_i %*% matrix(theta[t,,i]) )
  }
 probSigma2[t] <- dinvgamma(mean(na.omit(sigma2)), (m * n)/2 + a, ginv( 1/2 * sum(squares) + b^(-1)  ) ) 
}


# warunkowy rozkład a posteriori dla  Sigma
probSIGMA <- rep(NA, nSim)
for (t in 2: nSim){
  for (i in 1:n){
    thetaSquares[[i]] <- matrix(theta[t,,i] - apply(na.omit(theta0[,,]),2,mean) )  %*% t(matrix(theta[t,,i] - apply(na.omit(theta0[,,]),2,mean)  ) )
  }
  probSIGMA[t] <- diwish( apply(na.omit(SIGMA),2:3,mean), v = n + rho, S = Reduce("+", thetaSquares )  + rho * R )
}

probTheta2_Theta1 <- rep(NA, nSim)
for (t in 2: nSim){
  probTheta2_Theta1[t] <- probSIGMA[t] *  probSigma2[t] * probTheta0[t]
}

marginalLikelihood[model] <- prod(probTheta_i) * mean(probTheta2_Theta1[2:nSim])
}
