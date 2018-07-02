# amopunt of students (tests)
n <- 100
# number of simulations in the MCMc
nSim <- 1000
# number of test questions to answer
Q <- 10

# true probability of a succes in group 1 (randomizing students)
p1True <- 0.2
# true probability of a succes in group 2 (nonrandomizing students)
p2True <- 0.8
# probability of belonging to group 1
piTrue <- 0.05
# write a process to simulate the date
componentsTrue <- rep(NA, n)
# vector of true data
y <- rep(NA, n)

for (i in 1:n){
# first choose the component
  componentsTrue[i] <- rbinom(1,1,piTrue)
  # sample the number of correct answers within the specific component
  if (componentsTrue[i] == 1){
    y[i] <- rbinom(1,Q,p1True) 
  } else {
    y[i] <- rbinom(1,Q,p2True)
  }
}

# ------------- hyperparameters ---------------------------
a1 = a2 = a3 = b1 = b2 = b3 = 2

# probability of a succes in group 1 (randomizing students)
p1<-matrix(NA, nSim, 1)
# probability of a succes in group 2 (nonrandomizing students)
p2<-matrix(NA, nSim, 1)

pPi <- matrix(NA, nSim, 1)

tau <- matrix(NA, nSim, n)

# draw the first values from the prior
p1[1,1] <- rbeta(1, a2, b2)
p2[1,1] <- rbeta(1, a3, b3)
pPi[1,1] <-rbeta(1, a1, b1)

# sample first artifical adherence to components (pPi corresponds to the succes of getting in the randomizing)
tau[1,] <- rbinom(n,1,pPi[1,1])

# likelihood evaluations for n observations for component 1; across nSim simulations of MCMC
prob1Group <- matrix(NA, nSim, n)
# likelihood evaluations for n observations for component 2; across nSim simulations of MCMC
prob2Group <- matrix(NA, nSim, n)



for (i in 2:nSim){
  p1[i, 1] <- rbeta(1, a2 + sum(tau[i-1,] * y), b2 + sum(tau[i-1,] * (Q - y)) )
  p2[i, 1] <- rbeta(1, a3 + sum( (1 - tau[i-1,]) * y), b3 + sum( (1 - tau[i-1,]) * (Q - y)) )
  pPi[i, 1] <- rbeta(1, a1 + sum(tau[i-1,]), b1 + sum(1 - tau[i-1,]) )
  
  
  # condition to avoid the labeling problem
  if (p1[i, 1] > p2[i, 1]){
    pholder <- p1[i, 1] # not to forget the value
    p1[i, 1] <- p2[i, 1]
    p2[i, 1] <- pholder
    pPi[i, 1] <- 1- pPi[i, 1]
  }
  
  
  for (j in 1:n){
    prob1Group[i,j] <- pPi[i,1] * p1[i,1]^y[j] * (1 - p1[i,1])^(Q-y[j])
    prob2Group[i,j] <- (1 - pPi[i,1]) * p2[i,1]^y[j] * (1 - p2[i,1])^(Q-y[j])
    
    tau[i,j] <- rbinom(1, 1, prob1Group[i,j] / (prob1Group[i,j] + prob2Group[i,j]) )
  }
  
  
}



par(mfrow = c(3, 1), mar = c(4, 4, 1, 1))

hist(p1[(nSim/10):nSim, 1], col = "grey", breaks = 25, xlim = c(0,1), xlab = "", 
     main = expression(p_1), fre = F)
abline(v = p1True, lwd = 5, col = 'red')

hist(p2[(nSim/10):nSim, 1], col = "grey", breaks = 25, xlim = c(0,1), xlab = "", 
     main = expression(p_2), fre = F)
abline(v = p2True, lwd = 5, col = 'red')

hist(pPi[(nSim/10):nSim, 1], col = "grey", breaks = 25, xlim = c(0,1), xlab = "", 
     main = expression(pi), fre = F)
abline(v = piTrue, lwd = 5, col = 'red')


