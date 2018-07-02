h = function(x){(cos(50*x)+ sin(20*x))^2}
x<-matrix(h(runif(200*10^4)), ncol=200)
estint<-apply(x,2,cumsum)/(1:(10^4))

plot(estint[,1], type="l", ylim=c(0.8, 1.2))
y<-apply(estint,1,quantile,c(.025,.975))
polygon(c(1:10^4,10^4:1), c(y[1,], rev(y[2,])), col="red")

analyticalIntegral = integrate(h,0,1)


# alternatively we can obtain the bootstrapped replicas of the original vector

boot<- matrix(sample(x[,1],200*10^4, rep=T), nrow = 10^4, ncol = 200)
bootit<-apply(boot,2,cumsum)/(1:10^4)
bootup<-apply(bootit,1,quantile, .975 )
bootdown<-apply(bootit,1,quantile, .025 )

plot(bootit[,1], type="l", ylim=c(0.8, 1.2))
polygon(c(1:10^4,10^4:1), c(bootdown, rev(bootup)), col="wheat")

# example 4.2

# draw from normal with mean 2.5
norma<-matrix(rnorm(500*10^4), ncol=500) + 2.5
weit<-1/(1+norma^2)
# posterior expectation of the mean 
esti<-apply(norma*weit,2,cumsum)/apply(weit,2,cumsum)
plot(esti[,1], type="l", col="white",ylim=c(1.7,1.9))
band<-apply(esti,1,quantile, c(.025,.975))
polygon(c(1:10^4,10^4:1), c(band[1,], rev(band[2,])), col="wheat")

vare<-cumsum(weit[,1]*norma[,1]^2) / cumsum(weit[,1]) - esti[,1]^2
lines(esti[,1]+2*sqrt(vare/(1:10^4)))
lines(esti[,1]-2*sqrt(vare/(1:10^4)))

wis<-weit/apply(weit,2,cumsum)
sumwis2<-apply(wis^2,2,cumsum)
ess2<-(sumwis2)^(-1)
essbo2<-apply(ess2,1,quantile, c(.025,0.975))
plot(ess2[,1], type="l", col="white")
polygon(c(1:10^4,10^4:1), c(essbo2[1,], rev(essbo2[2,])), col="wheat")


ess<-apply(weit,2,cumsum)^2/apply(weit^2,2,cumsum)
essbo<-apply(ess,1,quantile, c(.025,0.975))

plot(ess[,1], type="l", col="white")
polygon(c(1:10^4,10^4:1), c(essbo[1,], rev(essbo[2,])), col="wheat")

# Example 4.5

Nsim<-10^4
norma<-rnorm(Nsim)+2.5
hnorm<-norma * dcauchy(norma)
munorm<-mean(hnorm)
sdnorm<-sd(hnorm)
f<-function(x)(cumsum(hnorm))[round(Nsim*x)]/round(x*Nsim)
f(2)
