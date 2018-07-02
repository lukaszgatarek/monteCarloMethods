f=function(a,b){exp(2*(lgamma(a+b) - lgamma(a) - lgamma(b)) + a * log(.3) + b * log(.2))}
aa=1:150
bb=1:100
posterior = outer(aa,bb,f)
image(aa,bb,posterior, xlab=expression(alpha), ylab=" ")
contour(aa,bb,posterior, add=T)