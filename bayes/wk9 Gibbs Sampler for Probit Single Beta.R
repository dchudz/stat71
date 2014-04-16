# Single beta 

library(truncnorm)
library(MASS)
library(coda)

x = c(rep(1,500),rep(2,500))
true_beta = 1
p = pnorm(x*true_beta)
y = rbinom(1000,1,p)

plot(1,type="n",xlab="i",ylab="beta",xlim=c(1,1000),ylim=c(0,2))
color = c("red","purple","black","green","blue")

mcmcoutput = NA
for (iter in 1:5) {
bvals = rep(NA,1000)
bvals[1] = runif(1,-3,3)

for (i in 1:999) {
  beta = bvals[i]
  wvals = rep(NA,1000)
  for (j in 1:1000) {
    if (y[j]==1) {
      wvals[j] = rtruncnorm(1, a=0, mean= x[j]*beta)
    } else {
      wvals[j] = rtruncnorm(1, b=0, mean= x[j]*beta)
    }
  }
  bvals[i+1] = rnorm(1, solve(t(x)%*%x) %*% t(x) %*% wvals, solve(t(x) %*% x))
}
points(c(1:1000), bvals,type="l",col=color[iter])

if (iter==1) {  mcmcoutput = bvals } 
else {  mcmcoutput = cbind(mcmcoutput, bvals) }
}

result = mcmc.list(as.mcmc(mcmcoutput[,1]),as.mcmc(mcmcoutput[,2]),
as.mcmc(mcmcoutput[,3]),as.mcmc(mcmcoutput[,4]),
as.mcmc(mcmcoutput[,5]))
gelman.diag(result)
dev.new()
gelman.plot(result)
