# two betas

library(truncnorm)
library(MASS)
library(coda)

chain_length = 2000

x = cbind( c(rep(1.5,500), rep(1,500)), runif(1000,-5,5) )
true_beta = c(1,2)
p = pnorm(x %*% true_beta)
y = rbinom(1000,1,p)

plot(1,type="n",xlab="i",ylab="beta1",xlim=c(1,chain_length),ylim=c(0,2))
dev.new()
plot(1, type="n",xlab="i",ylab="beta2",xlim=c(1,chain_length),ylim=c(0,3))
color = c("red","purple","black","green","blue")
mcmcoutput = NA

for (iter in 1:5) {
bvals = rbind(rep(NA,chain_length), rep(NA,chain_length))
#bvals[,1] = solve(t(x)%*%x) %*% (t(x) %*% y)
bvals[,1] = c(runif(1,0,3), runif(1,0,3))
for (i in 1:(chain_length-1)) {
  beta = bvals[,i]
  wvals = rep(NA,1000)
  for (j in 1:1000) {
    if (y[j]==1) {
      wvals[j] = rtruncnorm(1, a=0, mean=t(x[j,]) %*% beta)
    } else {
      wvals[j] = rtruncnorm(1, b=0, mean=t(x[j,]) %*% beta)
    }
  }
  bvals[,i+1] = mvrnorm(1, solve(t(x)%*%x) %*% t(x) %*% wvals, solve(t(x) %*% x))
}
dev.set(2)
points(c(1:chain_length), bvals[1,],type="l",col=color[iter])
dev.set(3)
points(c(1:chain_length), bvals[2,],type="l",col=color[iter])
if (iter==1) {  mcmcoutput = t(bvals) } 
else {  mcmcoutput = cbind(mcmcoutput, t(bvals)) }
}

result = mcmc.list(as.mcmc(mcmcoutput[,c(1,2)]),as.mcmc(mcmcoutput[,c(3,4)]),
as.mcmc(mcmcoutput[,c(5,6)]),as.mcmc(mcmcoutput[,c(7,8)]),
as.mcmc(mcmcoutput[,c(9,10)]))
gelman.diag(result)
