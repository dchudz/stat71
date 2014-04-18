x = runif(1000,-3,3) #set these values for x

true_beta = .5
y = rbinom(1000, 1, pnorm(x*true_beta))

beta_0 = solve(t(x)%*%x)%*%t(x)%*%y #least squares guess 

beta = beta_0
for(j in 1:20) { #right now, just 20 gibbs samples, seems to converge quickly
  print(beta)
  W = c() #vector of conditionals, given y's and our current beta
  for(i in 1:1000){
    if(y[i]==1) {
      w = rtruncnorm(1, a=0, mean = x[i]*beta) #constrained to be > 0
    }
    else {
      w = rtruncnorm(1, b=0, mean = x[i]*beta) #constrained to be < 0
    }
    W = c(W, w)
  }
  
  #new guess for beta
  new_beta = rnorm(1, solve(t(x)%*%x)%*%t(x)%*%cbind(W), solve(t(x)%*%x))
  beta = new_beta
}
