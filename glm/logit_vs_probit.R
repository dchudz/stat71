logitCDF <- function(t) {
  return(exp(t)/(exp(t)+1))
}

n = 100  			
beta=1				# Actual beta, TODO: 0.5, -1
betaList = seq(-5, 5, .1)	# Range of beta to plot
x = runif(n, -5, 5)		# Set random x vals between -5 and 5

N <- 100
lambda_vec <- c()
for(i in 1:N) {
  #Simulate Yi's using model 1
  u = runif(n)
  Y1 = rep(0,n)
  T1 = log(u/(1-u))

  for (i in 1:n) {
    if (x[i]*beta + T1[i]>0) {Y1[i]=1}
  }
  
  #calculate MLEs using glm
  logit <- glm(Y1~x, family=binomial()) 
  probit <- glm(Y1~x, family=binomial(link="probit")) 
  
  mle_logit <- logit$coefficients[-1]
  mle_probit <- probit$coefficients[-1]
  
  #calculate log-likelihoods for MLEs
  lik_logit <- sum(Y1*log(1 - logitCDF(-x*mle_logit))
                   + (1-Y1)*log(logitCDF(-x*mle_logit)))
  
  lik_probit <- sum(Y1*log(1 - pnorm(-x*mle_probit))
                    + (1-Y1)*log(pnorm(-x*mle_probit)))
  
  lambda <- lik_logit - lik_probit
  
  lambda_vec <- c(lambda_vec, lambda)
  
}

hist(lambda_vec)
