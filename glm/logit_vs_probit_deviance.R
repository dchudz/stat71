#generate data from probit
set.seed(1)

probLower = vector(length=1000)
devLogit = vector(length=1000)
devProbit = vector(length=1000)

for(i in 1:1000){      
  x = rnorm(1000)
  y = rbinom(n=1000, size=1, prob=pnorm(x))
  
  logitModel  = glm(y~x, family=binomial(link="logit"))
  probitModel = glm(y~x, family=binomial(link="probit"))
  
  probLower[i] = deviance(probitModel)<deviance(logitModel)
  devLogit[i]  = deviance(logitModel)
  devProbit[i] = deviance(probitModel)
}

sum(probLower)/1000 

plot(devProbit,devLogit)

################################################################################################################

#generate data from logit
set.seed(1)
probLower = vector(length=1000)
devLogit = vector(length=1000)
devProbit = vector(length=1000)

n <- 1000
beta0 <- 0
beta1 <- 1
pi_x <- exp(beta0 + beta1 * x) / (1 + exp(beta0 + beta1 * x)) 

for(i in 1:1000){      
  x <- runif(n=n, min=-5, max=5)
  y <- rbinom(n=length(x), size=1, prob=pi_x)
  
  logitModel  = glm(y~x, family=binomial(link="logit"))
  probitModel = glm(y~x, family=binomial(link="probit"))
  
  probLower[i] = deviance(logitModel)<deviance(probitModel)
  devLogit[i]  = deviance(logitModel)
  devProbit[i] = deviance(probitModel)
}

sum(probLower)/1000 

plot(devProbit,devLogit)
