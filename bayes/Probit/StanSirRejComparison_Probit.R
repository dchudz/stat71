# Compare Stan, SIR and Rejection Sampling for Probit

#{r results='hide'}
library(rstan)
library(ggplot2)
library(coda)
library(plyr)
source("StanSIRRej/drawbeta_logit.R")
source("StanSIRRej/Logit_rejection.R")
set.seed(908952)


## Comparing Stan and SIR in 1D Logit Example

#First we simulate some data (one predictor logit, no intercept):


nStanIterations <- 3000
n=1000
beta <- .5

# sample from the logit model via latent data
x=rnorm(n,0,1)
theta = x*beta
prob = exp(theta)/(1+exp(theta))
y = rbinom(n, 1, prob)


#Now fit via Stan and SIR:
s71logit_dat <- list (Y=y, x=x,n=1000) 
stanFitLogit1D <- stan("StanSIRRej/Logit1D.Stan", data = s71logit_dat, iter = nStanIterations, chains = 4)
stanFitLogit1D

betaSamplesStan <- extract(stanFitLogit1D)$beta1
betaSamplesSir <- drawbeta(5000, matrix(x), y, int=FALSE)
betaSamplesRej <- c(logit.rej(5000, matrix(x), y))		#could set prnt=F here

samplesDF <- rbind(data.frame(beta = betaSamplesStan, sampler = "STAN"),
                   data.frame(beta = betaSamplesSir, sampler = "SIR"),
			 data.frame(beta = betaSamplesRej, sampler = "Rejection"))


#Here's a histogram of the samples:


ggplot(samplesDF) + 
  geom_histogram(aes(x = beta, y = ..density..)) +
  facet_grid(sampler ~ .) +
  ggtitle("Histogram of Samples Under Each Method")


#Here's a plot of the quantiles. 

probs <- seq(0,1,by=.01)
quantilesDF <- ddply(samplesDF, "sampler", function(df) {
  return(data.frame(prob=probs, quantile=quantile(df$beta, probs=probs)))
  })

ggplot(quantilesDF) + 
  geom_point(aes(x=prob, y=quantile, color=sampler), alpha=.5, size=3) +
  ggtitle("Quantiles of the Samples Under Each Method")


## Logit with Intercept

#Create simulated data:

n=1000
x1 = rep(1,n) # intercept
x2 = runif(n, min=-2, max=2)
X <- cbind(x1,x2)
beta <- c(.5, 1)
theta = X %*% beta
prob = exp(theta)/(1+exp(theta))
y = rbinom(n, 1, prob)


#Now fit via Stan and SIR:


stanFitLogitMultiple <- stan("StanSIRRej/LogitMultipleRegression.Stan",
                              data = list (Y=y, X=X,n=n, m=2), 
                              iter = nStanIterations, chains = 4)

stanFitLogitMultiple
betaSamplesStan <- extract(stanFitLogitMultiple)$beta
beta1SamplesStan <- betaSamplesStan[,2]

betaSamplesSir <- drawbeta(5000, X, y)
beta1SamplesSir <- betaSamplesSir[2,]

betaSamplesRej <- c(logit.rej(5000, X, y)[2,])

samplesDF <- rbind(data.frame(beta1 = beta1SamplesStan, sampler = "STAN"),
                   data.frame(beta1 = beta1SamplesSir, sampler = "SIR"),
			 data.frame(beta1 = betaSamplesRej, sampler = "Rejection"))


#I'm only comparing the posteriors for the intercept parameter for now:


ggplot(samplesDF) + 
  geom_histogram(aes(x = beta1, y = ..density..)) +
  facet_grid(sampler ~ .) +
  ggtitle("Histogram of Samples Under Each Method")


#Here's a plot of the quantiles. Still pretty similar.

probs <- seq(0,1,by=.01)
quantilesDF <- ddply(samplesDF, "sampler", function(df) {
  return(data.frame(prob=probs, quantile=quantile(df$beta1, probs=probs)))
  })

ggplot(quantilesDF) + 
  geom_point(aes(x=prob, y=quantile, color=sampler), alpha=.5, size=3) +
  ggtitle("Quantiles of the Samples Under Each Method")

