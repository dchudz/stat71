library(rstan)
set.seed(987)
n <- 1000
x <- runif(n, min=-2, max=2)

nu.true <- 3
beta.true <- 0.5
p0 = pt(beta.true * x, nu.true)
y = rbinom(n, 1, p0)


stanFitProbit1D <- stan(model_code='
data {
 int<lower=0> n; // sample size
 int<lower=0,upper=1> Y[n]; // response variable
 real x[n]; // explanatory variable
 real nu;
}
parameters { 
 real beta;
}
model{
 for (i in 1:n)
  Y[i] ~ bernoulli(student_t_cdf(beta * x[i], nu, 0, 1));
}', data = list (Y=y, x=x,n=1000, nu=nu.true) , 
                        iter = 1000, chains = 4)

hist(extract(stanFitProbit1D)$beta)

