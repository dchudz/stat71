data {
 int<lower=0> n; // sample size
 int<lower=0> m; // number of features
 int<lower=0,upper=1> Y[n]; // response variable
 matrix[n,m] X; // explanatory variable
}
parameters { 
 vector[m] beta;
}
model{
 vector[n] theta;
 theta <- X * beta;
 for (i in 1:n)
  Y[i] ~ bernoulli(inv_logit(theta[i]));
}