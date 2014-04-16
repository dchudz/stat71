data {
 int<lower=0> n; // sample size
 int<lower=0,upper=1> Y[n]; // response variable
 real x[n]; // explanatory variable
}
parameters { 
 real beta1;
}
model{
 for (i in 1:n)
  Y[i] ~ bernoulli(Phi(beta1 * x[i]));
}