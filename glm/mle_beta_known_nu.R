# Derivative of t density function
fprime <- function(x, nu) {
  return (dt(x,nu)*(-(nu+1)*x/(nu+x^2)))
}

# Calculate log likelihood 
llik <- function(x,y,p,nu,beta) {
  return (sum( y*log(p) + (1-y)*log(1-p) ))
}

# First derivative of llik
lprime <- function(x, y, p, nu, beta ){
  return (sum(
              y/p*dt(x*beta,nu)*x + (1-y)/(1-p)*(-dt(x*beta,nu)*x)))
}

# Second derivative of llik
ldoubleprime <- function(x, y, p, nu, beta) {
  return (sum(
              y*( 1/p* fprime(x*beta,nu)*x^2 - dt(x*beta,nu)^2*x^2/(p^2) ) +
              (1-y)* (1/(1-p)*(-fprime(x*beta,nu)*x^2)- dt(x*beta,nu)^2*x^2/((1-p)^2))
              )
         )
}


# Simulate x, y values
x = c(rep(1,500), rep(2,500))
p = pt(.3*x, 5)			# Can try other p's here
y = rbinom(1000, 1, p)
nu = 5


beta = NA
oldbeta = 0
repeat {
  p = pt(x*oldbeta, nu)
  print(c(oldbeta, llik(x,y,p,nu,oldbeta)))
  newbeta = oldbeta - 
            lprime(x,y,p,nu,oldbeta)/ldoubleprime(x,y,p,nu,oldbeta)
  beta = c(beta, newbeta)
  if (abs(oldbeta-newbeta) <=.001) {
     print(c(newbeta, llik(x,y,p,nu,newbeta)))
     break}
  oldbeta = newbeta
}

llikmax = -100000000000
betahat = NA
for (beta in seq(0,1,.001)) {
  p = pt(x*beta,nu)
  #print(c(beta, llik(x,y,p,nu,beta)))
  if (beta==0) {plot(beta, llik(x,y,p,nu,beta))}
  points(beta, llik(x,y,p,nu,beta))
  if (llik(x,y,p,nu,beta)>llikmax) {
     llikmax = llik(x,y,p,nu,beta)
     betahat = beta}
}


