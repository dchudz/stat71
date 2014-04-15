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
p = pt(.3*x, 5)  		# Can try other p's here
y = rbinom(1000, 1, p)

for(j in 1:250){
  beta = NA
  oldbeta = 0
  i = 0
  result = matrix(NA, nrow = 29, ncol = 3)
  
  for(nu in 2:30){
    repeat{
      p = pt(x*oldbeta, nu)
      #print(c(oldbeta, llik(x,y,p,nu,oldbeta)))
      newbeta = oldbeta - 
        lprime(x,y,p,nu,oldbeta)/ldoubleprime(x,y,p,nu,oldbeta)
      beta = c(beta, newbeta)
      if (abs(oldbeta-newbeta) <=.001) {
        i = i + 1
        #print(c(newbeta, llik(x,y,p,nu,newbeta)))
        result[i,1] = nu
        result[i,2] = newbeta
        result[i,3] = llik(x,y,p,nu,newbeta)
        break
      }
      oldbeta = newbeta
    }
  }

if(j == 1){
  plot(result[,1],result[,3])
} else {
  points(result[,1],result[,3])
}

}
 



