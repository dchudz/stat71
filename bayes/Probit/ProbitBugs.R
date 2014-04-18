
#setting up the problem

n=1000
beta<-1
# set a grid of x and error values
x=rnorm(n,0,1)
error<-rnorm(n, 0, 1)
# calculate Z
z<-x*beta+error
#set Y
Y<-rep(0,n)
for (i in 1:n){
if(z[i]>0){
Y[i]<-1}
}


#BUGS CODE FOR STAT 71
library(arm)

data <- list ("Y", "x","n")
inits <- function() {list (beta1=runif(1,-10,10))}
parameters <- c("beta1")
stat71probit.sim <- bugs (data, inits, parameters, "Stat71ProbitModel2.bug", n.chains=3, n.iter=10000,debug=TRUE)


