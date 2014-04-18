# Everson, Schofield, Chudzicki, Ning, Oh, Sahai and Wang


lfrobit = function(beta, nu, X,y){
# returns log-likelihood function of beta for robit regression model.
# beta is a p-vector of regression coefficients.
# X is nxp and includes a column of 1's if you want an intercept. 
# y is an n-vector of 0's and 1's. 
 mu = X%*%beta
 ltheta = pt(mu,nu, log=T)
 lmtheta = pt(mu,nu, log=T, lower = F)
 sum(y*ltheta + (1-y)*lmtheta)
 }


d1lfrobit = function(beta,nu,X,y){
# returns 1st derivative of the log-likelihood function for
# the probit regression model. 
# beta is a p-vector of regression coefficients.
# X is nxp and includes a column of 1's if you want an intercept. 
# y is an n-vector of 0's and 1's. 
 mu = c(X%*%beta)
 theta = pt(mu,nu)
 v = dt(mu,nu)*(y-theta)/(theta*(1-theta))
 t(X)%*%v
 }


d1t = function(t,nu){
 # work out 1st derivative for the t(nu) pdf at t
 lc = lgamma((nu+1)/2) - lgamma(nu/2) - log(pi*nu)/2
 -t*exp(lc + log((nu+1)/nu) - (nu+3)*log(1+t^2/nu)/2)
 }


d2lfrobit = function(beta,nu, X,y){
# returns pxp 2nd derivative matrix of the log-likelihood function for
# the robit regression model with df nu. 
# beta is a p-vector of regression coefficients.
# X is nxp and includes a column of 1's if you want an intercept. 
# y is an n-vector of 0's and 1's. 
 mu = c(X%*%beta)
 theta = pt(mu,nu)
 D = diag(dt(mu,nu)^2*( y/theta^2 + (1-y)/(1-theta)^2) - 
   d1t(mu,nu)*(y-theta)/(theta*(1-theta)))
#
-t(X)%*%D%*%X
 }


lmvt = function(t,nu=1, mu=NA, Siginv=NA){
# returns the log of the multi-variate t density
 p=length(t)
 if(missing(mu)) mu = rep(0,p)
 if(missing(Siginv)) Siginv=diag(rep(1,p))
  -c((nu+p)*log(1+t(t-mu)%*%Siginv%*%(t-mu)/nu)/2)
 }

d1lmvt =  function(t,nu=1, mu=NA, Siginv=NA){
# returns the first derivative vector of the log of the multi-variate t density
 p=length(t)
 if(missing(mu)) mu = rep(0,p)
 if(missing(Siginv)) Siginv=diag(rep(1,p))
 u = c(1+t(t-mu)%*%Siginv%*%(t-mu)/nu)
 -Siginv%*%(t-mu)*c( (nu+p)/(nu*u))
}


d2lmvt =  function(t,nu=1, mu=NA, Siginv=NA){
# returns the second derivative vector of the log of the multi-variate t density
 p=length(t)
 if(missing(mu)) mu = rep(0,p)
 if(missing(Siginv)) Siginv=diag(rep(1,p))
 u = c(1+t(t-mu)%*%Siginv%*%(t-mu)/nu)
 v = Siginv%*%(t-mu)
 -Siginv*c( (nu+p)/(nu*u)) + c( 2*(nu+p)/(nu^2*u^2))*v%*%t(v)
}




robit.rej = function(N, X, y, nu=5, prnt=T){
# returns N independent posterior draws of beta for a robit regression
# model P(Y=1) = pt(X%*%beta, nu)
# output is a pxN matrix with beta values as columns.
#
# beta is a p-vector of regression coefficients.
# X is nxp matrix and includes a column of 1's if you want an intercept. 
# y is an n-vector of 0's and 1's. 
#
# note the dimension of beta:
 p = dim(X)[2]
# use glm to find a starting estimate for the posterior mode of beta:
 bhat = glm.fit(X,y, family=binomial(probit))$coef
 lf = lfrobit(bhat,nu, X,y)
# Use Newton-Raphson to find the posterior mode:
 cnt=1
 while(cnt==1){
  d1lf = d1lfrobit(bhat,nu,X,y)
  d2mat = d2lfrobit(bhat,nu, X,y) 
  bhat = bhat - solve(d2mat,d1lf)
  lf1= lfrobit(bhat,nu, X,y)
  if(lf1-lf <= .0000001) cnt=0
  lf=lf1
#  if(prnt) print(lf)
  }
#
# set inverse covariance matrix:
 d2mat = d2lfrobit(bhat,nu, X,y) 
 Vinv = -d2mat/(p+1)
 V = solve(Vinv)
# find matrix square root of the covariance matrix:
 rtV = t(chol(V)) 
#
# Use Newton-Raphson from bhat +/- rtV%*%1 to locate maximum density ratio
 maxlr = lfrobit(bhat,nu, X,y)- lmvt(bhat, mu=bhat, Siginv=Vinv)
 incr = rtV%*%rep(1,p)
#
# first start from bhat + rtV%*%1
#
 beta0 = bhat+incr
 lr0 = lfrobit(beta0,nu,X,y)- lmvt(beta0, mu=bhat, Siginv=Vinv)
 lr1 = lr0-1
 cnt=1
 while(cnt==1){
  d1lr = d1lfrobit(beta0,nu,X,y) - d1lmvt(beta0, mu=bhat, Siginv=Vinv)
  d2lr = d2lfrobit(beta0,nu,X,y) - d2lmvt(beta0, mu=bhat, Siginv=Vinv)
# take a Newton-Raphson step
  beta0 = beta0 - solve(d2lr,d1lr)
  lr1 = lfrobit(beta0,nu,X,y) - lmvt(beta0, mu=bhat, Siginv=Vinv)
# check for convergence:
  if(abs(lr1-lr0) < .000000001) 
   cnt=0
  lr0=lr1
#  if(prnt) print(lr0)
  }
# check if log-ratio is greater than at bhat
 if(lr0 >= maxlr)
  maxlr = lr0
#
# now start at bhat-rtV%*%1
 beta0 = bhat-incr
 lr0 = lfrobit(beta0,nu,X,y)- lmvt(beta0, mu=bhat, Siginv=Vinv)
 lr1 = lr0-1
 cnt=1
 while(cnt==1){
  d1lr = d1lfrobit(beta0,nu,X,y) - d1lmvt(beta0, mu=bhat, Siginv=Vinv)
  d2lr = d2lfrobit(beta0,nu,X,y) - d2lmvt(beta0, mu=bhat, Siginv=Vinv)
# take a Newton-Raphson step
  beta0 = beta0 - solve(d2lr,d1lr)
  lr1 = lfrobit(beta0,nu,X,y) - lmvt(beta0, mu=bhat, Siginv=Vinv)
  if(abs(lr1-lr0) < .0000001) 
   cnt=0
  lr0=lr1
#  if(prnt) print(lr0)
  }
# check if this is the larger than the other value:
 if(lr0 >= maxlr)
  maxlr = lr0
# assume we have found the max ratio to be maxlr
# Create space to store N accepted beta draws
 betamat = matrix(0, nrow=p,ncol=N)
#
 Nacc = 0; Nrem = N
 while(Nacc < N){
# make Nrem candidate draws; repeat until Nrem = 0 and Nacc=N
  betavals = c(bhat) + rtV%*%matrix(rt(p*Nrem, 1),nrow=p)
# evaluate log likelihood and approximating mv-t density for each draw
  lf = lf0 = rep(0,Nrem)
  for(iter in 1:Nrem){
   b = betavals[,iter]
   lf[iter] = lfrobit(b,nu,X,y)
   lf0[iter] = -((p+1)/2)*log(1+t(b-bhat)%*%Vinv%*%(b-bhat))
   }
# compute acceptance probabilities for each candidate draw:
  r = exp(lf-lf0-maxlr)
# generate uniforms to decide accept/reject
  u = runif(Nrem)
 betavals = betavals[,u <= r]
# record the number accepted on this pass:
  nacc = sum(u<=r)
# store new accepted draws
  if(nacc>0)
   betamat[,(Nacc+1):(Nacc+nacc)] = betavals
# update counts
  Nacc = Nacc + nacc;  Nrem = Nrem - nacc
  if(prnt) print(paste(Nacc, "accepted"),quote=F)
  } # end while(Nacc < N).  We should now have N accepted draws.
 betamat
 }


### example data and "by hand" computations:

n=20
x = seq(-2,2,length.out=n)
X=cbind(1,x)
p=2
prnt =T
cnt=T

mu = x
theta = pt(mu, 5)
y = rbinom(20,1,theta)

nu=5
 p = dim(X)[2]
# use glm to find a starting estimate for the posterior mode of beta:
 bhat = glm.fit(X,y, family=binomial(probit))$coef
 lf = lfrobit(bhat,nu, X,y)
# Use Newton-Raphson to find the posterior mode:
 cnt=1
 while(cnt==1){
  d1lf = d1lfrobit(bhat,nu,X,y)
  d2mat = d2lfrobit(bhat,nu, X,y) 
  bhat = bhat - solve(d2mat,d1lf)
  lf1= lfrobit(bhat,nu, X,y)
  if(lf1-lf <= .0000001) cnt=0
  lf=lf1
  points(bhat[1],bhat[2],col="blue")
  print(lf)
  }
#
# set inverse covariance matrix:
 d2mat = d2lfrobit(bhat,nu, X,y) 
 Vinv = -d2mat/(p+1)
 V = solve(Vinv)
# find matrix square root of the covariance matrix:
 rtV = t(chol(V)) 


beta0vals = seq(-1,1,length.out=50)
beta1vals = seq(0,2.5,length.out=50)

lfbmat = lfb0mat =matrix(0,ncol=50,nrow=50)

bmax = c(0,0)
lmax=-999999999999
for(i in 1:50)
 for(j in 1:50){
  b = c(beta0vals[i],beta1vals[j])
  theta = pt(beta0vals[i] + beta1vals[j]*x, nu)
  lfbmat[i,j] = sum(y*log(theta) + (1-y)*log(1-theta))
  if(lfbmat[i,j]>lmax){
   bmax=b
   lmax=lfbmat[i,j]
   }
  }

contour(beta0vals, beta1vals, exp(lfbmat-max(lfbmat)))

## having fit the mode and 2nd derivative, evaluate approximating density

for(i in 1:50)
 for(j in 1:50){
  b = c(beta0vals[i],beta1vals[j])
  lfb0mat[i,j] = -1.5*log(1+t(b-bhat)%*%Vinv%*%(b-bhat))
  }

par(new=T)
contour(beta0vals, beta1vals, exp(lfb0mat-max(lfb0mat)),col="red")

# look at ratio
contour(beta0vals, beta1vals, exp(lfbmat-max(lfbmat)))
par(new=T)
contour(beta0vals, beta1vals, exp(lfbmat-lfb0mat-max(lfbmat-lfb0mat)),col="blue")


betamat = robit.rej(10000,X,y,5)
points(betamat[1,],betamat[2,],pch=".")
