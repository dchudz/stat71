# Rejection Sampling for Probit
# (From Phil on April 11, 2014)

lfbeta = function(beta, X,y){
# returns log-likelihood function of beta for probit regression model.
# beta is a p-vector of regression coefficients.
# X is nxp and includes a column of 1's if you want an intercept. 
# y is an n-vector of 0's and 1's. 
 mu = X%*%beta
 ltheta = pnorm(mu,log=T)
 lmtheta = pnorm(mu, log=T, lower = F)
 sum(y*ltheta + (1-y)*lmtheta)
 }


d1lfbeta = function(beta,X,y){
# returns 1st derivative of the log-likelihood function for
# the probit regression model. 
# beta is a p-vector of regression coefficients.
# X is nxp and includes a column of 1's if you want an intercept. 
# y is an n-vector of 0's and 1's. 
 mu = c(X%*%beta)
 theta = pnorm(mu)
 v = dnorm(mu)*(y-theta)/(theta*(1-theta))
 t(X)%*%v
 }


d2lfbeta = function(beta,X,y){
# returns pxp 2nd derivative matrix of the log-likelihood function for
# the probit regression model. 
# beta is a p-vector of regression coefficients.
# X is nxp and includes a column of 1's if you want an intercept. 
# y is an n-vector of 0's and 1's. 
 mu = c(X%*%beta)
 theta = pnorm(mu)
 D = diag(dnorm(mu)^2*( y/theta^2 + (1-y)/(1-theta)^2) + 
   dnorm(mu)*mu*(y-theta)/(theta*(1-theta)))
#
-t(X)%*%D%*%X
 }


lmvt = function(t,nu=1, mu=NA, Siginv=NA){
 p=length(t)
 if(missing(mu)) mu = rep(0,p)
 if(missing(Siginv)) Siginv=diag(rep(1,p))
  -c((nu+p)*log(1+t(t-mu)%*%Siginv%*%(t-mu)/nu)/2)
 }

d1lmvt =  function(t,nu=1, mu=NA, Siginv=NA){
 p=length(t)
 if(missing(mu)) mu = rep(0,p)
 if(missing(Siginv)) Siginv=diag(rep(1,p))
 u = c(1+t(t-mu)%*%Siginv%*%(t-mu)/nu)
 -Siginv%*%(t-mu)*c( (nu+p)/(nu*u))
}


d2lmvt =  function(t,nu=1, mu=NA, Siginv=NA){
 p=length(t)
 if(missing(mu)) mu = rep(0,p)
 if(missing(Siginv)) Siginv=diag(rep(1,p))
 u = c(1+t(t-mu)%*%Siginv%*%(t-mu)/nu)
 v = Siginv%*%(t-mu)
 -Siginv*c( (nu+p)/(nu*u)) + c( 2*(nu+p)/(nu^2*u^2))*v%*%t(v)
}



probit.rej = function(N, X, y, prnt=T){
# returns N independent posterior draws of beta for a probit regression
# model P(Y=1) = pnorm(X%*%beta)
# output is a pxN matrix with beta values as columns.
#
# beta is a p-vector of regression coefficients.
# X is nxp matrix and includes a column of 1's if you want an intercept. 
# y is an n-vector of 0's and 1's. 
#
# note the dimension of beta:
 p = dim(X)[2]
# use glm to find MLE for beta:
 bhat = glm.fit(X,y, family=binomial(probit))$coef
 d2mat = d2lfbeta(bhat, X,y)
# set inverse covariance matrix:
 Vinv = -d2mat/(p+1)
 V = solve(Vinv)
# find matrix square root of the covariance matrix:
 rtV = t(chol(V))
#
# Use Newton-Raphson from bhat +/- rtV%*%1 to locate maximum density ratio
 maxlr = lfbeta(bhat,X,y)- lmvt(bhat, mu=bhat, Siginv=Vinv)
 incr = rtV%*%rep(1,p)
#
# first start from bhat + rtV%*%1
#
 beta0 = bhat+incr
 lr0 = lfbeta(beta0,X,y)- lmvt(beta0, mu=bhat, Siginv=Vinv)
 lr1 = lr0-1
 cnt=1
 while(cnt==1){
  d1lr = d1lfbeta(beta0,X,y) - d1lmvt(beta0, mu=bhat, Siginv=Vinv)
  d2lr = d2lfbeta(beta0,X,y) - d2lmvt(beta0, mu=bhat, Siginv=Vinv)
# take a Newton-Raphson step
  beta0 = beta0 - solve(d2lr,d1lr)
  lr1 = lfbeta(beta0,X,y) - lmvt(beta0, mu=bhat, Siginv=Vinv)
# check for convergence:
  if(lr1-lr0 < .0000001) 
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
 lr0 = lfbeta(beta0,X,y)- lmvt(beta0, mu=bhat, Siginv=Vinv)
 lr1 = lr0-1
 cnt=1
 while(cnt==1){
  d1lr = d1lfbeta(beta0,X,y) - d1lmvt(beta0, mu=bhat, Siginv=Vinv)
  d2lr = d2lfbeta(beta0,X,y) - d2lmvt(beta0, mu=bhat, Siginv=Vinv)
# take a Newton-Raphson step
  beta0 = beta0 - solve(d2lr,d1lr)
  lr1 = lfbeta(beta0,X,y) - lmvt(beta0, mu=bhat, Siginv=Vinv)
  if(lr1-lr0 < .0000001) 
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
  betavals = bhat + rtV%*%matrix(rt(p*Nrem, 1),nrow=p)
# evaluate log likelihood and approximating mv-t density for each draw
  lf = lf0 = rep(0,Nrem)
  for(iter in 1:Nrem){
   b = betavals[,iter]
   lf[iter] = lfbeta(b,X,y)
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
