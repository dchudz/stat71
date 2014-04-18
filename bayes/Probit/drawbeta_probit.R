# SIR for Probit
## (From Phil on April 8, 2014)

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

#?glm


drawbeta = function(N, X, y,int=T, prnt=T){
  # returns N independent posterior draws of beta for a probit regression
  # model P(Y=1) = pnorm(X%*%beta)
  # output is a pxN matrix with beta values as columns.
  #
  # beta is a p-vector of regression coefficients.
  # X is nxp matrix and includes a column of 1's if you want an intercept. 
  # y is an n-vector of 0's and 1's. 
  #
  # use glm to find MLE for beta:
  p = dim(X)[2]
  # check for intercept

  if(int)
    beta0 = glm(y~X[,-1], family=binomial(probit))$coef
  if(!int)
    beta0 = glm(y~ -1+X, family=binomial(probit))$coef 
 
  d2mat = d2lfbeta(beta0, X,y)
  # set inverse covariance matrix
  Vinv = -d2mat/(p+1)
  V = solve(Vinv)
  # find square root of the covariance matrix
  rtV = t(chol(V))
  #
  # make 10 times as many candidate draws as desired
  N0 = 10*N
  # generate a pxN0 matrix of candidate beta draws centered to have
  # mode beta0 and scaled to have 2nd derivative matrix Vinv at the mode:
  betavals1 = beta0 + rtV%*%matrix(rt(p*N0, 1),nrow=p)
  #
  # evaluate log likelihood and approximating mv-t density for each draw
  lf = lf0 = rep(0,N0)
  for(iter in 1:N0){
    if(prnt)
      if(iter/10000 == round(iter/10000)) cat(iter,", ")
    b = betavals1[,iter]
    lf[iter] = lfbeta(b,X,y)
    lf0[iter] = -((p+1)/2)*log(1+t(b-beta0)%*%Vinv%*%(b-beta0))
  }
  lr = lf-lf0
  r = exp(lr-max(lr))
  # drop values where r==0
  N0 = sum(r>0)
  betavals1 = betavals1[,r>0]
  r = r[r>0]
  indx = sample(1:N0, N, prob=r, replace=F)
  if(p==1)
    out = betavals1[indx]
  if(p > 1)
    out = betavals1[,indx]
  out
}

