#Higher dimension for SIR Logit

lfbeta = function(beta,X,y){
# returns log-likelihood function of beta for logit regression model.
# beta is a p-vector of regression coefficients.
# X is nxp and includes a column of 1's if you want an intercept. 
# y is an n-vector of 0's and 1's. 
 mu = X %*% beta
 ltheta = ifelse(mu<200, mu-log(1+exp(mu)), 0)
 lmtheta = ifelse(mu<200, -log(1+exp(mu)), -mu)
 sum( y*ltheta + (1-y)*lmtheta)
 }

d1lfbeta = function(beta,X,y){
# evaluate first derivative of logit log-likelihood
 mu = c(X %*% beta)
 theta = exp(mu)/(1+exp(mu))
 #d1theta = theta - theta^2	#derivative of theta w.r.t mu
 #v = d1theta*(y-theta)/(theta*(1-theta))
 v = y - theta
 t(X) %*% v
 }


d2lfbeta = function(beta,X,y){
# evaluate pxp 2nd derivative matrix of logit log-likelihood
 mu = c(X %*% beta)
 theta = exp(mu)/(1+exp(mu))
 #d1theta = theta-theta^2
 #d2theta = 2*theta^3 - 3*theta^2 + theta	     #2nd derivative of theta wrt mu
 D = -diag(theta*(1-theta))
 t(X) %*% D %*% X
}



drawbeta = function(N, X, y,int=T, prnt=T){
# make N independent approximate posterior draws for beta based on a logit model
# model P(Y=1) = logit(X %*% beta)
# output is a pxN matrix with beta values as columns

# use glm to find MLE for beta:
p = dim(X)[2]
# check for intercept
if(int)
 beta0 = glm(y~X[,-1], family=binomial(logit))$coef
if(!int)
 beta0 = glm(y~ -1+X, family=binomial(logit))$coef
d2mat = d2lfbeta(beta0, X,y)

# Find inverse covariance matrix
Vinv = -d2mat/(p+1)
V = solve(Vinv)
rtV = t(chol(V))

# make 10 times as many candidate draws as requested.
N0 = 10*N
# generate a pxN0 matrix of candidate beta draws centered to have
# mode beta0 and scaled to have 2nd derivative matrix Vinv at the mode:
betavals1 = beta0 + rtV %*% matrix(rt(p*N0, 1),nrow=p)

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
out}


## p=2 example:

n=1000
x = seq(-2,2,length.out=n)
X = cbind(1,x)
p=2

# beta = c(0,1)
theta0 = exp(x)/(1+exp(x))
y = rbinom(n,1,theta0)

beta.sample = drawbeta(10000, X,y)
par(mfrow=c(2,1))
hist(beta.sample[1,])
hist(beta.sample[2,])
par(mfrow=c(1,1))

out = glm(y~X[,-1], family=binomial(logit))
bhat = out$coef

apply(beta.sample, 1, mean)
apply(beta.sample, 1, sd)

summary(out)

# make contour plots
beta0vals = seq(-.25,.25,length.out=30)
beta1vals = seq(0.75,1.25,length.out=30)
d2mat = d2lfbeta(bhat,X,y)
Vinv = -d2mat/(p+1)
V= solve(Vinv)

lfbmat = matrix(0,ncol=30,nrow=30)
lf0mat = lfbmat
for(i in 1:30)
 for(j in 1:30){
  b = c(beta0vals[i],beta1vals[j])
  mu = beta0vals[i] + beta1vals[j]*x
  theta = exp(mu)/(1+exp(mu))
  lfbmat[i,j] = sum(y*log(theta) + (1-y)*log(1-theta))
  lf0mat[i,j] = -((p+1)/2)*log(1+t(b-bhat)%*%Vinv%*%(b-bhat))
  }

contour(beta0vals, beta1vals, exp(lfbmat-max(lfbmat)), xlab=expression(beta[0]),ylab=expression(beta[1]))
points(bhat[1],bhat[2], pch=19)
par(new=T)
contour(beta0vals, beta1vals, exp(lf0mat-max(lf0mat)),col="red")
points(bhat[1], bhat[2], pch=19)


## p=3 example

n=1000
p=3
x = seq(-2,2,length.out=n)
x2 = rnorm(n)

X = cbind(1,x,x2)

# beta = c(0,1,1)
theta0 = exp(x+x2)/(1+exp(x+x2))
y = rbinom(n,1,theta0)

out=glm(y~X[,-1], family=binomial(logit))
bhat = out$coef
summary(out)

beta.sample = drawbeta(10000, X,y)
apply(beta.sample,1,mean)
apply(beta.sample,1,sd)

par(mfrow=c(3,1))
hist(beta.sample[1,])
hist(beta.sample[2,])
hist(beta.sample[3,])


## Example with no intercept:
n=1000
x1 = runif(n, min=-2, max=2)
x2 = runif(n, min=-2, max=2)
X = cbind(x1,x2)
betaTrue <- c(.5,1)
theta0 = exp(X %*% betaTrue)/(1+exp(X %*% betaTrue))
y = rbinom(n,1,theta0)
beta.sample = drawbeta(10000, X,y, int=FALSE)
apply(beta.sample,1,mean)
apply(beta.sample,1,sd)

par(mfrow=c(2,1))
hist(beta.sample[1,])
hist(beta.sample[2,])








