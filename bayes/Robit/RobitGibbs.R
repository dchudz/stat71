
##################################
set.seed(021309)
library(mvtnorm)




rmvn=function(mu,Sigma){
     mu+crossprod(rnorm(length(mu)),chol(Sigma))
}




rcnorm = function(c){
     # generate n standard Normal deviates constrained to be greater than the elements
     # of the n-vector c. 
     n=length(c)
     z = rep(0, n)
     indx0 = c(1:n)[c < -2]
     c0 = c[c< -2]
     # identify cases to use rnorm
     n0 = length(indx0)
     if(n0 > 0){
          z0=rep(0,n0)
          indx0sub=1:n0
          cnt=1
          while(cnt==1){
               z0[indx0sub] = rnorm(n0)
               indx0sub = indx0sub[z0[indx0sub] < c0[indx0sub] ]
               n0 = length(indx0sub)
               if(n0==0) cnt = 0
          }
          z[indx0] = z0
     }
     #
     indx1 = c(1:n)[c >= -2 & c < 1];   c1 = c[c >= -2 & c<1]
     # identify easy cases where we can use qnorm
     n1 = length(indx1)
     if(n1 > 0){
          z1 = qnorm(runif(n1, pnorm(c1),1) )
          z[indx1] = z1
     }
     #
     indx2 = c(1:n)[c >= 1];  c2 = c[indx2]
     # identify hard cases for rejection sampling
     n2 = length(indx2)
     if(n2 > 0){
          z2 = rep(0,n2)
          indx2sub = 1:n2
          lambda = (c2 + sqrt(c2^2 + 4))/2
          logM = dnorm(lambda,log=T) - pnorm(-c2,log=T) - dgamma(lambda-c2,1,lambda,log=T)
          cnt=1
          while(cnt==1){
               z2[indx2sub] = rgamma(n2,1,lambda)+c2[indx2sub]
               p = exp( dnorm(z2[indx2sub],log=T) - pnorm(-c2[indx2sub],log=T) - 
                             dgamma(z2[indx2sub] - c2[indx2sub],1,lambda[indx2sub],log=T) - logM[indx2sub] )
               u = runif(n2)
               indx2sub = indx2sub[u > p]
               n2 = length(indx2sub)
               if(n2==0) cnt=0
          }
          z[indx2] =z2
     }
     z
}







robit_MCMC <- function(nu){
      
      # Simulate data
      n<-1000			                        # Sample Size
      x<-sample(0:4,n,replace=T)	          # Covariates
      b0<- -2.25				                    # True Intercept
      b1<-  1.75				                    # True Value of Beta
      eta<-b0+b1*x 			                  # Linear Predictor
      p<-exp(eta)/(1+exp(eta))	            # Logistic Prob. Vector
      y<-rbinom(n,1,p)				              # Observations
      X<-matrix(c(rep(1,n),x), ncol=2)     # Design matrix
      k<-2				                    # Number of parms
      
      #Initial Values
      #Beta<-c(fit$coef*.634)
      Betavals <- c(-2,-1,0,1,2)
      lambda<-rep(1,n)		                # Weights
      #nuvals <- c(1:1000)                    # t df (8 ~ logistic) 
      nuval <- nu
      
      # Store Results
      nsim<-1000  			                    # Number of Iterations of Gibbs Sampler
      z<-rep(0,n)					                # Latent Normal Variable
      Betamat<-matrix(0,nrow=nsim,ncol=k)  # Store Results
      
      #plot(1:nsim,Betamat[,2]/.634, type="l",ylim = c(0,3.5), col="lightgreen")
      #abline(h=mean.beta[2],col="black")
      colors = c("red","purple","black","green","blue")
      col = 1
      
      nuMat <- matrix(,nrow=5,ncol=3)
      j = 1 
      ###################
      # GIBBS SAMPLER	#
      ###################
      for (beta in Betavals){
            Beta = c(-beta,beta)
            for (i in 2:nsim) {
                  
                  # UPDATE Z
                  muz<-X%*%Beta			# Update Mean of Z
                  #z[y==0]<-qnorm(runif(n,0,pnorm(0,muz,sqrt(1/lambda))),muz,sqrt(1/lambda))[y==0]
                  #z[y==1]<-qnorm(runif(n,pnorm(0,muz,sqrt(1/lambda)),1),muz,sqrt(1/lambda))[y==1]
                  z <- ifelse(y==1, muz+rcnorm(-muz), muz-rcnorm(muz))
                  
                  #UPDATE BETA
                  vbeta<-solve(crossprod(X*sqrt(lambda)))
                  betahat <- vbeta%*%(t(X*lambda)%*%z)
                  #Beta<-Betamat[i,]<-c(rmvnorm(1,betahat,vbeta),method="chol")
                  Beta <- Betamat[i,] <- c(rmvn(as.vector(betahat),vbeta))
                  
                  # UPDATE LAMBDA
                  lambda<-rgamma(n,(nu+1)/2,scale=2/(nu+(z-X%*%Beta)^2))
                  if (i%%100==0) print(i)
            }
            #points(1:nsim,Betamat[,2]/.634, type="l",col=colors[col])
            col = col+1
            
            #Results
            mean.beta = colMeans(Betamat/.634) # Correction factor is 1/.634
            #print(mean.beta)
            nuMat[j,]<- rbind(c(nu,mean.beta[1],mean.beta[2]))
            j = j + 1 
            
      }
  nuMat
}

################
# EXAMPLE DATA #
################

bignuMat <- matrix(0,nrow=25,ncol=3)
colnames(bignuMat) <- c("nu","betamat[,1]","betamat[,2]")
count = 1
for(i in 3:7){
      
      nui <- robit_MCMC(i)
      for(j in 1:5){
            bignuMat[count,] <- nui[j,]
            count = count + 1
      }

}

bignuMat






