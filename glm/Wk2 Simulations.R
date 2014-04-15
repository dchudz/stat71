sim = 1000				# How many simulations to run
n = 100				# How many points to estimate

betaList = seq(-2,2,.2)
result = matrix(NA, nrow=length(betaList)*5, ncol = 4,
dimnames = list(NULL, c("xrange","beta","betaLogit/betaProbit","likLogit>likProbit")))

row = 1
for (beta in betaList)
{

for (xrange in 1:5) 
{

# Simulate x values
# x = runif(n, -xrange, xrange)		#Symmetrically pick x
# x = runif(n, 0, xrange)			#Just positive x
x = sample(1:xrange, n, replace=TRUE)	#Only positive integers

betaLogit = betaProbit = likLogit = likProbit = rep(NA, sim)

for (i in 1:sim) { 
  u = runif(n)
  T = log(u/(1-u))
  #T = qnorm(u)
  theta = x*beta + T

  # use theta to find Y values
  y = rep(0, n)
  y[which(theta>0)] = 1

  fitLogit = glm(y ~ x-1, family = binomial(logit))
  betaLogit[i] = fitLogit$coefficients
  likLogit[i] = sum( y*log(1/(1+exp(-x*betaLogit[i]))) + 
                (1-y)*log(exp(-x*betaLogit[i])/(1+exp(-x*betaLogit[i]))))

  fitProbit = glm(y ~ x-1, family=binomial(probit))
  betaProbit[i]=fitProbit$coefficients
  likProbit[i] = sum( y*log(pnorm(x*betaProbit[i])) + 
                     (1-y)*log(pnorm(-x*betaProbit[i])) )

}

# Out of curiosity: check the relationship between Logit and Probit MLE estimate
fit = lm(betaLogit ~ betaProbit)

result[row, "xrange"] = xrange
result[row, "beta"] = beta
result[row, "betaLogit/betaProbit"] = fit$coefficients[2]
result[row, "likLogit>likProbit"] = length(which(likLogit>likProbit))
row = row+1
} #end of xrange loop
} #end of beta loop



# save the results in an Excel file just in case we need it
library(XLConnect)

# If simulation was based on Logit:
workbook = loadWorkbook("result2.xlsx", create=TRUE)
createSheet(workbook, "simulate Logit pos integer x")
writeWorksheet(workbook, result, 2, startRow=1, startCol=1, header=TRUE)
saveWorkbook(workbook)

pdf("Simulate with positive integer x Logit__hist likLog big.pdf")
hist(result[,"likLogit>likProbit"])
dev.off()


# If simulation was based on Probit:
workbook = loadWorkbook("result.xlsx", create=TRUE)
createSheet(workbook, "simulate Probit positive integer x")
writeWorksheet(workbook, result, 6, startRow=1, startCol=1, header=TRUE)
saveWorkbook(workbook)

pdf("Simulate with positive integer x Probit__hist likLog big.pdf")
hist(result[,"likLogit>likProbit"])
dev.off()




