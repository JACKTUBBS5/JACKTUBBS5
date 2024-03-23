# Semi-parametric code -- based on comparative paper
# Data generated from sas (Alvarez_sim1.sas) same as (cao-PV_sim2)
#setwd("C:/Users/Sarah/Desktop/Box Sync/Box Sync")

library(quantreg)
library(foreign)
library(gtools)
binorm <- read.csv("binorm.csv")  # generated so that we have 50 diseased and 50 non-diseased

summary(binorm)
FPR <- c((1:99)/100) 

binorm.ref <- subset(binorm, LA == 0)
binorm.dis <- subset(binorm, LA == 1)

# quantile regression
qr1 <- rq(y ~ x, data=binorm.ref, tau = FPR )
coef(qr1)
pred1 <- predict.rq(qr1)

pred1means <- colMeans(pred1)

##################################################
# Calculating B.hat
#################################################
nd <- nrow(binorm.dis)
nq <- length(FPR)
np <- nq
dat <- matrix(, nrow = nd, ncol = nq)

B.hat <- NULL

for (j in 1:nq){
  for(i in 1:nd){
    dat[i,j] <- (binorm.dis$y[i] >= pred1means[j])*1
    
  }
}

for(i in 1:nd){
  B.hat <- append(B.hat, dat[i,])
}

#######################################################
# calculating placement values
#######################################################
pv <- rep(.9999, nd)  # if diseased point larger than predicted, default to .9999


# If working as intended, pv[1] = .16 for this example
for (i in 1:np){
  for(j in 1:(nq - 1)){
    if( binorm.dis$y[i] > pred1means[j] && binorm.dis$y[i] <= pred1means[j+1]){
      pv[i] <- 1 - .01*( j)
    }
    if(binorm.dis$y[i] < pred1means[1]){
      pv[i] <- .0001
    }
    
  } #end j 
} # end i

temp <-data.frame(cbind(binorm.dis, pv))

# Calculating pairwise differences of placement values
# "combn" function calculates row 2 - row 1, so take the negative to get row 1 - row2; store in pvDiff
pvComb <- combn(pv,2)
pvDiff <- -combn(pv, 2, diff)

# Calculating covariate differences
xDiff <- -combn(temp$x, 2, diff)

BinPV <- (pvDiff <= 0)*1  # if the difference is less than or equal to 0, we have a 1

# verifying that this works how we think it does
pvComb[,1:10]
pvDiff[1:10]
BinPV[1:10]


# storing binary and covariate difference vectors in a new data frame 
temp2 <- data.frame(cbind(BinPV, xDiff))
temp3 <- data.frame(cbind(B.hat = B.hat, xPbeta = rep(beta.hat*binorm.dis$x, 99)))

# First call to GLM to estimate the Betas

logitMod1 <- glm(BinPV ~ xDiff - 1 , family = binomial(link = "logit"), data = temp2) # No intercept
beta.hat <- logitMod1$coefficients

# Second call to GLM to estimate h0(t)
logitMod2 <- glm(B.hat ~ offset(xPbeta), family = binomial(link = "logit"), data = temp3)

h.hat <- logitMod2$coefficients

