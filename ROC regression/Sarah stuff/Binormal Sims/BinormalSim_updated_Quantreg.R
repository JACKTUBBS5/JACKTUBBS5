#  Binormal Simulation -- parametric (Alonzo), Beta

library(quantreg)
library(foreign)
library(betareg)



set.seed(123)
z <- runif(200,0,1)

  # normal data generation
  yDis <- sapply(z, function(x) rnorm(1, 2 + 4*x, (1.5)))        # generating diseased points
  yRef <- sapply(z, function(x) rnorm(1, 1.5 + 3*x,( 1.5)))      # generating reference points
  
  yInd <- rep(c(1,0), each = length(z))  # disease indicator -- 1 = diseased, 0 = reference
  yResp <- c(yDis, yRef)
  
  data421 <- data.frame(cbind(yInd, yResp, z ))
  
  binorm.ref <- data.frame(cbind("y" = yRef, "x" = z))  # response and covariate for reference
  binorm.dis <- data.frame(cbind("y" = yDis,"x" = z))   # response and covariate for diseased
  
  # specifying set of false positive rates
  FPR <- seq(.01960784, .980392, by = .01960784)  #(1:99/100)  #--- original FPR set
  
  # quantile regression to estimate reference survival
  qr1 <- rq(y ~ x, data=binorm.ref, tau = rev(FPR))#data421[101:200,], tau = rev(FPR) )
  
  #covariate adjusted survival for diseased observations
  pred1 <- predict.rq(qr1, newdata = binorm.dis)
  
  Inv.t <- qnorm(FPR)   # calculating inverse normal of FPRs
  nq <- length(FPR)     # number of quantiles
  nd <- nrow(binorm.dis)  # number of diseased points
  
  
  # reshaping the data
  trans.pred1 <- t(pred1)          # transposing to obtain 50 X 200 vector
  col.pred1 <- c(trans.pred1)      # concatenating into a single column vector
  col.t <- rep(Inv.t, nd)          # creating long vector of Inv.t  
  col.ydis <- rep(binorm.dis$y, each = 50)
  col.x <- rep(binorm.dis$x, each = 50) 
  col.uit <- as.numeric(col.ydis >= col.pred1)  
  
  probitData <- data.frame( "fdbar" = col.pred1, "phiInv" =  col.t, "disRes" = col.ydis, "covX" = col.x, "uit" = col.uit)
  
  probitMod <- glm(uit ~ phiInv + covX, family = binomial(link = "probit"), data = probitData)
  summary(probitMod)
  
  
  #####  coefficients from probit model used in calculation of ROC
  
  
  
  ################################################
  #   Beta Method
  ################################################
  
  # normal data generation -- see above in Alonzo section
  #############################################################
  # calculating placement values for each diseased point
  
  pv <- rep(.0001, nd)  # if diseased point larger than predicted, default to .9999
  #pred1means2 <- rev(pred1means)
  
  # If working as intended, pv[1] = .6862746 for this example because binorm.dis$y[1] = 2.0847
  for (i in 1:nd){
    for(j in 1:(nq - 1)){
      if( binorm.dis$y[i] > rev(pred1[i, ])[j] && binorm.dis$y[i] <= rev(pred1[i,])[j+1]){
        pv[i] <- 1 - FPR[j]
      }
      if(binorm.dis$y[i] < min(pred1)){
        pv[i] <- .9999
      }
      
    } #end j 
  } # end i
  
  ################################ using B.hat instead of dat (as of 8/11/2016)
  
  
  temp <-data.frame(cbind(binorm.dis, pv))
  
  BetaModel <- betareg(pv ~ x, data = temp, link.phi = "identity", link = "logit")
  summary(BetaModel)
  
  # coefficients from beta model used in calculation of the ROC
  
  