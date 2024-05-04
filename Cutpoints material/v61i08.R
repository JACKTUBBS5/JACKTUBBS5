## options
options(prompt = "R> ", width = 76, digits = 6)

## package
library("OptimalCutpoints")

## data
data("elas")
summary(elas)

###############
# Section 2.1 #
###############
                                                           
# "MinValueSe" criterion with Se >= 0.95:
cutpoint <- optimal.cutpoints(X = "elas", status = "status", tag.healthy = 0, methods = "MinValueSe", data = elas, 
  pop.prev = NULL, control = control.cutpoints(valueSe=0.95), ci.fit = TRUE)
summary(cutpoint)

# "SpEqualSe" method:
cutpoint <- optimal.cutpoints(X = "elas", status = "status", tag.healthy = 0, methods = "SpEqualSe", data = elas, 
  pop.prev = NULL, control = control.cutpoints(), ci.fit = TRUE)
summary(cutpoint)

# "Youden" method:
cutpoint <- optimal.cutpoints(X = "elas", status = "status", tag.healthy = 0, methods = "Youden", data = elas, 
  pop.prev = NULL, control = control.cutpoints(), ci.fit = TRUE)
summary(cutpoint)

# "MinValueSp" method with Sp >= 0.95::
cutpoint <- optimal.cutpoints(X = "elas", status = "status", tag.healthy = 0, methods = "MinValueSp", data = elas,
  pop.prev = NULL, control = control.cutpoints(valueSp=0.95), ci.fit = TRUE)
summary(cutpoint)

# "MinValueNPV" method with NPV >= 0.95:
cutpoint <- optimal.cutpoints(X = "elas", status = "status", tag.healthy = 0, methods = "MinValueNPV", data = elas, 
  pop.prev = NULL, control = control.cutpoints(valueNPV=0.95), ci.fit = TRUE)
summary(cutpoint)

# "ValueDLR.Positive" method with a DLR+ = 2:
cutpoint <- optimal.cutpoints(X = "elas", status = "status", tag.healthy = 0, methods = "ValueDLR.Positive", data = elas, 
  pop.prev = NULL, control = control.cutpoints(), ci.fit = TRUE)
summary(cutpoint)


###############
# Section 2.2 #
###############

# "MCT" method with CFP = 1, CFN = 3  
cutpoint <- optimal.cutpoints(X = "elas", status = "status", tag.healthy = 0, methods = "MCT", data = elas, 
  pop.prev = NULL, control = control.cutpoints(CFN=3), ci.fit = TRUE)
summary(cutpoint)


###############
## Section 4 ##
###############

## first example
cutpoint1 <- optimal.cutpoints(X = "elas", status = "status",
  tag.healthy = 0, methods = c("Youden", "SpEqualSe"), data = elas,
  categorical.cov = "gender", pop.prev = NULL,
  control = control.cutpoints(), ci.fit = TRUE)

names(cutpoint1)  
names(cutpoint1$Youden$Male)
names(cutpoint1$Youden$Male$measures.acc)
names(cutpoint1$Youden$Male$optimal.cutoff)

cutpoint1$Youden$Male$criterion
cutpoint1$Youden$Male$optimal.criterion
cutpoint1$Youden$Male$measures.acc$Se
cutpoint1$SpEqualSe$Male$measures.acc$Se


#################
## Section 4.1 ##
#################

options(digits = 7)
summary(cutpoint1)

cutpoint2 <- optimal.cutpoints(X = "elas", status = "status",
  tag.healthy = 0, methods = c("Youden","SpEqualSe"), data = elas,
  pop.prev = NULL, categorical.cov = "gender",
  control = control.cutpoints(ci.SeSp = "AgrestiCoull"), ci.fit = TRUE)
summary(cutpoint2)

cutpoint3 <- optimal.cutpoints(X = "elas", status = "status",
  tag.healthy = 0, methods = c("Youden","SpEqualSe"), data = elas,
  pop.prev = NULL, categorical.cov = "gender",
  control = control.cutpoints(costs.benefits.Youden = TRUE),
  ci.fit = TRUE)
summary(cutpoint3)

cutpoint4 <- optimal.cutpoints(X = "elas", status = "status",
  tag.healthy = 0, methods = c("Youden", "SpEqualSe"), data = elas,
  pop.prev = NULL, categorical.cov = "gender", control =
    control.cutpoints(generalized.Youden = TRUE, CFP = 1, CFN = 3),
  ci.fit = TRUE) 
summary(cutpoint4)

cutpoint5 <- optimal.cutpoints(X = "elas", status = "status",
  tag.healthy = 0, methods = c("Youden", "SpEqualSe"), data = elas,
  pop.prev = 0.5, categorical.cov = "gender",
  control = control.cutpoints(generalized.Youden = TRUE),
  ci.fit = TRUE)
summary(cutpoint5)

cutpoint5 <- optimal.cutpoints(X = "elas", status = "status",
  tag.healthy = 0, methods = c("Youden", "SpEqualSe"), data = elas,
  pop.prev = c(0.5, 0.5), categorical.cov = "gender",
  control = control.cutpoints(generalized.Youden = TRUE),
  ci.fit = TRUE) 
summary(cutpoint5)


#################
## Section 4.1 ##
#################

plot(cutpoint1)
plot(cutpoint1, which = 3, ylim = c(0, 1))

# In males:
plot(cutpoint1$SpEqualSe$Male$measures.acc$cutoffs,
  cutpoint1$SpEqualSe$Male$measures.acc$Se[,1], xlab = "Cutpoint",
  ylab = "sensitivity and specificity", type = "l",lty = 2,
  main = "sensitivity and specificity \n Males")

lines(cutpoint1$SpEqualSe$Male$measures.acc$cutoffs,
  cutpoint1$SpEqualSe$Male$measures.acc$Sp[,1], xlab = "Cutpoint",
  ylab = "sensitivity and specificity", type = "l")
legend("topright", legend = c("Se","Sp"), lty = c(2,1), bty = "n")

# In females:
plot(cutpoint1$SpEqualSe$Female$measures.acc$cutoffs,
  cutpoint1$SpEqualSe$Female$measures.acc$Se[,1], xlab = "Cutpoint",
  ylab = "sensitivity and specificity", type = "l",lty = 2,
  main = "sensitivity and specificity \n Females")

lines(cutpoint1$SpEqualSe$Female$measures.acc$cutoffs,
  cutpoint1$SpEqualSe$Female$measures.acc$Sp[,1], xlab = "Cutpoint",
  ylab = "sensitivity and specificity", type = "l")
legend("topright", legend = c("Se","Sp"), lty = c(2,1), bty = "n")
