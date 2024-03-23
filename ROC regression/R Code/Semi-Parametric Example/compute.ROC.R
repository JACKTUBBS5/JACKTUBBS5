compute.ROC <-
function(data.ROC, set.p, ROCreg) {

emp.surv.function<-function(sequence,values) {
  sapply(sequence,function(s,values)mean(values>=s),values=values)
 }        
quantile.surv.function<-function(sequence,values) {
 surv<-emp.surv.function(values,values)
  sapply(sequence,function(s,surv){min(values[surv<=s])},surv=surv)
}
 compute.NM <- function(data.ROC, set.p, ROC.fit) {
 binormal <- function(x.d, x.h, se.d, se.h, p) pnorm((x.d - x.h)/se.d + se.h/se.d*qnorm(p))
 
preROC.h <- model.matrix(delete.response(ROC.fit$terms$h), model.frame(delete.response(ROC.fit$terms$h),data.ROC))%*%ROC.fit$beta$h
preROC.d <- model.matrix(delete.response(ROC.fit$terms$d), model.frame(delete.response(ROC.fit$terms$d),data.ROC))%*%ROC.fit$beta$d

ROC <- sapply(set.p, function(s) binormal(preROC.d, preROC.h, ROC.fit$sd$d, ROC.fit$sd$h, s))
if(is.vector(ROC)) ROC <- matrix(ROC, nrow=1)
if(length(ROC) > 1) colnames(ROC) <- paste("p", round(set.p, 3), sep = "")
ROC
 }
compute.SM <- function(data.ROC, set.p, ROC.fit) { 
preROC.h <- model.matrix(delete.response(ROC.fit$terms$h), model.frame(delete.response(ROC.fit$terms$h),data.ROC))%*%ROC.fit$beta$h
preROC.d <- model.matrix(delete.response(ROC.fit$terms$d), model.frame(delete.response(ROC.fit$terms$d),data.ROC))%*%ROC.fit$beta$d
quantiles<-quantile.surv.function(set.p,ROC.fit$residuals$h/ROC.fit$sd$h)
  ROC<-sapply(quantiles, function(q) emp.surv.function((preROC.h-preROC.d)/ROC.fit$sd$d+(ROC.fit$sd$h/ROC.fit$sd$d)*q, ROC.fit$residuals$d/ROC.fit$sd$d)) 
if(is.vector(ROC)) ROC <- matrix(ROC, nrow=1)
if(length(ROC) > 1) colnames(ROC) <- paste("p", round(set.p, 3), sep = "")
ROC
}
compute.PROCGLM <- function(data.ROC, ROC.model, set.p, ROC.fit) {
    q.f <- switch(ROC.model, "binormal" = "qnorm", "logistic" = "qlogis")
link <- switch(ROC.model, "binormal" = "pnorm", "logistic" = "qnorm" )
h <- eval(parse(text = q.f))(set.p)
   ROC <- sapply(1:length(h), function(s) {
      eval(parse(text = link))(model.matrix(delete.response(ROC.fit$terms), model.frame(delete.response(ROC.fit$terms),cbind(h = h[s], data.ROC)))%*%ROC.fit$beta.ROC)})
   if(is.vector(ROC)) ROC<-matrix(ROC,nrow=1)           
if(length(ROC)>1) colnames(ROC) <- paste("p", round(set.p, 3), sep = "")      
    if (match("p0", paste("p", round(set.p, 20), sep = ""), -1) != -1) ROC[,"p0"] <- 0
    if (match("p1", paste("p", round(set.p, 20), sep = ""), -1) != -1) ROC[,"p1"] <- 1
    ROC                                    
}
compute.SROCGLM <- function(data.ROC, set.p, ROC.fit) {
  X <- model.matrix(ROC.fit$terms, model.frame(ROC.fit$terms, data.ROC))
    h <- approxfun(ROC.fit$p.h, ROC.fit$h)(set.p)
   ROC <- sapply(1:length(set.p), function(s) pnorm(h[s] + X[, -1, drop = F] %*% ROC.fit$beta.ROC))
   if(is.vector(ROC)) ROC <- matrix(ROC, nrow = 1)
if(length(ROC) > 1) colnames(ROC) <- paste("p", round(set.p, 3), sep = "")
ROC
}   
  ROC <- switch(ROCreg$method,
"NM" = compute.NM(data.ROC = data.ROC, set.p = set.p, ROC.fit = ROCreg$ROC.fit),
"SM" = compute.SM(data.ROC = data.ROC, set.p = set.p, ROC.fit = ROCreg$ROC.fit),
    "PROCGLM" = compute.PROCGLM(data.ROC = data.ROC, ROCreg$control$ROC.model, set.p = set.p, ROC.fit = ROCreg$ROC.fit),
 "SROCGLM" = compute.SROCGLM(data.ROC = data.ROC, set.p = set.p, ROC.fit = ROCreg$ROC.fit))
}

