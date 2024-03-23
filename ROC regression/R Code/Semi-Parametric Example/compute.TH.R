compute.TH <-
function(data.ROC, p.opt, ROCreg) {
 quantile.surv.function<-function(sequence,values) {    
  sapply(sequence,function(s,surv){min(values[surv<=s])},surv=sapply(values,function(s,values)mean(values>=s),values=values))
 }
compute.TH.NM <- function(data, p.opt, ROC.fit = ROC.fit) {
 x.h <- model.matrix(delete.response(ROC.fit$terms$h), model.frame(delete.response(ROC.fit$terms$h),data))%*%ROC.fit$beta$h
 qnorm(1 - p.opt, x.h, ROC.fit$sd$h)
 }
 compute.TH.SM <- function(data, p.opt, ROC.fit = ROC.fit) {         
  x.h <- model.matrix(delete.response(ROC.fit$terms$h), model.frame(delete.response(ROC.fit$terms$h),data))%*%ROC.fit$beta$h   
res <- as.vector(x.h + ROC.fit$sd$h*quantile.surv.function(p.opt,ROC.fit$residuals$h/ROC.fit$sd$h))
res       
}
compute.TH.GLM <- function(data, p.opt, est.surv , ROC.fit = ROC.fit) {            
  sd.h <- sqrt(sum((ROC.fit$fit.h$residuals)^2/ROC.fit$fit.h$df.residual))
   if(est.surv =="empirical") {
res<-predict(ROC.fit$fit.h, newdata = data) + sd.h*quantile.surv.function(p.opt,ROC.fit$fit.h$residuals/sd.h)      
} else if(est.surv =="normal") {    
 res<-predict(ROC.fit$fit.h, newdata = data) + sd.h*qnorm(1-p.opt)
}
res       
 }
  TH <- switch(ROCreg$method,
      "NM" = compute.TH.NM(data.ROC, p.opt, ROC.fit = ROCreg$ROC.fit),
      "SM" = compute.TH.SM(data.ROC, p.opt, ROC.fit = ROCreg$ROC.fit),
"PROCGLM" = ,
"SROCGLM" = compute.TH.GLM(data.ROC, p.opt, est.surv = ROCreg$control$est.surv, ROC.fit = ROCreg$ROC.fit))
}

