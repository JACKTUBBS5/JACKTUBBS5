ROCregAUCci <-
function(
   method,
   model,
   group,
   tag.healthy,
   marker,
   data,
   control = controlROCreg(),
   newdata = NULL,
   conf = 0.95) {

   bootstrap.sample <- function(data, group, method) {
      if(method=="coutcome")
      res<-do.call("rbind",lapply(split(data,data[,group]), 
         function(x)x[sample(nrow(x), replace=TRUE),]))
      else
      res<-data[sample(nrow(data), replace=TRUE),]
    res
   }
   if (method=="NM") {
fit.ROC<-ROCreg(method = method, model=model, group=group, tag.healthy=tag.healthy, marker=marker, data=data, control=control, se.fit=TRUE)
    if(is.null(newdata))
      newdata <- ROCregData(fit.ROC)
    if(sum(is.na(match(fit.ROC$names.cov$g, names(newdata)))))
      stop("Not all needed variables are supplied in newdata")
      
   AUC <-predict(fit.ROC,accuracy="AUC",newdata=newdata)$AUC
    m.h <- model.matrix(delete.response(fit.ROC$ROC.fit$terms$h), model.frame(delete.response(fit.ROC$ROC.fit$terms$h), newdata))
    m.d <- model.matrix(delete.response(fit.ROC$ROC.fit$terms$d), model.frame(delete.response(fit.ROC$ROC.fit$terms$d), newdata))
    x.h <-  m.h%*% fit.ROC$ROC.fit$beta$h
    x.d <-  m.d%*% fit.ROC$ROC.fit$beta$d
    delta<- (x.d - x.h)/ (sqrt(fit.ROC$ROC.fit$sd$d^2 + fit.ROC$ROC.fit$sd$h^2))
    f <- 2*((fit.ROC$ROC.fit$sd$h^2 + fit.ROC$ROC.fit$sd$d^2)/(fit.ROC$ROC.fit$sd$h^4/(fit.ROC$ROC.fit$df$h) + fit.ROC$ROC.fit$sd$d^4/(fit.ROC$ROC.fit$df$d)))
    a1 <- diag((m.h%*%fit.ROC$ROC.fit$AUCic$h%*%t(m.h))*fit.ROC$ROC.fit$sd$h^2)
    a2 <- diag((m.d%*%fit.ROC$ROC.fit$AUCic$d%*%t(m.d))*fit.ROC$ROC.fit$sd$d^2)
    M <- ((a1+a2)^(-1))*(fit.ROC$ROC.fit$sd$h^2+fit.ROC$ROC.fit$sd$d^2)
    AUC.ll<-pnorm(delta-sqrt((M^(-1)+((delta^2)/f)))*qnorm(1-(1-conf)/2))
    AUC.ul<-pnorm(delta+sqrt((M^(-1)+((delta^2)/f)))*qnorm(1-(1-conf)/2))
  } else {
      fit.ROC<-ROCreg(method = method, model=model, group=group, tag.healthy=tag.healthy, marker=marker, data=data, control=control, se.fit=TRUE)
      if(is.null(newdata))
          newdata <- ROCregData(fit.ROC)
      if(sum(is.na(match(fit.ROC$names.cov$g, names(newdata)))))
          stop("Not all needed variables are supplied in newdata")
      AUC <-predict(fit.ROC, accuracy="AUC", newdata=newdata)$AUC      
      AUC.b <- array(dim=c(dim(newdata)[1], control$nboot), dimnames=list(1:dim(newdata)[1],c(1:(control$nboot))))
      for (i in 1:control$nboot) {
         data.boot <- bootstrap.sample(data, group, control$resample.m)
         fit.ROC <- ROCreg(method=method, model=model, group=group, tag.healthy=tag.healthy, marker=marker, data=data.boot, control=control)
         AUC.b[,i]<-predict(fit.ROC, accuracy="AUC", newdata=newdata)$AUC
      }
      AUC.ll <- apply(AUC.b,1,function(x) quantile(x,(1-conf)/2))
      AUC.ul <- apply(AUC.b,1,function(x) quantile(x,1-(1-conf)/2))
}
   res<-list(data=newdata, conf=conf, AUC = AUC, AUC.ll = AUC.ll, AUC.ul = AUC.ul)
   invisible(res)
}

