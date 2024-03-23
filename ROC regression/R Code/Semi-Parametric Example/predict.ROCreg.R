predict.ROCreg <-
function(object, newdata = NULL, step.p=0.02, accuracy = NULL, AUC.calc=c("simpson","integrate"), control = controlROCregData(), ...) {
   if(is.null(newdata))
      newdata <- ROCregData(object, control = control)
   if(sum(is.na(match(object$names.cov$g, names(newdata)))))
      stop("Not all needed variables are supplied in newdata")        
   p<-seq(0, 1, by = step.p)
   res<-list()
   
   res$data<-newdata
   res$ROC<-compute.ROC(newdata, p, object)
   set.accuracy <- c("AUC", "YI", "TH")
   ind.accuracy <- is.element(set.accuracy, accuracy)
   if (ind.accuracy[1]) {
   AUC.calc<-match.arg(AUC.calc)    
   compute.AUC<-switch(AUC.calc,"simpson"="compute.AUC.simpson","integrate"="compute.AUC.integrate")  
      res$AUC <- eval(parse(text = compute.AUC))(res, p, object)
   }
   if(any(ind.accuracy[2], ind.accuracy[3])) {
      aux <- compute.YI(newdata, p, object, TH = ind.accuracy[3])
if(!is.null(res)) {
   if(ind.accuracy[2]) res$YI <- aux$YI
   if(ind.accuracy[3]) res$TH <- aux$TH
}
   }
   res$step.p<-step.p
  res$control<-control
   class(res) <- c("predict.ROCreg")
   res
}

