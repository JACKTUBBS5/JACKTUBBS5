compute.AUC.integrate <-
function(data.ROC, set.p, ROCreg) {

integrate.ROC <-function(x,ROCreg) {
obj <- function(t,x,ROCreg) {
unlist(compute.ROC(data.ROC=x, set.p=t, ROCreg=ROCreg))
}
integrate(obj, upper=1, lower=0, x=x, ROCreg=ROCreg)$value
}
compute.AUC.NM <- function(data, ROCreg) {    
      x.h <- model.matrix(delete.response(ROCreg$ROC.fit$terms$h), model.frame(delete.response(ROCreg$ROC.fit$terms$h), data)) %*% ROCreg$ROC.fit$beta$h
      x.d <- model.matrix(delete.response(ROCreg$ROC.fit$terms$d), model.frame(delete.response(ROCreg$ROC.fit$terms$d), data)) %*% ROCreg$ROC.fit$beta$d       
      res<-as.vector(pnorm((x.d - x.h)/(sqrt(ROCreg$ROC.fit$sd$d^2 + ROCreg$ROC.fit$sd$h^2))))
      res
   }

compute.AUC.GLM <- function(data, ROCreg)
sapply(1:dim(data)[1], 
function(i, data,ROCreg) {         
x <- data[i, , drop = FALSE]
integrate.ROC(x, ROCreg)
},data=data, ROCreg=ROCreg)
      
switch(ROCreg$method,
  "NM" = compute.AUC.NM(data = data.ROC$data, ROCreg=ROCreg),
     "SM"=,
  "PROCGLM" = ,
  "SROCGLM" = compute.AUC.GLM(data = data.ROC$data, ROCreg=ROCreg))
}

