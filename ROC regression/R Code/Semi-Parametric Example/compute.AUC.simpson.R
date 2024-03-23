compute.AUC.simpson <-
function(data.ROC, set.p, ROCreg) {
   simpson <- function(ROC, set.p) {
  l.set.p <- length(set.p)
  integral <- (set.p[l.set.p] - set.p[1])/(l.set.p - 1)/3*(ROC[1] + ROC[l.set.p] + 4*sum(ROC[seq(2,l.set.p - 1, by = 2)]) + 2*sum(ROC[seq(3, l.set.p - 2, by = 2)]))
   }
      
   compute.AUC.NM <- function(data, covar, set.p, ROC.fit) { 
      x.h <- model.matrix(delete.response(ROC.fit$terms$h), model.frame(delete.response(ROC.fit$terms$h), covar)) %*% ROC.fit$beta$h
      x.d <- model.matrix(delete.response(ROC.fit$terms$d), model.frame(delete.response(ROC.fit$terms$d), covar)) %*% ROC.fit$beta$d 
      res<-as.vector(pnorm((x.d - x.h)/(sqrt(ROC.fit$sd$d^2 + ROC.fit$sd$h^2))))
      res      
   }
   
   compute.AUC.GLM <- function(data, set.p)
      apply(data, 1, simpson, set.p = set.p)


   switch(ROCreg$method,
      "NM" = compute.AUC.NM(data = data.ROC$ROC, covar = data.ROC$data, set.p = set.p, ROC.fit = ROCreg$ROC.fit),
      "SM"=,
"PROCGLM" = ,
"SROCGLM" = compute.AUC.GLM(data = data.ROC$ROC, set.p = set.p))
}

