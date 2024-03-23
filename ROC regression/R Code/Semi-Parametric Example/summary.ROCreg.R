summary.ROCreg <-
function(object, conf = 0.95, ndigits = 4, ...) {
    beta.ROC <- object$ROC.fit$beta.ROC
    se.beta <- object$ROC.fit$se.beta
    d <- object$ROC.fit$df$d
    h <- object$ROC.fit$df$h
    
    coeff <- formatC(beta.ROC, format = "f", digits = ndigits)
        
   ci.legend <- paste(conf*100, "% Conf. Interval", sep = "")
   
  if(object$se.fit) {
 if(object$method %in%c("SM","PROCGLM","SROCGLM")) {   
se <- formatC(se.beta, format = "f", digits = ndigits)
q.alpha <- qnorm((1-(1-conf)/2))
       ci <- paste("(", formatC(beta.ROC - q.alpha*se.beta, digits = ndigits, format = "f"), ",",
               formatC(beta.ROC + q.alpha*se.beta, format = "f", digits = ndigits),")", sep = "")
       pvalue <- formatC((1 - pnorm(abs(beta.ROC/se.beta)))*2,digits = ndigits, format = "f")
       pvalue <- replace(pvalue, pvalue == format(0, nsmall = ndigits), "<0.0001")
       ci <- as.matrix(ci, ncol = 1)
       colnames(ci) <- ci.legend
              
       p.table <- cbind(coeff, se, ci, pvalue)
            dimnames(p.table) <- list(names(beta.ROC), c("Estimate","Std. Error", ci.legend, "p-value"))
            
  } else
   if (object$method == "NM") {
   M.hat <- object$ROC.fit$M  
     reg.coefs <- beta.ROC[!(names(beta.ROC) == "h")]
     delta.M <- reg.coefs*M.hat   
     pvalue <- formatC((1 - pt(abs(delta.M), df = d))*2, digits = ndigits, format = "f")
     pvalue <- replace(pvalue, pvalue == format(0, nsmall = ndigits), "<0.0001")   
     ci <- vector()      
for(i in 1:length(reg.coefs)) {
conf.int <- conf.limits.nct(delta.M[i], df = d, conf)
ci[i] <- paste("(", formatC(conf.int$Lower.Limit/M.hat[i], digits = ndigits, format = "f"), ",",
            formatC(conf.int$Upper.Limit/M.hat[i],digits = ndigits, format = "f"), ")", sep = "")          
      }
      h.coef <- reg.coefs <- beta.ROC[names(beta.ROC) == "h"]
         ci[i+1] <- paste("(", formatC(sqrt(h.coef * qf((1 - conf)/2, d, h)), digits = ndigits, format = "f"), ",",
                  formatC(sqrt(h.coef * qf((1 + conf)/2, d, h)), digits = ndigits, format = "f"), ")", sep = "")
        pvalue[length(pvalue)+1] <- formatC(2 * min(pf(h.coef^2, h, d), pf(h.coef^2, h, d, lower.tail = FALSE)), digits = ndigits, format = "f")
         ci <- as.matrix(ci, ncol = 1)
         
         p.table <- cbind(coeff, ci, pvalue)
            dimnames(p.table) <- list(names(beta.ROC), c("Estimate",ci.legend, "p-value"))            
}
   } else {
   p.table<-matrix(coeff, ncol=1)   
   dimnames(p.table)<-list(names(beta.ROC), "Estimate")
   
   }
   object$p.table <- p.table
   class(object) <- "summary.ROCreg"
   object
}

