print.ROCreg <-
function(x, ...) {
cat("\nCall:\n")
   print(x$call)    
   cat("\n\n*************************************************\n")
   method <- switch(x$method,
     "NM" = cat("Normal Method"),
     "SM" = cat("Semipar Method"),
     "PROCGLM" = cat("Parametric ROC-GLM Method"),
     "SROCGLM" = cat("Semiparametric ROC-GLM Method"))
   cat("\n*************************************************\n\n")
   cat("ROC Coefficients:\n")
   cat("----------------\n\n")
   print.default(format(x$ROC.fit$beta.ROC, digits = 4), print.gap = 2, quote = FALSE)
   cat("\n")
   invisible(x)
}

