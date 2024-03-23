print.summary.ROCreg <-
function(x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), ...) {   
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
   print(x$p.table, quote = FALSE, justify = "right")   
   cat("\n\n")   
   invisible(x)
}

