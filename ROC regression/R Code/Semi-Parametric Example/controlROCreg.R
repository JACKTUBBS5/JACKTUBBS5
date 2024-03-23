controlROCreg <-
function(
   step.h = 0.02,     
   card.T = 10,         
   ROC.model = c("binormal", "logistic"),   
   nboot = 500,
   resample.m =c("coutcome","ncoutcome"),
   FPFint = NULL,  
   est.surv=c("normal","empirical")
   )
   list(step.h = step.h, card.T = card.T, ROC.model = match.arg(ROC.model), nboot = nboot, resample.m = match.arg(resample.m),FPFint = FPFint, est.surv=match.arg(est.surv))

