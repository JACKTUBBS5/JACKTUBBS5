compute.EQ <-
function(data.ROC, set.p, ROCreg, TH = FALSE) {
tEQ <- function(t, x, ROCreg) {
obj <- compute.ROC(x,t,ROCreg)
    as.numeric(abs(obj[length(obj)] - 1+t))
 }
compute.EQ.sub <- function(data, set.p, ROCreg, TH) {
res <- t(sapply(1:dim(data)[1], 
function(i, data, range.FP, ROCreg=ROCreg) {
x <- data[i, , drop = FALSE]
unlist(optimize(tEQ, interval = c(0, 1), maximum = FALSE, x = x, ROCreg=ROCreg))
}, data = data, range.FP = range(set.p), ROCreg=ROCreg))

EQ <- data.frame(EQ = 1-res[, "minimum"])  

if(TH)
EQ$TH <- compute.TH(data, p.opt = res[ , "minimum"], ROCreg)

EQ
}
  
   EQ <- switch(ROCreg$method,
    "NM" = ,
"SM" =,
"PROCGLM" =,
"SROCGLM" = compute.EQ.sub(data=data.ROC, set.p=set.p, ROCreg=ROCreg, TH = TH))

}

