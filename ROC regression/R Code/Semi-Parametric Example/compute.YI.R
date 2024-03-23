compute.YI <-
function(data.ROC, set.p, ROCreg, TH = FALSE) {

tYI <- function(t, x, ROCreg) {
obj <- compute.ROC(x,t,ROCreg)
 as.numeric(abs(obj[length(obj)] - t))
}

compute.YI.sub <- function(data, set.p, ROCreg, TH) {
res <- t(sapply(1:dim(data)[1], 
function(i, data, range.FP, ROCreg) {
x <- data[i, , drop = FALSE]
unlist(optimize(tYI, interval = c(0, 1), maximum = TRUE, x = x, ROCreg=ROCreg))
}, data = data, range.FP = range(set.p), ROCreg=ROCreg))
YI <- data.frame(YI = res[, "objective"])  
if(TH)
 YI$TH <- compute.TH(data, p.opt = res[ , "maximum"], ROCreg) 
YI
}  
YI <- switch(ROCreg$method,
    "NM" = ,
"SM" =,
"PROCGLM" =,
"SROCGLM" = compute.YI.sub(data=data.ROC, set.p=set.p, ROCreg=ROCreg, TH = TH))
}

