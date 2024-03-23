plot.ROCreg <- function(x, data = NULL, step.p=0.02, accuracy = NULL, AUC.calc=c("simpson","integrate"), control = controlROCregData(),ask = TRUE, ...) {
     
   	change.ROC.format <- function(p, ROC) {
      	temp <- reshape(ROC, varying = paste("p", round(p, 3), sep = ""), sep = "",
      	v.names = "ROC", timevar = "p", times = p, idvar = "comb", direction = "long")
      	temp[order(temp$comb),]
	}	
	plot.accuracy <- function(ROC, names.cat, n.cat, n.levels, names.cont, exp.cat, dim.exp.cat, range.marker, accuracy, dots, ask) {
   		
   		if(ask)
   			readline("Press return for next page....")
        
   	    par(mfrow = c(1, 1))            	
		plot(0, 0, xlab = names.cont, ylab = accuracy, xlim = range(ROC[ , names.cont]), ylim = if(accuracy == "TH") range.marker else c(0,1), type = "n")
      	if (n.cat == 0)
	   		lines(ROC[ , names.cont], ROC[ , accuracy])
      	else
      		if (n.cat == 1) {
   	      		for(i in 1:dim.exp.cat)
	         		lines(ROC[ROC[, names.cat] == exp.cat[i, ], names.cont], ROC[ROC[, names.cat] == exp.cat[i, ], accuracy], lty = i)
	      			legend(if(!is.null(dots$pos.legend)) dots$pos.legend else "bottomright", legend = paste(names.cat, "=", exp.cat[, , drop = TRUE]),
               		lty = 1:dim.exp.cat, cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 1, y.intersp = if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 1)
         	} else {
            for (i in 1:dim.exp.cat) {
               	ind <- apply(t(apply(ROC[, names.cat], 1, function(x) x == exp.cat[i, ])), 1, all)
	         	lines(ROC[ind, names.cont], ROC[ind, accuracy], lty = i)
            }
      		legend(if(!is.null(dots$pos.legend)) dots$pos.legend else "bottomright", legend = sapply(1:dim.exp.cat, function(s) paste(paste(names.cat,"=",as.matrix(exp.cat)[s,]), collapse = ", ")),
       			lty = 1:dim.exp.cat, cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 1, y.intersp = if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 1)
       	}
		if (accuracy == "AUC") abline(h = 0.5, col = "grey")
   	}

   	dots <- list(...)   	
   	set.accuracy <- c("AUC", "YI", "TH") 
   	
   	if (is.null(data) | !inherits(data, "predict.ROCreg")) {
   		data <- predict(x, newdata = data, step.p=step.p, accuracy = accuracy, AUC.calc = AUC.calc, control = control)
   	}
   	
   	ROC<-cbind(data$data, data$ROC)
   
   	ind.accuracy <- is.element(set.accuracy, names(data))
   	if(any(ind.accuracy)) 
   		accuracy <- set.accuracy[ind.accuracy]   			
		else
			accuracy=NULL
			
		if(!is.null(accuracy)) {
			for (i in 1:length(accuracy)){
				aux<-names(ROC)
				ROC<-cbind(ROC, data[[accuracy[i]]])
				names(ROC)<-c(aux,accuracy[i]) 
			}
		}
		
   	step.p <- data$step.p
   	p <- seq(0, 1, by = step.p) 
   	n.p <- length(p)
   	range.marker <- if(ind.accuracy[3]) range(x$data[, x$marker])
   	names.cov <- names(ROC[, 1:(dim(ROC)[2] - n.p - sum(ind.accuracy)), drop = FALSE])
   	ind.cat <- unlist(lapply(ROC[ ,names.cov, drop = FALSE], is.factor))
   	names(ind.cat) <- names.cov     
   	names.cont <- names.cov[!ind.cat]        
   	n.cont <- length(names.cont)
   	n.cat <- sum(ind.cat)
   	names.cat <- if(n.cat > 0) names.cov[ind.cat]
   	if (n.cont > 1) ROC[, names.cont] <- apply(round(ROC[ , names.cont, drop = FALSE], 3), 2, factor)
   	if (n.cat > 0) {
 		exp.cat <- unique(ROC[, names.cat, drop = FALSE])
      	exp.cat.matrix <- as.matrix(exp.cat)
      	dim.exp.cat <- nrow(exp.cat)
      	levels.cat <- if(n.cat > 0) lapply(ROC[, names.cat, drop = FALSE], levels)          	
      	n.levels <- as.numeric(unlist(lapply(levels.cat, length))) 	      	
      	if(n.cont == 0) {
	   		ROC.long <- change.ROC.format(p, ROC)
	        print(xyplot(as.formula(paste("ROC ~ p |", paste(names.cat, collapse = "+"))),
	        data = ROC.long,
	        xlab = "FPF", ylab = "TPF",
	        ylim = c(-0.1,1.05),
	        strip = strip.custom(strip.names = TRUE, strip.levels = TRUE, sep = " = ",
	    	par.strip.text = list(cex = if(!is.null(dots$cex.par.strip.text)) dots$cex.par.strip.text else 0.75)),
	        panel = function(x, y, subscripts) {
            	panel.xyplot(x, y, type = "l")
				for (i in 1:3)
   		   			if(ind.accuracy[i]) 
                    	ltext(0.99, 0.01 + (if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 0.14) * (i-1),
                    	labels = paste(set.accuracy[i],"=", round(unique(ROC.long[subscripts, set.accuracy[i]]),2)), adj = c(1,0.5),
                    	cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 0.5)
	   		}, page = function(page){readline("Press return for next page....")}))      	 	      	
      	} else {
         	if(n.cont == 1) {
            	if (length(names.cat) == 1) {
           			par(mfrow = c(1, dim.exp.cat))            	
	         		for(i in 1:dim.exp.cat) {
		            	persp(p, ROC[ROC[ , names.cat] == exp.cat[i, ], names.cont], 
	                    t(as.matrix(ROC[ROC[ , names.cat] == exp.cat[i, ], -(c(1:(1 + n.cat), if(!is.null(accuracy)) ncol(ROC):(ncol(ROC) + 1 - length(accuracy))))])),
		   	         	xlab = "FPF", ylab = names.cont, zlab = "TPF",   	  	
		           	   	sub = exp.cat.matrix[i, ],    	  	
		               	theta = if (!is.null(dots$theta))dots$theta else 20,
			         	phi = if (!is.null(dots$phi))dots$phi else 30,
		          	   	col = if(!is.null(dots$col))dots$col else "white",
		           	   	shade = if(!is.null(dots$shade))dots$shade else 0.5, ticktype = "detailed",
	            	   	cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.sub = dots$cex.sub, cex = dots$cex)	            	            
		         	}
	         		if(any(ind.accuracy))
	            		for(i in (1:3)[ind.accuracy])
                  			plot.accuracy(ROC, names.cat, n.cat, n.levels, names.cont, exp.cat, dim.exp.cat, range.marker, set.accuracy[i], dots, ask)
     			} else {
	           		par(mfrow = n.levels[1:2])   
	               	for (i in 1:(dim.exp.cat/prod(n.levels[1:2]))) {
	               		if(i > 1) {
	               			if(ask)
   								readline("Press return for next page....")
        				}
	                  	k <- 0
	                  	for (j in 1:(n.levels[1]*n.levels[2])) { 
		                	ind <- apply(t(apply(ROC[,names.cat], 1, function(x) x == exp.cat[(i-1)*prod(n.levels[1:2]) + j,])), 1, all)        	
			               	persp(p, ROC[ind, names.cont],
		                 	t(as.matrix(ROC[ind, -(c(1:(n.cont + n.cat),if(!is.null(accuracy)) ncol(ROC):(ncol(ROC) + 1 - length(accuracy))))])),
		   	      	      	xlab = "FPF", ylab = names.cont, zlab="TPF",
		                    sub = paste(paste(names.cat, "=", c(exp.cat.matrix[j,1:2], exp.cat.matrix[1+(i-1)*prod(n.levels[1:2]),-(1:2)])),collapse = ", "),
			                theta = if (!is.null(dots$theta))dots$theta else 20,
			   	            phi = if (!is.null(dots$phi))dots$phi else 30,
			                col = if(!is.null(dots$col))dots$col else "white",
			                shade = if(!is.null(dots$shade))dots$shade else 0.5, ticktype = "detailed",
			                cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.sub = dots$cex.sub,cex = dots$cex)              			
	          			}
	       	  		}
	      			if(any(ind.accuracy))
	         			for(i in (1:3)[ind.accuracy])
                  			plot.accuracy(ROC, names.cat, n.cat, n.levels, names.cont, exp.cat, dim.exp.cat, range.marker, set.accuracy[i], dots, ask)
   				}
			} else {
   				ROC.long <- change.ROC.format(p, ROC)
	            cat.cont <- vector("list", dim.exp.cat)                  	
	            for(i in 1:dim.exp.cat) {
	         		cat.cont[[i]] <- vector("list", n.cont)
	     			for(j in 1:n.cont) {
	     	      		ind <- t(apply(ROC[ , names.cat, drop = F], 1, function(x) x == exp.cat[i, ]))         			
	     	      		if(dim(ind)[1] == 1) ind <- t(ind)   			
	       	      		cat.cont[[i]][[j]] <- unique(ROC[apply(ind, 1, all), names.cont[j]])
	     	   		}
	            }         
            	n.comb <- c(0, as.numeric(cumsum(unlist(lapply(cat.cont, function(x)cumprod(lapply(x, length))[n.cont]))))*n.p)         	
            	ROC.long <- change.ROC.format(p, ROC)
            	for (i in 1:dim.exp.cat) {
            	
            		if(i > 1 && ask) {
            			readline("Press return for next page....")            		
            		}            	 
           			print(xyplot(as.formula(paste("ROC ~ p |", paste(names.cont, collapse = "+"))),
              		data = ROC.long,
	              	ylim = c(-0.1,1.05),
	              	xlab = "FPF", ylab = "TPF",
	               	subset = (1 + n.comb[i]):n.comb[i + 1],
	               	strip = strip.custom(style = 3, strip.names = TRUE, strip.levels = TRUE, sep = " = ",
	            	par.strip.text = list(cex = if(!is.null(dots$par.strip.text)) dots$par.strip.text else 0.75)),
	              	panel = function(x, y, subscripts) {
				   	panel.xyplot(x, y, type = "l")
				   		for (i in 1:3)
	   			      		if(ind.accuracy[i]) 
	                           	ltext(0.99, 0.01 + (if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 0.14) * (i-1),
	                           	labels = paste(set.accuracy[i],"=", round(unique(ROC.long[subscripts, set.accuracy[i]]),2)), adj = c(1, 0.5),
	                           	cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 0.5)
			      	},
               		main = paste(names.cat, "=", exp.cat.matrix[i, ]),
               		page = function(page){readline("Press return for next page....")}))
      			}
			}
 		}
	} else {   
		if(n.cont == 1) {  	
    		persp(p, ROC[ , names.cont], t(as.matrix(ROC[ ,-(c(1:(1 + n.cat), if(!is.null(accuracy)) ncol(ROC):(ncol(ROC) + 1 - length(accuracy))))])),
           	xlab = "FPF", ylab = names.cont, zlab = "TPF",	
          	theta = if (!is.null(dots$theta))dots$theta else 20,
           	phi = if (!is.null(dots$phi))dots$phi else 30,
           	col = if(!is.null(dots$col))dots$col else "white",
           	shade = if(!is.null(dots$shade))dots$shade else 0.5, ticktype = "detailed",
           	cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.sub = dots$cex.sub,cex = dots$cex)
	      	if(any(ind.accuracy))
	      		for(i in (1:3)[ind.accuracy])
               		plot.accuracy(ROC, names.cat, n.cat, n.levels, names.cont, exp.cat, dim.exp.cat, range.marker, set.accuracy[i], dots, ask)
   	   } else {
            ROC.long <- change.ROC.format(p, ROC)
            print(xyplot(as.formula(paste("ROC ~ p |", paste(names.cont, collapse = "+"))),
          	data = ROC.long,
           	ylim = c(-0.1,1.05),
           	xlab = "FPF", ylab = "TPF",
           	strip = strip.custom(style = 1, strip.names = TRUE, strip.levels = TRUE, sep = " = ",
	        par.strip.text = list(cex = if(!is.null(dots$par.strip.text)) dots$par.strip.text else 0.75)),
            panel = function(x, y, subscripts) {
  		    	panel.xyplot(x, y, type = "l")
				for (i in 1:3)
   					if(ind.accuracy[i]) 
                		ltext(0.99, 0.01 + (if(!is.null(dots$y.intersp.legend)) dots$y.intersp.legend else 0.14) * (i-1),
                    	labels = paste(set.accuracy[i],"=", round(unique(ROC.long[subscripts, set.accuracy[i]]),2)), adj = c(1,0.5),
                        cex = if(!is.null(dots$cex.legend)) dots$cex.legend else 0.5)
	         },
	         page = function(page){readline("Press return for next page....")}
	      	))
		}                  
   	}
   	invisible(data)
}