title 'Alvarez Comparative Simulation';
title2 'Alonzo-Pepe ROC Regression';
title3 'Input Parameters';

proc options option = RLANG;
run;
* this is a modified version of Alvarez_sim1.sas;
* corresponding R-file is AP_method_datagen_sas.R;

data inparm; seed0=12345; seed1=67891;

    do n = 1 to 50;
    x0 = ranuni(seed0);
    x1 = ranuni(seed1);

    LA=0;
    m0 = 0.0;
    s0 = 1.0;
	   y0 = m0 + 3*x0 + sqrt(s0)*rannor(seed0);

    m1 = 0;
    s1 = 1.0;
	   y1 = m1 + 4*x1 + sqrt(s1)*rannor(seed1);

	la1A = (m1-m0+(x1))/s1;
    la2A = s0**2/s1**2;
    AUC_true = cdf("normal", (la1A/sqrt(1+la2A)), 0, 1);
    output;
    end;
run;
proc print data=inparm; where n <= 10; run;

data binorm; set inparm;
   y=y0; x=x0; LA=0; output;  *0 for reference;
   y=y1; x=x1; LA=1; output;  *1 for diseased;
   keep LA y x;
run;
title3 'Frequentist AUC';
proc logistic data=binorm plots=roc;
*ods select ROCAssociation;
   model LA(event='0') = y x / nofit;
   roc  y;
run;

*++++++++++++++++++++++++++++++++++++++++++++++++;


proc quantreg data=binorm ; where LA=0;
ods select ParameterEstimates;
*class gender;
model y = x/*age gender*/ /quantile=(.01 to .99 by .01) nosummary;
*model y = x /*age gender*/ /quantile=(.01 .05 .10 .15 .20 .25 .30 .35 .40 .45 .50  .55 .60 .65 .70 .75 .80 .85 .90 .95 .99) nosummary;
 output out=e p=p;
 run;
 
 proc means data=e noprint; var p1-p99; output out=a;
 
  
 data b; set a; where _STAT_='MEAN'; keep p1-p99; run;
 
data c; set binorm; where LA=1; run;

proc iml;
    use b; 
    read all into quant;
    FPR= { 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.20 0.21 0.22 0.23 0.24  0.25 0.26 0.27 0.28 0.29 0.30 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.40 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48  0.49 0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99};
    *FPR = {.01 .05 .10 .15 .20 .25 .30 .35 .40 .45 .50  .55 .60 .65 .70 .75 .80 .85 .90 .95 .99};
    
	qq=quant;
    
    nq     = ncol(quant);
    do j = 1 to nq;
    	qq[j] = quantile("normal", FPR[j]);
    end;
 *   print quant qq;
	use c;
    read all into data;
 * print data;
    y   = data[,2];
    x   = data[,3];
    np 	   = nrow(data);
    
  dat=j(np,nq,0);
  
    do j = 1 to nq;
        do i = 1 to np;
  	      dat[i,j] = (y[i] >= quant[j]);
      end;  
    end;
* print sdat;  

*#############################*
* Calculating placement values
*##############################;

pv = j(np,1,.9999);

  do i = 1 to np;
  
    do j = 1 to nq-1;
        if y[i] > quant[j] & y[i] <= quant[j+1] then pv[i] = 1 - j*.01 ;
        if y[i] < quant[1] then pv[i] = .0001;
	end; 
  end;	
  
 

new = dat||t||data[,2:3]||pv;;
*print new;
*varname={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21, age, bmiz, abmiz, nob, adult, t} ;
create f from new;*[colname= varname];
  append from new;
 
 
 *#################################
 * Calculating Pairwise Differences
 *#################################;
 
/** compute matrix of differences **/
/** 1. use matrices **/ 
start FullMatDiff(x); /** x is a column vector **/
   n = nrow(x);   
   m = shape(x, n, n); /** duplicate data **/
   return( m` - m );
finish;
 
D = FullMatDiff(pv);
*print D[c= 'placementVal' r='placementVal'];
 
 
start TriDiff(x);
   n = nrow(x);   
   diff = j(n, n, .);
   xt = T(x); 
   do i = 1 to n-1;
      cols = i+1:n;
      diff[i,cols] = x[i] - xt[,cols];
   end;
   return( diff );
finish;

PVdiff = TriDiff(pv);
Xdiff = TriDiff(x);
print PVdiff Xdiff;

/*  All pairwise PV differences into a column vector */

r1 = row(PVdiff);
c1 = col(PVdiff);
upperTri = loc(r1 < c1); /* upper tri indices in row major order */
v1 = PVdiff[upperTri];    /* vector contains n*(n-1)/2 upper triangular corr */
*print v1;

/* Converting pairwise PV diff to binary */

BinPV =(v1 <= 0);

/* All covariate differences into a column vector */
r2 = row(Xdiff);
c2 = col(Xdiff);
upperTri2 = loc(r2 < c2); /* upper tri indices in row major order */
v2 = Xdiff[upperTri2];    /* vector contains n*(n-1)/2 upper triangular corr */

/*/  Data set for first genmod call /*/
temp2 = v1||BinPV||v2;;
*print temp2;
*print v1 BinPV v2;
create g from temp2;*[colname= varname];
  append from temp2;
 
/*/ sending the "dat" matrix to R for second glm call /*/
run ExportMatrixToR(dat, "data1"); 
run ExportMatrixToR(BinPV, "BinPV");
run ExportMatrixtoR(v2, "xDiff");
run ExportMatrixtoR(x, "xVec");  /* sending the covariate vector*/
run ExportMatrixtoR(FPR, "FPRvec");

submit / R;
   	library(quantreg)
	library(foreign)
	library(gtools)
   probitMod1 <- glm(BinPV ~ xDiff - 1 , family = binomial(link = "probit")) # No intercept
   beta.hat <- probitMod1$coefficients
   xPbeta <- beta.hat*xVec
   h0 <- apply((1 - data1),2,function(s){
   glm(s ~ + offset(xPbeta), family = binomial(link = "probit"))$coefficient[1]})
   
    cov.data = seq(0,1, by = .2)
	tVec = seq(0.01,.99, by = .005)
	xseq = seq(0,1, by = .2)
	
   		xPbetaROC <- cov.data * beta.hat
  
  		h <- approxfun(FPRvec, h0)(tVec)  
  		ROC <- sapply(1:length(tVec), function(s) pnorm(h[s] + 
                                                   xPbetaROC))
  		
	newROC <- t(ROC)
	newROC.long <- c(newROC)
	t.long <- rep(tVec, length(xseq))
	x.long <- rep(xseq, each = length(tVec))
	
	ROCdata <- as.matrix(data.frame(cbind(ROC = newROC.long, t = t.long, x = x.long, factor.x = as.factor(x.long))))
	
	
endsubmit;
 
run ImportMatrixFromR(baselineH, "h0");
run ImportMatrixFromR(approxH, "h");
run ImportMatrixFromR(betaHat, "beta.hat");

run ImportMatrixFromR(ROC, "ROCdata");
*print baselineH, betaHat, ROC;

ROCdata = ROC;
create ROCdata from ROC;
append from ROC;

run;
quit; 



/*/  #########  End of Proc IML ########################### /*/

data temp2; set g; BinPV=col2; xDiff=col3;drop col1-col3; run;
proc sort data = temp2 out = sorted;
by BinPV;
run;

/*/  #########  Estimating the Betas -- did this above in IML, but can verify ################### /*/
proc probit data=sorted plot=predpplot ;
ods select ParameterEstimates;
   ods output ParameterEstimates=parms;
model BinPV = xDiff /noint ;
run;


/*/################# plotting the ROC for various values of the covariate x ####################### /*/
data temp3; set ROCdata; 
	ROC = col1;
	tvec = col2;
	x = col3;

keep ROC tvec x ;
run;


proc sgpanel data=temp3;
panelby x;
 series y=ROC x=tvec;
 run;
 
 /*/#################### calculating the AUC  ######################################/*/
 %macro loop(dsn= , cov=, title= );
 data temp; set &dsn; where x=&cov; keep x tvec ROC;run;
 title &title;
 proc sgplot data=temp;
 series y=roc x=tvec;
 run;
  
 proc iml;
	use temp;
    read all into data;
    x   = data[1,3];
    t   = data[,2];
    roc = data[,1];
    np 	   = nrow(data);
    dt = t;
    auc=roc; meanroc=roc; pAUC=AUC;
   dt[1] = t[1];
   meanROC[1] = ROC[1]/2;
   pAUC[1] = dt[1]*meanROC[1];
   AUC[1]= pAUC[1];

do i=2 to np by 1;
   dt[i] = t[i] - t[i-1];
   meanROC[i] = (ROC[i] + ROC[i-1])/2;
   pAUC[i] = dt[i]*meanROC[i];
   AUC[i] = AUC[i-1] + pAUC[i]; 
end;
AUC_PVest = AUC[np];
print x AUC_PVest;
quit;
%mend;
%loop(dsn = temp3, cov = 0.0, title = 'x = 0.0');
%loop(dsn = temp3, cov = 0.2, title = 'x = 0.2');
%loop(dsn = temp3, cov = 0.4, title = 'x = 0.4');
%loop(dsn = temp3, cov = 0.6, title = 'x = 0.6');  /*/ Error /*/
%loop(dsn = temp3, cov = 0.8, title = 'x = 0.8');
%loop(dsn = temp3, cov = 1.0, title = 'x = 1.0');


