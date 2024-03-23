title 'Alvarez Comparative Simulation';
title2 'Alonzo-Pepe ROC Regression';
title3 'Input Parameters';

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


proc quantreg data=binorm; where LA=0;
ods select ParameterEstimates;
*class gender;
model y = x/*age gender*/ /quantile=(.01 to .99 by .01) nosummary ;
*model y = x /*age gender*/ /quantile=(.01 .05 .10 .15 .20 .25 .30 .35 .40 .45 .50  .55 .60 .65 .70 .75 .80 .85 .90 .95 .99) nosummary;
 output out=e p=p ;
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
  t = j(np,nq,0);
    do j = 1 to nq;
  *      qq[j] = quantile("normal", quant[j]);
      do i = 1 to np;
        dat[i,j] = (y[i] >= quant[j]);
        t[i,j] = qq[j];
      end;  
    end;
* print sdat;   

new = dat||t||data[,2:3];;
*print new;
*varname={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21, age, bmiz, abmiz, nob, adult, t} ;
create f from new;*[colname= varname];
  append from new;
  
run;
quit; 



data temp; set f;
y=col199; 
x=col200;

array U[99] COL1-COL99; 
array phi[99] COL100 - COL198; 
	do i = 1 to 99;
	uit = U[i];
	phi_inv = phi[i];
	output; end;
	
keep y uit phi_inv x;
run;
/* proc print data = temp; run; */
/*
data temp; set temp; uit=1-uit; run;
*/


proc probit data=temp plot=predpplot;
ods select ParameterEstimates;
   ods output ParameterEstimates=parms;
model uit = phi_inv x ;
run;

proc transpose data=parms out=out_parms; run;

data out_parms; set out_parms;
	if _NAME_ = 'Estimate';run;

data temp1; set out_parms; 
	alpha_hat = col1;
	beta_hat = col2;
	theta = col3;

do x = 0 to 1 by .2;
	do s = 0.001 to 0.999 by .005;
		quant = quantile('normal', s, 0, 1);
		ROC = cdf('normal', +1*alpha_hat + beta_hat*quant + theta*x, 0,1);
		output;
end;end;

keep alpha_hat beta_hat s theta x ROC ;
run;
*proc print data = temp1; *run;

proc sgpanel data=temp1;
panelby x;
 series y=ROC x=s;
 run;

%macro loop(dsn= , cov=, title= );
 data temp; set &dsn; where x=&cov; keep x s roc;run;
 title &title;
 *proc sgplot data=temp;
 *series y=roc x=s;
 *run;
  
 proc iml;
	use temp;
    read all into data;
    x   = data[1,1];
    t   = data[,2];
    roc = data[,3];
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

%loop(dsn = temp1, cov = 0, title = 'x = 0');
/*%loop(dsn = temp1, cov = 0.1, title = 'x = 0.1');*/
%loop(dsn = temp1, cov = 0.2, title = 'x = 0.2');
*%loop(dsn = temp1, cov = 0.3, title = 'x = 0.3'); *error -- can't find 0.3;
%loop(dsn = temp1, cov = 0.4, title = 'x = 0.4');
*%loop(dsn = temp1, cov = 0.5, title = 'x = 0.5');
%loop(dsn = temp1, cov = 0.6, title = 'x = 0.6');
*%loop(dsn = temp1, cov = 0.7, title = 'x = 0.7');
%loop(dsn = temp1, cov = 0.8, title = 'x = 0.8'); * error -- can't find 0.8;
%loop(dsn = temp1, cov = 1.0, title = 'x = 1.0');
run;
quit;

quit;
 
