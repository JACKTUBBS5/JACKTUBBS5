title 'Cai-Pepe Placement Value ROC Regression';
title2 'Input Parameters';


data inparm; /*m0=0; s0=1; alpha=0.55; beta= 0.90; m1 = alpha/beta; s1=1/beta;*/ seed0=12345; seed1=67891;
/*data inpram; m0=0; s0=1; alpha=1.5; beta= 0.85; m1 = alpha/beta; s1=1/beta; seed0=12345; seed1=67891;*/
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
   y=y0; x=x0; LA=0; output;
   y=y1; x=x1; LA=1; output;

	keep LA y x;
run;


title3 'Frequentist AUC';
proc logistic data=binorm plots=roc;
ods select ROCAssociation;
   model LA(event='0') = y x / nofit;
   roc  y;
run;

*++++++++++++++++++++++++++++++++++++++++++++++++;

*ods listing close;
proc quantreg data=binorm ; where LA=0;
ods select ParameterEstimates;
*class gender;
model y = x/*age gender*/ /quantile=(.01 to .99 by .01) nosummary;
*model y = /*age gender*/ /quantile=(.01 .05 .10 .15 .20 .25 .30 .35 .40 .45 .50  .55 .60 .65 .70 .75 .80 .85 .90 .95 .99) nosummary;
 output out=e p=p;
 run;
 
 proc means data=e noprint; var p1-p99; output out=a;
 
  
 data b; set a; where _STAT_='MEAN'; keep p1-p99; run;
 
 data c; set binorm; where LA=1; run;

proc iml;
    use b; 
    read all into quant;
 
    nq     = ncol(quant);
 	use c;
    read all into data;
  
    y   = data[,2];
    x   = data[,3];
    np 	   = nrow(data);

pv = j(np,1,.0001);

  do i = 1 to np;
  
    do j = 1 to nq-1;
        if y[i] > quant[j] & y[i] <= quant[j+1] then pv[i] = 1 - j*.01;
        if y[i] < quant[1] then pv[i] = .9999;
	end; 
  end;	
  

new = y||x||pv;
*print new;
create f from new;
  append from new;
  
run;
quit; 

data temp; set f; y=col1; x=col2; pv=col3; phi_inv=col4 ;drop col1-col4; run;

proc glimmix data=temp;* plot=predpplot;
ods select ParameterEstimates;
*ods output ParameterEstimates=parm_est;
   ods output ParameterEstimates=parms;
model pv =  x/dist=beta solution ;
run;

proc transpose data=parms out = mn;
run;
*proc print data = mn; run;
 
data mn (keep= Intercept e_A1 /*e_A2 e_A3*/ scale mu /* mu0 mu1 /*mu_A1 alpha_A1 beta_A1 mu_A2
               alpha_A2 beta_A2 /*mu_A3 alpha_A3 beta_A3 mu_laser alpha_laser beta_laser*/ 
               omega  tau x roc s); 
	set mn;
if _NAME_='Estimate'  
	then do;
 		Intercept = col1; 
 		e_A1 = col2; 
  		scale = col3;
 *		mu0    = 1/(1+exp(-Intercept));
do x = 0 to 1.0 by .2;
* do x = 0.0 to 0.5 by .1;		
 		mu    = 1/(1+exp(Intercept + e_A1*x));
 		omega = mu*scale;
 		tau = (1 - mu)*scale;
 	do s = 0.001 to 0.999 by .005;
 	    ROC = cdf("beta",s, omega, tau);

 		output; 
 	 end;	
 end;
end; 
run;
*proc print data = mn; *run;

proc sgpanel data=mn;
panelby x;
 series y=roc x=s;
 run;
quit;


data mn2; set mn; keep  x s roc; run;
  
%macro loop(dsn= , cov=, title= );
 data temp; set &dsn; where x=&cov; keep x s roc;run;
 title &title;
 proc sgplot data=temp;
 series y=roc x=s;
 run;
  
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

%loop(dsn = mn2, cov = 0, title = 'x = 0');
*%loop(dsn = mn2, cov = 0.1, title = 'x = 0.1');
%loop(dsn = mn2, cov = 0.2, title = 'x = 0.2');
*%loop(dsn = mn2, cov = 0.3, title = 'x = 0.3'); *error -- can't find 0.3;
%loop(dsn = mn2, cov = 0.4, title = 'x = 0.4');
*%loop(dsn = mn2, cov = 0.5, title = 'x = 0.5'); 

*%loop(dsn = mn2, cov = 0.5, title = 'x = 0.5'); 
%loop(dsn = mn2, cov = 0.6, title = 'x = 0.6'); 
*%loop(dsn = mn2, cov = 0.7, title = 'x = 0.7'); 
%loop(dsn = mn2, cov = 0.8, title = 'x = 0.8'); *error -- can't find 0.8;
*%loop(dsn = mn2, cov = 0.9, title = 'x = 0.9'); *error -- can't find 0.9;
%loop(dsn = mn2, cov = 1, title = 'x = 1'); *error -- can't find 1;
   
run;
quit; 


