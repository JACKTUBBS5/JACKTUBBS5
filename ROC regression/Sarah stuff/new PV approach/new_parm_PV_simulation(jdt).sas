title 'Alonzo-Pepe ROC Regression';
title2 'Input Parameters';

libname sarah  '/folders/myfolders/Sarah_S/SAS stuff' ;
/*
/*data inparm; m0=0; s0=1; alpha=0.55; beta= 0.90; m1 = alpha/beta; s1=1/beta; seed0=12345; seed1=67891;*
data inparm; m0=0; s0=1; alpha=1.5; beta= 0.85; m1 = alpha/beta; s1=1/beta; seed0=12345; seed1=67891;
la1A = (m1-m0)/s1;
la2A = s0**2/s1**2;
AUC_true = cdf("normal", (la1A/sqrt(1+la2A)), 0, 1);
run;
proc print data=inparm; run;
*/	
data inparm; alpha0=0; alpha1=1.0;beta1 = 0.5;beta2 = 0.7; seed0=12345; seed1=67891; run;
data binorm; set inparm;
	do n = 1 to 100;
	LA=0;
	censor=1;
	e=rannor(seed0);
	z1=ranbin(seed0, 1,0.5);
	z2=ranuni(seed0);
	y = 2 + 0.5*z2 + e;
	if y < 0 then y=0.001;
	output; end;
	do n = 1 to 100;
	e=rannor(seed1);
	z1=ranbin(seed1, 1,0.5);
	z2=ranuni(seed1);
	LA =1;
	y = 2 + (alpha0 + beta1*z1 +(beta2 + 0.5*alpha1)*z2 + e)/alpha1;
	if y < 0 then y=0.001;
	output; end;	
	keep LA y z1 z2 censor;
run;
title3 'Frequentist AUC';
proc logistic data=binorm plots=roc;
*ods select ROCAssociation;
   model LA(event='0') = y / nofit;
   roc  y;
run;

proc sgplot data=binorm;
 histogram y /group=LA;
 density y / group=LA;
 run;

data in_data; set binorm; where LA=1; run;

*++++++++++++++++++++++++++++++++++++++++++++++++;

proc lifereg data=binorm plots=probplot; where LA=0;
   model y*censor(0) = z2/ d=weibull;
   output out=tempa  xbeta=xbeta cres=cres cdf=cdf;* censored=censored;
      probplot
      ppout;
*      npintervals=simul ;
      inset;
   run;
 
data tempb; set tempa; surv=exp(-1*cres); cdf2=1-surv; 
            keep la y z2 surv cdf; run; 
*proc print; run;
data temp3; merge tempb binorm; by LA; keep LA y z1 z2 surv; run;
proc sort data=temp3; by y; run;
proc print;
run;
/*
proc phreg data=binorm plots(overlay)=survival; where LA=0;
   model y*censor(0) = z2;
   output out=temp2 survival=surv xbeta=xbeta logsurv = lsurv;
 *  baseline covariates=in_data out=sur_pred survival=surv;
run;   

data temp2; set temp2; e_xbeta=exp(xbeta); s_0  = exp(lsurv/e_xbeta); run;
proc print data=temp2; run;

data temp3; merge temp2 binorm; by LA; keep LA y z1 z2 surv; run;
proc sort data=temp3; by y; run;
proc print data=temp3; run;

proc print data=tempc; run;
proc print data=temp3; run;
*/
 proc iml;
 
    use temp3; 
    read all into tc;
    LA=tc[,1];
    z1=tc[,5];
    z2 =tc[,2];
    y = tc[,3];
    surv = tc[,4];
    n    = nrow(tc);
 
    do j = 1 to n;
    if LA[j]=1 then surv[j]=surv[j-1];
    end;
new = LA||z1||z2||y||surv;
print new;
create f from new;
  append from new;
*/
run;
quit; 

data f; set f; LA=col1; z1=col2; z2=col3; y=col4; PV=col5; keep LA z1 z2 y PV;

proc sgplot data=f;
 histogram PV /group=LA;
 density PV / group=LA;
 run;


proc glimmix data=f;* plot=predpplot;where LA=1;  
     model PV = z1 z2 /dist=beta solution ; 
run; 


proc glimmix data=f;* plot=predpplot;where LA=1; 
ods select ParameterEstimates; 
ods output ParameterEstimates=parms; 
    model PV = z1 z2 /dist=beta solution ; 
run; 
 
proc transpose data=parms out = mn; 
run; 

data mn (keep= Intercept e_A1 e_A2 /*e_A3*/ scale /*mu_A1 alpha_A1 beta_A1 mu_A2 
               alpha_A2 beta_A2 mu_A3 alpha_A3 beta_A3 mu_laser alpha_laser beta_laser*/ mu omega tau roc s z1 z2);  
	set mn; 
if _NAME_='Estimate'   
	then do; 
 		Intercept = col1;  
 		e_A1 = col2;  
 		e_A2 = col3; 
 *		e_A3 = col4; 
 *		e_laser = col5; 
 		scale = col4; 
 *		mu_A1 = 1/(1+exp(-Intercept-e_A1));  
        do z1 = 0 to 1 by 1;
        do z2 = 0.2 to 0.8 by 0.6;
 		mu= 1/(1+exp(-Intercept - e_A1*z1 - e_A2*z2));  
 		omega = mu*scale; 
 		tau = (1 - mu)*scale; 
 	do s = 0.001 to 0.999 by .005; 
 	    ROC = cdf("beta",s, omega, tau); 
  *		mu_A2 = 1/(1+exp(-Intercept-e_A2)); 
 *		mu_A3 = 1/(1+exp(-Intercept-e_A3)); 
 *		mu_laser = 1/(1+exp(-Intercept-e_laser)); 
 *		alpha_A1=mu_A1*scale; 
 *		beta_A1=(1-mu_A1)*scale; 
 *		alpha_A2=mu_A2*scale; 
 *		beta_A2=(1-mu_A2)*scale; 
 *		alpha_A3=mu_A3*scale; 
 *		beta_A3=(1-mu_A3)*scale; 
 *		alpha_laser=mu_laser*scale; 
 *		beta_laser=(1-mu_laser)*scale; 
 		output;  
	end; 
	end;
	end;
	end;
run; 

data temp4; set mn; if z1=1; run;

proc sgplot data=temp4;
*panelby x;
 series y=roc x=s/ group=z2;
 run;
 
data temp5; set mn; if z2 = .2; run; 
 
 proc sgplot data=temp5;
*panelby x;
 series y=roc x=s/ group=z1;
 run;
 
 

run;
quit; 

