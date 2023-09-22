libname katie  '/home/jacktubbs/my_shared_file_links/jacktubbs/myfolders/ROC research' ;
options nodate center nonumber ps=200 ls=80 formdlim=' ';


/* Simplified LaTeX output that uses plain LaTeX tables  *
ods latex path='/home/jacktubbs/my_shared_file_links/jacktubbs/LaTeX/clean'
 file='borrow_5_sas.tex' 
stylesheet="sas.sty"(url="sas");

/*
http://support.sas.com/rnd/base/ods/odsmarkup/latex.html 
*
ods graphics / reset width=3.5in outputfmt=png
  antialias=on;
*/


title 'Simulation of Historical Borrowing Issues';
title2 'Emulating Thompson .et .al Example';
%let seed=2023;
%let mu0=100; %let sd0=20; %let n=100; 				/*  Parameters for Historical control */;
* case = 1;/*
title3 'Case 1';
%let d_mu = 1;%let d_sig = 1; %let d_n = .8;	    /*  ratio of current to historical   */;

*case = 2;/*
title3 'Case 2';
%let d_mu = 1.15;%let d_sig = 1; %let d_n = .8;	    /*  ratio of current to historical   */;

*case = 3;/*
title3 'Case 3';
%let d_mu = 1;%let d_sig = .6; %let d_n = .8;	    /*  ratio of current to historical   */;

*case = 4;/*
title3 'Case 4';
%let d_mu = 1.1;%let d_sig = .6; %let d_n = .8;	    /*  ratio of current to historical   */;

*case = 5;
title3 'Case 5';
%let d_mu = 1.1;%let d_sig = 1.6; %let d_n = .8;	    /*  ratio of current to historical   */;

*case = 6;/*
%let d_mu = 1.08;%let d_sig = 1.5; %let d_n = 1;	    /*  ratio of current to historical   */;


data D0;			/*  historical control*/; 
LA = 0;   	
label
 	LA = '0 Historical';			
call rannor(seed = &seed,z);
do i = 1 to &n;
ODI = rand('normal', &mu0,&sd0);
output;
end;
drop i seed z;
run;

data D1;			/*  new study control*/; 
LA = 1;   
label
 	LA = '1 Current';				
mu1 = &mu0*&d_mu; sd1 = &sd0*&d_sig; n1 = &n*&d_n;;
do i = 1 to n1;
ODI = rand('normal', mu1,sd1);
output;
end;
drop i;
run;

data borrow ; set D0 D1;
label
 	LA = '0 Historical 1 Current';
run;

proc sort data=borrow; by LA; run;
 
proc sgplot data=borrow;
density ODI/ type=kernel group=LA; run;

/*************************************************************
**************************************************************
*********Computing KS Statistic*******************************
**************************************************************
*************************************************************/
;
proc ttest data=borrow plots=none;
class LA; var ODI;
run;

proc npar1way edf data=borrow;
class LA;
var ODI;
ods select KS2Stats;
run;


/*************************************************************
**************************************************************
*********Using Frequentist Logistic Model*********************
**************************************************************
*************************************************************/
;

title3 'Frequentist AUC';
proc logistic data=borrow plots(only)=roc;
   model LA(event='1') = ODI ;
 *  ROC ODI;
   ods select  Association ROCcurve OddsRatios;
run;


/*************************************************************
**************************************************************
*********Using Beta Regression ROC Model**********************
**************************************************************
*************************************************************/
;

data binorm; set borrow;
*y=ldl;
y=ODI;
  run;

title3 'Placement Value Method with Beta Distribution';
proc lifereg data=binorm plots=probplot; where LA=0;
   model y = / d=weibull;
   ods select ParameterEstimates; 
run;
proc lifereg data=binorm noprint  plots=probplot; where LA=0;
   model y = / d=weibull;
   output out=tempa  xbeta=xbeta cres=cres cdf=cdf;* censored=censored;
      probplot
      ppout;
      inset;
run;
 
data tempb; set tempa; surv=exp(-1*cres); cdf2=1-surv; 
            keep  LA y  surv cdf; run; 

proc sort data=tempb; by LA; run;
proc sort data=binorm; by LA; run;

data temp3; merge tempb binorm; by LA;  keep   y  LA  surv; run;
proc sort data=temp3; by y; run;


proc iml;
 
    use temp3; 
    read all into tc;
    LA=tc[,1];
    y = tc[,2];
    surv = tc[,3];
    n    = nrow(tc);
    la[1]=0;surv[1]=1;
    do j = 2 to n;
    if LA[j]=1 then surv[j]=surv[j-1];
    end;
new = LA||y||surv;
create temp from new;
  append from new;

run;
quit; 

data f1; set temp; LA=col1; y=col2; PV=col3; keep LA y PV;

proc sort data=f1; by la; run;

proc sgplot data=f1;
 density PV / group=LA type=kernel ;
 run;
 
proc sgplot data=f1;
 hbox PV / group=LA;
 run;

*  Beta regression;
proc glimmix data=f1;where LA=1; 
ods select ParameterEstimates; 
ods output ParameterEstimates=parms; 
    model PV =  /dist=beta solution ; 
run; 
 
proc transpose data=parms out = mn; 
run; 


data mn (keep=  mu omega tau roc s  auc youden);  
	set mn; 
if _NAME_='Estimate'   
	then do; 
 		Intercept = col1;  
 		scale = col2;  
 		mu= 1/(1+exp(-Intercept));  
 		omega = mu*scale; 
 		tau = (1 - mu)*scale;
 		auc = tau / (omega + tau); 
 	   do s = 0.001 to 0.999 by .005; 
 	    ROC = cdf("beta",s, omega, tau); 
 	    youden=ROC - s;
 		output;  
	   end; 	 
    end;
run; 

proc means data=mn mean max; var auc youden;* by visit; run;

* ROC for the Beta Regression;
proc sgplot data=mn;
 series y=s x=s;
 series y=roc x=s;
 run;
 
