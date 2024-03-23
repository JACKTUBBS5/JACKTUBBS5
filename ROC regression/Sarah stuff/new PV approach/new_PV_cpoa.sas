title 'Cai-Pepe Placement Value ROC Regression';
title2 'Input Parameters';


libname mfpt "/folders/myfolders/JVZ/MFPT";
data cpao; set mfpt.cpao; nob=_N_; adult=(adultobi='normal');sex=0; if gender='male' then sex=1; drop studyid adultobi ;run;
proc print data=cpao; where nob < 20; run;
data cpao; set cpao; gender=sex; drop sex; run;
proc sort data=cpao; by adult gender; run;

proc freq data=cpao;
   tables adult*gender;
run;

/*Take a random sample of size n */

proc surveyselect data=cpao m=srs n=400
                  seed=1953 out=new_cpao;
run;

proc means data=new_cpao; var age bmiz; run;

proc freq data=new_cpao;
   tables adult*gender;
run;

proc sort data=new_cpao; by adult; run;
proc means data=new_cpao; var bmiz abmiz age; by adult; run;
data new_cpao; set new_cpao; censor=1; bmiz = bmiz + 2; if bmiz < 0 then bmiz=0; *sex=(gender='male'); run;

title3 'Frequentist AUC';
proc logistic data=new_cpao plots=roc;
*ods select ROCAssociation;
   model adult(event='0') = bmiz / nofit;
   roc  bmiz;
run;

proc sgplot data=new_cpao;
 histogram bmiz /group=adult;
 density bmiz / group=adult;
 run;

*++++++++++++++++++++++++++++++++++++++++++++++++;



proc phreg data=new_cpao plots(overlay)=survival; where adult=1;
class gender;
   model bmiz*censor(0) = age gender;
   output out=temp2 survival=surv;* xbeta=xbeta logsurv = lsurv;
 *  baseline covariates=in_data out=sur_pred survival=surv;
run; 


data temp3; merge temp2 new_cpao; by adult; keep adult age gender bmiz abmiz surv; run;
proc sort data=temp3; by bmiz; run;
*proc print data=temp3; run;

 proc iml;
 
    use temp3; 
    read all into tc;
*    age=tc[,1];
*    bmiz=tc[,2];
*    abmiz=tc[,3];
*    gender=tc[,4];
*    adult =tc[,5];
*    surv = tc[,6];
    n    = nrow(tc);
 
    do j = 1 to n;
    if tc[j,4]=0 then tc[j,5]=tc[j-1,5];
    end;
*new = adult||age||bmiz||abmiz||gender||surv;
print tc;
create f from tc;
  append from tc;
*/
run;
quit; 

data f; set f; age=col1; bmiz=col2; abmiz=col3; adult=col4; PV=col5; drop col1-col5;

proc sgplot data=f;
 histogram PV /group=adult;
 density PV / group=adult;
 run;


proc glimmix data=f;* plot=predpplot;where adult=0;  
     model PV = age abmiz /dist=beta solution ; 
run; 


proc glimmix data=f;* plot=predpplot;where adult=0; 
ods select ParameterEstimates; 
ods output ParameterEstimates=parms; 
    model PV = age abmiz /dist=beta solution ; 
run; 

proc transpose data=parms out = mn; 
run; 

data mn (keep= Intercept e_A1 e_A2 scale 
               mu omega tau roc s age abmiz);
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
        do abmiz = 3 to 6 by 3;
        do age = 5 to 15 by 5;
 		mu= 1/(1+exp(-Intercept - e_A1*age - e_A2*abmiz));  
 		omega = mu*scale; 
 		tau = (1 - mu)*scale; 
 	do s = 0.001 to 0.999 by .005; 
 	    ROC = cdf("beta",s, omega, tau); 
 		output;  
	end; 
	end;
	end;
	end;
run; 

data temp4; set mn; if abmiz=3; run;

proc sgplot data=temp4;
*panelby x;
 series y=roc x=s/ group=age;
 run;
 

data temp5; set mn; if age = 10; run; 
 
 proc sgplot data=temp5;
*panelby x;
 series y=roc x=s/ group=abmiz;
 run;
 

  
 %macro loop(dsn= , age= , roc = , title = );
 title &title;
 data tp; set &dsn; where age=&age; run;
 proc sgplot data=tp;
 series y=&roc x=s;
 run;
 
 proc iml;
/*
    use b; 
    read all into quant;
    qq=quant;
    nq     = ncol(quant);
    do j = 1 to nq;
    qq[j] = quantile("normal", quant[j]);
    end;
 *   print quant qq;
*/
	use tp;
    read all into data;
       x   = data[1,1];
       t   = data[,2];
       roc = data[,3];
       np  = nrow(data);
       dt  = t;
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
print &age AUC_PVest;
quit;
%mend;
data tempx; set temp4; keep age s roc; run;
%loop(dsn = tempx, age = 5, roc = ROC, title = 'Abmiz=3 Age 5');
%loop(dsn = tempx, age = 10, roc = ROC, title = 'Abmiz=3 Age 10');
%loop(dsn = tempx, age = 15, roc = ROC, title = 'Abmiz=3 Age 15');

data temp4; set mn; if abmiz=6; run;
data tempx; set temp4; keep age s roc; run;
%loop(dsn = tempx, age = 5, roc = ROC, title = 'Abmiz=6 Age 5');
%loop(dsn = tempx, age = 10, roc = ROC, title = 'Abmiz=6 Age 10');
%loop(dsn = tempx, age = 15, roc = ROC, title = 'Abmiz=6 Age 15');

/*
%loop(dsn = mn2, age = 5, roc = ROC_m, title = 'Males Age 5');
%loop(dsn = mn2, age = 10, roc = ROC_m, title = 'Males Age 10');
%loop(dsn = mn2, age = 5, roc=ROC_f, title = "Females Age 5"); 
%loop(dsn = mn2, age = 10, roc=ROC_f, title = "Females Age 10");
*/
run;

quit;