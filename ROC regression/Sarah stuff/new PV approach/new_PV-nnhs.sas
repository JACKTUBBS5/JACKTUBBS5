libname sarah "/folders/myfolders/Sarah_S/SAS stuff/";
data nnhs; set sarah.nnhs2; nob=_N_; *adult=(adultobi='normal');gender = 2 - gender; drop id ear sitenum ;run;
proc print data=nnhs; where nob < 20; run;

proc sort data=nnhs; by gender; run;

proc freq data=nnhs;
   tables d*gender;
run;

/*

proc surveyselect data=nnhs m=srs n=500
                  seed=1953 out=new_nnhs;
*   size 200;
  * strata d ;
run;


proc freq data=new_nnhs;
   tables d*gender;
run;
*/
proc sort data=nnhs; by d; run;
proc means data=nnhs; var y1 currage; by d; run;
data nnhs; set nnhs; censor=1; y1 = y1 + 40; if y1 < 0 then y1=0; run;

title3 'Frequentist AUC';
proc logistic data=nnhs plots=roc;
*ods select ROCAssociation;
   model d(event='0') = y1 / nofit;
   roc  y1;
run;

proc sgplot data=nnhs;
 histogram y1 /group=d;
 density y1 / group=d;
 run;

*++++++++++++++++++++++++++++++++++++++++++++++++;



proc phreg data=nnhs plots(overlay)=survival; where d=0;
   model y1*censor(0) = currage gender;
   output out=temp2 survival=surv;* xbeta=xbeta logsurv = lsurv;
 *  baseline covariates=in_data out=sur_pred survival=surv;
run; 


data temp3; merge temp2 nnhs; by d; keep d currage gender y1 surv; run;
proc sort data=temp3; by y1; run;
*proc print data=temp3; run;

 proc iml;
 
    use temp3; 
    read all into tc;
    currage=tc[,1];
    gender=tc[,2];
    d =tc[,3];
    y1 = tc[,4];
    surv = tc[,5];
    n    = nrow(tc);
 
    do j = 1 to n;
    if d[j]=1 then surv[j]=surv[j-1];
    end;
new = d||currage||gender||y1||surv;
*print new;
create f from new;
  append from new;
*/
run;
quit; 

data f; set f; d=col1; currage=col2; gender=col3; y1=col4; PV=col5; drop col1-col5;

proc sgplot data=f;
 histogram PV /group=d;
 density PV / group=d;
 run;


proc glimmix data=f;* plot=predpplot;where d=1;  
     model PV = currage gender /dist=beta solution ; 
run; 


proc glimmix data=f;* plot=predpplot;where d=1; 
ods select ParameterEstimates; 
ods output ParameterEstimates=parms; 
    model PV = currage gender /dist=beta solution ; 
run; 

proc transpose data=parms out = mn; 
run; 

data mn (keep= Intercept e_A1 e_A2 scale 
               mu omega tau roc s currage gender);
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
        do gender = 0 to 1 by 1;
        do currage = 30 to 50 by 10;
 		mu= 1/(1+exp(-Intercept - e_A1*currage - e_A2*gender));  
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

data temp4; set mn; if gender=1; run;

proc sgplot data=temp4;
*panelby x;
 series y=roc x=s/ group=currage;
 run;
 
data temp5; set mn; if currage = 40; run; 
 
 proc sgplot data=temp5;
*panelby x;
 series y=roc x=s/ group=gender;
 run;
 
 
 %macro loop(dsn= , currage=, roc = , title= );
 data temp; set &dsn; where currage=&currage; keep currage s &roc;run;
 title &title;
 proc sgplot data=temp;
 series y=&roc x=s;
 run;
 *proc print data=temp;
 proc iml;
	use temp;
    read all into data;
    x   = data[1,1];
    t   = data[,2];
    roc = data[,3];
    np 	   = nrow(data);
    dt = t;
*    print np dt x t roc;
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
data tempx; set temp4; keep currage s roc; run;
%loop(dsn = tempx, currage = 30, roc = ROC, title = 'Males Age 30');
%loop(dsn = tempx, currage = 40, roc = ROC, title = 'Males Age 40');
%loop(dsn = tempx, currage = 50, roc = ROC, title = 'Males Age 50');

run;
quit; 
