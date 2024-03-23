* beta approach to adjusted cpao data;

data cpao; set sasuser.adj_cpao; /* nob=_N_; adult=(adultobi='normal');
gen = (gender = "female"); drop studyid adultobi ;*/ run;
 
*proc print data=cpao; *run;
proc sort data=cpao; by adult id_gender; run;

proc freq data=cpao;
   tables adult*id_gender;
run;


data  new_cpao; set cpao; run;
/*proc surveyselect data=cpao m=srs n=400
                  seed=1953 out=new_cpao;
run;

*/


proc quantreg data=new_cpao ; where adult=1;
class id_gender;
model bmiz = age id_gender /quantile=(.01 to .99 by .01) nosummary;
   /* model bmiz = age id_gender /quantile=(.01 .05 .10 .15 .20 .25 .30 .35 .40 .45 .50  .55 .60 .65 .70 .75 .80 .85 .90 .95 .99) nosummary;
   */
 output out=e p=p;
 run;
 
proc means data=e; var p1-p99; output out=a;
 /*proc means data=e; var p1-p21; output out=a; */
  
data b; set a; where _STAT_='MEAN'; keep p1-p99; run;
 /* data b; set a; where _STAT_='MEAN'; keep p1-p21; run; */
 
 
data c; set new_cpao; where adult=0; keep age bmiz id_gender abmiz adult; run;
/*proc print data = c; run; */


 
  proc iml;
    use b; 
    read all into quant;
    qq=quant;
    nq     = ncol(quant); * 99 cols;

   *print quant;
	use c;
    read all into cpao;
    *print cpao;
    age    = cpao[,1];  
    bmiz   = cpao[,2];
    abmiz  = cpao[,3];
    gender = cpao[,5];
    np 	   = nrow(cpao);  
    pv = j(np,1,.0001);
    
do i = 1 to np;  
  do j = 1 to nq-1;
    if bmiz[i] > quant[j] & bmiz[i] <= quant[j+1] then pv[i] = 1 - .01*j;
    if bmiz[i] < quant[1] then pv[i] = .9999;
  end; 
end;


  
new = pv||bmiz||age||abmiz||gender;
create f from new;*[colname= varname];
  append from new;
  *print new;

run;


data temp; set f; pv = col1; bmiz=col2; age = col3; abmiz = col4; gender = col5; drop col1-col5; run;

/* proc print data = temp; run; */

proc glimmix data=temp;* plot=predpplot;
ods select ParameterEstimates;
   ods output ParameterEstimates=parms;
model pv =  age gender/dist=beta solution ;
run;


	
proc transpose data=parms out = mn;
run;
 proc print data = parms; run;
 proc print data = mn; run;
 	
data mn (keep= Intercept e_A1 /*e_A2 e_A3*/ scale mu_f mu_m /* mu0 mu1 /*mu_A1 alpha_A1 beta_A1 mu_A2
               alpha_A2 beta_A2 /*mu_A3 alpha_A3 beta_A3 mu_laser alpha_laser beta_laser*/ 
               omega_f  tau_f omega_m tau_m age roc_f roc_m s); 
	set mn;

	
if _NAME_='Estimate'  
	then do;
 		Intercept = col1;  * estimate for intercept; 
 		e_A1 = col2;  * estimate for age;
 		e_A2 = col3;  * estimate for gender;	

 		scale = col4; * estimate for scale
 *		mu0    = 1/(1+exp(-Intercept));
        gender=1;
    do age = 5 to 20 by 5;	* for ages 5, 10, 15, 20 we calculate...;	
 		mu_f    = 1/(1+exp(-Intercept-e_A1*age));
 		mu_m    = 1/(1+exp(-Intercept-e_A1*age-e_A2*gender));
 		omega_f = mu_f*scale;
 		tau_f = (1 - mu_f)*scale;
 		omega_m = mu_m*scale;
 		tau_m = (1 - mu_m)*scale;
 	do s = 0.001 to 0.999 by .005;
 	    ROC_f = cdf("beta",s, omega_f, tau_f); * shouldn't we use cdf here?;
 	    ROC_m = cdf("beta",s, omega_m, tau_m);
 
 		output; 
 	 end;	
 end;
end; 
run;


data mn2; set mn; keep  age s roc_m roc_f; run;
  
/* proc print data = mn2; where age = 5; run;  */
%macro loop(dsn= , age= , roc = , title = );
 title &title;
 data tp; set &dsn; where age=&age; run; *restricting by age group;
 proc sgplot data=tp;
 series y=&roc x=s;  * specifying ROC for males or females;
 run;
 
 proc iml;
	use tp;
    read all into data;
       x   = data[1,1];
       t   = data[,2];
       roc = data[,3];  * note that we restricted to ROC_f/ROC_M above;
       np  = nrow(data);
       dt  = t;
   auc=roc; meanroc=roc; pAUC=AUC;
   dt[1] = t[1];  * so in our case t[1] = .001;
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

%loop(dsn = mn2, age = 5, roc = ROC_m, title = 'Males Age 5');
%loop(dsn = mn2, age = 10, roc = ROC_m, title = 'Males Age 10');
%loop(dsn = mn2, age = 15, roc = ROC_m, title = 'Males Age 15');
%loop(dsn = mn2, age = 5, roc=ROC_f, title = "Females Age 5"); 
%loop(dsn = mn2, age = 10, roc=ROC_f, title = "Females Age 10");
%loop(dsn = mn2, age = 15, roc=ROC_f, title = "Females Age 15");
run;
quit;