/*
 This SAS macro is used fit the PH model using the power prior via MCMC. 
 The macro requires the following parameters:
  
  MACRO VARIABLES REQUIRED:
   [1] dsNew           = dataset of simulated new trials after reduction to sufficient statistics; This dataset needs to have been 
                         created using the %COLLAPSE macro or needs to be in the same structure as the dataset produced by that macro;
   [2] dsHist          = historical dataset after reduction to sufficient statistics; This dataset needs to have been 
                         created using the %COLLAPSE macro or needs to be in the same structure as the dataset produced by that macro;
                         If there is no historical dataset, this dataset should simply have zero observations;
   [3] noTreatment     = indicator for whether or not the treatment effect covariate is always zero;
   [4] a0              = value of a0 for the power prior;
   [5] dsOutParmEst    = dataset storing results from PROC MCMC;
   [6] mcmcOptions     = PROC statement options for PROC MCMC; This option can be used to set the MCMC seed;
   [7] lambdaPrior     = prior for baseline hazard parameters using PROC MCMC syntax;
   [8] regPrior        = prior for the hazard ratio regression parameters using PROC MCMC syntax;
   [9] sampMethod      = type of MCMC sampler to use (should be a valid option for the PARM statement in PROC MCMC);
*/

%macro fit_mcmc;

 proc sql noprint;
  select max(simStudy) into :numSimStudies
  from &dsNew.;
 quit;

proc sql noprint;
  select max(riskStratum) into :numRiskStratum
  from &dsNew.;
 quit;
 
data _null_;
 set &dsNew.(obs=1);
 array x(*) x:;
 
 call symput('p0',strip(put(dim(x),best.)));
run;
 
 

data __tempHist__;
 set &dsHist.(in=a);

  do simStudy = 1 to &numSimStudies.;
   output;
  end;
run; 

data __tempData__;
 set &dsNew. __tempHist__(in=a);

 if a then study = 0; else study  = 1;
 if a then wgt = &a0.; else wgt = 1.0;
run;

proc sort data = __tempData__ out = __hazCounter__(keep = riskStratum interval) nodupkey;
 by riskStratum interval;
run;

data __hazCounter__;
 set __hazCounter__;
  
   
  retain haz_id 0;
  haz_id+1;
  
  length Parameter $25.;
  PARAMETER = "LAMBDA"||strip(put(haz_id,best.));
run;

proc SQL noprint undo_policy=none;
 create table __tempData__ as 
 select a.*,b.haz_id
 from __tempData__ as a left join __hazCounter__ as b
  on a.riskStratum=b.riskStratum and a.interval=b.interval
  order by simStudy, riskStratum, interval;
quit;

proc sql noprint;
 select max(haz_id) into :numHazComp from __tempData__;
quit;
%let numHazComp = %sysfunc(compress(&numHazComp.));

sasfile __tempData__ load;
ods _all_ close;
ods output PostSummaries = _summary PostIntervals  = _intervals;
proc mcmc data = __tempData__ plots=(none) 
          monitor = (lambda1-lambda&numHazComp. %if &noTreatment. ^=1 %then %do; gamma HR postProb %end;  %do p = 1 %to %eval(&p0.-1); beta&p. %end; ) 
		  statistics = ALL &mcmcOptions. propcov=Quanew;
 by simStudy;

 array lambda[&numHazComp.];
 array logLambda[&numHazComp.];
 
 array x[&p0.] x1-x&p0.;
 array beta[&p0.] gamma %do p = 1 %to %eval(&p0.-1); beta&p. %end;;
 
 beginnodata;

  do k = 1 to dim(lambda);
   logLambda[k] = log(lambda[k]);
  end;
  
   %if &noTreatment. = 1 %then %do; gamma = 0; %end;
   HR = exp(gamma);
   postProb = (HR<1.0);
 endnodata;

 parms lambda: / &sampMethod. ;
 parms %if &noTreatment. ^= 1 %then %do; gamma %end; beta: / &sampMethod. ;

 prior lambda: ~ &lambdaPrior.;
 prior %if &noTreatment. ^= 1 %then %do; gamma %end; beta: ~ &regPrior.; 


 
  linPred = 0;
  do p = 1 to &p0.;
    linPred = linPred + x[p]*beta[p];
  end;

  
 

 likelihood = wgt*(
                    numEvents*linPred + numEvents*logLambda[haz_id] - riskTime*lambda[haz_id]*exp(linPred)
                  );
 model general(likelihood);
run;
sasfile __tempData__ close;

data &dsOutParmEst;
  merge _summary _intervals;
  parameter = upcase(Parameter);
 run;


data &dsOutParmEst;
 length Parameter $25.;
  set &dsOutParmEst;
  parameter = upcase(Parameter);
 run;
 
 proc sql noprint undo_policy=none;
  create table &dsOutParmEst as
  select a.*,b.riskStratum,b.interval
  from &dsOutParmEst as a left join __hazCounter__ as b
   on a.parameter = b.parameter;
 quit;
 
 
 data &dsOutParmEst;
  set &dsOutParmEst;
    if riskStratum >. then parameter = 'LAMBDA_'||strip(put(riskStratum,best.))||"_"||strip(put(interval,best.));
    keep simStudy N mean stdDev p50 parameter HPD:;
 run; proc sort; by simStudy parameter; run;


 proc datasets library=work noprint;
  delete __tempData__ __tempHist__ __hazCounter__ _intervals _summary;
 run;
 quit;  

%mend fit_mcmc;