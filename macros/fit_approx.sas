/*
 This SAS macro is used fit the PH model using the power prior via weighted maximum likelihood approximation. 
 The macro requires the following parameters:
  
  MACRO VARIABLES REQUIRED:
   [1] dsNew           = dataset of simulated new trials after reduction to sufficient statistics; This dataset needs to have been 
                         created using the %COLLAPSE macro or needs to be in the same structure as the dataset produced by that macro.
   [2] dsHist          = historical dataset after reduction to sufficient statistics; This dataset needs to have been 
                         created using the %COLLAPSE macro or needs to be in the same structure as the dataset produced by that macro.
                         If there is no historical dataset, this dataset should simply have zero observations;
   [3] shareParms      = indicator variable for whether or not nuisance parameters are shared between studies; A value of 0 indicates
                         that nuisance parameters are not shared; A value of 1 indicates they are shared;
   [4] nCov            = number of hazard ratio regression parameters;
   [5] a0              = value of a0 for the power prior;
   [6] dsOutParmEst    = dataset storing parameter estimate results from PROC GENMOD;
   [7] dsOutModFit     = dataset storing model fit statistics from PROC GENMOD;
*/

%macro fit_Approx;

 proc sql noprint;
  select max(simStudy) into :numSimStudies
  from &dsNew.;
 quit;

 proc sql noprint;
  select min(simStudy) into :startSimStudies
  from &dsNew.;
 quit;


proc sql noprint;
  select max(riskStratum) into :numRiskStratum
  from &dsNew.;
 quit; 
 
 
data __tempHist__;
 set &dsHist.(in=a);
 
  %if &shareParms. = 0 %then %do;
   riskStratum = riskStratum + &numRiskStratum.;
   rename %do i = 2 %to &nCov.; x&i. = x%eval(&nCov.+ &i -1) %end;;
  %end;
 
 do simStudy = &startSimStudies. to &numSimStudies.;
  output;
 end;
run; proc sort; by simStudy; run;

data __tempData__;
 set &dsNew. __tempHist__(in=a);
 by simStudy;
 if a then wgt = &a0.; else wgt = 1.0;
 
 
 array x(*) x:;
  do j = 1 to dim(x);
    if missing(x(j)) then x(j) = 0;
  end;
 
 ** Fudge Factors;
  *if riskTime = 0 then riskTime = 1e-7;
  *if numEvents = 0 then numEvents = 1e-7;
  
  logRiskTime = log(riskTime);

 riskInt = strip(put(riskStratum,z3.))||'-'||strip(put(interval,z3.));
run;

 ods select none;
 ods output ParameterEstimates = &dsOutParmEst. Modelfit = &dsOutModFit.;
 proc genmod data = __tempData__;
  by simStudy;
  weight wgt;
  class riskInt(ref='001-001');
  model numEvents = riskInt x: / dist=poisson link=log offset=logRiskTime;
 run;
 quit;

 
 data &dsOutParmEst.;
  length parameter $25.;
  set &dsOutParmEst.;

   parameter = upcase(parameter);
  
  retain int;
  if parameter = 'INTERCEPT' then do; int = estimate; end;
  if find(parameter,'RISKINT','i') then do;
      estimate = exp(estimate + int);
	  hazCounter+1;
	  parameter = 'LAMBDA_'||strip(put(input(scan(Level1,1,'-'),best.),best.))||"_"||strip(put(input(scan(Level1,2,'-'),best.),best.));
  end;
  

  
  if parameter = 'X1' then PARAMETER = 'GAMMA';
  if parameter not in ("SCALE" "INTERCEPT") then output;
   
  if PARAMETER = 'GAMMA' and df > 0;
  
    est = estimate;
    std = stderr;
  
    parameter= 'HR';
    estimate = exp(est);
    stderr   = sqrt(std**2*est**2);
    output;
  
    parameter= 'POSTPROB';
    estimate  = CDF('normal',(log(1.3)-est)/std);
    stderr = .;
    output;
  
 * keep simStudy estimate stdErr parameter;
  
 run; proc sort; by simStudy parameter; run;

 proc datasets library=work noprint;
  delete __tempData__ __tempHist__;
 run;
 quit; 

%mend fit_Approx;